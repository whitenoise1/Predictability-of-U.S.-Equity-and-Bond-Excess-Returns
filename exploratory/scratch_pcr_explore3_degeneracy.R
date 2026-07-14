# scratch_pcr_explore3_degeneracy.R
# -----------------------------------------------------------------------------
# EXPLORATION 3 — verify the degeneracy claim.
#
# Claim (METHOD_AUDIT_GS_PCA.md follow-up): the shared selector's PCA-grouping +
# max-min DISTANCE screen was designed for COLLINEAR RAW inputs. Fed the already-
# orthogonal PC SCORES (what the Research method / pcr_core hands it), the
# internal prcomp(scale.=TRUE) sees a correlation matrix == I -> flat eigenvalues
# -> all pairwise loading-distances equal -> the distance screen is INERT (its
# picks are tie/PCA-order-driven, not signal-driven). Net: on scores the real
# selection collapses onto the GS supervised t-gate over the top-k variance PCs.
#
# We verify this on the bond_5y h=1 cell by reproducing pcr_core steps 1-2 to get
# the scores Z, then instrumenting exactly what hybrid_feature_selection does
# internally (Step 1 PCA, Step 2 greedy max-min) on Z vs the RAW collinear X.
# -----------------------------------------------------------------------------
source("working_paper/scratch_pcr_props_setup.R")
suppressMessages(library(grid))

cell <- build_cell("bond_5y", 1L)
X <- cell$X_train                       # 375 x 70 raw z-scored series (collinear)
cat(sprintf("\n== bond_5y h=1 | T=%s | X_train %d x %d ==\n",
            cell$T_date, nrow(X), ncol(X)))

# ---- reproduce pcr_core steps 1-2: PIT-standardize -> SVD -> top-k scores Z ----
ctr <- colMeans(X); scl <- apply(X, 2L, sd)
Zs  <- sweep(sweep(X, 2L, ctr, "-"), 2L, scl, "/")
sv  <- svd(Zs); ev <- sv$d^2; cumv <- cumsum(ev) / sum(ev)
k   <- min(which(cumv >= PCR_MIN_VAR)[1L], min(nrow(X) - 1L, ncol(X)))
V   <- sv$v[, seq_len(k), drop = FALSE]
Z   <- Zs %*% V                         # 375 x k  ORTHOGONAL PC scores fed to selector
colnames(Z) <- paste0("PC", seq_len(k))
cat(sprintf("   PCR retains k=%d PCs at >=%.0f%% cum var (pcr_core step 1)\n",
            k, 100 * PCR_MIN_VAR))

# ---- helper: instrument selector Steps 1-2 on a given design matrix ------------
# Returns the diagnostics that decide the distance screen: correlation off-diag,
# internal-PCA eigenvalue spread, n_pcs at min_var, and the greedy max-min trace.
# NB: the OLD shipped defaults (10L / 0.5) are hardcoded here on purpose -- this
# script DEMONSTRATES that degenerate behaviour; PCR_N_GROUPS/PCR_MIN_DISTANCE now
# default to Inf/0 (the fix), so reading the constants would erase the contrast.
instrument <- function(M, min_var = 0.80, n_groups = 10L,
                       min_distance = 0.5) {
  cm <- cor(M); offd <- abs(cm[upper.tri(cm)])
  pr <- prcomp(M, scale. = TRUE, center = TRUE)
  vr <- pr$sdev^2 / sum(pr$sdev^2); cu <- cumsum(vr)
  npc <- which(cu >= min_var)[1L]; if (is.na(npc)) npc <- length(vr)
  L <- pr$rotation[, seq_len(npc), drop = FALSE]          # Step-2 loadings
  D <- as.matrix(dist(L)); od <- D[upper.tri(D)]          # pairwise feature distances
  # greedy max-min selection trace (Step 2), recording the max-min gap each step
  feats <- colnames(M); sel <- names(which.max(abs(L[, 1])))
  rem <- setdiff(feats, sel); gaps <- numeric(0)
  while (length(sel) < min(n_groups, length(feats)) && length(rem) > 0) {
    md <- vapply(rem, function(f) min(D[f, sel]), numeric(1))
    g  <- max(md); gaps <- c(gaps, g)
    if (g < min_distance && length(sel) >= 2L) break
    nf <- names(which.max(md)); sel <- c(sel, nf); rem <- setdiff(rem, nf)
  }
  list(offdiag = offd, eig = vr, npc = npc, cumvar_at_npc = cu[npc],
       dist = od, sel = sel, gaps = gaps, lam_ratio = vr[1] / vr[length(vr)])
}

cv <- function(x) sd(x) / mean(x)                         # coefficient of variation
rawI <- instrument(X)                                     # raw collinear input
zI   <- instrument(Z)                                     # orthogonal scores input

# ---- report table --------------------------------------------------------------
fmt <- function(a, b, lab, dgt = 3) cat(sprintf("  %-34s %12s %12s\n", lab,
            formatC(a, format = "g", digits = dgt), formatC(b, format = "g", digits = dgt)))
cat("\n-- selector internals: RAW collinear X  vs  orthogonal scores Z --\n")
cat(sprintf("  %-34s %12s %12s\n", "", "RAW X", "scores Z"))
fmt(ncol(X), ncol(Z),                       "candidate columns")
fmt(mean(rawI$offdiag), mean(zI$offdiag),   "mean |corr| off-diagonal")
fmt(max(rawI$offdiag),  max(zI$offdiag),    "max  |corr| off-diagonal")
fmt(rawI$lam_ratio, zI$lam_ratio,           "internal-PCA eig ratio lam1/lamL")
fmt(cv(rawI$eig), cv(zI$eig),               "eigenvalue spread (CV)")
fmt(rawI$npc, zI$npc,                        "n_pcs at 80% var (Step 1)")
fmt(cv(rawI$dist), cv(zI$dist),             "pairwise loading-dist (CV)  <= INERT iff ~0")
fmt(sd(rawI$gaps), sd(zI$gaps),             "greedy max-min gap SD across steps")

# Is the distance pick on Z just the PCA order (PC1..PCn_groups)? -> screen inert.
z_first_n  <- paste0("PC", seq_len(length(zI$sel)))
cat(sprintf("\n  Z distance-survivors == first %d PCs in PCA order?  %s\n",
            length(zI$sel), identical(sort(zI$sel), sort(z_first_n))))
cat(sprintf("  Z survivors: %s\n", paste(zI$sel, collapse = ", ")))

# ---- confirm via the ACTUAL shared selector (byte-identical file) --------------
# Run hybrid_feature_selection on Z three times under COLUMN REORDERINGS. If the
# distance screen carried signal, the survivor SET would be reorder-invariant.
# If it is inert (tie/order-driven), reordering changes which survive.
set.seed(1)
run_sel <- function(perm) {
  Zp <- Z[, perm, drop = FALSE]
  df <- as.data.frame(Zp); df$y <- cell$y_train
  invisible(capture.output(
    hs <- hybrid_feature_selection(df, "y", colnames(Zp),
            n_groups = min(10L, k), min_var_explained = PCR_MIN_VAR,
            min_distance = 0.5, min_t_value = 0, verbose = FALSE)))   # OLD defaults (the broken case)
  sort(hs$gs_result$summary_stats$Predictor)              # distance-survivors (pre t-gate)
}
ord  <- run_sel(seq_len(k))
perm1 <- run_sel(sample(k)); perm2 <- run_sel(sample(k))
cat(sprintf("\n  Shared-selector distance-survivors reorder-invariant?  %s\n",
            identical(ord, perm1) && identical(ord, perm2)))
cat(sprintf("    natural order : %s\n", paste(ord,  collapse = ", ")))
cat(sprintf("    permutation 1 : %s\n", paste(perm1, collapse = ", ")))

# ---- consequence: does the INERT screen discard PREDICTIVE PCs? ----------------
# n_groups=10 < k=16, so the distance screen drops 6 of the 16 retained PCs BEFORE
# the supervised t-gate sees them. If which-6 is arbitrary (above) AND some dropped
# PC carries target signal, the screen is silently discarding predictivity.
uni_t <- function(z) { f <- summary(lm(cell$y_train ~ z))$coefficients
                       if (nrow(f) < 2) NA_real_ else f[2, "t value"] }
t_all   <- vapply(seq_len(k), function(j) uni_t(Z[, j]), numeric(1))
names(t_all) <- colnames(Z)
dropped <- setdiff(colnames(Z), ord)                     # 6 PCs the screen discarded
cat(sprintf("\n  univariate |t| of the %d screen-DROPPED PCs (vs target):\n", length(dropped)))
cat(sprintf("    %s\n", paste(sprintf("%s=%.2f", dropped, abs(t_all[dropped])), collapse = "  ")))
cat(sprintf("    -> dropped PCs with |t|>=1.645 (signal discarded arbitrarily): %s\n",
            paste(dropped[abs(t_all[dropped]) >= 1.645], collapse = ", ")))

pr <- pcr_core(X, cell$y_train, cell$x_T)
cat(sprintf("  pcr_core on this cell: k=%d retained -> %d survive distance screen -> %d pass t-gate\n",
            pr$k, length(ord), pr$n_selected))
cat(sprintf("    final selected PCs: %s\n", paste(pr$selected_pcs, collapse = ", ")))

# ---- figure: scree (concentrated vs flat) + loading-distance spread ------------
png("working_paper/figures/scratch_pcr_degeneracy.png", 1100, 480, res = 130)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
matplot(cbind(sort(rawI$eig, TRUE), c(sort(zI$eig, TRUE), rep(NA, length(rawI$eig) - length(zI$eig)))),
        type = "b", pch = c(16, 1), lty = 1, col = c("firebrick", "steelblue"),
        xlab = "component", ylab = "variance share",
        main = "Internal-PCA scree")
legend("topright", c("RAW X (concentrated)", "scores Z (flat -> degenerate)"),
       col = c("firebrick", "steelblue"), pch = c(16, 1), lty = 1, bty = "n", cex = 0.85)
plot(density(rawI$dist), col = "firebrick", lwd = 2, xlab = "pairwise loading-distance",
     main = "Distance screen: dispersed vs inert",
     xlim = range(c(rawI$dist, zI$dist)))
lines(density(zI$dist), col = "steelblue", lwd = 2)
legend("topright", c("RAW X (informative)", "scores Z (all equal)"),
       col = c("firebrick", "steelblue"), lwd = 2, bty = "n", cex = 0.85)
par(op); invisible(dev.off())
cat("\n[fig] working_paper/figures/scratch_pcr_degeneracy.png\n")
