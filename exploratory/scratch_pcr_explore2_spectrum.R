# scratch_pcr_explore2_spectrum.R
# -----------------------------------------------------------------------------
# EXPLORATION 2 — position of the PCR (Research) method on the shrinkage
# spectrum, via the SVD filter-factor lens (ESL Fig 3.17).
#
# In the SVD basis of the standardized training design Zs = U D V', every linear
# shrinkage estimator shrinks the OLS PC-coordinate coefficient by a filter
# factor f_i, and effective df = sum_i f_i:
#   OLS         f_i = 1            (df = p)
#   ridge       f_i = d_i^2/(d_i^2 + lambda)            (smooth)
#   PCR-select  f_i = 1{component i selected}, else 0   (HARD subset on the
#                                                        variance-ordered basis)
# We place pcr_core (as-shipped n_groups=10, and the n_groups>=k fix from
# Exploration 3) against the CV-tuned ridge/PLS the paper already ships, compute
# PCR's IMPLIED ridge-lambda (the lambda whose ridge df equals PCR's selected
# count), and trace the regularization path over (min_var_explained, min_t_value).
# bond_5y h=1 cell.
# -----------------------------------------------------------------------------
source("working_paper/scratch_pcr_props_setup.R")

cell <- build_cell("bond_5y", 1L)
X <- cell$X_train; y <- cell$y_train; xT <- cell$x_T
p <- ncol(X)
cat(sprintf("\n== bond_5y h=1 | X_train %d x %d ==\n", nrow(X), p))

# ---- SVD of the standardized design (the common basis for all filters) --------
ctr <- colMeans(X); scl <- apply(X, 2L, sd)
Zs  <- sweep(sweep(X, 2L, ctr, "-"), 2L, scl, "/")
d2  <- svd(Zs)$d^2                                   # length-p singular values^2

edf_ridge <- function(lambda) sum(d2 / (d2 + lambda))           # df(lambda)
implied_lambda <- function(df_target) {
  if (df_target >= p) return(0)
  uniroot(function(l) edf_ridge(l) - df_target, c(1e-8, 1e10))$root
}

# ---- the PCR method: as-shipped vs the n_groups>=k fix ------------------------
run_pcr <- function(ng) {
  pr <- pcr_core(X, y, xT, n_groups = ng)
  idx <- as.integer(sub("^PC", "", pr$selected_pcs))             # variance-order indices
  list(pr = pr, idx = idx, df = length(idx))
}
pcr_ship <- run_pcr(10L)              # the OLD default (n_groups=10 < k) -> arbitrary cap active
                                      # (PCR_N_GROUPS now defaults to Inf == the fix; see METHODS_NOTE §4.1)
pcr_fix  <- run_pcr(p)                # n_groups >= k        -> all retained PCs t-gated
k_ret <- pcr_ship$pr$k

# ---- CV-tuned shrinkage the paper already ships (via the harness) -------------
rg <- pit_expanding_cv(X, y, h = 1L, seed = 3L,
                       fold_predict = ridge_fold_predict(RIDGE_LAMBDA_GRID),
                       n_candidates = length(RIDGE_LAMBDA_GRID))
lam_cv <- RIDGE_LAMBDA_GRID[rg]; ridge_cv <- ridge_fit(X, y, lam_cv)
m_max <- min(3L, p)
pg <- pit_expanding_cv(X, y, h = 1L, seed = 3L,
                       fold_predict = pls_fold_predict(m_max), n_candidates = m_max)
pls_m <- pg

# ---- spectrum summary table ---------------------------------------------------
row <- function(method, df, lam, note)
  cat(sprintf("  %-22s  df=%6.2f  lambda=%-11s  %s\n", method, df,
              if (is.na(lam)) "--" else formatC(lam, format = "g", digits = 3), note))
cat(sprintf("\n-- shrinkage spectrum (p=%d, k_retained=%d) --\n", p, k_ret))
row("OLS (no shrink)",        p,              0,                  "full design, df=p")
row("ridge (CV-tuned)",       ridge_cv$edf,  lam_cv,             "smooth filter; paper's C-RIDGE")
row("PLS (CV-tuned)",         pls_m,         NA,                 sprintf("m=%d components (df~m proxy)", pls_m))
row("PCR-select (shipped)",   pcr_ship$df,   implied_lambda(pcr_ship$df),
    sprintf("PCs {%s}; lambda=implied-ridge", paste(pcr_ship$idx, collapse=",")))
row("PCR-select (n_groups>=k)", pcr_fix$df,  implied_lambda(pcr_fix$df),
    sprintf("PCs {%s}; lambda=implied-ridge", paste(pcr_fix$idx, collapse=",")))
cat(sprintf("  [scale] PCR truncation alone (top-%d) -> df=%d; ridge-CV df=%.1f sits at ~%.0f%% of p\n",
            k_ret, k_ret, ridge_cv$edf, 100 * ridge_cv$edf / p))

# ---- regularization path: (min_var_explained, min_t_value) -> (k, n_selected) -
cat("\n-- regularization path (n_groups>=k; pure truncation x t-gate) --\n")
mv_grid <- c(0.50, 0.70, 0.85, 0.95, 0.99)
mt_grid <- c(0.0, 1.0, 1.645, 2.0, 2.5)
cat(sprintf("  %-8s", "min_t\\mv")); for (mv in mv_grid) cat(sprintf(" %6.2f", mv)); cat("   (cell=n_selected; k in [])\n")
for (mt in mt_grid) {
  cat(sprintf("  %-8.3f", mt))
  for (mv in mv_grid) {
    pr <- tryCatch(pcr_core(X, y, xT, min_var_explained = mv, n_groups = p, min_t_value = mt),
                   error = function(e) NULL)
    cat(sprintf(" %6s", if (is.null(pr)) "  0" else sprintf("%d[%d]", pr$n_selected, pr$k)))
  }
  cat("\n")
}

# ---- figure: SVD filter factors (ESL Fig 3.17 for this cell) ------------------
f_ridge_cv  <- d2 / (d2 + lam_cv)
f_ridge_imp <- d2 / (d2 + implied_lambda(pcr_fix$df))
f_pcr_ship  <- as.numeric(seq_len(p) %in% pcr_ship$idx)
f_pcr_fix   <- as.numeric(seq_len(p) %in% pcr_fix$idx)
png("working_paper/figures/scratch_pcr_filterfactors.png", 1100, 520, res = 130)
op <- par(mar = c(4.2, 4.2, 3, 1))
plot(seq_len(p), f_ridge_cv, type = "l", lwd = 2, col = "darkgreen", ylim = c(0, 1.05),
     xlab = "principal component (variance-ordered)", ylab = "shrinkage filter factor  f_i",
     main = "bond_5y h=1 : where PCR-select sits on the shrinkage spectrum")
abline(h = 1, col = "grey70", lty = 3)                       # OLS reference
lines(seq_len(p), f_ridge_imp, lwd = 2, col = "orange", lty = 2)
abline(v = k_ret, col = "grey60", lty = 2)
points(pcr_fix$idx,  f_pcr_fix[pcr_fix$idx],   pch = 1,  col = "steelblue", cex = 1.6, lwd = 2)
points(pcr_ship$idx, rep(1.0, length(pcr_ship$idx)), pch = 16, col = "firebrick", cex = 1.3)
legend("topright", bty = "n", cex = 0.8,
       legend = c("OLS (f=1)", sprintf("ridge CV (lambda=%.2g, df=%.1f)", lam_cv, ridge_cv$edf),
                  sprintf("ridge @ PCR-implied lambda (df=%d)", pcr_fix$df),
                  sprintf("PCR-select n_groups>=k {%s}", paste(pcr_fix$idx, collapse=",")),
                  sprintf("PCR-select shipped {%s}", paste(pcr_ship$idx, collapse=",")),
                  sprintf("k retained = %d", k_ret)),
       col = c("grey70","darkgreen","orange","steelblue","firebrick","grey60"),
       lty = c(3,1,2,NA,NA,2), pch = c(NA,NA,NA,1,16,NA), lwd = 2)
par(op); invisible(dev.off())
cat("\n[fig] working_paper/figures/scratch_pcr_filterfactors.png\n")
