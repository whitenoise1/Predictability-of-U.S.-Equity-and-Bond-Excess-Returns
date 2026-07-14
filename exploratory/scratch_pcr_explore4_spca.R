# scratch_pcr_explore4_spca.R
# -----------------------------------------------------------------------------
# Option 2 prototype — SUPERVISED PCA (Bair, Hastie, Paul & Tibshirani 2006) as
# the variance-CONDITIONED repair of the GS+PCA selection problem.
#
# The defect (Explorations 2-3): the t-gate selects PCs on in-sample correlation
# alone (t_k ~ corr(y,Z_k)*sqrt(n) -- variance cancels), so it grabs fragile
# LOW-variance directions {PC12,13,16}. Supervised PCA injects the missing
# variance condition by SCREENING RAW FEATURES first, then taking the LEADING
# (high-variance) PCs of the screened block -- components that are predictive AND
# stable by construction:
#   1. screen: keep raw features with univariate |t| vs y >= theta
#   2. PCA the screened block; take leading m components (PIT: train-only basis)
#   3. lm(y ~ leading scores); project X_pred through the training basis
#   (theta, m) chosen by the SAME pit_expanding_cv engine ridge/PLS use.
#
# Question: does variance-conditioned select-then-fit beat the marginal-screened
# PCR and close the gap to PLS? bond_5y h=1 cell.
# -----------------------------------------------------------------------------
source("working_paper/scratch_pcr_props_setup.R")

# ---- spca_core ---------------------------------------------------------------
# X_train/X_pred: raw (z-scored-panel) matrices WITHOUT intercept. Returns NULL
# when fewer than 2 features clear the screen (caller falls back to benchmark --
# surfaced, never silent).
spca_core <- function(X_train, y_train, X_pred, theta, m) {
  n <- nrow(X_train)
  r  <- as.numeric(cor(X_train, y_train))                  # univariate corr (Bair score)
  tt <- abs(r) * sqrt((n - 2) / pmax(1 - r^2, 1e-12))      # corr t-stat
  tt[!is.finite(tt)] <- 0                                  # constant col -> no info
  keep <- which(tt >= theta)
  if (length(keep) < 2L) return(NULL)
  ctr <- colMeans(X_train)[keep]; scl <- apply(X_train[, keep, drop = FALSE], 2L, sd)
  scl[!is.finite(scl) | scl == 0] <- 1
  Zs <- sweep(sweep(X_train[, keep, drop = FALSE], 2L, ctr, "-"), 2L, scl, "/")
  sv <- svd(Zs); mm <- min(as.integer(m), ncol(Zs), n - 1L)
  V  <- sv$v[, seq_len(mm), drop = FALSE]
  S  <- Zs %*% V
  fit <- lm.fit(cbind(1, S), y_train); b <- fit$coefficients
  Zp  <- sweep(sweep(X_pred[, keep, drop = FALSE], 2L, ctr, "-"), 2L, scl, "/")
  y_hat <- as.numeric(cbind(1, Zp %*% V) %*% b)
  beta_x <- numeric(ncol(X_train)); beta_x[keep] <- as.numeric(V %*% b[-1L])
  names(beta_x) <- colnames(X_train)
  list(y_hat = y_hat, n_screened = length(keep), m = mm, df = mm,
       keep = keep, V = V, var_share = sv$d[seq_len(mm)]^2 / sum(sv$d^2),
       beta_x = beta_x)
}

# Lean prototype: fix theta at Bair's canonical screen, CV only the component
# count m (3 candidates). Full (theta x m) CV is a refinement, not the question.
SPCA_THETA <- 1.645; SPCA_M <- 1:3
SPCA_CANDS <- lapply(SPCA_M, function(mm) list(theta = SPCA_THETA, m = as.integer(mm)))
spca_fold_predict <- function(cands) function(Xtr, ytr, xval)
  vapply(cands, function(cd) {
    r <- spca_core(Xtr, ytr, xval, cd$theta, cd$m)
    if (is.null(r)) mean(ytr) else r$y_hat                 # benchmark fallback (finite)
  }, numeric(1))
spca_cv <- function(X, y, xT, h = 1L) {
  cands <- SPCA_CANDS
  best <- pit_expanding_cv(X, y, h = h, seed = h + 2L,
            fold_predict = spca_fold_predict(cands), n_candidates = length(cands))
  if (is.null(best)) return(NULL)
  cd <- cands[[best]]; fit <- spca_core(X, y, xT, cd$theta, cd$m)
  if (!is.null(fit)) { fit$theta <- cd$theta }
  fit
}

# =============================================================================
# A. Final-cell diagnostic: where does variance-conditioned selection land?
# =============================================================================
cell <- build_cell("bond_5y", 1L); X <- cell$X_train; y <- cell$y_train
cat(sprintf("\n== bond_5y h=1 | X_train %d x %d ==\n", nrow(X), ncol(X)))
sp <- spca_cv(X, y, cell$x_T)
pr <- pcr_core(X, y, cell$x_T)                             # marginal-screened PCR (fixed defaults)

# full-design variance rank of the directions each method actually uses
d2_full <- svd(sweep(sweep(X,2,colMeans(X),"-"),2,apply(X,2,sd),"/"))$d^2
var_share_full <- d2_full / sum(d2_full)
pcr_idx <- as.integer(sub("^PC","",pr$selected_pcs))
cat("\n-- SUPERVISED PCA (Bair) --\n")
cat(sprintf("  CV-chosen theta=%.3f, m=%d  (df=%d)\n", sp$theta, sp$m, sp$df))
cat(sprintf("  screened %d of %d raw features clear |t|>=%.2f\n", sp$n_screened, ncol(X), sp$theta))
top_scr <- cell$cols[sp$keep][order(-abs(cor(X[,sp$keep,drop=FALSE], y)))][1:min(6,sp$n_screened)]
cat(sprintf("  top screened features: %s\n", paste(top_scr, collapse=", ")))
cat(sprintf("  leading-component variance share WITHIN screened block: %s\n",
            paste(sprintf("%.1f%%", 100*sp$var_share), collapse=", ")))
cat("\n-- contrast: marginal-screened PCR --\n")
cat(sprintf("  selected PCs {%s} -> full-design variance rank %s (of %d)\n",
            paste(pr$selected_pcs,collapse=","), paste(pcr_idx,collapse=","), ncol(X)))
cat(sprintf("  those PCs' full-design variance share: %s  <= LOW-variance\n",
            paste(sprintf("%.1f%%", 100*var_share_full[pcr_idx]), collapse=", ")))
cat(sprintf("  SPCA leading dir lives in the HIGH-variance region: screened-block PC1 var=%.0f%% vs PCR picks ~%.1f%%\n",
            100*sp$var_share[1], 100*mean(var_share_full[pcr_idx])))

# =============================================================================
# B. Walk-forward: SPCA vs PCR vs PLS vs RIDGE on IDENTICAL per-as-of rectangles
# =============================================================================
MAT <- "bond_5y"; H <- 1L
yv  <- targets_bd[[MAT]][[as.character(H)]]; ref <- panel_xfm$reference_date
t_start <- t_start_for(S2_DEEP_COLS)
t_eval  <- seq(t_start, by = "120 months", length.out = 2L)[2L]
pool <- S2_POOL_COLS
fit_lo0 <- which(ref >= t_start)[1L]
eval_T  <- which(ref >= t_eval & ref <= AS_OF & !is.na(yv))
eval_T  <- eval_T[seq(1L, length(eval_T), by = 3L)]         # quarterly subsample (lean)
cat(sprintf("\n[wf] %d candidate as-ofs, quarterly (%s .. %s) ...\n", length(eval_T), ref[min(eval_T)], ref[max(eval_T)]))

f_pls <- function(X,y,xT){ mm<-min(3L,ncol(X)); b<-pit_expanding_cv(X,y,h=H,seed=H+2L,
            fold_predict=pls_fold_predict(mm),n_candidates=mm); if(is.null(b)) NA_real_
            else pls_predict(nipals_pls(X,y,b),xT,b) }
f_ridge <- function(X,y,xT){ b<-pit_expanding_cv(X,y,h=H,seed=H+2L,
            fold_predict=ridge_fold_predict(RIDGE_LAMBDA_GRID),n_candidates=length(RIDGE_LAMBDA_GRID))
            if(is.null(b)) NA_real_ else ridge_fit(X,y,RIDGE_LAMBDA_GRID[b])$predict(xT) }
f_pcr <- function(X,y,xT){ r<-pcr_core(X,y,xT); if(is.null(r)) NA_real_ else r$y_hat }
f_spca <- function(X,y,xT){ r<-spca_cv(X,y,xT,H); if(is.null(r)) NA_real_ else r$y_hat }
fns <- list(SPCA=f_spca, PCR=f_pcr, PLS=f_pls, RIDGE=f_ridge)

rows <- list()
for (T_idx in eval_T) {
  fit_hi <- T_idx - H; if (fit_hi < fit_lo0) next
  rr <- fit_lo0:fit_hi
  xT <- panel_xfm[T_idx, pool, drop=FALSE]; win <- panel_xfm[rr, pool, drop=FALSE]
  ok <- vapply(pool, function(cc) !is.na(xT[[cc]]) && !anyNA(win[[cc]]), logical(1))
  cols <- pool[ok]; if (length(cols) < 2L) next
  Xa <- as.matrix(panel_xfm[rr, cols, drop=FALSE]); yt <- yv[rr]
  keep <- stats::complete.cases(Xa) & !is.na(yt)
  if (sum(keep) < 30L) next
  Xtr <- Xa[keep,,drop=FALSE]; ytr <- yt[keep]; xTm <- as.matrix(panel_xfm[T_idx, cols, drop=FALSE])
  yh <- vapply(fns, function(f) tryCatch(as.numeric(f(Xtr,ytr,xTm)), error=function(e) NA_real_), numeric(1))
  rows[[length(rows)+1L]] <- data.frame(as_of=ref[T_idx], y_real=yv[T_idx], y_bench=mean(ytr),
                                        t(yh))
  if (length(rows) %% 10L == 0L) message(sprintf("    ... %d as-ofs done", length(rows)))
}
D <- do.call(rbind, rows); meth <- names(fns)
J <- D[stats::complete.cases(D[, meth]), ]
cat(sprintf("[align] %d common OOS as-ofs\n", nrow(J)))

cm <- cor(as.matrix(J[, meth]))
cat("\n-- (a) OOS forecast correlation --\n        "); cat(sprintf("%8s", meth),"\n")
for (i in meth){ cat(sprintf("  %-6s", i)); cat(sprintf("%8.3f", cm[i,])); cat("\n") }
cat(sprintf("\n  >>> corr(SPCA,PLS)=%.3f  corr(SPCA,PCR)=%.3f  corr(PCR,PLS)=%.3f\n",
            cm["SPCA","PLS"], cm["SPCA","PCR"], cm["PCR","PLS"]))
cat("\n-- (b) OOS skill --\n"); cat(sprintf("  %-7s %9s %9s\n","method","R2_oos","sharpe"))
for (mth in meth)
  cat(sprintf("  %-7s %9.4f %9.3f\n", mth,
              r2_oos(J[[mth]], J$y_real, J$y_bench), signal_sharpe(J[[mth]], J$y_real, h=H)))

# ---- figure ------------------------------------------------------------------
png("working_paper/figures/scratch_spca_compare.png", 1150, 460, res=130)
op <- par(mfrow=c(1,2), mar=c(4,4,3,1))
lim <- range(as.matrix(J[,c("SPCA","PLS","PCR")]))
plot(J$PLS, J$SPCA, pch=16, col="#22884499", xlim=lim, ylim=lim,
     xlab="C-PLS forecast", ylab="forecast", main=sprintf("SPCA & PCR vs PLS"))
points(J$PLS, J$PCR, pch=1, col="#bb333399"); abline(0,1,col="grey60",lty=2)
legend("topleft", c(sprintf("SPCA (r=%.2f)",cm["SPCA","PLS"]),
                    sprintf("PCR (r=%.2f)",cm["PCR","PLS"]),"45 deg"),
       pch=c(16,1,NA), lty=c(NA,NA,2), col=c("#228844","#bb3333","grey60"), bty="n", cex=0.8)
ord <- order(J$as_of)
matplot(J$as_of[ord], as.matrix(J[ord,c("SPCA","PCR","PLS")]), type="l", lty=1, lwd=1.4,
        col=c("seagreen","firebrick","steelblue"), xlab="as-of", ylab="forecast",
        main="OOS forecast paths")
legend("topright", c("SPCA","PCR","PLS"), lwd=1.4, col=c("seagreen","firebrick","steelblue"), bty="n", cex=0.8)
par(op); invisible(dev.off())
cat("\n[fig] working_paper/figures/scratch_spca_compare.png\n")
