# scratch_pcr_battery.R
# -----------------------------------------------------------------------------
# Widen the PCR shrinkage-property battery beyond bond_5y to SHORT-END BONDS and
# EQUITY S2, to test whether the three findings generalize:
#   (1) marginal-screened PCR selects LOW-variance (fragile) PCs;
#   (2) supervised PCA (variance-conditioned) closes the gap to the shrinkage class;
#   (3) the near-null holds — no method (incl. PCR/SPCA) beats the benchmark OOS.
# Equity uses the artifact-free SP500TR month-end target with {SP500TR,SHILLER_PRICE}
# purged from its pool (build_cell(target="equity")). All cells h=1, expanding,
# lean quarterly walk-forward (matches scratch_pcr_explore4 protocol).
# -----------------------------------------------------------------------------
source("working_paper/scratch_pcr_props_setup.R")

CELLS <- list(
  list(key = "bond_1y", target = "bond",   maturity = "bond_1y"),
  list(key = "bond_2y", target = "bond",   maturity = "bond_2y"),
  list(key = "bond_5y", target = "bond",   maturity = "bond_5y"),
  list(key = "equity",  target = "equity", maturity = NA_character_))
H <- 1L

# per-as-of predictors (identical rectangles for all methods)
f_pls <- function(X,y,xT){ mm<-min(3L,ncol(X)); b<-pit_expanding_cv(X,y,h=H,seed=H+2L,
            fold_predict=pls_fold_predict(mm),n_candidates=mm); if(is.null(b)) NA_real_ else pls_predict(nipals_pls(X,y,b),xT,b) }
f_ridge <- function(X,y,xT){ b<-pit_expanding_cv(X,y,h=H,seed=H+2L,
            fold_predict=ridge_fold_predict(RIDGE_LAMBDA_GRID),n_candidates=length(RIDGE_LAMBDA_GRID))
            if(is.null(b)) NA_real_ else ridge_fit(X,y,RIDGE_LAMBDA_GRID[b])$predict(xT) }
f_pcr <- function(X,y,xT){ r<-pcr_core(X,y,xT); if(is.null(r)) NA_real_ else r$y_hat }
f_spca <- function(X,y,xT){ r<-spca_cv(X,y,xT,H); if(is.null(r)) NA_real_ else r$y_hat }
FNS <- list(SPCA=f_spca, PCR=f_pcr, PLS=f_pls, RIDGE=f_ridge)

run_wf <- function(target, maturity) {
  if (target == "equity") { yv <- targets_eq[[as.character(H)]]; pool <- S2_POOL_COLS_EQ; anchor <- S2_DEEP_COLS_EQ }
  else { yv <- targets_bd[[maturity]][[as.character(H)]]; pool <- S2_POOL_COLS; anchor <- S2_DEEP_COLS }
  ref <- panel_xfm$reference_date
  t_start <- t_start_for(anchor); t_eval <- seq(t_start, by="120 months", length.out=2L)[2L]
  fit_lo0 <- which(ref >= t_start)[1L]
  eval_T  <- which(ref >= t_eval & ref <= AS_OF & !is.na(yv)); eval_T <- eval_T[seq(1L,length(eval_T),by=6L)]  # semiannual (lean)
  rows <- list()
  for (T_idx in eval_T) {
    fit_hi <- T_idx - H; if (fit_hi < fit_lo0) next
    rr <- fit_lo0:fit_hi
    xT <- panel_xfm[T_idx, pool, drop=FALSE]; win <- panel_xfm[rr, pool, drop=FALSE]
    ok <- vapply(pool, function(cc) !is.na(xT[[cc]]) && !anyNA(win[[cc]]), logical(1))
    cols <- pool[ok]; if (length(cols) < 2L) next
    Xa <- as.matrix(panel_xfm[rr, cols, drop=FALSE]); yt <- yv[rr]
    keep <- stats::complete.cases(Xa) & !is.na(yt); if (sum(keep) < 30L) next
    Xtr <- Xa[keep,,drop=FALSE]; ytr <- yt[keep]; xTm <- as.matrix(panel_xfm[T_idx, cols, drop=FALSE])
    yh <- vapply(FNS, function(f) tryCatch(as.numeric(f(Xtr,ytr,xTm)), error=function(e) NA_real_), numeric(1))
    rows[[length(rows)+1L]] <- data.frame(as_of=ref[T_idx], y_real=yv[T_idx], y_bench=mean(ytr), t(yh))
  }
  D <- do.call(rbind, rows); D[stats::complete.cases(D[, names(FNS)]), ]
}

spec_rows <- list(); red_rows <- list()
for (cl in CELLS) {
  cat(sprintf("\n######## %s ########\n", cl$key))
  cell <- build_cell(cl$maturity, H, target = cl$target)
  X <- cell$X_train; y <- cell$y_train
  d2 <- svd(sweep(sweep(X,2,colMeans(X),"-"),2,apply(X,2,sd),"/"))$d^2; vs <- d2/sum(d2)
  pr <- pcr_core(X, y, cell$x_T)
  sp <- spca_cv(X, y, cell$x_T, H)
  lam <- RIDGE_LAMBDA_GRID[pit_expanding_cv(X,y,h=H,seed=H+2L,
            fold_predict=ridge_fold_predict(RIDGE_LAMBDA_GRID),n_candidates=length(RIDGE_LAMBDA_GRID))]
  redf <- ridge_fit(X,y,lam)$edf
  plsm <- pit_expanding_cv(X,y,h=H,seed=H+2L,fold_predict=pls_fold_predict(min(3L,ncol(X))),n_candidates=min(3L,ncol(X)))
  pcr_idx <- if (is.null(pr)) integer(0) else as.integer(sub("^PC","",pr$selected_pcs))
  spec_rows[[cl$key]] <- data.frame(cell=cl$key, p=ncol(X), k=if(is.null(pr)) NA else pr$k,
    pcr_nsel=if(is.null(pr)) 0L else pr$n_selected,
    pcr_selvar_pct=if(length(pcr_idx)) round(100*mean(vs[pcr_idx]),2) else NA,
    spca_nscr=if(is.null(sp)) 0L else sp$n_screened, spca_m=if(is.null(sp)) NA else sp$m,
    spca_leadvar_pct=if(is.null(sp)) NA else round(100*sp$var_share[1],1),
    ridge_edf=round(redf,1), pls_m=plsm)
  cat(sprintf("  PCR picks {%s} var%%=%.2f | SPCA screens %s -> lead var%%=%.1f | ridge edf=%.1f pls m=%d\n",
      paste(pr$selected_pcs,collapse=","), if(length(pcr_idx)) 100*mean(vs[pcr_idx]) else NA_real_,
      if(is.null(sp)) "0" else sp$n_screened, if(is.null(sp)) NA_real_ else 100*sp$var_share[1], redf, plsm))

  J <- run_wf(cl$target, cl$maturity); meth <- names(FNS)
  cm <- cor(as.matrix(J[, meth]))
  red_rows[[cl$key]] <- data.frame(cell=cl$key, n_oos=nrow(J),
    r2_PCR=r2_oos(J$PCR,J$y_real,J$y_bench), r2_SPCA=r2_oos(J$SPCA,J$y_real,J$y_bench),
    r2_PLS=r2_oos(J$PLS,J$y_real,J$y_bench), r2_RIDGE=r2_oos(J$RIDGE,J$y_real,J$y_bench),
    corr_PCR_PLS=cm["PCR","PLS"], corr_SPCA_PLS=cm["SPCA","PLS"])
  cat(sprintf("  [wf] n_oos=%d | R2: PCR=%.3f SPCA=%.3f PLS=%.3f RIDGE=%.3f | corr-to-PLS: PCR=%.2f SPCA=%.2f\n",
      nrow(J), red_rows[[cl$key]]$r2_PCR, red_rows[[cl$key]]$r2_SPCA, red_rows[[cl$key]]$r2_PLS,
      red_rows[[cl$key]]$r2_RIDGE, cm["PCR","PLS"], cm["SPCA","PLS"]))
}

SPEC <- do.call(rbind, spec_rows); RED <- do.call(rbind, red_rows)
cat("\n\n================ SPECTRUM POSITION (per cell) ================\n"); print(SPEC, row.names=FALSE)
cat("\n================ REDUNDANCY + SKILL (per cell) ===============\n"); print(RED, row.names=FALSE, digits=3)
saveRDS(list(spec=SPEC, red=RED), "working_paper/.scratch_pcr_battery_results.rds")

# ---- summary figure ----------------------------------------------------------
png("working_paper/figures/scratch_pcr_battery.png", 1200, 480, res=130)
op <- par(mfrow=c(1,2), mar=c(5,4,3,1))
r2m <- t(as.matrix(RED[,c("r2_PCR","r2_SPCA","r2_PLS","r2_RIDGE")])); colnames(r2m) <- RED$cell
barplot(r2m, beside=TRUE, col=c("firebrick","seagreen","steelblue","darkorange"),
        ylab="OOS R^2", main="OOS skill by method x cell", las=2, cex.names=0.8)
abline(h=0, col="grey40"); legend("bottomleft", c("PCR","SPCA","PLS","RIDGE"),
       fill=c("firebrick","seagreen","steelblue","darkorange"), bty="n", cex=0.75)
cm2 <- t(as.matrix(RED[,c("corr_PCR_PLS","corr_SPCA_PLS")])); colnames(cm2) <- RED$cell
barplot(cm2, beside=TRUE, col=c("firebrick","seagreen"), ylim=c(0,1),
        ylab="corr with PLS", main="Forecast corr to PLS: PCR vs SPCA", las=2, cex.names=0.8)
legend("topright", c("PCR","SPCA (variance-conditioned)"), fill=c("firebrick","seagreen"), bty="n", cex=0.75)
par(op); invisible(dev.off())
cat("\n[fig] working_paper/figures/scratch_pcr_battery.png\n")
