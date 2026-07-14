# scratch_pcr_explore1_redundancy.R
# -----------------------------------------------------------------------------
# EXPLORATION 1 — is the PCR (Research) method redundant with the shrinkage
# class the paper already ships? (audit METHOD_AUDIT_GS_PCA.md s7: A vs B.)
#
# Walk-forward OOS over the bond_5y h=1 cell, identical per-as-of rectangles, for
# C-PCR vs C-RIDGE / C-PLS / C-ENET (+ C-OLS deep-cohort reference). Then:
#   (a) cross-method OOS forecast correlation  (PCR vs C-PLS is the key cell);
#   (b) OOS R^2 / signal-Sharpe side by side;
#   (c) forecast-DIRECTION cosine at the final as-of (beta vectors in the shared
#       standardized predictor basis) + PCR-retained vs PLS principal angles.
# -----------------------------------------------------------------------------
source("working_paper/scratch_pcr_props_setup.R")

MAT <- "bond_5y"; H <- 1L
y   <- targets_bd[[MAT]][[as.character(H)]]
ref <- panel_xfm$reference_date
t_start <- t_start_for(S2_DEEP_COLS)
t_eval  <- seq(t_start, by = "120 months", length.out = 2L)[2L]   # chapter 2 t_eval_from
cat(sprintf("\n== %s h=%d | t_start=%s t_eval=%s t_end=%s ==\n",
            MAT, H, t_start, t_eval, AS_OF))

panel_shrink <- panel_xfm |> dplyr::select(reference_date, dplyr::all_of(S2_POOL_COLS))
panel_deep   <- panel_xfm |> dplyr::select(reference_date, dplyr::all_of(S2_DEEP_COLS))

wf <- function(cfg, pool_arg) {
  args <- list(y = y, panel = panel_shrink, horizon = H, window_kind = "expanding",
               t_start = t_start, t_eval_start = t_eval, t_end = AS_OF,
               model_config = cfg, target = "bond")
  out <- do.call(walk_forward_fs, c(args, pool_arg))
  out$forecasts |> dplyr::transmute(as_of_date, y_realized, y_bench, !!cfg := y_hat)
}
cat("[wf] running PCR / RIDGE / PLS / ENET ...\n")
f_pcr   <- wf("C-PCR",   list(pcr_pool    = S2_POOL_COLS))
f_ridge <- wf("C-RIDGE", list(shrink_pool = S2_POOL_COLS))
f_pls   <- wf("C-PLS",   list(shrink_pool = S2_POOL_COLS))
f_enet  <- wf("C-ENET",  list(enet_pool   = S2_POOL_COLS))
# C-OLS reference on the deep cohort (needs n>p complete cases; fixed features).
f_ols <- do.call(walk_forward_fs, list(
  y = y, panel = panel_deep, horizon = H, window_kind = "expanding",
  t_start = t_start, t_eval_start = t_eval, t_end = AS_OF,
  model_config = "C-OLS", target = "bond", features = S2_DEEP_COLS))$forecasts |>
  dplyr::transmute(as_of_date, `C-OLS` = y_hat)

# ---- align on common as-of dates (all methods non-NA) -------------------------
key <- f_pcr |> dplyr::select(as_of_date, y_realized, y_bench)
J <- key |>
  dplyr::inner_join(f_pcr   |> dplyr::select(as_of_date, `C-PCR`),   by = "as_of_date") |>
  dplyr::inner_join(f_ridge |> dplyr::select(as_of_date, `C-RIDGE`), by = "as_of_date") |>
  dplyr::inner_join(f_pls   |> dplyr::select(as_of_date, `C-PLS`),   by = "as_of_date") |>
  dplyr::inner_join(f_enet  |> dplyr::select(as_of_date, `C-ENET`),  by = "as_of_date") |>
  dplyr::inner_join(f_ols,  by = "as_of_date") |>
  na.omit()
methods <- c("C-PCR","C-RIDGE","C-PLS","C-ENET","C-OLS")
cat(sprintf("[align] %d common OOS as-of dates\n", nrow(J)))

# ---- (a) forecast correlation matrix -----------------------------------------
M  <- as.matrix(J[, methods])
cm <- cor(M)
cat("\n-- (a) OOS forecast correlation --\n      ")
cat(sprintf("%8s", methods), "\n")
for (i in methods) { cat(sprintf("  %-7s", i)); cat(sprintf("%8.3f", cm[i, ])); cat("\n") }
cat(sprintf("\n  >>> corr(C-PCR, C-PLS) = %.3f   corr(C-PCR, C-RIDGE) = %.3f\n",
            cm["C-PCR","C-PLS"], cm["C-PCR","C-RIDGE"]))

# ---- (b) OOS R^2 / Sharpe -----------------------------------------------------
cat("\n-- (b) OOS skill (common support) --\n")
cat(sprintf("  %-8s %9s %9s\n", "method", "R2_oos", "sharpe"))
for (mth in methods) {
  r2 <- r2_oos(J[[mth]], J$y_realized, J$y_bench)
  sh <- signal_sharpe(J[[mth]], J$y_realized, h = H)
  cat(sprintf("  %-8s %9.4f %9.3f\n", mth, r2, sh))
}

# ---- (c) forecast-direction cosine at the final as-of ------------------------
cell <- build_cell(MAT, H)
Xc <- cell$X_train; yc <- cell$y_train
pr <- pcr_core(Xc, yc, cell$x_T)                           # fixed defaults (no pre-cap)
lam_cv <- RIDGE_LAMBDA_GRID[pit_expanding_cv(Xc, yc, h = H, seed = H + 2L,
              fold_predict = ridge_fold_predict(RIDGE_LAMBDA_GRID),
              n_candidates = length(RIDGE_LAMBDA_GRID))]
rd <- ridge_fit(Xc, yc, lam_cv)
m_max <- min(3L, ncol(Xc))
pm <- pit_expanding_cv(Xc, yc, h = H, seed = H + 2L,
          fold_predict = pls_fold_predict(m_max), n_candidates = m_max)
plf <- nipals_pls(Xc, yc, pm)
B_pls <- as.numeric(plf$W %*% solve(crossprod(plf$P, plf$W), plf$Q))   # std-X coef
names(B_pls) <- colnames(Xc)

cosang <- function(a, b) sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
dirs <- list(PCR = pr$beta_x, RIDGE = rd$beta_std, PLS = B_pls)
cat(sprintf("\n-- (c) forecast-direction cosine at final as-of (%s) --\n", cell$T_date))
nm <- names(dirs)
for (i in seq_along(nm)) for (j in seq_along(nm)) if (j > i)
  cat(sprintf("  cos(%s, %s) = %+.3f\n", nm[i], nm[j], cosang(dirs[[nm[i]]], dirs[[nm[j]]])))

# PCR-retained subspace vs PLS-component subspace: principal angles.
Vk <- pr$loadings                                          # p x k retained PCA basis
Wp <- qr.Q(qr(plf$W))                                      # p x m PLS direction basis
sing <- svd(crossprod(Vk, Wp))$d
cat(sprintf("  principal angles PCR-retained(k=%d) vs PLS(m=%d): %s deg\n",
            ncol(Vk), ncol(Wp), paste(sprintf("%.1f", acos(pmin(sing,1)) * 180/pi), collapse=", ")))

cat(sprintf("\n  selected: PCR PCs {%s} (df=%d) | PLS m=%d | ridge df=%.1f\n",
            paste(pr$selected_pcs, collapse=","), pr$n_selected, pm, rd$edf))

# ---- figure: forecast scatter PCR vs PLS / RIDGE + time series ----------------
png("working_paper/figures/scratch_pcr_redundancy.png", 1150, 460, res = 130)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
lim <- range(M[, c("C-PCR","C-PLS","C-RIDGE")])
plot(J$`C-PLS`, J$`C-PCR`, pch = 16, col = "#33667799", xlim = lim, ylim = lim,
     xlab = "C-PLS forecast", ylab = "C-PCR forecast",
     main = sprintf("PCR vs PLS forecasts (r=%.2f)", cm["C-PCR","C-PLS"]))
abline(0, 1, col = "grey60", lty = 2)
points(J$`C-RIDGE`, J$`C-PCR`, pch = 1, col = "#cc663399")
legend("topleft", c("vs C-PLS","vs C-RIDGE","45 deg"), pch = c(16,1,NA), lty = c(NA,NA,2),
       col = c("#336677","#cc6633","grey60"), bty = "n", cex = 0.8)
ord <- order(J$as_of_date)
matplot(J$as_of_date[ord], M[ord, c("C-PCR","C-PLS","C-RIDGE")], type = "l", lty = 1, lwd = 1.5,
        col = c("firebrick","steelblue","darkgreen"), xlab = "as-of", ylab = "forecast",
        main = "OOS forecast paths")
legend("topright", c("C-PCR","C-PLS","C-RIDGE"), lwd = 1.5,
       col = c("firebrick","steelblue","darkgreen"), bty = "n", cex = 0.8)
par(op); invisible(dev.off())
cat("\n[fig] working_paper/figures/scratch_pcr_redundancy.png\n")
