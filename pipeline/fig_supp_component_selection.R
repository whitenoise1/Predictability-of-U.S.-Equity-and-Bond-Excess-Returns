# fig_supp_component_selection.R
# -----------------------------------------------------------------------------
# Canonical generator for the four Supplementary-Methods figures (S1-S4) of the
# component-selection cautionary note. De-scratched companion to the exploratory
# scratch_pcr_explore{2,3,4}/battery scripts (which remain the working record);
# this one re-plots the SAME quantities in the main-paper ggplot aesthetic.
#
#   S1 fig_supp_pcr_degeneracy.png  scree + loading-distance density (explore3)
#   S2 fig_supp_pcr_spectrum.png    SVD filter factors / shrinkage spectrum (explore2)
#   S3 fig_supp_pcr_battery.png     OOS R2 by method + forecast-corr to PLS (battery cache)
#   S4 fig_supp_pcr_spca.png        SPCA & PCR vs PLS scatter + OOS paths (explore4 wf)
#
# Runs from the repo ROOT (props_setup sources R/* and the cached panel from root).
# All cell-level work is the bond_5y h=1 cell, the note's running specimen.
# -----------------------------------------------------------------------------
if (basename(getwd()) == "working_paper") setwd("..")
suppressMessages({library(dplyr); library(tibble); library(tidyr)
                  library(ggplot2); library(patchwork)})
source("working_paper/exploratory/scratch_pcr_props_setup.R")

FIGDIR  <- "working_paper/figures"
SPCACHE <- "working_paper/.scratch_pcr_spca_paths.rds"   # cache the slow explore4 walk-forward

## ---- shared aesthetic (main-paper match) -----------------------------------
raw_col <- "#b2182b"; sco_col <- "#2c6fbb"
meth_col <- c(PCR = "#b2182b", SPCA = "#1b7837", PLS = "#2c6fbb", RIDGE = "#e08214")
g_theme <- theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "grey80", linewidth = 0.3),
        legend.position = "top", legend.title = element_blank(),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 7.5),
        plot.title = element_blank(),
        plot.subtitle = element_text(face = "bold", size = 8.5, hjust = 0.5),
        axis.title = element_text(size = 8))
save_fig <- function(p, file, w = 6.6, h = 2.7)
  ggsave(file.path(FIGDIR, file), p, width = w, height = h, dpi = 200)

cell <- build_cell("bond_5y", 1L)
X <- cell$X_train; y <- cell$y_train; xT <- cell$x_T; p <- ncol(X)
ctr <- colMeans(X); scl <- apply(X, 2L, sd)
Zs  <- sweep(sweep(X, 2L, ctr, "-"), 2L, scl, "/")
sv  <- svd(Zs); d2 <- sv$d^2; cumv <- cumsum(d2) / sum(d2)
k   <- min(which(cumv >= PCR_MIN_VAR)[1L], min(nrow(X) - 1L, p))

## ===== S1 — degeneracy: scree + distance density ============================
instrument_fig <- function(M, min_var = 0.80) {
  pr <- prcomp(M, scale. = TRUE, center = TRUE)
  vr <- pr$sdev^2 / sum(pr$sdev^2); cu <- cumsum(vr)
  npc <- which(cu >= min_var)[1L]; if (is.na(npc)) npc <- length(vr)
  L <- pr$rotation[, seq_len(npc), drop = FALSE]
  D <- as.matrix(dist(L)); list(eig = vr, dist = D[upper.tri(D)])
}
Z <- Zs %*% sv$v[, seq_len(k), drop = FALSE]; colnames(Z) <- paste0("PC", seq_len(k))
rawI <- instrument_fig(X); zI <- instrument_fig(Z)
lev <- c("raw inputs (X)", "component scores (Z)")
scree <- bind_rows(
  tibble(component = seq_along(rawI$eig), share = sort(rawI$eig, TRUE), set = lev[1]),
  tibble(component = seq_along(zI$eig),   share = sort(zI$eig, TRUE),   set = lev[2])) |>
  mutate(set = factor(set, lev))
dens <- bind_rows(tibble(d = rawI$dist, set = lev[1]),
                  tibble(d = zI$dist,   set = lev[2])) |> mutate(set = factor(set, lev))
p_scree <- ggplot(scree, aes(component, share)) +
  geom_line(aes(colour = set), linewidth = 0.5) +
  geom_point(aes(colour = set, shape = set, fill = set), size = 1.5, stroke = 0.5) +
  scale_colour_manual(values = setNames(c(raw_col, sco_col), lev)) +
  scale_fill_manual(values = setNames(c(raw_col, "white"), lev)) +
  scale_shape_manual(values = setNames(c(16, 21), lev)) +
  labs(subtitle = "Internal-PCA scree", x = "component", y = "variance share") + g_theme
p_dens <- ggplot(dens, aes(d, colour = set)) + geom_density(linewidth = 0.7) +
  scale_colour_manual(values = setNames(c(raw_col, sco_col), lev)) +
  labs(subtitle = "Distance screen: dispersed vs inert",
       x = "pairwise loading-distance", y = "density") + g_theme +
  theme(legend.position = "none")
save_fig(p_scree + plot_spacer() + p_dens + plot_layout(widths = c(1, 0.15, 1)),
         "fig_supp_pcr_degeneracy.png", h = 3.4)

## ===== S2 — shrinkage spectrum (SVD filter factors) =========================
edf_ridge <- function(lambda) sum(d2 / (d2 + lambda))
implied_lambda <- function(df_target) {
  if (df_target >= p) return(0)
  uniroot(function(l) edf_ridge(l) - df_target, c(1e-8, 1e10))$root
}
run_pcr <- function(ng) { pr <- pcr_core(X, y, xT, n_groups = ng)
  list(idx = as.integer(sub("^PC", "", pr$selected_pcs)), kk = pr$k) }
pcr_ship <- run_pcr(10L); pcr_fix <- run_pcr(p); k_ret <- pcr_ship$kk
rg <- pit_expanding_cv(X, y, h = 1L, seed = 3L,
        fold_predict = ridge_fold_predict(RIDGE_LAMBDA_GRID), n_candidates = length(RIDGE_LAMBDA_GRID))
lam_cv <- RIDGE_LAMBDA_GRID[rg]; edf_cv <- edf_ridge(lam_cv)
df_fix <- length(pcr_fix$idx)
curves <- bind_rows(
  tibble(component = seq_len(p), f = d2 / (d2 + lam_cv),
         series = sprintf("ridge, CV-tuned (df=%.1f)", edf_cv)),
  tibble(component = seq_len(p), f = d2 / (d2 + implied_lambda(df_fix)),
         series = sprintf("ridge at PCR-implied lambda (df=%d)", df_fix)))
clev <- unique(curves$series)
pts <- bind_rows(
  tibble(component = pcr_fix$idx, f = 1, sel = sprintf("PCR-select, n_groups>=k {%s}", paste(sort(pcr_fix$idx), collapse = ","))),
  tibble(component = pcr_ship$idx, f = 1, sel = sprintf("PCR-select, as shipped {%s}", paste(sort(pcr_ship$idx), collapse = ","))))
slev <- unique(pts$sel)
p_spec <- ggplot() +
  geom_hline(yintercept = 1, colour = "grey70", linetype = 3) +
  geom_vline(xintercept = k_ret, colour = "grey75", linetype = 2, linewidth = 0.3) +
  annotate("text", x = k_ret + 1.5, y = 0.9, label = sprintf("k retained = %d", k_ret),
           hjust = 0, size = 2.5, colour = "grey45") +
  geom_line(data = curves, aes(component, f, colour = series, linetype = series), linewidth = 0.6) +
  geom_point(data = pts, aes(component, f, shape = sel, fill = sel),
             colour = raw_col, size = 2, stroke = 0.6) +
  scale_colour_manual(values = setNames(c(sco_col, "#7fb0d6"), clev)) +
  scale_linetype_manual(values = setNames(c(1, 2), clev)) +
  scale_shape_manual(values = setNames(c(21, 16), slev)) +
  scale_fill_manual(values = setNames(c("white", raw_col), slev)) +
  labs(subtitle = "Where PCR-selection sits on the shrinkage spectrum (bond 5y, h=1)",
       x = "principal component (variance-ordered)", y = expression("shrinkage filter factor  " * f[i])) +
  guides(colour = guide_legend(order = 1, nrow = 2), linetype = guide_legend(order = 1, nrow = 2),
         shape = guide_legend(order = 2, nrow = 2), fill = guide_legend(order = 2, nrow = 2)) +
  g_theme + theme(legend.text = element_text(size = 6.8))
save_fig(p_spec, "fig_supp_pcr_spectrum.png", w = 5.2, h = 5.0)

## ===== S3 — cross-cell battery: OOS R2 + corr-to-PLS ========================
b <- readRDS("working_paper/.scratch_pcr_battery_results.rds")$red
clv <- c("bond_1y", "bond_2y", "bond_5y", "equity"); mlv <- c("PCR", "SPCA", "PLS", "RIDGE")
r2 <- b |> select(cell, PCR = r2_PCR, SPCA = r2_SPCA, PLS = r2_PLS, RIDGE = r2_RIDGE) |>
  pivot_longer(-cell, names_to = "method", values_to = "r2") |>
  mutate(cell = factor(cell, clv), method = factor(method, mlv))
p_r2 <- ggplot(r2, aes(cell, r2, fill = method)) +
  geom_col(width = 0.74, alpha = 0.6, position = position_dodge2(padding = 0.1)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  scale_fill_manual(values = meth_col) +
  labs(subtitle = "Out-of-sample R² by method", x = NULL, y = expression(OOS~R^2)) + g_theme
cr <- b |> select(cell, PCR = corr_PCR_PLS, SPCA = corr_SPCA_PLS) |>
  pivot_longer(-cell, names_to = "method", values_to = "corr") |>
  mutate(cell = factor(cell, clv), method = factor(method, c("PCR", "SPCA")))
p_cr <- ggplot(cr, aes(cell, corr, fill = method)) +
  geom_col(width = 0.66, alpha = 0.6, position = position_dodge2(padding = 0.12)) +
  scale_fill_manual(values = meth_col[c("PCR", "SPCA")]) + ylim(0, 1) +
  labs(subtitle = "Forecast correlation with PLS", x = NULL, y = "corr with PLS") + g_theme
save_fig(p_r2 + plot_spacer() + p_cr + plot_layout(widths = c(1, 0.15, 1)),
         "fig_supp_pcr_battery.png", h = 3.4)

## ===== S4 — SPCA vs PLS: scatter + OOS paths (explore4 walk-forward) ========
if (file.exists(SPCACHE)) { J <- readRDS(SPCACHE) } else {
  MAT <- "bond_5y"; H <- 1L
  yv <- targets_bd[[MAT]][[as.character(H)]]; ref <- panel_xfm$reference_date
  t_start <- t_start_for(S2_DEEP_COLS)
  t_eval  <- seq(t_start, by = "120 months", length.out = 2L)[2L]
  pool <- S2_POOL_COLS; fit_lo0 <- which(ref >= t_start)[1L]
  eval_T <- which(ref >= t_eval & ref <= AS_OF & !is.na(yv)); eval_T <- eval_T[seq(1L, length(eval_T), by = 3L)]
  f_pls <- function(X, y, xT) { mm <- min(3L, ncol(X)); b <- pit_expanding_cv(X, y, h = H, seed = H + 2L,
              fold_predict = pls_fold_predict(mm), n_candidates = mm); if (is.null(b)) NA_real_ else pls_predict(nipals_pls(X, y, b), xT, b) }
  f_ridge <- function(X, y, xT) { b <- pit_expanding_cv(X, y, h = H, seed = H + 2L,
              fold_predict = ridge_fold_predict(RIDGE_LAMBDA_GRID), n_candidates = length(RIDGE_LAMBDA_GRID))
              if (is.null(b)) NA_real_ else ridge_fit(X, y, RIDGE_LAMBDA_GRID[b])$predict(xT) }
  f_pcr <- function(X, y, xT) { r <- pcr_core(X, y, xT); if (is.null(r)) NA_real_ else r$y_hat }
  f_spca <- function(X, y, xT) { r <- spca_cv(X, y, xT, H); if (is.null(r)) NA_real_ else r$y_hat }
  fns <- list(SPCA = f_spca, PCR = f_pcr, PLS = f_pls, RIDGE = f_ridge); rows <- list()
  for (T_idx in eval_T) {
    fit_hi <- T_idx - H; if (fit_hi < fit_lo0) next
    rr <- fit_lo0:fit_hi; xTr <- panel_xfm[T_idx, pool, drop = FALSE]; win <- panel_xfm[rr, pool, drop = FALSE]
    ok <- vapply(pool, function(cc) !is.na(xTr[[cc]]) && !anyNA(win[[cc]]), logical(1)); cols <- pool[ok]
    if (length(cols) < 2L) next
    Xa <- as.matrix(panel_xfm[rr, cols, drop = FALSE]); yt <- yv[rr]
    keep <- stats::complete.cases(Xa) & !is.na(yt); if (sum(keep) < 30L) next
    Xtr <- Xa[keep, , drop = FALSE]; ytr <- yt[keep]; xTm <- as.matrix(panel_xfm[T_idx, cols, drop = FALSE])
    yh <- vapply(fns, function(f) tryCatch(as.numeric(f(Xtr, ytr, xTm)), error = function(e) NA_real_), numeric(1))
    rows[[length(rows) + 1L]] <- data.frame(as_of = ref[T_idx], y_real = yv[T_idx], y_bench = mean(ytr), t(yh))
  }
  D <- do.call(rbind, rows); J <- D[stats::complete.cases(D[, names(fns)]), ]; saveRDS(J, SPCACHE)
}
cmJ <- cor(as.matrix(J[, c("SPCA", "PCR", "PLS")]))
sc <- bind_rows(
  tibble(pls = J$PLS, fc = J$SPCA, m = sprintf("SPCA  (r=%.2f)", cmJ["SPCA", "PLS"])),
  tibble(pls = J$PLS, fc = J$PCR,  m = sprintf("PCR  (r=%.2f)",  cmJ["PCR", "PLS"])))
slv <- unique(sc$m)
p_sc <- ggplot(sc, aes(pls, fc, colour = m, shape = m, fill = m)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2, linewidth = 0.3) +
  geom_point(size = 1.5, stroke = 0.5, alpha = 0.85) +
  scale_colour_manual(values = setNames(c(meth_col["SPCA"], meth_col["PCR"]), slv)) +
  scale_fill_manual(values = setNames(c(meth_col["SPCA"], "white"), slv)) +
  scale_shape_manual(values = setNames(c(16, 21), slv)) +
  coord_equal() +
  labs(subtitle = "Forecasts vs C-PLS forecast", x = "C-PLS forecast", y = "forecast") + g_theme
pa <- J |> select(as_of, SPCA, PCR, PLS) |> pivot_longer(-as_of, names_to = "method", values_to = "fc") |>
  mutate(method = factor(method, c("SPCA", "PCR", "PLS")))
p_path <- ggplot(pa, aes(as_of, fc, colour = method)) + geom_line(linewidth = 0.5) +
  scale_colour_manual(values = meth_col[c("SPCA", "PCR", "PLS")]) +
  labs(subtitle = "Out-of-sample forecast paths", x = "as-of date", y = "forecast") + g_theme
save_fig(p_sc + plot_spacer() + p_path + plot_layout(widths = c(1, 0.15, 1)),
         "fig_supp_pcr_spca.png", h = 3.4)

cat("[fig] wrote S1-S4 supplement figures to", FIGDIR, "\n")
