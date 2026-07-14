# working_paper/ws4_inference_selftest.R
# WS4 Stage-B module self-tests (G3 LRV/CW, G4 FDR, G5 SPA). Synthetic data only;
# no dependency on the residual re-run. RUN FROM REPO ROOT. Fail-loud throughout.
suppressMessages({library(dplyr); library(tibble)})
source("R/predictive_metrics.R"); source("R/predictive_regression.R")
source("R/overlap_robust_inference.R"); source("R/spa_test.R")
set.seed(42)
ok <- function(cond, msg) if (!isTRUE(cond)) stop("FAIL: ", msg) else cat("  ok:", msg, "\n")

cat("[G3] overlap LRV / CW reductions\n")
f <- rnorm(300)
# L=0 NW LRV == gamma_0 = sum((x-xbar)^2)/n
ok(abs(newey_west_lrv_scalar(f, 0L) - sum((f-mean(f))^2)/length(f)) < 1e-12,
   "NW LRV at L=0 equals gamma_0")
# cw_overlap_robust at h=1 -> L_ovl=0; cw_pub uses L=floor(1.5)=1; both finite, n_eff=n
yh <- rnorm(200); yr <- 0.3*yh + rnorm(200); yb <- rep(0, 200)
cro1 <- cw_overlap_robust(yh, yr, yb, h = 1L)
ok(cro1$n_eff == 200L, "h=1 effective n == n")
ok(is.finite(cro1$cw_pub) && is.finite(cro1$cw_overlap), "h=1 CW stats finite")
# cw_pub must equal the canonical clark_west_stat at L=floor(1.5h)
h <- 12L; yh <- rnorm(250); yr <- rnorm(250); yb <- rep(mean(yr), 250)
fT <- clark_west_pointwise(y_hat=yh, y_realized=yr, y_bench=yb)
ref <- clark_west_stat(fT, L = floor(1.5*h))$stat
cro <- cw_overlap_robust(yh, yr, yb, h)
ok(abs(cro$cw_pub - ref) < 1e-10, "cw_pub reproduces clark_west_stat exactly")
ok(cro$n_eff == floor(250/12), "h=12 effective n == floor(n/h)")
ok(cro$p_t_eff > cro$p_normal, "t(n_eff) reference is more conservative than N(0,1) for cw>0")

cat("[G3b] Hodrick 1B vs NW-HAC MZ — overlap inflates SE with h\n")
n <- 400; x <- rnorm(n); r1 <- 0.0*x + rnorm(n)           # one-period returns, no predictability
mk_h <- function(h) { yh <- x; yr <- as.numeric(stats::filter(r1, rep(1,h), sides=1)); # h-sum
                      keep <- !is.na(yr); list(yr=yr[keep], yh=yh[keep], r1=r1[keep]) }
se3  <- hodrick_1b_mz(mk_h(3)$yr,  mk_h(3)$yh,  mk_h(3)$r1,  3L)$se_hodrick
se12 <- hodrick_1b_mz(mk_h(12)$yr, mk_h(12)$yh, mk_h(12)$r1, 12L)$se_hodrick
ok(se12 > se3, "Hodrick 1B SE grows with horizon (more overlap)")
ok(is.finite(nw_hac_mz(mk_h(12)$yr, mk_h(12)$yh, 12L)$se_nwhac), "NW-HAC MZ SE finite at h=12")

cat("[G4] FDR wrappers equal stats::p.adjust\n")
p <- runif(50)
ok(identical(bh_fdr(p), p.adjust(p, "BH")), "bh_fdr == p.adjust BH")
ok(identical(by_fdr(p), p.adjust(p, "BY")), "by_fdr == p.adjust BY")

cat("[G5] SPA behaviour: noise family vs planted signal\n")
nT <- 240; m <- 60
# Pure-noise family: loss differentials mean 0 -> SPA should NOT reject (p large).
Dnoise <- matrix(rnorm(nT*m, mean = 0, sd = 1), nT, m)
sp_noise <- spa_test(Dnoise, q = 1/12, B = 500L, seed = 7L)
ok(sp_noise$p_consistent > 0.10, sprintf("noise family: p_consistent=%.3f not significant", sp_noise$p_consistent))
# Plant one genuinely superior model (positive-mean loss differential).
Dsig <- Dnoise; Dsig[,1] <- rnorm(nT, mean = 0.4, sd = 1)
sp_sig <- spa_test(Dsig, q = 1/12, B = 500L, seed = 7L)
ok(sp_sig$p_consistent < 0.05, sprintf("planted signal: p_consistent=%.3f significant", sp_sig$p_consistent))
ok(sp_sig$p_lower <= sp_sig$p_consistent && sp_sig$p_consistent <= sp_sig$p_upper,
   "p ordering lower <= consistent <= upper")

cat("\nALL SELF-TESTS PASSED\n")
