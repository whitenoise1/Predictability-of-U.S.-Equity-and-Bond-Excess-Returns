# working_paper/verify_equity_target_fix.R
#
# Side-by-side: equity standalone predictability under the OLD target
# (SHILLER_PRICE = monthly AVERAGE, price-only) vs the FIXED target
# (SP500TR = month-end close, total return; frozen Yahoo ^SP500TR snapshot).
# Reuses construct_equity_excess_return() via its sp_col arg (no targets.R edit).
# RUN FROM REPO ROOT. Read-only; writes nothing to data/.

suppressMessages({library(arrow);library(dplyr);library(tibble);library(readr)})
source("R/pit_query.R");source("R/panel.R");source("R/transforms.R")
source("R/targets.R");source("R/predictive_regression.R")
source("R/predictive_metrics.R");source("R/walk_forward.R")

AS_OF <- as.Date("2026-04-30"); W_Z <- 60L; H <- 1L
LEVEL_FACTORS <- c("AAII_BULLBEAR","SKEW","VIXCLS","ANFCI","NFCICREDIT","NFCI",
                   "NFCIRISK","NFCILEVERAGE","GACDFSA066MSFRBPHI",
                   "CORESTICKM159SFRBATL","AAA","ADS_INDEX")
RATIO_FACTORS <- "CPIAUCSL"; FACTORS <- c(LEVEL_FACTORS, RATIO_FACTORS)

tr <- read_csv("sources/equity/sp500tr_monthly_yahoo_2026-05-31.csv",
               show_col_types = FALSE) |>
  mutate(reference_date = as.Date(reference_date))

panel <- build_panel(AS_OF, unique(c(FACTORS,"SHILLER_PRICE","DGS3MO")),
                     "M","forward") |>
  left_join(tr, by = "reference_date") |>
  rename(SP500TR = tr_close)

xfm <- panel |>
  apply_transform(yoy_log_diff,   factor_ids = RATIO_FACTORS) |>
  apply_transform(rolling_zscore, factor_ids = RATIO_FACTORS, window = W_Z, min_obs = W_Z) |>
  apply_transform(rolling_zscore, factor_ids = LEVEL_FACTORS, window = W_Z, min_obs = W_Z)

y_old <- construct_equity_excess_return(panel, h = H, sp_col = "SHILLER_PRICE") # avg, price-only
y_new <- construct_equity_excess_return(panel, h = H, sp_col = "SP500TR")        # month-end, TR

ar1 <- function(v){v<-v[is.finite(v)];round(acf(v,1,plot=FALSE)$acf[2],3)}
cat(sprintf("AR(1) equity return  OLD(Shiller-avg)=%.3f   NEW(SP500TR month-end)=%.3f\n",
            ar1(y_old), ar1(y_new)))
cat("(true month-end returns ~0; the OLD +0.23 was the averaging artifact)\n\n")

scan_one <- function(fac, y){
  ref <- xfm$reference_date
  i0 <- which(!is.na(xfm[[fac]]))[1L]; if(is.na(i0)) return(NULL)
  ts <- ref[i0]; te <- seq(ts, by="120 months", length.out=2L)[2L]
  if (is.na(te) || te >= AS_OF) return(NULL)
  pc <- xfm |> select(reference_date, all_of(fac))
  wf <- tryCatch(walk_forward_predict(y=y, panel=pc, horizon=H, window_kind="expanding",
                 t_start=ts, t_eval_start=te, t_end=AS_OF), error=function(e) NULL)
  if (is.null(wf) || nrow(wf) < 24L) return(NULL)
  f_T <- clark_west_pointwise(wf$y_hat, wf$y_realized, wf$y_bench)
  tibble(r2=r2_oos(wf$y_hat,wf$y_realized,wf$y_bench),
         cw=clark_west_stat(f_T,L=floor(1.5*H))$stat,
         sh=signal_sharpe(wf$y_hat,wf$y_realized,h=H),
         n=nrow(wf))
}

cmp <- bind_rows(lapply(FACTORS, function(f){
  o <- scan_one(f, y_old); n <- scan_one(f, y_new)
  if (is.null(o)||is.null(n)) return(NULL)
  tibble(factor=f,
         r2_old=o$r2, r2_new=n$r2,
         cw_old=o$cw, cw_new=n$cw,
         sh_old=o$sh, sh_new=n$sh,
         acc_old=accept_decision(o$r2,o$cw,o$sh),
         acc_new=accept_decision(n$r2,n$cw,n$sh))
})) |> arrange(desc(r2_old))

cat("================ OLD (Shiller avg) vs NEW (SP500TR month-end) — equity h=1 ================\n")
print(cmp |> mutate(across(where(is.numeric), \(x) round(x,3))), n=Inf, width=Inf)
cat(sprintf("\nStandalone accepters: OLD=%d   NEW=%d   (of %d factors)\n",
            sum(cmp$acc_old,na.rm=TRUE), sum(cmp$acc_new,na.rm=TRUE), nrow(cmp)))

# skip-month IS R^2 on the NEW target — confirm the overlap artifact is gone
lead <- function(x,k) c(x[(k+1):length(x)],rep(NA,k))
isr2 <- function(xz,yv){ok<-is.finite(xz)&is.finite(yv); if(sum(ok)<50) NA else round(cor(xz[ok],yv[ok])^2,3)}
cat("\n================ skip-month IS R2 on the NEW month-end TR target ================\n")
for(f in c("AAII_BULLBEAR","VIXCLS")){
  xz <- xfm[[f]]
  cat(sprintf("%-14s r_{T,T+1}=%.3f  r_{T+1,T+2}=%.3f  r_{T+2,T+3}=%.3f\n",
              f, isr2(xz,y_new), isr2(xz,lead(y_new,1)), isr2(xz,lead(y_new,2))))
}
cat("(if flat across the gap -> genuine 1-mo-ahead signal; no averaging overlap)\n[done]\n")
