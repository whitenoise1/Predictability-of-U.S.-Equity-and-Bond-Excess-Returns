#!/usr/bin/env Rscript
# run_revision_audit.R — corpus-wide data-revision audit (Annex exhibit + §3 footnote).
#
# Quantifies how much the panel's GENUINELY-REVISED series get revised between first
# release and latest vintage, across EVERY such series — not one hand-picked example.
# This is the evidence that vintage leakage is a pervasive, distributed phenomenon.
#
# Universe: the 132 VINTAGE-tier, non-forecast series (vintage_source ALFRED or RTDSM;
# family != FORECASTS). SPF survey rounds are excluded — they are successive forecasts,
# not revisions of a realized value.
#
# Metric: revisions are measured on the YEAR-OVER-YEAR log-growth (×100), i.e. the
# transform the model actually consumes (§2.3) and the growth-rate convention of the
# real-time-data literature (Croushore & Stark 2001; Aruoba 2008). Growth cancels the
# base-year rebasing that contaminates chained-$ LEVEL revisions. For each reference
# month t:
#   first-release growth  g0(t) = 100*[log v_adv(t) - log v_adv(t-12)]   (advance vintage of t)
#   latest growth         g1(t) = 100*[log v_lat(t) - log v_lat(t-12)]   (latest vintage)
#   revision              r(t)  = g1(t) - g0(t)
# Per series we report the noise-to-signal ratio NS = sd(r)/sd(g1) (Aruoba), the median
# |r| (pp), the mean signed revision, the share of upward revisions, and counts.
#
# Left-censoring control: a reference month is used only if its advance vintage was
# published within 365 days of the reference month (so g0 is a true first release, not a
# mid-life vintage from before the series' archive began). Series whose growth transform
# is undefined (any non-positive value in range) are EXCLUDED and listed — never silently
# dropped (no-silent-fails).
#
# Run from the repo ROOT. Reads data/raw (in-sandbox; no FRED). Writes slim committed
# artifacts to data/audit/ so the manuscript render needs no data/raw.
#
#   Rscript working_paper/pipeline/run_revision_audit.R

suppressMessages({library(arrow); library(dplyr); library(lubridate)})
source("R/pit_query.R")                       # parse_vintage_to_date, get_factor_path, get_pit_value
meta <- read_parquet(DEFAULT_META_PATH)

universe <- meta %>%
  filter(pit_quality == "VINTAGE", family != "FORECASTS",
         vintage_source %in% c("ALFRED", "RTDSM"))
stopifnot(nrow(universe) == 132L)

CENSOR_DAYS <- 365L
MIN_REF     <- 20L                       # need a usable revision sample to report a series
AS_OF       <- as.Date("2026-04-30")     # "today" — the latest leg is what the panel shows now

one_series <- function(fid, fam, freq) {
  d <- tryCatch(read_parquet(get_factor_path(fid, meta)), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d$pub <- parse_vintage_to_date(d$vintage_id)
  # PIT-valid publication window only: release on/after the reference period and on/before
  # AS_OF. This also defuses a raw-store quirk where some RTDSM monthly vintage_ids have a
  # 2-digit year mis-expanded into the 21st century (e.g. 1964 -> "2064m12"); those parse to
  # a far-future pub and would otherwise hijack the "latest" leg. (Production PIT queries are
  # unaffected: get_pit_value(as_of <= today) already filters future-dated vintages out.)
  d <- d[!is.na(d$pub) & !is.na(d$value) & d$pub >= d$reference_date & d$pub <= AS_OF, ]
  if (!nrow(d)) return(NULL)

  val_key <- d %>% distinct(vintage_id, reference_date, .keep_all = TRUE) %>%
    select(vintage_id, reference_date, value)
  vlookup <- function(vid, rd) {
    val_key$value[match(paste(vid, rd), paste(val_key$vintage_id, val_key$reference_date))]
  }
  adv <- d %>% group_by(reference_date) %>%
    slice_min(pub, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(reference_date, vid_adv = vintage_id, pub_adv = pub, v_adv = value) %>%
    filter(pub_adv <= reference_date + days(CENSOR_DAYS))    # left-censoring control
  if (nrow(adv) < MIN_REF) return(NULL)

  # latest leg via the production accessor itself: per reference month, the most recent
  # vintage with pub <= AS_OF (so the future-dated bad ids are excluded automatically).
  lat <- get_pit_value(fid, AS_OF, meta)
  lat_v <- setNames(lat$value, as.character(lat$reference_date))
  rd_lag <- adv$reference_date %m-% months(12L)
  v_adv_lag <- vlookup(adv$vid_adv, rd_lag)                 # t-12 in the SAME advance vintage
  v_lat_t   <- unname(lat_v[as.character(adv$reference_date)])
  v_lat_lag <- unname(lat_v[as.character(rd_lag)])

  g0 <- 100 * (log(adv$v_adv) - log(v_adv_lag))            # first-release YoY growth
  g1 <- 100 * (log(v_lat_t)  - log(v_lat_lag))            # latest YoY growth
  ok <- is.finite(g0) & is.finite(g1)                      # finite => all inputs > 0
  n_undef <- sum(!ok)
  g0 <- g0[ok]; g1 <- g1[ok]; r <- g1 - g0
  if (length(r) < MIN_REF || sd(g1) == 0) return(NULL)

  data.frame(
    factor_id = fid, family = fam, native_frequency = freq,
    n_ref = length(r), n_vint = length(unique(d$vintage_id)),
    n_undef = n_undef,
    med_abs_rev = median(abs(r)),          # pp of YoY growth
    mean_rev = mean(r),                    # signed (≈0 => unbiased / two-sided)
    ns_ratio = sd(r) / sd(g1),             # Aruoba noise-to-signal
    pct_up = mean(r > 0))
}

rows <- Map(one_series, universe$factor_id, universe$family, universe$native_frequency)
tab  <- bind_rows(rows[!vapply(rows, is.null, logical(1))]) %>% arrange(desc(ns_ratio))

excluded <- setdiff(universe$factor_id, tab$factor_id)     # disclosed, not hidden
cat(sprintf("revision audit: %d/%d series reportable (%d excluded: %s)\n",
            nrow(tab), nrow(universe), length(excluded),
            if (length(excluded)) paste(excluded, collapse = ", ") else "none"))

# --- the named striking example: 2008-Q4 real-GDP q/q annualized (the famous number) ---
gex <- get_pit_value("GDPC1", as.Date("2009-02-15"), meta)
gex2 <- get_pit_value("GDPC1", as.Date("2026-04-30"), meta)
qg <- function(v) { a <- v$value[v$reference_date == as.Date("2008-07-01")]
                    b <- v$value[v$reference_date == as.Date("2008-10-01")]; ((b/a)^4 - 1) * 100 }
g_adv <- qg(gex); g_rev <- qg(gex2)
stopifnot(abs(g_adv - (-3.80)) < 0.05, abs(g_rev - (-8.47)) < 0.05)

dt <- Sys.Date()
write_parquet(tab, file.path("data/audit", sprintf("revision_audit_%s.parquet", dt)))
write_parquet(data.frame(g_adv = g_adv, g_rev = g_rev),
              file.path("data/audit", sprintf("revision_gdp_example_%s.parquet", dt)))

cat(sprintf("\nAGGREGATE (n=%d reportable series):\n", nrow(tab)))
cat(sprintf("  median noise-to-signal ratio : %.2f\n", median(tab$ns_ratio)))
cat(sprintf("  series with NS > 0.5         : %d (%.0f%%)\n",
            sum(tab$ns_ratio > 0.5), 100*mean(tab$ns_ratio > 0.5)))
cat(sprintf("  median |revision| (pp YoY)   : %.2f\n", median(tab$med_abs_rev)))
cat(sprintf("  pooled share upward          : %.0f%% (mean of per-series pct_up)\n",
            100*mean(tab$pct_up)))
cat(sprintf("  GDP 2008-Q4 q/q ann example  : %.1f%% -> %.1f%%\n", g_adv, g_rev))
print(head(tab, 12)); print(tail(tab, 4))
