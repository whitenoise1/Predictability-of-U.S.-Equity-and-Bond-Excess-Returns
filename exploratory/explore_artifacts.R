#!/usr/bin/env Rscript
# explore_artifacts.R — interactive loader for every generated artifact.
#
# Purpose: pull every committed artifact into a live R session so you can poke at
# any object directly. Loads all data/audit/*.parquet into a named list `A`, binds
# the canonical (latest-dated) ones to short names, exposes a raw-factor loader,
# and prints an index. Nothing here mutates anything on disk.
#
# Usage (from an R session, working dir anywhere inside the repo):
#   source("working_paper/explore_artifacts.R")
#   ls_artifacts()              # index of everything loaded
#   peek(p10)                   # dim + head of an object
#   cols(eqfc)                  # column names
#   get_raw("UNRATE")           # load one raw PIT factor parquet
#   src_R()                     # (optional) source R/ modules so functions are callable
#
# Short names bound: meta, p10, fdr, spa, eqi, s2man, m3, usedv, eqfc

suppressMessages({library(arrow); library(dplyr); library(tibble)})

# --- locate project root (the dir that contains data/audit) -------------------
.find_root <- function() {
  d <- normalizePath(getwd())
  for (i in 1:6) {
    if (dir.exists(file.path(d, "data", "audit"))) return(d)
    d <- dirname(d)
  }
  stop("explore_artifacts: no data/audit found at or above ", getwd())
}
ROOT <- .find_root()

# --- load every audit parquet into the named list A ---------------------------
.audit_dir <- file.path(ROOT, "data", "audit")
.files <- list.files(.audit_dir, pattern = "\\.parquet$", full.names = TRUE)
A <- setNames(lapply(.files, function(f)
                tryCatch(read_parquet(f), error = function(e) {
                  message("  ! skipped (read error): ", basename(f)); NULL })),
              sub("\\.parquet$", "", basename(.files)))
A <- A[!vapply(A, is.null, logical(1))]

# --- canonical (latest-dated) objects, bound to short names -------------------
.latest <- function(stem) {
  hits <- grep(paste0("^", stem), names(A), value = TRUE)
  if (!length(hits)) return(NULL)
  A[[sort(hits, decreasing = TRUE)[1]]]
}
meta  <- read_parquet(file.path(ROOT, "data", "metadata", "factor_metadata.parquet"))
p10   <- .latest("p10_rescored")              # all trials: r2_oos / cw / sharpe / accept_p10 ...
fdr   <- .latest("stage4_fdr")                # universe-wide BH/BY-FDR (q_*_self/pool/univ)
spa   <- .latest("stage4_spa")                # Hansen SPA per (horizon x sample)
eqi   <- .latest("stage4_equity_inference")   # equity-S2 overlap-robust CW + Hodrick + FDR
s2man <- .latest("s2_pool_manifest")          # U6 114-feature S2 pool composition
m3    <- .latest("m3_attribution")            # native shrinkage attribution (ridge/PLS)
usedv <- .latest("annex_used_variables")      # audited 119-feature table
eqfc  <- .latest("equity_s2_forecasts")       # pointwise y_hat / y_realized / y_bench / as_of_date

# --- raw factor loader (single-series PIT parquet) ----------------------------
get_raw <- function(factor_id) {
  row <- meta[meta$factor_id == factor_id, , drop = FALSE]
  if (!nrow(row)) stop("unknown factor_id: ", factor_id)
  read_parquet(file.path(ROOT, "data", "raw", row$family[1L],
                         paste0(tolower(factor_id), ".parquet")))
}

# --- helpers ------------------------------------------------------------------
ls_artifacts <- function() {
  tibble(object = names(A),
         rows   = vapply(A, nrow, 0L),
         cols   = vapply(A, ncol, 0L)) |>
    arrange(object) |> as.data.frame()
}
peek <- function(x, n = 6L) {
  cat("dim:", paste(dim(x), collapse = " x "), "\n")
  print(utils::head(as.data.frame(x), n)); invisible(x)
}
cols <- function(x) names(x)
src_R <- function() {                          # opt-in: side-effecting
  fs <- list.files(file.path(ROOT, "R"), pattern = "\\.R$", full.names = TRUE)
  ok <- vapply(fs, function(f) isTRUE(tryCatch({source(f); TRUE},
                error = function(e) {message("  ! ", basename(f), ": ", conditionMessage(e)); FALSE})),
               logical(1))
  cat(sprintf("sourced %d/%d R/ modules\n", sum(ok), length(fs)))
}

cat(sprintf("explore_artifacts: root=%s\n", ROOT))
cat(sprintf("  %d audit parquets in list `A`; metadata in `meta` (%d series)\n",
            length(A), nrow(meta)))
cat("  short names: p10, fdr, spa, eqi, s2man, m3, usedv, eqfc\n")
cat("  helpers: ls_artifacts(), peek(x), cols(x), get_raw(id), src_R()\n")
