#!/usr/bin/env Rscript
# build_annex_table.R — R5/M5 annex: the audited table of variables EFFECTIVELY USED
# in a model under the methodology upgrade (U6). One row per variable, with academic
# source (RESEARCH.md §2 sources S1-S21 / "Standard practitioner").
#
# Used set = the 114-var U6 Scenario-2 pool (the canonical oracle
# data/audit/s2_pool_manifest_2026-06-06.parquet) UNION the S1 theory pools, the CP
# forward-rate inputs, and the three Part-II SPF surprises. The S2 pool is the
# single source of truth for membership/family/formulation; the manifest carries
# `family` and `ratio_or_level` so name/formulation need no re-derivation.
#
# Audit gate (manuscript-completeness standard): every used variable is listed and
# every listed variable is used (0 stale / 0 missing), and every row carries a
# non-empty academic source and pool membership. Output:
#   data/audit/annex_used_variables_<today>.parquet
suppressMessages({library(arrow); library(dplyr)})

meta  <- read_parquet("data/metadata/factor_metadata.parquet")
s2man <- read_parquet("data/audit/s2_pool_manifest_2026-06-06.parquet")  # 114 vars
TODAY <- Sys.Date()

# ---- pool membership ---------------------------------------------------------
s2_vars <- s2man$factor_id                                       # 114 (U6 S2 pool)
dash13  <- s2man$factor_id[s2man$subset == "dashboard"]          # 13 practitioner constructs (M1)
in_s1_eq <- c("CPIAUCSL","INDPRO","PCEPI","DGS3MO","DGS10","TERM_SPREAD",
              "DEFAULT_SPREAD","VIXCLS","EXPINF10YR")            # Welch-Goyal ∪ Chen-Roll-Ross
in_ff    <- c("TERM_SPREAD","DEFAULT_SPREAD")                    # Fama-French 1989
in_cp    <- c("DGS1","DGS2","DGS3","DGS5","DGS10")              # Cochrane-Piazzesi forward inputs
surprises<- c("CPI_SURPRISE","INDPRO_SURPRISE","UNRATE_SURPRISE")
derived  <- c("TERM_SPREAD","DEFAULT_SPREAD", surprises)         # not raw factor_ids in meta

pools_of <- function(id) {
  p <- character(0)
  if (id %in% s2_vars)   p <- c(p, "S2")
  if (id %in% in_s1_eq)  p <- c(p, "S1-eq")
  if (id %in% in_ff)     p <- c(p, "S1-bond/FF")
  if (id %in% in_cp)     p <- c(p, "S1-bond/CP")
  if (id %in% surprises) p <- c(p, "Ch3-surprise")
  if (id %in% dash13)    p <- c(p, "dashboard")
  paste(p, collapse = ", ")
}

# ---- academic source: RESEARCH.md §2 (S1-S21) family anchors + per-code overrides
fam_src <- c(
  CREDIT_SPREADS = "Fama & French (1989)",
  YIELD_CURVE = "Adrian, Crump & Moench (2013)",
  MONETARY_POLICY_RATES = "Welch & Goyal (2008); FRB H.15",
  INFLATION = "Chen, Roll & Ross (1986)",
  GROWTH = "Chen, Roll & Ross (1986)",
  LABOR_MARKET = "Boyd, Hu & Jagannathan (2005)",
  HOUSING = "Standard practitioner (Census / S&P / FHFA)",
  CONSUMPTION_WEALTH = "Lettau & Ludvigson (2001)",
  EQUITY_VALUATION = "Campbell & Shiller (1988)",
  EQUITY_VOLATILITY = "CBOE (practitioner)",
  EQUITY_SENTIMENT = "Baker & Wurgler (2006); Hull & Qiao (2017)",
  FINANCIAL_CONDITIONS_COMPOSITE = "Brave & Butters (2011)",
  FX_CURRENCY = "Standard practitioner (FRB H.10)",
  COMMODITIES = "Bernanke (2016); standard practitioner",
  POLICY_UNCERTAINTY_GEOPOLITICAL = "Baker, Bloom & Davis (2016); Caldara & Iacoviello (2022)",
  FORECASTS = "Philadelphia Fed SPF (Croushore 1993)")
override <- c(
  # precise per-series anchors (RESEARCH.md §2)
  AAA="Fama & French (1989)", BAA="Fama & French (1989)",
  DEFAULT_SPREAD="Fama & French (1989)", TERM_SPREAD="Fama & French (1989)",
  AAII_BULLBEAR="Hull & Qiao (2017)",
  WURGLER_SENT="Baker & Wurgler (2006); Hull & Qiao (2017)",
  SHILLER_CAPE="Campbell & Shiller (1988)", SHILLER_PRICE="Campbell & Shiller (1988)",
  SP500TR="Standard practitioner (S&P)",
  SKEW="CBOE (practitioner)", VIXCLS="CBOE (practitioner)",
  GVZCLS="CBOE (practitioner)", OVXCLS="CBOE (practitioner)",
  VVIXCLS="CBOE (practitioner)", VXEEMCLS="CBOE (practitioner)",
  ANFCI="Brave & Butters (2011)", NFCI="Brave & Butters (2011)",
  NFCICREDIT="Brave & Butters (2011)", NFCILEVERAGE="Brave & Butters (2011)",
  NFCIRISK="Brave & Butters (2011)",
  KCFSI="Hakkio & Keeton (2009)", STLFSI4="Kliesen & Smith (2010)",
  ADS_INDEX="Aruoba, Diebold & Scotti (2009)",
  CFNAI="Brave, Butters & Justiniano (2019)",
  GACDFSA066MSFRBPHI="Philadelphia Fed (practitioner)",
  GACDISA066MSFRBNY="NY Fed Empire State (practitioner)",
  INDPRO="Chen, Roll & Ross (1986)",
  INDPRO_SURPRISE="Chen, Roll & Ross (1986) MP channel; SPF",
  HOUST="U.S. Census (practitioner)", PERMIT="U.S. Census (practitioner)",
  HSN1F="U.S. Census (practitioner)", CSUSHPINSA="Case & Shiller; S&P (practitioner)",
  USSTHPI="FHFA (practitioner)", MORTGAGE30US="Freddie Mac PMMS (practitioner)",
  CORESTICKM159SFRBATL="Bryan & Meyer (2010)",
  CPIAUCSL="Chen, Roll & Ross (1986); BLS", CPILFESL="Chen, Roll & Ross (1986); BLS",
  CPI_SURPRISE="Chen, Roll & Ross (1986) UI channel; SPF",
  PCEPILFE="Chen, Roll & Ross (1986); BEA",
  EXPINF10YR="Haubrich, Pennacchi & Ritchken (2012)",
  DFII10="Gürkaynak, Sack & Wright (2010)",
  PCEPI="Chen, Roll & Ross (1986); BEA", PCEC96="Lettau & Ludvigson (2001)",
  A229RX0="Lettau & Ludvigson (2001)",
  PAYEMS="BLS (practitioner)", UNRATE="BLS (practitioner)",
  UNRATE_SURPRISE="Boyd, Hu & Jagannathan (2005); SPF",
  TCU="Standard practitioner (FRB G.17)",
  OECDPRINTO01GYSAM="OECD (practitioner)",
  GDP="BEA; Chen, Roll & Ross (1986)", GDPC1="BEA; Chen, Roll & Ross (1986)",
  GDPNOW="Atlanta Fed GDPNow (Higgins 2014)",
  DFF="Booth & Booth (1997); FRB H.15",
  DGS1="FRB H.15 (Treasury yield)", DGS2="FRB H.15 (Treasury yield)",
  DGS3="FRB H.15 (Treasury yield)", DGS5="FRB H.15 (Treasury yield)",
  DGS10="Welch & Goyal (2008); FRB H.15", DGS3MO="Welch & Goyal (2008); FRB H.15",
  NYFED_ACMRNY10="Adrian, Crump & Moench (2013)",
  NYFED_ACMTP10="Adrian, Crump & Moench (2013)",
  RECPROUSM156N="Estrella & Mishkin (1998)",
  EMVOVERALLEMV="Baker, Bloom, Davis & Kost (2019)",
  IACO_GPR="Caldara & Iacoviello (2022)", IACO_GPRA="Caldara & Iacoviello (2022)",
  IACO_GPRT="Caldara & Iacoviello (2022)",
  DCOILBRENTEU="Bernanke (2016); EIA", GOLD="Standard practitioner",
  AGG="Standard practitioner (iShares)", HYG="Standard practitioner (iShares)",
  LQD="Standard practitioner (iShares)", EEM="Standard practitioner (iShares)",
  EFA="Standard practitioner (iShares)",
  DEXJPUS="FRB H.10 (practitioner)", DEXUSEU="FRB H.10 (practitioner)",
  DTWEXBGS="FRB H.10 (practitioner)", DTWEXEMEGS="FRB H.10 (practitioner)",
  "00XEFDEZ19M086NEST"="Eurostat (practitioner)")
# BBD EPU/disagreement family -> Baker, Bloom & Davis (2016)
acad_of <- function(id, fam) {
  if (id %in% names(override)) return(unname(override[id]))
  if (grepl("^BBD_", id))      return("Baker, Bloom & Davis (2016)")
  if (id %in% dash13)          return("T. Rowe Price (practitioner dashboard)")
  if (fam %in% names(fam_src)) return(unname(fam_src[fam]))
  ""
}

# ---- assemble rows -----------------------------------------------------------
src_short <- function(u) ifelse(is.na(u) | u == "", "constructed",
  sub("^https?://(www[.])?([^/]+).*", "\\2", u))
raw_ids <- setdiff(unique(c(s2_vars, in_s1_eq, in_cp)), derived)   # raw factor_ids (in meta)
raw_ids <- intersect(raw_ids, meta$factor_id)
fam_lookup <- setNames(s2man$family, s2man$factor_id)
rl_lookup  <- setNames(s2man$ratio_or_level, s2man$factor_id)

raw <- meta %>% filter(factor_id %in% raw_ids) %>%
  transmute(code = factor_id, name = factor_name,
            family = ifelse(factor_id %in% names(fam_lookup), fam_lookup[factor_id], family),
            data_source = src_short(source_url),
            formulation = ifelse(code %in% names(rl_lookup) & rl_lookup[code] == "ratio",
                                 "Y/Y log difference, then trailing 60-mo z-score",
                                 "level/spread, trailing 60-mo z-score"))

# derived constructs (dashboard + spreads + surprises) — not raw meta rows
dash_rows <- tibble(code = dash13,
  name = paste0("Practitioner dashboard construct (", dash13, ")"),
  family = unname(fam_lookup[dash13]),
  data_source = "constructed (TRP dashboard)",
  formulation = ifelse(rl_lookup[dash13] == "ratio",
                       "relative-value / Y/Y construct, 60-mo z-score",
                       "level construct, 60-mo z-score"))
spread_rows <- tibble(
  code = c("TERM_SPREAD","DEFAULT_SPREAD","CPI_SURPRISE","INDPRO_SURPRISE","UNRATE_SURPRISE"),
  name = c("Term spread (10Y minus 3M Treasury)","Default spread (Baa minus Aaa)",
           "CPI inflation surprise (realized minus SPF)",
           "Industrial-production surprise (realized minus SPF)",
           "Unemployment surprise (realized minus SPF)"),
  family = c("YIELD_CURVE","CREDIT_SPREADS","INFLATION","GROWTH","LABOR_MARKET"),
  data_source = c("FRED (DGS10, DGS3MO)","FRED (BAA, AAA)",
                  "FRED ALFRED + Philadelphia Fed SPF","FRED ALFRED + Philadelphia Fed SPF",
                  "FRED ALFRED + Philadelphia Fed SPF"),
  formulation = c("derived spread (level); 60-mo z-score","derived spread (level); 60-mo z-score",
                  "realized first-vintage minus SPF consensus nowcast; 60-mo z-score",
                  "realized first-vintage minus SPF consensus nowcast; 60-mo z-score",
                  "realized first-vintage minus SPF consensus nowcast; 60-mo z-score"))

tab <- bind_rows(raw, dash_rows, spread_rows) %>%
  mutate(pools = vapply(code, pools_of, ""),
         academic_source = mapply(acad_of, code, family)) %>%
  arrange(family, code) %>%
  mutate(No = row_number()) %>%
  select(No, code, name, family, data_source, formulation, academic_source, pools)

# ---- AUDIT: 0 stale / 0 missing ---------------------------------------------
target_set <- sort(unique(c(s2_vars, in_s1_eq, in_ff, in_cp, surprises, dash13)))
listed_set <- sort(unique(tab$code))
missing <- setdiff(target_set, listed_set)   # used but not listed
stale   <- setdiff(listed_set, target_set)   # listed but not used
if (length(missing)) cat("MISSING:", paste(missing, collapse=", "), "\n")
if (length(stale))   cat("STALE:",   paste(stale,   collapse=", "), "\n")
stopifnot("0 missing (every used variable is listed)" = length(missing) == 0L,
          "0 stale (every listed variable is used)"   = length(stale)   == 0L,
          "every row has an academic source" = all(nzchar(tab$academic_source)),
          "every row has a pool"             = all(nzchar(tab$pools)))

out <- sprintf("data/audit/annex_used_variables_%s.parquet", TODAY)
write_parquet(tab, out)
cat(sprintf("OK: %d used variables (%d S2 + S1/derived/surprise); 0 stale / 0 missing → %s\n",
            nrow(tab), length(s2_vars), out))
print(as.data.frame(tab %>% count(family, name = "n_vars")), row.names = FALSE)
