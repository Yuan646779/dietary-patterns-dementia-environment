# ============================================================================
# Build the final master database: v3
#
# Revision:
#   Dementia outcomes are read directly from the pre-calculated country-year
#   60+ internally age-standardized file:
#   GBD_dementia_ASR60plus_country_year.csv
#
# This script does not re-aggregate age-specific dementia data and does not
# derive population weights from Number/Rate. Existing v2 output is replaced.
# ============================================================================

required_packages <- c("data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "))
}
library(data.table)

script_version <- "v3_ASR60plus_direct_country_year_merge"
analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
master_dir <- file.path(analysis_dir, "02_processed", "Master")
gbd_file <- file.path(
  analysis_dir, "01_raw", "GBD",
  "GBD_dementia_ASR60plus_country_year.csv"
)
master_file_in <- file.path(
  master_dir, "GDD_final_master_database_185_country_year_v1.csv"
)
master_file_out <- file.path(
  master_dir, "GDD_final_master_database_185_country_year_v2.csv"
)

dir.create(master_dir, recursive = TRUE, showWarnings = FALSE)

missingness_file <- file.path(
  master_dir, "QC_GDD_final_master_database_v2_missingness.csv"
)
key_qc_file <- file.path(
  master_dir, "QC_GDD_final_master_database_v2_key.csv"
)
health_qc_file <- file.path(
  master_dir, "QC_GBD_dementia_ASR60plus_direct_merge_v3.csv"
)
negative_file <- file.path(
  master_dir, "QC_GDD_final_master_v2_negative_competing_mortality.csv"
)
manifest_file <- file.path(
  master_dir, "QC_GDD_final_master_v2_input_manifest.csv"
)

study_years <- c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L)
expected_countries <- 185L
expected_rows <- expected_countries * length(study_years)

file_required <- function(path) {
  if (!file.exists(path)) stop("Required input not found: ", path)
  path
}

normalise_names <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

# ---- 1. Read and validate the existing non-dementia panel -----------------

master <- fread(
  file_required(master_file_in),
  na.strings = c("", "NA", "NaN", "."),
  encoding = "UTF-8"
)
setnames(master, names(master), normalise_names(names(master)))

required_master <- c(
  "iso3", "year", "superregion2",
  "allcause_death_rate_60plus", "dementia_death_rate_60plus"
)
missing_master <- setdiff(required_master, names(master))
if (length(missing_master)) {
  stop("Existing master is missing: ", paste(missing_master, collapse = ", "))
}

master[, `:=`(
  iso3 = toupper(trimws(as.character(iso3))),
  year = as.integer(year)
)]
master <- master[year %in% study_years]

if (anyDuplicated(master[, .(iso3, year)])) {
  stop("Duplicate existing-master iso3-year keys detected.")
}
if (nrow(master) != expected_rows || uniqueN(master$iso3) != expected_countries) {
  stop("Existing master must contain exactly 185 countries and 1,295 country-year rows.")
}
if (!identical(sort(unique(master$year)), study_years)) {
  stop("Existing master years do not match the seven study years.")
}

# Preserve the predefined competitive-mortality definition.
master[, competitive_mortality_rate_60plus :=
         allcause_death_rate_60plus - dementia_death_rate_60plus]
master[, competitive_mortality_definition :=
         "All-cause mortality rate minus dementia mortality rate"]

# ---- 2. Read the direct country-year ASR60plus outcomes --------------------

asr <- fread(
  file_required(gbd_file),
  na.strings = c("", "NA", "NaN", "."),
  encoding = "UTF-8"
)
setnames(asr, names(asr), normalise_names(names(asr)))

required_asr <- c(
  "iso3", "year",
  "dementia_daly_asr60plus",
  "dementia_daly_asr60plus_lower",
  "dementia_daly_asr60plus_upper",
  "dementia_prevalence_asr60plus",
  "dementia_prevalence_asr60plus_lower",
  "dementia_prevalence_asr60plus_upper"
)
missing_asr <- setdiff(required_asr, names(asr))
if (length(missing_asr)) {
  stop(
    "ASR60plus file is missing: ",
    paste(missing_asr, collapse = ", ")
  )
}

asr[, `:=`(
  iso3 = toupper(trimws(as.character(iso3))),
  year = as.integer(year)
)]
asr <- asr[year %in% study_years]

asr_duplicates <- asr[, .N, by = .(iso3, year)][N > 1L]
if (nrow(asr_duplicates)) {
  fwrite(asr_duplicates, health_qc_file, bom = TRUE)
  stop("Duplicate ASR60plus iso3-year keys detected.")
}

asr_missing_keys <- fsetdiff(
  master[, .(iso3, year)], asr[, .(iso3, year)],
  all = TRUE
)
asr_extra_keys <- fsetdiff(
  asr[, .(iso3, year)], master[, .(iso3, year)],
  all = TRUE
)

if (nrow(asr_missing_keys) || nrow(asr_extra_keys)) {
  fwrite(
    rbindlist(
      list(
        asr_missing_keys[, key_status := "missing_from_ASR"],
        asr_extra_keys[, key_status := "extra_in_ASR"]
      ),
      fill = TRUE
    ),
    health_qc_file,
    bom = TRUE
  )
  stop("ASR60plus country-year keys do not exactly match the master panel.")
}

asr_values <- asr[, c("iso3", "year", setdiff(required_asr, c("iso3", "year"))), with = FALSE]
asr_qc <- asr_values[, .(
  n_rows = .N,
  n_countries = uniqueN(iso3),
  n_years = uniqueN(year),
  n_missing_daly = sum(is.na(dementia_daly_asr60plus)),
  n_missing_prevalence = sum(is.na(dementia_prevalence_asr60plus)),
  n_nonpositive_daly = sum(!is.na(dementia_daly_asr60plus) & dementia_daly_asr60plus <= 0),
  n_nonpositive_prevalence = sum(!is.na(dementia_prevalence_asr60plus) & dementia_prevalence_asr60plus <= 0)
)]
fwrite(asr_qc, health_qc_file, bom = TRUE)

outcome_columns <- setdiff(required_asr, c("iso3", "year"))
if (any(asr_values[, lapply(.SD, function(z) anyNA(z)), .SDcols = outcome_columns])) {
  stop("Missing dementia outcomes remain in the ASR60plus file.")
}
if (any(asr_values[, lapply(.SD, function(z) any(!is.finite(z))), .SDcols = outcome_columns])) {
  stop("Non-finite dementia outcomes found in the ASR60plus file.")
}
if (any(asr_values$dementia_daly_asr60plus <= 0) ||
    any(asr_values$dementia_prevalence_asr60plus <= 0)) {
  stop("Non-positive dementia DALY or prevalence rates found.")
}

# ---- 3. Merge ASR60plus outcomes into the master panel --------------------

# Remove any stale dementia outcome columns if the input v1 panel already
# contains them; the direct ASR60plus file is the sole outcome source here.
stale_outcome_columns <- intersect(outcome_columns, names(master))
if (length(stale_outcome_columns)) {
  master[, (stale_outcome_columns) := NULL]
}
master <- merge(master, asr_values, by = c("iso3", "year"), all.x = TRUE)
setorder(master, iso3, year)

if (nrow(master) != expected_rows || anyDuplicated(master[, .(iso3, year)])) {
  stop("Final master key QC failed after ASR60plus merge.")
}
if (anyNA(master[, ..outcome_columns])) {
  stop("Missing ASR60plus outcomes remain after merging.")
}

negative_competing <- master[
  competitive_mortality_rate_60plus < -1e-10,
  .(
    iso3, year, allcause_death_rate_60plus,
    dementia_death_rate_60plus,
    competitive_mortality_rate_60plus
  )
]

key_qc <- master[, .(
  n_rows = .N,
  n_countries = uniqueN(iso3),
  n_years = uniqueN(year),
  min_year = min(year),
  max_year = max(year)
), by = .(iso3, year)]

if (nrow(key_qc) != expected_rows || any(key_qc$n_rows != 1L)) {
  stop("Final key QC failed.")
}

qc_cols <- setdiff(names(master), c("iso3", "year"))
missingness_qc <- data.table(
  variable = qc_cols,
  n_missing = vapply(master[, ..qc_cols], function(z) sum(is.na(z)), integer(1)),
  missing_percent = vapply(master[, ..qc_cols], function(z) 100 * mean(is.na(z)), numeric(1))
)

# ---- 4. Write outputs; v2 is intentionally overwritten --------------------

fwrite(master, master_file_out, bom = TRUE, na = "NA")
fwrite(missingness_qc, missingness_file, bom = TRUE)
fwrite(key_qc, key_qc_file, bom = TRUE)
fwrite(negative_competing, negative_file, bom = TRUE)

manifest <- data.table(
  component = c("Existing master input", "Direct dementia ASR60plus input"),
  file = c(basename(master_file_in), basename(gbd_file)),
  directory = c(dirname(master_file_in), dirname(gbd_file)),
  role = c(
    "Dietary, environmental, covariate, and mortality country-year panel",
    "Pre-calculated country-year dementia DALY and prevalence rates internally age-standardized among adults aged 60+"
  ),
  processing_note = c(
    "Validated before merge",
    "Merged directly; no re-aggregation from age-specific Number/Rate data"
  )
)
fwrite(manifest, manifest_file, bom = TRUE)

message("Completed: final master database rebuilt with direct ASR60plus outcomes")
message("Script version: ", script_version)
message("Countries: ", uniqueN(master$iso3))
message("Years: ", paste(sort(unique(master$year)), collapse = ", "))
message("Rows: ", format(nrow(master), big.mark = ","))
message("Dementia source: ", basename(gbd_file))
message("Dementia outcomes: direct country-year ASR60plus values")
message("Competitive mortality: all-cause mortality rate minus dementia mortality rate")
message("Negative competitive-mortality records: ", nrow(negative_competing))
message("Master output overwritten: ", master_file_out)
message("Missingness QC: ", missingness_file)
message("Key QC: ", key_qc_file)
message("Health merge QC: ", health_qc_file)
message("Input manifest: ", manifest_file)
