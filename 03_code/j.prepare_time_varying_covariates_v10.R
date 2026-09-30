# j.prepare_time_varying_covariates_v1.R
# Prepare harmonized time-varying covariates for the 185-country GDD analysis
# Study years: 1990, 1995, 2000, 2005, 2010, 2015, 2018

required_packages <- c("data.table", "readxl", "countrycode")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages: ", paste(missing_packages, collapse = ", "))
library(data.table)
library(readxl)

analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
raw_dir <- file.path(analysis_dir, "01_raw")
processed_dir <- file.path(analysis_dir, "02_processed")
gdd_processed_dir <- file.path(processed_dir, "GDD")
output_dir <- file.path(processed_dir, "Covariates")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

study_years <- c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L)

file_required <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path)
  path
}

clean_name <- function(x) {
  x <- trimws(as.character(x))
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- tolower(gsub("[^a-z0-9]", "", x))
  x
}

make_skeleton <- function() {
  base_file <- file.path(gdd_processed_dir, "GDD_four_dietary_patterns_country_year_v1.csv")
  if (!file.exists(base_file)) stop("GDD dietary-pattern base file not found: ", base_file)
  base <- fread(base_file, encoding = "UTF-8")
  base[, .(iso3 = as.character(iso3), year = as.integer(year))][
    , .(iso3, year)
  ][, .SD[year %in% study_years]][, unique(.SD)]
}

skeleton <- make_skeleton()
if (nrow(skeleton) != 185L * length(study_years)) {
  stop("Unexpected country-year skeleton: ", nrow(skeleton), " rows; expected ", 185L * length(study_years), ".")
}

add_covariate <- function(dt, value_col, covariate_name, source, unit, years = study_years) {
  dt <- as.data.table(dt)
  keep <- intersect(c("iso3", "year", value_col), names(dt))
  if (!all(c("iso3", "year", value_col) %in% keep)) stop("Cannot standardize ", covariate_name)
  out <- dt[, .(iso3 = as.character(iso3), year = as.integer(year), value = as.numeric(get(value_col)))][year %in% years]
  out <- out[!is.na(iso3) & !is.na(year)]
  if (anyDuplicated(out, by = c("iso3", "year"))) stop("Duplicate iso3-year records in ", covariate_name)
  out[, `:=`(covariate = covariate_name, source = source, unit = unit)]
  out[, .(iso3, year, covariate, value, unit, source)]
}

gbd_crosswalk_file <- file.path(gdd_processed_dir, "GBD_location_to_ISO3_crosswalk.csv")
crosswalk <- fread(file_required(gbd_crosswalk_file), encoding = "UTF-8")
setnames(crosswalk, names(crosswalk), tolower(names(crosswalk)))
cw_id_col <- names(crosswalk)[tolower(names(crosswalk)) %in% c("location_id", "gbd_location_id")][1L]
iso_col <- names(crosswalk)[tolower(names(crosswalk)) %in% c("iso3", "iso3c", "iso3_code")][1L]
name_col <- names(crosswalk)[tolower(names(crosswalk)) %in% c("location_name", "gbd_location_name", "location")][1L]
if (is.na(iso_col) || is.na(name_col)) stop("Crosswalk must contain ISO3 and location-name columns.")
if (!is.na(cw_id_col)) {
  crosswalk <- unique(crosswalk[, .(location_id = as.integer(get(cw_id_col)), iso3 = toupper(get(iso_col)), location_name = as.character(get(name_col)))])
} else {
  crosswalk <- unique(crosswalk[, .(iso3 = toupper(get(iso_col)), location_name = as.character(get(name_col)))])
}

# ---- 1. GBD SDI ------------------------------------------------------------
# Prefer the already harmonized country-year SDI file used in the earlier
# analysis. The newer 1950-2023 file is retained as a fallback.
sdi_candidates <- c(
  file.path(raw_dir, "GBD", "SDI_country_year_1990_2018.csv"),
  file.path(raw_dir, "GBD", "SDI", "IHME_GBD_2023_SDI_1950_2023_Y2025M10D12.csv")
)
sdi_existing <- sdi_candidates[file.exists(sdi_candidates)]
if (!length(sdi_existing)) stop("No SDI input found. Expected SDI_country_year_1990_2018.csv or the IHME GBD SDI file.")
sdi_file <- sdi_existing[1L]
sdi <- fread(file_required(sdi_file), encoding = "UTF-8")
sdi_names <- tolower(names(sdi))
setnames(sdi, names(sdi), sdi_names)
if ("year_id" %in% names(sdi)) {
  sdi[, year := as.integer(year_id)]
} else if ("year" %in% names(sdi)) {
  sdi[, year := as.integer(year)]
} else {
  stop("SDI file must contain year_id or year.")
}
sdi <- sdi[age_group_name %in% c("All Ages", "All ages") & sex %in% c("Both", "Both sexes") & location_name != "Global"]
sdi <- if ("location_id" %in% names(sdi) && "location_id" %in% names(crosswalk)) {
  merge(sdi, crosswalk[, .(location_id, iso3)], by = "location_id", all.x = TRUE)
} else {
  merge(sdi, crosswalk[, .(location_name, iso3)], by = "location_name", all.x = TRUE)
}
sdi <- sdi[iso3 %in% skeleton$iso3]
sdi <- sdi[, .(mean_value = mean(mean_value, na.rm = TRUE)), by = .(iso3, year)]
sdi_long <- add_covariate(sdi, "mean_value", "sdi", "IHME GBD 2023 SDI", "SDI", years = study_years)

# ---- 2. GBD competing mortality and diabetes prevalence --------------------
read_gbd_outcome <- function(path, covariate_name, unit) {
  x <- fread(file_required(path), encoding = "UTF-8")
  setnames(x, names(x), tolower(names(x)))
  x <- x[age_name %in% c("60+ years", "60 plus years", "60+ years and over") &
    sex_name %in% c("Both", "Both sexes")]
  x <- if ("location_id" %in% names(x) && "location_id" %in% names(crosswalk)) {
    merge(x, crosswalk[, .(location_id, iso3)], by = "location_id", all.x = TRUE)
  } else {
    merge(x, crosswalk[, .(location_name, iso3)], by = "location_name", all.x = TRUE)
  }
  add_covariate(x, "val", covariate_name, "IHME GBD 2023", unit)
}

allcause_long <- read_gbd_outcome(
  file.path(raw_dir, "GBD", "GBD_allcause_death_rate_60plus_country_year.csv"),
  "allcause_death_rate_60plus", "deaths per 100,000 population"
)
dementia_long <- read_gbd_outcome(
  file.path(raw_dir, "GBD", "GBD_dementia_death_rate_60plus_country_year.csv"),
  "dementia_death_rate_60plus", "deaths per 100,000 population"
)
diabetes_long <- read_gbd_outcome(
  file.path(raw_dir, "GBD", "GBD2023_diabetes_prevalence_60plus_percent_raw.csv"),
  "diabetes_prevalence_60plus", "percent"
)

# ---- 3. UNDP mean years of schooling ---------------------------------------
mys_file <- file.path(raw_dir, "UNDP", "UNDP_mean_years_schooling_raw.xlsx")
mys <- as.data.table(read_excel(file_required(mys_file), sheet = 1))
setnames(mys, names(mys), tolower(names(mys)))
if ("countryisocode" %in% names(mys)) {
  mys[, iso3 := toupper(trimws(as.character(countryisocode)))]
} else if ("iso3" %in% names(mys)) {
  mys[, iso3 := toupper(trimws(as.character(iso3)))]
} else {
  stop("UNDP schooling file must contain countryIsoCode or iso3.")
}
mys <- mys[year %in% study_years]
mys_long <- add_covariate(mys, "actualvalue", "mean_years_schooling", "UNDP Human Development Index", "years")

# ---- 4. WHO indicators -----------------------------------------------------
who_m49_to_iso3 <- function(x) {
  code <- suppressWarnings(as.numeric(gsub("[^0-9]", "", as.character(x))))
  # Use the broadly supported UN numeric-code dictionary name. Some older
  # countrycode versions expose this field as `un`, not `un_m49`.
  toupper(countrycode::countrycode(code, origin = "un", destination = "iso3c", warn = FALSE))
}

read_who_age_standardized <- function(path, covariate_name, unit, both_sex = TRUE) {
  x <- fread(file_required(path), encoding = "UTF-8")
  setnames(x, names(x), toupper(names(x)))
  if ("DIM_SEX" %in% names(x) && both_sex) x <- x[DIM_SEX %in% c("TOTAL", "BOTH", "BOTH SEXES")]
  x[, iso3 := who_m49_to_iso3(DIM_GEO_CODE_M49)]
  x[, year := as.integer(DIM_TIME)]
  x[, value := as.numeric(RATE_PER_100_N)]
  add_covariate(x, "value", covariate_name, "WHO Data Platform", unit)
}

obesity_long <- read_who_age_standardized(
  file.path(raw_dir, "WHO", "WHO_adult_obesity_age_standardized_raw.csv"),
  "adult_obesity", "percent"
)
hypertension_long <- read_who_age_standardized(
  file.path(raw_dir, "WHO", "WHO_hypertension_prevalence_raw.csv"),
  "hypertension", "percent"
)

read_who_physical_activity <- function(path) {
  x <- fread(file_required(path), encoding = "UTF-8")
  x <- x[Dim1ValueCode %in% c("SEX_BTSX", "BTSX") | Dim1 %in% c("Both sexes", "Both")]
  x[, iso3 := toupper(SpatialDimValueCode)]
  x[, year := as.integer(Period)]
  x[, value := as.numeric(FactValueNumeric)]
  add_covariate(x, "value", "insufficient_physical_activity", "WHO Data Platform", "percent")
}
physical_activity_long <- read_who_physical_activity(file.path(raw_dir, "WHO", "WHO_insufficient_physical_activity_raw.csv"))

# ---- 5. World Bank indicators ----------------------------------------------
read_wdi_wide <- function(path, covariate_name, indicator_code, unit) {
  target_indicator_code <- indicator_code
  # The current World Bank files have already had metadata rows removed;
  # the first line is the actual CSV header.
  x <- fread(file_required(path), encoding = "UTF-8", check.names = FALSE,
             header = TRUE)
  normalized <- gsub("[^a-z0-9]", "", tolower(names(x)))
  country_name_col <- names(x)[normalized == "countryname"][1L]
  country_code_col <- names(x)[normalized == "countrycode"][1L]
  indicator_name_col <- names(x)[normalized == "indicatorname"][1L]
  indicator_code_col <- names(x)[normalized == "indicatorcode"][1L]
  if (any(is.na(c(country_name_col, country_code_col, indicator_name_col, indicator_code_col)))) {
    stop("World Bank header columns not found in ", path, ". Detected: ", paste(names(x), collapse = ", "))
  }
  setnames(x, c(country_name_col, country_code_col, indicator_name_col, indicator_code_col),
           c("country_name", "country_code", "indicator_name", "indicator_code"))
  year_cols <- intersect(as.character(study_years), names(x))
  if (!length(year_cols)) stop("No study-year columns found in ", path)
  out <- melt(x, id.vars = c("country_name", "country_code", "indicator_name", "indicator_code"),
              measure.vars = year_cols, variable.name = "year", value.name = "value", variable.factor = FALSE)
  out[, `:=`(iso3 = toupper(as.character(country_code)), year = as.integer(as.character(year)), value = as.numeric(value))]
  out <- out[indicator_code == target_indicator_code]
  add_covariate(out, "value", covariate_name, "World Bank World Development Indicators", unit)
}

uhc_long <- read_wdi_wide(file.path(raw_dir, "WorldBank", "UHC service coverage index.csv"), "uhc_service_coverage", "SH_UHC_SCI", "index, 0-100")
smoking_long <- read_wdi_wide(file.path(raw_dir, "WorldBank", "WorldBank_SH_PRV_SMOK_raw.csv"), "smoking_prevalence", "SH.PRV.SMOK", "percent")
pm25_long <- read_wdi_wide(file.path(raw_dir, "WorldBank", "WorldBank_PM25_raw.csv"), "pm25", "EN.ATM.PM25.MC.M3", "micrograms per cubic meter")

# ---- 6. Harmonize to the complete 185-country x 7-year panel --------------
covariate_long <- rbindlist(list(sdi_long, allcause_long, dementia_long, diabetes_long,
  mys_long, obesity_long, hypertension_long, physical_activity_long,
  uhc_long, smoking_long, pm25_long), use.names = TRUE, fill = TRUE)
covariate_long[, iso3 := toupper(trimws(iso3))]
covariate_long <- merge(skeleton, covariate_long, by = c("iso3", "year"), all.x = TRUE, allow.cartesian = FALSE)
setorder(covariate_long, iso3, year, covariate)

duplicate_keys <- covariate_long[, .N, by = .(iso3, year, covariate)][N > 1L]
if (nrow(duplicate_keys)) stop("Duplicate country-year-covariate keys after harmonization.")

coverage_qc <- covariate_long[, .(
  n_expected = .N,
  n_available = sum(!is.na(value)),
  n_missing = sum(is.na(value)),
  coverage_percent = 100 * mean(!is.na(value)),
  min_year_available = if (any(!is.na(value))) min(year[!is.na(value)]) else NA_integer_,
  max_year_available = if (any(!is.na(value))) max(year[!is.na(value)]) else NA_integer_
), by = covariate]

covariate_wide <- dcast(covariate_long, iso3 + year ~ covariate, value.var = "value")
source_wide <- dcast(unique(covariate_long[, .(iso3, year, covariate, source)]), iso3 + year ~ covariate, value.var = "source")
unit_wide <- dcast(unique(covariate_long[, .(iso3, year, covariate, unit)]), iso3 + year ~ covariate, value.var = "unit")

long_file <- file.path(output_dir, "time_varying_covariates_185_country_year_long_v1.csv")
wide_file <- file.path(output_dir, "time_varying_covariates_185_country_year_wide_v1.csv")
source_file <- file.path(output_dir, "time_varying_covariates_185_country_year_sources_v1.csv")
unit_file <- file.path(output_dir, "time_varying_covariates_185_country_year_units_v1.csv")
qc_file <- file.path(output_dir, "QC_time_varying_covariates_coverage_v1.csv")
manifest_file <- file.path(output_dir, "QC_time_varying_covariates_file_manifest_v1.csv")

fwrite(covariate_long, long_file, bom = TRUE)
fwrite(covariate_wide, wide_file, bom = TRUE)
fwrite(source_wide, source_file, bom = TRUE)
fwrite(unit_wide, unit_file, bom = TRUE)
fwrite(coverage_qc, qc_file, bom = TRUE)

manifest <- data.table(
  covariate = c("sdi", "allcause_death_rate_60plus", "dementia_death_rate_60plus", "diabetes_prevalence_60plus", "mean_years_schooling", "adult_obesity", "hypertension", "insufficient_physical_activity", "uhc_service_coverage", "smoking_prevalence", "pm25"),
  source_type = c("GBD", "GBD", "GBD", "GBD", "UNDP", "WHO", "WHO", "WHO", "World Bank", "World Bank", "World Bank"),
  years_requested = paste(study_years, collapse = ", "),
  notes = c("Country-level SDI; All Ages, Both sexes; Global excluded", "60+ years, Both sexes, Rate", "60+ years, Both sexes, Rate", "60+ years, Both sexes, Prevalence, Percent", "actualValue used", "Age-standardized adults 18+, TOTAL", "Adults 30-79, TOTAL", "Age-standardized adults 18+, Both sexes", "SH_UHC_SCI", "SH.PRV.SMOK", "EN.ATM.PM25.MC.M3")
)
fwrite(manifest, manifest_file, bom = TRUE)

message("Completed: time-varying covariate harmonization")
message("Countries: ", uniqueN(covariate_wide$iso3))
message("Years: ", paste(sort(unique(covariate_wide$year)), collapse = ", "))
message("Covariates: ", uniqueN(covariate_long$covariate))
message("Long rows: ", format(nrow(covariate_long), big.mark = ","))
message("Wide rows: ", format(nrow(covariate_wide), big.mark = ","))
message("Long output: ", long_file)
message("Wide output: ", wide_file)
message("Coverage QC: ", qc_file)
message("Next: review QC_time_varying_covariates_coverage_v1.csv before merging covariates into the analysis panel.")
