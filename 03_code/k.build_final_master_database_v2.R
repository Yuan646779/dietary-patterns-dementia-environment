# k.build_final_master_database_v1.R
# Build the final 185-country x 7-year analysis master database.
# Competitive mortality is generated here as:
#   all-cause mortality rate - dementia mortality rate.

required_packages <- c("data.table")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required package: data.table")
library(data.table)

analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
gdd_dir <- file.path(analysis_dir, "02_processed", "GDD")
env_dir <- file.path(analysis_dir, "02_processed", "Environment", "LCA_processed")
cov_dir <- file.path(analysis_dir, "02_processed", "Covariates")
output_dir <- file.path(analysis_dir, "02_processed", "Master")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

study_years <- c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L)
expected_countries <- 185L
expected_rows <- expected_countries * length(study_years)

required_file <- function(path) {
  if (!file.exists(path)) stop("Required input not found: ", path)
  path
}

read_keyed <- function(path, label) {
  x <- fread(required_file(path), encoding = "UTF-8")
  if (!all(c("iso3", "year") %in% names(x))) stop(label, " must contain iso3 and year.")
  x[, `:=`(iso3 = toupper(trimws(as.character(iso3))), year = as.integer(year))]
  x <- x[year %in% study_years]
  if (anyDuplicated(x[, .(iso3, year)])) stop(label, " contains duplicate iso3-year keys.")
  x
}

rename_if_present <- function(x, old, new) {
  if (old %in% names(x) && !(new %in% names(x))) setnames(x, old, new)
  x
}

# ---- 1. GDD outcomes and population-weighted estimates --------------------
gdd_long <- fread(required_file(file.path(gdd_dir, "GDD_country_60plus_all47_country_year.csv")), encoding = "UTF-8")
setnames(gdd_long, names(gdd_long), tolower(names(gdd_long)))
gdd_long[, `:=`(iso3 = toupper(trimws(as.character(iso3))), year = as.integer(year))]
gdd_long <- gdd_long[year %in% study_years]
if (anyDuplicated(gdd_long[, .(iso3, year, variable_code)])) stop("GDD base contains duplicate country-year-variable keys.")
if (!all(c("estimate", "variable_code") %in% names(gdd_long))) stop("GDD base must contain estimate and variable_code.")

gdd_core <- dcast(gdd_long, iso3 + year ~ variable_code, value.var = "estimate")
gdd_core_cols <- setdiff(names(gdd_core), c("iso3", "year"))
setnames(gdd_core, gdd_core_cols, paste0("gdd_", gdd_core_cols))
gdd_superregion <- unique(gdd_long[, .(iso3, year, superregion2)], by = c("iso3", "year"))

# ---- 2. Four dietary-pattern scores ----------------------------------------
diet <- read_keyed(file.path(gdd_dir, "GDD_four_dietary_patterns_country_year_v1.csv"), "Dietary-pattern file")
diet_cols <- setdiff(names(diet), c("iso3", "year"))
setnames(diet, diet_cols, paste0("diet_", diet_cols))

# ---- 3. Fifteen environmental footprints -----------------------------------
foot <- read_keyed(file.path(env_dir, "GDD15_environmental_footprints_country_year_v2.csv"), "Environmental-footprint file")
foot_cols <- setdiff(names(foot), c("iso3", "year"))
setnames(foot, foot_cols, paste0("env_", foot_cols))

# ---- 4. Harmonized time-varying covariates ---------------------------------
covariates <- read_keyed(file.path(cov_dir, "time_varying_covariates_185_country_year_wide_v1.csv"), "Covariate file")

# ---- 5. Merge all components on the complete country-year key --------------
master <- Reduce(function(x, y) merge(x, y, by = c("iso3", "year"), all = TRUE),
                 list(gdd_core, gdd_superregion, diet, foot, covariates))
setorder(master, iso3, year)

if (nrow(master) != expected_rows) {
  stop("Final master database has ", nrow(master), " rows; expected ", expected_rows, ".")
}
if (uniqueN(master$iso3) != expected_countries) stop("Final master database does not contain 185 countries.")
if (!identical(sort(unique(master$year)), study_years)) stop("Final master years do not match the study years.")

# ---- 6. Generate competitive mortality -------------------------------------
required_mortality <- c("allcause_death_rate_60plus", "dementia_death_rate_60plus")
missing_mortality <- setdiff(required_mortality, names(master))
if (length(missing_mortality)) stop("Cannot calculate competitive mortality; missing: ", paste(missing_mortality, collapse = ", "))

master[, competitive_mortality_rate_60plus := allcause_death_rate_60plus - dementia_death_rate_60plus]
master[, competitive_mortality_definition := "All-cause mortality rate minus dementia mortality rate"]
negative_competing <- master[competitive_mortality_rate_60plus < -1e-10,
                            .(iso3, year, allcause_death_rate_60plus, dementia_death_rate_60plus,
                              competitive_mortality_rate_60plus)]

# Keep only the core analysis fields plus complete provenance/derived fields.
setcolorder(master, c("iso3", "year", "superregion2", "competitive_mortality_rate_60plus",
                      setdiff(names(master), c("iso3", "year", "superregion2", "competitive_mortality_rate_60plus"))))

# ---- 7. Quality control ----------------------------------------------------
missingness_qc <- data.table(
  variable = setdiff(names(master), c("iso3", "year", "competitive_mortality_definition")),
  n_missing = vapply(master[, setdiff(names(master), c("iso3", "year", "competitive_mortality_definition")), with = FALSE],
                     function(z) sum(is.na(z)), integer(1)),
  missing_percent = vapply(master[, setdiff(names(master), c("iso3", "year", "competitive_mortality_definition")), with = FALSE],
                           function(z) 100 * mean(is.na(z)), numeric(1))
)

key_qc <- master[, .(n_rows = .N, n_countries = uniqueN(iso3), n_years = uniqueN(year)), by = .(iso3, year)]
if (nrow(key_qc) != expected_rows || any(key_qc$n_rows != 1L)) stop("Final master key QC failed.")

master_file <- file.path(output_dir, "GDD_final_master_database_185_country_year_v1.csv")
qc_file <- file.path(output_dir, "QC_GDD_final_master_database_v1.csv")
negative_file <- file.path(output_dir, "QC_GDD_final_master_negative_competing_mortality_v1.csv")
manifest_file <- file.path(output_dir, "QC_GDD_final_master_input_manifest_v1.csv")

fwrite(master, master_file, bom = TRUE, na = "NA")
fwrite(missingness_qc, qc_file, bom = TRUE)
fwrite(negative_competing, negative_file, bom = TRUE)
manifest <- data.table(
  component = c("GDD base", "Dietary patterns", "Environmental footprints", "Time-varying covariates"),
  file = c("GDD_country_60plus_all47_country_year.csv", "GDD_four_dietary_patterns_country_year_v1.csv", "GDD15_environmental_footprints_country_year_v2.csv", "time_varying_covariates_185_country_year_wide_v1.csv"),
  directory = c(gdd_dir, gdd_dir, env_dir, cov_dir),
  role = c("47 GDD country-year estimates", "mPHDI, MED, DASH, AHEI", "15-food-group environmental footprints", "11 time-varying covariates")
)
fwrite(manifest, manifest_file, bom = TRUE)

message("Completed: final master database")
message("Countries: ", uniqueN(master$iso3))
message("Years: ", paste(sort(unique(master$year)), collapse = ", "))
message("Rows: ", format(nrow(master), big.mark = ","))
message("Competitive mortality: all-cause minus dementia mortality")
message("Negative competitive-mortality records: ", nrow(negative_competing))
message("Master output: ", master_file)
message("Missingness QC: ", qc_file)
message("Input manifest: ", manifest_file)
message("Next: review the master missingness QC and negative competitive-mortality QC before analysis models.")
