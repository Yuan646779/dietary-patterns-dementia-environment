# ============================================================================
# Primary TWFE analyses of dietary-pattern scores and dementia burden
#
# Study design:
#   185 countries/territories; seven years: 1990, 1995, 2000, 2005, 2010,
#   2015, and 2018; adults aged >=60 years.
#
# Primary health analysis:
#   Exact 10-year exposure lag; dementia DALY rate is the primary outcome and
#   dementia prevalence rate is complementary.
#
# Model sequence:
#   Model 0: country FE + year FE
#   Model 1: Model 0 + SDI (primary prespecified model)
#   Model 2: Model 1 + mean years of schooling
#   Model 3: Model 2 + smoking, obesity, hypertension, diabetes, and physical
#            activity
#   Model 4: Model 3 + PM2.5 + UHC service coverage
#   Model 5: Model 4 + competitive mortality rate
#
# Missing data:
#   Each model uses its own complete-case sample. No simple imputation is used.
#   Structural low-coverage covariates are explicitly reported in QC files.
#
# Lag construction:
#   Only exact observed year pairs are used. No interpolation is performed.
#   10-year lag is primary; 5-year and 15-year lags are sensitivity analyses.
#
# Inference:
#   Country-clustered standard errors; BH-FDR is applied within each model,
#   lag, and outcome family across the four dietary exposures.
#
# Interpretation:
#   This is a country-level ecological association analysis and is not an
#   individual-level causal effect estimate.
# ============================================================================

# ---- 1. Packages, paths, and options ---------------------------------------
required_packages <- c("data.table", "fixest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install required packages before running: ",
       paste(missing_packages, collapse = ", "))
}

library(data.table)
library(fixest)

script_version <- "health_TWFE_v5_ASR60plus_exact_lags_complete_case"
message("Running: ", script_version)

analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
input_file <- file.path(
  analysis_dir, "02_processed", "Master",
  "GDD_final_master_database_185_country_year_v2.csv"
)

output_dir <- file.path(
  analysis_dir, "TWFE_health_dietary_patterns_v5"
)
table_dir <- file.path(output_dir, "tables")
data_dir <- file.path(output_dir, "model_data")
model_dir <- file.path(output_dir, "model_objects")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

study_years <- c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L)
health_lags <- c(0L, 5L, 10L, 15L)
primary_lag <- 10L

if (!file.exists(input_file)) stop("Input file not found: ", input_file)

ssc_setting <- fixest::ssc(
  adj = TRUE,
  cluster.adj = TRUE,
  cluster.df = "min",
  t.df = "min",
  fixef.K = "nested"
)

# ---- 2. Read and validate the final master database -----------------------
dt <- fread(input_file, na.strings = c("", "NA", "NaN"), encoding = "UTF-8")
setnames(dt, names(dt), tolower(names(dt)))
dt[, iso3 := toupper(trimws(as.character(iso3)))]
dt[, year := as.integer(year)]

diet_vars <- c(
  "diet_mphdi_score", "diet_med_score", "diet_dash_score", "diet_ahei_score"
)
diet_labels <- c(
  diet_mphdi_score = "mPHDI",
  diet_med_score = "MED",
  diet_dash_score = "DASH",
  diet_ahei_score = "AHEI"
)
diet_increments <- c(
  diet_mphdi_score = 10,
  diet_med_score = 1,
  diet_dash_score = 5,
  diet_ahei_score = 10
)

health_outcomes <- data.table(
  outcome_code = c("daly", "prevalence"),
  outcome_variable = c(
    "dementia_daly_asr60plus", "dementia_prevalence_asr60plus"
  ),
  outcome_label = c(
    "Dementia DALY rate",
    "Dementia prevalence rate"
  ),
  primary_outcome = c(TRUE, FALSE)
)

covariates <- c(
  "sdi", "mean_years_schooling", "smoking_prevalence", "adult_obesity",
  "hypertension", "diabetes_prevalence_60plus",
  "insufficient_physical_activity", "pm25", "uhc_service_coverage",
  "competitive_mortality_rate_60plus"
)

# The final master v2 may not retain the population field. It is not needed
# for the unweighted primary health models. Identify an available equivalent
# only for the optional population-weighted sensitivity analysis.
population_candidates <- c(
  "population_60plus_from_gbd", "population_60plus_full",
  "population_60plus_represented", "population_60plus"
)
population_candidates <- population_candidates[population_candidates %in% names(dt)]
has_population_weight <- length(population_candidates) > 0L
population_weight_source <- if (has_population_weight) population_candidates[1L] else NA_character_

required_columns <- c(
  "iso3", "year", "superregion2", diet_vars,
  health_outcomes$outcome_variable, covariates
)
missing_columns <- setdiff(required_columns, names(dt))
if (length(missing_columns) > 0L) {
  stop("Final master is missing: ", paste(missing_columns, collapse = ", "))
}

# CSV imports can represent numeric fields as character strings when a small
# number of missing/non-standard entries is present. Convert only variables
# used by the models and preserve identifiers/labels as character fields.
numeric_model_columns <- c(
  diet_vars, health_outcomes$outcome_variable, covariates
)
for (v in numeric_model_columns) {
  dt[, (v) := suppressWarnings(as.numeric(get(v)))]
}

duplicate_keys <- dt[, .N, by = .(iso3, year)][N > 1L]
if (nrow(duplicate_keys) > 0L) {
  fwrite(duplicate_keys, file.path(table_dir, "QC_duplicate_country_year.csv"))
  stop("Duplicate iso3-year records detected.")
}

if (uniqueN(dt$iso3) != 185L || nrow(dt) != 1295L) {
  stop("Expected 185 countries and 1,295 country-year rows; observed ",
       uniqueN(dt$iso3), " countries and ", nrow(dt), " rows.")
}
if (!identical(sort(unique(dt$year)), study_years)) {
  stop("Unexpected study years: ", paste(sort(unique(dt$year)), collapse = ", "))
}

if (any(vapply(health_outcomes$outcome_variable, function(v) {
  any(dt[[v]] <= 0, na.rm = TRUE)
}, logical(1)))) {
  stop("Non-positive dementia outcome values cannot be log-transformed.")
}
dt[, `:=`(
  ln_dementia_daly = log(dementia_daly_asr60plus),
  ln_dementia_prevalence = log(dementia_prevalence_asr60plus)
)]

# ---- 3. Missingness and data-coverage QC ----------------------------------
missing_qc <- rbindlist(lapply(c(diet_vars, health_outcomes$outcome_variable,
                                 covariates), function(v) {
  data.table(
    variable = v,
    n_rows = nrow(dt),
    n_missing = sum(is.na(dt[[v]])),
    missing_percent = 100 * mean(is.na(dt[[v]])),
    n_countries_with_data = uniqueN(dt[!is.na(get(v)), iso3]),
    observed_years = paste(sort(unique(dt[!is.na(get(v)), year])), collapse = ", ")
  )
}))
fwrite(missing_qc, file.path(table_dir, "QC_health_covariate_missingness_v1.csv"), bom = TRUE)

coverage_by_year <- rbindlist(lapply(c(diet_vars, health_outcomes$outcome_variable,
                                       covariates), function(v) {
  dt[, .(
    variable = v,
    n_rows = .N,
    n_observed = sum(!is.na(get(v))),
    n_missing = sum(is.na(get(v))),
    observed_percent = 100 * mean(!is.na(get(v)))
  ), by = year]
}))
fwrite(coverage_by_year, file.path(table_dir, "QC_health_covariate_coverage_by_year_v1.csv"), bom = TRUE)

# ---- 4. Exact lag data construction ---------------------------------------
make_exact_lag_data <- function(data, lag_years) {
  exposure <- data[, c(
    list(iso3 = iso3, exposure_year = year),
    setNames(.SD, paste0("x_", names(.SD)))
  ), .SDcols = c(diet_vars, covariates)]
  exposure[, outcome_year := exposure_year + lag_years]

  outcome <- data[, .(
    iso3,
    outcome_year = year,
    superregion2,
    dementia_daly_asr60plus,
    dementia_prevalence_asr60plus
  )]
  if (has_population_weight) {
    outcome[, population_60plus_weight := data[[population_weight_source]]]
  }
  out <- merge(outcome, exposure, by = c("iso3", "outcome_year"), all = FALSE)
  out[, lag_years := lag_years]
  out[, `:=`(
    ln_dementia_daly = log(dementia_daly_asr60plus),
    ln_dementia_prevalence = log(dementia_prevalence_asr60plus)
  )]
  out[]
}

lag_data <- setNames(
  lapply(health_lags, function(z) make_exact_lag_data(dt, z)),
  paste0("lag", health_lags)
)

lag_coverage <- rbindlist(lapply(lag_data, function(x) {
  x[, .(
    lag_years = first(lag_years),
    country_year_rows = .N,
    countries = uniqueN(iso3),
    exposure_years = paste(sort(unique(exposure_year)), collapse = ", "),
    outcome_years = paste(sort(unique(outcome_year)), collapse = ", "),
    exact_year_pairs = paste(
      sort(unique(paste(exposure_year, outcome_year, sep = " -> "))),
      collapse = "; "
    )
  )]
}))
setorder(lag_coverage, lag_years)
fwrite(lag_coverage, file.path(table_dir, "QC_health_exact_lag_coverage_v1.csv"), bom = TRUE)

for (z in health_lags) {
  fwrite(lag_data[[paste0("lag", z)]],
         file.path(data_dir, paste0("health_exact_", z, "year_lag_data_v1.csv")),
         bom = TRUE)
}

# ---- 5. Model helpers ------------------------------------------------------
model_specs <- list(
  M0 = list(label = "Model 0", vars = character(0), role = "Basic FE model"),
  M1 = list(label = "Model 1", vars = c("x_sdi"), role = "Primary prespecified model"),
  M2 = list(label = "Model 2", vars = c("x_sdi", "x_mean_years_schooling"), role = "Education-adjusted extension"),
  M3 = list(label = "Model 3", vars = c(
    "x_sdi", "x_mean_years_schooling", "x_smoking_prevalence",
    "x_adult_obesity", "x_hypertension", "x_diabetes_prevalence_60plus",
    "x_insufficient_physical_activity"
  ), role = "Lifestyle and metabolic extension"),
  M4 = list(label = "Model 4", vars = c(
    "x_sdi", "x_mean_years_schooling", "x_smoking_prevalence",
    "x_adult_obesity", "x_hypertension", "x_diabetes_prevalence_60plus",
    "x_insufficient_physical_activity", "x_pm25", "x_uhc_service_coverage"
  ), role = "Environment and health-system extension"),
  M5 = list(vars = c(
    "x_sdi", "x_mean_years_schooling", "x_smoking_prevalence",
    "x_adult_obesity", "x_hypertension", "x_diabetes_prevalence_60plus",
    "x_insufficient_physical_activity", "x_pm25", "x_uhc_service_coverage",
    "x_competitive_mortality_rate_60plus"
  ), label = "Model 5", role = "Competing-mortality extension")
)

exposure_terms <- c(
  diet_mphdi_score = "x_diet_mphdi_score_per_increment",
  diet_med_score = "x_diet_med_score_per_increment",
  diet_dash_score = "x_diet_dash_score_per_increment",
  diet_ahei_score = "x_diet_ahei_score_per_increment"
)

prepare_model_data <- function(x, spec) {
  vars <- c("iso3", "outcome_year",
            "ln_dementia_daly", "ln_dementia_prevalence",
            unname(exposure_terms), spec$vars,
            if (has_population_weight) "population_60plus_weight")
  y <- copy(x[, ..vars])
  y[, complete_case := complete.cases(y)]
  y[complete_case == TRUE]
}

fit_one <- function(data, outcome_variable, exposure_term,
                    covariate_terms = character(0), weighted = FALSE) {
  rhs <- paste(c(exposure_term, covariate_terms), collapse = " + ")
  fml <- as.formula(paste0(outcome_variable, " ~ ", rhs,
                           " | iso3 + outcome_year"))
  if (weighted) {
    feols(fml, data = data, weights = ~population_60plus_weight,
          vcov = ~iso3, ssc = ssc_setting, notes = FALSE)
  } else {
    feols(fml, data = data, vcov = ~iso3,
          ssc = ssc_setting, notes = FALSE)
  }
}

extract_result <- function(model, metadata, term) {
  ct <- coeftable(model)
  if (!term %in% rownames(ct)) stop("Exposure term not found: ", term)
  ci <- confint(model, parm = term, level = 0.95)
  fs <- tryCatch(fitstat(model, "wr2")[[1]], error = function(e) NA_real_)
  beta_hat <- as.numeric(ct[term, "Estimate"])
  se_hat <- as.numeric(ct[term, "Std. Error"])
  p_hat <- as.numeric(ct[term, "Pr(>|t|)"])
  beta_low <- as.numeric(ci[1, 1])
  beta_high <- as.numeric(ci[1, 2])
  out <- copy(metadata)
  out[, `:=`(
    beta = beta_hat,
    standard_error = se_hat,
    p_value = p_hat,
    beta_lower_95 = beta_low,
    beta_upper_95 = beta_high,
    percent_difference = 100 * (exp(beta_hat) - 1),
    percent_difference_lower_95 = 100 * (exp(beta_low) - 1),
    percent_difference_upper_95 = 100 * (exp(beta_high) - 1),
    observations = nobs(model),
    countries = tryCatch(length(fixef(model)$iso3), error = function(e) NA_integer_),
    within_r_squared = as.numeric(fs)
  )]
  out[]
}

# ---- 6. Primary and lagged health models ----------------------------------
results <- list()
models <- list()
model_qc <- list()
counter <- 0L

run_specification <- function(x, lag_value, spec_id, spec, weighted = FALSE,
                              analysis_type = "Health model") {
  prepared <- prepare_model_data(x, spec)
  current_covariates <<- spec$vars

  initial_qc_index <- length(model_qc) + 1L
  model_qc[[initial_qc_index]] <<- data.table(
    model_spec = spec_id,
    model_label = spec$label,
    lag_years = lag_value,
    analysis_type = analysis_type,
    n_available_before_complete_case = nrow(x),
    n_complete_case = nrow(prepared),
    n_countries = uniqueN(prepared$iso3),
    n_outcome_years = uniqueN(prepared$outcome_year),
    complete_panel = nrow(prepared) == 185L * uniqueN(prepared$outcome_year),
    adjustment_variables = if (length(spec$vars)) paste(spec$vars, collapse = "; ") else "None",
    fit_status = if (nrow(prepared) == 0L) "SKIPPED_NO_COMPLETE_CASE"
                 else if (uniqueN(prepared$iso3) < 10L) "SKIPPED_TOO_FEW_COUNTRIES"
                 else if (uniqueN(prepared$outcome_year) < 2L) "SKIPPED_TOO_FEW_YEARS"
                 else "PENDING"
  )

  if (nrow(prepared) == 0L || uniqueN(prepared$iso3) < 10L ||
      uniqueN(prepared$outcome_year) < 2L) return(invisible(NULL))

  for (d in diet_vars) {
    term <- exposure_terms[[d]]
    for (j in seq_len(nrow(health_outcomes))) {
      outcome_code <- health_outcomes$outcome_code[j]
      outcome_var <- c(
        daly = "ln_dementia_daly",
        prevalence = "ln_dementia_prevalence"
      )[[outcome_code]]
      model_id <- paste("H", outcome_code, d, spec_id,
                        paste0("lag", lag_value), sep = "_")
      fit <- tryCatch(
        fit_one(prepared, outcome_var, term, covariate_terms = spec$vars,
                weighted = weighted),
        error = function(e) e
      )
      if (inherits(fit, "error")) {
        model_qc[[initial_qc_index]][, fit_status := "ERROR"]
        model_qc[[length(model_qc) + 1L]] <<- data.table(
          model_spec = spec_id, model_label = spec$label,
          lag_years = lag_value, analysis_type = analysis_type,
          n_available_before_complete_case = nrow(x), n_complete_case = nrow(prepared),
          n_countries = uniqueN(prepared$iso3), n_outcome_years = uniqueN(prepared$outcome_year),
          complete_panel = FALSE, adjustment_variables = paste(spec$vars, collapse = "; "),
          fit_status = "ERROR", fit_error = conditionMessage(fit)
        )
        next
      }
      model_qc[[initial_qc_index]][, fit_status := "SUCCESS"]
      counter <<- counter + 1L
      models[[model_id]] <<- fit
      saveRDS(fit, file.path(model_dir, paste0(model_id, ".rds")))
      metadata <- data.table(
        model_id = model_id,
        model_spec = spec_id,
        model_label = spec$label,
        model_role = spec$role,
        analysis_type = analysis_type,
        diet_variable = d,
        diet_label = unname(diet_labels[[d]]),
        score_increment = unname(diet_increments[[d]]),
        outcome_code = outcome_code,
        outcome_label = health_outcomes$outcome_label[j],
        lag_years = lag_value,
        exposure_years = paste(sort(unique(x$exposure_year)), collapse = ", "),
        outcome_years = paste(sort(unique(x$outcome_year)), collapse = ", "),
        weighting = if (weighted) "Population aged 60+" else "Unweighted",
        fixed_effects = "Country and outcome year",
        clustered_by = "Country",
        complete_case = TRUE,
        primary_model = spec_id == "M1" && lag_value == primary_lag && !weighted,
        exposure_term = term
      )
      results[[counter]] <<- extract_result(fit, metadata, term)
    }
  }
  invisible(NULL)
}

for (nm in names(lag_data)) {
  lag_data[[nm]][, `:=`(
    x_diet_mphdi_score_per_increment = x_diet_mphdi_score / diet_increments[["diet_mphdi_score"]],
    x_diet_med_score_per_increment = x_diet_med_score / diet_increments[["diet_med_score"]],
    x_diet_dash_score_per_increment = x_diet_dash_score / diet_increments[["diet_dash_score"]],
    x_diet_ahei_score_per_increment = x_diet_ahei_score / diet_increments[["diet_ahei_score"]]
  )]
}

for (z in health_lags[health_lags %in% c(0L, 5L, 10L, 15L)]) {
  analysis_type <- if (z == primary_lag) "Primary 10-year lag" else "Lag sensitivity"
  run_specification(lag_data[[paste0("lag", z)]], z, "M0", model_specs$M0,
                    weighted = FALSE, analysis_type = analysis_type)
  run_specification(lag_data[[paste0("lag", z)]], z, "M1", model_specs$M1,
                    weighted = FALSE, analysis_type = analysis_type)
}
for (m in c("M2", "M3", "M4", "M5")) {
  run_specification(lag_data[[paste0("lag", primary_lag)]], primary_lag,
                    m, model_specs[[m]], weighted = FALSE,
                    analysis_type = "10-year covariate extension")
}

if (has_population_weight) {
  run_specification(lag_data[[paste0("lag", primary_lag)]], primary_lag, "M1",
                    model_specs$M1, weighted = TRUE,
                    analysis_type = "10-year population-weighted sensitivity")
} else {
  message("Population-weighted sensitivity skipped: no 60+ population column in the final master v2.")
}

# ---- 7. BH-FDR correction and output tables -------------------------------
if (!length(results)) {
  model_qc_dt <- rbindlist(model_qc, fill = TRUE)
  fwrite(model_qc_dt,
         file.path(table_dir, "QC_health_model_samples_and_fit_status_v1.csv"),
         bom = TRUE)
  error_text <- model_qc_dt[
    fit_status == "ERROR" & !is.na(fit_error) & nzchar(fit_error),
    fit_error
  ][1L]
  stop("No health models were successfully fitted. First recorded fit error: ",
       ifelse(is.na(error_text), "none; inspect the QC file for skipped models", error_text))
}
results_dt <- rbindlist(results, fill = TRUE)

results_dt[, q_value_bh := p.adjust(p_value, method = "BH"),
           by = .(model_spec, lag_years, outcome_code, analysis_type, weighting)]
results_dt[, fdr_family := paste(model_spec, paste0("lag", lag_years),
                                 outcome_code, weighting, sep = "_")]

fwrite(results_dt, file.path(table_dir, "TABLE_health_TWFE_all_results_v1.csv"), bom = TRUE)

primary_results <- results_dt[
  primary_model == TRUE & lag_years == primary_lag & weighting == "Unweighted"
]
fwrite(primary_results, file.path(table_dir, "TABLE_health_primary_10year_lag_Model1_v1.csv"), bom = TRUE)

model_comparison <- results_dt[
  lag_years == primary_lag & weighting == "Unweighted" &
    model_spec %in% names(model_specs)
]
fwrite(model_comparison, file.path(table_dir, "TABLE_health_Model0_to_Model5_10year_v1.csv"), bom = TRUE)

lag_sensitivity <- results_dt[
  model_spec == "M1" & weighting == "Unweighted" & lag_years %in% c(5L, 15L)
]
fwrite(lag_sensitivity, file.path(table_dir, "TABLE_health_5year_15year_lag_sensitivity_v1.csv"), bom = TRUE)

population_weighted <- results_dt[
  analysis_type == "10-year population-weighted sensitivity"
]
fwrite(population_weighted, file.path(table_dir, "TABLE_health_population_weighted_sensitivity_v1.csv"), bom = TRUE)

fwrite(rbindlist(model_qc, fill = TRUE),
       file.path(table_dir, "QC_health_model_samples_and_fit_status_v1.csv"), bom = TRUE)

manifest <- data.table(
  script_version = script_version,
  input_file = input_file,
  analysis_type = "Dietary-pattern scores and dementia burden",
  countries_expected = 185L,
  rows_expected = 1295L,
  study_years = paste(study_years, collapse = ", "),
  primary_lag_years = primary_lag,
  sensitivity_lags = paste(c(5L, 15L), collapse = ", "),
  primary_dementia_outcome = "dementia_daly_asr60plus",
  complementary_outcome = "dementia_prevalence_asr60plus",
  fixed_effects = "iso3 + outcome_year",
  standard_errors = "Country-clustered",
  missing_data = "Model-specific complete-case analysis; no simple imputation",
  fdr_method = "Benjamini-Hochberg within model-lag-outcome-weighting family",
  interpolation_used = FALSE,
  population_weight_source = population_weight_source,
  environment_models_included = FALSE
)
fwrite(manifest, file.path(table_dir, "MANIFEST_health_TWFE_analysis_v1.csv"), bom = TRUE)

message("\nCompleted: dietary-pattern and dementia health TWFE analyses")
message("Countries: ", uniqueN(dt$iso3))
message("Rows: ", format(nrow(dt), big.mark = ","))
message("Primary lag: ", primary_lag, " years; exact observed year pairs only")
message("Successful fitted models: ", length(models))
message("Primary results: ", file.path(table_dir, "TABLE_health_primary_10year_lag_Model1_v1.csv"))
message("Model comparison: ", file.path(table_dir, "TABLE_health_Model0_to_Model5_10year_v1.csv"))
message("Lag sensitivity: ", file.path(table_dir, "TABLE_health_5year_15year_lag_sensitivity_v1.csv"))
message("QC: ", file.path(table_dir, "QC_health_model_samples_and_fit_status_v1.csv"))
message("No environmental outcome models were fitted in this script.")

