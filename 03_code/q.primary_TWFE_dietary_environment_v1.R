# ============================================================================
# q.primary_TWFE_dietary_environment_v1.R
# Concurrent TWFE associations between dietary-pattern scores and environmental
# burdens in the final 185-country x 7-year master panel.
# Environmental values cover the 15 parameterized GDD food groups, not total diet.
# ============================================================================

required_packages <- c("data.table", "fixest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required packages before running: ",
       paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

script_version <- "environment_TWFE_v1_concurrent_complete_case"
analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
input_file <- file.path(
  analysis_dir, "02_processed", "Master",
  "GDD_final_master_database_185_country_year_v2.csv"
)
population_file <- file.path(
  analysis_dir, "01_raw", "GBD",
  "GBD_age_specific_population_60plus.csv"
)
output_dir <- file.path(analysis_dir, "03_analysis", "Environment_TWFE_v1")
table_dir <- file.path(output_dir, "tables")
model_dir <- file.path(output_dir, "model_objects")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

study_years <- c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L)
expected_rows <- 1295L
expected_countries <- 185L

normalise_names <- function(x) {
  x <- sub("^\ufeff", "", x)
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}
as_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
}

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
if (!file.exists(population_file)) stop("Population file not found: ", population_file)

ssc_setting <- fixest::ssc(
  adj = TRUE, cluster.adj = TRUE, cluster.df = "min",
  t.df = "min", fixef.K = "nested"
)

# ---- 1. Read and validate the final master -------------------------------
dt <- fread(input_file, na.strings = c("", "NA", "NaN", ".", "NULL"),
            encoding = "UTF-8")
setnames(dt, names(dt), normalise_names(names(dt)))
if (!all(c("iso3", "year") %in% names(dt))) {
  stop("Final master must contain iso3 and year.")
}
dt[, c("iso3", "year") := list(
  toupper(trimws(as.character(iso3))), as.integer(year)
)]

diet_vars <- c(
  "diet_mphdi_score", "diet_med_score",
  "diet_dash_score", "diet_ahei_score"
)
diet_labels <- c(
  diet_mphdi_score = "mPHDI", diet_med_score = "MED",
  diet_dash_score = "DASH", diet_ahei_score = "AHEI"
)
diet_increments <- c(
  diet_mphdi_score = 10, diet_med_score = 1,
  diet_dash_score = 5, diet_ahei_score = 10
)
exposure_terms <- c(
  diet_mphdi_score = "x_diet_mphdi_score",
  diet_med_score = "x_diet_med_score",
  diet_dash_score = "x_diet_dash_score",
  diet_ahei_score = "x_diet_ahei_score"
)

environment_vars <- c(
  "env_ghg_kg_co2e_person_day",
  "env_land_m2_person_day",
  "env_water_scarcity_l_equivalent_person_day",
  "env_eutrophication_g_po4e_person_day"
)
environment_info <- data.table(
  outcome_code = c("ghg", "land", "water", "eutrophication"),
  outcome_term = c("ln_ghg", "ln_land", "ln_water", "ln_eutrophication"),
  raw_variable = environment_vars,
  outcome_label = c(
    "GHG emissions", "Land use", "Scarcity-weighted water use",
    "Eutrophication potential"
  ),
  unit = c(
    "kg CO2-eq/person/day", "m2/person/day", "L-eq/person/day",
    "g PO4-eq/person/day"
  )
)

required <- c(
  "iso3", "year", "superregion2", diet_vars, environment_vars,
  "sdi", "env_coverage_proportion"
)
missing_columns <- setdiff(required, names(dt))
if (length(missing_columns)) {
  stop("Final master is missing: ", paste(missing_columns, collapse = ", "))
}
if (nrow(dt) != expected_rows || uniqueN(dt$iso3) != expected_countries ||
    !identical(sort(unique(dt$year)), study_years)) {
  stop("Expected 185 countries x 7 years (1,295 rows); observed ",
       uniqueN(dt$iso3), " countries and ", nrow(dt), " rows.")
}
if (anyDuplicated(dt[, .(iso3, year)])) {
  stop("Duplicate iso3-year keys detected.")
}

for (v in c(diet_vars, environment_vars, "sdi", "env_coverage_proportion")) {
  dt[, (v) := as_numeric_safe(get(v))]
}

# Restore population aged 60+ when it is absent from master v2.
population_candidates <- c(
  "population_60plus_from_gbd", "population_60plus",
  "population_60_plus", "population_60plus_full"
)
population_column <- intersect(population_candidates, names(dt))[1L]
if (is.na(population_column)) {
  pop <- fread(population_file, na.strings = c("", "NA", "NaN", ".", "NULL"),
               encoding = "UTF-8")
  setnames(pop, names(pop), normalise_names(names(pop)))
  pop_value_candidates <- c("population", "val", "value", "estimate")
  pop_value_column <- intersect(pop_value_candidates, names(pop))[1L]
  if (!all(c("iso3", "year") %in% names(pop)) ||
      is.na(pop_value_column)) {
    stop("Cannot identify iso3, year, and population columns in: ",
         population_file)
  }
  pop[, c("iso3", "year", "population_value") := list(
    toupper(trimws(as.character(iso3))), as.integer(year),
    as_numeric_safe(get(pop_value_column))
  )]
  pop <- pop[iso3 %in% unique(dt$iso3) & year %in% study_years]
  pop <- pop[, .(population_60plus_from_gbd =
                   sum(population_value, na.rm = TRUE)), by = .(iso3, year)]
  dt <- merge(dt, pop, by = c("iso3", "year"), all.x = TRUE)
  population_source <- basename(population_file)
} else {
  if (population_column != "population_60plus_from_gbd") {
    setnames(dt, population_column, "population_60plus_from_gbd")
  }
  population_source <- "Final master population column"
}
if (anyNA(dt$population_60plus_from_gbd) ||
    any(!is.finite(dt$population_60plus_from_gbd)) ||
    any(dt$population_60plus_from_gbd <= 0)) {
  stop("Invalid population aged 60+ values.")
}

if (any(!is.na(dt$env_coverage_proportion) &
        (dt$env_coverage_proportion < 0 |
         dt$env_coverage_proportion > 1))) {
  stop("env_coverage_proportion must be between 0 and 1.")
}
dt[, c("coverage_per_10pct", "coverage_ge80") := list(
  env_coverage_proportion / 0.10,
  !is.na(env_coverage_proportion) & env_coverage_proportion >= 0.80
)]
for (v in environment_vars) {
  if (any(dt[[v]] <= 0, na.rm = TRUE)) {
    stop("Non-positive environmental values cannot be log-transformed: ", v)
  }
}
dt[, c("ln_ghg", "ln_land", "ln_water", "ln_eutrophication",
       "x_diet_mphdi_score", "x_diet_med_score",
       "x_diet_dash_score", "x_diet_ahei_score") := list(
  log(env_ghg_kg_co2e_person_day),
  log(env_land_m2_person_day),
  log(env_water_scarcity_l_equivalent_person_day),
  log(env_eutrophication_g_po4e_person_day),
  diet_mphdi_score / diet_increments[["diet_mphdi_score"]],
  diet_med_score / diet_increments[["diet_med_score"]],
  diet_dash_score / diet_increments[["diet_dash_score"]],
  diet_ahei_score / diet_increments[["diet_ahei_score"]]
)]

# ---- 2. QC -----------------------------------------------------------------
qc_vars <- c(diet_vars, environment_vars, "sdi",
             "env_coverage_proportion", "population_60plus_from_gbd")
missingness_qc <- rbindlist(lapply(qc_vars, function(v) data.table(
  variable = v, n_rows = nrow(dt), n_missing = sum(is.na(dt[[v]])),
  missing_percent = 100 * mean(is.na(dt[[v]])),
  observed_years = paste(sort(unique(dt[!is.na(get(v)), year])), collapse = ", ")
)), fill = TRUE)
coverage_qc <- dt[, .(
  n_rows = .N, n_observed = sum(!is.na(env_coverage_proportion)),
  observed_percent = 100 * mean(!is.na(env_coverage_proportion)),
  n_coverage_ge80 = sum(coverage_ge80, na.rm = TRUE)
), by = year]
fwrite(missingness_qc,
       file.path(table_dir, "QC_environment_missingness_v1.csv"), bom = TRUE)
fwrite(coverage_qc,
       file.path(table_dir, "QC_environment_coverage_by_year_v1.csv"), bom = TRUE)

# ---- 3. TWFE helpers ------------------------------------------------------
fit_twfe <- function(data, outcome_term, exposure_term,
                     adjustment_terms = character(0), weighted = FALSE) {
  rhs <- paste(c(exposure_term, adjustment_terms), collapse = " + ")
  fml <- as.formula(paste0(outcome_term, " ~ ", rhs, " | iso3 + year"))
  if (weighted) {
    feols(fml, data = data, weights = ~population_60plus_from_gbd,
          vcov = ~iso3, ssc = ssc_setting, notes = FALSE)
  } else {
    feols(fml, data = data, vcov = ~iso3,
          ssc = ssc_setting, notes = FALSE)
  }
}
fitstat_safe <- function(model, stat) {
  tryCatch(as.numeric(fitstat(model, stat)[[1L]]),
           error = function(e) NA_real_)
}
extract_result <- function(model, metadata, term) {
  ct <- coeftable(model)
  if (!term %in% rownames(ct)) stop("Exposure term not found: ", term)
  ci <- confint(model, parm = term, level = 0.95)
  beta_hat <- as.numeric(ct[term, "Estimate"])
  se_hat <- as.numeric(ct[term, "Std. Error"])
  p_hat <- as.numeric(ct[term, "Pr(>|t|)"])
  beta_low <- as.numeric(ci[1L, 1L])
  beta_high <- as.numeric(ci[1L, 2L])
  out <- copy(metadata)
  out[, c("beta", "standard_error", "p_value",
          "beta_lower_95", "beta_upper_95",
          "percent_difference", "percent_difference_lower_95",
          "percent_difference_upper_95", "observations",
          "countries", "within_r_squared") := list(
    beta_hat, se_hat, p_hat, beta_low, beta_high,
    100 * (exp(beta_hat) - 1),
    100 * (exp(beta_low) - 1),
    100 * (exp(beta_high) - 1),
    nobs(model), length(fixef(model)$iso3),
    fitstat_safe(model, "wr2")
  )]
  out[]
}

# Each model uses a specification-specific complete-case sample.
specifications <- list(
  list(id = "unadjusted", label = "Unadjusted",
       filter = function(x) x, adjust = character(0),
       adjust_label = "No time-varying covariates",
       sample_label = "All available country-years",
       weighted = FALSE, primary = FALSE),
  list(id = "primary", label = "Primary analysis",
       filter = function(x) x, adjust = c("sdi", "coverage_per_10pct"),
       adjust_label = "SDI and continuous environmental intake coverage",
       sample_label = "All available country-years",
       weighted = FALSE, primary = TRUE),
  list(id = "coverage80", label = "Coverage >=80%",
       filter = function(x) x[coverage_ge80 == TRUE], adjust = "sdi",
       adjust_label = "SDI",
       sample_label = "Country-years with environmental coverage >=80%",
       weighted = FALSE, primary = FALSE),
  list(id = "coverage80_adjusted",
       label = "Coverage >=80% with additional adjustment for coverage",
       filter = function(x) x[coverage_ge80 == TRUE],
       adjust = c("sdi", "coverage_per_10pct"),
       adjust_label = "SDI and continuous environmental intake coverage",
       sample_label = "Country-years with environmental coverage >=80%",
       weighted = FALSE, primary = FALSE),
  list(id = "population_weighted", label = "Population weighted",
       filter = function(x) x, adjust = c("sdi", "coverage_per_10pct"),
       adjust_label = "SDI and continuous environmental intake coverage",
       sample_label = "All available country-years",
       weighted = TRUE, primary = FALSE)
)

all_models <- list()
all_results <- list()
model_qc <- list()

for (diet_code in names(exposure_terms)) {
  exposure_term <- unname(exposure_terms[[diet_code]])
  diet_label <- unname(diet_labels[[diet_code]])
  exposure_increment <- unname(diet_increments[[diet_code]])
  exposure_unit <- paste0("per ", exposure_increment,
                          "-point higher ", diet_label)

  for (j in seq_len(nrow(environment_info))) {
    outcome <- environment_info[j]
    for (spec in specifications) {
      dat <- spec$filter(dt)
      needed <- c(outcome$outcome_term, exposure_term, spec$adjust)
      if (spec$weighted) needed <- c(needed, "population_60plus_from_gbd")
      dat <- dat[complete.cases(dat[, needed, with = FALSE])]

      model_id <- paste("E", outcome$outcome_code, diet_code, spec$id, sep = "_")
      qc_row <- data.table(
        model_id = model_id, specification = spec$label,
        n_available = nrow(spec$filter(dt)), n_complete_case = nrow(dat),
        countries = uniqueN(dat$iso3), years = uniqueN(dat$year),
        fit_status = "SKIPPED"
      )
      if (nrow(dat) == 0L || uniqueN(dat$iso3) < 10L ||
          uniqueN(dat$year) < 2L) {
        model_qc[[length(model_qc) + 1L]] <- qc_row
        next
      }

      fit <- tryCatch(
        fit_twfe(dat, outcome$outcome_term, exposure_term,
                 adjustment_terms = spec$adjust, weighted = spec$weighted),
        error = function(e) e
      )
      if (inherits(fit, "error")) {
        qc_row[, c("fit_status", "fit_error") :=
                 list("ERROR", conditionMessage(fit))]
        model_qc[[length(model_qc) + 1L]] <- qc_row
        next
      }

      qc_row[, fit_status := "SUCCESS"]
      model_qc[[length(model_qc) + 1L]] <- qc_row
      metadata <- data.table(
        model_id = model_id, family = "Environment",
        diet_code = diet_code, diet_pattern = diet_label,
        exposure_unit = exposure_unit,
        outcome_code = outcome$outcome_code,
        outcome = outcome$outcome_label, unit = outcome$unit,
        specification = spec$label, specification_id = spec$id,
        sample_definition = spec$sample_label,
        adjustment = spec$adjust_label,
        weighting = if (spec$weighted) "Population aged 60+" else "Unweighted",
        fixed_effects = "Country and calendar year",
        clustered_by = "Country", temporal_relation = "Concurrent",
        food_scope = "15 parameterized GDD food groups; not total diet",
        primary_model = spec$primary, term = exposure_term
      )
      all_models[[model_id]] <- fit
      all_results[[length(all_results) + 1L]] <-
        extract_result(fit, metadata, exposure_term)
      saveRDS(fit, file.path(model_dir, paste0(model_id, ".rds")))
    }
  }
}

if (!length(all_results)) {
  qc_file <- file.path(table_dir,
                       "QC_environment_model_samples_and_fit_status_v1.csv")
  fwrite(rbindlist(model_qc, fill = TRUE), qc_file, bom = TRUE)
  stop("No environmental models were successfully fitted. Inspect: ", qc_file)
}

results_dt <- rbindlist(all_results, fill = TRUE)
results_dt[, q_value_bh := p.adjust(p_value, method = "BH"),
           by = .(outcome_code, specification_id, weighting)]
results_dt[, fdr_family := paste(outcome_code, specification_id,
                                 weighting, sep = "_")]

primary_results <- results_dt[primary_model == TRUE]
setorder(primary_results, outcome_code, diet_code)
results_display <- copy(results_dt)
results_display[, estimate_95ci := sprintf(
  "%.2f%% (%.2f%% to %.2f%%)",
  percent_difference, percent_difference_lower_95,
  percent_difference_upper_95
)]
results_display[, p_value_display := format.pval(
  p_value, digits = 3, eps = 0.001
)]
results_display[, q_value_bh_display := format.pval(
  q_value_bh, digits = 3, eps = 0.001
)]

fwrite(results_dt,
       file.path(table_dir, "TABLE_environment_TWFE_all_results_v1.csv"),
       bom = TRUE)
fwrite(primary_results,
       file.path(table_dir, "TABLE_environment_TWFE_primary_models_v1.csv"),
       bom = TRUE)
fwrite(results_display,
       file.path(table_dir, "TABLE_environment_TWFE_all_results_display_v1.csv"),
       bom = TRUE)
fwrite(results_display[specification_id != "primary"],
       file.path(table_dir, "TABLE_environment_TWFE_sensitivity_models_v1.csv"),
       bom = TRUE)

primary_display <- results_display[primary_model == TRUE, .(
  diet_pattern, exposure_unit, outcome, unit, estimate_95ci,
  p_value = p_value_display, q_value_bh = q_value_bh_display,
  observations, countries, adjustment, fixed_effects, clustered_by,
  temporal_relation, food_scope
)]
fwrite(primary_display,
       file.path(table_dir, "TABLE_2_environment_primary_manuscript_display_v1.csv"),
       bom = TRUE)

fwrite(rbindlist(model_qc, fill = TRUE),
       file.path(table_dir, "QC_environment_model_samples_and_fit_status_v1.csv"),
       bom = TRUE)

manifest <- data.table(
  script_version = script_version, input_file = input_file,
  output_dir = output_dir, countries_expected = expected_countries,
  rows_expected = expected_rows,
  study_years = paste(study_years, collapse = ", "),
  exposure_variables = paste(diet_vars, collapse = "; "),
  exposure_increments = "mPHDI 10; MED 1; DASH 5; AHEI 10",
  environmental_outcomes = paste(environment_vars, collapse = "; "),
  temporal_relation = "Concurrent",
  fixed_effects = "Country and calendar year",
  standard_errors = "Country-clustered",
  primary_adjustment = "SDI + continuous environmental intake coverage",
  sensitivity_analyses = paste(
    c("Unadjusted", "Coverage >=80%",
      "Coverage >=80% with additional coverage adjustment",
      "Population weighted"), collapse = "; "
  ),
  missing_data = "Specification-specific complete-case analysis; no simple imputation",
  fdr_method = "BH q-values retained separately; raw P values retained",
  population_weight_source = population_source,
  food_scope = "15 parameterized GDD food groups; not total diet"
)
fwrite(manifest,
       file.path(table_dir, "MANIFEST_environment_TWFE_analysis_v1.csv"),
       bom = TRUE)

message("\nCompleted: dietary-pattern and environmental-burden TWFE analyses")
message("Countries: ", uniqueN(dt$iso3))
message("Rows: ", nrow(dt))
message("Successful environmental models: ", length(all_models))
message("Primary table: ",
        file.path(table_dir,
                  "TABLE_2_environment_primary_manuscript_display_v1.csv"))
message("Sensitivity table: ",
        file.path(table_dir, "TABLE_environment_TWFE_sensitivity_models_v1.csv"))
message("QC: ",
        file.path(table_dir,
                  "QC_environment_model_samples_and_fit_status_v1.csv"))
