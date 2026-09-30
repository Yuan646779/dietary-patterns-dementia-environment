# Same-sample comparison: refit Models 0-2 on the Model 3-5 analytic sample.
# Standalone add-on to health_TWFE_v5_ASR60plus_exact_lags_complete_case.
# Only the exact observed 10-year exposure lag is used. No imputation.
# Existing primary-model estimates and files are not changed.

required_packages <- c("data.table", "fixest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}
library(data.table)
library(fixest)

# Match the paths of the original health-analysis script. Edit if needed.
analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
input_file <- file.path(
  analysis_dir, "02_processed", "Master",
  "GDD_final_master_database_185_country_year_v2.csv"
)
output_dir <- file.path(
  analysis_dir, "TWFE_health_dietary_patterns_v5",
  "same_sample_comparison_CC306"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

diet_vars <- c("diet_mphdi_score", "diet_med_score", "diet_dash_score",
               "diet_ahei_score")
diet_labels <- c("mPHDI", "MED", "DASH", "AHEI")
increments <- c(10, 1, 5, 10)
outcome_vars <- c("dementia_daly_asr60plus",
                  "dementia_prevalence_asr60plus")
outcome_codes <- c("daly", "prevalence")

# Model 5 includes every time-varying covariate from Models 1-5.
# Consequently, its joint complete cases define the common sample.
covariates_m5 <- c(
  "sdi", "mean_years_schooling", "smoking_prevalence", "adult_obesity",
  "hypertension", "diabetes_prevalence_60plus",
  "insufficient_physical_activity", "pm25", "uhc_service_coverage",
  "competitive_mortality_rate_60plus"
)

dt <- fread(input_file, na.strings = c("", "NA", "NaN"), encoding = "UTF-8")
setnames(dt, names(dt), tolower(names(dt)))
required <- c("iso3", "year", diet_vars, outcome_vars, covariates_m5)
missing_columns <- setdiff(required, names(dt))
if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}
dt[, iso3 := toupper(trimws(as.character(iso3)))]
dt[, year := as.integer(year)]
population_candidates <- c(
  "population_60plus_from_gbd", "population_60plus_full",
  "population_60plus_represented", "population_60plus"
)
population_candidates <- population_candidates[
  population_candidates %in% names(dt)
]
population_weight_source <- if (length(population_candidates)) {
  population_candidates[1L]
} else {
  NA_character_
}
for (v in c(diet_vars, outcome_vars, covariates_m5)) {
  dt[, (v) := suppressWarnings(as.numeric(get(v)))]
}
if (dt[, anyDuplicated(paste(iso3, year))] > 0L) {
  stop("Duplicate country-year records in the master dataset.")
}
if (nrow(dt) != 1295L || uniqueN(dt$iso3) != 185L ||
    !identical(sort(unique(dt$year)),
               c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L))) {
  stop("Master dataset does not match the original 185-country, 7-year panel.")
}

# Mirror the original exact-lag merge: exposures/covariates at exposure year,
# outcomes at outcome year. In particular, the 2018 outcome is not included.
exposure <- dt[, c(
  list(iso3 = iso3, exposure_year = year),
  setNames(.SD, paste0("x_", names(.SD)))
), .SDcols = c(diet_vars, covariates_m5)]
exposure[, outcome_year := exposure_year + 10L]
outcome <- dt[, .(
  iso3, outcome_year = year,
  dementia_daly_asr60plus, dementia_prevalence_asr60plus
)]
if (!is.na(population_weight_source)) {
  # The original unweighted prepare_model_data also checks this field when
  # present. Keep its complete-case sample identical to that of Model 5.
  outcome[, population_60plus_weight := dt[[population_weight_source]]]
}
lag10 <- merge(outcome, exposure, by = c("iso3", "outcome_year"),
               all = FALSE)
if (nrow(lag10) != 740L) {
  stop("Expected 740 exact 10-year exposure-outcome pairs; found ",
       nrow(lag10), ".")
}
if (any(lag10[, c(outcome_vars), with = FALSE] <= 0, na.rm = TRUE)) {
  stop("Dementia rates must be positive before log transformation.")
}
lag10[, `:=`(
  ln_dementia_daly = log(dementia_daly_asr60plus),
  ln_dementia_prevalence = log(dementia_prevalence_asr60plus)
)]
for (i in seq_along(diet_vars)) {
  lag10[, (paste0("x_", diet_vars[i], "_per_increment")) :=
          get(paste0("x_", diet_vars[i])) / increments[i]]
}

exposure_terms <- paste0("x_", diet_vars, "_per_increment")
all_model_terms <- c(
  "iso3", "outcome_year", "ln_dementia_daly", "ln_dementia_prevalence",
  exposure_terms, paste0("x_", covariates_m5),
  if (!is.na(population_weight_source)) "population_60plus_weight"
)
# This reproduces the unweighted prepare_model_data(lag10, Model 5) rule in
# the original script: both outcomes, four dietary scores, and all M5
# covariates. It yields 312 complete records from 159 countries.
common_complete <- copy(
  lag10[complete.cases(lag10[, ..all_model_terms])]
)
if (nrow(common_complete) != 312L ||
    uniqueN(common_complete$iso3) != 159L) {
  stop("Expected 312 complete records in 159 countries before fixed-effect ",
       "singleton removal; found ", nrow(common_complete), " rows in ",
       uniqueN(common_complete$iso3), " countries. Check the master data.")
}

# fixest removes countries represented by only one record because their
# observations are fixed-effect singletons. The original Model 3-5 estimates
# therefore use 306 observations from 153 countries. Remove the same six
# singleton records explicitly so Models 0-2 use the identical analytic sample.
country_counts <- common_complete[, .N, by = iso3]
singleton_countries <- country_counts[N == 1L, iso3]
common <- common_complete[!iso3 %in% singleton_countries]
if (length(singleton_countries) != 6L || nrow(common) != 306L ||
    uniqueN(common$iso3) != 153L) {
  stop("Expected removal of 6 singleton countries to leave 306 observations ",
       "in 153 countries; found ", length(singleton_countries),
       " singleton countries and ", nrow(common), " rows in ",
       uniqueN(common$iso3), " countries.")
}
setorder(common, iso3, outcome_year)

fwrite(
  common_complete[, .(
    iso3, exposure_year, outcome_year,
    retained_in_estimation = !iso3 %in% singleton_countries
  )],
  file.path(output_dir,
            "QC_exact10_complete_records_and_singletons_CC306.csv"),
  bom = TRUE
)
fwrite(common[, .(iso3, exposure_year, outcome_year)],
       file.path(output_dir, "QC_exact10_common_analytic_sample_keys_CC306.csv"),
       bom = TRUE)

ssc_setting <- fixest::ssc(
  adj = TRUE, cluster.adj = TRUE, cluster.df = "min", t.df = "min",
  fixef.K = "nested"
)
model_covariates <- list(
  M0 = character(0),
  M1 = "x_sdi",
  M2 = c("x_sdi", "x_mean_years_schooling")
)

results <- vector("list", length(model_covariates) * length(diet_vars) *
                         length(outcome_codes))
k <- 0L
for (model_id in names(model_covariates)) {
  for (i in seq_along(diet_vars)) {
    for (j in seq_along(outcome_codes)) {
      term <- exposure_terms[i]
      outcome_term <- c("ln_dementia_daly", "ln_dementia_prevalence")[j]
      rhs <- paste(c(term, model_covariates[[model_id]]), collapse = " + ")
      fml <- as.formula(paste0(outcome_term, " ~ ", rhs,
                               " | iso3 + outcome_year"))
      fit <- feols(fml, data = common, vcov = ~iso3,
                   ssc = ssc_setting, notes = FALSE)
      if (nobs(fit) != nrow(common)) {
        stop(model_id, "/", diet_labels[i], "/", outcome_codes[j],
             " dropped observations during estimation: ", nobs(fit),
             " instead of ", nrow(common), ".")
      }
      ct <- coeftable(fit)
      ci <- confint(fit, parm = term, level = 0.95)
      beta <- unname(ct[term, "Estimate"])
      low <- unname(ci[1, 1])
      high <- unname(ci[1, 2])
      k <- k + 1L
      results[[k]] <- data.table(
        analysis = "Same analytic sample as Models 3-5",
        model = paste0(model_id, "_CC306"),
        diet_pattern = diet_labels[i], score_increment = increments[i],
        outcome = outcome_codes[j], lag_years = 10L,
        beta = beta, standard_error = unname(ct[term, "Std. Error"]),
        percent_difference = 100 * (exp(beta) - 1),
        lower_95 = 100 * (exp(low) - 1),
        upper_95 = 100 * (exp(high) - 1),
        p_value = unname(ct[term, "Pr(>|t|)"]),
        country_years = nobs(fit),
        countries = uniqueN(common$iso3),
        exposure_years = paste(sort(unique(common$exposure_year)),
                               collapse = ", "),
        outcome_years = paste(sort(unique(common$outcome_year)),
                              collapse = ", ")
      )
    }
  }
}

result_dt <- rbindlist(results)
# Same FDR family as the original: four dietary patterns within each model
# and outcome (all models here have the same exact 10-year lag and weighting).
result_dt[, q_value_bh := p.adjust(p_value, method = "BH"),
          by = .(model, outcome)]
setorder(result_dt, outcome, diet_pattern, model)
fwrite(result_dt,
       file.path(output_dir,
                 "TABLE_health_Model0to2_exact10_same_sample_CC306.csv"),
       bom = TRUE)
message("Completed same-sample comparison: ", nrow(common),
        " country-years; ", uniqueN(common$iso3), " countries; ",
        nrow(result_dt), " estimates (Models 0-2 only).")
print(result_dt[diet_pattern == "mPHDI" & outcome == "daly",
                .(model, percent_difference, lower_95, upper_95,
                  p_value, q_value_bh, country_years)])
