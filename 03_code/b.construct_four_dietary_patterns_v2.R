# ============================================================
# b.construct_four_dietary_patterns_v1.R
#
# Purpose:
#   Construct four dietary-pattern scores for adults aged 60 years
#   and older at the GDD country-year level:
#   mPHDI, AHEI, DASH, and MED.
#
# Input:
#   02_processed/GDD/GDD_country_60plus_all47_country_year.csv
#
# Important methodological note:
#   The input from step a is already aggregated across the eight
#   older-adult age groups. Therefore this script does not apply the
#   former age-specific 2,000/1,700-kcal normalization. It uses the
#   country-year GDD estimates on their reported scale and records this
#   decision explicitly in the output. Exact age-specific kcal scaling
#   requires the age-specific GDD file and must be performed before
#   aggregation.
#
# Key correction:
#   Tubers/potatoes = v03 Potatoes + v04 Other starchy vegetables.
#   Refined grains (v07) are retained in the input but are not silently
#   substituted for whole grains (v08); the four prespecified scores use
#   whole grains where applicable.
# ============================================================

required_packages <- c("dplyr", "readr", "tidyr", "tibble")
missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]
if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}
invisible(lapply(required_packages, library, character.only = TRUE))

project_path <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
input_dir <- file.path(project_path, "02_processed", "GDD")
output_dir <- input_dir
input_file <- file.path(
  input_dir,
  "GDD_country_60plus_all47_country_year.csv"
)

output_file <- file.path(
  output_dir,
  "GDD_four_dietary_patterns_country_year_v1.csv"
)
output_file_complete <- file.path(
  output_dir,
  "GDD_four_dietary_patterns_country_year_complete_v1.csv"
)
qc_summary_file <- file.path(
  output_dir,
  "QC_four_dietary_patterns_summary_v1.csv"
)
qc_component_file <- file.path(
  output_dir,
  "QC_four_dietary_patterns_component_summary_v1.csv"
)
qc_input_file <- file.path(
  output_dir,
  "QC_four_dietary_patterns_input_v1.csv"
)

analysis_years <- c(1990, 1995, 2000, 2005, 2010, 2015, 2018)
expected_countries <- 185L
expected_variables <- 47L

if (!file.exists(input_file)) {
  stop("Input file was not found: ", input_file)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gdd_long <- readr::read_csv(input_file, show_col_types = FALSE)

required_long <- c(
  "iso3", "year", "superregion2", "variable_code", "estimate",
  "population_coverage"
)
missing_long <- setdiff(required_long, names(gdd_long))
if (length(missing_long) > 0) {
  stop("Missing required input columns: ", paste(missing_long, collapse = ", "))
}

gdd_long <- gdd_long %>%
  mutate(
    iso3 = as.character(iso3),
    year = as.integer(year),
    variable_code = as.character(variable_code),
    estimate = as.numeric(estimate),
    population_coverage = as.numeric(population_coverage)
  )

expected_key_n <- n_distinct(gdd_long$iso3) *
  n_distinct(gdd_long$year) * n_distinct(gdd_long$variable_code)
duplicate_keys <- gdd_long %>%
  count(iso3, year, variable_code, name = "n") %>%
  filter(n != 1L)
if (nrow(duplicate_keys) > 0) {
  stop("Duplicate country-year-variable keys detected.")
}

observed_years <- sort(unique(as.numeric(gdd_long$year)))
if (length(observed_years) != length(analysis_years) ||
    any(observed_years != analysis_years)) {
  stop(
    "Years do not match the study years. Observed: ",
    paste(observed_years, collapse = ", ")
  )
}
if (n_distinct(gdd_long$iso3) != expected_countries) {
  warning("Expected 185 countries; observed ", n_distinct(gdd_long$iso3))
}
if (n_distinct(gdd_long$variable_code) != expected_variables) {
  warning("Expected 47 variables; observed ", n_distinct(gdd_long$variable_code))
}
if (anyNA(gdd_long$estimate)) {
  stop("The country-year GDD input contains missing estimates.")
}

food_codes <- c(
  "v01", "v02", "v03", "v04", "v05", "v06", "v07", "v08",
  "v09", "v10", "v11", "v12", "v13", "v14", "v15", "v16", "v57"
)
nutrient_codes <- c("v27", "v28", "v29", "v30", "v31", "v35", "v37")
scoring_codes <- c(food_codes, nutrient_codes)
missing_scoring_codes <- setdiff(scoring_codes, unique(gdd_long$variable_code))
if (length(missing_scoring_codes) > 0) {
  stop("Missing scoring variables: ", paste(missing_scoring_codes, collapse = ", "))
}

country_year_meta <- gdd_long %>%
  group_by(iso3, year) %>%
  summarise(
    superregion2 = first(na.omit(superregion2)),
    population_coverage = min(population_coverage, na.rm = TRUE),
    .groups = "drop"
  )

gdd_wide <- gdd_long %>%
  select(iso3, year, variable_code, estimate) %>%
  tidyr::pivot_wider(
    id_cols = c(iso3, year),
    names_from = variable_code,
    values_from = estimate
  ) %>%
  left_join(country_year_meta, by = c("iso3", "year"))

if (nrow(gdd_wide) != expected_countries * length(analysis_years)) {
  stop(
    "Unexpected country-year row count: ", nrow(gdd_wide),
    ". Expected ", expected_countries * length(analysis_years), "."
  )
}
if (anyNA(gdd_wide[scoring_codes])) {
  missing_by_var <- vapply(gdd_wide[scoring_codes], function(x) sum(is.na(x)), integer(1))
  stop(
    "Missing scoring values after reshaping: ",
    paste(names(missing_by_var)[missing_by_var > 0], collapse = ", ")
  )
}

clamp <- function(x, lower, upper) {
  pmax(lower, pmin(upper, x))
}
score_positive <- function(x, min_val, max_val, max_score = 10) {
  clamp((x - min_val) / (max_val - min_val) * max_score, 0, max_score)
}
score_reverse <- function(x, best_val, worst_val, max_score = 10) {
  clamp((worst_val - x) / (worst_val - best_val) * max_score, 0, max_score)
}
score_dash_positive <- function(x, cutpoints) {
  dplyr::case_when(
    is.na(x) ~ NA_real_,
    x <= cutpoints[1] ~ 1,
    x <= cutpoints[2] ~ 2,
    x <= cutpoints[3] ~ 3,
    x <= cutpoints[4] ~ 4,
    TRUE ~ 5
  )
}
score_dash_reverse <- function(x, cutpoints) {
  dplyr::case_when(
    is.na(x) ~ NA_real_,
    x <= cutpoints[1] ~ 5,
    x <= cutpoints[2] ~ 4,
    x <= cutpoints[3] ~ 3,
    x <= cutpoints[4] ~ 2,
    TRUE ~ 1
  )
}
score_med_positive <- function(x, cutoff) {
  dplyr::case_when(is.na(x) ~ NA_real_, x >= cutoff ~ 1, TRUE ~ 0)
}
score_med_reverse <- function(x, cutoff) {
  dplyr::case_when(is.na(x) ~ NA_real_, x <= cutoff ~ 1, TRUE ~ 0)
}

serving_g <- list(
  fruit = 110,
  nonstarchy_vegetables = 40,
  beans_legumes = 86.5,
  nuts_seeds = 29.75,
  whole_grains = 48.975,
  processed_meat = 53.705,
  unprocessed_red_meat = 85,
  seafood = 85.78,
  ssb = 368,
  fruit_juice = 209.25,
  dairy = 245
)

ahei_whole_grains_max_g <- mean(c(75, 90))
dash_cutpoints <- list(
  fruit = colMeans(rbind(c(0.71, 0.98, 1.26, 1.75), c(0.65, 0.86, 1.11, 1.52))),
  nonstarchy_vegetables = colMeans(rbind(c(0.94, 1.35, 1.74, 2.29), c(0.86, 1.24, 1.59, 2.13))),
  legumes_nuts = colMeans(rbind(c(0.25, 0.41, 0.60, 0.93), c(0.26, 0.41, 0.60, 0.95))),
  whole_grains = colMeans(rbind(c(0.27, 0.47, 0.77, 1.30), c(0.26, 0.46, 0.75, 1.29))),
  low_fat_dairy = colMeans(rbind(c(0.34, 0.66, 1.12, 1.71), c(0.31, 0.62, 1.07, 1.64))),
  red_processed_meat = colMeans(rbind(c(0.37, 0.62, 0.97, 1.49), c(0.38, 0.66, 1.10, 1.66))),
  ssb = colMeans(rbind(c(0.17, 0.37, 0.64, 1.08), c(0.19, 0.40, 0.68, 1.14))),
  sodium = colMeans(rbind(c(1857.58, 2233.79, 2615.16, 3037.78), c(1974.49, 2374.67, 2814.69, 3228.36)))
)
med_cutpoints <- list(
  fruit_nuts = mean(c(1.45, 1.29)),
  nonstarchy_vegetables = mean(c(1.54, 1.41)),
  legumes = mean(c(0.22, 0.23)),
  whole_grains = mean(c(0.61, 0.59)),
  dairy = mean(c(0.89, 0.84)),
  red_processed_meat = mean(c(0.80, 0.86)),
  seafood = mean(c(0.31, 0.30)),
  mufa_sfa = mean(c(0.98, 1.05))
)

gdd_scores <- gdd_wide %>%
  mutate(
    energy_normalization = "not_applied_after_country_year_aggregation",
    v07_refined_grains_g_day = v07,
    v08_whole_grains_g_day = v08,
    tubers_g_day = v03 + v04,
    fruits_g_day = v01,
    nonstarchy_vegetables_g_day = v02,
    beans_legumes_g_day = v05,
    nuts_seeds_g_day = v06,
    whole_grains_g_day = v08,
    processed_meat_g_day = v09,
    unprocessed_red_meat_g_day = v10,
    seafood_g_day = v11,
    eggs_g_day = v12,
    dairy_foods_g_day = v13 + v14 + v57,
    ssb_g_day = v15,
    fruit_juice_g_day = v16,
    seafood_omega3_pct = (v30 / 1000) * 9 / 2000 * 100,
    plant_omega3_pct = (v31 / 1000) * 9 / 2000 * 100,
    unsaturated_oils_pct = v28 + v29 + seafood_omega3_pct + plant_omega3_pct,
    added_sugar_pct = v35,
    ahei_fruit_serv = v01 / serving_g$fruit,
    ahei_nonstarchy_vegetables_serv = v02 / serving_g$nonstarchy_vegetables,
    ahei_whole_grains_g = v08,
    ahei_ssb_fruit_juice_serv = v15 / serving_g$ssb + v16 / serving_g$fruit_juice,
    ahei_legumes_nuts_serv = v05 / serving_g$beans_legumes + v06 / serving_g$nuts_seeds,
    ahei_red_processed_meat_serv = v10 / serving_g$unprocessed_red_meat + v09 / serving_g$processed_meat,
    ahei_seafood_omega3_mg = v30,
    ahei_pufa_pct = v29 + seafood_omega3_pct + plant_omega3_pct,
    ahei_sodium_mg = v37,
    dash_fruit_serv = v01 / serving_g$fruit,
    dash_nonstarchy_vegetables_serv = v02 / serving_g$nonstarchy_vegetables,
    dash_legumes_nuts_serv = v05 / serving_g$beans_legumes + v06 / serving_g$nuts_seeds,
    dash_whole_grains_serv = v08 / serving_g$whole_grains,
    dash_low_fat_dairy_serv = (v57 + v14 + v13) / serving_g$dairy,
    dash_red_processed_meat_serv = v10 / serving_g$unprocessed_red_meat + v09 / serving_g$processed_meat,
    dash_ssb_serv = v15 / serving_g$ssb,
    dash_sodium_mg = v37,
    med_fruit_nuts_serv = v01 / serving_g$fruit + v06 / serving_g$nuts_seeds,
    med_nonstarchy_vegetables_serv = v02 / serving_g$nonstarchy_vegetables,
    med_legumes_serv = v05 / serving_g$beans_legumes,
    med_whole_grains_serv = v08 / serving_g$whole_grains,
    med_dairy_serv = (v57 + v14 + v13) / serving_g$dairy,
    med_red_processed_meat_serv = v10 / serving_g$unprocessed_red_meat + v09 / serving_g$processed_meat,
    med_seafood_serv = v11 / serving_g$seafood,
    med_mufa_sfa_ratio = if_else(v27 > 0, v28 / v27, NA_real_)
  ) %>%
  mutate(
    mphdi_whole_grains_score = score_positive(whole_grains_g_day, 0, 82.5),
    mphdi_tubers_score = score_reverse(tubers_g_day, 50, 200),
    mphdi_nonstarchy_vegetables_score = score_positive(nonstarchy_vegetables_g_day, 0, 300),
    mphdi_fruits_score = score_positive(fruits_g_day, 0, 200),
    mphdi_dairy_foods_score = score_reverse(dairy_foods_g_day, 250, 1000),
    mphdi_red_processed_meat_score = score_reverse(processed_meat_g_day + unprocessed_red_meat_g_day, 14, 100),
    mphdi_eggs_score = score_reverse(eggs_g_day, 13, 120),
    mphdi_total_seafood_score = score_positive(seafood_g_day, 0, 28),
    mphdi_nuts_seeds_score = score_positive(nuts_seeds_g_day, 0, 50),
    mphdi_beans_legumes_score = score_positive(beans_legumes_g_day, 0, 150),
    mphdi_unsaturated_oils_score = score_positive(unsaturated_oils_pct, 3.5, 21),
    mphdi_added_sugar_score = score_reverse(added_sugar_pct, 5, 25),
    ahei_fruit_score = score_positive(ahei_fruit_serv, 0, 4),
    ahei_nonstarchy_vegetables_score = score_positive(ahei_nonstarchy_vegetables_serv, 0, 5),
    ahei_whole_grains_score = score_positive(ahei_whole_grains_g, 0, ahei_whole_grains_max_g),
    ahei_ssb_fruit_juice_score = score_reverse(ahei_ssb_fruit_juice_serv, 0, 1),
    ahei_legumes_nuts_score = score_positive(ahei_legumes_nuts_serv, 0, 1),
    ahei_red_processed_meat_score = score_reverse(ahei_red_processed_meat_serv, 0, 1.5),
    ahei_seafood_omega3_score = score_positive(ahei_seafood_omega3_mg, 0, 250),
    ahei_pufa_score = score_positive(ahei_pufa_pct, 2, 10),
    ahei_sodium_score = score_reverse(ahei_sodium_mg, 1657.16, 10292.65),
    dash_fruit_score = score_dash_positive(dash_fruit_serv, dash_cutpoints$fruit),
    dash_nonstarchy_vegetables_score = score_dash_positive(dash_nonstarchy_vegetables_serv, dash_cutpoints$nonstarchy_vegetables),
    dash_legumes_nuts_score = score_dash_positive(dash_legumes_nuts_serv, dash_cutpoints$legumes_nuts),
    dash_whole_grains_score = score_dash_positive(dash_whole_grains_serv, dash_cutpoints$whole_grains),
    dash_low_fat_dairy_score = score_dash_positive(dash_low_fat_dairy_serv, dash_cutpoints$low_fat_dairy),
    dash_red_processed_meat_score = score_dash_reverse(dash_red_processed_meat_serv, dash_cutpoints$red_processed_meat),
    dash_ssb_score = score_dash_reverse(dash_ssb_serv, dash_cutpoints$ssb),
    dash_sodium_score = score_dash_reverse(dash_sodium_mg, dash_cutpoints$sodium),
    med_fruit_nuts_score = score_med_positive(med_fruit_nuts_serv, med_cutpoints$fruit_nuts),
    med_nonstarchy_vegetables_score = score_med_positive(med_nonstarchy_vegetables_serv, med_cutpoints$nonstarchy_vegetables),
    med_legumes_score = score_med_positive(med_legumes_serv, med_cutpoints$legumes),
    med_whole_grains_score = score_med_positive(med_whole_grains_serv, med_cutpoints$whole_grains),
    med_dairy_score = score_med_reverse(med_dairy_serv, med_cutpoints$dairy),
    med_red_processed_meat_score = score_med_reverse(med_red_processed_meat_serv, med_cutpoints$red_processed_meat),
    med_seafood_score = score_med_positive(med_seafood_serv, med_cutpoints$seafood),
    med_mufa_sfa_score = score_med_positive(med_mufa_sfa_ratio, med_cutpoints$mufa_sfa)
  )

mphdi_score_vars <- c(
  "mphdi_whole_grains_score", "mphdi_tubers_score",
  "mphdi_nonstarchy_vegetables_score", "mphdi_fruits_score",
  "mphdi_dairy_foods_score", "mphdi_red_processed_meat_score",
  "mphdi_eggs_score", "mphdi_total_seafood_score",
  "mphdi_nuts_seeds_score", "mphdi_beans_legumes_score",
  "mphdi_unsaturated_oils_score", "mphdi_added_sugar_score"
)
ahei_score_vars <- c(
  "ahei_fruit_score", "ahei_nonstarchy_vegetables_score",
  "ahei_whole_grains_score", "ahei_ssb_fruit_juice_score",
  "ahei_legumes_nuts_score", "ahei_red_processed_meat_score",
  "ahei_seafood_omega3_score", "ahei_pufa_score", "ahei_sodium_score"
)
dash_score_vars <- c(
  "dash_fruit_score", "dash_nonstarchy_vegetables_score",
  "dash_legumes_nuts_score", "dash_whole_grains_score",
  "dash_low_fat_dairy_score", "dash_red_processed_meat_score",
  "dash_ssb_score", "dash_sodium_score"
)
med_score_vars <- c(
  "med_fruit_nuts_score", "med_nonstarchy_vegetables_score",
  "med_legumes_score", "med_whole_grains_score", "med_dairy_score",
  "med_red_processed_meat_score", "med_seafood_score", "med_mufa_sfa_score"
)

gdd_scores <- gdd_scores %>%
  mutate(
    mPHDI_score = rowSums(across(all_of(mphdi_score_vars)), na.rm = FALSE),
    AHEI_score_raw_0_90 = rowSums(across(all_of(ahei_score_vars)), na.rm = FALSE),
    AHEI_score = AHEI_score_raw_0_90 / 90 * 100,
    DASH_score = rowSums(across(all_of(dash_score_vars)), na.rm = FALSE),
    MED_score = rowSums(across(all_of(med_score_vars)), na.rm = FALSE),
    mPHDI_complete = !is.na(mPHDI_score),
    AHEI_complete = !is.na(AHEI_score),
    DASH_complete = !is.na(DASH_score),
    MED_complete = !is.na(MED_score),
    all_diet_scores_complete = mPHDI_complete & AHEI_complete & DASH_complete & MED_complete
  ) %>%
  arrange(iso3, year)

check_range <- function(data, variable, min_value, max_value) {
  x <- data[[variable]]
  x <- x[is.finite(x)]
  if (length(x) == 0 || min(x) < min_value || max(x) > max_value) {
    stop("Unexpected range for ", variable, ".")
  }
}
check_range(gdd_scores, "mPHDI_score", 0, 120)
check_range(gdd_scores, "AHEI_score_raw_0_90", 0, 90)
check_range(gdd_scores, "AHEI_score", 0, 100)
check_range(gdd_scores, "DASH_score", 8, 40)
check_range(gdd_scores, "MED_score", 0, 8)

component_vars <- c(mphdi_score_vars, ahei_score_vars, dash_score_vars, med_score_vars)
component_summary <- dplyr::bind_rows(lapply(component_vars, function(v) {
  x <- gdd_scores[[v]]
  tibble(
    component = v,
    minimum = min(x, na.rm = TRUE),
    maximum = max(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    missing_n = sum(is.na(x))
  )
}))

final_summary <- tibble(
  dietary_pattern = c("mPHDI", "AHEI", "DASH", "MED", "All four"),
  score_min = c(
    min(gdd_scores$mPHDI_score), min(gdd_scores$AHEI_score),
    min(gdd_scores$DASH_score), min(gdd_scores$MED_score), NA_real_
  ),
  score_max = c(
    max(gdd_scores$mPHDI_score), max(gdd_scores$AHEI_score),
    max(gdd_scores$DASH_score), max(gdd_scores$MED_score), NA_real_
  ),
  score_mean = c(
    mean(gdd_scores$mPHDI_score), mean(gdd_scores$AHEI_score),
    mean(gdd_scores$DASH_score), mean(gdd_scores$MED_score), NA_real_
  ),
  complete_records = c(
    sum(gdd_scores$mPHDI_complete), sum(gdd_scores$AHEI_complete),
    sum(gdd_scores$DASH_complete), sum(gdd_scores$MED_complete),
    sum(gdd_scores$all_diet_scores_complete)
  ),
  missing_records = c(
    sum(!gdd_scores$mPHDI_complete), sum(!gdd_scores$AHEI_complete),
    sum(!gdd_scores$DASH_complete), sum(!gdd_scores$MED_complete),
    sum(!gdd_scores$all_diet_scores_complete)
  )
)

input_qc <- tibble(
  input_file = input_file,
  countries = n_distinct(gdd_long$iso3),
  years = n_distinct(gdd_long$year),
  variables = n_distinct(gdd_long$variable_code),
  country_year_rows = nrow(gdd_wide),
  long_rows = nrow(gdd_long),
  duplicate_key_n = nrow(duplicate_keys),
  estimate_missing_n = sum(is.na(gdd_long$estimate)),
  minimum_population_coverage = min(gdd_long$population_coverage, na.rm = TRUE),
  maximum_population_coverage = max(gdd_long$population_coverage, na.rm = TRUE),
  tubers_definition = "v03 Potatoes + v04 Other starchy vegetables",
  refined_grains_handling = "v07 retained; not substituted for v08 whole grains",
  energy_normalization = "not applied after country-year aggregation"
)

readr::write_csv(gdd_scores, output_file)
readr::write_csv(gdd_scores %>% filter(all_diet_scores_complete), output_file_complete)
readr::write_csv(final_summary, qc_summary_file)
readr::write_csv(component_summary, qc_component_file)
readr::write_csv(input_qc, qc_input_file)

message("Completed: four dietary-pattern construction")
message("Rows: ", nrow(gdd_scores))
message("Countries: ", n_distinct(gdd_scores$iso3))
message("Years: ", paste(sort(unique(gdd_scores$year)), collapse = ", "))
message("Tubers: v03 + v04")
message("Energy normalization: not applied after country-year aggregation")
message("Output: ", output_file)
message("Complete output: ", output_file_complete)
message("QC summary: ", qc_summary_file)
message("QC component summary: ", qc_component_file)
message("QC input: ", qc_input_file)
