# ============================================================
# c.rebuild_dietary_pattern_qc_v1.R
#
# Purpose:
#   Rebuild the quality-control results for the four country-year
#   dietary-pattern scores generated in step b.
#
# Input:
#   GDD_four_dietary_patterns_country_year_v1.csv
#
# This script is QC only. It does not calculate age-standardized
# dementia outcomes, environmental footprints, or regression models.
# ============================================================

required_packages <- c("dplyr", "readr", "tidyr", "tibble", "stringr")
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
data_dir <- file.path(project_path, "02_processed", "GDD")
input_file <- file.path(
  data_dir,
  "GDD_four_dietary_patterns_country_year_v1.csv"
)

overall_output <- file.path(
  data_dir,
  "QC_dietary_patterns_overall_v1.csv"
)
year_output <- file.path(
  data_dir,
  "QC_dietary_patterns_by_year_v1.csv"
)
region_output <- file.path(
  data_dir,
  "QC_dietary_patterns_by_region_v1.csv"
)
component_output <- file.path(
  data_dir,
  "QC_dietary_pattern_components_v1.csv"
)
outlier_output <- file.path(
  data_dir,
  "QC_dietary_pattern_outliers_v1.csv"
)
correlation_output <- file.path(
  data_dir,
  "QC_dietary_pattern_correlations_v1.csv"
)
check_output <- file.path(
  data_dir,
  "QC_dietary_pattern_checks_v1.csv"
)

analysis_years <- c(1990, 1995, 2000, 2005, 2010, 2015, 2018)
expected_countries <- 185L
expected_rows <- expected_countries * length(analysis_years)

total_score_vars <- c("mPHDI_score", "AHEI_score", "DASH_score", "MED_score")

component_groups <- list(
  mPHDI = c(
    "mphdi_whole_grains_score", "mphdi_tubers_score",
    "mphdi_nonstarchy_vegetables_score", "mphdi_fruits_score",
    "mphdi_dairy_foods_score", "mphdi_red_processed_meat_score",
    "mphdi_eggs_score", "mphdi_total_seafood_score",
    "mphdi_nuts_seeds_score", "mphdi_beans_legumes_score",
    "mphdi_unsaturated_oils_score", "mphdi_added_sugar_score"
  ),
  AHEI = c(
    "ahei_fruit_score", "ahei_nonstarchy_vegetables_score",
    "ahei_whole_grains_score", "ahei_ssb_fruit_juice_score",
    "ahei_legumes_nuts_score", "ahei_red_processed_meat_score",
    "ahei_seafood_omega3_score", "ahei_pufa_score", "ahei_sodium_score"
  ),
  DASH = c(
    "dash_fruit_score", "dash_nonstarchy_vegetables_score",
    "dash_legumes_nuts_score", "dash_whole_grains_score",
    "dash_low_fat_dairy_score", "dash_red_processed_meat_score",
    "dash_ssb_score", "dash_sodium_score"
  ),
  MED = c(
    "med_fruit_nuts_score", "med_nonstarchy_vegetables_score",
    "med_legumes_score", "med_whole_grains_score", "med_dairy_score",
    "med_red_processed_meat_score", "med_seafood_score", "med_mufa_sfa_score"
  )
)
component_score_vars <- unname(unlist(component_groups, use.names = FALSE))

if (!file.exists(input_file)) {
  stop("Input file was not found: ", input_file)
}

diet <- readr::read_csv(input_file, show_col_types = FALSE)

required_vars <- c(
  "iso3", "year", "superregion2", "population_coverage",
  "v03", "v04", "tubers_g_day", "v07_refined_grains_g_day",
  "v08_whole_grains_g_day", total_score_vars,
  "mPHDI_complete", "AHEI_complete", "DASH_complete", "MED_complete",
  "all_diet_scores_complete", component_score_vars
)
missing_vars <- setdiff(required_vars, names(diet))
if (length(missing_vars) > 0) {
  stop("Missing variables in dietary-pattern output: ", paste(missing_vars, collapse = ", "))
}

diet <- diet %>%
  mutate(
    iso3 = as.character(iso3),
    year = as.numeric(year),
    superregion2 = as.character(superregion2)
  )

duplicate_keys <- diet %>%
  count(iso3, year, name = "n") %>%
  filter(n != 1L)

observed_years <- sort(unique(diet$year))
duplicate_tubers_n <- sum(abs(diet$tubers_g_day - (diet$v03 + diet$v04)) > 1e-10)
all_complete <- all(diet$all_diet_scores_complete)

check_table <- tibble(
  check = c(
    "row_count", "country_count", "year_count", "year_set",
    "unique_country_year_keys", "duplicate_country_year_keys",
    "all_diet_scores_complete", "tubers_definition",
    "population_coverage_range", "energy_normalization_flag"
  ),
  observed = c(
    as.character(nrow(diet)),
    as.character(n_distinct(diet$iso3)),
    as.character(n_distinct(diet$year)),
    paste(observed_years, collapse = ", "),
    as.character(nrow(diet) - sum(duplicate_keys$n - 1L)),
    as.character(nrow(duplicate_keys)),
    as.character(all_complete),
    as.character(duplicate_tubers_n == 0),
    paste(round(range(diet$population_coverage, na.rm = TRUE), 6), collapse = " to "),
    paste(unique(diet$energy_normalization), collapse = "; ")
  ),
  expected = c(
    as.character(expected_rows), as.character(expected_countries), "7",
    paste(analysis_years, collapse = ", "), "1295", "0", "TRUE", "TRUE",
    "approximately 0.994 to 1", "not applied after country-year aggregation"
  ),
  status = c(
    ifelse(nrow(diet) == expected_rows, "PASS", "FAIL"),
    ifelse(n_distinct(diet$iso3) == expected_countries, "PASS", "FAIL"),
    ifelse(n_distinct(diet$year) == length(analysis_years), "PASS", "FAIL"),
    ifelse(identical(observed_years, analysis_years), "PASS", "FAIL"),
    ifelse(nrow(duplicate_keys) == 0, "PASS", "FAIL"),
    ifelse(nrow(duplicate_keys) == 0, "PASS", "FAIL"),
    ifelse(all_complete, "PASS", "FAIL"),
    ifelse(duplicate_tubers_n == 0, "PASS", "FAIL"),
    ifelse(all(diet$population_coverage > 0 & diet$population_coverage <= 1.000001), "PASS", "REVIEW"),
    "DOCUMENTED"
  )
)

if (any(check_table$status == "FAIL")) {
  print(check_table %>% filter(status == "FAIL"))
  stop("One or more essential dietary-pattern QC checks failed.")
}

diet_long <- diet %>%
  select(iso3, year, all_of(total_score_vars)) %>%
  pivot_longer(
    cols = all_of(total_score_vars),
    names_to = "dietary_pattern",
    values_to = "score"
  )

summarise_scores <- function(data) {
  data %>%
    group_by(dietary_pattern) %>%
    summarise(
      records = n(),
      missing_n = sum(is.na(score)),
      minimum = min(score, na.rm = TRUE),
      p01 = quantile(score, 0.01, na.rm = TRUE),
      p25 = quantile(score, 0.25, na.rm = TRUE),
      median = median(score, na.rm = TRUE),
      mean = mean(score, na.rm = TRUE),
      p75 = quantile(score, 0.75, na.rm = TRUE),
      p99 = quantile(score, 0.99, na.rm = TRUE),
      maximum = max(score, na.rm = TRUE),
      sd = sd(score, na.rm = TRUE),
      .groups = "drop"
    )
}

overall_summary <- summarise_scores(diet_long)
year_summary <- diet_long %>%
  group_by(year) %>%
  group_modify(~ summarise_scores(.x)) %>%
  ungroup()

map_region <- function(x) {
  upper <- stringr::str_to_upper(stringr::str_trim(x))
  dplyr::case_when(
    stringr::str_detect(upper, "SAARC|SOUTH ASIA") ~ "South Asia",
    upper == "ASIA" | stringr::str_detect(upper, "EAST.*SOUTHEAST|EAST.*SOUTH.*ASIA") ~ "East and Southeast Asia",
    stringr::str_detect(upper, "SSA|SUB.?SAHARAN") ~ "Sub-Saharan Africa",
    stringr::str_detect(upper, "MENA|MIDDLE EAST|NORTH AFRICA") ~ "Middle East and North Africa",
    stringr::str_detect(upper, "LAC|LATIN AMERICA|CARIBBEAN") ~ "Latin America and the Caribbean",
    stringr::str_detect(upper, "HIC|HIGH.?INCOME") ~ "High-income countries",
    upper == "FSU" | stringr::str_detect(upper, "CEE|CENTRAL.*EASTERN|EASTERN EUROPE|CENTRAL ASIA") ~ "Central/Eastern Europe and Central Asia",
    TRUE ~ NA_character_
  )
}

diet_long_region <- diet %>%
  mutate(super_region = map_region(superregion2)) %>%
  select(super_region, year, all_of(total_score_vars)) %>%
  pivot_longer(
    cols = all_of(total_score_vars),
    names_to = "dietary_pattern",
    values_to = "score"
  )

if (anyNA(diet_long_region$super_region)) {
  stop("Some superregion2 values could not be mapped to the seven analysis regions.")
}

region_summary <- diet_long_region %>%
  group_by(super_region, dietary_pattern) %>%
  summarise(
    records = n(),
    missing_n = sum(is.na(score)),
    minimum = min(score, na.rm = TRUE),
    median = median(score, na.rm = TRUE),
    mean = mean(score, na.rm = TRUE),
    maximum = max(score, na.rm = TRUE),
    sd = sd(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(super_region, dietary_pattern)

component_summary <- bind_rows(lapply(names(component_groups), function(pattern) {
  vars <- component_groups[[pattern]]
  bind_rows(lapply(vars, function(v) {
    x <- diet[[v]]
    tibble(
      dietary_pattern = pattern,
      component = v,
      records = length(x),
      missing_n = sum(is.na(x)),
      minimum = min(x, na.rm = TRUE),
      maximum = max(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      expected_minimum = 0,
      expected_maximum = 10,
      range_status = ifelse(
        all(x >= 0 & x <= 10, na.rm = TRUE), "PASS", "FAIL"
      )
    )
  }))
}))

outlier_summary <- diet_long %>%
  group_by(dietary_pattern) %>%
  mutate(
    q1 = quantile(score, 0.25, na.rm = TRUE),
    q3 = quantile(score, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower_bound = q1 - 3 * iqr,
    upper_bound = q3 + 3 * iqr,
    outlier_flag = score < lower_bound | score > upper_bound
  ) %>%
  filter(outlier_flag) %>%
  select(iso3, year, dietary_pattern, score, lower_bound, upper_bound) %>%
  arrange(dietary_pattern, desc(score))

correlation_summary <- as.data.frame(
  cor(diet[, total_score_vars], use = "pairwise.complete.obs")
) %>%
  tibble::rownames_to_column("dietary_pattern") %>%
  pivot_longer(
    cols = -dietary_pattern,
    names_to = "comparison_pattern",
    values_to = "correlation"
  )

readr::write_csv(check_table, check_output)
readr::write_csv(overall_summary, overall_output)
readr::write_csv(year_summary, year_output)
readr::write_csv(region_summary, region_output)
readr::write_csv(component_summary, component_output)
readr::write_csv(outlier_summary, outlier_output)
readr::write_csv(correlation_summary, correlation_output)

message("Completed: dietary-pattern QC")
message("Input rows: ", nrow(diet))
message("Countries: ", n_distinct(diet$iso3))
message("Years: ", paste(observed_years, collapse = ", "))
message("All four scores complete: ", all_complete)
message("Tubers definition verified: v03 + v04")
message("QC outputs saved under: ", data_dir)
message("Next: review QC_dietary_patterns_overall_v1.csv, QC_dietary_patterns_by_year_v1.csv, and QC_dietary_pattern_components_v1.csv.")
