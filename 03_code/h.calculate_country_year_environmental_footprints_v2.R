# =============================================================================
# h.calculate_country_year_environmental_footprints_v2.R
#
# Purpose
#   Calculate country-year environmental footprints from 15 GDD food groups
#   and the country-specific fixed composite LCA coefficients.
#
#   Intake is converted from g/person/day to kg/person/day before multiplying
#   by LCA coefficients per kg food. The four outputs are:
#   GHG, land use, scarcity-weighted water, and eutrophication.
#
#   Scope: selected 15 parameterized GDD food groups. This is not a total-diet
#   footprint because coffee, tea, fruit juice, and sugar-sweetened beverages
#   remain outside the current LCA parameterization.
# =============================================================================

rm(list = ls())

required_packages <- c("data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install required package(s) first: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages(library(data.table))

analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
gdd_dir <- file.path(analysis_dir, "02_processed", "GDD")
environment_dir <- file.path(
  analysis_dir, "02_processed", "Environment", "LCA_processed"
)
output_dir <- environment_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gdd_food_file <- file.path(
  output_dir,
  "GDD_country_60plus_foods_country_year_environment_v1.csv"
)
coefficient_file <- file.path(
  environment_dir,
  "GDD15_country_fixed_composite_LCA_coefficients_v3.csv"
)
if (!file.exists(gdd_food_file)) stop("GDD food input not found: ", gdd_food_file)
if (!file.exists(coefficient_file)) stop("LCA coefficient input not found: ", coefficient_file)

indicator_columns <- c(
  "ghg_kg_co2e_per_kg_food",
  "land_m2_per_kg_food",
  "scarcity_water_l_per_kg_food",
  "eutrophication_g_po4e_per_kg_food"
)
indicator_short <- c(
  ghg_kg_co2e_per_kg_food = "ghg",
  land_m2_per_kg_food = "land",
  scarcity_water_l_per_kg_food = "water_scarcity",
  eutrophication_g_po4e_per_kg_food = "eutrophication"
)

food <- fread(gdd_food_file, encoding = "UTF-8")
coefficients <- fread(coefficient_file, encoding = "UTF-8")

required_food_columns <- c(
  "superregion2", "iso3", "year", "variable_code", "variable_label",
  "unit", "estimate", "mc_lower_95", "mc_upper_95"
)
missing_food_columns <- setdiff(required_food_columns, names(food))
if (length(missing_food_columns)) {
  stop("GDD food file is missing: ", paste(missing_food_columns, collapse = ", "))
}
required_coefficient_columns <- c(
  "superregion2", "iso3", "variable_code", "variable_label",
  indicator_columns
)
missing_coefficient_columns <- setdiff(required_coefficient_columns, names(coefficients))
if (length(missing_coefficient_columns)) {
  stop("LCA coefficient file is missing: ", paste(missing_coefficient_columns, collapse = ", "))
}

parameterized_codes <- c(
  "v01", "v02", "v03", "v04", "v05", "v06", "v07", "v08",
  "v09", "v10", "v11", "v12", "v13", "v14", "v57"
)
unresolved_codes <- c("v15", "v16", "v17", "v18")
study_years <- c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L)

food[, `:=`(
  iso3 = as.character(iso3),
  year = as.integer(year),
  variable_code = as.character(variable_code),
  estimate = as.numeric(estimate),
  mc_lower_95 = as.numeric(mc_lower_95),
  mc_upper_95 = as.numeric(mc_upper_95)
)]
coefficients[, `:=`(
  iso3 = as.character(iso3),
  variable_code = as.character(variable_code)
)]

food <- food[year %in% study_years]
food_parameterized <- food[variable_code %chin% parameterized_codes]
if (nrow(food_parameterized) == 0L) stop("No parameterized GDD food records found.")

food_keys <- food[, .N, by = .(iso3, year, variable_code)][N > 1L]
if (nrow(food_keys)) stop("Duplicate GDD country-year-food keys detected.")
coefficient_keys <- coefficients[, .N, by = .(iso3, variable_code)][N > 1L]
if (nrow(coefficient_keys)) stop("Duplicate country-food LCA coefficients detected.")

coefficients <- coefficients[variable_code %chin% parameterized_codes]
if (uniqueN(coefficients$iso3) != 185L || uniqueN(coefficients$variable_code) != 15L) {
  stop("Expected 185 countries and 15 parameterized food categories in LCA coefficients.")
}

joined <- merge(
  food_parameterized,
  coefficients,
  by = c("superregion2", "iso3", "variable_code", "variable_label"),
  all.x = TRUE,
  suffixes = c("_gdd", "_lca")
)
if (nrow(joined) != nrow(food_parameterized)) {
  stop("LCA merge changed the number of GDD food records.")
}
if (anyNA(joined[, ..indicator_columns])) {
  stop("Missing LCA coefficients after joining to GDD food records.")
}

# Convert g/day to kg/day before applying per-kg LCA coefficients.
joined[, `:=`(
  intake_kg_day = estimate / 1000,
  intake_lower_kg_day = mc_lower_95 / 1000,
  intake_upper_kg_day = mc_upper_95 / 1000
)]

joined[, `:=`(
  ghg_contribution = intake_kg_day * ghg_kg_co2e_per_kg_food,
  land_contribution = intake_kg_day * land_m2_per_kg_food,
  water_scarcity_contribution = intake_kg_day * scarcity_water_l_per_kg_food,
  eutrophication_contribution = intake_kg_day * eutrophication_g_po4e_per_kg_food,
  ghg_lower = intake_lower_kg_day * ghg_kg_co2e_per_kg_food,
  ghg_upper = intake_upper_kg_day * ghg_kg_co2e_per_kg_food,
  land_lower = intake_lower_kg_day * land_m2_per_kg_food,
  land_upper = intake_upper_kg_day * land_m2_per_kg_food,
  water_scarcity_lower = intake_lower_kg_day * scarcity_water_l_per_kg_food,
  water_scarcity_upper = intake_upper_kg_day * scarcity_water_l_per_kg_food,
  eutrophication_lower = intake_lower_kg_day * eutrophication_g_po4e_per_kg_food,
  eutrophication_upper = intake_upper_kg_day * eutrophication_g_po4e_per_kg_food
)]

contribution_file <- file.path(
  output_dir, "GDD15_food_group_environmental_contributions_country_year_v2.csv"
)
footprint_file <- file.path(
  output_dir, "GDD15_environmental_footprints_country_year_v2.csv"
)
coverage_file <- file.path(
  output_dir, "QC_GDD15_environmental_coverage_country_year_v2.csv"
)
summary_file <- file.path(
  output_dir, "QC_GDD15_environmental_footprint_summary_v2.csv"
)

footprints <- joined[
  , .(
    covered_intake_g_day = sum(estimate, na.rm = TRUE),
    ghg_kg_co2e_person_day = sum(ghg_contribution, na.rm = TRUE),
    land_m2_person_day = sum(land_contribution, na.rm = TRUE),
    water_scarcity_l_equivalent_person_day = sum(water_scarcity_contribution, na.rm = TRUE),
    eutrophication_g_po4e_person_day = sum(eutrophication_contribution, na.rm = TRUE),
    ghg_lower = sum(ghg_lower, na.rm = TRUE), ghg_upper = sum(ghg_upper, na.rm = TRUE),
    land_lower = sum(land_lower, na.rm = TRUE), land_upper = sum(land_upper, na.rm = TRUE),
    water_scarcity_lower = sum(water_scarcity_lower, na.rm = TRUE),
    water_scarcity_upper = sum(water_scarcity_upper, na.rm = TRUE),
    eutrophication_lower = sum(eutrophication_lower, na.rm = TRUE),
    eutrophication_upper = sum(eutrophication_upper, na.rm = TRUE),
    n_food_groups = uniqueN(variable_code)
  ),
  by = .(superregion2, iso3, year)
]

# Comparable intake denominator excludes cup-based coffee and tea. It includes
# all available GDD g/day categories, including the four currently unresolved
# LCA categories where applicable.
comparable_food <- food[unit == "g/day"]
comparable_intake <- comparable_food[
  , .(total_comparable_intake_g_day = sum(estimate, na.rm = TRUE)),
  by = .(iso3, year)
]
coverage <- merge(
  footprints[, .(iso3, year, covered_intake_g_day)],
  comparable_intake,
  by = c("iso3", "year"), all.x = TRUE
)
coverage[, coverage_proportion := fifelse(
  total_comparable_intake_g_day > 0,
  covered_intake_g_day / total_comparable_intake_g_day,
  NA_real_
)]
coverage[, coverage_definition := "15 parameterized GDD groups / all GDD g/day groups; cup-based coffee and tea excluded"]
footprints <- merge(footprints, coverage, by = c("iso3", "year"), all.x = TRUE)
setorder(footprints, iso3, year)

footprint_summary <- footprints[
  , .(
    countries = uniqueN(iso3), years = uniqueN(year), rows = .N,
    minimum_coverage = min(coverage_proportion, na.rm = TRUE),
    median_coverage = median(coverage_proportion, na.rm = TRUE),
    maximum_coverage = max(coverage_proportion, na.rm = TRUE),
    ghg_median = median(ghg_kg_co2e_person_day, na.rm = TRUE),
    land_median = median(land_m2_person_day, na.rm = TRUE),
    water_scarcity_median = median(water_scarcity_l_equivalent_person_day, na.rm = TRUE),
    eutrophication_median = median(eutrophication_g_po4e_person_day, na.rm = TRUE)
  )
]

fwrite(joined, contribution_file, bom = TRUE)
fwrite(footprints, footprint_file, bom = TRUE)
fwrite(coverage, coverage_file, bom = TRUE)
fwrite(footprint_summary, summary_file, bom = TRUE)

message("Completed: country-year environmental footprints")
message("Countries: ", uniqueN(footprints$iso3))
message("Years: ", paste(sort(unique(footprints$year)), collapse = ", "))
message("Rows: ", format(nrow(footprints), big.mark = ","))
message("Food groups per row: 15")
message("Output: ", footprint_file)
message("Contributions: ", contribution_file)
message("Coverage QC: ", coverage_file)
message("Summary QC: ", summary_file)
message("Scope: selected 15 GDD food groups, not total diet")
