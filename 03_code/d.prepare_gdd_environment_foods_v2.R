# =============================================================================
# d.prepare_gdd_environment_foods_v1.R
#
# Purpose
#   Prepare the 19 GDD food/beverage categories required for the environmental
#   footprint pipeline and create an auditable GDD-to-LCA mapping template.
#
# Input
#   02_processed/GDD/GDD_country_60plus_all47_country_year.csv
#
# Output directory
#   02_processed/Environment/LCA_processed
#
# This script does not download or modify LCA/FAOSTAT files. It only prepares
# the GDD food table. The GDD lower/upper estimates are retained as approximate
# marginal bounds for later sensitivity calculations; they are not Monte Carlo
# intervals.
# =============================================================================

rm(list = ls())

required_packages <- c("data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install required package(s) first: ",
    paste(missing_packages, collapse = ", ")
  )
}
suppressPackageStartupMessages(library(data.table))

analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
gdd_dir <- file.path(analysis_dir, "02_processed", "GDD")
environment_dir <- file.path(analysis_dir, "02_processed", "Environment")
lca_processed_dir <- file.path(environment_dir, "LCA_processed")
dir.create(lca_processed_dir, recursive = TRUE, showWarnings = FALSE)

gdd_long_file <- file.path(
  gdd_dir,
  "GDD_country_60plus_all47_country_year.csv"
)

food_output_file <- file.path(
  lca_processed_dir,
  "GDD_country_60plus_foods_country_year_environment_v1.csv"
)
template_output_file <- file.path(
  lca_processed_dir,
  "GDD19_LCA_coefficient_template_v1.csv"
)
qc_output_file <- file.path(
  lca_processed_dir,
  "QC_GDD19_environment_food_table_v1.csv"
)

study_years <- c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L)
expected_countries <- 185L

if (!file.exists(gdd_long_file)) {
  stop("Required GDD input file was not found: ", gdd_long_file)
}

gdd_long <- fread(gdd_long_file, encoding = "UTF-8")

required_columns <- c(
  "superregion2", "iso3", "year", "variable_code", "variable_label",
  "unit", "estimate", "lowerci_95", "upperci_95"
)
missing_columns <- setdiff(required_columns, names(gdd_long))
if (length(missing_columns) > 0L) {
  stop("GDD input is missing required columns: ", paste(missing_columns, collapse = ", "))
}

gdd_dictionary <- data.table(
  variable_code = c(
    "v01", "v02", "v03", "v04", "v05", "v06", "v07", "v08",
    "v09", "v10", "v11", "v12", "v13", "v14", "v15", "v16",
    "v17", "v18", "v57"
  ),
  variable_label = c(
    "Fruits", "Non-starchy vegetables", "Potatoes",
    "Other starchy vegetables", "Beans and legumes", "Nuts and seeds",
    "Refined grains", "Whole grains", "Total processed meats",
    "Unprocessed red meats", "Total seafoods", "Eggs", "Cheese",
    "Yoghurt (including fermented milk)", "Sugar-sweetened beverages",
    "Fruit juices", "Coffee", "Tea", "Total milk"
  ),
  unit = c(rep("g/day", 16L), "cups/day", "cups/day", "g/day")
)

food_codes <- gdd_dictionary$variable_code
missing_food_codes <- setdiff(food_codes, unique(gdd_long$variable_code))
if (length(missing_food_codes) > 0L) {
  stop("The GDD input is missing food codes: ", paste(missing_food_codes, collapse = ", "))
}

gdd_long[, `:=`(
  iso3 = as.character(iso3),
  year = as.integer(year),
  variable_code = as.character(variable_code),
  estimate = as.numeric(estimate),
  lowerci_95 = as.numeric(lowerci_95),
  upperci_95 = as.numeric(upperci_95)
)]

gdd_food_raw <- gdd_long[
  variable_code %chin% food_codes,
  .(
    superregion2,
    iso3,
    year,
    variable_code,
    estimate,
    mc_lower_95 = lowerci_95,
    mc_median = estimate,
    mc_upper_95 = upperci_95,
    source_interval_type = "GDD 95% CI; approximate marginal bound, not Monte Carlo"
  )
]

# Add authoritative labels and units from the locked dictionary.
gdd_food <- merge(
  gdd_dictionary,
  gdd_food_raw,
  by = "variable_code",
  all.y = TRUE,
  sort = FALSE
)
setcolorder(
  gdd_food,
  c(
    "superregion2", "iso3", "year", "variable_code", "variable_label",
    "unit", "estimate", "mc_lower_95", "mc_median", "mc_upper_95",
    "source_interval_type"
  )
)
setorder(gdd_food, iso3, year, variable_code)

expected_food_rows <- expected_countries * length(study_years) * length(food_codes)
duplicate_keys <- gdd_food[, .N, by = .(iso3, year, variable_code)][N > 1L]
interval_errors <- gdd_food[
  mc_lower_95 > mc_median + 1e-10 | mc_median > mc_upper_95 + 1e-10
]
missing_values <- gdd_food[
  is.na(estimate) | is.na(mc_lower_95) |
    is.na(mc_median) | is.na(mc_upper_95)
]

if (nrow(gdd_food) != expected_food_rows) {
  stop(
    "The reconstructed food table contains ", nrow(gdd_food),
    " rows; expected ", expected_food_rows, "."
  )
}
if (uniqueN(gdd_food$iso3) != expected_countries) {
  stop("The reconstructed food table does not contain 185 countries.")
}
if (!all(sort(unique(gdd_food$year)) == study_years)) {
  stop("The reconstructed food table does not contain the seven study years.")
}
if (nrow(duplicate_keys) > 0L) {
  stop("Duplicate country-year-food keys were detected.")
}
if (nrow(missing_values) > 0L) {
  stop("Missing GDD food estimates or interval limits were detected.")
}
if (nrow(interval_errors) > 0L) {
  stop("GDD food interval ordering errors were detected.")
}

# Candidate mapping is deliberately transparent. It is not a final coefficient
# assignment; FAOSTAT composition weights and LCA product coverage are applied
# in later scripts.
candidate_mapping <- data.table(
  variable_code = c(
    rep("v01", 5L), rep("v02", 5L), "v03", "v04",
    rep("v05", 3L), rep("v06", 2L), rep("v07", 4L), rep("v08", 4L),
    rep("v09", 3L), rep("v10", 4L), rep("v11", 2L), "v12", "v13", "v14",
    rep("v15", 2L), "v16", "v17", "v18", "v57"
  ),
  candidate_lca_product = c(
    "Apples", "Bananas", "Berries & Grapes", "Citrus Fruit", "Other Fruit",
    "Tomatoes", "Brassicas", "Onions & Leeks", "Root Vegetables", "Other Vegetables",
    "Potatoes", "Cassava", "Peas", "Other Pulses", "Tofu", "Nuts", "Groundnuts",
    "Wheat & Rye", "Maize", "Rice", "Oatmeal", "Wheat & Rye", "Maize", "Rice", "Oatmeal",
    "Pig Meat", "Poultry Meat", "Beef (beef herd)", "Beef (beef herd)", "Beef (dairy herd)",
    "Lamb & Mutton", "Pig Meat", "Fish (farmed)", "Prawns (farmed)", "Eggs", "Cheese", "Milk",
    "Cane Sugar", "Beet Sugar", "Fruit juice", "Coffee", "Tea", "Milk"
  ),
  mapping_role = c(
    rep("constituent", 5L), rep("constituent", 5L), "direct", "partial_constituent_only",
    rep("constituent", 3L), rep("constituent", 2L), rep("constituent", 4L), rep("constituent", 4L),
    rep("proxy_constituent", 3L), rep("constituent", 4L), rep("constituent", 2L),
    "direct", "direct", "provisional_proxy", rep("ingredient_only", 2L),
    "unmatched_expected", "direct_if_available", "unmatched_expected", "direct"
  )
)

candidate_mapping <- merge(
  candidate_mapping,
  gdd_dictionary,
  by = "variable_code",
  all.x = TRUE,
  sort = FALSE
)
setcolorder(
  candidate_mapping,
  c("variable_code", "variable_label", "unit", "candidate_lca_product", "mapping_role")
)
candidate_mapping[, `:=`(
  coefficient_status = "TO_BE_FILLED_AFTER_LCA_PRODUCT_COVERAGE_REVIEW",
  composition_weight_status = "TO_BE_FILLED_FROM_FAOSTAT",
  coefficient_note = "Do not use candidate product directly for aggregated GDD categories without documented weighting."
)]

qc_summary <- data.table(
  check = c(
    "countries", "years", "food_categories", "rows", "duplicate_keys",
    "missing_values", "interval_order_errors", "v03_present", "v04_present",
    "coffee_tea_units", "interval_interpretation"
  ),
  observed = c(
    uniqueN(gdd_food$iso3),
    paste(sort(unique(gdd_food$year)), collapse = ", "),
    uniqueN(gdd_food$variable_code),
    nrow(gdd_food),
    nrow(duplicate_keys),
    nrow(missing_values),
    nrow(interval_errors),
    "TRUE", "TRUE",
    paste(unique(gdd_food[variable_code %chin% c("v17", "v18"), unit]), collapse = "; "),
    unique(gdd_food$source_interval_type)
  ),
  expected = c(
    185L,
    paste(study_years, collapse = ", "),
    19L,
    expected_food_rows,
    0L, 0L, 0L, "TRUE", "TRUE", "cups/day",
    "GDD 95% CI retained; not Monte Carlo"
  ),
  status = c(
    ifelse(uniqueN(gdd_food$iso3) == expected_countries, "PASS", "FAIL"),
    ifelse(all(sort(unique(gdd_food$year)) == study_years), "PASS", "FAIL"),
    ifelse(uniqueN(gdd_food$variable_code) == 19L, "PASS", "FAIL"),
    ifelse(nrow(gdd_food) == expected_food_rows, "PASS", "FAIL"),
    ifelse(nrow(duplicate_keys) == 0L, "PASS", "FAIL"),
    ifelse(nrow(missing_values) == 0L, "PASS", "FAIL"),
    ifelse(nrow(interval_errors) == 0L, "PASS", "FAIL"),
    "PASS", "PASS", "PASS", "DOCUMENTED"
  )
)

fwrite(gdd_food, food_output_file, bom = TRUE)
fwrite(candidate_mapping, template_output_file, bom = TRUE)
fwrite(qc_summary, qc_output_file, bom = TRUE)

message("Completed: GDD environmental-food preparation")
message("Countries: ", uniqueN(gdd_food$iso3))
message("Years: ", paste(sort(unique(gdd_food$year)), collapse = ", "))
message("Food categories: ", uniqueN(gdd_food$variable_code))
message("Rows: ", format(nrow(gdd_food), big.mark = ","))
message("Tubers remain separate for LCA mapping: v03 Potatoes and v04 Other starchy vegetables")
message("Food table: ", food_output_file)
message("LCA template: ", template_output_file)
message("QC file: ", qc_output_file)
