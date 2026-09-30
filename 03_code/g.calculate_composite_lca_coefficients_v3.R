# =============================================================================
# g.calculate_composite_lca_coefficients_v3.R
# Calculate country-specific fixed composite LCA coefficients for GDD foods
#
# Inputs:
#   environment analysis/FAOSTAT_processed/
#     FAOSTAT_2010_fixed_composition_weights_GDD185.csv
#   environment analysis/LCA_processed/
#     OWID_Poore_LCA_four_indicators_wide.csv
#
# Outputs:
#   environment analysis/LCA_processed/
#
# This script calculates coefficients only. It DOES NOT multiply coefficients
# by GDD intake and DOES NOT calculate dietary environmental footprints.
# Poore/OWID provides global mean point estimates, not coefficient-level UIs.
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
})

# ---- 1. Paths ------------------------------------------------------------

analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
faostat_dir <- file.path(
  analysis_dir, "02_processed", "Environment", "FAOSTAT_processed"
)
lca_dir <- file.path(
  analysis_dir, "02_processed", "Environment", "LCA_processed"
)

weights_file <- file.path(
  faostat_dir,
  "FAOSTAT_2010_fixed_composition_weights_GDD185.csv"
)
lca_candidates <- file.path(
  lca_dir,
  c(
    "OWID_Poore_LCA_four_indicators_wide_v2.csv",
    "OWID_Poore_LCA_four_indicators_wide_v1.csv"
  )
)
lca_file <- lca_candidates[file.exists(lca_candidates)][1L]
gdd_wide_file <- file.path(
  analysis_dir, "02_processed", "GDD",
  "GDD_country_60plus_all47_country_year.csv"
)

required_files <- c(weights_file, lca_file, gdd_wide_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"))
}

dir.create(lca_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Read and validate inputs ----------------------------------------

weights <- fread(weights_file, encoding = "UTF-8")
lca_products <- fread(lca_file, encoding = "UTF-8")
gdd_locations <- unique(
  fread(
    gdd_wide_file,
    select = c("superregion2", "iso3"),
    encoding = "UTF-8"
  )
)

indicator_columns <- c(
  "ghg_kg_co2e_per_kg_food",
  "land_m2_per_kg_food",
  "scarcity_water_l_per_kg_food",
  "eutrophication_g_po4e_per_kg_food"
)

required_weight_columns <- c(
  "superregion2", "iso3", "variable_code", "variable_label",
  "lca_product_key", "lca_product_original", "final_weight",
  "weight_source", indicator_columns
)

if (length(setdiff(required_weight_columns, names(weights))) > 0L) {
  stop("Fixed composition-weight file is missing required columns.")
}

if (nrow(gdd_locations) != 185L || uniqueN(gdd_locations$iso3) != 185L) {
  stop("GDD location table does not contain exactly 185 ISO3 locations.")
}

if (uniqueN(weights$iso3) != 185L ||
    uniqueN(weights$variable_code) != 15L) {
  stop("Fixed composition-weight file has unexpected country/category coverage.")
}

if (anyNA(weights$final_weight) || any(weights$final_weight < 0)) {
  stop("Missing or negative composition weights detected.")
}

weight_sums <- weights[
  , .(weight_sum = sum(final_weight)), by = .(iso3, variable_code)
]
if (any(abs(weight_sums$weight_sum - 1) > 1e-10)) {
  stop("Composition weights do not sum to 1 for every country-category.")
}

if (anyNA(weights[, ..indicator_columns])) {
  stop("Missing Poore indicators detected in the weighted input table.")
}

# ---- 3. Metadata for the 10 composition-weighted categories -------------

weighted_category_metadata <- data.table(
  variable_code = c(
    "v01", "v02", "v04", "v05", "v06",
    "v07", "v08", "v09", "v10", "v11"
  ),
  coefficient_method = c(
    "FAOSTAT_2010_fixed_composition_weighted",
    "FAOSTAT_2010_fixed_composition_weighted_coarse_vegetable_groups",
    "FAOSTAT_2010_composition_all_mapped_to_cassava_proxy",
    "FAOSTAT_2010_fixed_composition_with_tofu_proxy_for_soy",
    "FAOSTAT_2010_fixed_composition_weighted",
    "FAOSTAT_2010_fixed_commodity_composition_processing_not_distinguished",
    "FAOSTAT_2010_fixed_commodity_composition_processing_not_distinguished",
    "FAOSTAT_2010_meat_species_composition_processed_form_not_observed",
    "FAOSTAT_2010_fixed_species_composition_beef_system_proxy",
    "FAOSTAT_2010_seafood_group_composition_production_system_proxy"
  ),
  evidence_tier = c(
    "B", "B", "C", "B", "B", "B", "B", "C", "B", "C"
  ),
  primary_limitation = c(
    "FAOSTAT supply composition is a population-level proxy for older-adult intake composition.",
    "FAOSTAT does not separately identify all Poore vegetable subgroups.",
    "Sweet potato, yam, other roots and plantain use cassava LCA as a proxy.",
    "Soyabean supply uses tofu as the nearest available Poore product proxy.",
    "FAOSTAT nuts and groundnuts are used as composition proxies.",
    "FAOSTAT does not distinguish refined from whole-grain commodity shares or processing impacts.",
    "FAOSTAT does not distinguish refined from whole-grain commodity shares or processing impacts.",
    "FAOSTAT meat supply does not identify the processed fraction or processing burden.",
    "Bovine and other meat categories require simplified Poore production-system proxies.",
    "Wild capture and detailed aquaculture systems are not represented by the two Poore seafood proxies."
  )
)

# The v7 fixed-weight file contains 15 categories, including five direct
# mappings. Calculate weighted coefficients from the ten weighted categories
# only; direct mappings are appended once in the next section.
weights <- weights[
  variable_code %chin% weighted_category_metadata$variable_code
]

# ---- 4. Calculate weighted composite coefficients -----------------------

weighted_coefficients <- weights[
  ,
  .(
    ghg_kg_co2e_per_kg_food = sum(
      final_weight * ghg_kg_co2e_per_kg_food
    ),
    land_m2_per_kg_food = sum(
      final_weight * land_m2_per_kg_food
    ),
    scarcity_water_l_per_kg_food = sum(
      final_weight * scarcity_water_l_per_kg_food
    ),
    eutrophication_g_po4e_per_kg_food = sum(
      final_weight * eutrophication_g_po4e_per_kg_food
    ),
    n_lca_products = .N,
    effective_lca_products = sum(final_weight > 0),
    maximum_component_weight = max(final_weight),
    dominant_lca_product = lca_product_original[which.max(final_weight)],
    weight_source = paste(sort(unique(weight_source)), collapse = "; ")
  ),
  by = .(superregion2, iso3, variable_code, variable_label)
]

weighted_coefficients <- weighted_category_metadata[
  weighted_coefficients,
  on = "variable_code"
]
weighted_coefficients[, coefficient_reference_year := 2010L]
weighted_coefficients[, coefficient_geography := "country composition; global Poore product intensities"]
weighted_coefficients[, coefficient_uncertainty_available := FALSE]

# ---- 5. Add direct and explicitly provisional single-product matches ----

direct_mapping <- data.table(
  variable_code = c("v03", "v12", "v13", "v14", "v57"),
  variable_label = c(
    "Potatoes", "Eggs", "Cheese",
    "Yoghurt (including fermented milk)", "Total milk"
  ),
  lca_product_original = c("Potatoes", "Eggs", "Cheese", "Milk", "Milk"),
  coefficient_method = c(
    "direct_single_product",
    "direct_single_product",
    "direct_single_product",
    "provisional_milk_proxy_for_yoghurt",
    "direct_single_product"
  ),
  evidence_tier = c("A", "A", "A", "C", "A"),
  primary_limitation = c(
    "Direct Poore product match.",
    "Direct Poore product match.",
    "Direct Poore product match.",
    "Milk coefficient excludes yoghurt-specific fermentation, processing and packaging differences.",
    "Direct Poore product match."
  )
)

normalize_food_name <- function(x) {
  x <- tolower(trimws(enc2utf8(x)))
  x <- gsub("&", " and ", x, fixed = TRUE)
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

direct_mapping[, lca_product_key := normalize_food_name(lca_product_original)]

direct_coefficients_base <- lca_products[
  direct_mapping,
  on = "lca_product_key",
  nomatch = 0L
]

if (nrow(direct_coefficients_base) != nrow(direct_mapping)) {
  stop("One or more direct/provisional LCA products could not be matched.")
}

gdd_locations_cross <- copy(gdd_locations)
direct_coefficients_cross <- direct_coefficients_base[
    ,
    c(
      "variable_code", "variable_label", "coefficient_method",
      "evidence_tier", "primary_limitation", "lca_product_original",
      indicator_columns
    ),
    with = FALSE
  ]
gdd_locations_cross[, temporary_cross_key := 1L]
direct_coefficients_cross[, temporary_cross_key := 1L]

direct_coefficients <- merge(
  gdd_locations_cross,
  direct_coefficients_cross,
  by = "temporary_cross_key",
  allow.cartesian = TRUE
)
direct_coefficients[, temporary_cross_key := NULL]
direct_coefficients[, n_lca_products := 1L]
direct_coefficients[, effective_lca_products := 1L]
direct_coefficients[, maximum_component_weight := 1]
direct_coefficients[, dominant_lca_product := lca_product_original]
direct_coefficients[, weight_source := "not_applicable_direct_or_proxy"]
direct_coefficients[, coefficient_reference_year := 2010L]
direct_coefficients[, coefficient_geography := "global Poore product intensity"]
direct_coefficients[, coefficient_uncertainty_available := FALSE]

# ---- 6. Combine 15 currently parameterized GDD categories ---------------

common_columns <- c(
  "superregion2", "iso3", "variable_code", "variable_label",
  indicator_columns,
  "n_lca_products", "effective_lca_products", "maximum_component_weight",
  "dominant_lca_product", "weight_source", "coefficient_method",
  "evidence_tier", "primary_limitation", "coefficient_reference_year",
  "coefficient_geography", "coefficient_uncertainty_available"
)

country_coefficients <- rbindlist(
  list(
    weighted_coefficients[, ..common_columns],
    direct_coefficients[, ..common_columns]
  ),
  use.names = TRUE
)
setorder(country_coefficients, iso3, variable_code)

# Four beverage categories remain unresolved and are intentionally excluded.
unresolved_categories <- data.table(
  variable_code = c("v15", "v16", "v17", "v18"),
  variable_label = c(
    "Sugar-sweetened beverages", "Fruit juices", "Coffee", "Tea"
  ),
  gdd_unit = c("g/day", "g/day", "cups/day", "cups/day"),
  exclusion_reason = c(
    "Sugar ingredients do not capture beverage processing, dilution and packaging.",
    "No directly compatible Poore fruit-juice product coefficient.",
    "GDD is cups/day whereas Poore Coffee is a mass-based product coefficient; brewing conversion required.",
    "No Poore Tea coefficient and GDD is cups/day; supplementary LCA and brewing conversion required."
  ),
  required_next_input = c(
    "Beverage LCA plus sugar concentration/recipe assumptions",
    "Fruit-juice LCA compatible with beverage mass",
    "Dry coffee grams per 8-oz cup plus system-boundary-compatible LCA",
    "Dry tea grams per 8-oz cup plus system-boundary-compatible LCA"
  )
)

# ---- 7. Strict QC --------------------------------------------------------

expected_rows <- 185L * 15L
duplicate_keys <- country_coefficients[
  , .N, by = .(iso3, variable_code)
][N > 1L]

if (nrow(country_coefficients) != expected_rows) {
  stop(
    "Country coefficient table has ", nrow(country_coefficients),
    " rows; expected ", expected_rows, "."
  )
}
if (uniqueN(country_coefficients$iso3) != 185L ||
    uniqueN(country_coefficients$variable_code) != 15L) {
  stop("Unexpected country or GDD-category coverage in coefficient table.")
}
if (nrow(duplicate_keys) > 0L) {
  stop("Duplicate country-category coefficient keys detected.")
}
if (anyNA(country_coefficients[, ..indicator_columns])) {
  stop("Missing composite LCA coefficients detected.")
}
if (any(unlist(country_coefficients[, ..indicator_columns]) < 0)) {
  stop("Negative composite LCA coefficients detected.")
}

coefficient_summary <- country_coefficients[
  ,
  .(
    countries = uniqueN(iso3),
    evidence_tier = paste(sort(unique(evidence_tier)), collapse = "; "),
    coefficient_method = paste(sort(unique(coefficient_method)), collapse = "; "),
    ghg_minimum = min(ghg_kg_co2e_per_kg_food),
    ghg_median = median(ghg_kg_co2e_per_kg_food),
    ghg_maximum = max(ghg_kg_co2e_per_kg_food),
    land_minimum = min(land_m2_per_kg_food),
    land_median = median(land_m2_per_kg_food),
    land_maximum = max(land_m2_per_kg_food),
    scarcity_water_minimum = min(scarcity_water_l_per_kg_food),
    scarcity_water_median = median(scarcity_water_l_per_kg_food),
    scarcity_water_maximum = max(scarcity_water_l_per_kg_food),
    eutrophication_minimum = min(eutrophication_g_po4e_per_kg_food),
    eutrophication_median = median(eutrophication_g_po4e_per_kg_food),
    eutrophication_maximum = max(eutrophication_g_po4e_per_kg_food)
  ),
  by = .(variable_code, variable_label)
]

evidence_summary <- unique(
  country_coefficients[
    , .(variable_code, variable_label, evidence_tier, coefficient_method)
  ]
)[
  , .(categories = .N), by = evidence_tier
][order(evidence_tier)]

# ---- 8. Export -----------------------------------------------------------

fwrite(
  country_coefficients,
  file.path(lca_dir, "GDD15_country_fixed_composite_LCA_coefficients_v3.csv"),
  bom = TRUE
)
fwrite(
  unresolved_categories,
  file.path(lca_dir, "GDD4_unresolved_beverage_LCA_categories_v3.csv"),
  bom = TRUE
)
fwrite(
  coefficient_summary,
  file.path(lca_dir, "QC_GDD15_composite_LCA_coefficient_summary_v3.csv"),
  bom = TRUE
)
fwrite(
  evidence_summary,
  file.path(lca_dir, "QC_GDD15_LCA_evidence_tier_summary_v3.csv"),
  bom = TRUE
)

# ---- 9. Final report -----------------------------------------------------

message("\nCountry fixed composite LCA coefficients completed successfully.")
message("Analysis directory: ", analysis_dir)
message("FAOSTAT input directory: ", faostat_dir)
message("LCA input/output directory: ", lca_dir)
message("Countries: ", uniqueN(country_coefficients$iso3))
message("Parameterized GDD categories: ", uniqueN(country_coefficients$variable_code))
message("Coefficient rows: ", nrow(country_coefficients))
message("Duplicate keys: ", nrow(duplicate_keys))
message("Missing coefficient values: ",
        sum(is.na(country_coefficients[, ..indicator_columns])))
message("Unresolved beverage categories: ", nrow(unresolved_categories))
message("Evidence-tier distribution:")
print(evidence_summary)
message("No dietary environmental footprints were calculated.")
