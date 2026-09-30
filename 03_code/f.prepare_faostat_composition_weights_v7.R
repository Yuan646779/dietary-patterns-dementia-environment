# =============================================================================
# f.prepare_faostat_composition_weights_v7.R
# Build fixed 2010 FAOSTAT composition weights for GDD-to-Poore mapping
#
# Input directories:
#   revised_analysis/02_processed/Environment/FAOSTAT_processed
#   revised_analysis/02_processed/Environment/LCA_processed
#
# Output directory:
#   revised_analysis/02_processed/Environment/FAOSTAT_processed
#
# Primary hierarchy for each country and GDD category:
#   1. FAOSTAT FBS current-method 2010 composition
#   2. FAOSTAT FBSH historical-method 2010 composition
#   3. Median composition among countries in the same GDD superregion
#   4. Global median composition
#
# This script creates composition weights only. It DOES NOT calculate final
# LCA coefficients or dietary environmental footprints.
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

faostat_long_file <- file.path(
  faostat_dir,
  "FAOSTAT_FBS_GDD185_study_and_overlap_years_long_v3.csv"
)
gdd_wide_file <- file.path(
  analysis_dir,
  "02_processed", "GDD", "GDD_country_60plus_all47_country_year.csv"
)
lca_candidates <- file.path(
  lca_dir,
  c(
    "OWID_Poore_LCA_four_indicators_wide_v2.csv",
    "OWID_Poore_LCA_four_indicators_wide_v1.csv"
  )
)
lca_product_file <- lca_candidates[file.exists(lca_candidates)][1L]

required_files <- c(faostat_long_file, gdd_wide_file, lca_product_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"))
}

# ---- 2. Read inputs ------------------------------------------------------

fbs <- fread(faostat_long_file, encoding = "UTF-8")
gdd_locations <- unique(
  fread(
    gdd_wide_file,
    select = c("superregion2", "iso3"),
    encoding = "UTF-8"
  )
)
lca_products <- fread(lca_product_file, encoding = "UTF-8")

if (nrow(gdd_locations) != 185L || uniqueN(gdd_locations$iso3) != 185L) {
  stop("GDD location table does not contain exactly 185 ISO3 locations.")
}

required_fbs_columns <- c(
  "source_code", "iso3", "year", "fao_item_code", "fao_item", "value"
)
if (length(setdiff(required_fbs_columns, names(fbs))) > 0L) {
  stop("FAOSTAT long table is missing required columns.")
}

required_lca_columns <- c(
  "lca_product_key", "lca_product_original",
  "ghg_kg_co2e_per_kg_food", "land_m2_per_kg_food",
  "scarcity_water_l_per_kg_food",
  "eutrophication_g_po4e_per_kg_food"
)
if (length(setdiff(required_lca_columns, names(lca_products))) > 0L) {
  stop("Poore LCA product table is missing required columns.")
}

normalize_food_name <- function(x) {
  x <- tolower(trimws(enc2utf8(x)))
  x <- gsub("&", " and ", x, fixed = TRUE)
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

fbs[, fao_item_numeric := gsub("[^0-9]", "", fao_item_code)]
fbs[, value := as.numeric(value)]

# ---- 3. Auditable FAOSTAT-to-GDD-to-Poore mapping -----------------------
#
# Aggregate FAOSTAT item codes beginning with 29 are intentionally excluded
# to avoid double counting their more detailed component items.
# Production-stage Poore coefficients do not distinguish refined from whole
# grain processing; v07 and v08 therefore use the same commodity composition.

mapping <- data.table(
  fao_item_numeric = c(
    # Fruits
    "2617", "2615", "2620", "2611", "2612", "2613", "2614",
    "2618", "2619", "2625",
    # Non-starchy vegetables
    "2601", "2602", "2605",
    # Other starchy vegetables
    "2532", "2533", "2534", "2535", "2616",
    # Beans and legumes
    "2546", "2547", "2549", "2555",
    # Nuts and seeds
    "2552", "2551",
    # Refined grains
    "2511", "2515", "2514", "2807", "2516", "2513", "2517", "2518", "2520",
    # Whole grains: same FAOSTAT composition, separate GDD category
    "2511", "2515", "2514", "2807", "2516", "2513", "2517", "2518", "2520",
    # Processed meat proxy composition
    "2731", "2732", "2733", "2734", "2735",
    # Unprocessed red meat
    "2731", "2732", "2733", "2735",
    # Seafood
    "2761", "2762", "2763", "2764", "2765", "2766", "2767", "2769"
  ),
  variable_code = c(
    rep("v01", 10), rep("v02", 3), rep("v04", 5), rep("v05", 4),
    rep("v06", 2), rep("v07", 9), rep("v08", 9), rep("v09", 5),
    rep("v10", 4), rep("v11", 8)
  ),
  variable_label = c(
    rep("Fruits", 10), rep("Non-starchy vegetables", 3),
    rep("Other starchy vegetables", 5), rep("Beans and legumes", 4),
    rep("Nuts and seeds", 2), rep("Refined grains", 9),
    rep("Whole grains", 9), rep("Total processed meats", 5),
    rep("Unprocessed red meats", 4), rep("Total seafoods", 8)
  ),
  lca_product_original = c(
    # Fruits
    "Apples", "Bananas", "Berries & Grapes",
    rep("Citrus Fruit", 4), rep("Other Fruit", 3),
    # Vegetables
    "Tomatoes", "Onions & Leeks", "Other Vegetables",
    # Other starchy vegetables: cassava is the available Poore proxy
    rep("Cassava", 5),
    # Legumes
    "Other Pulses", "Peas", "Other Pulses", "Tofu",
    # Nuts
    "Groundnuts", "Nuts",
    # Refined grains
    "Wheat & Rye", "Wheat & Rye", "Maize", "Rice", "Oatmeal",
    "Barley", "Maize", "Maize", "Wheat & Rye",
    # Whole grains
    "Wheat & Rye", "Wheat & Rye", "Maize", "Rice", "Oatmeal",
    "Barley", "Maize", "Maize", "Wheat & Rye",
    # Processed meats
    "Beef (beef herd)", "Lamb & Mutton", "Pig Meat",
    "Poultry Meat", "Beef (beef herd)",
    # Unprocessed red meats
    "Beef (beef herd)", "Lamb & Mutton", "Pig Meat", "Beef (beef herd)",
    # Seafood
    rep("Fish (farmed)", 4), "Prawns (farmed)",
    rep("Fish (farmed)", 3)
  ),
  mapping_quality = c(
    rep("direct_or_group_match", 10),
    "direct", "direct", "group_match",
    "direct", rep("proxy_cassava", 4),
    "group_match", "direct", "group_match", "proxy_tofu_for_soy",
    "direct", "group_match",
    rep("direct_or_nearest_grain_proxy", 18),
    rep("processed_form_not_observed_proxy", 5),
    rep("direct_or_group_match", 4),
    rep("production_system_proxy", 8)
  )
)

mapping[, lca_product_key := normalize_food_name(lca_product_original)]

if (anyDuplicated(mapping[, .(fao_item_numeric, variable_code)])) {
  stop("Duplicate FAOSTAT-item/GDD-category rows in the mapping table.")
}

missing_lca_products <- setdiff(
  unique(mapping$lca_product_key),
  unique(lca_products$lca_product_key)
)
if (length(missing_lca_products) > 0L) {
  stop(
    "Mapped LCA products absent from Poore table: ",
    paste(missing_lca_products, collapse = ", ")
  )
}

# Attach observed FAOSTAT names for audit.
# The old and new FAOSTAT sources can differ only in capitalization (e.g.
# "Cereals, Other" versus "Cereals, other"). Collapse labels by numeric item
# code before joining so that supplies are never duplicated.
observed_item_names <- fbs[
  ,
  .(
    fao_item_mapping_label = paste(
      sort(unique(fao_item)),
      collapse = " / "
    )
  ),
  by = fao_item_numeric
]
mapping <- observed_item_names[mapping, on = "fao_item_numeric"]

if (anyNA(mapping$fao_item_mapping_label)) {
  stop(
    "Some mapped FAOSTAT item codes were not found: ",
    paste(
      mapping[is.na(fao_item_mapping_label), fao_item_numeric],
      collapse = ", "
    )
  )
}

duplicate_mapping_keys_after_label_join <- mapping[
  , .N, by = .(fao_item_numeric, variable_code)
][N > 1L]

if (nrow(duplicate_mapping_keys_after_label_join) > 0L) {
  stop(
    "Duplicate FAOSTAT-item/GDD-category rows appeared after attaching ",
    "item labels."
  )
}

if (nrow(mapping) != 59L) {
  stop("The audited mapping contains ", nrow(mapping), " rows; expected 59.")
}

# ---- 4. Aggregate mapped 2010 supplies by source ------------------------

fbs_2010 <- fbs[year == 2010L & source_code %chin% c("FBS", "FBSH")]

mapped_supply <- merge(
  fbs_2010,
  mapping,
  by = "fao_item_numeric",
  allow.cartesian = TRUE
)

mapped_supply <- mapped_supply[
  !is.na(iso3) & iso3 %chin% gdd_locations$iso3
]

supply_product <- mapped_supply[
  ,
  .(
    supply_kg_capita_year = sum(value, na.rm = TRUE),
    contributing_fao_items = paste(sort(unique(fao_item)), collapse = "; "),
    mapping_quality = paste(sort(unique(mapping_quality)), collapse = "; ")
  ),
  by = .(
    source_code, iso3, variable_code, variable_label,
    lca_product_key, lca_product_original
  )
]

category_totals <- supply_product[
  ,
  .(category_supply = sum(supply_kg_capita_year, na.rm = TRUE)),
  by = .(source_code, iso3, variable_code)
]

# Select current FBS when the category has positive supply; otherwise FBSH.
source_choice <- category_totals[
  category_supply > 0,
  .(
    selected_source = if ("FBS" %chin% source_code) "FBS" else "FBSH"
  ),
  by = .(iso3, variable_code)
]

selected_supply <- source_choice[
  supply_product,
  on = .(iso3, variable_code, selected_source = source_code),
  nomatch = 0L
]

selected_supply[
  ,
  country_weight := supply_kg_capita_year / sum(supply_kg_capita_year),
  by = .(iso3, variable_code)
]

selected_supply <- gdd_locations[
  selected_supply,
  on = "iso3"
]

# ---- 5. Regional and global fallback compositions -----------------------

# First create explicit zero-weight rows for products absent within otherwise
# observed country-category compositions.
category_products <- unique(
  mapping[
    , .(variable_code, variable_label, lca_product_key, lca_product_original)
  ]
)

observed_country_categories <- unique(
  selected_supply[, .(iso3, superregion2, variable_code, selected_source)]
)

complete_country_grid <- merge(
  observed_country_categories,
  category_products,
  by = "variable_code",
  allow.cartesian = TRUE
)

country_weights_complete <- merge(
  complete_country_grid,
  selected_supply[
    , .(iso3, variable_code, lca_product_key, country_weight)
  ],
  by = c("iso3", "variable_code", "lca_product_key"),
  all.x = TRUE
)
country_weights_complete[is.na(country_weight), country_weight := 0]

regional_weights <- country_weights_complete[
  ,
  .(regional_raw = median(country_weight, na.rm = TRUE)),
  by = .(superregion2, variable_code, lca_product_key)
]
regional_weights[
  , regional_weight := regional_raw / sum(regional_raw),
  by = .(superregion2, variable_code)
]

global_weights <- country_weights_complete[
  ,
  .(global_raw = median(country_weight, na.rm = TRUE)),
  by = .(variable_code, lca_product_key)
]
global_weights[
  , global_weight := global_raw / sum(global_raw),
  by = variable_code
]

if (anyNA(regional_weights$regional_weight) || anyNA(global_weights$global_weight)) {
  stop("Regional or global fallback weights could not be normalized.")
}

# ---- 6. Create fixed weights for all 185 countries ----------------------

gdd_locations_cross <- copy(gdd_locations)
category_products_cross <- copy(category_products)
gdd_locations_cross[, temporary_cross_key := 1L]
category_products_cross[, temporary_cross_key := 1L]

target_grid <- merge(
  gdd_locations_cross,
  category_products_cross,
  by = "temporary_cross_key",
  allow.cartesian = TRUE
)
target_grid[, temporary_cross_key := NULL]

target_grid <- merge(
  target_grid,
  country_weights_complete[
    , .(iso3, variable_code, lca_product_key, country_weight, selected_source)
  ],
  by = c("iso3", "variable_code", "lca_product_key"),
  all.x = TRUE
)
target_grid <- merge(
  target_grid,
  regional_weights[
    , .(superregion2, variable_code, lca_product_key, regional_weight)
  ],
  by = c("superregion2", "variable_code", "lca_product_key"),
  all.x = TRUE
)
target_grid <- merge(
  target_grid,
  global_weights[, .(variable_code, lca_product_key, global_weight)],
  by = c("variable_code", "lca_product_key"),
  all.x = TRUE
)

target_grid[
  ,
  `:=`(
    final_weight = fcase(
      !is.na(country_weight), country_weight,
      !is.na(regional_weight), regional_weight,
      default = global_weight
    ),
    weight_source = fcase(
      !is.na(country_weight) & selected_source == "FBS", "country_FBS_2010",
      !is.na(country_weight) & selected_source == "FBSH", "country_FBSH_2010",
      is.na(country_weight) & !is.na(regional_weight), "superregion_median_2010",
      default = "global_median_2010"
    )
  )
]

# Re-normalize after fallback as a final numerical safeguard.
target_grid[
  , final_weight := final_weight / sum(final_weight),
  by = .(iso3, variable_code)
]

# Add the five direct one-to-one mappings for every GDD country. They are not
# composition-weighted and therefore do not require FAOSTAT supply shares.
direct_categories <- data.table(
  variable_code = c("v03", "v12", "v13", "v14", "v57"),
  variable_label = c(
    "Potatoes", "Eggs", "Cheese",
    "Yoghurt (including fermented milk)", "Total milk"
  ),
  lca_product_original = c("Potatoes", "Eggs", "Cheese", "Milk", "Milk")
)
direct_categories[, lca_product_key := normalize_food_name(lca_product_original)]
direct_locations <- copy(gdd_locations[, .(iso3, superregion2)])
direct_locations[, temporary_cross_key := 1L]
direct_categories_cross <- copy(direct_categories)
direct_categories_cross[, temporary_cross_key := 1L]
direct_target <- merge(
  direct_locations,
  direct_categories_cross,
  by = "temporary_cross_key",
  allow.cartesian = TRUE
)
direct_target[, temporary_cross_key := NULL]
direct_target[, `:=`(
  country_weight = NA_real_,
  regional_weight = NA_real_,
  global_weight = NA_real_,
  selected_source = "DIRECT",
  final_weight = 1,
  weight_source = "direct_one_to_one_2010"
)]
target_grid <- rbindlist(list(target_grid, direct_target), use.names = TRUE, fill = TRUE)

# ---- 7. Strict QC --------------------------------------------------------

weight_sums <- target_grid[
  ,
  .(
    weight_sum = sum(final_weight),
    products = .N,
    weight_source = paste(sort(unique(weight_source)), collapse = "; ")
  ),
  by = .(iso3, superregion2, variable_code, variable_label)
]

weight_sum_errors <- weight_sums[abs(weight_sum - 1) > 1e-10]

if (nrow(weight_sums) != 185L * (uniqueN(mapping$variable_code) + nrow(direct_categories))) {
  stop("Unexpected number of country-category weight groups.")
}
if (nrow(weight_sum_errors) > 0L) {
  stop("Some fixed composition weights do not sum to 1.")
}
if (anyNA(target_grid$final_weight) || any(target_grid$final_weight < 0)) {
  stop("Missing or negative final composition weights detected.")
}

source_summary <- unique(
  target_grid[, .(iso3, variable_code, weight_source)]
)[
  , .(country_category_groups = .N), by = weight_source
][order(weight_source)]

category_summary <- weight_sums[
  ,
  .(
    countries = uniqueN(iso3),
    minimum_weight_sum = min(weight_sum),
    maximum_weight_sum = max(weight_sum),
    source_types = paste(sort(unique(weight_source)), collapse = "; ")
  ),
  by = .(variable_code, variable_label)
][order(variable_code)]

# Attach LCA values for auditing only; do not calculate composite coefficients.
# Remove the duplicate display label before joining; the authoritative product
# label comes from the Poore LCA table.
target_grid_for_lca <- copy(target_grid)
target_grid_for_lca[, lca_product_original := NULL]

fixed_weights <- lca_products[
  target_grid_for_lca,
  on = "lca_product_key"
]

# ---- 8. Export -----------------------------------------------------------

fwrite(
  mapping,
  file.path(faostat_dir, "FAOSTAT_to_GDD_to_Poore_mapping_v7.csv"),
  bom = TRUE
)
fwrite(
  fixed_weights,
  file.path(faostat_dir, "FAOSTAT_2010_fixed_composition_weights_GDD185.csv"),
  bom = TRUE
)
fwrite(
  weight_sums,
  file.path(faostat_dir, "QC_FAOSTAT_2010_fixed_weight_sums_v7.csv"),
  bom = TRUE
)
fwrite(
  source_summary,
  file.path(faostat_dir, "QC_FAOSTAT_2010_weight_source_summary_v7.csv"),
  bom = TRUE
)
fwrite(
  category_summary,
  file.path(faostat_dir, "QC_FAOSTAT_2010_weight_category_summary_v7.csv"),
  bom = TRUE
)
fwrite(
  selected_supply,
  file.path(faostat_dir, "QC_FAOSTAT_2010_selected_country_supply_v7.csv"),
  bom = TRUE
)

# ---- 9. Final report -----------------------------------------------------

message("\nFixed 2010 FAOSTAT composition weights completed successfully.")
message("Analysis directory: ", analysis_dir)
message("FAOSTAT input/output directory: ", faostat_dir)
message("LCA input directory: ", lca_dir)
message("Countries: ", uniqueN(fixed_weights$iso3))
message("Weighted GDD categories: ", uniqueN(fixed_weights$variable_code))
message("Country-category groups: ", nrow(weight_sums))
message("Weight-sum errors: ", nrow(weight_sum_errors))
message("Missing final weights: ", sum(is.na(fixed_weights$final_weight)))
message("Negative final weights: ", sum(fixed_weights$final_weight < 0))
message("Weight-source distribution:")
print(source_summary)
message(
  "No final composite LCA coefficients or dietary footprints were calculated."
)
message(
  "Review FAOSTAT_to_GDD_to_Poore_mapping_v7.csv and all ",
  "QC_FAOSTAT_2010_* files next."
)
