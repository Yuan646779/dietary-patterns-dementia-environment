# =============================================================================
# e.prepare_lca_coefficients_v1.R
#
# Purpose
#   Read the already copied OWID/Poore–Nemecek raw LCA files, standardize the
#   four environmental indicators, and assess their coverage for the 19 GDD
#   food/beverage categories prepared in step d.
#
# No download is performed in this script.
#
# Raw input directory:
#   01_raw/LCA_raw
#
# Processed output directory:
#   02_processed/Environment/LCA_processed
#
# This script does not calculate final composite coefficients or dietary
# footprints. Those are subsequent steps after the candidate mappings and
# FAOSTAT composition weights have been reviewed.
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
raw_dir <- file.path(analysis_dir, "01_raw", "LCA_raw")
processed_dir <- file.path(
  analysis_dir,
  "02_processed",
  "Environment",
  "LCA_processed"
)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

template_file <- file.path(
  processed_dir,
  "GDD19_LCA_coefficient_template_v1.csv"
)

indicator_wide_output <- file.path(
  processed_dir,
  "OWID_Poore_LCA_four_indicators_wide_v1.csv"
)
indicator_long_output <- file.path(
  processed_dir,
  "OWID_Poore_LCA_four_indicators_long_v1.csv"
)
coverage_output <- file.path(
  processed_dir,
  "QC_OWID_Poore_product_indicator_coverage_v1.csv"
)
indicator_summary_output <- file.path(
  processed_dir,
  "QC_OWID_Poore_indicator_summary_v1.csv"
)
mapping_output <- file.path(
  processed_dir,
  "GDD19_to_OWID_Poore_candidate_mapping_v1.csv"
)
mapping_status_output <- file.path(
  processed_dir,
  "GDD19_LCA_mapping_status_v1.csv"
)
unresolved_output <- file.path(
  processed_dir,
  "GDD19_unresolved_LCA_categories_v1.csv"
)
overall_qc_output <- file.path(
  processed_dir,
  "QC_OWID_Poore_GDD19_preparation_summary_v1.csv"
)

source_catalogue <- data.table(
  indicator_code = c("ghg", "land", "scarcity_water", "eutrophication"),
  indicator_label = c(
    "Greenhouse gas emissions",
    "Land use",
    "Scarcity-weighted water use",
    "Eutrophying emissions"
  ),
  owid_slug = c(
    "ghg-per-kg-poore",
    "land-use-per-kg-poore",
    "scarcity-water-per-kg-poore",
    "eutrophying-emissions-per-kg-poore"
  ),
  standardized_column = c(
    "ghg_kg_co2e_per_kg_food",
    "land_m2_per_kg_food",
    "scarcity_water_l_per_kg_food",
    "eutrophication_g_po4e_per_kg_food"
  ),
  standardized_unit = c(
    "kg CO2-eq/kg food",
    "m2/kg food",
    "L scarcity-weighted water/kg food",
    "g PO4-eq/kg food"
  )
)

if (!file.exists(template_file)) {
  stop("Step d template was not found: ", template_file)
}

normalize_food_name <- function(x) {
  x <- enc2utf8(as.character(x))
  x <- tolower(trimws(x))
  x <- gsub("&", " and ", x, fixed = TRUE)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

find_raw_indicator_file <- function(slug) {
  candidates <- list.files(
    raw_dir,
    pattern = paste0("^", slug, "_full\\.csv$"),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(candidates) != 1L) {
    stop(
      "Expected exactly one raw LCA file for ", slug,
      "; found ", length(candidates),
      ". Expected filename: ", paste0(slug, "_full.csv"),
      " under ", raw_dir
    )
  }
  candidates[[1L]]
}

read_one_indicator <- function(catalogue_row) {
  raw_file <- find_raw_indicator_file(catalogue_row$owid_slug)
  dt <- fread(raw_file, encoding = "UTF-8")

  entity_column <- intersect(c("Entity", "entity"), names(dt))
  year_column <- intersect(c("Year", "year"), names(dt))
  if (length(entity_column) != 1L || length(year_column) != 1L) {
    stop(
      "Could not uniquely identify Entity and Year in ", basename(raw_file),
      ". Columns: ", paste(names(dt), collapse = ", ")
    )
  }

  identifier_columns <- intersect(
    c("Entity", "entity", "Code", "code", "Year", "year"),
    names(dt)
  )
  value_columns <- setdiff(names(dt), identifier_columns)
  if (length(value_columns) != 1L) {
    stop(
      "Expected exactly one environmental value column in ",
      basename(raw_file), "; found: ", paste(value_columns, collapse = ", ")
    )
  }

  result <- dt[
    , .(
      lca_product_original = as.character(get(entity_column)),
      reference_year = as.integer(get(year_column)),
      indicator_value = as.numeric(get(value_columns))
    )
  ]
  result[, `:=`(
    lca_product_key = normalize_food_name(lca_product_original),
    indicator_code = catalogue_row$indicator_code,
    indicator_label = catalogue_row$indicator_label,
    standardized_unit = catalogue_row$standardized_unit,
    original_value_column = value_columns,
    raw_file = basename(raw_file)
  )]

  if (anyNA(result$lca_product_original) || anyNA(result$indicator_value)) {
    stop("Missing product names or values found in ", basename(raw_file), ".")
  }
  if (any(result$indicator_value < 0)) {
    stop("Negative environmental coefficients found in ", basename(raw_file), ".")
  }

  duplicate_keys <- result[
    , .N, by = .(lca_product_key, reference_year, indicator_code)
  ][N > 1L]
  if (nrow(duplicate_keys) > 0L) {
    stop("Duplicate product-year rows found in ", basename(raw_file), ".")
  }
  result[]
}

indicator_long <- rbindlist(
  lapply(seq_len(nrow(source_catalogue)), function(i) {
    read_one_indicator(source_catalogue[i])
  }),
  use.names = TRUE,
  fill = TRUE
)

product_name_crosswalk <- unique(
  indicator_long[, .(lca_product_key, lca_product_original, reference_year)]
)
name_conflicts <- product_name_crosswalk[
  , .(n_names = uniqueN(lca_product_original)),
  by = .(lca_product_key, reference_year)
][n_names > 1L]
if (nrow(name_conflicts) > 0L) {
  stop("Normalized LCA food-name collisions were detected.")
}

indicator_wide <- dcast(
  indicator_long,
  lca_product_key + lca_product_original + reference_year ~ indicator_code,
  value.var = "indicator_value"
)
missing_indicator_columns <- setdiff(
  source_catalogue$indicator_code,
  names(indicator_wide)
)
if (length(missing_indicator_columns) > 0L) {
  stop("Missing standardized indicator columns: ", paste(missing_indicator_columns, collapse = ", "))
}
setnames(
  indicator_wide,
  old = source_catalogue$indicator_code,
  new = source_catalogue$standardized_column
)
indicator_columns <- source_catalogue$standardized_column
indicator_wide[, n_indicators_available := rowSums(!is.na(.SD)), .SDcols = indicator_columns]
indicator_wide[, all_four_indicators_available := n_indicators_available == 4L]
setorder(indicator_wide, lca_product_original, reference_year)

product_coverage <- indicator_wide[
  , c(
    list(
      lca_product_key = lca_product_key,
      lca_product_original = lca_product_original,
      reference_year = reference_year
    ),
    setNames(
      lapply(.SD, function(x) !is.na(x)),
      paste0(names(.SD), "_available")
    ),
    list(
      n_indicators_available = n_indicators_available,
      all_four_indicators_available = all_four_indicators_available
    )
  ),
  .SDcols = indicator_columns
]

candidate_mapping <- fread(template_file, encoding = "UTF-8")
required_template_columns <- c(
  "variable_code", "variable_label", "unit",
  "candidate_lca_product", "mapping_role"
)
missing_template_columns <- setdiff(required_template_columns, names(candidate_mapping))
if (length(missing_template_columns) > 0L) {
  stop("The step d LCA template is missing: ", paste(missing_template_columns, collapse = ", "))
}

candidate_mapping[, candidate_lca_product_key := normalize_food_name(candidate_lca_product)]
candidate_mapping <- merge(
  candidate_mapping,
  product_coverage,
  by.x = "candidate_lca_product_key",
  by.y = "lca_product_key",
  all.x = TRUE,
  sort = FALSE
)
candidate_mapping[, candidate_found_in_owid := !is.na(lca_product_original)]
candidate_mapping[is.na(n_indicators_available), n_indicators_available := 0L]
candidate_mapping[is.na(all_four_indicators_available), all_four_indicators_available := FALSE]

availability_columns <- paste0(indicator_columns, "_available")
for (v in availability_columns) {
  if (!v %in% names(candidate_mapping)) candidate_mapping[, (v) := FALSE]
  candidate_mapping[, (v) := fifelse(is.na(get(v)), FALSE, get(v))]
}

category_mapping_status <- candidate_mapping[
  , .(
    proposed_candidates = .N,
    candidates_found_in_owid = sum(candidate_found_in_owid),
    candidates_with_all_four = sum(all_four_indicators_available),
    found_candidate_names = paste(sort(unique(lca_product_original[!is.na(lca_product_original)])), collapse = "; "),
    missing_candidate_names = paste(sort(unique(candidate_lca_product[!candidate_found_in_owid])), collapse = "; ")
  ),
  by = .(variable_code, variable_label, unit)
]
category_mapping_status[
  , mapping_status := fcase(
    variable_code %chin% c("v03", "v12", "v13", "v57") & candidates_with_all_four >= 1L,
    "DIRECT_SINGLE_PRODUCT_AVAILABLE",
    variable_code %chin% c("v01", "v02", "v05", "v06", "v07", "v08", "v10", "v11") & candidates_found_in_owid >= 1L,
    "COMPOSITION_WEIGHTS_REQUIRED",
    variable_code == "v09" & candidates_found_in_owid >= 1L,
    "PROXY_COMPOSITION_REQUIRES_REVIEW",
    variable_code == "v14" & candidates_found_in_owid >= 1L,
    "PROVISIONAL_MILK_PROXY_REQUIRES_SUPPLEMENT",
    variable_code == "v04" & candidates_found_in_owid >= 1L,
    "PARTIAL_COVERAGE_REQUIRES_SUPPLEMENT",
    variable_code %chin% c("v15", "v16", "v17", "v18"),
    "SUPPLEMENTARY_LCA_SOURCE_REQUIRED",
    default = "MANUAL_REVIEW_REQUIRED"
  )
]
category_mapping_status[, ready_for_final_coefficient := mapping_status == "DIRECT_SINGLE_PRODUCT_AVAILABLE"]
unresolved_categories <- category_mapping_status[ready_for_final_coefficient == FALSE]

indicator_summary <- indicator_long[
  , .(
    products = uniqueN(lca_product_key),
    years = paste(sort(unique(reference_year)), collapse = ", "),
    missing_values = sum(is.na(indicator_value)),
    minimum = min(indicator_value),
    median = median(indicator_value),
    maximum = max(indicator_value),
    original_value_column = paste(unique(original_value_column), collapse = "; "),
    unit = paste(unique(standardized_unit), collapse = "; ")
  ),
  by = .(indicator_code, indicator_label)
][order(indicator_code)]

overall_qc <- data.table(
  metric = c(
    "raw_lca_files_read", "unique_lca_products", "products_with_all_four_indicators",
    "gdd_food_categories", "candidate_mapping_rows", "categories_directly_ready",
    "categories_requiring_weights_review_or_supplement", "raw_directory"
  ),
  value = c(
    nrow(source_catalogue), uniqueN(indicator_wide$lca_product_key),
    sum(indicator_wide$all_four_indicators_available),
    uniqueN(candidate_mapping$variable_code), nrow(candidate_mapping),
    sum(category_mapping_status$ready_for_final_coefficient),
    nrow(unresolved_categories), raw_dir
  )
)

fwrite(indicator_long, indicator_long_output, bom = TRUE)
fwrite(indicator_wide, indicator_wide_output, bom = TRUE)
fwrite(product_coverage, coverage_output, bom = TRUE)
fwrite(indicator_summary, indicator_summary_output, bom = TRUE)
fwrite(candidate_mapping, mapping_output, bom = TRUE)
fwrite(category_mapping_status, mapping_status_output, bom = TRUE)
fwrite(unresolved_categories, unresolved_output, bom = TRUE)
fwrite(overall_qc, overall_qc_output, bom = TRUE)

message("Completed: OWID/Poore LCA coefficient preparation")
message("Raw directory: ", raw_dir)
message("Processed directory: ", processed_dir)
message("LCA files read: ", nrow(source_catalogue))
message("Unique LCA products: ", uniqueN(indicator_wide$lca_product_key))
message("Products with all four indicators: ", sum(indicator_wide$all_four_indicators_available))
message("GDD categories: ", uniqueN(candidate_mapping$variable_code))
message("Candidate mappings: ", nrow(candidate_mapping))
message("Directly ready categories: ", sum(category_mapping_status$ready_for_final_coefficient))
message("Categories requiring weights/review/supplement: ", nrow(unresolved_categories))
message("Next: review GDD19_LCA_mapping_status_v1.csv before constructing FAOSTAT composition weights.")
