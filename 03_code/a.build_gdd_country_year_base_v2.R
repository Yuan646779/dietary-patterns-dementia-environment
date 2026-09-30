# ============================================================
# a.build_gdd_country_year_base.R
# Build the GDD country-year base dataset for adults aged >=60 years.
#
# This script performs three operations in one reproducible step:
#   1) reads the 47 GDD country-level files;
#   2) links the eight older-age strata to GBD 2023 population weights;
#   3) aggregates age-specific GDD estimates to one country-year estimate.
#
# Raw files are never modified. Outputs are written to 02_processed/GDD.
# ============================================================

required_packages <- c("data.table", "countrycode")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "))
}

library(data.table)
library(countrycode)

# ---- 1. Paths ------------------------------------------------------------
project_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
gdd_input_dir <- "F:/A-博士阶段/GDD公开数据/GDD_FinalEstimates_01102022/Country-level estimates"
gbd_input_file <- file.path(project_dir, "01_raw", "GBD", "IHME-GBD_2023_DATA-47a8000f-1.csv")
output_dir <- file.path(project_dir, "02_processed", "GDD")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(gdd_input_dir)) stop("GDD input directory not found: ", gdd_input_dir)
if (!file.exists(gbd_input_file)) stop("GBD population file not found: ", gbd_input_file)

# ---- 2. Prespecified study dimensions -----------------------------------
study_years <- c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L)
older_ages <- c(62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 92.5, 97.5)

expected_codes <- c(
  sprintf("v%02d", 1:18), "v22", "v23", sprintf("v%02d", 27:31),
  sprintf("v%02d", 33:43), sprintf("v%02d", 45:54), "v57"
)

# The first 18 variables are the food/beverage variables used later for LCA.
food_codes <- c(sprintf("v%02d", 1:18), "v57")

variable_dictionary <- data.table(
  variable_code = expected_codes,
  variable_label = c(
    "Fruits", "Non-starchy vegetables", "Potatoes", "Other starchy vegetables",
    "Beans and legumes", "Nuts and seeds", "Refined grains", "Whole grains",
    "Total processed meats", "Unprocessed red meats", "Total seafoods", "Eggs",
    "Cheese", "Yoghurt (including fermented milk)", "Sugar-sweetened beverages",
    "Fruit juices", "Coffee", "Tea", "Total carbohydrates", "Total protein",
    "Saturated fat", "Monounsaturated fatty acids", "Total omega-6 fat",
    "Seafood omega-3 fat", "Plant omega-3 fat", "Dietary cholesterol",
    "Dietary fiber", "Added sugars", "Calcium", "Dietary sodium", "Iodine",
    "Iron", "Magnesium", "Potassium", "Selenium", "Vitamin A with supplements",
    "Vitamin B1", "Vitamin B2", "Vitamin B3", "Vitamin B6", "Vitamin B9 (Folate)",
    "Vitamin B12", "Vitamin C", "Vitamin D", "Vitamin E", "Zinc", "Total milk"
  ),
  unit = c(
    rep("g/day", 16), "cups/day", "cups/day",
    "% total kcal/day", "g/day",
    "% total kcal/day", "% total kcal/day", "% total kcal/day",
    "mg/day", "mg/day", "mg/day", "g/day", "% total kcal/day",
    "mg/day", "mg/day", "ug/day", "mg/day", "mg/day", "mg/day",
    "ug/day", "ug RAE/day", "mg/day", "mg/day", "mg/day", "mg/day",
    "ug DFE/day", "ug/day", "mg/day", "ug/day", "mg/day", "mg/day",
    "g/day"
  ),
  variable_class = c(rep("food_or_beverage", 18), rep("nutrient", 28), "food_or_beverage")
)

# ---- 3. Locate exactly one source file for each GDD variable -------------
source_files <- list.files(
  gdd_input_dir, pattern = "^v[0-9]{2}_cnty\\.csv$",
  full.names = TRUE, ignore.case = TRUE
)
manifest <- data.table(
  source_file = source_files,
  variable_code = sub("_cnty\\.csv$", "", basename(source_files), ignore.case = TRUE)
)

missing_files <- setdiff(expected_codes, manifest$variable_code)
if (length(missing_files) > 0L) {
  stop("Missing GDD source file(s): ", paste(missing_files, collapse = ", "))
}
manifest <- manifest[variable_code %chin% expected_codes]
manifest <- manifest[match(expected_codes, variable_code)]
if (anyDuplicated(manifest$variable_code)) stop("Duplicate GDD variable files detected.")

# ---- 4. Read and standardize one GDD file -------------------------------
read_gdd_file <- function(file_path, variable_code) {
  header <- names(fread(file_path, nrows = 0L, showProgress = FALSE))
  id_cols <- c("superregion2", "iso3", "age", "female", "urban", "edu", "year")
  value_cols <- c("median", "lowerci_95", "upperci_95", "serving", "s_lowerci_95", "s_upperci_95")
  required <- c(id_cols, "median", "lowerci_95", "upperci_95")
  absent <- setdiff(required, header)
  if (length(absent) > 0L) stop(basename(file_path), " missing: ", paste(absent, collapse = ", "))

  dt <- fread(
    file_path,
    select = intersect(c(id_cols, value_cols), header),
    na.strings = c("NA", "", "."), showProgress = TRUE
  )
  dt <- dt[
    age %in% older_ages & female == 999 & urban == 999 & edu == 999 &
      year %in% study_years
  ]
  for (x in setdiff(value_cols, names(dt))) dt[, (x) := NA_real_]
  dt[, `:=`(variable_code = variable_code, source_file = basename(file_path))]

  duplicate_keys <- dt[, .N, by = .(iso3, year, age)][N > 1L]
  if (nrow(duplicate_keys) > 0L) stop("Duplicate GDD keys in ", basename(file_path))
  dt[]
}

message("Reading GDD country-level files...")
gdd_list <- lapply(seq_len(nrow(manifest)), function(i) {
  message(sprintf("[%02d/%02d] %s", i, nrow(manifest), basename(manifest$source_file[i])))
  read_gdd_file(manifest$source_file[i], manifest$variable_code[i])
})
gdd <- rbindlist(gdd_list, use.names = TRUE, fill = TRUE)
rm(gdd_list)

gdd <- variable_dictionary[gdd, on = "variable_code"]
setorder(gdd, iso3, year, age, variable_code)

# ---- 5. Read GBD 2023 population file ----------------------------------
message("Reading GBD 2023 population file...")
population_raw <- fread(gbd_input_file, na.strings = c("NA", "", "."), showProgress = TRUE)

required_population <- c(
  "population_group_id", "measure_id", "location_id", "location_name",
  "sex_id", "age_id", "age_name", "metric_id", "year", "val", "lower", "upper"
)
absent_population <- setdiff(required_population, names(population_raw))
if (length(absent_population) > 0L) {
  stop("GBD population file missing: ", paste(absent_population, collapse = ", "))
}

# Expected query: All population, Population, Both sexes, Number.
if (!all(unique(population_raw$population_group_id) == 1L) ||
    !all(unique(population_raw$measure_id) == 44L) ||
    !all(unique(population_raw$sex_id) == 3L) ||
    !all(unique(population_raw$metric_id) == 1L)) {
  stop("GBD population file is not restricted to All population/Both sexes/Number.")
}

age_crosswalk <- data.table(
  age_id = c(17L, 18L, 19L, 20L, 30L, 31L, 32L, 235L),
  age_name_expected = c(
    "60-64 years", "65-69 years", "70-74 years", "75-79 years",
    "80-84 years", "85-89 years", "90-94 years", "95+ years"
  ),
  age = older_ages
)

custom_match <- c(
  "Bahamas" = "BHS", "Bolivia (Plurinational State of)" = "BOL",
  "Brunei Darussalam" = "BRN", "Cabo Verde" = "CPV", "Congo" = "COG",
  "Côte d'Ivoire" = "CIV", "Czechia" = "CZE",
  "Democratic People's Republic of Korea" = "PRK",
  "Democratic Republic of the Congo" = "COD", "Eswatini" = "SWZ",
  "Gambia" = "GMB", "Iran (Islamic Republic of)" = "IRN",
  "Lao People's Democratic Republic" = "LAO",
  "Micronesia (Federated States of)" = "FSM", "North Macedonia" = "MKD",
  "Palestine" = "PSE", "Republic of Korea" = "KOR",
  "Republic of Moldova" = "MDA", "Russian Federation" = "RUS",
  "Sao Tome and Principe" = "STP", "Syrian Arab Republic" = "SYR",
  "Taiwan" = "TWN", "Türkiye" = "TUR",
  "United Republic of Tanzania" = "TZA",
  "United States of America" = "USA",
  "Venezuela (Bolivarian Republic of)" = "VEN", "Viet Nam" = "VNM"
)

crosswalk <- unique(population_raw[, .(location_id, location_name)])
crosswalk[, iso3 := countrycode(
  location_name, origin = "country.name", destination = "iso3c",
  custom_match = custom_match, warn = TRUE
)]
if (anyDuplicated(crosswalk[!is.na(iso3), iso3])) stop("Duplicate ISO3 mappings in GBD locations.")

gdd_countries <- sort(unique(gdd$iso3))
if (length(gdd_countries) != 185L) stop("Expected 185 GDD countries; found ", length(gdd_countries))

pop <- population_raw[
  year %in% study_years & age_id %in% age_crosswalk$age_id
]
pop <- crosswalk[pop, on = c("location_id", "location_name")]
pop <- age_crosswalk[pop, on = "age_id"]
pop <- pop[iso3 %chin% gdd_countries]

if (any(pop$age_name != pop$age_name_expected)) stop("GBD age-name crosswalk failed.")
setnames(pop, c("val", "lower", "upper"), c("population", "population_lower", "population_upper"))
pop <- pop[, .(location_id, location_name, iso3, year, age_id, age_name, age,
               population, population_lower, population_upper)]
if (nrow(pop) != 185L * length(study_years) * length(older_ages)) {
  stop("Unexpected number of GBD population country-year-age rows: ", nrow(pop))
}
if (any(pop$population <= 0, na.rm = TRUE) || anyNA(pop$population)) stop("Invalid GBD population estimates.")
if (any(pop$population_lower > pop$population | pop$population > pop$population_upper, na.rm = TRUE)) {
  stop("Invalid GBD population uncertainty intervals.")
}
pop[, population_60plus := sum(population), by = .(iso3, year)]
pop[, population_weight_full := population / population_60plus]

# ---- 6. Merge and calculate available-age weights ----------------------
gdd_pop <- merge(gdd, pop, by = c("iso3", "year", "age"), all.x = TRUE, sort = FALSE)
if (nrow(gdd_pop) != nrow(gdd) || anyNA(gdd_pop$population)) stop("GDD-population merge failed.")

gdd_pop[, population_available := sum(population[!is.na(median)]), by = .(iso3, year, variable_code)]
gdd_pop[, population_weight_available := fifelse(
  !is.na(median), population / population_available, NA_real_
)]

# ---- 7. Aggregate age-specific estimates to country-year ----------------
# Point estimates are population-weighted across available age strata.
# Published lower/upper limits are weighted in the same transparent manner.
# These UI columns are retained for documentation; regression models use the
# point estimate only unless a later analysis explicitly specifies otherwise.
gdd_available <- gdd_pop[!is.na(median)]
gdd_country_year <- gdd_available[, .(
  superregion2 = unique(superregion2)[1L],
  variable_label = unique(variable_label)[1L],
  variable_class = unique(variable_class)[1L],
  unit = unique(unit)[1L],
  estimate = sum(median * population_weight_available),
  lowerci_95 = sum(lowerci_95 * population_weight_available),
  upperci_95 = sum(upperci_95 * population_weight_available),
  n_available_ages = .N,
  population_covered = sum(population),
  population_60plus = unique(population_60plus)[1L],
  population_coverage = sum(population) / unique(population_60plus)[1L]
), by = .(iso3, year, variable_code)]

gdd_country_year <- variable_dictionary[gdd_country_year, on = "variable_code"]
gdd_country_year[, variable_order__ := match(variable_code, expected_codes)]
setorder(gdd_country_year, iso3, year, variable_order__)
gdd_country_year[, variable_order__ := NULL]

# ---- 8. Quality control --------------------------------------------------
expected_rows <- 185L * length(study_years) * length(expected_codes)
if (nrow(gdd_country_year) != expected_rows) {
  stop("Expected ", expected_rows, " country-year-variable rows; found ", nrow(gdd_country_year))
}
if (anyDuplicated(gdd_country_year[, .(iso3, year, variable_code)])) stop("Duplicate output keys detected.")
if (any(gdd_country_year$lowerci_95 > gdd_country_year$estimate |
        gdd_country_year$estimate > gdd_country_year$upperci_95, na.rm = TRUE)) {
  stop("Aggregated uncertainty interval order error detected.")
}

qc_summary <- gdd_country_year[, .(
  n_rows = .N, n_countries = uniqueN(iso3), n_years = uniqueN(year),
  n_missing_estimate = sum(is.na(estimate)),
  minimum_estimate = min(estimate, na.rm = TRUE), maximum_estimate = max(estimate, na.rm = TRUE),
  minimum_population_coverage = min(population_coverage, na.rm = TRUE),
  maximum_population_coverage = max(population_coverage, na.rm = TRUE)
), by = variable_code]
qc_missing <- gdd_pop[is.na(median), .(
  iso3, year, age, variable_code, variable_label, population, population_60plus
)]

# ---- 9. Write processed outputs -----------------------------------------
fwrite(gdd_country_year, file.path(output_dir, "GDD_country_60plus_all47_country_year.csv"), na = "NA")
fwrite(gdd_country_year[variable_code %chin% food_codes],
       file.path(output_dir, "GDD_country_60plus_foods_country_year.csv"), na = "NA")
fwrite(pop, file.path(output_dir, "GBD2023_population_GDD185_60plus_age_specific.csv"), na = "NA")
fwrite(crosswalk, file.path(output_dir, "GBD_location_to_ISO3_crosswalk.csv"), na = "NA")
fwrite(variable_dictionary, file.path(output_dir, "GDD_variable_dictionary.csv"), na = "NA")
fwrite(qc_summary, file.path(output_dir, "QC_GDD_country_year_variable_summary.csv"), na = "NA")
fwrite(qc_missing, file.path(output_dir, "QC_GDD_missing_age_records.csv"), na = "NA")
fwrite(manifest, file.path(output_dir, "QC_GDD_source_file_manifest.csv"), na = "NA")

message("\nCompleted: GDD country-year base dataset")
message("Countries: ", uniqueN(gdd_country_year$iso3))
message("Years: ", paste(sort(unique(gdd_country_year$year)), collapse = ", "))
message("Variables: ", uniqueN(gdd_country_year$variable_code))
message("Rows: ", format(nrow(gdd_country_year), big.mark = ","))
message("Output directory: ", output_dir)
message("Next: review QC_GDD_country_year_variable_summary.csv before constructing dietary indices.")

