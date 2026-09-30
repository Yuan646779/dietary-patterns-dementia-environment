# Descriptive tables and panel QC for final master v2.
# Version v5: restores the missing 60+ population denominator from the raw GBD
# age-specific population file when it is absent from the final master.

required_packages <- c("data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "))
}
library(data.table)

analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
input_file <- file.path(
  analysis_dir, "02_processed", "Master",
  "GDD_final_master_database_185_country_year_v2.csv"
)
population_file <- file.path(
  analysis_dir, "01_raw", "GBD",
  "GBD_age_specific_population_60plus.csv"
)
output_dir <- file.path(analysis_dir, "03_analysis", "Descriptive")
table_dir <- file.path(output_dir, "tables")
qc_dir <- file.path(output_dir, "qc")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

study_years <- c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L)
expected_rows <- 1295L
expected_countries <- 185L

file_required <- function(x) {
  if (!file.exists(x)) stop("Input not found: ", x)
  x
}
normalise_names <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

# ---- 1. Read final master --------------------------------------------------

dt <- fread(
  file_required(input_file),
  na.strings = c("", "NA", "NaN", "."),
  encoding = "UTF-8"
)
setnames(dt, names(dt), normalise_names(names(dt)))
dt[, `:=`(iso3 = toupper(trimws(as.character(iso3))), year = as.integer(year))]

diet_vars <- c(
  "diet_mphdi_score", "diet_med_score", "diet_dash_score", "diet_ahei_score"
)
health_vars <- c(
  "dementia_daly_asr60plus", "dementia_prevalence_asr60plus"
)
environment_vars <- c(
  "env_ghg_kg_co2e_person_day", "env_land_m2_person_day",
  "env_water_scarcity_l_equivalent_person_day",
  "env_eutrophication_g_po4e_person_day"
)
other_vars <- c("sdi", "env_coverage_proportion")
analysis_vars_without_population <- c(diet_vars, health_vars, environment_vars, other_vars)

required <- c("iso3", "year", "superregion2", analysis_vars_without_population)
missing <- setdiff(required, names(dt))
if (length(missing)) {
  stop("Final master is missing: ", paste(missing, collapse = ", "))
}
if (nrow(dt) != expected_rows || uniqueN(dt$iso3) != expected_countries ||
    !identical(sort(unique(dt$year)), study_years)) {
  stop("Final master is not the expected 185-country x 7-year panel.")
}
if (anyDuplicated(dt[, .(iso3, year)])) {
  stop("Duplicate country-year keys detected.")
}

# ---- 2. Restore a valid 60+ population weight -----------------------------

population_candidates <- c(
  "population_60plus_from_gbd", "population_60plus", "population_60_plus"
)
population_column <- intersect(population_candidates, names(dt))[1L]

if (!is.na(population_column)) {
  if (population_column != "population_60plus_from_gbd") {
    setnames(dt, population_column, "population_60plus_from_gbd")
  }
  population_source <- "final master population column"
} else {
  population <- fread(
    file_required(population_file),
    na.strings = c("", "NA", "NaN", "."),
    encoding = "UTF-8"
  )
  setnames(population, names(population), normalise_names(names(population)))
  population_candidates_raw <- c("population", "val", "value", "estimate")
  population_column_raw <- intersect(population_candidates_raw, names(population))[1L]
  required_population <- c("iso3", "year")
  missing_population <- setdiff(required_population, names(population))
  if (length(missing_population) || is.na(population_column_raw)) {
    stop(
      "Cannot identify the 60+ population source. Required iso3/year and a population column are missing from: ",
      population_file
    )
  }
  population[, `:=`(
    iso3 = toupper(trimws(as.character(iso3))),
    year = as.integer(year),
    population_value = as.numeric(get(population_column_raw))
  )]
  population <- population[iso3 %in% unique(dt$iso3) & year %in% study_years]
  population <- population[, .(population_60plus_from_gbd = sum(population_value, na.rm = TRUE)),
                           by = .(iso3, year)]
  if (any(population$population_60plus_from_gbd <= 0) ||
      any(!is.finite(population$population_60plus_from_gbd))) {
    stop("Raw GBD 60+ population contains non-positive or non-finite totals.")
  }
  dt <- merge(dt, population, by = c("iso3", "year"), all.x = TRUE)
  population_source <- basename(population_file)
}

if (anyNA(dt$population_60plus_from_gbd) ||
    any(!is.finite(dt$population_60plus_from_gbd)) ||
    any(dt$population_60plus_from_gbd <= 0)) {
  stop("60+ population weights are missing, non-finite, or non-positive.")
}

all_analysis_vars <- c(analysis_vars_without_population, "population_60plus_from_gbd")

# Fixed region names and order for all descriptive outputs.
region_map <- c(
  HIC = "High-income countries", FSU = "Central/Eastern Europe and Central Asia",
  CEECA = "Central/Eastern Europe and Central Asia",
  Asia = "East and Southeast Asia", LAC = "Latin America and the Caribbean",
  MENA = "Middle East and North Africa", SAARC = "South Asia",
  SSA = "Sub-Saharan Africa"
)
region_order <- c(
  "High-income countries", "Central/Eastern Europe and Central Asia",
  "East and Southeast Asia", "Latin America and the Caribbean",
  "Middle East and North Africa", "South Asia", "Sub-Saharan Africa"
)
dt[, geographical_region := unname(region_map[as.character(superregion2)])]
dt[is.na(geographical_region), geographical_region := as.character(superregion2)]

# ---- 3. Variable dictionary and descriptive statistics --------------------

variable_dictionary <- data.table(
  variable = all_analysis_vars,
  label = c(
    "mPHDI score", "Mediterranean diet score", "DASH score", "AHEI score",
    "Dementia DALY rate", "Dementia prevalence rate", "GHG emissions",
    "Land use", "Scarcity-weighted water use", "Eutrophication",
    "Socio-demographic Index", "Environmental intake coverage",
    "60+ population"
  ),
  unit = c(
    rep("score", 4L), rep("per 100,000", 2L),
    "kg CO2-eq/person/day", "m2/person/day", "L-eq/person/day",
    "g PO4-eq/person/day", "index", "proportion", "population"
  )
)
fwrite(variable_dictionary,
       file.path(table_dir, "Table_variable_dictionary_v5.csv"), bom = TRUE)

describe <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(data.table(n = 0L, mean = NA_real_, sd = NA_real_, median = NA_real_,
                      q25 = NA_real_, q75 = NA_real_, minimum = NA_real_, maximum = NA_real_))
  }
  data.table(
    n = length(x), mean = mean(x), sd = sd(x), median = median(x),
    q25 = as.numeric(quantile(x, .25)), q75 = as.numeric(quantile(x, .75)),
    minimum = min(x), maximum = max(x)
  )
}
table1 <- rbindlist(lapply(all_analysis_vars, function(v) {
  z <- describe(dt[[v]])
  z[, `:=`(variable = v,
           label = variable_dictionary[variable == v, label],
           unit = variable_dictionary[variable == v, unit])]
  z
}), fill = TRUE)
setcolorder(table1, c("variable", "label", "unit", "n", "mean", "sd",
                      "median", "q25", "q75", "minimum", "maximum"))
fwrite(table1, file.path(table_dir, "Table_1_overall_descriptive_statistics_v5.csv"), bom = TRUE)

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}
weighted_summary <- function(x) {
  out <- x[, lapply(.SD, weighted_mean, w = population_60plus_from_gbd),
           by = geographical_level, .SDcols = all_analysis_vars]
  melt(out, id.vars = "geographical_level",
       variable.name = "variable", value.name = "population_weighted_mean")
}

dt2018 <- dt[year == 2018L]
global2018 <- copy(dt2018)
global2018[, geographical_level := "Global"]
regional2018 <- dt2018[, .SD, by = geographical_region]
setnames(regional2018, "geographical_region", "geographical_level")
table_s1 <- rbindlist(list(weighted_summary(global2018), weighted_summary(regional2018)), fill = TRUE)
table_s1 <- merge(table_s1, variable_dictionary, by = "variable", all.x = TRUE)
table_s1[, geographical_level := factor(geographical_level, levels = c("Global", region_order))]
table_s1[, variable_order := match(variable, all_analysis_vars)]
setorder(table_s1, geographical_level, variable_order)
table_s1[, c("geographical_level", "variable_order") := list(as.character(geographical_level), NULL)]
setcolorder(table_s1, c("geographical_level", "variable", "label", "unit", "population_weighted_mean"))
fwrite(table_s1, file.path(table_dir, "Supplementary_Table_S1_global_regional_distributions_2018_v5.csv"), bom = TRUE)

weighted_changes <- function(x) {
  first <- x[year == 1990L, lapply(.SD, weighted_mean, w = population_60plus_from_gbd), .SDcols = all_analysis_vars]
  last <- x[year == 2018L, lapply(.SD, weighted_mean, w = population_60plus_from_gbd), .SDcols = all_analysis_vars]
  out <- data.table(variable = all_analysis_vars, estimate_1990 = as.numeric(first), estimate_2018 = as.numeric(last))
  out[, `:=`(
    absolute_change = estimate_2018 - estimate_1990,
    relative_change_percent = 100 * (estimate_2018 - estimate_1990) / abs(estimate_1990)
  )]
  out
}
change_list <- list(cbind(data.table(geographical_level = "Global"), weighted_changes(dt)))
for (r in region_order) {
  change_list[[length(change_list) + 1L]] <- cbind(
    data.table(geographical_level = r), weighted_changes(dt[geographical_region == r])
  )
}
table_s2 <- rbindlist(change_list, fill = TRUE)
table_s2 <- merge(table_s2, variable_dictionary, by = "variable", all.x = TRUE)
table_s2[, geographical_level := factor(geographical_level, levels = c("Global", region_order))]
table_s2[, variable_order := match(variable, all_analysis_vars)]
setorder(table_s2, geographical_level, variable_order)
table_s2[, c("geographical_level", "variable_order") := list(as.character(geographical_level), NULL)]
setcolorder(table_s2, c("geographical_level", "variable", "label", "unit",
                        "estimate_1990", "estimate_2018", "absolute_change",
                        "relative_change_percent"))
fwrite(table_s2, file.path(table_dir, "Supplementary_Table_S2_global_regional_changes_1990_2018_v5.csv"), bom = TRUE)

# ---- 4. QC -----------------------------------------------------------------

missingness <- data.table(
  variable = all_analysis_vars,
  n_missing = vapply(dt[, ..all_analysis_vars], function(z) sum(is.na(z)), integer(1)),
  missing_percent = vapply(dt[, ..all_analysis_vars], function(z) 100 * mean(is.na(z)), numeric(1)),
  n_countries_available = vapply(dt[, ..all_analysis_vars], function(z) uniqueN(dt$iso3[!is.na(z)]), integer(1)),
  n_years_available = vapply(dt[, ..all_analysis_vars], function(z) uniqueN(dt$year[!is.na(z)]), integer(1))
)
missingness <- merge(missingness, variable_dictionary, by = "variable", all.x = TRUE)
fwrite(missingness, file.path(qc_dir, "QC_descriptive_missingness_by_variable_v5.csv"), bom = TRUE)

panel_qc <- dt[, .(
  n_rows = .N, n_countries = uniqueN(iso3), n_years = uniqueN(year),
  min_year = min(year), max_year = max(year),
  n_complete_diet_scores = sum(complete.cases(.SD))
), by = year, .SDcols = diet_vars]
panel_qc <- merge(panel_qc,
  dt[, .(n_complete_health_outcomes = sum(complete.cases(.SD))), by = year, .SDcols = health_vars], by = "year")
panel_qc <- merge(panel_qc,
  dt[, .(n_complete_environment = sum(complete.cases(.SD))), by = year, .SDcols = environment_vars], by = "year")
fwrite(panel_qc, file.path(qc_dir, "QC_descriptive_panel_by_year_v5.csv"), bom = TRUE)

population_qc <- dt[, .(
  n_rows = .N, n_countries = uniqueN(iso3), n_years = uniqueN(year),
  min_population_60plus = min(population_60plus_from_gbd),
  max_population_60plus = max(population_60plus_from_gbd)
)]
population_qc[, population_source := population_source]
fwrite(population_qc, file.path(qc_dir, "QC_descriptive_population_weights_v5.csv"), bom = TRUE)

message("Completed: descriptive tables and panel QC v5")
message("Population source: ", population_source)
message("Rows: ", nrow(dt), "; countries: ", uniqueN(dt$iso3),
        "; years: ", paste(study_years, collapse = ", "))
message("Table 1: ", file.path(table_dir, "Table_1_overall_descriptive_statistics_v5.csv"))
message("Supplementary Table S1: ", file.path(table_dir, "Supplementary_Table_S1_global_regional_distributions_2018_v5.csv"))
message("Supplementary Table S2: ", file.path(table_dir, "Supplementary_Table_S2_global_regional_changes_1990_2018_v5.csv"))
message("QC directory: ", qc_dir)
