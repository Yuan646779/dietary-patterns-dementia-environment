# i.generate_reviewer2_tables_v1.R
rm(list = ls())
if (!requireNamespace("data.table", quietly = TRUE)) stop("Install data.table first.")
suppressPackageStartupMessages(library(data.table))

analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
lca_dir <- file.path(analysis_dir, "02_processed", "Environment", "LCA_processed")
coef_file <- file.path(lca_dir, "GDD15_country_fixed_composite_LCA_coefficients_v3.csv")
weight_file <- file.path(analysis_dir, "02_processed", "Environment", "FAOSTAT_processed", "FAOSTAT_2010_fixed_composition_weights_GDD185.csv")
if (!file.exists(coef_file)) stop("Missing coefficient file: ", coef_file)
if (!file.exists(weight_file)) stop("Missing weight file: ", weight_file)
coef <- fread(coef_file, encoding = "UTF-8")
weights <- fread(weight_file, encoding = "UTF-8")

table1 <- data.table(
  issue = c("LCA product intensity", "Country specificity", "Regional fallback", "Time specificity", "Imported food attribution", "Water-scarcity spatial scale", "System boundary", "Uncertainty", "Unresolved groups"),
  implementation = c(
    "Poore and Nemecek product-level intensities were used for GHG, land use, scarcity-weighted water, and eutrophication.",
    "Aggregated GDD categories use country-specific 2010 FAOSTAT food-supply composition weights.",
    "If a country lacks a positive 2010 composition, the same GDD superregion median composition is used; global fallback is used only when regional data are unavailable.",
    "Composition weights are fixed at 2010 and applied to 1990, 1995, 2000, 2005, 2010, 2015, and 2018; product intensities are not year-specific.",
    "FAOSTAT food-balance quantities are attributed to the reporting country. Import origin and production-country attribution are not separately identifiable.",
    "The Poore scarcity-weighted water coefficient is retained as reported; no additional country-year water-scarcity weighting is applied.",
    "The Poore and Nemecek LCA system boundary is retained; no separate country-specific transport model is added.",
    "LCA coefficient-level uncertainty intervals were unavailable. GDD intake bounds are propagated to contribution bounds.",
    "Coffee, tea, fruit juice, and sugar-sweetened beverages are excluded because compatible mass-based LCA coefficients or conversions are unavailable."
  ),
  quantitative_detail = c(
    "Four indicators per kg food",
    paste0(uniqueN(coef$iso3), " country-specific coefficient sets"),
    paste0(sum(unique(weights[, .(iso3, variable_code, weight_source)])$weight_source == "superregion_median_2010"), " country-category groups"),
    "2010 fixed composition reference year",
    "Country-level food availability; origin not modelled",
    "Scarcity-weighted freshwater, L-equivalents per kg food",
    "Poore and Nemecek source boundary",
    "No LCA UI; GDD bounds only",
    "Four GDD groups excluded"
  )
)

table1_category_summary <- coef[, .(
  countries = uniqueN(iso3),
  coefficient_methods = paste(sort(unique(coefficient_method)), collapse = "; "),
  evidence_tiers = paste(sort(unique(evidence_tier)), collapse = "; "),
  reference_years = paste(sort(unique(coefficient_reference_year)), collapse = "; "),
  geographies = paste(sort(unique(coefficient_geography)), collapse = "; ")
), by = .(variable_code, variable_label)]

table2 <- data.table(
  food_group_code = c("v01","v02","v03","v04","v05","v06","v07","v08","v09","v10","v11","v12","v13","v14","v57"),
  environmental_food_group = c("Fruits","Non-starchy vegetables","Potatoes","Other starchy vegetables","Beans and legumes","Nuts and seeds","Refined grains","Whole grains","Total processed meats","Unprocessed red meats","Total seafoods","Eggs","Cheese","Yoghurt (including fermented milk)","Total milk"),
  lca_status = c(rep("Composition-weighted",2),"Direct single product","Composition-weighted",rep("Composition-weighted",7),rep("Direct single product",5)),
  PHDI_mPHDI = c("Positive","Positive","Reverse (combined tubers)","Reverse (combined tubers)","Positive","Positive","Not separate","Positive","Reverse (red/processed meat)","Reverse (red/processed meat)","Positive","Reverse","Reverse (dairy)","Reverse (dairy)","Reverse (dairy)"),
  AHEI = c("Positive","Positive","Not separate","Not separate","Positive (legumes/nuts)","Positive (legumes/nuts)","Not separate","Positive","Reverse (red/processed meat)","Reverse (red/processed meat)","Positive (seafood omega-3)","Not separate","Not separate","Not separate","Not separate"),
  DASH = c("Positive","Positive","Not separate","Not separate","Positive (legumes/nuts)","Positive (legumes/nuts)","Not separate","Positive","Reverse (red/processed meat)","Reverse (red/processed meat)","Not separate","Not separate","Positive (low-fat dairy)","Positive (low-fat dairy)","Positive (low-fat dairy)"),
  MED = c("Positive (fruit/nuts)","Positive","Not separate","Not separate","Positive","Positive (fruit/nuts)","Not separate","Positive","Reverse (red/processed meat)","Reverse (red/processed meat)","Positive","Not separate","Reverse (dairy)","Reverse (dairy)","Reverse (dairy)"),
  included_in_environmental_total = "Yes"
)

notes <- data.table(note = c(
  "mPHDI tubers are v03 Potatoes plus v04 Other starchy vegetables.",
  "The environmental total uses all 15 listed groups regardless of whether each is a separate component of a particular index.",
  "v07 and v08 use the same commodity composition because compatible processing-specific separation is unavailable.",
  "v15 to v18 (sugar-sweetened beverages, fruit juices, coffee, and tea) are excluded from the current environmental total."
))

fwrite(table1, file.path(lca_dir, "Reviewer2_Table1_LCA_specificity_attribution_v1.csv"), bom = TRUE)
fwrite(table1_category_summary, file.path(lca_dir, "Reviewer2_Table1_LCA_category_summary_v1.csv"), bom = TRUE)
fwrite(table2, file.path(lca_dir, "Reviewer2_Table2_dietary_index_to_15_food_group_mapping_v1.csv"), bom = TRUE)
fwrite(notes, file.path(lca_dir, "Reviewer2_Table2_mapping_notes_v1.csv"), bom = TRUE)
message("Completed: Reviewer 2 audit tables")
