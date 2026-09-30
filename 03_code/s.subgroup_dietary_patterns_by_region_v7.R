# Subgroup TWFE analyses and table-style forest plot by super-region.
# Final v2 master; dementia uses an exact 10-year lag and environmental
# outcomes are concurrent. Confidence intervals omit percent signs internally.

pkgs <- c("data.table", "fixest", "ggplot2", "patchwork", "scales")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install required package(s): ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table); library(fixest); library(ggplot2)
  library(patchwork); library(scales)
})

analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
input_file <- file.path(analysis_dir, "02_processed", "Master",
  "GDD_final_master_database_185_country_year_v2.csv")
output_dir <- file.path(analysis_dir, "03_analysis",
  "Subgroup_dietary_patterns_by_region_v7")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
qc_dir <- file.path(output_dir, "qc")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

years <- c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L)
regions <- c("High-income countries",
  "Central/Eastern Europe and Central Asia", "East and Southeast Asia",
  "Latin America and the Caribbean", "Middle East and North Africa",
  "South Asia", "Sub-Saharan Africa")
diets <- c(diet_mphdi_score = "mPHDI", diet_med_score = "MED",
  diet_dash_score = "DASH", diet_ahei_score = "AHEI")
increments <- c(diet_mphdi_score = 10, diet_med_score = 1,
  diet_dash_score = 5, diet_ahei_score = 10)

resolve <- function(x, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% names(x)][1L]
  if (is.na(hit) || !length(hit)) {
    if (required) stop("Cannot resolve ", label, ". Tried: ",
      paste(candidates, collapse = ", "))
    return(NA_character_)
  }
  hit
}

dt <- fread(input_file, encoding = "UTF-8", na.strings = c("", "NA", "."))
setnames(dt, sub("^\\ufeff", "", names(dt)))
dt[, iso3 := toupper(trimws(as.character(iso3)))]
dt[, year := as.integer(year)]
dt <- dt[year %in% years]
if (!all(names(diets) %in% names(dt)))
  stop("Missing dietary columns: ",
    paste(setdiff(names(diets), names(dt)), collapse = ", "))

sdi <- resolve(dt, c("sdi", "SDI"), "SDI")
daly <- resolve(dt, c("dementia_daly_asr60plus", "dementia_daly_rate_60plus"),
  "dementia DALY rate")
prev <- resolve(dt, c("dementia_prevalence_asr60plus",
  "dementia_prevalence_rate_60plus"), "dementia prevalence rate")
env <- c(
  "GHG emissions" = resolve(dt, c("env_ghg_kg_co2e_person_day", "ghg_point",
    "ghg_emissions", "GHG_emissions", "ghg"), "GHG emissions"),
  "Land use" = resolve(dt, c("env_land_m2_person_day", "land_point",
    "land_use", "Land_use", "land"), "Land use"),
  "Water-scarcity" = resolve(dt, c(
    "env_water_scarcity_l_equivalent_person_day", "scarcity_water_point",
    "scarcity_weighted_water_use", "water_scarcity", "Water_scarcity",
    "water"), "Water-scarcity"),
  "Eutrophication" = resolve(dt, c(
    "env_eutrophication_g_po4e_person_day", "eutrophication_point",
    "eutrophication_potential", "Eutrophication", "eutrophication"),
    "Eutrophication")
)
coverage <- resolve(dt, c("env_coverage_proportion",
  "included_mass_intake_coverage", "environmental_intake_coverage",
  "food_intake_coverage"), "environmental intake coverage", FALSE)
if (!"superregion2" %in% names(dt)) stop("Missing superregion2.")
region_map <- c(HIC = regions[1], FSU = regions[2], CEECA = regions[2],
  Asia = regions[3], LAC = regions[4], MENA = regions[5],
  SAARC = regions[6], SSA = regions[7])
dt[, region := unname(region_map[as.character(superregion2)])]
dt[is.na(region), region := trimws(as.character(superregion2))]
dt[, region := factor(region, levels = regions)]
if (anyNA(dt$region)) stop("Some rows have no valid super-region.")
for (v in unique(c(names(diets), sdi, daly, prev, unname(env), coverage))) {
  if (!is.na(v)) dt[, (v) := as.numeric(get(v))]
}
if (anyDuplicated(dt[, .(iso3, year)])) stop("Duplicate iso3-year keys.")

# Exact 10-year lag: exposure years 1990, 1995, 2000, and 2005.
ex <- dt[, c(list(iso3 = iso3, exposure_year = year),
  setNames(.SD, paste0("x_", names(.SD)))), .SDcols = c(names(diets), sdi)]
ex[, outcome_year := exposure_year + 10L]
oy <- dt[, .(iso3, outcome_year = year, region,
  log_daly = log(get(daly)), log_prevalence = log(get(prev)))]
health <- merge(ex, oy, by = c("iso3", "outcome_year"), all = FALSE)
health[, year := outcome_year]
for (d in names(diets))
  health[, (paste0("x_", d)) := get(paste0("x_", d)) / increments[[d]]]
health[, x_sdi := .SD[[1L]], .SDcols = paste0("x_", sdi)]
health <- health[is.finite(log_daly) & is.finite(log_prevalence)]

environment <- copy(dt)
environment[, log_ghg := log(get(env[["GHG emissions"]]))]
environment[, log_land := log(get(env[["Land use"]]))]
environment[, log_water := log(get(env[["Water-scarcity"]]))]
environment[, log_eutrophication := log(get(env[["Eutrophication"]]))]
environment[, x_sdi := .SD[[1L]], .SDcols = sdi]
for (d in names(diets))
  environment[, (paste0("x_", d)) := get(d) / increments[[d]]]
if (!is.na(coverage)) environment[, x_coverage := get(coverage) / 0.10]

fit_one <- function(x, y, term, adjust) {
  needed <- c("iso3", "year", y, term, adjust)
  z <- x[complete.cases(x[, ..needed])]
  if (nrow(z) < 10L || uniqueN(z$iso3) < 3L || uniqueN(z$year) < 2L)
    return(list(model = NULL, data = z, error = "Insufficient complete-case data"))
  rhs <- paste(c(term, adjust), collapse = " + ")
  m <- tryCatch(feols(as.formula(paste0(y, " ~ ", rhs, " | iso3 + year")),
    data = z, vcov = ~iso3, notes = FALSE), error = function(e) e)
  if (inherits(m, "error"))
    return(list(model = NULL, data = z, error = conditionMessage(m)))
  list(model = m, data = z, error = NA_character_)
}
extract <- function(m, meta, term) {
  ct <- coeftable(m); ci <- confint(m, parm = term, level = .95)
  b <- as.numeric(ct[term, "Estimate"])
  out <- copy(meta)
  out[, estimate_percent := 100 * (exp(b) - 1)]
  out[, lower_percent := 100 * (exp(as.numeric(ci[1, 1])) - 1)]
  out[, upper_percent := 100 * (exp(as.numeric(ci[1, 2])) - 1)]
  out[, p_value := as.numeric(ct[term, "Pr(>|t|)"])]
  out[, observations := nobs(m)]
  out[, countries := tryCatch(length(fixef(m)$iso3),
                              error = function(e) uniqueN(m$model$iso3))]
  out[]
}
interaction_p <- function(x, y, term, adjust) {
  needed <- c("iso3", "year", "region", y, term, adjust)
  z <- x[complete.cases(x[, ..needed])]
  if (nrow(z) < 20L) return(NA_real_)
  rhs <- paste(c(paste0(term, " * i(region)"), adjust), collapse = " + ")
  m <- tryCatch(feols(as.formula(paste0(y, " ~ ", rhs, " | iso3 + year")),
    data = z, vcov = ~iso3, notes = FALSE), error = function(e) NULL)
  if (is.null(m)) return(NA_real_)
  cn <- names(coef(m))
  int <- grep(paste0("^", term, ":region::|^region::.*:", term),
    cn, value = TRUE)
  if (!length(int)) return(NA_real_)
  w <- tryCatch(fixest::wald(m, keep = int), error = function(e) NULL)
  if (is.null(w)) return(NA_real_)
  as.numeric(if (!is.null(w$p.value)) w$p.value else w$p)[1L]
}

health_outcomes <- c("Dementia DALY rate" = "log_daly",
  "Dementia prevalence rate" = "log_prevalence")
environment_outcomes <- c("GHG emissions" = "log_ghg",
  "Land use" = "log_land", "Water-scarcity" = "log_water",
  "Eutrophication" = "log_eutrophication")
env_adjust <- c("x_sdi", if (!is.na(coverage)) "x_coverage" else character(0))
results <- list(); qc <- list()
run_family <- function(x, family, outcomes, adjust, temporal) {
  for (o in names(outcomes)) for (d in names(diets)) {
    term <- paste0("x_", d)
    pint <- interaction_p(x, outcomes[[o]], term, adjust)
    for (r in regions) {
      z <- fit_one(x[as.character(region) == r], outcomes[[o]], term, adjust)
      qc[[length(qc) + 1L]] <<- data.table(
        family = family, outcome = o, dietary_pattern = unname(diets[[d]]),
        region = r, observations = nrow(z$data),
        countries = uniqueN(z$data$iso3), years = uniqueN(z$data$year),
        fit_status = ifelse(is.null(z$model), "FAILED", "SUCCESS"),
        fit_error = z$error)
      if (!is.null(z$model)) results[[length(results) + 1L]] <<- extract(
        z$model, data.table(family = family, outcome = o,
          dietary_pattern = unname(diets[[d]]), region = r,
          exposure_increment = increments[[d]], p_interaction = pint,
          temporal_relation = temporal), term)
    }
  }
}
run_family(health, "Dementia", health_outcomes, "x_sdi", "10-year lag")
run_family(environment, "Environment", environment_outcomes, env_adjust, "Concurrent")
res <- rbindlist(results, fill = TRUE); qc <- rbindlist(qc, fill = TRUE)
fwrite(qc, file.path(qc_dir, "QC_subgroup_model_samples_and_fit_status_v7.csv"), bom = TRUE)
if (!nrow(res)) stop("No subgroup models were successfully fitted.")

fmt_p <- function(x) ifelse(is.na(x), "", ifelse(x < .001, "<0.001", sprintf("%.3f", x)))
fmt_effect <- function(a, l, h) ifelse(is.na(a), "",
  sprintf("%+.1f%% (%+.1f, %+.1f)", a, l, h))
res[, estimate_95ci := fmt_effect(estimate_percent, lower_percent, upper_percent)]
res[, p_interaction_text := fmt_p(p_interaction)]
res[, region_order := match(region, regions)]
res[, diet_order := match(dietary_pattern, unname(diets))]
setorder(res, family, outcome, region_order, diet_order)
res[, c("region_order", "diet_order") := NULL]
fwrite(res, file.path(table_dir, "TABLE_subgroup_dietary_patterns_by_region_v7.csv"), bom = TRUE)

# Forest plot: row spacing is increased to 5 units and text is enlarged.
diet_order <- unname(diets)
diet_cols <- c(mPHDI = "#A10C89", MED = "#4B3A99",
  DASH = "#268E86", AHEI = "#2E86BD")
make_panel <- function(x, selected) {
  x <- x[outcome %in% selected]
  x[, region := factor(region, levels = regions)]
  x[, dietary_pattern := factor(dietary_pattern, levels = diet_order)]
  base <- CJ(region = factor(regions, levels = regions),
    dietary_pattern = factor(diet_order, levels = diet_order))
  base[, row_id := .I]; base[, y := (nrow(base) - row_id + 1L) * 5]
  base[, region_center := mean(y), by = region]
  base[, region_label := ifelse(dietary_pattern == diet_order[1L],
    as.character(region), "")]
  z <- merge(base, x, by = c("region", "dietary_pattern"), all.x = TRUE, sort = FALSE)
  z[, low := min(lower_percent, na.rm = TRUE), by = outcome]
  z[, high := max(upper_percent, na.rm = TRUE), by = outcome]
  # Each outcome block has a separate text column and a wider forest column.
  # The expanded x-range prevents long CI strings from entering the forest.
  layout <- data.table(outcome = selected, text_x = c(.420, .920, 1.420),
    left = c(.610, 1.110, 1.610), right = c(.840, 1.340, 1.840),
    header_x = c(.725, 1.225, 1.725))
  z <- merge(z, layout, by = "outcome", all.x = TRUE)
  z[, x := left + (estimate_percent - low) / (high - low) * (right - left)]
  z[, xl := left + (lower_percent - low) / (high - low) * (right - left)]
  z[, xh := left + (upper_percent - low) / (high - low) * (right - left)]
  z[, xr := left + (0 - low) / (high - low) * (right - left)]
  z[, xl := pmax(xl, left)]; z[, xh := pmin(xh, right)]
  z[, x := pmax(pmin(x, right), left)]
  h <- unique(x[, .(outcome, dietary_pattern, p_interaction_text)])
  h <- h[order(match(dietary_pattern, diet_order)), .(
    p_header = paste0("P-interaction: ",
      paste0(as.character(dietary_pattern), "=", p_interaction_text,
             collapse = "; "))), by = outcome]
  h <- merge(layout, h, by = "outcome", all.x = TRUE)
  h[is.na(p_header), p_header := "P-interaction: "]
  n <- 28L; step <- 5; ymax <- n * step + 12L
  band <- data.table(i = seq_along(regions))
  band[, top := (n - (i - 1L) * 4L) * step + step / 2]
  band[, bottom := (n - i * 4L + 1L) * step - step / 2]
  ggplot() +
    geom_rect(data = band[i %% 2L == 1L],
      aes(xmin = 0, xmax = 1.875, ymin = bottom, ymax = top),
      fill = "#F3F6F6", colour = NA) +
    geom_text(data = unique(base[, .(region, region_label, region_center)]),
      aes(x = .015, y = region_center, label = region_label),
      hjust = 0, size = 3.75, fontface = "bold") +
    geom_text(data = base, aes(.275, y, label = dietary_pattern),
      hjust = 0, size = 3.45, fontface = "bold") +
    geom_text(data = z[is.finite(estimate_percent)],
      aes(text_x, y, label = estimate_95ci), size = 3.15, fontface = "bold") +
    geom_segment(data = z[is.finite(estimate_percent)],
      aes(x = xr, xend = xr, y = 2.5, yend = ymax - 3),
      colour = "grey45", linetype = "dashed", linewidth = .35) +
    geom_segment(data = z[is.finite(estimate_percent)],
      aes(x = xl, xend = xh, y = y, yend = y,
          colour = dietary_pattern), linewidth = .60) +
    geom_point(data = z[is.finite(estimate_percent)],
      aes(x = x, y = y, colour = dietary_pattern), size = 2.0) +
    geom_text(data = h, aes(header_x, ymax - 1.0, label = outcome),
      size = 3.55, fontface = "bold") +
    geom_text(data = h, aes(header_x, ymax - 3.2, label = p_header),
      size = 2.45, colour = "grey30") +
    annotate("text", .015, ymax - 1, label = "Super region",
      hjust = 0, size = 3.55, fontface = "bold") +
    annotate("text", .275, ymax - 1.0, label = "Dietary pattern",
      hjust = 0, size = 3.55, fontface = "bold") +
    annotate("text", .400, ymax - 5.2, label = "Percent difference (95% CI)",
      size = 2.85, fontface = "bold") +
    annotate("text", .810, ymax - 5.2, label = "Percent difference (95% CI)",
      size = 2.85, fontface = "bold") +
    annotate("text", 1.220, ymax - 5.2, label = "Percent difference (95% CI)",
      size = 2.85, fontface = "bold") +
    geom_hline(yintercept = ymax - 5.8, colour = "black", linewidth = .45) +
    scale_colour_manual(values = diet_cols, drop = FALSE) +
    coord_cartesian(xlim = c(0, 1.875), ylim = c(2.5, ymax), clip = "off") +
    theme_void(base_family = "Arial") +
    theme(legend.position = "none", plot.margin = margin(5, 10, 5, 10))
}
p1 <- make_panel(res[outcome %in% c("Dementia DALY rate",
  "Dementia prevalence rate", "GHG emissions")],
  c("Dementia DALY rate", "Dementia prevalence rate", "GHG emissions"))
p2 <- make_panel(res[outcome %in% c("Land use", "Water-scarcity",
  "Eutrophication")], c("Land use", "Water-scarcity", "Eutrophication"))
fig <- p1 / p2 + plot_layout(heights = c(1, 1))
ggsave(file.path(figure_dir, "Figure_3_subgroup_dietary_patterns_by_region_v7.png"),
  fig, width = 28, height = 15.5, units = "in", dpi = 600,
  limitsize = FALSE, bg = "white")
ggsave(file.path(figure_dir, "Figure_3_subgroup_dietary_patterns_by_region_v7.pdf"),
  fig, width = 28, height = 15.5, units = "in",
  device = grDevices::cairo_pdf, limitsize = FALSE, bg = "white")
message("Completed subgroup analyses and figure generation.")
