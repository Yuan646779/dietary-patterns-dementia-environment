# Figure 2: vertical absolute changes, reference-style colours and 2/4/4 layout.
# v19: abbreviated bottom-row x labels rotated 45 degrees to prevent overlap,
# with equal panel heights, equal two-line strips, and enlarged figure fonts.
# Data calculations retain the supplied script's population-weighted summaries.
# Edit analysis_dir below, then source this complete script in R.
# Global dashed lines represent 2018 minus 1990, not the 2018 level.
# No uncertainty/target shading is inferred or drawn.
# Output: PDF + 600-dpi PNG and the underlying plot-data CSV.

required_packages <- c("data.table", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
input_file <- file.path(analysis_dir, "02_processed", "Master", "GDD_final_master_database_185_country_year_v2.csv")
population_file <- file.path(analysis_dir, "01_raw", "GBD", "GBD_age_specific_population_60plus.csv")
output_dir <- file.path(analysis_dir, "03_analysis", "Figures", "Global_regional_trends_v19_angled_xlabels_244")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

study_years <- c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L)
if (!file.exists(input_file)) stop("Input file not found: ", input_file)
dt <- fread(input_file, encoding = "UTF-8", na.strings = c("", "NA", "."))
dt[, `:=`(iso3 = toupper(trimws(as.character(iso3))), year = as.integer(year))]

# The v2 master may not retain the 60+ denominator. Restore it from the raw
# GBD age-specific population file for population-weighted summaries.
population_candidates <- c("population_60plus_from_gbd", "population_60plus", "population_60_plus")
population_column <- intersect(population_candidates, names(dt))[1L]
if (is.na(population_column)) {
  pop <- fread(population_file, encoding = "UTF-8", na.strings = c("", "NA", "."))
  names(pop) <- tolower(gsub("[^a-zA-Z0-9]+", "_", sub("^\\ufeff", "", names(pop))))
  pop_value <- intersect(c("population", "val", "value", "estimate"), names(pop))[1L]
  if (!all(c("iso3", "year") %in% names(pop)) || is.na(pop_value)) {
    stop("Cannot identify iso3, year, and population columns in: ", population_file)
  }
  pop[, `:=`(
    iso3 = toupper(trimws(as.character(iso3))),
    year = as.integer(year),
    population_value = as.numeric(get(pop_value))
  )]
  pop <- pop[iso3 %in% unique(dt$iso3) & year %in% study_years,
             .(population_60plus_from_gbd = sum(population_value, na.rm = TRUE)),
             by = .(iso3, year)]
  dt <- merge(dt, pop, by = c("iso3", "year"), all.x = TRUE)
} else if (population_column != "population_60plus_from_gbd") {
  setnames(dt, population_column, "population_60plus_from_gbd")
}
if (anyNA(dt$population_60plus_from_gbd) || any(dt$population_60plus_from_gbd <= 0)) {
  stop("Invalid 60+ population weights in the master or raw GBD population file.")
}
if (nrow(dt) != 1295L || uniqueN(dt$iso3) != 185L || !identical(sort(unique(dt$year)), study_years)) {
  stop("Unexpected final panel size or study years.")
}
if (anyDuplicated(dt[, .(iso3, year)])) stop("Duplicate iso3-year keys detected.")

vars <- data.table(
  variable = c(
    "diet_mphdi_score", "diet_med_score", "diet_dash_score", "diet_ahei_score",
    "dementia_daly_asr60plus", "dementia_prevalence_asr60plus",
    "env_ghg_kg_co2e_person_day", "env_land_m2_person_day",
    "env_water_scarcity_l_equivalent_person_day", "env_eutrophication_g_po4e_person_day"
  ),
  panel = LETTERS[1:10],
  title = c("mPHDI", "MED", "DASH", "AHEI", "Dementia DALY rate", "Dementia prevalence rate",
            "GHG emissions", "Land use", "Water scarcity", "Eutrophication"),
  y_label = c("Population-weighted mPHDI score", "Population-weighted MED score",
              "Population-weighted DASH score", "Population-weighted AHEI score",
              "DALY rate (per 100,000 adults aged 60+)", "Prevalence rate (per 100,000 adults aged 60+)",
              "kg CO2-eq/person/day", "m2/person/day", "L-eq/person/day", "g PO4-eq/person/day"),
  digits = c(1L, 2L, 1L, 1L, 1L, 1L, 2L, 2L, 0L, 2L)
)

wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  weighted.mean(x[ok], w[ok])
}
summarise_geo <- function(x, geo) {
  x[, geographical_level := geo]
  x[, lapply(.SD, wmean, w = population_60plus_from_gbd),
    by = .(geographical_level, year), .SDcols = vars$variable]
}

global <- summarise_geo(dt, "Global")
region_map <- c(HIC = "High-income countries", FSU = "Central and Eastern Europe and Central Asia",
                CEECA = "Central and Eastern Europe and Central Asia", Asia = "East and Southeast Asia",
                LAC = "Latin America and the Caribbean", MENA = "Middle East and North Africa",
                SAARC = "South Asia", SSA = "Sub-Saharan Africa")
region_order <- c("Central and Eastern Europe and Central Asia",
                  "East and Southeast Asia", "High-income countries", "Latin America and the Caribbean",
                  "Middle East and North Africa", "South Asia", "Sub-Saharan Africa")
dt[, geographical_region := unname(region_map[as.character(superregion2)])]
dt[is.na(geographical_region), geographical_region := as.character(superregion2)]
dt[geographical_region == "Central/Eastern Europe and Central Asia",
   geographical_region := "Central and Eastern Europe and Central Asia"]
if (anyNA(dt$geographical_region) ||
    any(!dt$geographical_region %in% region_order)) {
  stop("Unrecognised superregion2 values: ",
       paste(unique(dt[is.na(geographical_region) |
                         !geographical_region %in% region_order, superregion2]), collapse = ", "))
}
if (!all(region_order %in% dt$geographical_region)) stop("Missing required region(s).")
regional <- rbindlist(lapply(region_order, function(r) {
  summarise_geo(dt[geographical_region == r], r)
}), fill = TRUE)
global_long <- melt(global, id.vars = c("geographical_level", "year"), variable.name = "variable", value.name = "value")
regional_long <- melt(regional, id.vars = c("geographical_level", "year"), variable.name = "variable", value.name = "value")
global_long <- merge(global_long, vars, by = "variable", all.x = TRUE)
regional_long <- merge(regional_long, vars, by = "variable", all.x = TRUE)
fwrite(global_long, file.path(table_dir, "Table_global_population_weighted_trends_v2_style_fixed_v1.csv"), bom = TRUE)
fwrite(regional_long, file.path(table_dir, "Table_regional_population_weighted_trends_v2_style_fixed_v1.csv"), bom = TRUE)

changes <- global_long[, .(
  estimate_1990 = value[year == 1990L][1L],
  estimate_2018 = value[year == 2018L][1L]
), by = .(geographical_level, variable, panel,
          label = title, unit = y_label, digits)]
regional_changes <- regional_long[, .(
  estimate_1990 = value[year == 1990L][1L],
  estimate_2018 = value[year == 2018L][1L]
), by = .(geographical_level, variable, panel,
          label = title, unit = y_label, digits)]
change_table <- rbindlist(list(changes, regional_changes), fill = TRUE)
change_table[, absolute_change := estimate_2018 - estimate_1990]
change_table[, relative_change_percent := 100 * absolute_change / abs(estimate_1990)]
fwrite(change_table, file.path(table_dir, "Supplementary_Table_S2_global_regional_changes_1990_2018_v2.csv"), bom = TRUE)

navy <- "#1F3A56"
salmon <- "#C86E5B"
grid_grey <- "#D9D9D9"

format_value <- function(x, digits) {
  scales::number(x, accuracy = 10^(-digits), big.mark = ",", trim = TRUE)
}

make_global_plot <- function(v) {
  info <- vars[variable == v]
  x <- global_long[variable == v]
  x[, value_label := format_value(value, info$digits)]
  x[, label_vjust := -0.75]
  x[variable == "dementia_prevalence_asr60plus" & year == 2015L,
    label_vjust := -1.25]
  x[variable == "dementia_prevalence_asr60plus" & year == 2018L,
    label_vjust := 1.50]
  ggplot(x, aes(x = year, y = value, group = 1L)) +
    geom_line(colour = navy, linewidth = 0.9) +
    geom_point(shape = 21, fill = "white", colour = navy, size = 2.7, stroke = 0.9) +
    geom_text(aes(label = value_label, vjust = label_vjust),
              colour = salmon, size = 3.55,
              family = "Arial", lineheight = 0.9) +
    scale_x_continuous(breaks = study_years, expand = expansion(mult = c(0.06, 0.04))) +
    scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(mult = c(0.08, 0.14))) +
    labs(title = paste0(info$panel, " ", info$title), x = "Year", y = info$y_label) +
    theme_classic(base_family = "Arial", base_size = 9.5) +
    theme(
      plot.title = element_text(size = 10.5, face = "bold", hjust = 0, colour = "black", margin = margin(b = 5)),
      axis.title.x = element_text(size = 9, colour = "black", margin = margin(t = 4)),
      axis.title.y = element_text(size = 9, colour = "black", margin = margin(r = 5)),
      axis.text = element_text(size = 8.2, colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      panel.grid.major.y = element_line(colour = grid_grey, linewidth = 0.3),
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(7, 7, 7, 7)
    )
}

# Supplementary Figure S2: global temporal trends; separate from Figure 2.
EXPORT_SUPPLEMENTARY_TRENDS <- TRUE
if (EXPORT_SUPPLEMENTARY_TRENDS) {
plots <- lapply(vars$variable, make_global_plot)
combined_figure <- wrap_plots(plots, ncol = 5, nrow = 2, guides = "keep")
ggsave(file.path(figure_dir, "Supplementary_Figure_S2_global_temporal_trends_v2.png"),
       combined_figure, width = 18, height = 9, units = "in", dpi = 600, limitsize = FALSE, bg = "white")
ggsave(file.path(figure_dir, "Supplementary_Figure_S2_global_temporal_trends_v2.pdf"),
       combined_figure, width = 18, height = 9, units = "in", limitsize = FALSE, bg = "white")

message("Completed: manuscript-style global temporal trends")
message("Input: ", input_file)
message("Figure: ", file.path(figure_dir, "Supplementary_Figure_S2_global_temporal_trends_v2.png"))
message("Tables: ", table_dir)

}

# ---- Figure 2: reference-style horizontal bars ---------------------------
FONT_FAMILY <- "Arial"
AXIS_SIZE <- 14.5
LABEL_SIZE <- 15.5
LEGEND_TEXT_SIZE <- LABEL_SIZE + 1.5
FIGURE_WIDTH <- 20
FIGURE_HEIGHT <- 22   # portrait-oriented output with equal-height panel rows
EXPORT_INDIVIDUAL_PANELS <- FALSE

location_order <- c("Global", region_order)
location_axis_labels <- c(
  "Global" = "Global",
  "Central and Eastern Europe and Central Asia" = "CEECA",
  "East and Southeast Asia" = "ESEA",
  "High-income countries" = "HIC",
  "Latin America and the Caribbean" = "LAC",
  "Middle East and North Africa" = "MENA",
  "South Asia" = "SA",
  "Sub-Saharan Africa" = "SSA"
)
location_legend_labels <- c(
  "Global" = "Global",
  "Central and Eastern Europe and Central Asia" =
    "Central/Eastern Europe and Central Asia (CEECA)",
  "East and Southeast Asia" = "East and Southeast Asia (ESEA)",
  "High-income countries" = "High-income countries (HIC)",
  "Latin America and the Caribbean" = "Latin America and the Caribbean (LAC)",
  "Middle East and North Africa" = "Middle East and North Africa (MENA)",
  "South Asia" = "South Asia (SA)",
  "Sub-Saharan Africa" = "Sub-Saharan Africa (SSA)"
)
# Muted, darker red-yellow-green palette. The sequence follows the reference
# palette while reducing brightness and saturation for publication use.
region_colors <- c(
  "Global" = "#B83A32",
  "Central and Eastern Europe and Central Asia" = "#CF5B3E",
  "East and Southeast Asia" = "#D98A45",
  "High-income countries" = "#D4B45B",
  "Latin America and the Caribbean" = "#AAA957",
  "Middle East and North Africa" = "#789858",
  "South Asia" = "#4F865D",
  "Sub-Saharan Africa" = "#276B52"
)

# Neutral styling shared by bars, legend keys, panel frames, and title strips.
bar_border_colour <- "#4A4A4A"
strip_fill_colour <- "#E3DED3"
panel_border_colour <- "#6F6B64"
# Line breaks are only in long titles, never in region names.
panel_titles <- c(
  "A. mPHDI\n(0–120)",
  "B. MED\n(0–8)",
  "C. DASH\n(8–40)",
  "D. AHEI\n(0–100)",
  "E. Dementia DALY rate\n(per 100,000 adults)",
  "F. Dementia prevalence rate\n(per 100,000 adults)",
  "G. GHG emissions\n(kg CO₂-equivalents/person/day)",
  "H. Land use\n(m²/person/day)",
  "I. Water scarcity\n(L-equivalents/person/day)",
  "J. Eutrophication\n(g PO₄-equivalents/person/day)"
)
names(panel_titles) <- LETTERS[1:10]

change_plot_data <- copy(change_table[, .(
  location = as.character(geographical_level), variable, panel,
  indicator = label, unit, estimate_1990, estimate_2018, absolute_change
)])
change_plot_data[, location := factor(location, levels = location_order)]
if (nrow(change_plot_data) != 80L || anyNA(change_plot_data$location) ||
    any(!is.finite(change_plot_data$absolute_change)) ||
    anyDuplicated(change_plot_data[, .(panel, location)])) {
  stop("Figure 2 requires exactly 10 indicators x 8 locations, with finite changes.")
}

make_change_plot <- function(panel_id, show_regions = panel_id %in% c("G", "H", "I", "J")) {
  d <- copy(change_plot_data[panel == panel_id])
  global_change <- d[location == "Global", absolute_change]
  if (length(global_change) != 1L) stop("Expected one Global value in panel ", panel_id)
  d[, strip_title := unname(panel_titles[panel_id])]
  ggplot(d, aes(x = location, y = absolute_change, fill = location)) +
    geom_col(width = 0.76, colour = bar_border_colour, linewidth = 0.28) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
    geom_hline(yintercept = global_change, colour = "grey20",
               linetype = "dashed", linewidth = 0.55) +
    scale_fill_manual(values = region_colors, limits = location_order, drop = FALSE) +
    scale_x_discrete(limits = location_order,
                     labels = location_axis_labels,
                     drop = FALSE,
                     expand = expansion(add = 0.50)) +
    scale_y_continuous(labels = scales::label_number(big.mark = ","),
                       breaks = scales::breaks_pretty(n = 4),
                       expand = expansion(mult = c(0.06, 0.08))) +
    facet_wrap(~strip_title, ncol = 1) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_family = FONT_FAMILY, base_size = AXIS_SIZE) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(size = AXIS_SIZE, face = "plain", colour = "black"),
      axis.text.x = if (show_regions)
        element_text(size = AXIS_SIZE, face = "plain", colour = "black",
                     angle = 45, hjust = 1, vjust = 1,
                     margin = margin(t = 7)) else element_blank(),
      axis.text.y = element_text(size = AXIS_SIZE, face = "plain",
                                 colour = "black", angle = 0),
      axis.title.x = element_blank(),
      axis.ticks = element_line(colour = "grey40", linewidth = 0.35),
      axis.ticks.x = if (show_regions) element_line(colour = "grey40", linewidth = 0.35)
                     else element_blank(),
      panel.grid.major = element_line(colour = "#E8E8E8", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = 0.85),
      strip.background = element_rect(fill = strip_fill_colour,
                                      colour = panel_border_colour, linewidth = 0.85),
      strip.text = element_text(size = LABEL_SIZE, face = "plain", colour = "black",
                                hjust = 0.5, vjust = 0.5, lineheight = 1.05,
                                margin = margin(7, 5, 7, 5)),
      legend.position = "none",
      plot.margin = if (show_regions) margin(10, 12, 16, 8)
                    else margin(10, 12, 8, 8)
    )
}

# The custom grid legend uses text 1.5 pt larger than the panel labels.
# Its row pitch is two character heights, adding half a character of spacing
# relative to the previous 1.5-character setting.
make_label_legend <- function() {
  row_pitch <- 2.0 * LEGEND_TEXT_SIZE
  text_gp <- grid::gpar(fontfamily = FONT_FAMILY,
                        fontsize = LEGEND_TEXT_SIZE, col = "black")
  block_height <- 32 + 7 * row_pitch + 16 + 30
  top <- grid::unit(0.5, "npc") +
    grid::unit(block_height / 2 + 12, "pt")
  # The legend spans columns 3 and 4 in the top row. One quarter of this
  # width aligns its start with the centreline of the third plot column.
  left <- grid::unit(0.25, "npc")
  children <- list(grid::textGrob("Label", x = left, y = top,
                                 just = c("left", "top"), gp = text_gp))
  for (i in seq_along(location_order)) {
    yy <- top - grid::unit(32 + (i - 1) * row_pitch, "pt")
    children[[length(children) + 1L]] <- grid::rectGrob(
      x = left + grid::unit(7, "pt"), y = yy,
      width = grid::unit(14, "pt"), height = grid::unit(14, "pt"),
      gp = grid::gpar(fill = unname(region_colors[location_order[i]]),
                      col = bar_border_colour, lwd = 0.7))
    children[[length(children) + 1L]] <- grid::textGrob(
      unname(location_legend_labels[location_order[i]]),
      x = left + grid::unit(23, "pt"), y = yy,
      just = "left", gp = text_gp)
  }
  # Explain the reference line within the legend, with no figure note.
  yy <- top - grid::unit(32 + 7 * row_pitch + 30, "pt")
  children[[length(children) + 1L]] <- grid::segmentsGrob(
    x0 = left, x1 = left + grid::unit(16, "pt"), y0 = yy, y1 = yy,
    gp = grid::gpar(col = "grey20", lty = "dashed", lwd = 1.5))
  children[[length(children) + 1L]] <- grid::textGrob(
    "Global absolute change", x = left + grid::unit(23, "pt"), y = yy,
    just = "left", gp = text_gp)
  grid::gTree(children = do.call(grid::gList, children))
}

change_plots <- lapply(LETTERS[1:10], make_change_plot)
names(change_plots) <- LETTERS[1:10]
legend_panel <- patchwork::wrap_elements(full = make_label_legend(), clip = FALSE)
figure_2 <- patchwork::wrap_plots(
  c(change_plots, list(K = legend_panel)),
  design = "ABKK\nCDEF\nGHIJ",
  widths = rep(1, 4), heights = rep(1, 3)
)
# All ten plotting panels use the same width and height.
# No figure title, subtitle, numeric annotations, or note.
save_plot_both <- function(plot_object, filename, width, height) {
  ggsave(file.path(figure_dir, paste0(filename, ".png")), plot = plot_object,
         width = width, height = height, units = "in", dpi = 600,
         limitsize = FALSE, bg = "white")
  if (!capabilities("cairo")) stop("Cairo support is required for Unicode PDF labels.")
  ggsave(file.path(figure_dir, paste0(filename, ".pdf")), plot = plot_object,
         width = width, height = height, units = "in", device = grDevices::cairo_pdf,
         limitsize = FALSE, bg = "white")
}
output_stem <- "Figure_2_v22_enlarged_raised_legend"
save_plot_both(figure_2, output_stem, FIGURE_WIDTH, FIGURE_HEIGHT)
if (EXPORT_INDIVIDUAL_PANELS) {
  for (panel_id in LETTERS[1:10]) {
    save_plot_both(make_change_plot(panel_id, show_regions = TRUE),
                   paste0("Figure_2_panel_", panel_id), width = 10, height = 5)
  }
}
setorder(change_plot_data, panel, location)
fwrite(change_plot_data, file.path(table_dir, paste0(output_stem, "_plot_data.csv")), bom = TRUE)
message("Figure 2 saved: ", file.path(figure_dir, paste0(output_stem, ".pdf")))
message("Figure 2 saved: ", file.path(figure_dir, paste0(output_stem, ".png")))
