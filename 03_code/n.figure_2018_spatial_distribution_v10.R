# 2018 spatial distributions for the final v2 master database.
# Consolidates the former V12/V14 plus direct-ASR60plus update map scripts using the current field names.

required_packages <- c("data.table", "ggplot2", "sf", "rnaturalearth", "rnaturalearthdata", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required package(s): ", paste(missing_packages, collapse = ", "))
library(data.table); library(ggplot2); library(sf); library(rnaturalearth); library(rnaturalearthdata); library(patchwork); library(scales)

analysis_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
input_file <- file.path(analysis_dir, "02_processed", "Master", "GDD_final_master_database_185_country_year_v2.csv")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}
output_dir <- file.path(analysis_dir, "03_analysis", "Figures", "Figure_2018_spatial_distribution")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dt <- fread(input_file, encoding = "UTF-8")
dt <- dt[year == 2018L]
if (nrow(dt) != 185L) stop("The 2018 panel must contain 185 countries.")

map_vars <- data.table(
  variable = c("diet_mphdi_score", "diet_med_score", "diet_dash_score", "diet_ahei_score",
               "dementia_daly_asr60plus", "dementia_prevalence_asr60plus",
               "env_ghg_kg_co2e_person_day", "env_land_m2_person_day",
               "env_water_scarcity_l_equivalent_person_day", "env_eutrophication_g_po4e_person_day"),
  panel = LETTERS[1:10],
  label = c("mPHDI", "MED", "DASH", "AHEI", "Dementia DALY rate", "Dementia prevalence rate",
            "GHG emissions", "Land use", "Scarcity-weighted water use", "Eutrophication"),
  unit = c("score", "score", "score", "score", "per 100,000", "per 100,000",
           "kg CO2-eq/person/day", "m2/person/day", "L-eq/person/day", "g PO4-eq/person/day")
)
world <- ne_countries(scale = "small", type = "countries", returnclass = "sf")
world <- world[world$continent != "Antarctica", ]
world$iso3 <- toupper(as.character(world$adm0_a3))
if (any(is.na(world$iso3) | world$iso3 %in% c("-99", ""))) {
  fallback <- intersect(c("iso_a3", "gu_a3", "sov_a3"), names(world))
  if (length(fallback)) {
    bad <- is.na(world$iso3) | world$iso3 %in% c("-99", "")
    world$iso3[bad] <- toupper(as.character(world[[fallback[1L]]][bad]))
  }
}
world <- world[, c("iso3", "geometry")]

# -----------------------------------------------------------------------------
# Reference-style appearance ONLY. Input, 2018 filter, variables, ISO matching,
# mapping QC, CSV tables and data paths are unchanged.
# Two columns x five rows; headings and colour bars are outside the map frames.
# New figure filenames use _v10_rectangular_orange_purple_blue; CSV filenames are unchanged.
# No Antarctica is added: the original geography and mapping QC are preserved.
# -----------------------------------------------------------------------------
FONT_FAMILY <- if (capabilities("cairo")) "Arial" else "sans"
SIZE_TITLE <- 14  # pt: panel letters and full titles, including units
SIZE_DETAIL <- 12 # pt: colour-bar numbers; no third font size
FIG_WIDTH <- 14
FIG_HEIGHT <- 19.5 # inches; space for external headings and larger type
EXPORT_DPI <- 600
COUNTRY_LINE_MM <- 0.16
FRAME_LINE_WIDTH <- 0.65
NA_COLOUR <- "#E6E6E6"

# Reference palettes selected from the supplied figures.
# A–D: upper orange sequence from figure 1; E–F: middle purple sequence from
# figure 2; G–J: blue sequence from figure 1.
# Darker colour means a larger value; each variable keeps its own scale.
diet_palette <- c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59",
                  "#EF6548", "#D7301F", "#B30000", "#7F0000")
dementia_palette <- c("#F7F7FB", "#ECE7F2", "#D0D1E6", "#A6BDDB", "#9E9AC8",
                      "#807DBA", "#6A51A3", "#54278F", "#3F007D")
environment_palette <- c("#F7FCF0", "#E0F3DB", "#CCEBC5", "#A8DDB5", "#7BCCC4",
                         "#4EB3D3", "#2B8CBE", "#0868AC", "#084081")

# Figure-only headings: preserve the input data and original CSV tables.
# Ranges below are the requested index labels, not colour-scale limits.
# All titles are retained on one line at the same title font size.
display_titles <- setNames(c(
  "mPHDI (0\u2013120)",
  "MED (0\u20138)",
  "DASH (8\u201340)",
  "AHEI (0\u2013100)",
  "Dementia DALY rate (per 100,000 adults)",
  "Dementia prevalence rate (per 100,000 adults)",
  "GHG emissions (kg CO\u2082-equivalents/person/day)",
  "Land use (m\u00b2/person/day)",
  "Water scarcity (L-equivalents/person/day)",
  "Eutrophication (g PO\u2084-equivalents/person/day)"
), map_vars$variable)
if (!capabilities("cairo")) {
  # Standard PDF devices may not support these Unicode characters.
  for (pair in list(c("\u2013", "-"), c("\u2265", ">="),
                    c("\u2082", "2"), c("\u2084", "4"), c("\u00b2", "2"))) {
    display_titles <- gsub(pair[1], pair[2], display_titles, fixed = TRUE)
  }
  names(display_titles) <- map_vars$variable
}

theme_map_manuscript <- theme_void(base_family = FONT_FAMILY,
                                    base_size = SIZE_DETAIL) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.25),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

# Draw a continuous bar just below the map frame, with no unit/title above it.
# This avoids version-dependent inset-legend positioning and overlapping maps.
# Map and bar use the SAME palette, limits, Lab interpolation and transform.
# See https://scales.r-lib.org/reference/pal_gradient_n.html
colourbar_grob <- function(palette, limits, breaks, labels,
                           log_scale = FALSE) {
  tf <- if (log_scale) log10 else identity
  tick_fraction <- (tf(breaks) - tf(limits[1])) / diff(tf(limits))
  bar_left <- 0.20
  bar_width <- 0.60
  bar_bottom <- 0.105
  bar_height <- 0.040
  tick_x <- bar_left + bar_width * tick_fraction
  bar_colours <- scales::gradient_n_pal(palette, space = "Lab")(
    seq(0, 1, length.out = 1024)
  )
  grid::grobTree(
    grid::rasterGrob(as.raster(matrix(bar_colours, nrow = 1)),
      x = bar_left, y = bar_bottom, width = bar_width, height = bar_height,
      just = c("left", "bottom"), interpolate = TRUE),
    grid::rectGrob(x = bar_left, y = bar_bottom,
      width = bar_width, height = bar_height, just = c("left", "bottom"),
      gp = grid::gpar(fill = NA, col = "black", lwd = FRAME_LINE_WIDTH)),
    grid::segmentsGrob(x0 = tick_x, x1 = tick_x,
      y0 = bar_bottom, y1 = bar_bottom - 0.017,
      gp = grid::gpar(col = "black", lwd = FRAME_LINE_WIDTH)),
    grid::textGrob(labels, x = tick_x, y = 0.060,
      gp = grid::gpar(fontfamily = FONT_FAMILY, fontsize = SIZE_DETAIL,
                      col = "black", lineheight = 1))
  )
}

map_one <- function(v) {
  info <- map_vars[variable == v]
  x <- merge(world, dt[, .(iso3, value = get(v))], by = "iso3", all.x = TRUE)
  is_diet <- v %in% map_vars$variable[1:4]
  is_dementia <- v %in% map_vars$variable[5:6]
  is_environment <- !is_diet && !is_dementia
  palette <- if (is_diet) diet_palette else if (is_dementia) {
    dementia_palette
  } else environment_palette

  # Retain original scaling: dietary/dementia = linear; environment = log10.
  valid <- x$value[is.finite(x$value)]
  if (!length(valid)) stop("No finite mapped values for: ", v)
  if (is_environment && any(valid <= 0)) {
    stop("The original log10 colour scale requires positive values: ", v)
  }
  limits <- range(valid)
  if (diff(limits) == 0) {
    limits <- if (is_environment) limits * c(1 / 1.01, 1.01) else {
      limits + c(-1, 1) * max(abs(limits[1]) * 0.01, 0.1)
    }
  }
  breaks <- if (is_environment) scales::breaks_log(n = 4)(limits) else {
    scales::breaks_extended(n = 5)(limits)
  }
  breaks <- sort(unique(breaks[is.finite(breaks) &
                               breaks >= limits[1] & breaks <= limits[2]]))
  if (length(breaks) < 2L) breaks <- limits
  if (length(breaks) > 5L) {
    breaks <- breaks[unique(round(seq(1, length(breaks), length.out = 5)))]
  }
  # Keep labels readable, with thousands separators and no unnecessary .0.
  step <- min(diff(breaks))
  accuracy <- if (step >= 1) 1 else 10^floor(log10(step))
  labels <- scales::label_number(accuracy = accuracy, big.mark = ",",
                                 decimal.mark = ".", trim = TRUE)(breaks)

  p <- ggplot(x) +
    geom_sf(aes(fill = value), colour = "black", linewidth = COUNTRY_LINE_MM) +
    # Rectangular equirectangular (Plate Carree) view approximating reference 1.
    # The original screenshot does not establish an exact CRS.
    # Limits are in degrees via default_crs, without latitude-axis stretching.
    coord_sf(crs = "+proj=eqc +lat_ts=0 +lon_0=0 +datum=WGS84 +units=m +no_defs",
             default_crs = sf::st_crs(4326),
             xlim = c(-180, 180), ylim = c(-60, 85),
             datum = NA, expand = FALSE, clip = "on") +
    scale_fill_gradientn(colours = palette, space = "Lab", limits = limits,
      na.value = NA_COLOUR,
      trans = if (is_environment) "log10" else "identity",
      guide = "none") + theme_map_manuscript

  # Fixed map/heading/bar strips. Only the geographic panel has a frame.
  map_grob <- grid::grobTree(ggplotGrob(p),
    vp = grid::viewport(x = 0.5, y = 0.525, width = 0.95, height = 0.70))
  title_text <- unname(display_titles[v])
  two_lines <- grepl("\n", title_text, fixed = TRUE)
  panel_grob <- grid::grobTree(
    grid::rectGrob(gp = grid::gpar(fill = "white", col = NA)),
    map_grob,
    grid::textGrob(paste0(info$panel, "."), x = 0.025,
      y = if (two_lines) 0.966 else 0.940, just = "left",
      gp = grid::gpar(fontfamily = FONT_FAMILY, fontsize = SIZE_TITLE,
                      fontface = "bold", col = "black", lineheight = 1)),
    grid::textGrob(title_text, x = 0.078, y = 0.940, just = "left",
      gp = grid::gpar(fontfamily = FONT_FAMILY, fontsize = SIZE_TITLE,
                      fontface = "plain", col = "black", lineheight = 1.05)),
    colourbar_grob(palette, limits, breaks, labels, log_scale = is_environment)
  )
  ggsave(file.path(output_dir, paste0("Figure_1_2018_spatial_panel_",
    info$panel, "_", v, "_v10_rectangular_orange_purple_blue.png")), panel_grob,
    width = FIG_WIDTH / 2, height = FIG_HEIGHT / 5, units = "in",
    dpi = EXPORT_DPI, bg = "white")
  panel_grob
}

plots <- lapply(map_vars$variable, map_one)
# Row-major sequence: A/B, C/D, E/F, G/H, I/J.
# Grid retains fixed font sizes and identical map/legend spacing in every cell.
placed_panels <- lapply(seq_along(plots), function(i) {
  grid::grobTree(plots[[i]], vp = grid::viewport(
    layout.pos.row = (i - 1L) %/% 2L + 1L,
    layout.pos.col = (i - 1L) %% 2L + 1L))
})
figure <- grid::gTree(children = do.call(grid::gList, placed_panels),
  vp = grid::viewport(layout = grid::grid.layout(nrow = 5, ncol = 2)))

ggsave(file.path(output_dir, "Figure_1_2018_ten_panel_spatial_distribution_v10_rectangular_orange_purple_blue.png"),
  figure, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in",
  dpi = EXPORT_DPI, limitsize = FALSE, bg = "white")
# Cairo supports Arial and Unicode unit labels in the PDF on typical Windows R.
# On builds without Cairo, use the generic sans family and standard PDF device.
ggsave(file.path(output_dir, "Figure_1_2018_ten_panel_spatial_distribution_v10_rectangular_orange_purple_blue.pdf"),
  figure, device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf,
  width = FIG_WIDTH, height = FIG_HEIGHT, units = "in",
  limitsize = FALSE, bg = "white")

mapping_qc <- data.table(variable = map_vars$variable, label = map_vars$label,
                         n_2018_records = vapply(map_vars$variable, function(v) sum(!is.na(dt[[v]])), integer(1)),
                         n_mapped_polygons = vapply(map_vars$variable, function(v) sum(!is.na(merge(world, dt[, .(iso3, value = get(v))], by = "iso3", all.x = TRUE)$value)), integer(1)))
fwrite(mapping_qc, file.path(output_dir, "QC_2018_spatial_mapping_summary_v4.csv"), bom = TRUE)
fwrite(map_vars, file.path(output_dir, "Table_2018_spatial_variable_summary_v4.csv"), bom = TRUE)
message("Completed: 2018 spatial distribution maps")
message("Output directory: ", output_dir)
