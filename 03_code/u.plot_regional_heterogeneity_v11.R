# Figure 3 v9: top labels extend upward from ticks without covering cells.
# Each large cell is split into four dietary-pattern quadrants, clockwise:
# mPHDI (top-left), MED (top-right), DASH (bottom-right), AHEI (bottom-left).
# Tile shade: omnibus regional P interaction. Numbers/stars: subgroup P values.
# P interaction is shared across regions for each outcome-pattern combination.
# Grey macro-cell surrounds and white internal crosses; close inset legends.

pkgs <- c('data.table', 'ggplot2', 'patchwork')
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop('Please install: ', paste(missing_pkgs, collapse = ', '))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# -------------------- User settings --------------------
analysis_dir <- 'F:/A-博士阶段/GDD公开数据/revised_analysis'
result_file <- file.path(analysis_dir, '03_analysis',
  'Subgroup_dietary_patterns_by_region_v7', 'tables',
  'TABLE_subgroup_dietary_patterns_by_region_v7.csv')
figure_dir <- file.path(analysis_dir, '03_analysis',
  'Subgroup_dietary_patterns_by_region_v7', 'figures')
output_name <- 'Figure_3_heatmap_v9_top_labels_clear_panel'
font_family <- 'Arial' # Set to an installed font if Arial is unavailable.
figure_width <- 22
figure_height <- 14
jpeg_dpi <- 600

regions <- c('High-income countries',
  'Central/Eastern Europe and Central Asia', 'East and Southeast Asia',
  'Latin America and the Caribbean', 'Middle East and North Africa',
  'South Asia', 'Sub-Saharan Africa')
diets <- c('mPHDI', 'MED', 'DASH', 'AHEI')
outcomes <- c('Dementia DALY rate', 'Dementia prevalence rate',
  'GHG emissions', 'Land use', 'Water-scarcity', 'Eutrophication')

# Adjacent colour families requested by the user. First shade is significant;
# second shade is nonsignificant, based on the regional interaction P value.
dark_cols <- c(mPHDI = '#ec9565', MED = '#ffca85',
               DASH = '#c1d69c', AHEI = '#7bb384')
light_cols <- c(mPHDI = '#e6a481', MED = '#fde091',
                DASH = '#e3e7b5', AHEI = '#9ac6a0')

# -------------------- Load and validate --------------------
if (!file.exists(result_file)) stop('Result file not found: ', result_file,
  '\nRun the existing v7 subgroup analysis first, or edit result_file.')
res <- fread(result_file, encoding = 'UTF-8')
setnames(res, sub('^\ufeff', '', names(res)))
required <- c('outcome', 'dietary_pattern', 'region', 'p_value', 'p_interaction')
if (!all(required %in% names(res))) stop('Missing columns: ',
  paste(setdiff(required, names(res)), collapse = ', '))
for (v in c('outcome', 'dietary_pattern', 'region'))
  set(res, j = v, value = trimws(as.character(res[[v]])))
# Accept the equivalent full spelling of this region.
res[region == 'Central and Eastern Europe and Central Asia',
    region := 'Central/Eastern Europe and Central Asia']
res <- res[outcome %in% outcomes & dietary_pattern %in% diets]
if (any(!res$region %in% regions)) stop('Unexpected region labels: ',
  paste(unique(res$region[!res$region %in% regions]), collapse = ', '))
keys <- c('outcome', 'dietary_pattern', 'region')
if (anyDuplicated(res[, ..keys])) stop('Duplicate outcome-pattern-region rows.')
expected <- CJ(outcome = outcomes, dietary_pattern = diets, region = regions)
missing_rows <- fsetdiff(expected, res[, ..keys])
if (nrow(missing_rows)) {
  print(missing_rows)
  stop('Missing subgroup results; expected 6 x 4 x 7 = 168 rows.')
}
# Parse both numeric P values and strings such as '<0.001'.
parse_p <- function(x) {
  s <- trimws(as.character(x))
  out <- suppressWarnings(as.numeric(s))
  lt <- grepl('^<\\s*0?\\.001$', s)
  out[lt] <- 0.0005
  out[is.na(out) & !lt] <- NA_real_
  out
}
res[, p_subgroup := parse_p(p_value)]
if (anyNA(res$p_subgroup) || any(!is.finite(res$p_subgroup)) ||
    any(res$p_subgroup < 0 | res$p_subgroup > 1))
  stop('Invalid or missing subgroup P values.')
res[, p_int := parse_p(p_interaction)]
if (anyNA(res$p_int) || any(!is.finite(res$p_int)) ||
    any(res$p_int < 0 | res$p_int > 1))
  stop('Invalid or missing interaction P values.')
interaction_check <- res[, .(n_p = uniqueN(p_int)), by = .(outcome, dietary_pattern)]
if (any(interaction_check$n_p != 1L))
  stop('Expected one omnibus interaction P per outcome-pattern combination.')

# Quadrants arranged clockwise: top-left, top-right, bottom-right, bottom-left.
quadrants <- data.table(dietary_pattern = diets,
  dx = c(-0.25, 0.25, 0.25, -0.25),
  dy = c(0.25, 0.25, -0.25, -0.25))
res <- merge(res, quadrants, by = 'dietary_pattern', all.x = TRUE, sort = FALSE)
res[, x := match(region, regions) + dx]
res[, y := length(outcomes) + 1L - match(outcome, outcomes) + dy]
res[, fill_colour := ifelse(p_int < 0.05,
  unname(dark_cols[dietary_pattern]), unname(light_cols[dietary_pattern]))]
res[, star := fifelse(p_subgroup < 0.01, '**',
  fifelse(p_subgroup < 0.05, '*', ''))]
res[, p_label := ifelse(p_subgroup < 0.001, '<0.001', sprintf('%.3f', p_subgroup))]
res[, p_label := paste0(p_label, star)]
res[, dietary_pattern := factor(dietary_pattern, levels = diets)]

region_labels <- regions
macro_cells <- CJ(x = seq_along(regions), y = seq_along(outcomes))
main_plot <- ggplot() +
  # Grey surround, white inset base, and four coloured tiles.
  # This separates the grey outer frame from the white internal cross.
  geom_tile(data = macro_cells, aes(x, y), width = 1, height = 1,
    fill = '#C8C8C8', colour = NA) +
  geom_tile(data = macro_cells, aes(x, y), width = 0.978, height = 0.978,
    fill = 'white', colour = NA) +
  geom_tile(data = res, aes(x = x, y = y, fill = fill_colour),
            width = 0.466, height = 0.466, colour = NA) +
  geom_text(data = res, aes(x = x, y = y, label = p_label),
            family = font_family, size = 12.5 / (72.27 / 25.4),
            colour = '#303030', fontface = 'plain') +
  scale_fill_identity() +
  scale_x_continuous(breaks = seq_along(regions), labels = region_labels,
    position = 'top',
    expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(outcomes), labels = rev(outcomes),
    expand = c(0, 0)) +
  coord_fixed(ratio = 1, xlim = c(0.5, 7.5), ylim = c(0.5, 6.5),
              clip = 'off') +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = font_family, base_size = 14) +
  theme(
    text = element_text(size = 14, colour = 'black', face = 'plain'),
    panel.background = element_rect(fill = 'white', colour = NA),
    plot.background = element_rect(fill = 'white', colour = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.text.x.top = element_text(size = 16.5, colour = 'black', face = 'plain',
                                   angle = 35, hjust = 0, vjust = 0,
                                   margin = margin(b = 10)),
    axis.ticks.x = element_line(colour = 'black', linewidth = 0.45),
    axis.ticks.x.top = element_line(colour = 'black', linewidth = 0.45),
    axis.text.y = element_text(size = 16.5, colour = 'black', face = 'plain',
                               margin = margin(r = 8)),
    axis.ticks = element_line(colour = 'black', linewidth = 0.45),
    axis.ticks.length = grid::unit(0.10, 'cm'),
    axis.title = element_blank(),
    legend.position = 'none',
    plot.margin = margin(28, 165, 10, 12)
  )

# Right-side note: the same clockwise positions as the main plot.
# >= includes the boundary value 0.05 in the nonsignificant group.
legend_centres <- c(3.9, 2.0)
legend_data <- rbindlist(lapply(seq_along(legend_centres), function(i) {
  data.table(dietary_pattern = diets,
    x = 1 + quadrants$dx,
    y = legend_centres[i] + quadrants$dy,
    fill_colour = unname(if (i == 1) light_cols[diets] else dark_cols[diets]))
}))
legend_boxes <- data.table(y = legend_centres)
legend_plot <- ggplot() +
  geom_tile(data = legend_boxes, aes(x = 1, y = y),
    width = 1, height = 1, fill = '#C8C8C8', colour = NA) +
  geom_tile(data = legend_boxes, aes(x = 1, y = y),
    width = 0.978, height = 0.978, fill = 'white', colour = NA) +
  geom_tile(data = legend_data, aes(x, y, fill = fill_colour),
    width = 0.466, height = 0.466, colour = NA) +
  geom_text(data = legend_data, aes(x, y, label = dietary_pattern),
    family = font_family, fontface = 'plain', colour = '#303030',
    size = 11 / (72.27 / 25.4)) +
  annotate('text', x = 1, y = 5.05, label = 'Note:',
    family = font_family, size = 12 / (72.27 / 25.4)) +
  annotate('text', x = 1, y = 4.65, label = 'P interaction >= 0.05',
    family = font_family, size = 11 / (72.27 / 25.4)) +
  annotate('text', x = 1, y = 2.75, label = 'P interaction < 0.05',
    family = font_family, size = 11 / (72.27 / 25.4)) +
  scale_fill_identity() +
  coord_fixed(xlim = c(0.35, 1.65), ylim = c(1.2, 5.4), clip = 'off') +
  theme_void(base_family = font_family) +
  theme(legend.position = 'none', plot.margin = margin(0, 0, 0, 0))

# Position relative to the heatmap panel to avoid an extra layout column gap.
fig <- main_plot + inset_element(legend_plot,
  left = 1.008, right = 1.195, bottom = 0.15, top = 0.85,
  align_to = 'panel', clip = FALSE)

# -------------------- Export --------------------
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figure_dir, paste0(output_name, '.jpeg')), fig,
       width = figure_width, height = figure_height, units = 'in',
       device = 'jpeg', dpi = jpeg_dpi, quality = 100,
       bg = 'white', limitsize = FALSE)
pdf_device <- if (capabilities('cairo')) grDevices::cairo_pdf else grDevices::pdf
ggsave(file.path(figure_dir, paste0(output_name, '.pdf')), fig,
       width = figure_width, height = figure_height, units = 'in',
       device = pdf_device, bg = 'white', limitsize = FALSE)
message('Saved: ', file.path(figure_dir, paste0(output_name, '.jpeg')))
message('Saved: ', file.path(figure_dir, paste0(output_name, '.pdf')))
