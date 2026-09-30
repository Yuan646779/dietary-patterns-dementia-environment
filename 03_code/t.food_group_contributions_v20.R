# Figure 4: transposed rank heatmap for 15 food groups and 6 outcomes.
# Plotting only; no statistical model is refitted.
#
# Rows    = outcomes.
# Columns = food groups.
# Rank 1  = highest signed health estimate or highest environmental contribution.
# Colours represent the original estimates, not the ranks:
#   * dementia outcomes: two negative and two positive intervals, split at zero;
#   * environmental outcomes: four intervals centred on the rank-8 value;
#     each side is divided at its midpoint, and rank 8 enters the third colour.
# Ranks are descriptive and should not be interpreted as causal harm rankings.
# No figure title or note is drawn. Requires only base R and grid.

project_dir <- "F:/A-博士阶段/GDD公开数据/revised_analysis"

# Optional: enter the full path to the completed 90-row results CSV here.
input_file <- ""

output_dir <- file.path(
  project_dir,
  "03_analysis",
  "food_group_rank_heatmap_v23"
)

if (!nzchar(input_file)) {
  candidates <- c(
    file.path(
      project_dir,
      "03_analysis",
      "food_group_contributions_v11",
      "tables",
      "TABLE_S7_food_group_contributions_all_6_outcomes_v11.csv"
    ),
    file.path(getwd(), "TABLE_S7_food_group_contributions_all_6_outcomes_v11.csv"),
    file.path(getwd(), "06df202d-98d4-4b78-ae72-762a376e3e66.csv")
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits)) {
    input_file <- hits[1L]
  } else if (interactive()) {
    input_file <- file.choose()
  } else {
    stop("Set input_file to the completed 90-row six-outcome CSV.")
  }
}

if (!file.exists(input_file)) stop("Input not found: ", input_file)

food_codes <- c(sprintf("v%02d", 1:14), "v57")
food_names <- c(
  "Fruits",
  "Non-starchy vegetables",
  "Potatoes",
  "Other starchy vegetables",
  "Beans and legumes",
  "Nuts and seeds",
  "Refined grains",
  "Whole grains",
  "Total processed meats",
  "Unprocessed red meats",
  "Total seafoods",
  "Eggs",
  "Cheese",
  "Yoghurt",
  "Total milk"
)

outcome_codes <- c(
  "dementia_daly",
  "dementia_prevalence",
  "ghg",
  "land",
  "water",
  "eutrophication"
)

outcome_names <- c(
  "Dementia DALY rate",
  "Dementia prevalence rate",
  "GHG emissions",
  "Land use",
  "Water-scarcity",
  "Eutrophication"
)

# The first two outcome rows use the supplied yellow and orange colour cards.
# Four well-separated levels (10, 40, 70, and 100) are retained from each
# 10-level card. The remaining environmental rows keep the prior reference
# palettes. Colours run from lower to higher estimates within each row.
palettes <- rbind(
  c("#f8e5db", "#FFCa85", "#e6a481", "#ec9565"),
  c("#f8e5db", "#FFCa85", "#e6a481", "#ec9565"),
  c("#fff9e9", "#f6f2cb", "#fee9b2", "#fde091"),
  c("#fff9e9", "#f6f2cb", "#fee9b2", "#fde091"),
  c("#e3e7b5", "#c1d69c", "#629f6a", "#7bb384"),
  c("#e3e7b5", "#c1d69c", "#629f6a", "#7bb384")
)
rownames(palettes) <- outcome_codes

d <- read.csv(
  input_file,
  fileEncoding = "UTF-8-BOM",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required <- c("variable_code", "indicator_code", "estimate_display")
if (!all(required %in% names(d))) {
  stop("Missing columns: ", paste(setdiff(required, names(d)), collapse = ", "))
}

d <- d[
  d$variable_code %in% food_codes & d$indicator_code %in% outcome_codes,
  ,
  drop = FALSE
]

if (nrow(d) != 90L) stop("Expected exactly 90 rows (15 foods x 6 outcomes).")

keys <- paste(d$variable_code, d$indicator_code, sep = "|")
if (anyDuplicated(keys)) stop("Duplicate food-outcome pairs in input.")

if ("fit_status" %in% names(d) &&
    any(is.na(d$fit_status) | d$fit_status != "OK")) {
  stop("At least one input result has a missing or unsuccessful fit status.")
}

d$estimate_display <- suppressWarnings(as.numeric(d$estimate_display))
if (any(!is.finite(d$estimate_display))) {
  stop("Missing or nonfinite estimates; ranking stopped.")
}

values <- ranks <- bins <- matrix(
  NA_real_,
  nrow = length(outcome_codes),
  ncol = length(food_codes),
  dimnames = list(outcome_names, food_names)
)

breaks <- vector("list", length(outcome_codes))
legend_labels <- vector("list", length(outcome_codes))
audit <- vector("list", length(outcome_codes))

format_breaks <- function(x, initial_digits) {
  digits <- initial_digits
  repeat {
    formatted <- formatC(x, format = "f", digits = digits, big.mark = ",")
    if (!anyDuplicated(formatted) || digits >= 10L) return(formatted)
    digits <- digits + 1L
  }
}

for (j in seq_along(outcome_codes)) {
  z <- d[d$indicator_code == outcome_codes[j], , drop = FALSE]
  z <- z[match(food_codes, z$variable_code), , drop = FALSE]
  if (anyNA(z$variable_code)) stop("Incomplete food groups for ", outcome_codes[j])

  v <- z$estimate_display
  if (j > 2L && any(v < 0)) {
    stop("Negative environmental contribution detected for ", outcome_names[j], ".")
  }

  values[j, ] <- v

  # Signed descending ordering. Exact ties follow the fixed food-group order.
  if (anyDuplicated(v)) {
    warning(
      "Exact ties for ", outcome_names[j],
      "; fixed food order was used to assign distinct ranks 1-15."
    )
  }
  ranks[j, ] <- rank(-v, ties.method = "first")

  lo <- min(v)
  hi <- max(v)
  if (lo == hi) stop("Constant outcome cannot be binned: ", outcome_names[j])

  if (j <= 2L) {
    if (!(lo < 0 && hi > 0)) {
      stop("Dementia outcomes must contain both negative and positive estimates.")
    }

    # Two equally spaced negative intervals and two equally spaced positive
    # intervals. Zero is the shared boundary and is assigned to the second bin.
    negative_breaks <- seq(lo, 0, length.out = 3L)
    positive_breaks <- seq(0, hi, length.out = 3L)
    b <- c(negative_breaks, positive_breaks[-1L])

    current_bins <- integer(length(v))
    current_bins[v < 0] <- as.integer(cut(
      v[v < 0],
      breaks = negative_breaks,
      include.lowest = TRUE,
      right = TRUE
    ))
    current_bins[v == 0] <- 2L
    current_bins[v > 0] <- 2L + as.integer(cut(
      v[v > 0],
      breaks = positive_breaks,
      include.lowest = FALSE,
      right = TRUE
    ))
  } else {
    # Use the estimate ranked eighth as the central boundary. Divide the lower
    # and upper ranges at their respective midpoints to obtain four colours.
    # Values equal to the rank-8 estimate are assigned to the third colour.
    rank8_value <- v[which(ranks[j, ] == 8L)]
    lower_midpoint <- mean(c(lo, rank8_value))
    upper_midpoint <- mean(c(rank8_value, hi))
    b <- c(lo, lower_midpoint, rank8_value, upper_midpoint, hi)

    if (anyDuplicated(b)) {
      stop(
        "Environmental colour boundaries are not unique for ",
        outcome_names[j], "."
      )
    }

    current_bins <- ifelse(
      v < lower_midpoint,
      1L,
      ifelse(v < rank8_value, 2L, ifelse(v < upper_midpoint, 3L, 4L))
    )
  }

  if (anyNA(current_bins) || any(!current_bins %in% 1:4)) {
    stop("Binning failed for ", outcome_names[j], ".")
  }

  breaks[[j]] <- b
  bins[j, ] <- current_bins

  initial_digits <- if (outcome_codes[j] == "water") 1L else 2L
  f <- format_breaks(b, initial_digits)
  if (j <= 2L) {
    legend_labels[[j]] <- vapply(
      1:4,
      function(k) {
        paste0(if (k == 1L) "[" else "(", f[k], ", ", f[k + 1L], "]")
      },
      character(1)
    )
  } else {
    # The rank-8 boundary is left-closed in the third interval.
    legend_labels[[j]] <- c(
      paste0("[", f[1], ", ", f[2], ")"),
      paste0("[", f[2], ", ", f[3], ")"),
      paste0("[", f[3], ", ", f[4], ")"),
      paste0("[", f[4], ", ", f[5], "]")
    )
  }

  z$food_label <- food_names
  z$outcome_label <- outcome_names[j]
  z$rank_within_outcome <- as.integer(ranks[j, ])
  z$exact_tie <- duplicated(v) | duplicated(v, fromLast = TRUE)
  z$colour_bin <- current_bins
  z$bin_lower <- b[current_bins]
  z$bin_upper <- b[current_bins + 1L]
  z$colour_hex <- palettes[j, current_bins]
  audit[[j]] <- z

  stopifnot(
    identical(sort(as.integer(ranks[j, ])), 1:15),
    !anyNA(bins[j, ])
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  do.call(rbind, audit),
  file.path(output_dir, "Figure_4_v21_rank_and_bin_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

library(grid)

# Transposed matrix: 15 columns x 6 rectangular rows.
# All matrix labels and cell ranks are 14-point type.
# Both cell dimensions are 25% larger than in version 17.
cell_width <- 0.55 * 1.25
cell_height <- 0.50 * 1.25
matrix_width <- length(food_names) * cell_width
matrix_height <- length(outcome_names) * cell_height
matrix_left <- 3.25
# Move the heatmap down slightly so the rotated food-group labels are not clipped.
matrix_bottom <- 3.15

figure_width <- 14.0
figure_height <- 8.8

# Compact legend: six rows, four intervals per row, 11-point type.
# Legend labels begin at the left edge of the third heatmap cell. Subsequent
# positions are calculated from the rendered text widths to avoid overlap.
legend_label_x <- matrix_left + 2 * cell_width
legend_row_height <- 0.30
legend_top <- 2.85
legend_font_size <- 11
legend_swatch_width <- 0.13
legend_swatch_height <- 0.13
legend_label_gap <- 0.14
legend_swatch_text_gap <- 0.05
legend_entry_gap <- 0.14

draw_figure <- function() {
  grid.newpage()

  txt <- function(
    label,
    x,
    y,
    size = 14,
    just = "centre",
    rot = 0,
    col = "#111111",
    fontface = "plain"
  ) {
    grid.text(
      label,
      x = unit(x, "in"),
      y = unit(y, "in"),
      just = just,
      rot = rot,
      gp = gpar(
        fontfamily = "sans",
        fontface = fontface,
        fontsize = size,
        col = col
      )
    )
  }

  # Measure text in inches using the same font and size as the legend.
  text_width_in <- function(label, size = legend_font_size) {
    convertWidth(
      grobWidth(textGrob(
        label,
        gp = gpar(fontfamily = "sans", fontsize = size)
      )),
      unitTo = "in",
      valueOnly = TRUE
    )
  }

  # Cells and ranks.
  for (j in seq_along(outcome_names)) {
    for (i in seq_along(food_names)) {
      x <- matrix_left + (i - 0.5) * cell_width
      y <- matrix_bottom + (length(outcome_names) - j + 0.5) * cell_height
      colour <- palettes[j, bins[j, i]]

      grid.rect(
        x = unit(x, "in"),
        y = unit(y, "in"),
        width = unit(cell_width, "in"),
        height = unit(cell_height, "in"),
        gp = gpar(fill = colour, col = "#222222", lwd = 0.85)
      )

      rgb <- col2rgb(colour) / 255
      luminance <- sum(c(0.2126, 0.7152, 0.0722) * rgb)
      txt(
        as.integer(ranks[j, i]),
        x,
        y,
        size = 14,
        col = if (luminance < 0.47) "white" else "#111111"
      )
    }
  }

  # Outcome labels on the y axis.
  for (j in seq_along(outcome_names)) {
    y <- matrix_bottom + (length(outcome_names) - j + 0.5) * cell_height
    txt(outcome_names[j], matrix_left - 0.14, y, size = 14, just = "right")
  }

  # Food-group labels on the upper x axis, rotated 45 degrees.
  matrix_top <- matrix_bottom + matrix_height
  for (i in seq_along(food_names)) {
    x <- matrix_left + (i - 0.5) * cell_width
    grid.lines(
      x = unit(c(x, x), "in"),
      y = unit(c(matrix_top, matrix_top + 0.08), "in"),
      gp = gpar(lwd = 0.8, col = "#222222")
    )
    txt(
      food_names[i],
      x,
      matrix_top + 0.11,
      size = 14,
      just = "left",
      rot = 45
    )
  }

  # Slightly stronger separator between health and environmental outcomes.
  separator_y <- matrix_bottom + 4 * cell_height
  grid.lines(
    x = unit(c(matrix_left, matrix_left + matrix_width), "in"),
    y = unit(c(separator_y, separator_y), "in"),
    gp = gpar(col = "#111111", lwd = 1.3)
  )

  # Six compact legend rows. Column positions use the widest text in each
  # interval column. The 0.14-inch gap is approximately two English characters.
  legend_label_width <- max(vapply(
    outcome_names,
    text_width_in,
    numeric(1),
    size = legend_font_size
  ))
  legend_interval_widths <- vapply(
    1:4,
    function(k) max(vapply(
      legend_labels,
      function(x) text_width_in(x[k], legend_font_size),
      numeric(1)
    )),
    numeric(1)
  )

  legend_swatch_x <- numeric(4)
  legend_swatch_x[1] <- legend_label_x + legend_label_width +
    legend_label_gap + legend_swatch_width / 2
  for (k in 2:4) {
    previous_text_start <- legend_swatch_x[k - 1] +
      legend_swatch_width / 2 + legend_swatch_text_gap
    legend_swatch_x[k] <- previous_text_start +
      legend_interval_widths[k - 1] + legend_entry_gap +
      legend_swatch_width / 2
  }

  for (j in seq_along(outcome_names)) {
    y <- legend_top - (j - 1) * legend_row_height
    txt(
      outcome_names[j],
      legend_label_x,
      y,
      size = legend_font_size,
      just = "left"
    )

    for (k in 1:4) {
      x <- legend_swatch_x[k]
      grid.rect(
        x = unit(x, "in"),
        y = unit(y, "in"),
        width = unit(legend_swatch_width, "in"),
        height = unit(legend_swatch_height, "in"),
        gp = gpar(fill = palettes[j, k], col = "#333333", lwd = 0.6)
      )
      txt(
        legend_labels[[j]][k],
        x + legend_swatch_width / 2 + legend_swatch_text_gap,
        y,
        size = legend_font_size,
        just = "left"
      )
    }
  }

  invisible(NULL)
}

base_name <- file.path(
  output_dir,
  "Figure_4_v21_yellow_orange_health_rows"
)

png(
  paste0(base_name, ".png"),
  width = figure_width,
  height = figure_height,
  units = "in",
  res = 600,
  bg = "white"
)
tryCatch(draw_figure(), finally = dev.off())

if (capabilities("cairo")) {
  cairo_pdf(
    paste0(base_name, ".pdf"),
    width = figure_width,
    height = figure_height
  )
} else {
  pdf(
    paste0(base_name, ".pdf"),
    width = figure_width,
    height = figure_height,
    useDingbats = FALSE
  )
}
tryCatch(draw_figure(), finally = dev.off())

message("Saved PNG, PDF, and rank/bin audit to: ", output_dir)
