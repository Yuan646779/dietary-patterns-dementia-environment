# ==============================================================================
# r.sensitivity_loco_nonlinear_interpolation_v3.R
# Leave-one-country-out, nonlinear, and interpolation-based exploratory analyses
# for GDD_final_master_database_185_country_year_v2.csv
# ==============================================================================
install.packages("ggplot2")
library(ggplot2)
pkgs <- c("data.table", "fixest", "ggplot2")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
library(data.table)
library(fixest)

analysis_root <- "F:/A-博士阶段/GDD公开数据/revised_analysis"
input_file <- file.path(analysis_root, "02_processed", "Master",
                         "GDD_final_master_database_185_country_year_v2.csv")
output_root <- file.path(analysis_root, "03_analysis",
                         "r_sensitivity_loco_nonlinear_interpolation_v3")
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")
data_dir <- file.path(output_root, "analysis_data")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

fixest::setFixest_nthreads(max(1L, parallel::detectCores(logical = FALSE) - 1L))
ssc_setting <- fixest::ssc(adj = TRUE, cluster.adj = TRUE,
                           cluster.df = "min", t.df = "min",
                           fixef.K = "nested")
study_years <- c(1990L, 1995L, 2000L, 2005L, 2010L, 2015L, 2018L)
health_lags <- c(0L, 5L, 10L, 15L)
primary_health_lag <- 10L

dt <- fread(input_file, na.strings = c("", "NA", "NaN"))
dt[, iso3 := toupper(trimws(as.character(iso3)))]
dt[, year := as.integer(year)]
setorder(dt, iso3, year)

resolve <- function(x, label, required = TRUE) {
  z <- x[x %in% names(dt)]
  if (length(z)) return(z[1L])
  if (required) stop("Cannot resolve ", label, ". Tried: ",
                     paste(x, collapse = ", "))
  NA_character_
}

diet <- data.table(
  diet = c("mPHDI", "MED", "DASH", "AHEI"),
  increment = c(10, 1, 5, 10),
  increment_label = c("per 10-point higher mPHDI",
                      "per 1-point higher MED",
                      "per 5-point higher DASH",
                      "per 10-point higher AHEI"),
  score = c(
    resolve(c("diet_mPHDI_score", "diet_mphdi_score", "mphdi_total", "mPHDI"), "mPHDI"),
    resolve(c("diet_MED_score", "diet_med_score", "med_total", "MED"), "MED"),
    resolve(c("diet_DASH_score", "diet_dash_score", "dash_total", "DASH"), "DASH"),
    resolve(c("diet_AHEI_score", "diet_ahei_score", "ahei_total", "AHEI"), "AHEI")
  )
)
sdi_col <- resolve(c("sdi", "SDI"), "SDI")
daly_col <- resolve(c("dementia_daly_asr60plus", "dementia_daly_rate_60plus"),
                    "dementia DALY rate")
prev_col <- resolve(c("dementia_prevalence_asr60plus",
                      "dementia_prevalence_rate_60plus"),
                    "dementia prevalence rate")
coverage_col <- resolve(c("included_mass_intake_coverage",
                          "env_coverage_proportion",
                          "environmental_intake_coverage",
                          "food_intake_coverage"),
                        "environmental coverage", required = FALSE)

env_raw <- c(
  ghg = resolve(c("ln_ghg", "env_ghg_kg_co2e_person_day",
                  "ghg_emissions", "GHG_emissions", "ghg"), "GHG"),
  land = resolve(c("ln_land", "env_land_m2_person_day",
                   "land_use", "Land_use", "land"), "land use"),
  water = resolve(c("ln_scarcity_water",
                    "env_water_scarcity_l_equivalent_person_day",
                    "scarcity_weighted_water_use", "water_scarcity",
                    "Water_scarcity", "water"),
                  "water scarcity"),
  eutrophication = resolve(c("ln_eutrophication",
                             "env_eutrophication_g_po4e_person_day",
                             "eutrophication_potential", "Eutrophication",
                             "eutrophication"),
                            "eutrophication")
)

check_cols <- c("iso3", "year", sdi_col, daly_col, prev_col,
                diet$score, unname(env_raw))
if (anyNA(dt[, ..check_cols]))
  stop("Missing values found in required columns of the master database.")
if (nrow(dt) != 1295L || uniqueN(dt$iso3) != 185L)
  stop("Expected 1,295 rows and 185 countries.")
if (!identical(sort(unique(dt$year)), study_years))
  stop("Unexpected study years detected.")
if (nrow(dt[, .N, by = .(iso3, year)][N > 1L]))
  stop("Duplicate country-year keys detected.")

for (z in c(sdi_col, daly_col, prev_col, diet$score, unname(env_raw)))
  set(dt, j = z, value = as.numeric(dt[[z]]))
dt[, ln_daly := log(get(daly_col))]
dt[, ln_prevalence := log(get(prev_col))]

env_log <- character(length(env_raw))
names(env_log) <- names(env_raw)
for (i in seq_along(env_raw)) {
  z <- env_raw[i]
  if (startsWith(z, "ln_")) {
    env_log[i] <- z
  } else {
    new_z <- paste0("ln_", names(env_raw)[i])
    if (any(dt[[z]] <= 0)) stop("Non-positive values in ", z)
    dt[, (new_z) := log(get(z))]
    env_log[i] <- new_z
  }
}
for (i in seq_len(nrow(diet))) {
  new_z <- paste0("x_", tolower(diet$diet[i]), "_per_increment")
  dt[, (new_z) := get(diet$score[i]) / diet$increment[i]]
  diet[i, exposure := new_z]
}
if (!is.na(coverage_col)) dt[, coverage_per_10pct := get(coverage_col) / 0.10]
fwrite(dt, file.path(data_dir, "analysis_master_checked_v2.csv"), bom = TRUE)

fit_twfe <- function(data, outcome, exposure, covars = character(),
                     time_var = "year") {
  rhs <- paste(c(exposure, covars), collapse = " + ")
  fml <- as.formula(paste0(outcome, " ~ ", rhs, " | iso3 + ", time_var))
  feols(fml, data = data, vcov = ~iso3, ssc = ssc_setting, notes = FALSE)
}
extract_effect <- function(m, term, n_countries = NA_integer_) {
  ct <- coeftable(m)
  if (!term %in% rownames(ct)) stop("Term not estimated: ", term)
  ci <- confint(m, parm = term, level = .95)
  b <- as.numeric(ct[term, "Estimate"])
  lo <- as.numeric(ci[1L, 1L]); hi <- as.numeric(ci[1L, 2L])
  data.table(beta = b, standard_error = as.numeric(ct[term, "Std. Error"]),
             p_value = as.numeric(ct[term, "Pr(>|t|)"]),
             beta_lower_95 = lo, beta_upper_95 = hi,
             percent_difference = 100 * (exp(b) - 1),
             percent_difference_lower_95 = 100 * (exp(lo) - 1),
             percent_difference_upper_95 = 100 * (exp(hi) - 1),
             observations = nobs(m), countries = as.integer(n_countries))
}
safe_effect <- function(data, outcome, exposure, covars, time_var) {
  tryCatch(extract_effect(fit_twfe(data, outcome, exposure, covars, time_var),
                          exposure, uniqueN(data$iso3)),
           error = function(e)
             data.table(fit_status = "ERROR", fit_error = conditionMessage(e)))
}

make_health_lag <- function(lag_years, interpolate = FALSE) {
  outcomes <- dt[, .(iso3, outcome_year = year, ln_daly, ln_prevalence)]
  source <- dt[, .(iso3, exposure_year = year,
                   sdi_value = get(sdi_col),
                   x_mphdi_per_increment, x_med_per_increment,
                   x_dash_per_increment, x_ahei_per_increment)]
  if (!interpolate) {
    source[, outcome_year := exposure_year + lag_years]
    source[, exposure_year := NULL]
    return(merge(outcomes, source,
                 by = c("iso3", "outcome_year"), all = FALSE))
  }
  wanted <- outcomes[, .(iso3, outcome_year,
                         exposure_year = outcome_year - lag_years)]
  vars <- c("sdi_value", "x_mphdi_per_increment", "x_med_per_increment",
            "x_dash_per_increment", "x_ahei_per_increment")
  filled <- source[, {
    xy <- source[iso3 == .BY$iso3]
    q <- wanted[iso3 == .BY$iso3, .(outcome_year, exposure_year)]
    for (v in vars) {
      q[, (v) := approx(xy$exposure_year, xy[[v]],
                        xout = exposure_year, rule = 1, ties = mean)$y]
    }
    q
  }, by = iso3]
  merge(outcomes, filled, by = c("iso3", "outcome_year"), all = FALSE)
}

health10 <- make_health_lag(primary_health_lag, FALSE)
if (nrow(health10) != 740L)
  stop("Exact 10-year health data should contain 740 rows.")
health_covars <- "sdi_value"
env_covars <- if (is.na(coverage_col)) sdi_col else c(sdi_col, "coverage_per_10pct")

outcomes <- data.table(
  outcome_code = c("daly", "prevalence", "ghg", "land", "water", "eutrophication"),
  outcome_label = c("Dementia DALY rate", "Dementia prevalence rate",
                    "GHG emissions", "Land use",
                    "Water-scarcity", "Eutrophication"),
  family = c("Health", "Health", rep("Environment", 4L)),
  outcome_column = c("ln_daly", "ln_prevalence", unname(env_log)),
  health = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
)

# ---- Leave-one-country-out sensitivity analyses -----------------------------
specs <- rbindlist(lapply(seq_len(nrow(outcomes)), function(i)
  rbindlist(lapply(seq_len(nrow(diet)), function(j)
    data.table(outcome_code = outcomes$outcome_code[i],
               outcome_label = outcomes$outcome_label[i],
               family = outcomes$family[i],
               outcome_column = outcomes$outcome_column[i],
               diet = diet$diet[j], exposure = diet$exposure[j],
               increment_label = diet$increment_label[j])))))

loco <- list(); qc <- list(); nr <- 0L; nq <- 0L
for (i in seq_len(nrow(specs))) {
  sp <- specs[i]
  source <- if (sp$family == "Health") health10 else dt
  time_var <- if (sp$family == "Health") "outcome_year" else "year"
  covars <- if (sp$family == "Health") health_covars else env_covars
  for (omit in sort(unique(dt$iso3))) {
    ans <- safe_effect(source[iso3 != omit], sp$outcome_column,
                       sp$exposure, covars, time_var)
    nq <- nq + 1L
    qc[[nq]] <- data.table(outcome_code = sp$outcome_code,
                           outcome_label = sp$outcome_label,
                           family = sp$family, diet = sp$diet,
                           omitted_iso3 = omit,
                           fit_status = if ("fit_status" %in% names(ans))
                             ans$fit_status else "SUCCESS",
                           fit_error = if ("fit_error" %in% names(ans))
                             ans$fit_error else NA_character_)
    if (!"fit_status" %in% names(ans)) {
      nr <- nr + 1L
      loco[[nr]] <- cbind(data.table(
        outcome_code = sp$outcome_code, outcome_label = sp$outcome_label,
        family = sp$family, diet = sp$diet,
        increment_label = sp$increment_label, omitted_iso3 = omit), ans)
    }
  }
}
loco_results <- rbindlist(loco, fill = TRUE)
loco_results[, q_value_bh := p.adjust(p_value, "BH"),
            by = .(outcome_code, family, diet)]
fwrite(loco_results, file.path(table_dir,
       "TABLE_leave_one_country_out_all_estimates_v3.csv"), bom = TRUE)
fwrite(rbindlist(qc, fill = TRUE), file.path(table_dir,
       "QC_leave_one_country_out_fit_status_v3.csv"), bom = TRUE)
loco_summary <- loco_results[, .(
  omission_models = .N,
  estimate_median = median(percent_difference),
  estimate_min = min(percent_difference),
  estimate_max = max(percent_difference),
  estimate_q025 = quantile(percent_difference, .025),
  estimate_q975 = quantile(percent_difference, .975),
  p_min = min(p_value), p_max = max(p_value),
  positive_n = sum(percent_difference > 0),
  negative_n = sum(percent_difference < 0)
), by = .(outcome_code, outcome_label, family, diet, increment_label)]
## Old manuscript-style summary: full-sample estimate and 95% CI, range after
## country omission, and maximum absolute change in percentage points.
full_results <- rbindlist(lapply(seq_len(nrow(specs)), function(i) {
  sp <- specs[i]
  source <- if (sp$family == "Health") health10 else dt
  time_var <- if (sp$family == "Health") "outcome_year" else "year"
  covars <- if (sp$family == "Health") health_covars else env_covars
  ans <- safe_effect(source, sp$outcome_column, sp$exposure, covars, time_var)
  if ("fit_status" %in% names(ans)) return(NULL)
  cbind(sp[, .(outcome_code, outcome_label, family, diet, increment_label)], ans)
}), fill = TRUE)
full_key <- full_results[, .(outcome_code, diet,
  full_estimate = percent_difference,
  full_lower_95 = percent_difference_lower_95,
  full_upper_95 = percent_difference_upper_95)]
loco_summary <- loco_results[, .(
  range_after_country_omission_min = min(percent_difference),
  range_after_country_omission_max = max(percent_difference)
), by = .(outcome_code, outcome_label, family, diet, increment_label)]
loco_summary <- merge(loco_summary, full_key,
                      by = c("outcome_code", "diet"), all.x = TRUE)
loco_summary[, maximum_absolute_change_percentage_points := pmax(
  abs(range_after_country_omission_min - full_estimate),
  abs(range_after_country_omission_max - full_estimate))]
setcolorder(loco_summary, c("outcome_code", "outcome_label", "family", "diet",
  "increment_label", "full_estimate", "full_lower_95", "full_upper_95",
  "range_after_country_omission_min", "range_after_country_omission_max",
  "maximum_absolute_change_percentage_points"))
fwrite(loco_summary, file.path(table_dir,
       "TABLE_leave_one_country_out_summary_v3.csv"), bom = TRUE)

# Previous-style LOCO stability figure: 2 columns x 3 rows.
plot_dt <- merge(loco_results, full_key[, .(outcome_code, diet, full_estimate)],
                 by = c("outcome_code", "diet"), all.x = TRUE)
plot_dt[, outcome_label := factor(outcome_label, levels = outcomes$outcome_label)]
plot_dt[, diet := factor(diet, levels = c("mPHDI", "MED", "DASH", "AHEI"))]
diet_cols <- c(mPHDI = "#C05AB7", MED = "#8B83B9",
               DASH = "#73B7B3", AHEI = "#75ADD0")
stability_plot <- ggplot(plot_dt, aes(x = percent_difference,
                                      fill = diet, colour = diet)) +
  geom_histogram(bins = 14, alpha = 0.78, position = "identity", linewidth = 0.15) +
  geom_vline(aes(xintercept = full_estimate, colour = diet),
             linetype = "dashed", linewidth = 0.65, show.legend = FALSE) +
  facet_wrap(~ outcome_label, ncol = 2, scales = "free") +
  scale_fill_manual(values = diet_cols, drop = FALSE) +
  scale_colour_manual(values = diet_cols, drop = FALSE) +
  labs(x = "Percent difference by dietary pattern exposure unit",
       y = "Number of country-omission models") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())
ggsave(file.path(figure_dir, "Supplementary_Figure_S3_leave_one_country_out_v3.png"),
       stability_plot, width = 10, height = 8, dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "Supplementary_Figure_S3_leave_one_country_out_v3.pdf"),
       stability_plot, width = 10, height = 8, device = cairo_pdf)

# ---- Tests for nonlinear associations: quadratic exposure term ---------------
nonlinear <- list(); nr <- 0L
for (i in seq_len(nrow(outcomes))) {
  oo <- outcomes[i]
  source <- if (oo$health) health10 else dt
  time_var <- if (oo$health) "outcome_year" else "year"
  covars <- if (oo$health) health_covars else env_covars
  for (j in seq_len(nrow(diet))) {
    dd <- diet[j]; x2 <- paste0(dd$exposure, "_squared")
    source[, (x2) := get(dd$exposure)^2]
    fit <- tryCatch(fit_twfe(source, oo$outcome_column, dd$exposure,
                             c(covars, x2), time_var),
                    error = function(e) e)
    nr <- nr + 1L
    if (inherits(fit, "error")) {
      nonlinear[[nr]] <- data.table(
        outcome_code = oo$outcome_code, outcome_label = oo$outcome_label,
        family = oo$family, diet = dd$diet,
        increment_label = dd$increment_label,
        nonlinear_term = x2, fit_status = "ERROR",
        fit_error = conditionMessage(fit))
    } else {
      ct <- coeftable(fit)
      nonlinear[[nr]] <- data.table(
        outcome_code = oo$outcome_code, outcome_label = oo$outcome_label,
        family = oo$family, diet = dd$diet,
        increment_label = dd$increment_label,
        nonlinear_term = x2, quadratic_beta = as.numeric(ct[x2, "Estimate"]),
        quadratic_se = as.numeric(ct[x2, "Std. Error"]),
        nonlinear_p_value = as.numeric(ct[x2, "Pr(>|t|)"]),
        observations = nobs(fit), countries = uniqueN(source$iso3),
        fit_status = "SUCCESS")
    }
  }
}
nonlinear_results <- rbindlist(nonlinear, fill = TRUE)
nonlinear_results[, nonlinear_q_value_bh := p.adjust(nonlinear_p_value, "BH"),
                  by = .(outcome_code, family)]
fwrite(nonlinear_results, file.path(table_dir,
       "TABLE_nonlinear_quadratic_tests_v3.csv"), bom = TRUE)

# ---- Interpolation-based exploratory analysis --------------------------------
interpolation <- list(); nr <- 0L
for (lag_i in health_lags) {
  lag_data <- make_health_lag(lag_i, TRUE)
  for (i in 1:2) {
    oo <- outcomes[i]
    for (j in seq_len(nrow(diet))) {
      dd <- diet[j]
      fit <- tryCatch(fit_twfe(lag_data, oo$outcome_column, dd$exposure,
                               health_covars, "outcome_year"),
                      error = function(e) e)
      nr <- nr + 1L
      if (inherits(fit, "error")) {
        interpolation[[nr]] <- data.table(
          analysis = "Interpolation-based exploratory analysis",
          lag_years = lag_i, outcome_code = oo$outcome_code,
          outcome_label = oo$outcome_label, diet = dd$diet,
          increment_label = dd$increment_label, fit_status = "ERROR",
          fit_error = conditionMessage(fit))
      } else {
        interpolation[[nr]] <- cbind(data.table(
          analysis = "Interpolation-based exploratory analysis",
          lag_years = lag_i, outcome_code = oo$outcome_code,
          outcome_label = oo$outcome_label, diet = dd$diet,
          increment_label = dd$increment_label, fit_status = "SUCCESS"),
          extract_effect(fit, dd$exposure, uniqueN(lag_data$iso3)))
      }
    }
  }
}
interpolation_results <- rbindlist(interpolation, fill = TRUE)
interpolation_results[, q_value_bh := p.adjust(p_value, "BH"),
                      by = .(lag_years, outcome_code)]
fwrite(interpolation_results, file.path(table_dir,
       "TABLE_interpolation_based_exploratory_health_v3.csv"), bom = TRUE)

target_audit <- rbindlist(lapply(health_lags, function(z)
  data.table(lag_years = z, outcome_year = study_years,
             exposure_year = study_years - z)))
target_audit[, observed_exact := exposure_year %in% study_years]
target_audit[, interpolation_used := !observed_exact]
fwrite(target_audit, file.path(table_dir,
       "QC_interpolation_target_years_v3.csv"), bom = TRUE)

message("Completed: ", output_root)
message("LOCO result rows: ", nrow(loco_results))
message("Nonlinear result rows: ", nrow(nonlinear_results))
message("Interpolation result rows: ", nrow(interpolation_results))
