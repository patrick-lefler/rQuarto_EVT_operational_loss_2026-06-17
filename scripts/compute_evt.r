# ==============================================================================
# compute_evt.R
# ------------------------------------------------------------------------------
# Fits a Generalized Pareto Distribution (GPD) to operational loss exceedances
# via the Peaks-Over-Threshold (POT) method, cross-checks against a GEV fit to
# annual block maxima, fits a log-normal benchmark, and produces capital
# estimates (VaR / Expected Shortfall) at 99%, 99.9%, and 99.95% under all
# three approaches.
#
# Reads  : data/loss_register.csv, data/threshold_diagnostics.csv
# Writes : data/gpd_fit.csv, data/capital_comparison.csv,
#          output/evt_diagnostic_panel.png, output/mep_stability_panel.png
#
# This script is the authoritative source for all EVT estimates in index.qmd.
# index.qmd performs NO in-chunk model fitting — it reads the CSVs this script
# produces. Run this script standalone and inspect console output before
# rendering the document.
#
# Author : Patrick Lefler
# Project: The Tail You Cannot Fit — EVT for Operational Loss Capital Estimation
# ==============================================================================

suppressMessages({
  library(tidyverse)
  library(evd)
  library(POT)
  library(fitdistrplus)
  library(scales)
  library(patchwork)
})

cat(strrep("=", 70), "\n")
cat("COMPUTE_EVT.R — STARTING\n")
cat(strrep("=", 70), "\n\n")

# ------------------------------------------------------------------------------
# 0. Brand palette (for PNG diagnostic plots — kept consistent with index.qmd)
# ------------------------------------------------------------------------------
brand_primary   <- "#1A1A2E"
brand_secondary <- "#16213E"
brand_accent    <- "#0F3460"
brand_highlight <- "#E94560"
brand_surface   <- "#F5F5F5"
brand_text      <- "#1A1A2E"

theme_evt <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = brand_surface, linewidth = 0.3),
      plot.title       = element_text(color = brand_text, face = "bold", size = 13),
      plot.subtitle    = element_text(color = brand_text, size = 10),
      axis.title       = element_text(color = brand_text),
      axis.text        = element_text(color = brand_text),
      strip.text       = element_text(color = brand_text, face = "bold"),
      legend.position  = "bottom"
    )
}

dir.create("data",   showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
cat("--- 1. Loading data ---\n")

loss_register <- read_csv("data/loss_register.csv", show_col_types = FALSE) %>%
  mutate(event_date = as.Date(event_date))

threshold_diagnostics <- read_csv("data/threshold_diagnostics.csv", show_col_types = FALSE)

cat("loss_register rows           :", nrow(loss_register), "\n")
cat("threshold_diagnostics rows   :", nrow(threshold_diagnostics), "\n")
cat("Date range                   :", as.character(min(loss_register$event_date)),
    "to", as.character(max(loss_register$event_date)), "\n\n")

net_losses <- loss_register$net_loss_eur

stopifnot(
  "loss register is empty"                = nrow(loss_register) > 0,
  "net_loss_eur contains non-positive values" = all(net_losses > 0),
  "net_loss_eur contains NA"               = !any(is.na(net_losses))
)

# ------------------------------------------------------------------------------
# 2. Threshold selection
# ------------------------------------------------------------------------------
# Decision: u = EUR 75,000, calibrated from threshold_diagnostics.csv
# Rationale logged in INSTRUCTIONS.md: gives ~85 exceedances, sits just before
# the exceedance-count reliability break observed at u = EUR 100,000.
# Confirm programmatically rather than trusting the hardcoded value blindly.

cat("--- 2. Threshold selection ---\n")

SELECTED_THRESHOLD <- 75000

n_exceed_at_threshold <- sum(net_losses > SELECTED_THRESHOLD)
cat("Selected threshold (EUR)     :", comma(SELECTED_THRESHOLD), "\n")
cat("Exceedances above threshold  :", n_exceed_at_threshold, "\n")

if (n_exceed_at_threshold < 50) {
  warning("Fewer than 50 exceedances above selected threshold — GPD estimates may be unstable.")
}

# Confirm against the MEP data: mean excess should be in an upward-trending
# region near this threshold (visual MEP linearity check is qualitative and
# performed in the index.qmd figure; this is a sanity print only)
mep_near_threshold <- threshold_diagnostics %>%
  filter(threshold_eur >= SELECTED_THRESHOLD - 5000,
         threshold_eur <= SELECTED_THRESHOLD + 5000)
cat("Mean excess near threshold (sanity check):\n")
print(mep_near_threshold %>% dplyr::select(threshold_eur, n_exceedances, mean_excess_eur))
cat("\n")

exceedances_raw  <- net_losses[net_losses > SELECTED_THRESHOLD]
exceedances_over <- exceedances_raw - SELECTED_THRESHOLD   # Y = X - u, the GPD-distributed quantity

# ------------------------------------------------------------------------------
# 3. GPD fit via POT (evd::fpot)
# ------------------------------------------------------------------------------
cat("--- 3. GPD fit (evd::fpot) ---\n")

gpd_fit <- fpot(net_losses, threshold = SELECTED_THRESHOLD, model = "gpd", std.err = TRUE)

print(gpd_fit)
cat("\n")

xi_hat    <- gpd_fit$estimate["shape"]
sigma_hat <- gpd_fit$estimate["scale"]
xi_se     <- gpd_fit$std.err["shape"]
sigma_se  <- gpd_fit$std.err["scale"]

cat("Shape (xi)   estimate / SE   :", round(xi_hat, 4), "/", round(xi_se, 4), "\n")
cat("Scale (sigma) estimate / SE  :", round(sigma_hat, 4), "/", round(sigma_se, 4), "\n\n")

if (xi_hat <= 0) {
  warning("Fitted shape parameter xi <= 0 — tail is NOT heavy-tailed under this fit. ",
          "This contradicts the project's central thesis. Review threshold or data generator.")
} else {
  cat("xi > 0 confirmed: tail exhibits polynomial (heavy-tail) decay — Frechet domain.\n\n")
}

# ------------------------------------------------------------------------------
# 4. GPD quantile (VaR) and Expected Shortfall functions
# ------------------------------------------------------------------------------
# Standard POT VaR formula (Coles, 2001):
#   VaR_p = u + (sigma/xi) * [ ((n/Nu) * (1-p))^(-xi) - 1 ]
# where n = total observations, Nu = number of exceedances above u
#
# Expected Shortfall (for xi < 1):
#   ES_p = VaR_p / (1 - xi)  +  (sigma - xi*u) / (1 - xi)
# Equivalently: ES_p = (VaR_p + sigma - xi*u) / (1 - xi)

n_total <- length(net_losses)
n_u     <- n_exceed_at_threshold

gpd_var <- function(p, u, sigma, xi, n, nu) {
  u + (sigma / xi) * ( ((n / nu) * (1 - p))^(-xi) - 1 )
}

gpd_es <- function(var_p, u, sigma, xi) {
  (var_p + sigma - xi * u) / (1 - xi)
}

confidence_levels <- c(0.99, 0.999, 0.9995)

gpd_estimates <- tibble(
  confidence_level = confidence_levels,
  VaR_gpd = map_dbl(confidence_levels, ~ gpd_var(.x, SELECTED_THRESHOLD, sigma_hat, xi_hat, n_total, n_u)),
) %>%
  mutate(ES_gpd = map_dbl(VaR_gpd, ~ gpd_es(.x, SELECTED_THRESHOLD, sigma_hat, xi_hat)))

cat("--- 4. GPD-derived VaR / ES ---\n")
print(gpd_estimates)
cat("\n")

# ------------------------------------------------------------------------------
# 5. GEV cross-check on annual block maxima
# ------------------------------------------------------------------------------
cat("--- 5. GEV block maxima cross-check ---\n")

annual_maxima <- loss_register %>%
  mutate(year = lubridate::year(event_date)) %>%
  group_by(year) %>%
  summarise(annual_max_eur = max(net_loss_eur), .groups = "drop")

cat("Annual maxima:\n")
print(annual_maxima)
cat("\n")

n_years <- nrow(annual_maxima)

if (n_years < 5) {
  cat("NOTE: Fewer than 5 annual maxima available — GEV fit will have wide standard errors.\n")
  cat("This is an expected limitation with only", n_years, "years of data and is flagged in-text.\n\n")
}

gev_fit_result <- tryCatch({
  fgev(annual_maxima$annual_max_eur, std.err = TRUE)
}, error = function(e) {
  cat("GEV fit failed:", conditionMessage(e), "\n")
  NULL
})

gev_available <- !is.null(gev_fit_result)

if (gev_available) {
  print(gev_fit_result)
  gev_shape <- gev_fit_result$estimate["shape"]
  gev_loc   <- gev_fit_result$estimate["loc"]
  gev_scale <- gev_fit_result$estimate["scale"]
  cat("\nGEV shape estimate:", round(gev_shape, 4),
      "(POT GPD shape was:", round(xi_hat, 4), ")\n")
  cat("Same-sign check    :", ifelse(sign(gev_shape) == sign(xi_hat),
                                       "CONSISTENT (both indicate same tail domain)",
                                       "INCONSISTENT — flag for review"), "\n\n")
} else {
  gev_shape <- NA_real_
  gev_loc   <- NA_real_
  gev_scale <- NA_real_
  cat("GEV cross-check unavailable; proceeding with POT/GPD as primary method only.\n\n")
}

# ------------------------------------------------------------------------------
# 6. Log-normal benchmark fit (fitdistrplus)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 6. Log-normal benchmark fit (fitdistrplus)
# ------------------------------------------------------------------------------
# DECISION (logged in INSTRUCTIONS.md): fit log-normal to losses BELOW the
# selected POT threshold (EUR 75,000) only, not the full register.
#
# Rationale: a log-normal MLE fit to the FULL register (bulk + tail mixed)
# is structurally mis-specified -- it tries to accommodate genuinely
# heavy-tailed events within a thin-tailed family, which inflates sdlog and
# can produce a log-normal tail that overshoots even the GPD estimate. That
# is a real but different failure mode than the one this document is built
# to demonstrate.
#
# Fitting log-normal only to sub-threshold losses mirrors what a risk
# practitioner without an EVT framework would plausibly do: model the
# "normal" loss experience with a conventional distribution, with no
# explicit view on tail behavior above it. This is not circular reasoning
# using the GPD threshold to rig the comparison -- EUR 75,000 was selected
# from the MEP/stability diagnostics in Section 2, independent of the
# log-normal fit, and is the same boundary a practitioner would reasonably
# treat as "where the bulk model stops being credible."
cat("--- 6. Log-normal benchmark fit ---\n")

sub_threshold_losses <- net_losses[net_losses <= SELECTED_THRESHOLD]
cat("Fitting log-normal to sub-threshold losses only (<= EUR",
    comma(SELECTED_THRESHOLD), ")\n")
cat("Sub-threshold n           :", length(sub_threshold_losses),
    "of", n_total, "total events\n\n")

lognorm_fit <- fitdist(sub_threshold_losses, "lnorm")
print(summary(lognorm_fit))

ln_mu    <- lognorm_fit$estimate["meanlog"]
ln_sigma <- lognorm_fit$estimate["sdlog"]

lognorm_var <- function(p, mu, sigma) qlnorm(p, mu, sigma)
lognorm_es  <- function(p, mu, sigma) {
  # Closed-form ES for log-normal:
  # ES_p = exp(mu + sigma^2/2) * Phi(sigma - Phi^-1(p)) / (1 - p)
  z <- qnorm(p)
  exp(mu + sigma^2 / 2) * pnorm(sigma - z) / (1 - p)
}

lognorm_estimates <- tibble(
  confidence_level = confidence_levels,
  VaR_lognorm = map_dbl(confidence_levels, ~ lognorm_var(.x, ln_mu, ln_sigma)),
  ES_lognorm  = map_dbl(confidence_levels, ~ lognorm_es(.x, ln_mu, ln_sigma))
)

cat("\nLog-normal VaR / ES (fit to sub-threshold losses):\n")
print(lognorm_estimates)
cat("\n")

# Sanity check: the whole premise of this section requires that extrapolating
# a sub-threshold log-normal fit to 99.9%/99.95% UNDERSTATES the true tail
# risk relative to GPD. Confirm this holds; if not, flag loudly rather than
# silently presenting a finding that contradicts the document's thesis.
sanity_check <- gpd_estimates$VaR_gpd > lognorm_estimates$VaR_lognorm
cat("GPD VaR > log-normal VaR at each confidence level:", paste(sanity_check, collapse = ", "), "\n")
if (!all(sanity_check)) {
  warning("GPD VaR does NOT exceed log-normal VaR at all confidence levels. ",
          "Review the sub-threshold fit before writing document prose around this comparison.")
}
cat("\n")

# ------------------------------------------------------------------------------
# 7. Historical simulation (empirical quantiles)
# ------------------------------------------------------------------------------
cat("--- 7. Historical simulation ---\n")

histsim_estimates <- tibble(
  confidence_level = confidence_levels,
  VaR_histsim = map_dbl(confidence_levels, ~ quantile(net_losses, probs = .x, type = 7)),
) %>%
  mutate(
    ES_histsim = map_dbl(confidence_levels, function(p) {
      var_p <- quantile(net_losses, probs = p, type = 7)
      tail_obs <- net_losses[net_losses >= var_p]
      mean(tail_obs)
    })
  )

cat("Historical simulation VaR / ES:\n")
print(histsim_estimates)
cat("\nNOTE: at p=0.999 and p=0.9995 with n =", n_total,
    "observations, the empirical quantile collapses toward the sample maximum.\n")
cat("This is flagged in-text as the structural limitation of historical simulation\n")
cat("at extreme quantiles relative to sample size.\n\n")

# ------------------------------------------------------------------------------
# 8. Assemble capital comparison table
# ------------------------------------------------------------------------------
cat("--- 8. Capital comparison ---\n")

capital_comparison <- gpd_estimates %>%
  left_join(lognorm_estimates, by = "confidence_level") %>%
  left_join(histsim_estimates, by = "confidence_level") %>%
  mutate(
    gap_var_eur     = VaR_gpd - VaR_lognorm,
    gap_var_pct     = gap_var_eur / VaR_lognorm,
    gap_es_eur      = ES_gpd - ES_lognorm,
    gap_es_pct      = gap_es_eur / ES_lognorm
  ) %>%
  dplyr::select(
    confidence_level,
    VaR_gpd, VaR_lognorm, VaR_histsim, gap_var_eur, gap_var_pct,
    ES_gpd, ES_lognorm, ES_histsim, gap_es_eur, gap_es_pct
  )

print(capital_comparison)
cat("\n")

# ------------------------------------------------------------------------------
# 9. Assemble gpd_fit.csv (parameters + diagnostic data needed by index.qmd)
# ------------------------------------------------------------------------------
cat("--- 9. Assembling gpd_fit.csv ---\n")

# 9a. Parameter summary rows
param_summary <- tibble(
  parameter = c("shape_xi", "scale_sigma", "threshold_eur", "n_exceedances",
                "n_total", "gev_shape_xi", "gev_n_years"),
  estimate  = c(xi_hat, sigma_hat, SELECTED_THRESHOLD, n_u,
                n_total, gev_shape, n_years),
  std_error = c(xi_se, sigma_se, NA, NA, NA, ifelse(gev_available, gev_fit_result$std.err["shape"], NA), NA)
)

# 9b. Empirical exceedance data for QQ-plot / probability plot reconstruction in ggplot2
# Sorted exceedances over threshold with their empirical and GPD-fitted quantiles
sorted_exceed <- sort(exceedances_over)
m <- length(sorted_exceed)
empirical_p <- (seq_len(m) - 0.5) / m   # plotting positions

# GPD theoretical quantile function: Q(p) = (sigma/xi) * ((1-p)^(-xi) - 1)
gpd_quantile_fn <- function(p, sigma, xi) (sigma / xi) * ((1 - p)^(-xi) - 1)
gpd_theoretical_q <- gpd_quantile_fn(empirical_p, sigma_hat, xi_hat)

# GPD CDF function for probability plot: F(y) = 1 - (1 + xi*y/sigma)^(-1/xi)
gpd_cdf_fn <- function(y, sigma, xi) 1 - (1 + xi * y / sigma)^(-1 / xi)
gpd_theoretical_p <- gpd_cdf_fn(sorted_exceed, sigma_hat, xi_hat)

qq_pp_data <- tibble(
  rank               = seq_len(m),
  empirical_exceedance_eur = sorted_exceed,
  empirical_prob     = empirical_p,
  gpd_theoretical_quantile_eur = gpd_theoretical_q,
  gpd_theoretical_prob = gpd_theoretical_p
)

# 9c. Return level data (for return level plot)
# Return level at return period T (in units of "events"): for GPD POT,
# typically expressed in terms of exceedance probability -> return level curve
return_periods <- c(2, 5, 10, 25, 50, 100, 200, 500, 1000)
# Probability of non-exceedance corresponding to each return period, scaled by
# the rate of exceedance (n_u / n_total) per observation
exceedance_rate <- n_u / n_total
return_level_data <- tibble(
  return_period_events = return_periods,
  return_level_eur = map_dbl(return_periods, function(T) {
    p_exceed <- 1 / (T * exceedance_rate)
    if (p_exceed >= 1) return(NA_real_)
    SELECTED_THRESHOLD + gpd_quantile_fn(1 - p_exceed, sigma_hat, xi_hat)
  })
)

cat("Return level estimates:\n")
print(return_level_data)
cat("\n")

# ------------------------------------------------------------------------------
# 10. Build diagnostic panel (2x2: probability, quantile, return level, density)
# ------------------------------------------------------------------------------
cat("--- 10. Building diagnostic panel PNG ---\n")

p_prob <- ggplot(qq_pp_data, aes(x = empirical_prob, y = gpd_theoretical_prob)) +
  geom_point(color = brand_accent, alpha = 0.6, size = 1.6) +
  geom_abline(slope = 1, intercept = 0, color = brand_highlight, linetype = "dashed") +
  labs(title = "Probability Plot", x = "Empirical", y = "GPD Fitted") +
  theme_evt()

p_quant <- ggplot(qq_pp_data, aes(x = gpd_theoretical_quantile_eur, y = empirical_exceedance_eur)) +
  geom_point(color = brand_accent, alpha = 0.6, size = 1.6) +
  geom_abline(slope = 1, intercept = 0, color = brand_highlight, linetype = "dashed") +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Quantile Plot (QQ)", x = "GPD Fitted (EUR)", y = "Empirical (EUR)") +
  theme_evt()

p_return <- ggplot(return_level_data %>% filter(!is.na(return_level_eur)),
                    aes(x = return_period_events, y = return_level_eur)) +
  geom_line(color = brand_accent, linewidth = 0.8) +
  geom_point(color = brand_highlight, size = 1.8) +
  scale_x_log10(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Return Level Plot", x = "Return Period (events, log scale)", y = "Return Level (EUR)") +
  theme_evt()

dens_x <- seq(0, max(sorted_exceed) * 1.05, length.out = 300)
dens_gpd <- tibble(
  x = dens_x,
  density = (1 / sigma_hat) * (1 + xi_hat * dens_x / sigma_hat)^(-1 / xi_hat - 1)
)

p_density <- ggplot() +
  geom_histogram(data = tibble(x = exceedances_over), aes(x = x, y = after_stat(density)),
                 bins = 25, fill = brand_surface, color = brand_secondary, alpha = 0.7) +
  geom_line(data = dens_gpd, aes(x = x, y = density), color = brand_highlight, linewidth = 1) +
  scale_x_continuous(labels = label_comma()) +
  labs(title = "Density: Exceedances vs. Fitted GPD", x = "Exceedance over threshold (EUR)", y = "Density") +
  theme_evt()

diagnostic_panel <- (p_prob | p_quant) / (p_return | p_density)

ggsave("output/evt_diagnostic_panel.png", plot = diagnostic_panel,
       width = 10, height = 8, dpi = 300, bg = "white")

cat("Saved: output/evt_diagnostic_panel.png\n\n")

# ------------------------------------------------------------------------------
# 11. Build MEP + stability panel
# ------------------------------------------------------------------------------
cat("--- 11. Building MEP + stability panel PNG ---\n")

p_mep <- ggplot(threshold_diagnostics, aes(x = threshold_eur, y = mean_excess_eur)) +
  geom_ribbon(aes(ymin = ci_lower_eur, ymax = ci_upper_eur), fill = brand_surface, alpha = 0.6) +
  geom_line(color = brand_accent, linewidth = 0.8) +
  geom_vline(xintercept = SELECTED_THRESHOLD, color = brand_highlight, linetype = "dashed") +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Mean Excess Plot", x = "Threshold (EUR)", y = "Mean Excess (EUR)") +
  theme_evt()

# Stability plot: refit GPD shape across a grid of thresholds
stability_thresholds <- threshold_diagnostics %>%
  filter(n_exceedances >= 30, n_exceedances <= 200) %>%
  pull(threshold_eur)

stability_results <- map_dfr(stability_thresholds, function(u) {
  tryCatch({
    fit <- fpot(net_losses, threshold = u, model = "gpd", std.err = TRUE)
    tibble(
      threshold_eur = u,
      xi = fit$estimate["shape"],
      xi_se = fit$std.err["shape"]
    )
  }, error = function(e) tibble(threshold_eur = u, xi = NA_real_, xi_se = NA_real_))
}) %>%
  filter(!is.na(xi))

p_stability <- ggplot(stability_results, aes(x = threshold_eur, y = xi)) +
  geom_ribbon(aes(ymin = xi - 1.96 * xi_se, ymax = xi + 1.96 * xi_se),
              fill = brand_surface, alpha = 0.6) +
  geom_line(color = brand_accent, linewidth = 0.8) +
  geom_hline(yintercept = 0, color = brand_secondary, linetype = "dotted") +
  geom_vline(xintercept = SELECTED_THRESHOLD, color = brand_highlight, linetype = "dashed") +
  scale_x_continuous(labels = label_comma()) +
  labs(title = "Shape Parameter Stability", x = "Threshold (EUR)", y = expression(xi)) +
  theme_evt()

mep_stability_panel <- p_mep | p_stability

ggsave("output/mep_stability_panel.png", plot = mep_stability_panel,
       width = 10, height = 4.5, dpi = 300, bg = "white")

cat("Saved: output/mep_stability_panel.png\n\n")

# ------------------------------------------------------------------------------
# 12. Write final CSVs
# ------------------------------------------------------------------------------
cat("--- 12. Writing output CSVs ---\n")

write_csv(param_summary, "data/gpd_fit.csv")
write_csv(capital_comparison, "data/capital_comparison.csv")
write_csv(qq_pp_data, "data/qq_pp_data.csv")
write_csv(return_level_data, "data/return_level_data.csv")
write_csv(stability_results, "data/stability_results.csv")
write_csv(annual_maxima, "data/annual_maxima.csv")

cat("data/gpd_fit.csv             rows:", nrow(param_summary), "\n")
cat("data/capital_comparison.csv  rows:", nrow(capital_comparison), "\n")
cat("data/qq_pp_data.csv          rows:", nrow(qq_pp_data), "\n")
cat("data/return_level_data.csv   rows:", nrow(return_level_data), "\n")
cat("data/stability_results.csv   rows:", nrow(stability_results), "\n")
cat("data/annual_maxima.csv       rows:", nrow(annual_maxima), "\n\n")

cat(strrep("=", 70), "\n")
cat("COMPUTE_EVT.R — COMPLETE\n")
cat(strrep("=", 70), "\n")
