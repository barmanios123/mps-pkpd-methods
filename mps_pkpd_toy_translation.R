############################################################
# mps_pkpd_toy_translation.R
#
# Purpose:
# Exploratory translational PK/PD prototype linking:
#   1) MPS APOC3 mRNA knockdown summary data
#   2) A toy clinical-like PK model
#   3) A simple Emax clinical PD model
#
# Notes:
# - This script is for methods illustration and hypothesis generation.
# - The MPS inputs below are manually encoded summary values.
# - Clinical PK and PD parameters are illustrative, not validated.
# - This is not a production model and is not integrated into the
#   main Python scoring pipeline.
#
# Suggested usage:
#   source("scripts/mps_pkpd_toy_translation.R")
#
# Outputs:
#   outputs/pkpd_demo/chip_pd_points.csv
#   outputs/pkpd_demo/chip_emax_params.csv
#   outputs/pkpd_demo/chip_pred_curve.csv
#   outputs/pkpd_demo/clinical_pkpd_timeseries.csv
#   outputs/pkpd_demo/chip_ec50_reference_points.csv
#   outputs/pkpd_demo/time_above_thresholds.csv
#   outputs/pkpd_demo/*.png
############################################################

suppressPackageStartupMessages({
  library(PKPDsim)
  library(ggplot2)
  library(dplyr)
  library(minpack.lm)
})

out_dir <- "outputs/pkpd_demo"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

theme_set(theme_minimal(base_size = 13))

############################################################
# 1. MPS summary data
############################################################
# Approximate APOC3 mRNA fractions from summarized Day 14 PHH+KC MPS data.
# Entered manually from summary visualization / extracted summary values.
# 0 nM      -> 1.00 mRNA fraction
# 100 nM    -> 0.16 mRNA fraction
# 10000 nM  -> 0.07 mRNA fraction

chip_pd <- data.frame(
  conc_nM       = c(0, 100, 10000),
  mRNA_fraction = c(1.00, 0.16, 0.07)
) %>%
  mutate(
    knockdown_pct = (1 - mRNA_fraction) * 100
  )

write.csv(chip_pd, file.path(out_dir, "chip_pd_points.csv"), row.names = FALSE)

############################################################
# 2. Fit chip Emax model
############################################################
# E(C) = Emax * C / (EC50 + C)

emax_fun <- function(C, Emax, EC50) {
  Emax * C / (EC50 + C)
}

fit_chip <- nlsLM(
  knockdown_pct ~ emax_fun(conc_nM, Emax, EC50),
  data = chip_pd,
  start = list(Emax = 95, EC50 = 10),
  control = nls.lm.control(maxiter = 200)
)

coef_chip <- coef(fit_chip)
Emax_chip <- unname(coef_chip["Emax"])
EC50_chip <- unname(coef_chip["EC50"])

chip_emax_params <- data.frame(
  parameter = c("Emax_chip", "EC50_chip_nM"),
  value = c(Emax_chip, EC50_chip)
)

write.csv(chip_emax_params, file.path(out_dir, "chip_emax_params.csv"), row.names = FALSE)

new_C <- 10^seq(-1, 5, length.out = 200)

pred_chip <- data.frame(
  conc_nM = new_C,
  knockdown_fit = emax_fun(new_C, Emax_chip, EC50_chip)
)

write.csv(pred_chip, file.path(out_dir, "chip_pred_curve.csv"), row.names = FALSE)

p_chip <- ggplot() +
  geom_point(
    data = chip_pd,
    aes(x = conc_nM, y = knockdown_pct),
    size = 3,
    color = "#b22222"
  ) +
  geom_line(
    data = pred_chip,
    aes(x = conc_nM, y = knockdown_fit),
    linewidth = 1.1,
    color = "black"
  ) +
  scale_x_log10() +
  labs(
    title = "MPS APOC3 knockdown vs concentration",
    subtitle = "Three-point summary-data Emax fit; exploratory only",
    x = "Concentration in MPS (nM, log scale)",
    y = "APOC3 knockdown (%)"
  )

ggsave(
  filename = file.path(out_dir, "mps_chip_emax_fit.png"),
  plot = p_chip,
  width = 8,
  height = 5,
  dpi = 300
)

############################################################
# 3. Toy clinical PK model
############################################################
# Simple 1-compartment oral/subQ-like model used as an exposure scaffold.
# Parameters below are illustrative and not fitted to trial data.

mod_pk <- new_ode_model("pk_1cmt_oral")

par_pk <- list(
  CL = 5,     # L/h
  V  = 50,    # L
  KA = 0.3    # 1/h
)

doses_mg <- c(10, 75, 225)

reg_list <- lapply(doses_mg, function(d) {
  new_regimen(
    amt = d,
    times = 0,
    type = "bolus",
    cmt = 1
  )
})
names(reg_list) <- paste0("Dose_", doses_mg, "mg")

t_obs <- seq(0, 2000, by = 12)

sim_pk_list <- lapply(seq_along(doses_mg), function(i) {
  out <- sim(
    ode = mod_pk,
    parameters = par_pk,
    regimen = reg_list[[i]],
    t_obs = t_obs
  )
  out$dose_mg <- doses_mg[i]
  out
})

sim_pk_all <- bind_rows(sim_pk_list)

############################################################
# 4. Toy clinical PD link
############################################################
# E(t) = Emax_eff * C(t) / (EC50_eff + C(t))
# These are illustrative effect parameters.

Emax_eff <- 95
EC50_eff <- 0.5

sim_pk_pd <- sim_pk_all %>%
  mutate(
    conc = y,
    effect_lpa = Emax_eff * conc / (EC50_eff + conc),
    lpa_pct_baseline = 100 - effect_lpa,
    multiple_EC50_eff = conc / EC50_eff
  )

write.csv(
  sim_pk_pd,
  file.path(out_dir, "clinical_pkpd_timeseries.csv"),
  row.names = FALSE
)

p_pk <- ggplot(sim_pk_pd, aes(x = t, y = conc, color = factor(dose_mg))) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Toy PK concentration-time profiles",
    subtitle = "Illustrative single-dose simulation",
    x = "Time (h)",
    y = "Concentration (arbitrary units)",
    color = "Dose (mg)"
  )

ggsave(
  filename = file.path(out_dir, "toy_pk_profiles.png"),
  plot = p_pk,
  width = 8,
  height = 5,
  dpi = 300
)

p_pd <- ggplot(sim_pk_pd, aes(x = t, y = lpa_pct_baseline, color = factor(dose_mg))) +
  geom_line(linewidth = 1.1) +
  ylim(0, 110) +
  labs(
    title = "Toy PD: Lp(a) percent of baseline vs time",
    subtitle = "Clinical PD effect derived from illustrative PK/PD parameters",
    x = "Time (h)",
    y = "Lp(a) (% baseline)",
    color = "Dose (mg)"
  )

ggsave(
  filename = file.path(out_dir, "toy_pd_profiles.png"),
  plot = p_pd,
  width = 8,
  height = 5,
  dpi = 300
)

############################################################
# 5. Compare EC50 multiples
############################################################

chip_points <- data.frame(
  conc_nM = c(100, 10000)
) %>%
  mutate(
    multiple_EC50_chip = conc_nM / EC50_chip,
    knockdown_pred = emax_fun(conc_nM, Emax_chip, EC50_chip)
  )

write.csv(
  chip_points,
  file.path(out_dir, "chip_ec50_reference_points.csv"),
  row.names = FALSE
)

p_norm <- ggplot(sim_pk_pd, aes(x = t, y = multiple_EC50_eff, color = factor(dose_mg))) +
  geom_line(linewidth = 1.1) +
  geom_hline(
    yintercept = chip_points$multiple_EC50_chip[1],
    linetype = "dashed",
    color = "darkgreen"
  ) +
  geom_hline(
    yintercept = chip_points$multiple_EC50_chip[2],
    linetype = "dotted",
    color = "darkred"
  ) +
  labs(
    title = "Clinical exposure as multiples of EC50_eff",
    subtitle = "Dashed = chip 100 nM reference; dotted = chip 10000 nM reference",
    x = "Time (h)",
    y = "C / EC50_eff",
    color = "Dose (mg)"
  )

ggsave(
  filename = file.path(out_dir, "ec50_multiple_comparison.png"),
  plot = p_norm,
  width = 8,
  height = 5,
  dpi = 300
)

############################################################
# 6. Time above EC50-based thresholds
############################################################

thresholds <- c(1, 3, 10)
dt <- diff(t_obs)[1]

time_above <- sim_pk_pd %>%
  group_by(dose_mg) %>%
  summarize(
    time_above_1x_h = sum(multiple_EC50_eff >= thresholds[1]) * dt,
    time_above_3x_h = sum(multiple_EC50_eff >= thresholds[2]) * dt,
    time_above_10x_h = sum(multiple_EC50_eff >= thresholds[3]) * dt,
    .groups = "drop"
  )

write.csv(
  time_above,
  file.path(out_dir, "time_above_thresholds.csv"),
  row.names = FALSE
)

message("Done. Outputs written to: ", out_dir)
message("Estimated chip Emax: ", round(Emax_chip, 1), "%")
message("Estimated chip EC50: ", signif(EC50_chip, 3), " nM")
