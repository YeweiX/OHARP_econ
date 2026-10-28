# Run from the 1028_code folder in RStudio (Session > Set Working Directory > To Source File Location).
# Fit7 posterior, 1,000 paired PSA simulations, seed 123.
library(rstan)
load("../updated_fit7_20261007/RData/fit7.RData")
posterior_samples <- as.data.frame(fit7, pars = c(
  "beta_hh", "beta_eh", "beta_hc", "beta_sc_ac",
  "beta_sc_ec", "beta_sh_sih", "delta_rcc_ric"
))
N_PSA <- 1000L
SEED <- 123L
Sys.setenv(N_PSA_SAMPLES = N_PSA, PSA_BASE_SEED = SEED, TOTAL_COST_N_CORES = 4L)
out <- "outputs/total_cost_fit7"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

summary_row <- function(x, outcome, unit = "million SGD", scale = 1e6) {
  stopifnot(length(x) > 0L, all(is.finite(x)))
  data.frame(outcome = outcome, unit = unit, n_draws = length(x),
    mean = mean(x) / scale,
    lower_95_UI = unname(quantile(x, 0.025, type = 7)) / scale,
    upper_95_UI = unname(quantile(x, 0.975, type = 7)) / scale)
}

write_result <- function(x, name, out) {
  write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
}

# Extracted from total_cost_0623.R, original lines 71-457.
library(deSolve)
library(dplyr)
library(foreach)
library(doSNOW)
library(purrr)

## ============================================================
## PART 1: GLOBAL SETTINGS AND PARAMETERS
## ============================================================

full_sim_years <- 20
economic_horizon_years <- 5
discount_rate <- 0.03
wtp_threshold <- 45000
life_expectancy_val <- 83.5
n_psa_samples <- as.integer(Sys.getenv("N_PSA_SAMPLES", unset = "1000"))
psa_base_seed <- as.integer(Sys.getenv("PSA_BASE_SEED", unset = "123"))
if (is.na(n_psa_samples) || n_psa_samples <= 0) {
  stop("N_PSA_SAMPLES must be a positive integer.")
}
if (is.na(psa_base_seed)) {
  stop("PSA_BASE_SEED must be an integer.")
}
set.seed(psa_base_seed)

fixed_params <- list(
  beta_hh = 0.00789,
  beta_sh_sih = 0.23953,
  delta_rch_rih = 0.03395,
  delta_rch_sih = 0.00085,
  beta_sc_sic = 7.92e-6,
  delta_rcc_ric = 8.56e-5,
  delta_rcc_sic = 0.00085,
  beta_eh = 0.00626,
  beta_hc = 0.00366,
  beta_sc_ac = 0.00068,
  beta_sc_ec = 0.00026,
  mu_sih = 0.00312,
  mu_rih = 0.00405,
  gamma_sih_sh = 0.14286,
  gamma_rch_sh = 0.00612,
  gamma_rih_rch = 0.11111,
  gamma_sic_sc = 0.25,
  gamma_rcc_sc = 0.009,
  gamma_ric_rcc = 0.17,
  mu_sic = 0.00323,
  mu_ric = 0.00363,
  alpha_adm = 0.00032,
  alpha_dis = 0.235,
  mu_b = 2.52e-5,
  mu_m = 2e-6,
  mu_d = 1.35e-5
)

mcmc_pars <- c(
  "beta_hh", "beta_eh", "beta_hc", "beta_sc_ac",
  "beta_sc_ec", "beta_sh_sih", "delta_rcc_ric"
)

if (!exists("posterior_samples")) {
  stop("posterior_samples is not found. Please load posterior_samples before running this script.")
}
missing_mcmc_cols <- setdiff(mcmc_pars, colnames(posterior_samples))
if (length(missing_mcmc_cols) > 0) {
  stop("posterior_samples is missing columns: ", paste(missing_mcmc_cols, collapse = ", "))
}
cat("posterior_samples rows:", nrow(posterior_samples), "\n")
stopifnot(!any(is.na(posterior_samples[, mcmc_pars])))

los_epi_ranges <- list(
  gamma_sih_sh  = list(min = 0.09090909, max = 0.20000000),
  gamma_rih_rch = list(min = 0.06666667, max = 0.16666667),
  mu_sih        = list(min = 0.00233333, max = 0.00390000),
  mu_rih        = list(min = 0.00254333, max = 0.00554667)
)

gamma_mom_from_mean_ci <- function(mean, lower, upper) {
  if (mean <= 0 || lower <= 0 || upper <= lower) {
    stop("Invalid mean/lower/upper in gamma_mom_from_mean_ci().")
  }
  sd <- (upper - lower) / (2 * 1.96)
  shape <- (mean / sd)^2
  scale <- sd^2 / mean
  list(alpha = shape, beta = scale)
}

gamma_mom_from_mean_sd <- function(mean, sd) {
  if (mean <= 0 || sd <= 0) {
    stop("Invalid mean/sd in gamma_mom_from_mean_sd().")
  }
  shape <- (mean / sd)^2
  scale <- sd^2 / mean
  list(alpha = shape, beta = scale)
}

cost_params_dist <- list(
  # Bed-day costs, 2024 SGD, fitted to the reported mean and SD.
  icu_bed_daily = gamma_mom_from_mean_sd(mean = 996.86, sd = 311.52),
  a_class_bed_daily = gamma_mom_from_mean_sd(mean = 454.06, sd = 133.15),
  b1_class_bed_daily = gamma_mom_from_mean_sd(mean = 297.02, sd = 92.72),
  b2_class_bed_daily = gamma_mom_from_mean_sd(mean = 219.56, sd = 90.19),
  c_class_bed_daily = gamma_mom_from_mean_sd(mean = 217.57, sd = 68.31),
  ward_bed_proportions = c(A = 0.068, B1 = 0.097, B2 = 0.373, C = 0.462),

  # Community antibiotic course costs, 2025 SGD, SGH-based.
  # ESBL: moxifloxacin 400 mg, $1.30/tab, 1 tab/day × 14 days = SGD 18.20.
  ric_antibiotic_course = gamma_mom_from_mean_ci(
    mean = 18.20, lower = 18.20 * 0.70, upper = 18.20 * 1.30
  ),

  # Non-ESBL: co-amoxiclav 1 g, $0.43/tab, 2 tabs/day × 14 days = SGD 12.04.
  sic_antibiotic_course = gamma_mom_from_mean_ci(
    mean = 12.04, lower = 12.04 * 0.70, upper = 12.04 * 1.30
  ),

  outpatient_visit = list(min = 17.5, max = 100)
)

required_cost_keys <- c(
  "icu_bed_daily", "a_class_bed_daily", "b1_class_bed_daily",
  "b2_class_bed_daily", "c_class_bed_daily", "ward_bed_proportions",
  "ric_antibiotic_course", "sic_antibiotic_course", "outpatient_visit"
)
missing_keys <- setdiff(required_cost_keys, names(cost_params_dist))
if (length(missing_keys) > 0) {
  stop("Missing cost_params_dist entries: ", paste(missing_keys, collapse = ", "))
}

cat(sprintf(
  "Community ESBL antibiotic course Gamma(shape = %.2f, scale = %.4f), mean = %.2f\n",
  cost_params_dist$ric_antibiotic_course$alpha,
  cost_params_dist$ric_antibiotic_course$beta,
  cost_params_dist$ric_antibiotic_course$alpha * cost_params_dist$ric_antibiotic_course$beta
))
cat(sprintf(
  "Community non-ESBL antibiotic course Gamma(shape = %.2f, scale = %.4f), mean = %.2f\n",
  cost_params_dist$sic_antibiotic_course$alpha,
  cost_params_dist$sic_antibiotic_course$beta,
  cost_params_dist$sic_antibiotic_course$alpha * cost_params_dist$sic_antibiotic_course$beta
))

icu_prop_params <- list(
  rih = gamma_mom_from_mean_ci(
    mean = 0.086, lower = 0.086 * 0.70, upper = 0.086 * 1.30
  ),
  sih = gamma_mom_from_mean_ci(
    mean = 0.040, lower = 0.040 * 0.70, upper = 0.040 * 1.30
  )
)
cat(sprintf(
  "RIH ICU proportion Gamma(shape = %.2f, scale = %.5f), mean = %.4f\n",
  icu_prop_params$rih$alpha,
  icu_prop_params$rih$beta,
  icu_prop_params$rih$alpha * icu_prop_params$rih$beta
))
cat(sprintf(
  "SIH ICU proportion Gamma(shape = %.2f, scale = %.5f), mean = %.4f\n",
  icu_prop_params$sih$alpha,
  icu_prop_params$sih$beta,
  icu_prop_params$sih$alpha * icu_prop_params$sih$beta
))

lp_params_dist <- list(
  mean_income_lt30 = 3549 * 12,
  mean_income_30_39 = 5660 * 12,
  mean_income_40_49 = 6376 * 12,
  mean_income_50_59 = 4437 * 12,
  mean_income_60 = 2543 * 12,
  mean_lfpr_20_24 = 0.557,
  mean_lfpr_25_29 = 0.850,
  mean_lfpr_30_34 = 0.898,
  mean_lfpr_35_39 = 0.882,
  mean_lfpr_40_44 = 0.868,
  mean_lfpr_45_49 = 0.835,
  mean_lfpr_50_54 = 0.817,
  mean_lfpr_55_59 = 0.747,
  mean_lfpr_60_64 = 0.635,
  mean_lfpr_65_69 = 0.490,
  mean_lfpr_70plus = 0.204,
  life_expectancy = life_expectancy_val
)

fit_lnorm_from_quantiles <- function(p, q) {
  opt <- optim(
    par = c(meanlog = log(median(q)), sdlog = 0.25),
    fn = function(par) {
      if (par[2] <= 0) return(Inf)
      pred <- qlnorm(p, meanlog = par[1], sdlog = par[2])
      sum((log(pred) - log(q))^2)
    },
    method = "L-BFGS-B",
    lower = c(-Inf, 1e-6),
    upper = c(Inf, Inf)
  )
  c(meanlog = opt$par[1], sdlog = opt$par[2])
}

fit_ri_death <- fit_lnorm_from_quantiles(c(0.25, 0.5, 0.75), c(53.0, 66.0, 74.0))
lognormal_params_ri_death <- list(meanlog = fit_ri_death[1], sdlog = fit_ri_death[2])
fit_si_death <- fit_lnorm_from_quantiles(c(0.25, 0.5, 0.75), c(56.1, 67.2, 76.5))
lognormal_params_si_death <- list(meanlog = fit_si_death[1], sdlog = fit_si_death[2])

beta_mom <- function(mu, lower, upper) {
  sigma <- (upper - lower) / (2 * 1.96)
  nu <- mu * (1 - mu) / sigma^2 - 1
  if (nu <= 0) stop("Invalid Beta moment parameters. Check mean and CI.")
  list(shape1 = mu * nu, shape2 = (1 - mu) * nu)
}
intervention_eff_params <- list(
  hosp_hh = beta_mom(0.5, 0.35, 0.65),
  comm_hh = beta_mom(0.458, 0.333, 0.619)
)

## ============================================================
## PART 2: CORE ODE MODEL AND HELPER FUNCTIONS
## ============================================================

annual_income_func <- function(age, lp) {
  a <- floor(age)
  if (a < 30) return(lp$mean_income_lt30)
  if (a < 40) return(lp$mean_income_30_39)
  if (a < 50) return(lp$mean_income_40_49)
  if (a < 60) return(lp$mean_income_50_59)
  return(lp$mean_income_60)
}

lfpr_func <- function(age, lp) {
  a <- floor(age)
  if (a < 20) return(0.0)
  if (a < 25) return(lp$mean_lfpr_20_24)
  if (a < 30) return(lp$mean_lfpr_25_29)
  if (a < 35) return(lp$mean_lfpr_30_34)
  if (a < 40) return(lp$mean_lfpr_35_39)
  if (a < 45) return(lp$mean_lfpr_40_44)
  if (a < 50) return(lp$mean_lfpr_45_49)
  if (a < 55) return(lp$mean_lfpr_50_54)
  if (a < 60) return(lp$mean_lfpr_55_59)
  if (a < 65) return(lp$mean_lfpr_60_64)
  if (a < 70) return(lp$mean_lfpr_65_69)
  return(lp$mean_lfpr_70plus)
}

calculate_lp_per_death <- function(age_at_death, lp, r = 0.03) {
  a <- max(0, age_at_death)
  le <- lp$life_expectancy
  if (a >= le) return(0)
  y0 <- floor(a)
  frac <- a - y0
  lp_total <- 0
  cur_share <- if (frac > 0) (1 - frac) else 1
  if (cur_share > 0) {
    lp_total <- lp_total + annual_income_func(y0, lp) * lfpr_func(y0, lp) * cur_share
  }
  le_floor <- floor(le)
  start_age <- y0 + 1
  if (start_age <= le_floor - 1) {
    for (t in start_age:(le_floor - 1)) {
      years_elapsed <- t - y0
      disc_factor <- 1 / (1 + r)^years_elapsed
      lp_total <- lp_total + annual_income_func(t, lp) * lfpr_func(t, lp) * disc_factor
    }
  }
  tail_share <- le - le_floor
  if (tail_share > 0) {
    years_elapsed <- le_floor - y0
    disc_factor <- 1 / (1 + r)^years_elapsed
    lp_total <- lp_total + annual_income_func(le_floor, lp) * lfpr_func(le_floor, lp) * tail_share * disc_factor
  }
  return(lp_total)
}

make_psa_seeds <- function(n_psa, base_seed = 123) {
  set.seed(base_seed)
  sample.int(.Machine$integer.max, n_psa)
}

constrain_flow <- function(flow_rate, population, dt = 1) {
  max_flow <- population / dt
  min(flow_rate, max_flow)
}

infection_model_ode <- function(t, y, params) {
  with(as.list(c(y, params)), {
    Sh <- max(0, y["Sh"]); RCh <- max(0, y["RCh"]); RIh <- max(0, y["RIh"]); SIh <- max(0, y["SIh"])
    Sc <- max(0, y["Sc"]); RCc <- max(0, y["RCc"]); RIc <- max(0, y["RIc"]); SIc <- max(0, y["SIc"])
    Nh_current <- max(1, Sh + RCh + RIh + SIh)
    Nc_current <- max(1, Sc + RCc + RIc + SIc)
    N_total_living <- Nh_current + Nc_current

    flow_beta_eh_Sh <- constrain_flow(beta_eh * Sh, Sh)
    flow_beta_hh_Sh <- constrain_flow(beta_hh * Sh * RCh / Nh_current, Sh)
    flow_beta_sh_sih_Sh <- constrain_flow(beta_sh_sih * Sh * SIh / Nh_current, Sh)
    flow_alpha_dis_Sh <- constrain_flow(alpha_dis * Sh, Sh)
    flow_mu_d_Sh <- constrain_flow(mu_d * Sh, Sh)

    flow_delta_rch_rih_RCh <- constrain_flow(delta_rch_rih * RCh, RCh)
    flow_gamma_rch_sh_RCh <- constrain_flow(gamma_rch_sh * RCh, RCh)
    flow_alpha_dis_RCh <- constrain_flow(alpha_dis * RCh, RCh)
    flow_delta_rch_sih_RCh <- constrain_flow(delta_rch_sih * RCh, RCh)
    flow_mu_d_RCh <- constrain_flow(mu_d * RCh, RCh)

    flow_gamma_rih_rch_RIh <- constrain_flow(gamma_rih_rch * RIh, RIh)
    flow_mu_rih_RIh <- constrain_flow(mu_rih * RIh, RIh)

    flow_gamma_sih_sh_SIh <- constrain_flow(gamma_sih_sh * SIh, SIh)
    flow_mu_sih_SIh <- constrain_flow(mu_sih * SIh, SIh)

    flow_beta_sc_ec_Sc <- constrain_flow(beta_sc_ec * Sc, Sc)
    flow_beta_sc_ac_Sc <- constrain_flow(beta_sc_ac * Sc * 0.33, Sc)
    flow_beta_hc_Sc <- constrain_flow(beta_hc * Sc * RCc / Nc_current, Sc)
    flow_beta_sc_sic_Sc <- constrain_flow(beta_sc_sic * Sc * SIc / Nc_current, Sc)
    flow_alpha_adm_Sc <- constrain_flow(alpha_adm * Sc, Sc)
    flow_mu_d_Sc <- constrain_flow(mu_d * Sc, Sc)

    flow_delta_rcc_ric_RCc <- constrain_flow(delta_rcc_ric * RCc, RCc)
    flow_delta_rcc_sic_RCc <- constrain_flow(delta_rcc_sic * RCc, RCc)
    flow_gamma_rcc_sc_RCc <- constrain_flow(gamma_rcc_sc * RCc, RCc)
    flow_alpha_adm_RCc <- constrain_flow(alpha_adm * RCc, RCc)
    flow_mu_d_RCc <- constrain_flow(mu_d * RCc, RCc)

    flow_gamma_ric_rcc_RIc <- constrain_flow(gamma_ric_rcc * RIc, RIc)
    flow_mu_ric_RIc <- constrain_flow(mu_ric * RIc, RIc)
    flow_alpha_adm_RIc <- constrain_flow(alpha_adm * RIc, RIc)

    flow_gamma_sic_sc_SIc <- constrain_flow(gamma_sic_sc * SIc, SIc)
    flow_mu_sic_SIc <- constrain_flow(mu_sic * SIc, SIc)
    flow_alpha_adm_SIc <- constrain_flow(alpha_adm * SIc, SIc)

    dSh_dt <- -flow_beta_eh_Sh - flow_beta_hh_Sh - flow_beta_sh_sih_Sh - flow_alpha_dis_Sh - flow_mu_d_Sh + flow_gamma_sih_sh_SIh + flow_gamma_rch_sh_RCh + flow_alpha_adm_Sc
    dRCh_dt <- -flow_delta_rch_rih_RCh - flow_gamma_rch_sh_RCh - flow_alpha_dis_RCh - flow_delta_rch_sih_RCh - flow_mu_d_RCh + flow_beta_eh_Sh + flow_beta_hh_Sh + flow_gamma_rih_rch_RIh + flow_alpha_adm_RCc
    dRIh_dt <- -flow_gamma_rih_rch_RIh - flow_mu_rih_RIh + flow_delta_rch_rih_RCh + flow_alpha_adm_RIc
    dSIh_dt <- -flow_gamma_sih_sh_SIh - flow_mu_sih_SIh + flow_beta_sh_sih_Sh + flow_delta_rch_sih_RCh + flow_alpha_adm_SIc
    dSc_dt <- -flow_beta_sc_ec_Sc - flow_beta_sc_ac_Sc - flow_beta_hc_Sc - flow_beta_sc_sic_Sc - flow_alpha_adm_Sc - flow_mu_d_Sc + flow_gamma_sic_sc_SIc + flow_gamma_rcc_sc_RCc + flow_alpha_dis_Sh + (mu_b + mu_m) * N_total_living
    dRCc_dt <- -flow_delta_rcc_ric_RCc - flow_delta_rcc_sic_RCc - flow_gamma_rcc_sc_RCc - flow_alpha_adm_RCc - flow_mu_d_RCc + flow_beta_sc_ec_Sc + flow_beta_sc_ac_Sc + flow_beta_hc_Sc + flow_gamma_ric_rcc_RIc + flow_alpha_dis_RCh
    dRIc_dt <- -flow_gamma_ric_rcc_RIc - flow_mu_ric_RIc - flow_alpha_adm_RIc + flow_delta_rcc_ric_RCc
    dSIc_dt <- -flow_gamma_sic_sc_SIc - flow_mu_sic_SIc - flow_alpha_adm_SIc + flow_beta_sc_sic_Sc + flow_delta_rcc_sic_RCc

    dDrh_dt <- flow_mu_rih_RIh
    dDsh_dt <- flow_mu_sih_SIh
    dDrc_dt <- flow_mu_ric_RIc
    dDsc_dt <- flow_mu_sic_SIc
    dC_RIh_dt <- flow_delta_rch_rih_RCh + flow_alpha_adm_RIc
    dC_SIh_dt <- flow_beta_sh_sih_Sh + flow_delta_rch_sih_RCh + flow_alpha_adm_SIc
    dC_RIc_dt <- flow_delta_rcc_ric_RCc
    dC_SIc_dt <- flow_beta_sc_sic_Sc + flow_delta_rcc_sic_RCc

    list(c(dSh_dt, dRCh_dt, dRIh_dt, dSIh_dt, dSc_dt, dRCc_dt, dRIc_dt, dSIc_dt,
           dDrh_dt, dDsh_dt, dDrc_dt, dDsc_dt, dC_RIh_dt, dC_SIh_dt, dC_RIc_dt, dC_SIc_dt))
  })
}

solve_model <- function(y_init, times, params) {
  out <- tryCatch(
    ode(y = y_init, times = times, func = infection_model_ode, parms = params,
        method = "radau", rtol = 1e-10, atol = 1e-12, hmax = 1),
    error = function(e) {
      ode(y = y_init, times = times, func = infection_model_ode, parms = params,
          method = "lsoda", rtol = 1e-8, atol = 1e-10, hmax = 1)
    }
  )
  as.data.frame(out)
}

calculate_annual_outcomes <- function(output_df, economic_start_year, economic_horizon_years) {
  if (is.null(output_df)) return(NULL)
  econ_start_time <- economic_start_year * 365
  econ_end_time <- (economic_start_year + economic_horizon_years) * 365
  econ_data <- output_df %>% filter(time >= econ_start_time & time <= econ_end_time)
  if (nrow(econ_data) < 2) return(NULL)
  annual_outcomes <- vector("list", nrow(econ_data) - 1)
  for (i in 2:nrow(econ_data)) {
    year_idx <- i - 1
    annual_outcomes[[i - 1]] <- data.frame(
      year = economic_start_year + year_idx,
      year_idx = year_idx,
      annual_res_deaths_h = econ_data$Drh[i] - econ_data$Drh[i - 1],
      annual_sus_deaths_h = econ_data$Dsh[i] - econ_data$Dsh[i - 1],
      annual_res_deaths_c = econ_data$Drc[i] - econ_data$Drc[i - 1],
      annual_sus_deaths_c = econ_data$Dsc[i] - econ_data$Dsc[i - 1],
      annual_res_cases_h = econ_data$C_RIh[i] - econ_data$C_RIh[i - 1],
      annual_sus_cases_h = econ_data$C_SIh[i] - econ_data$C_SIh[i - 1],
      annual_res_cases_c = econ_data$C_RIc[i] - econ_data$C_RIc[i - 1],
      annual_sus_cases_c = econ_data$C_SIc[i] - econ_data$C_SIc[i - 1]
    )
  }
  bind_rows(annual_outcomes)
}

# Adapted from total_cost_0623.R, original lines 459-719.
# Total burden and scenario comparisons only; attribution is in attributable_burden.R.
## ============================================================
## PART 3: PSA-CEA ENGINE
## ============================================================

env_n_cores <- Sys.getenv("TOTAL_COST_N_CORES", unset = NA_character_)
if (!is.na(env_n_cores) && nzchar(env_n_cores)) {
  n_cores <- as.integer(env_n_cores)
  if (is.na(n_cores) || n_cores < 1) {
    stop("TOTAL_COST_N_CORES must be a positive integer.")
  }
} else {
  detected_cores <- parallel::detectCores()
  if (is.na(detected_cores) || detected_cores < 2) {
    n_cores <- 1L
  } else {
    n_cores <- max(1L, detected_cores - 1L)
  }
}
if (n_cores > 1L) {
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)
  clusterEvalQ(cl, { library(deSolve); library(dplyr) })
} else {
  cl <- NULL
  registerDoSEQ()
  cat("Running PSA sequentially because n_cores = 1.\n")
}
cat("Parallel workers:", n_cores, "\n")

needed_objects <- c(
  "fixed_params", "full_sim_years", "economic_horizon_years", "discount_rate",
  "posterior_samples", "mcmc_pars", "los_epi_ranges", "cost_params_dist",
  "icu_prop_params", "lp_params_dist", "lognormal_params_ri_death",
  "lognormal_params_si_death", "intervention_eff_params", "solve_model",
  "infection_model_ode", "constrain_flow", "calculate_annual_outcomes",
  "calculate_lp_per_death", "annual_income_func", "lfpr_func"
)
psa_seeds <- make_psa_seeds(n_psa_samples, base_seed = psa_base_seed)
if (!is.null(cl)) {
  clusterExport(cl, c(needed_objects, "psa_seeds"))
}

pb <- txtProgressBar(max = n_psa_samples, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

psa_raw_list <- foreach(
  i = 1:n_psa_samples,
  .packages = c("deSolve", "dplyr"),
  .options.snow = opts,
  .errorhandling = "pass"
) %dopar% {
  set.seed(psa_seeds[i])

  p_iter <- fixed_params
  mcmc_row <- posterior_samples[sample(nrow(posterior_samples), 1), mcmc_pars, drop = FALSE]
  for (par in names(mcmc_row)) p_iter[[par]] <- as.numeric(mcmc_row[[par]])

  p_iter$gamma_rih_rch <- runif(1, los_epi_ranges$gamma_rih_rch$min, los_epi_ranges$gamma_rih_rch$max)
  p_iter$gamma_sih_sh <- runif(1, los_epi_ranges$gamma_sih_sh$min, los_epi_ranges$gamma_sih_sh$max)
  p_iter$mu_rih <- runif(1, los_epi_ranges$mu_rih$min, los_epi_ranges$mu_rih$max)
  p_iter$mu_sih <- runif(1, los_epi_ranges$mu_sih$min, los_epi_ranges$mu_sih$max)

  N <- 5180000
  Nh <- 37588 * 0.857
  Nc <- N - Nh
  RCh <- 0.124 * Nh
  RIh <- 0.119 * 0.314 * 0.377 * Nh
  SIh <- 0.119 * 0.314 * 0.56 * Nh
  Sh <- Nh - RCh - RIh - SIh
  RCc <- 0.257 * Nc
  RIc <- 0.001 * Nc
  SIc <- 0.01 * Nc
  Sc <- Nc - RCc - RIc - SIc
  y_init <- c(Sh = Sh, RCh = RCh, RIh = RIh, SIh = SIh, Sc = Sc, RCc = RCc, RIc = RIc, SIc = SIc,
              Drh = 0, Dsh = 0, Drc = 0, Dsc = 0, C_RIh = 0, C_SIh = 0, C_RIc = 0, C_SIc = 0)

  s_icu <- rgamma(1, shape = cost_params_dist$icu_bed_daily$alpha, scale = cost_params_dist$icu_bed_daily$beta)
  s_a <- rgamma(1, shape = cost_params_dist$a_class_bed_daily$alpha, scale = cost_params_dist$a_class_bed_daily$beta)
  s_b1 <- rgamma(1, shape = cost_params_dist$b1_class_bed_daily$alpha, scale = cost_params_dist$b1_class_bed_daily$beta)
  s_b2 <- rgamma(1, shape = cost_params_dist$b2_class_bed_daily$alpha, scale = cost_params_dist$b2_class_bed_daily$beta)
  s_c <- rgamma(1, shape = cost_params_dist$c_class_bed_daily$alpha, scale = cost_params_dist$c_class_bed_daily$beta)

  prop_rih_icu <- min(rgamma(1, shape = icu_prop_params$rih$alpha, scale = icu_prop_params$rih$beta), 1)
  prop_sih_icu <- min(rgamma(1, shape = icu_prop_params$sih$alpha, scale = icu_prop_params$sih$beta), 1)

  ward_bed_props <- cost_params_dist$ward_bed_proportions
  daily_general_ward <- s_a * ward_bed_props["A"] + s_b1 * ward_bed_props["B1"] +
    s_b2 * ward_bed_props["B2"] + s_c * ward_bed_props["C"]
  daily_h_ri <- s_icu * prop_rih_icu + daily_general_ward * (1 - prop_rih_icu)
  daily_h_si <- s_icu * prop_sih_icu + daily_general_ward * (1 - prop_sih_icu)

  s_comm_vis <- runif(1, cost_params_dist$outpatient_visit$min, cost_params_dist$outpatient_visit$max)
  s_comm_abx_r <- rgamma(1, shape = cost_params_dist$ric_antibiotic_course$alpha, scale = cost_params_dist$ric_antibiotic_course$beta)
  s_comm_abx_s <- rgamma(1, shape = cost_params_dist$sic_antibiotic_course$alpha, scale = cost_params_dist$sic_antibiotic_course$beta)
  s_cost_c_ri <- s_comm_vis + s_comm_abx_r
  s_cost_c_si <- s_comm_vis + s_comm_abx_s

  u_lp <- runif(1)
  age_death_ri <- qlnorm(u_lp, meanlog = lognormal_params_ri_death$meanlog, sdlog = lognormal_params_ri_death$sdlog)
  unit_lp_ri <- calculate_lp_per_death(age_death_ri, lp_params_dist, r = discount_rate)
  remaining_ly_ri <- max(0, lp_params_dist$life_expectancy - age_death_ri)
  ly_lost_val <- if (remaining_ly_ri > 0) {
    if (discount_rate == 0) remaining_ly_ri else (1 - (1 + discount_rate)^(-remaining_ly_ri)) / discount_rate
  } else 0

  age_death_si <- qlnorm(u_lp, meanlog = lognormal_params_si_death$meanlog, sdlog = lognormal_params_si_death$sdlog)
  unit_lp_si <- calculate_lp_per_death(age_death_si, lp_params_dist, r = discount_rate)
  remaining_ly_si <- max(0, lp_params_dist$life_expectancy - age_death_si)
  ly_lost_val_si <- if (remaining_ly_si > 0) {
    if (discount_rate == 0) remaining_ly_si else (1 - (1 + discount_rate)^(-remaining_ly_si)) / discount_rate
  } else 0

  draw_hosp_hh <- rbeta(1, intervention_eff_params$hosp_hh$shape1, intervention_eff_params$hosp_hh$shape2)
  draw_comm_hh <- rbeta(1, intervention_eff_params$comm_hh$shape1, intervention_eff_params$comm_hh$shape2)

  burnin_years <- full_sim_years - economic_horizon_years
  econ_horizon <- economic_horizon_years
  out_burnin <- tryCatch(solve_model(y_init, seq(0, 365 * burnin_years, by = 365), p_iter), error = function(e) NULL)
  if (is.null(out_burnin)) {
    NULL
  } else {
    final_row <- out_burnin[nrow(out_burnin), ]
    y_equilibrium <- c(
      Sh = as.numeric(final_row["Sh"]), RCh = as.numeric(final_row["RCh"]), RIh = as.numeric(final_row["RIh"]), SIh = as.numeric(final_row["SIh"]),
      Sc = as.numeric(final_row["Sc"]), RCc = as.numeric(final_row["RCc"]), RIc = as.numeric(final_row["RIc"]), SIc = as.numeric(final_row["SIc"]),
      Drh = 0, Dsh = 0, Drc = 0, Dsc = 0, C_RIh = 0, C_SIh = 0, C_RIc = 0, C_SIc = 0
    )
    scen_map <- list(
      Baseline = c(hosp_hh = 0, comm_hh = 0),
      S1_Hosp_HH = c(hosp_hh = draw_hosp_hh, comm_hh = 0),
      S2_Comm_HH = c(hosp_hh = 0, comm_hh = draw_comm_hh)
    )
    iter_results <- vector("list", length(scen_map))
    baseline_annual_list <- list()

    for (s_idx in seq_along(scen_map)) {
      sn <- names(scen_map)[s_idx]
      eff <- scen_map[[sn]]
      p_scen <- p_iter
      if (eff["hosp_hh"] > 0) {
        p_scen$beta_hh <- p_scen$beta_hh * (1 - eff["hosp_hh"])
        p_scen$beta_sh_sih <- p_scen$beta_sh_sih * (1 - eff["hosp_hh"])
      }
      if (eff["comm_hh"] > 0) {
        p_scen$beta_hc <- p_scen$beta_hc * (1 - eff["comm_hh"])
        p_scen$beta_sc_sic <- p_scen$beta_sc_sic * (1 - eff["comm_hh"])
      }
      out <- tryCatch(solve_model(y_equilibrium, seq(0, 365 * econ_horizon, by = 365), p_scen), error = function(e) NULL)
      if (is.null(out)) next
      epi <- tryCatch(calculate_annual_outcomes(out, 0, econ_horizon), error = function(e) NULL)
      if (is.null(epi) || nrow(epi) == 0) next

      los_rih <- 1 / (p_scen$gamma_rih_rch + p_scen$mu_rih)
      los_sih <- 1 / (p_scen$gamma_sih_sh + p_scen$mu_sih)

      epi_costed <- epi %>% mutate(
        disc = (1 + discount_rate)^(-(year - 1)),
        iter = i,
        scenario = sn,
        cost_hosp_ri_yr = annual_res_cases_h * (daily_h_ri * los_rih) * disc,
        cost_comm_ri_yr = annual_res_cases_c * s_cost_c_ri * disc,
        cost_hosp_si_yr = annual_sus_cases_h * (daily_h_si * los_sih) * disc,
        cost_comm_si_yr = annual_sus_cases_c * s_cost_c_si * disc,
        cost_lp_ri_h_yr = annual_res_deaths_h * unit_lp_ri * disc,
        cost_lp_ri_c_yr = annual_res_deaths_c * unit_lp_ri * disc,
        cost_lp_si_h_yr = annual_sus_deaths_h * unit_lp_si * disc,
        cost_lp_si_c_yr = annual_sus_deaths_c * unit_lp_si * disc
      )

      iter_results[[s_idx]] <- epi_costed %>% summarise(
        iter = i, scenario = sn,
        cost_hosp_ri = sum(cost_hosp_ri_yr),
        cost_comm_ri = sum(cost_comm_ri_yr),
        cost_lp_ri_h = sum(cost_lp_ri_h_yr),
        cost_lp_ri_c = sum(cost_lp_ri_c_yr),
        cost_lp_ri = cost_lp_ri_h + cost_lp_ri_c,
        cost_hosp_si = sum(cost_hosp_si_yr),
        cost_comm_si = sum(cost_comm_si_yr),
        cost_lp_si_h = sum(cost_lp_si_h_yr),
        cost_lp_si_c = sum(cost_lp_si_c_yr),
        cost_lp_si = cost_lp_si_h + cost_lp_si_c,
        cost_total = cost_hosp_ri + cost_comm_ri + cost_lp_ri + cost_hosp_si + cost_comm_si + cost_lp_si,
        ly_lost = sum((annual_res_deaths_h + annual_res_deaths_c) * ly_lost_val * disc +
                        (annual_sus_deaths_h + annual_sus_deaths_c) * ly_lost_val_si * disc),
        total_cases_rih = sum(annual_res_cases_h),
        total_cases_sih = sum(annual_sus_cases_h),
        total_cases_ric = sum(annual_res_cases_c),
        total_cases_sic = sum(annual_sus_cases_c),
        total_deaths_rh = sum(annual_res_deaths_h),
        total_deaths_sh = sum(annual_sus_deaths_h),
        total_deaths_rc = sum(annual_res_deaths_c),
        total_deaths_sc = sum(annual_sus_deaths_c),
        draw_beta_hh = p_iter$beta_hh,
        draw_beta_eh = p_iter$beta_eh,
        draw_beta_sh_sih = p_iter$beta_sh_sih,
        draw_beta_hc = p_iter$beta_hc,
        draw_beta_sc_ac = p_iter$beta_sc_ac,
        draw_beta_sc_ec = p_iter$beta_sc_ec,
        draw_delta_rcc_ric = p_iter$delta_rcc_ric,
        draw_gamma_rih_rch = p_iter$gamma_rih_rch,
        draw_gamma_sih_sh = p_iter$gamma_sih_sh,
        draw_mu_rih = p_iter$mu_rih,
        draw_mu_sih = p_iter$mu_sih,
        draw_los_rih = los_rih,
        draw_los_sih = los_sih,
        draw_prop_rih_icu = prop_rih_icu,
        draw_prop_sih_icu = prop_sih_icu,
        draw_daily_h_ri = daily_h_ri,
        draw_daily_h_si = daily_h_si,
        draw_s_comm_abx_r = s_comm_abx_r,
        draw_s_comm_abx_s = s_comm_abx_s,
        draw_s_cost_c = s_cost_c_ri,
        draw_s_cost_c_si = s_cost_c_si,
        draw_unit_lp = unit_lp_ri,
        draw_ly_lost_val = ly_lost_val,
        draw_unit_lp_si = unit_lp_si,
        draw_ly_lost_val_si = ly_lost_val_si,
        draw_hosp_hh = draw_hosp_hh,
        draw_comm_hh = draw_comm_hh,
        .groups = "drop"
      )

      if (sn == "Baseline") {
        baseline_annual_list[[length(baseline_annual_list) + 1]] <- epi_costed %>% select(
          iter, scenario, year, year_idx,
          annual_res_cases_h, annual_sus_cases_h, annual_res_cases_c, annual_sus_cases_c,
          annual_res_deaths_h, annual_sus_deaths_h, annual_res_deaths_c, annual_sus_deaths_c,
          cost_hosp_ri_yr, cost_comm_ri_yr, cost_lp_ri_h_yr, cost_lp_ri_c_yr,
          cost_hosp_si_yr, cost_comm_si_yr, cost_lp_si_h_yr, cost_lp_si_c_yr
        )
      }
    }
    list(summary = bind_rows(iter_results), baseline_annual = bind_rows(baseline_annual_list))
  }
}

close(pb)
if (!is.null(cl)) {
  stopCluster(cl)
}

psa_results <- bind_rows(lapply(psa_raw_list, function(x) if (is.list(x)) x$summary else NULL))
psa_baseline_annual <- bind_rows(lapply(psa_raw_list, function(x) if (is.list(x)) x$baseline_annual else NULL))


expected_n <- N_PSA
main <- as.data.frame(psa_results)
annual <- as.data.frame(psa_baseline_annual)
stopifnot(!any(startsWith(names(main), "attr_")), !any(startsWith(names(annual), "attr_")))
stopifnot(!anyDuplicated(main[c("iter", "scenario")]),
          !anyDuplicated(annual[c("iter", "year_idx")]),
          nrow(main) == 3L * expected_n, nrow(annual) == 5L * expected_n)
ids <- sort(main$iter[main$scenario == "Baseline"])
stopifnot(length(ids) == expected_n,
          setequal(unique(main$scenario), c("Baseline", "S1_Hosp_HH", "S2_Comm_HH")))
for (scenario_name in unique(main$scenario)) {
  stopifnot(identical(sort(main$iter[main$scenario == scenario_name]), ids))
}
stopifnot(identical(sort(unique(annual$iter)), ids), all(table(annual$iter) == 5L),
          setequal(unique(annual$year_idx), 1:5))
for (table in list(main, annual)) {
  stopifnot(all(vapply(table[vapply(table, is.numeric, logical(1))],
    function(x) all(is.finite(x)), logical(1))))
}
saveRDS(main, file.path(out, "psa_results.rds"), compress = "xz")
saveRDS(annual, file.path(out, "psa_baseline_annual.rds"), compress = "xz")
b <- main[main$scenario == "Baseline", ]
b <- b[order(b$iter), ]
baseline <- list(total_cost = b$cost_hosp_ri + b$cost_comm_ri + b$cost_lp_ri,
  hospital_direct = b$cost_hosp_ri, community_direct = b$cost_comm_ri,
  healthcare = b$cost_hosp_ri + b$cost_comm_ri, productivity = b$cost_lp_ri,
  hospital_infection_entries = b$total_cases_rih,
  community_infection_entries = b$total_cases_ric,
  infection_entries = b$total_cases_rih + b$total_cases_ric,
  deaths = b$total_deaths_rh + b$total_deaths_rc)
baseline_summary <- do.call(rbind, lapply(names(baseline), function(n) {
  counts <- grepl("entries|deaths", n)
  summary_row(baseline[[n]], n, if (counts) "count" else "million SGD", if (counts) 1 else 1e6)
}))
write_result(baseline_summary, "baseline_esbl_5yr", out)
ann <- annual
ann$calendar_year <- 2025L + ann$year_idx
ann$esbl_total_cost <- ann$cost_hosp_ri_yr + ann$cost_comm_ri_yr + ann$cost_lp_ri_h_yr + ann$cost_lp_ri_c_yr
ann$esbl_productivity <- ann$cost_lp_ri_h_yr + ann$cost_lp_ri_c_yr
ann_cols <- c("annual_res_cases_h", "annual_res_cases_c", "annual_sus_cases_h", "annual_sus_cases_c",
  "annual_res_deaths_h", "annual_res_deaths_c", "annual_sus_deaths_h", "annual_sus_deaths_c",
  "cost_hosp_ri_yr", "cost_comm_ri_yr", "esbl_productivity", "esbl_total_cost")
ann_summary <- do.call(rbind, lapply(split(ann, ann$calendar_year), function(z) {
  do.call(rbind, lapply(ann_cols, function(n) {
    counts <- grepl("cases|deaths", n)
    cbind(calendar_year = z$calendar_year[1], summary_row(z[[n]], n,
      if (counts) "count" else "million SGD", if (counts) 1 else 1e6))
  }))
}))
write_result(ann_summary, "baseline_annual", out)
delta_rows <- list()
for (sc in c("S1_Hosp_HH", "S2_Comm_HH")) {
  s <- main[main$scenario == sc, ]
  s <- s[match(b$iter, s$iter), ]
  for (org in c("RI", "SI", "All")) {
    values <- function(x) {
      ri <- list(cases = x$total_cases_rih + x$total_cases_ric, deaths = x$total_deaths_rh + x$total_deaths_rc,
        healthcare = x$cost_hosp_ri + x$cost_comm_ri, productivity = x$cost_lp_ri)
      si <- list(cases = x$total_cases_sih + x$total_cases_sic, deaths = x$total_deaths_sh + x$total_deaths_sc,
        healthcare = x$cost_hosp_si + x$cost_comm_si, productivity = x$cost_lp_si)
      v <- switch(org, RI = ri, SI = si, All = Map(`+`, ri, si))
      v$total <- v$healthcare + v$productivity
      v
    }
    diff <- Map(`-`, values(s), values(b))
    for (n in names(diff)) {
      counts <- n %in% c("cases", "deaths")
      delta_rows[[length(delta_rows) + 1L]] <- cbind(scenario = sc, phenotype = org,
        contrast = "intervention minus baseline", summary_row(diff[[n]], n,
          if (counts) "count" else "million SGD", if (counts) 1 else 1e6))
    }
  }
}
intervention_summary <- do.call(rbind, delta_rows)
write_result(intervention_summary, "intervention_deltas_paired", out)

# Print the same signs, rounding and columns as manuscript Table 3.
format_table3_result <- function(scenario_name, phenotype_name, outcome_name) {
  r <- intervention_summary[intervention_summary$scenario == scenario_name &
    intervention_summary$phenotype == phenotype_name &
    intervention_summary$outcome == outcome_name, ]
  stopifnot(nrow(r) == 1L)
  v <- c(r$mean, r$lower_95_UI, r$upper_95_UI)
  counts <- outcome_name %in% c("cases", "deaths")
  if (counts) v <- c(-v[1], -v[3], -v[2])
  s <- formatC(v, format = "f", digits = if (counts) 0 else 2, big.mark = ",")
  if (phenotype_name == "All") s[1] else paste0(s[1], " (", s[2], " to ", s[3], ")")
}
table3_outcomes <- c(cases = "Infection-state entries averted, n",
  deaths = "Deaths averted, n", total = "Change in total costs, million SGD",
  healthcare = "Change in healthcare costs, million SGD",
  productivity = "Change in productivity costs, million SGD")
table3_display <- do.call(rbind, lapply(c("S1_Hosp_HH", "S2_Comm_HH"), function(sc) {
  do.call(rbind, lapply(names(table3_outcomes), function(n) data.frame(
    Intervention = if (sc == "S1_Hosp_HH") "Hospital HH" else "Community HH",
    Outcome = unname(table3_outcomes[n]),
    ESBL = format_table3_result(sc, "RI", n),
    Non_ESBL = format_table3_result(sc, "SI", n),
    All_Enterobacterales = format_table3_result(sc, "All", n),
    check.names = FALSE)))
}))
cat("\n=== TABLE 3: PROJECTED EFFECTS OF INTERVENTIONS, 2026-2030 ===\n")
cat("Fresh R simulation; 1,000 paired PSA iterations using the fit7 joint posterior.\n")
cat("ESBL/non-ESBL: mean (95% uncertainty interval); All Enterobacterales: mean only.\n")
cat("Counts: baseline minus intervention (averted). Costs: intervention minus baseline; implementation costs excluded.\n")
cat("Outcomes span hospital and community; counts are undiscounted; costs are discounted at 3% annually, with 2026 undiscounted.\n")
print(table3_display, row.names = FALSE, right = FALSE, quote = FALSE)
cat("Full-precision paired differences: intervention_deltas_paired.csv\n")
cat("=== END TABLE 3 ===\n\n")


print(baseline_summary, row.names = FALSE)
psa_baseline_annual_print <- transform(psa_baseline_annual, calendar_year = 2025L + year_idx)
# Extracted from total_cost_0623.R, original lines 1256-1323.
## Hospital vs community cases, 2026-2030
## ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

ri_cases_annual_plot <- psa_baseline_annual_print %>%
  group_by(calendar_year) %>%
  summarise(
    Hospital = mean(annual_res_cases_h, na.rm = TRUE),
    Community = mean(annual_res_cases_c, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Total = Hospital + Community,
    Hospital_pct = Hospital / Total,
    Community_pct = Community / Total
  ) %>%
  pivot_longer(
    cols = c(Hospital, Community),
    names_to = "Setting",
    values_to = "Cases"
  ) %>%
  mutate(
    Setting = factor(Setting, levels = c("Community", "Hospital"))
  )

p_ri_bar <- ggplot(
  ri_cases_annual_plot,
  aes(x = factor(calendar_year), y = Cases, fill = Setting)
) +
  geom_col(
    color = "black",
    linewidth = 0.3,
    width = 0.65
  ) +
  geom_text(
    aes(label = comma(round(Cases))),
    position = position_stack(vjust = 0.5),
    size = 3.2
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_fill_manual(
    values = c(
      "Community" = "#f4a582",
      "Hospital" = "#ca0020"
    )
  ) +
  labs(
    title = "Annual Burden of ESBL-Producing Enterobacterales Infections",
    subtitle = "Projected baseline cases by setting, Singapore, 2026-2030",
    x = NULL,
    y = "Annual cases",
    fill = "Setting",
    caption = "Bars show mean projected annual ESBL-producing Enterobacterales infection cases."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

if (interactive()) print(p_ri_bar)

ggsave(file.path(out, "figure2_annual_esbl_entries.png"), p_ri_bar, width = 8, height = 5, dpi = 300, bg = "white")
ggsave(file.path(out, "figure2_annual_esbl_entries.pdf"), p_ri_bar, width = 8, height = 5, bg = "white")


save(psa_results, psa_baseline_annual, posterior_samples, baseline_summary,
     intervention_summary, table3_display, file = file.path(out, "economic_fit7.RData"))
cat("Results saved to:", normalizePath(out), "\n")
