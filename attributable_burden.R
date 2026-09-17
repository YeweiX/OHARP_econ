# Resistance-attributable economic burden: designated 317.29 route
# Based on econ_incremental_0514.R; no other R script needs to be sourced.
# Change RUN_MODE to "fresh" to rerun the transmission/economic PSA.
RUN_MODE <- Sys.getenv("OHARP_MODE", "reported")
N_PSA <- as.integer(Sys.getenv("OHARP_N", "1000"))
SEED <- 123L
OUTPUT_FOLDER <- file.path("outputs", paste0("attributable_burden_", RUN_MODE))
stopifnot(RUN_MODE %in% c("reported", "fresh"), !is.na(N_PSA), N_PSA > 0L)

script_directory <- function() {
  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile)) return(dirname(normalizePath(frame$ofile, mustWork = TRUE)))
  }
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(arg) == 1L) {
    path <- sub("^--file=", "", arg)
    if (!file.exists(path)) path <- gsub("~+~", " ", path, fixed = TRUE)
    return(dirname(normalizePath(path, mustWork = TRUE)))
  }
  if (dir.exists("data")) return(normalizePath(getwd()))
  stop("Use source() on this file, or run it with Rscript.")
}
root <- script_directory()
previous_directory <- getwd()
setwd(root)
options(digits = 17)
override_output <- Sys.getenv("OHARP_OUTPUT", "")
out <- if (nzchar(override_output)) override_output else OUTPUT_FOLDER
if (file.exists(file.path(out, "RUN_COMPLETE.txt"))) {
  stop("This output folder already contains a completed run. Change OUTPUT_FOLDER at the top of this script.")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)
require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install these packages first: ", paste(missing, collapse = ", "))
}
assert_close <- function(actual, expected, label, rtol = 1e-10, atol = 1e-5) {
  if (length(actual) != length(expected) || any(!is.finite(actual)) ||
      any(!is.finite(expected)) ||
      any(abs(actual - expected) > atol + rtol * abs(expected))) {
    stop("Validation failed: ", label)
  }
  invisible(TRUE)
}

summary_row <- function(x, outcome, unit = "million SGD", scale = 1e6) {
  stopifnot(length(x) > 0L, all(is.finite(x)))
  data.frame(outcome = outcome, unit = unit, n_draws = length(x),
    mean = mean(x) / scale,
    lower_95_UI = unname(quantile(x, 0.025, type = 7)) / scale,
    upper_95_UI = unname(quantile(x, 0.975, type = 7)) / scale)
}

sum_by_iteration <- function(x, cols) {
  aggregate(x[cols], by = list(iter = x$iter), FUN = sum)
}

write_result <- function(x, name, out) {
  write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
}

rebuild_attribution <- function(x, parameters) {
  # The designated positive-premium route from econ_incremental_0514.R.
  r <- parameters$discount_rate
  p <- parameters$fixed_params
  out <- x[c("iter", "year_idx")]
  out$calendar_year <- 2025L + out$year_idx
  dh <- pmax(0, x$draw_daily_h_ri / (x$draw_gamma_rih_rch + x$draw_mu_rih) -
               x$draw_daily_h_si / (x$draw_gamma_sih_sh + x$draw_mu_sih))
  dc <- pmax(0, x$draw_s_comm_abx_r - x$draw_s_comm_abx_s)
  ph <- pmax(0, x$draw_mu_rih / (x$draw_gamma_rih_rch + x$draw_mu_rih) -
               x$draw_mu_sih / (x$draw_gamma_sih_sh + x$draw_mu_sih))
  pc <- pmax(0, p$mu_ric / (p$mu_ric + p$gamma_ric_rcc + p$alpha_adm) -
               p$mu_sic / (p$mu_sic + p$gamma_sic_sc + p$alpha_adm))
  disc <- (1 + r)^(-(x$year_idx - 1L))
  out$direct_h <- x$annual_res_cases_h * dh * disc
  out$direct_c <- x$annual_res_cases_c * dc * disc
  out$deaths_h <- x$annual_res_cases_h * ph
  out$deaths_c <- x$annual_res_cases_c * pc
  out$productivity_h <- out$deaths_h * x$unit_lp_ri * disc
  out$productivity_c <- out$deaths_c * x$unit_lp_ri * disc
  out$direct <- out$direct_h + out$direct_c
  out$productivity <- out$productivity_h + out$productivity_c
  out$total <- out$direct + out$productivity
  out$deaths <- out$deaths_h + out$deaths_c
  checks <- list(direct_h = "attr_direct_h", direct_c = "attr_direct_c",
    deaths_h = "attr_deaths_h", deaths_c = "attr_deaths_c",
    productivity_h = "attr_lp_h", productivity_c = "attr_lp_c", total = "attr_total")
  for (n in names(checks)) assert_close(out[[n]], x[[checks[[n]]]], paste("attribution", n))
  assert_close(disc, x$disc_factor, "incremental discounting")
  out
}


parameters <- readRDS("data/reported_input_parameters.rds")
if (RUN_MODE == "fresh") {
  require_packages(c("deSolve", "dplyr", "foreach", "doSNOW", "purrr",
                     "tidyr", "ggplot2", "scales", "doParallel"))
  Sys.setenv(N_PSA_SAMPLES = N_PSA, PSA_BASE_SEED = SEED, TOTAL_COST_N_CORES = 1L)
  posterior_samples <- readRDS("data/posterior_samples.rds")
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

  for (key in names(parameters)) {
    stopifnot(isTRUE(all.equal(get(key), parameters[[key]], check.attributes = FALSE, tolerance = 1e-10)))
  }
  # Extracted from econ_incremental_0514.R, original lines 46-627.
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(doParallel)
  library(foreach)

  ## ============================================================
  ## PART 0: Helper distributions
  ## ============================================================

  gamma_mom_from_mean_ci <- function(mean, lower, upper) {
    if (mean <= 0 || lower <= 0 || upper <= lower) {
      stop("Invalid mean/lower/upper in gamma_mom_from_mean_ci().")
    }
    sd <- (upper - lower) / (2 * 1.96)
    shape <- (mean / sd)^2
    scale <- sd^2 / mean
    list(alpha = shape, beta = scale)
  }

  # Use revised LOS-related uncertainty ranges if not already defined.
  if (!exists("los_epi_ranges")) {
    los_epi_ranges <- list(
      gamma_sih_sh  = list(min = 0.09090909, max = 0.20000000),
      gamma_rih_rch = list(min = 0.06666667, max = 0.16666667),
      mu_sih        = list(min = 0.00233333, max = 0.00390000),
      mu_rih        = list(min = 0.00254333, max = 0.00554667)
    )
  }

  # Use revised ICU proportion distributions if not already defined.
  if (!exists("icu_prop_params")) {
    icu_prop_params <- list(
      rih = gamma_mom_from_mean_ci(
        mean  = 0.086,
        lower = 0.086 * 0.70,
        upper = 0.086 * 1.30
      ),
      sih = gamma_mom_from_mean_ci(
        mean  = 0.040,
        lower = 0.040 * 0.70,
        upper = 0.040 * 1.30
      )
    )
  }

  # Update community antibiotic distributions if the current cost object is old.
  # This overwrites only the community antibiotic course distributions.
  if (exists("cost_params_dist")) {
    cost_params_dist$ric_antibiotic_course <- gamma_mom_from_mean_ci(
      mean  = 18.20,
      lower = 18.20 * 0.70,
      upper = 18.20 * 1.30
    )

    cost_params_dist$sic_antibiotic_course <- gamma_mom_from_mean_ci(
      mean  = 12.04,
      lower = 12.04 * 0.70,
      upper = 12.04 * 1.30
    )

    # Standardise ward proportion naming for compatibility.
    if (!"ward_bed_proportions" %in% names(cost_params_dist) &&
        "sih_bed_proportions" %in% names(cost_params_dist)) {
      cost_params_dist$ward_bed_proportions <- cost_params_dist$sih_bed_proportions
    }
  }

  ## ============================================================
  ## PART 1: Parallel PSA — GRAM-Style Attributable Premium
  ## ============================================================

  run_incremental_psa <- function(n_psa_samples = 1000,
                                  posterior_samples,
                                  fixed_params,
                                  cost_params_dist,
                                  lp_params_dist,
                                  lognormal_params_ri_death,
                                  lognormal_params_si_death,
                                  los_epi_ranges,
                                  icu_prop_params,
                                  full_sim_years = 20,
                                  economic_horizon_years = 5,
                                  discount_rate = 0.03,
                                  base_seed = 123L,
                                  n_cores = 1L,
                                  mcmc_pars = c(
                                    "beta_hh",
                                    "beta_eh",
                                    "beta_hc",
                                    "beta_sc_ac",
                                    "beta_sc_ec",
                                    "beta_sh_sih",
                                    "delta_rcc_ric"
                                  )) {

    required_cost_keys <- c(
      "icu_bed_daily",
      "a_class_bed_daily",
      "b1_class_bed_daily",
      "b2_class_bed_daily",
      "c_class_bed_daily",
      "ward_bed_proportions",
      "ric_antibiotic_course",
      "sic_antibiotic_course",
      "outpatient_visit"
    )

    missing_cost_keys <- setdiff(required_cost_keys, names(cost_params_dist))
    if (length(missing_cost_keys) > 0) {
      stop("cost_params_dist is missing required keys: ",
           paste(missing_cost_keys, collapse = ", "))
    }

    missing_mcmc_cols <- setdiff(mcmc_pars, colnames(posterior_samples))
    if (length(missing_mcmc_cols) > 0) {
      stop("posterior_samples is missing columns: ",
           paste(missing_mcmc_cols, collapse = ", "))
    }

    seeds <- make_psa_seeds(n_psa_samples, base_seed = base_seed)

    # ── Parallel setup ──────────────────────────────────────────
    stopifnot(n_cores >= 1L, n_psa_samples >= 1L)
    cl <- NULL
    if (n_cores > 1L) {
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
    } else {
      foreach::registerDoSEQ()
    }

    on.exit({
      if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)

    cat(
      "Running revised GRAM-style incremental PSA (",
      n_psa_samples, " iterations, ", n_cores, " cores)...\n",
      sep = ""
    )

    if (!is.null(cl)) {
    parallel::clusterExport(
      cl,
      c(
        "fixed_params", "posterior_samples", "mcmc_pars", "los_epi_ranges",
        "icu_prop_params", "cost_params_dist", "lp_params_dist",
        "lognormal_params_ri_death", "lognormal_params_si_death",
        "full_sim_years", "economic_horizon_years", "discount_rate", "seeds",
        "solve_model", "infection_model_ode", "constrain_flow",
        "calculate_annual_outcomes", "make_psa_seeds"
      ),
      envir = environment()
    )

    parallel::clusterEvalQ(cl, {
      library(deSolve)
      library(dplyr)
    })
    }

    # ── Parallel loop ───────────────────────────────────────────
    results_list <- foreach(
      i = 1:n_psa_samples,
      .packages = c("deSolve", "dplyr"),
      .errorhandling = "pass"
    ) %dopar% {

      set.seed(seeds[i])

      # =========================================================
      # A. Epidemiological parameter sampling
      # =========================================================

      params_iter <- fixed_params

      # Sample calibrated posterior parameters.
      rand_idx <- sample(nrow(posterior_samples), 1)
      mcmc_row <- posterior_samples[rand_idx, mcmc_pars, drop = FALSE]

      for (par_name in names(mcmc_row)) {
        if (par_name %in% names(params_iter)) {
          params_iter[[par_name]] <- as.numeric(mcmc_row[[par_name]])
        }
      }

      # Sample LOS-related hospital recovery and mortality parameters BEFORE ODE.
      # These same values are reused in LOS/costing.
      params_iter$gamma_rih_rch <- runif(
        1,
        los_epi_ranges$gamma_rih_rch$min,
        los_epi_ranges$gamma_rih_rch$max
      )

      params_iter$gamma_sih_sh <- runif(
        1,
        los_epi_ranges$gamma_sih_sh$min,
        los_epi_ranges$gamma_sih_sh$max
      )

      params_iter$mu_rih <- runif(
        1,
        los_epi_ranges$mu_rih$min,
        los_epi_ranges$mu_rih$max
      )

      params_iter$mu_sih <- runif(
        1,
        los_epi_ranges$mu_sih$min,
        los_epi_ranges$mu_sih$max
      )

      # =========================================================
      # B. Run baseline ODE
      # =========================================================

      N <- 5.18e6
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

      y_init <- c(
        Sh = Sh, RCh = RCh, RIh = RIh, SIh = SIh,
        Sc = Sc, RCc = RCc, RIc = RIc, SIc = SIc,
        Drh = 0, Dsh = 0, Drc = 0, Dsc = 0,
        C_RIh = 0, C_SIh = 0, C_RIc = 0, C_SIc = 0
      )

      times <- seq(0, 365 * full_sim_years, by = 365)

      out_ode <- tryCatch(
        solve_model(y_init, times, params_iter),
        error = function(e) NULL
      )

      if (is.null(out_ode)) {
        return(NULL)
      }

      econ_start_year <- full_sim_years - economic_horizon_years

      epi_annual <- calculate_annual_outcomes(
        out_ode,
        econ_start_year,
        economic_horizon_years
      )

      if (is.null(epi_annual) || nrow(epi_annual) == 0) {
        return(NULL)
      }

      # =========================================================
      # C. Cost parameter sampling
      # =========================================================

      # ---------------------------------------------------------
      # C1. LOS
      # ---------------------------------------------------------
      # LOS is the model-implied residence time in hospital infection
      # compartments and uses the same sampled parameters as the ODE.
      los_ri_raw <- 1 / (params_iter$gamma_rih_rch + params_iter$mu_rih)
      los_si <- 1 / (params_iter$gamma_sih_sh + params_iter$mu_sih)

      # For GRAM-style attributable burden, RI direct premium is constrained
      # to be non-negative below, so we do not need to force LOS_RI >= LOS_SI.
      los_ri <- los_ri_raw

      # ---------------------------------------------------------
      # C2. Bed-day costs
      # ---------------------------------------------------------

      s_icu <- rgamma(
        1,
        shape = cost_params_dist$icu_bed_daily$alpha,
        scale = cost_params_dist$icu_bed_daily$beta
      )

      s_a <- rgamma(
        1,
        shape = cost_params_dist$a_class_bed_daily$alpha,
        scale = cost_params_dist$a_class_bed_daily$beta
      )

      s_b1 <- rgamma(
        1,
        shape = cost_params_dist$b1_class_bed_daily$alpha,
        scale = cost_params_dist$b1_class_bed_daily$beta
      )

      s_b2 <- rgamma(
        1,
        shape = cost_params_dist$b2_class_bed_daily$alpha,
        scale = cost_params_dist$b2_class_bed_daily$beta
      )

      s_c <- rgamma(
        1,
        shape = cost_params_dist$c_class_bed_daily$alpha,
        scale = cost_params_dist$c_class_bed_daily$beta
      )

      ward_props <- cost_params_dist$ward_bed_proportions

      s_ward_avg <-
        s_a  * ward_props["A"] +
        s_b1 * ward_props["B1"] +
        s_b2 * ward_props["B2"] +
        s_c  * ward_props["C"]

      # ICU proportions sampled as Gamma distributions using mean ±30% bounds.
      prop_rih_in_icu <- rgamma(
        1,
        shape = icu_prop_params$rih$alpha,
        scale = icu_prop_params$rih$beta
      )

      prop_sih_in_icu <- rgamma(
        1,
        shape = icu_prop_params$sih$alpha,
        scale = icu_prop_params$sih$beta
      )

      # Safety cap.
      prop_rih_in_icu <- min(prop_rih_in_icu, 1)
      prop_sih_in_icu <- min(prop_sih_in_icu, 1)

      # Hospital direct cost: bed-day cost only.
      s_bed_expected_ri <-
        s_icu * prop_rih_in_icu +
        s_ward_avg * (1 - prop_rih_in_icu)

      s_bed_expected_si <-
        s_icu * prop_sih_in_icu +
        s_ward_avg * (1 - prop_sih_in_icu)

      s_daily_h_ri <- s_bed_expected_ri
      s_daily_h_si <- s_bed_expected_si

      # ---------------------------------------------------------
      # C3. Community unit costs
      # ---------------------------------------------------------
      # Community direct cost: outpatient visit + antibiotic course only.
      # Community diagnostic cost removed.

      s_comm_vis <- runif(
        1,
        cost_params_dist$outpatient_visit$min,
        cost_params_dist$outpatient_visit$max
      )

      s_comm_abx_r <- rgamma(
        1,
        shape = cost_params_dist$ric_antibiotic_course$alpha,
        scale = cost_params_dist$ric_antibiotic_course$beta
      )

      s_comm_abx_s <- rgamma(
        1,
        shape = cost_params_dist$sic_antibiotic_course$alpha,
        scale = cost_params_dist$sic_antibiotic_course$beta
      )

      s_cost_c_ri <- s_comm_vis + s_comm_abx_r
      s_cost_c_si <- s_comm_vis + s_comm_abx_s

      # =========================================================
      # D. Productivity loss per death
      # =========================================================

      u_lp <- runif(1)

      age_ri <- qlnorm(
        u_lp,
        lognormal_params_ri_death$meanlog,
        lognormal_params_ri_death$sdlog
      )

      calc_lp_inline <- function(age_at_death, lp, r) {
        a <- max(0, age_at_death)
        le <- lp$life_expectancy

        if (a >= le) {
          return(0)
        }

        get_income <- function(age) {
          a2 <- floor(age)
          if (a2 < 30) return(lp$mean_income_lt30)
          if (a2 < 40) return(lp$mean_income_30_39)
          if (a2 < 50) return(lp$mean_income_40_49)
          if (a2 < 60) return(lp$mean_income_50_59)
          return(lp$mean_income_60)
        }

        get_lfpr <- function(age) {
          a2 <- floor(age)
          if (a2 < 20) return(0)
          if (a2 < 25) return(lp$mean_lfpr_20_24)
          if (a2 < 30) return(lp$mean_lfpr_25_29)
          if (a2 < 35) return(lp$mean_lfpr_30_34)
          if (a2 < 40) return(lp$mean_lfpr_35_39)
          if (a2 < 45) return(lp$mean_lfpr_40_44)
          if (a2 < 50) return(lp$mean_lfpr_45_49)
          if (a2 < 55) return(lp$mean_lfpr_50_54)
          if (a2 < 60) return(lp$mean_lfpr_55_59)
          if (a2 < 65) return(lp$mean_lfpr_60_64)
          if (a2 < 70) return(lp$mean_lfpr_65_69)
          return(lp$mean_lfpr_70plus)
        }

        y0 <- floor(a)
        frac <- a - y0
        lp_total <- 0

        cur_share <- if (frac > 0) (1 - frac) else 1

        if (cur_share > 0) {
          lp_total <- lp_total +
            get_income(y0) * get_lfpr(y0) * cur_share
        }

        le_floor <- floor(le)

        if (y0 + 1 <= le_floor - 1) {
          for (t in (y0 + 1):(le_floor - 1)) {
            yrs <- t - y0
            df <- if (r == 0) 1 else 1 / (1 + r)^yrs

            lp_total <- lp_total +
              get_income(t) * get_lfpr(t) * df
          }
        }

        tail_share <- le - le_floor

        if (tail_share > 0) {
          yrs <- le_floor - y0
          df <- if (r == 0) 1 else 1 / (1 + r)^yrs

          lp_total <- lp_total +
            get_income(le_floor) * get_lfpr(le_floor) * tail_share * df
        }

        return(lp_total)
      }

      unit_lp_ri <- calc_lp_inline(age_ri, lp_params_dist, discount_rate)

      # =========================================================
      # E. Per-case attributable direct cost premiums
      # =========================================================
      # Revised hospital unit cost = bed-day cost only.
      # No hospital antibiotics, diagnostics, routine clinical tests, or procedures.

      unit_cost_h_ri <- s_daily_h_ri * los_ri
      unit_cost_h_si <- s_daily_h_si * los_si

      premium_direct_h_raw <- unit_cost_h_ri - unit_cost_h_si
      premium_direct_c_raw <- s_cost_c_ri - s_cost_c_si

      # GRAM-style attributable burden should be non-negative.
      premium_direct_h <- max(0, premium_direct_h_raw)
      premium_direct_c <- max(0, premium_direct_c_raw)

      # =========================================================
      # F. Excess mortality probability attributable to resistance
      # =========================================================

      # Hospital death probabilities.
      # RIh/SIh leave through recovery or death.
      prob_death_rih <- params_iter$mu_rih /
        (params_iter$mu_rih + params_iter$gamma_rih_rch)

      prob_death_sih <- params_iter$mu_sih /
        (params_iter$mu_sih + params_iter$gamma_sih_sh)

      # Community death probabilities.
      # RIc/SIc leave through recovery, death, or admission.
      prob_death_ric <- params_iter$mu_ric /
        (params_iter$mu_ric + params_iter$gamma_ric_rcc + params_iter$alpha_adm)

      prob_death_sic <- params_iter$mu_sic /
        (params_iter$mu_sic + params_iter$gamma_sic_sc + params_iter$alpha_adm)

      excess_prob_death_h_raw <- prob_death_rih - prob_death_sih
      excess_prob_death_c_raw <- prob_death_ric - prob_death_sic

      # Attributable burden should not be negative.
      excess_prob_death_h <- max(0, excess_prob_death_h_raw)
      excess_prob_death_c <- max(0, excess_prob_death_c_raw)

      # =========================================================
      # G. Apply attributable premiums to annual RI cases
      # =========================================================

      psa_df <- epi_annual %>%
        mutate(
          iter = i,

          # year_idx = 1 has no additional discount.
          disc_factor = (1 + discount_rate)^(-(year_idx - 1)),

          # Direct attributable cost.
          attr_direct_h = annual_res_cases_h * premium_direct_h * disc_factor,
          attr_direct_c = annual_res_cases_c * premium_direct_c * disc_factor,

          # Excess deaths attributable to resistance.
          attr_deaths_h = annual_res_cases_h * excess_prob_death_h,
          attr_deaths_c = annual_res_cases_c * excess_prob_death_c,

          # Mortality-related excess productivity loss.
          attr_lp_h = attr_deaths_h * unit_lp_ri * disc_factor,
          attr_lp_c = attr_deaths_c * unit_lp_ri * disc_factor,

          # Total attributable burden.
          attr_direct = attr_direct_h + attr_direct_c,
          attr_lp = attr_lp_h + attr_lp_c,
          attr_total = attr_direct + attr_lp,

          # Epi counts for reference.
          new_ri_h = annual_res_cases_h,
          new_ri_c = annual_res_cases_c,
          new_ri = annual_res_cases_h + annual_res_cases_c,

          deaths_ri_h = annual_res_deaths_h,
          deaths_ri_c = annual_res_deaths_c,
          deaths_ri = annual_res_deaths_h + annual_res_deaths_c,

          attr_deaths = attr_deaths_h + attr_deaths_c,

          # Diagnostics.
          los_ri_raw = los_ri_raw,
          los_ri = los_ri,
          los_si = los_si,

          prop_rih_in_icu = prop_rih_in_icu,
          prop_sih_in_icu = prop_sih_in_icu,

          unit_cost_h_ri = unit_cost_h_ri,
          unit_cost_h_si = unit_cost_h_si,

          premium_direct_h_raw = premium_direct_h_raw,
          premium_direct_c_raw = premium_direct_c_raw,
          premium_direct_h = premium_direct_h,
          premium_direct_c = premium_direct_c,

          prob_death_rih = prob_death_rih,
          prob_death_sih = prob_death_sih,
          prob_death_ric = prob_death_ric,
          prob_death_sic = prob_death_sic,

          excess_prob_death_h_raw = excess_prob_death_h_raw,
          excess_prob_death_c_raw = excess_prob_death_c_raw,
          excess_prob_death_h = excess_prob_death_h,
          excess_prob_death_c = excess_prob_death_c,

          age_ri = age_ri,
          unit_lp_ri = unit_lp_ri,

          # Draws for diagnostics / EVPPI if needed.
          draw_gamma_rih_rch = params_iter$gamma_rih_rch,
          draw_gamma_sih_sh = params_iter$gamma_sih_sh,
          draw_mu_rih = params_iter$mu_rih,
          draw_mu_sih = params_iter$mu_sih,
          draw_s_comm_abx_r = s_comm_abx_r,
          draw_s_comm_abx_s = s_comm_abx_s,
          draw_daily_h_ri = s_daily_h_ri,
          draw_daily_h_si = s_daily_h_si
        )

      psa_df
    }

    if (!all(vapply(results_list, is.data.frame, logical(1)))) {
      stop("At least one incremental iteration failed; no incomplete summary exported.")
    }
    results_list <- Filter(function(x) is.data.frame(x), results_list)

    cat("\nBinding results from", length(results_list), "successful iterations...\n")

    bind_rows(results_list)
  }

  psa_inc_results <- run_incremental_psa(n_psa_samples = N_PSA,
    posterior_samples = posterior_samples, fixed_params = fixed_params,
    cost_params_dist = cost_params_dist, lp_params_dist = lp_params_dist,
    lognormal_params_ri_death = lognormal_params_ri_death,
    lognormal_params_si_death = lognormal_params_si_death,
    los_epi_ranges = los_epi_ranges, icu_prop_params = icu_prop_params,
    full_sim_years = full_sim_years, economic_horizon_years = economic_horizon_years,
    discount_rate = discount_rate, base_seed = SEED, n_cores = 1L)
  expected_n <- N_PSA
} else {
  psa_inc_results <- readRDS("data/psa_inc_results.rds")
  expected_n <- 1000L
}
incremental <- as.data.frame(psa_inc_results)
stopifnot(nrow(incremental) == 5L * expected_n,
          length(unique(incremental$iter)) == expected_n,
          !anyDuplicated(incremental[c("iter", "year_idx")]),
          all(table(incremental$iter) == 5L), setequal(unique(incremental$year_idx), 1:5))
stopifnot(all(vapply(incremental[vapply(incremental, is.numeric, logical(1))],
  function(x) all(is.finite(x)), logical(1))))
saveRDS(incremental, file.path(out, "psa_inc_results.rds"), compress = "xz")
attr <- rebuild_attribution(incremental, parameters)
attr5 <- sum_by_iteration(attr, c("direct", "productivity", "total", "deaths"))
write_result(attr5, "attribution_by_iteration", out)
attr_summary <- do.call(rbind, lapply(c("total", "direct", "productivity", "deaths"), function(n) {
  summary_row(attr5[[n]], n, if (n == "deaths") "deaths" else "million SGD", if (n == "deaths") 1 else 1e6)
}))
write_result(attr_summary, "attributable_burden_5yr", out)
attr_annual <- do.call(rbind, lapply(split(attr, attr$calendar_year), function(z) {
  do.call(rbind, lapply(c("total", "direct", "productivity", "deaths"), function(n) {
    cbind(calendar_year = z$calendar_year[1],
      summary_row(z[[n]], n, if (n == "deaths") "deaths" else "million SGD", if (n == "deaths") 1 else 1e6))
  }))
}))
write_result(attr_annual, "attributable_burden_annual", out)

if (RUN_MODE == "reported") stopifnot(round(attr_summary$mean[attr_summary$outcome == "total"], 2) == 317.29)
print(attr_summary, row.names = FALSE)
input_description <- paste(RUN_MODE, "dedicated positive-premium PSA")

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
writeLines(c("Completed", paste("Input:", input_description),
  paste("Time:", format(Sys.time(), tz = "UTC", usetz = TRUE))), file.path(out, "RUN_COMPLETE.txt"))
setwd(previous_directory)
cat("\nResults saved to:", out, "\n")
