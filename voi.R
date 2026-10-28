# Run after total_cost_1028.R, from the 1028_code folder in RStudio.
psa_results <- readRDS("outputs/total_cost_fit7/psa_results.rds")
set.seed(123)
output_dir <- "outputs/voi_fit7"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ============================================================
## BASELINE EVPPI FIGURE — REVISED 4-OUTCOME VERSION
## ESBL-producing Enterobacterales cost outcomes
## Compatible with revised PSA output
## ============================================================
##
## Outcomes:
##   1. Total ESBL cost:
##        cost_hosp_ri + cost_comm_ri + cost_lp_ri_h + cost_lp_ri_c
##
##   2. Hospital direct cost:
##        cost_hosp_ri
##
##   3. Community direct cost:
##        cost_comm_ri
##
##   4. Productivity loss:
##        cost_lp_ri_h + cost_lp_ri_c
##
## Y-axis: % of variance explained
## X-axis: top 5 parameters within each outcome
##
## Notes:
##   - Hospital direct cost is bed-day-only in the revised PSA.
##   - Community direct cost is outpatient visit + antibiotic course only.
##   - Hospital antibiotic, diagnostic, and routine-cost parameters are excluded
##     because they are no longer part of the revised costing model.
## ============================================================

library(voi)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(purrr)

## ── 1. Check required PSA columns ────────────────────────────

required_cols <- c(
  "scenario", "iter",
  "cost_hosp_ri", "cost_comm_ri", "cost_lp_ri_h", "cost_lp_ri_c"
)

missing_cols <- setdiff(required_cols, names(psa_results))
if (length(missing_cols) > 0) {
  stop("Missing required columns in psa_results: ", paste(missing_cols, collapse = ", "))
}

## ── 2. Create ESBL cost outcomes ─────────────────────────────

psa_results_granular <- psa_results %>%
  mutate(
    cost_esbl_total =
      cost_hosp_ri + cost_comm_ri + cost_lp_ri_h + cost_lp_ri_c,

    cost_esbl_hospital_direct =
      cost_hosp_ri,

    cost_esbl_community_direct =
      cost_comm_ri,

    cost_esbl_productivity =
      cost_lp_ri_h + cost_lp_ri_c
  )

## ── 3. Restrict to baseline scenario only ────────────────────

baseline_df <- psa_results_granular %>%
  filter(scenario == "Baseline") %>%
  arrange(iter)

if (nrow(baseline_df) == 0) {
  stop("No rows found for Baseline scenario.")
}

cat("Number of baseline PSA iterations:", nrow(baseline_df), "\n\n")

## ── 4. Define outcomes for the figure ────────────────────────

plot_outcomes <- list(
  list(
    type  = "ESBL_Total_Cost",
    label = "Total cost",
    col   = "cost_esbl_total"
  ),
  list(
    type  = "ESBL_Hospital_Direct_Cost",
    label = "Hospital direct cost",
    col   = "cost_esbl_hospital_direct"
  ),
  list(
    type  = "ESBL_Community_Direct_Cost",
    label = "Community direct cost",
    col   = "cost_esbl_community_direct"
  ),
  list(
    type  = "ESBL_Productivity_Loss",
    label = "Productivity loss",
    col   = "cost_esbl_productivity"
  )
)

## ── 5. Prepare input parameter matrix ────────────────────────

inputs_raw <- baseline_df %>%
  select(starts_with("draw_")) %>%
  as.data.frame()

if (ncol(inputs_raw) == 0) {
  stop("No draw_ parameter columns found in baseline_df.")
}

names(inputs_raw) <- gsub("^draw_", "", names(inputs_raw))

param_var <- apply(inputs_raw, 2, var, na.rm = TRUE)
variable_pars <- names(param_var[!is.na(param_var) & param_var > 1e-20])

## Baseline-only EVPPI: exclude intervention efficacy parameters.
baseline_excluded_pars <- c("hosp_hh", "comm_hh")
actionable_pars <- setdiff(variable_pars, baseline_excluded_pars)

if (length(actionable_pars) == 0) {
  stop("No actionable parameters available after exclusions.")
}

inputs_actionable <- inputs_raw[, actionable_pars, drop = FALSE]

cat("Number of variable parameters found:", length(variable_pars), "\n")
cat("Number of actionable parameters used in EVPPI:", length(actionable_pars), "\n\n")

## ── 6. Parameter labels and categories ───────────────────────

param_labels <- c(
  ## Transmission / epidemiology
  "beta_hh"       = "Hospital human acquisition, ESBL",
  "beta_sh_sih"   = "Hospital human acquisition, non-ESBL",
  "beta_eh"       = "Hospital environmental acquisition, ESBL",
  "beta_hc"       = "Community human acquisition, ESBL",
  "beta_sc_ac"    = "Community animal acquisition, ESBL",
  "beta_sc_ec"    = "Community environmental acquisition, ESBL",
  "delta_rcc_ric" = "Community progression, ESBL",
  
  ## LOS / clinical
  "gamma_rih_rch" = "Recovery rate, hospital ESBL",
  "gamma_sih_sh"  = "Recovery rate, hospital non-ESBL",
  "mu_rih"        = "Mortality rate, hospital ESBL",
  "mu_sih"        = "Mortality rate, hospital non-ESBL",
  "los_rih"       = "LOS, hospital ESBL",
  "los_sih"       = "LOS, hospital non-ESBL",
  
  ## ICU and hospital bed-day cost
  "prop_rih_icu"  = "ICU proportion, ESBL",
  "prop_sih_icu"  = "ICU proportion, non-ESBL",
  "daily_h_ri"    = "Daily hospital bed cost, ESBL",
  "daily_h_si"    = "Daily hospital bed cost, non-ESBL",
  
  ## Community cost
  "s_comm_abx_r"  = "Community antibiotic course cost, ESBL",
  "s_comm_abx_s"  = "Community antibiotic course cost, non-ESBL",
  "s_cost_c"      = "Community cost per infection episode, ESBL",
  "s_cost_c_si"   = "Community cost per infection episode, non-ESBL",
  
  ## Productivity loss
  "unit_lp"        = "Productivity loss per death, ESBL",
  "unit_lp_si"     = "Productivity loss per death, non-ESBL",
  "ly_lost_val"    = "Life-years lost per death, ESBL",
  "ly_lost_val_si" = "Life-years lost per death, non-ESBL"
)

clinical_epi_params <- c(
  "beta_hh", "beta_sh_sih", "beta_eh", "beta_hc",
  "beta_sc_ac", "beta_sc_ec", "delta_rcc_ric",
  "gamma_rih_rch", "gamma_sih_sh", "mu_rih", "mu_sih",
  "los_rih", "los_sih"
)

econ_params <- c(
  "prop_rih_icu", "prop_sih_icu",
  "daily_h_ri", "daily_h_si",
  "s_comm_abx_r", "s_comm_abx_s",
  "s_cost_c", "s_cost_c_si",
  "unit_lp", "unit_lp_si",
  "ly_lost_val", "ly_lost_val_si"
)

## ── 7. Helper function for parameter-level EVPPI ─────────────

run_single_evppi <- function(Y, inputs_actionable, actionable_pars) {
  valid_idx <- complete.cases(Y, inputs_actionable)

  Y_clean <- Y[valid_idx]
  X_clean <- inputs_actionable[valid_idx, , drop = FALSE]

  if (length(Y_clean) < 30) {
    stop("Too few complete PSA iterations available for EVPPI.")
  }

  if (var(Y_clean, na.rm = TRUE) <= 1e-20) {
    warning("Outcome has near-zero variance. EVPPI results may be uninformative.")
  }

  evppi_result <- voi::evppivar(
    outputs = Y_clean,
    inputs  = X_clean,
    pars    = as.list(actionable_pars),
    method  = "gam"
  )

  list(
    evppi_result = evppi_result,
    Y_clean = Y_clean
  )
}

## ── 8. Run EVPPI for the four outcomes ───────────────────────

evppi_results <- list()

for (oc in plot_outcomes) {
  cat("\nRunning EVPPI for:", oc$label, "\n")

  Y <- baseline_df[[oc$col]]
  evpi_var <- var(Y, na.rm = TRUE)

  if (is.na(evpi_var) || evpi_var <= 1e-20) {
    warning(paste0("Skipping ", oc$label, ": near-zero outcome variance."))
    next
  }

  evppi_run <- run_single_evppi(
    Y = Y,
    inputs_actionable = inputs_actionable,
    actionable_pars = actionable_pars
  )

  evppi_result <- evppi_run$evppi_result
  Y_clean <- evppi_run$Y_clean

  evpi_var_clean <- var(Y_clean, na.rm = TRUE)
  evpi_sd_clean  <- sd(Y_clean, na.rm = TRUE)
  mean_y_clean   <- mean(Y_clean, na.rm = TRUE)

  evppi_results[[oc$type]] <- data.frame(
    outcome_type = oc$type,
    outcome_label = oc$label,
    param = actionable_pars,
    evppi = as.numeric(evppi_result$evppi),
    evpi_var = evpi_var_clean,
    evpi_sd = evpi_sd_clean,
    mean_projected_cost = mean_y_clean,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      evppi = pmax(0, evppi),
      pct_variance = ifelse(evpi_var > 0, evppi / evpi_var * 100, NA_real_),
      Category = case_when(
        param %in% econ_params ~ "Economic / Costing",
        param %in% clinical_epi_params ~ "Clinical / Epidemiological",
        TRUE ~ "Other / Derived"
      ),
      param_label = ifelse(
        param %in% names(param_labels),
        param_labels[param],
        param
      )
    ) %>%
    arrange(desc(pct_variance))
}

evppi_results <- bind_rows(evppi_results)

if (nrow(evppi_results) == 0) {
  stop("No EVPPI results were generated.")
}

## ── 9. Select top 5 parameters within each outcome ───────────

evppi_top5 <- evppi_results %>%
  group_by(outcome_label) %>%
  slice_max(order_by = pct_variance, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    ## Raw EVPPI is variance reduction; sqrt(EVPPI) gives SD-equivalent cost scale
    evppi_sd_equiv_m_sgd = sqrt(evppi) / 1e6,
    evppi_label = paste0("$", round(evppi_sd_equiv_m_sgd, 1), "M")
  )

## Facet-specific x-axis ordering.
evppi_top5 <- evppi_top5 %>%
  group_by(outcome_label) %>%
  arrange(pct_variance, .by_group = TRUE) %>%
  mutate(
    param_label_facet = factor(
      paste(param_label, outcome_label, sep = "___"),
      levels = paste(param_label, outcome_label, sep = "___")
    ),
    outcome_label = factor(
      outcome_label,
      levels = c(
        "Total cost",
        "Hospital direct cost",
        "Community direct cost",
        "Productivity loss"
      )
    )
  ) %>%
  ungroup()

## ── 10. Make combined figure ─────────────────────────────────

p_evppi_top5 <- ggplot(
  evppi_top5,
  aes(x = param_label_facet, y = pct_variance, fill = Category)
) +
  geom_col(color = "black", linewidth = 0.3, width = 0.75) +
  geom_text(
    aes(label = evppi_label),
    vjust = -0.35,
    size = 3
  ) +
  facet_wrap(
    ~ outcome_label,
    scales = "free_x",
    nrow = 2
  ) +
  scale_x_discrete(
    labels = function(x) gsub("___.*$", "", x)
  ) +
  scale_y_continuous(
    labels = function(x) paste0(round(x, 1), "%"),
    expand = expansion(mult = c(0, 0.18))
  ) +
  scale_fill_manual(
    values = c(
      "Clinical / Epidemiological" = "#de2d26",
      "Economic / Costing" = "#3182bd",
      "Other / Derived" = "gray60"
    )
  ) +
  labs(
    title = "Cost Outcomes for ESBL-producing Enterobacterales Infections",
    subtitle = "Top 5 parameters by percentage of variance explained; bar labels show SD-equivalent EVPPI in million SGD",
    x = NULL,
    y = "% of variance explained",
    fill = NULL,
    caption = "Parameter-level estimation EVPPI for projected baseline cost outcomes."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", margin = margin(b = 5, t = 5)),
    strip.background = element_rect(fill = "gray95", color = NA),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    plot.margin = margin(12, 12, 12, 36),
    panel.grid.major.x = element_blank()
  )

if (interactive()) print(p_evppi_top5)


###EVPI
evpi_summary <- evppi_results %>%
  distinct(
    outcome_label,
    mean_projected_cost,
    evpi_var,
    evpi_sd
  ) %>%
  mutate(
    mean_projected_cost_M_SGD = mean_projected_cost / 1e6,
    evpi_sd_M_SGD = evpi_sd / 1e6
  ) %>%
  select(
    outcome_label,
    mean_projected_cost_M_SGD,
    evpi_var,
    evpi_sd_M_SGD
  )

cat("\nEVPI summary by outcome:\n")
print(as.data.frame(evpi_summary))
## ── 11. Save figure and table ────────────────────────────────

ggsave(
  filename = file.path(output_dir, "baseline_esbl_evppi_top5_4_cost_outcomes.png"),
  plot = p_evppi_top5,
  width = 13,
  height = 8,
  bg = "white",
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "baseline_esbl_evppi_top5_4_cost_outcomes.pdf"),
  plot = p_evppi_top5,
  width = 13,
  height = 8,
  bg = "white"
)

write.csv(
  evppi_results,
  file.path(output_dir, "baseline_esbl_evppi_4_cost_outcomes_all_results.csv"),
  row.names = FALSE
)

write.csv(
  evppi_top5,
  file.path(output_dir, "baseline_esbl_evppi_4_cost_outcomes_top5.csv"),
  row.names = FALSE
)

cat("\nSaved outputs:\n")
cat("  - baseline_esbl_evppi_top5_4_cost_outcomes.png\n")
cat("  - baseline_esbl_evppi_top5_4_cost_outcomes.pdf\n")
cat("  - baseline_esbl_evppi_4_cost_outcomes_all_results.csv\n")
cat("  - baseline_esbl_evppi_4_cost_outcomes_top5.csv\n")

write.csv(evpi_summary, file.path(output_dir, "evpi_summary.csv"), row.names = FALSE)


save(evpi_summary, evppi_results, evppi_top5, file = file.path(output_dir, "voi_fit7.RData"))
