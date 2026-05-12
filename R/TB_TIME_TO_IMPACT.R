################################################################################
#
#  TB VACCINATION MODEL — TIME-TO-IMPACT (PLAUSIBLE + OPTIMISTIC)
#
################################################################################
#
# PURPOSE:
#   Time-to-impact / transient dynamics under both scenarios. Extracts annual
#   cases prevented at years 1, 2, 3, 5, 10, 15, 20, 25, 30 and computes
#   cumulative cases prevented and NNV at each time horizon.
#
#   Plausible scenario:  VE = 50%, strategy-specific psi (6 strategies)
#   Optimistic scenario: VE = 70%, psi = 50%/yr for all strategies
#                        (7 strategies, includes FB Stratum Only)
#
#   Each scenario uses its own historical strategy set so existing CSV
#   outputs are reproduced bit-for-bit.
#
# REQUIRES:
#   - deSolve package
#   - TB_VACCINATION_MODEL_COMPLETE.R in working directory
#
# USAGE:
#   1. Set working directory to folder containing the model script
#   2. Source this file: source("TB_TIME_TO_IMPACT.R")
#   3. Results saved as CSVs in working directory
#
#   To run only one scenario, comment out the unwanted entry in the
#   `scenarios` list below.
#
# OUTPUT FILES:
#   Plausible:
#     - TB_TIME_TO_IMPACT_RESULTS.csv
#     - TB_TIME_TO_IMPACT_RESULTS.RData
#   Optimistic:
#     - TB_TIME_TO_IMPACT_RESULTS_OPTIMISTIC.csv
#     - TB_TIME_TO_IMPACT_OPTIMISTIC_RESULTS.RData
#
# EXPECTED RUNTIME: ~3-5 minutes total
#
# DATE: May 2026
#
################################################################################

suppressPackageStartupMessages({
  library(deSolve)
})

cat("================================================================\n")
cat("  TB VACCINATION MODEL — TIME-TO-IMPACT (BOTH SCENARIOS)\n")
cat("================================================================\n\n")

# --- Source the base model (runs calibration automatically) ---
cat("Loading base model and running calibration...\n")
invisible(capture.output(suppressMessages(
  source("TB_VACCINATION_MODEL_COMPLETE.R")
)))
cat("Base model loaded.\n")
cat(sprintf("  Baseline beta = %.4f\n", calibration$beta))
cat(sprintf("  Actual 2024 cases = %d\n\n", actual_2024_cases))

baseline_cal <- calibration

# ---- Common simulation parameters ----------------------------------------
sim_years   <- 30
time_points <- c(1, 2, 3, 5, 10, 15, 20, 25, 30)

# ---- Plausible strategies (6, heterogeneous psi) -------------------------
strategies_plaus <- list(
  all_high_risk = list(name = "All High-Risk",   psi = 0.05, theta = c(1, 1, 1, 0)),
  plwh_all_nusb = list(name = "PLWH + All NUSB", psi = 0.05, theta = c(1, 0.18, 1, 0)),
  universal     = list(name = "Universal",       psi = 0.02, theta = c(1, 1, 1, 1)),
  medical       = list(name = "Medical Only",    psi = 0.05, theta = c(0, 1, 0, 0)),
  hiv           = list(name = "HIV Only",        psi = 0.10, theta = c(1, 0, 0, 0)),
  hiv_medical   = list(name = "HIV + Medical",   psi = 0.05, theta = c(1, 1, 0, 0))
)

# ---- Optimistic strategies (7, all psi = 0.50) ---------------------------
OPT_PSI <- 0.50
strategies_opt <- list(
  all_high_risk = list(name = "All High-Risk (HIV+Med+FB)", psi = OPT_PSI, theta = c(1, 1, 1, 0)),
  plwh_nusb     = list(name = "PLWH + Non-U.S.-Born",       psi = OPT_PSI, theta = c(1, 0.18, 1, 0)),
  fb_stratum    = list(name = "FB Stratum Only",            psi = OPT_PSI, theta = c(0, 0, 1, 0)),
  universal     = list(name = "Universal",                  psi = OPT_PSI, theta = c(1, 1, 1, 1)),
  medical       = list(name = "Medical Only",               psi = OPT_PSI, theta = c(0, 1, 0, 0)),
  hiv           = list(name = "HIV Only",                   psi = OPT_PSI, theta = c(1, 0, 0, 0)),
  hiv_medical   = list(name = "HIV + Medical",              psi = OPT_PSI, theta = c(1, 1, 0, 0))
)

# ---- Define scenarios ----------------------------------------------------
scenarios <- list(
  plausible = list(
    label      = "PLAUSIBLE",
    strategies = strategies_plaus,
    VE         = 0.50,
    csv        = "TB_TIME_TO_IMPACT_RESULTS.csv",
    rdata      = "TB_TIME_TO_IMPACT_RESULTS.RData"
  ),
  optimistic = list(
    label      = "OPTIMISTIC",
    strategies = strategies_opt,
    VE         = 0.70,
    csv        = "TB_TIME_TO_IMPACT_RESULTS_OPTIMISTIC.csv",
    rdata      = "TB_TIME_TO_IMPACT_OPTIMISTIC_RESULTS.RData"
  )
)

# ---- Helper: extract outcomes at specific time points from trajectory ----
extract_at_times <- function(trajectory, params, t_pts) {
  out <- data.frame(year = numeric(), annual_cases = numeric(),
                    cases_HIV = numeric(), cases_Med = numeric(),
                    cases_FB = numeric(), cases_USB = numeric(),
                    annual_vax = numeric())
  for (t in t_pts) {
    row_idx <- which.min(abs(trajectory[, 1] - t))
    state   <- as.numeric(trajectory[row_idx, -1])
    cases   <- calc_annual_cases(state, params)
    vax     <- calc_vaccinations(state, params$psi, params$theta)
    out <- rbind(out, data.frame(
      year         = t,
      annual_cases = cases$total,
      cases_HIV    = cases$by_stratum[1],
      cases_Med    = cases$by_stratum[2],
      cases_FB     = cases$by_stratum[3],
      cases_USB    = cases$by_stratum[4],
      annual_vax   = vax$total
    ))
  }
  out
}

# ---- Helper: compute cumulative impact via trapezoidal integration ------
add_cumulative_impact <- function(impact, prev_vals, vax_vals, t_pts) {
  for (k in seq_len(nrow(impact))) {
    idx <- 1:k
    if (length(idx) >= 2) {
      yrs <- t_pts[idx]
      cp  <- cv <- 0
      for (j in 2:length(idx)) {
        dt <- yrs[j] - yrs[j - 1]
        cp <- cp + dt * (prev_vals[idx[j - 1]] + prev_vals[idx[j]]) / 2
        cv <- cv + dt * (vax_vals[idx[j - 1]]  + vax_vals[idx[j]])  / 2
      }
      impact$cumul_prevented[k] <- round(cp)
      impact$cumul_vax[k]       <- round(cv)
      impact$cumul_NNV[k]       <- if (cp > 0) round(cv / cp) else NA
    } else {
      # Year 1: rectangle approximation
      impact$cumul_prevented[k] <- round(prev_vals[1])
      impact$cumul_vax[k]       <- round(vax_vals[1])
      impact$cumul_NNV[k]       <- if (prev_vals[1] > 0)
                                      round(vax_vals[1] / prev_vals[1])
                                   else NA
    }
  }
  impact
}

################################################################################
# RUN EACH SCENARIO
################################################################################

for (scen_name in names(scenarios)) {

  scen <- scenarios[[scen_name]]

  cat("\n\n")
  cat("################################################################\n")
  cat(sprintf("#  SCENARIO: %s  (VE = %.0f%%)\n", scen$label, scen$VE * 100))
  cat("################################################################\n\n")

  # ---- Run baseline (no vaccination) for this scenario ---
  params_base <- list(
    beta = baseline_cal$beta, kappa = c(1, 1, 1, 1),
    psi = 0, theta = c(0, 0, 0, 0),
    mtb_prev = baseline_cal$mtb_prev,
    M = M_default, N = N_vec, VE = scen$VE
  )
  base_run <- run_model(params_base, years = sim_years)
  base_ts  <- extract_at_times(base_run$trajectory, params_base, time_points)

  # Scale factor: model equilibrium -> 2024 surveillance
  base_eq_cases <- base_ts$annual_cases[nrow(base_ts)]
  sf <- actual_2024_cases / base_eq_cases

  # ---- Run each strategy and compute time-varying impact ---
  time_impact_all <- list()

  for (s in names(scen$strategies)) {
    strat <- scen$strategies[[s]]
    cat(sprintf("  %-30s  (psi = %g/yr) ... ", strat$name, strat$psi))

    params_int <- list(
      beta = baseline_cal$beta, kappa = c(1, 1, 1, 1),
      psi = strat$psi, theta = strat$theta,
      mtb_prev = baseline_cal$mtb_prev,
      M = M_default, N = N_vec, VE = scen$VE
    )

    int_run <- run_model(params_int, years = sim_years)
    int_ts  <- extract_at_times(int_run$trajectory, params_int, time_points)

    # Compute impact at each time point
    impact <- data.frame(
      year     = time_points,
      strategy = strat$name,

      # Instantaneous annual impact
      baseline_cases     = round(base_ts$annual_cases * sf),
      intervention_cases = round(int_ts$annual_cases * sf),
      cases_prevented    = round((base_ts$annual_cases - int_ts$annual_cases) * sf),
      pct_reduction      = round((base_ts$annual_cases - int_ts$annual_cases) /
                                  base_ts$annual_cases * 100, 2),
      annual_vax         = round(int_ts$annual_vax),

      # Placeholders for cumulative values (filled below)
      cumul_prevented = NA_real_,
      cumul_vax       = NA_real_,
      cumul_NNV       = NA_real_,

      # Stratum-level prevented
      prev_HIV = round((base_ts$cases_HIV - int_ts$cases_HIV) * sf),
      prev_Med = round((base_ts$cases_Med - int_ts$cases_Med) * sf),
      prev_FB  = round((base_ts$cases_FB  - int_ts$cases_FB)  * sf),
      prev_USB = round((base_ts$cases_USB - int_ts$cases_USB) * sf)
    )

    # Cumulative via trapezoidal integration
    prev_vals <- (base_ts$annual_cases - int_ts$annual_cases) * sf
    vax_vals  <- int_ts$annual_vax
    impact    <- add_cumulative_impact(impact, prev_vals, vax_vals, time_points)

    time_impact_all[[s]] <- impact
    cat("done\n")
  }

  # Combine into one data frame
  time_impact_combined <- do.call(rbind, time_impact_all)
  rownames(time_impact_combined) <- NULL

  # ===========================================================================
  # Cross-strategy comparison summary
  # ===========================================================================
  cat(sprintf("\n  --- Cross-strategy comparison at key time horizons (%s) ---\n",
              scen$label))
  for (t_horizon in c(5, 10, 20, 30)) {
    cat(sprintf("\n  At Year %d:\n", t_horizon))
    cat(sprintf("    %-30s  %10s  %8s  %12s  %10s\n",
                "Strategy", "Prev/yr", "% Red", "Cumul Prev", "Cumul NNV"))
    cat("    ", paste(rep("-", 78), collapse = ""), "\n", sep = "")
    for (s in names(scen$strategies)) {
      df <- time_impact_all[[s]]
      row <- df[df$year == t_horizon, ]
      if (nrow(row) == 1) {
        cat(sprintf("    %-30s  %10d  %7.1f%%  %12s  %10s\n",
                    scen$strategies[[s]]$name,
                    row$cases_prevented,
                    row$pct_reduction,
                    format(row$cumul_prevented, big.mark = ","),
                    ifelse(is.na(row$cumul_NNV), "--",
                           format(row$cumul_NNV, big.mark = ","))))
      }
    }
  }

  # Write CSV and RData
  write.csv(time_impact_combined, scen$csv, row.names = FALSE)
  save(time_impact_all, time_impact_combined,
       baseline_cal,
       file = scen$rdata)

  cat(sprintf("\n  Outputs (%s):\n", scen$label))
  cat(sprintf("    %s\n", scen$csv))
  cat(sprintf("    %s\n", scen$rdata))
}

cat("\n\n================================================================\n")
cat("  ALL SCENARIOS COMPLETE\n")
cat("================================================================\n")
