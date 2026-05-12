################################################################################
#
#  TB VACCINATION MODEL — STRATEGY EXTRACTION (PLAUSIBLE + OPTIMISTIC)
#
################################################################################
#
# PURPOSE:
#   Runs each vaccination strategy under both scenarios and extracts:
#     1. Equilibrium compartment sizes by stratum
#     2. Strategy-level summary (cases prevented, NNV, % reduction)
#     3. Case sources (fast / slow / relapse breakdown)
#     4. LTBI pool changes by scenario
#     5. Annual transition flows
#
#   These outputs share the same set of ODE runs, so bundling them avoids
#   running each strategy 5 separate times.
#
#   Plausible scenario:  strategy-specific psi, VE = 50%
#   Optimistic scenario: psi = 50%/yr for all strategies, VE = 70%
#
# REQUIRES:
#   - deSolve package
#   - TB_VACCINATION_MODEL_COMPLETE.R in working directory
#
# USAGE:
#   1. Set working directory to folder containing the model script
#   2. Source this file: source("TB_EXTRACT_STRATEGIES.R")
#   3. Results saved as CSVs in working directory
#
#   To run only one scenario, comment out the unwanted entry in the
#   `scenarios` list below.
#
# OUTPUT FILES:
#   Plausible:
#     - TB_COMPARTMENT_SIZES.csv
#     - TB_STRATEGY_SUMMARY.csv
#     - TB_CASE_SOURCES.csv
#     - TB_LTBI_POOL_ANALYSIS.csv
#     - TB_ANNUAL_FLOWS.csv
#   Optimistic:
#     - TB_COMPARTMENT_SIZES_OPTIMISTIC.csv
#     - TB_STRATEGY_SUMMARY_OPTIMISTIC.csv
#     - TB_CASE_SOURCES_OPTIMISTIC.csv
#     - TB_LTBI_POOL_ANALYSIS_OPTIMISTIC.csv
#     - TB_ANNUAL_FLOWS_OPTIMISTIC.csv
#
# NOTE:
#   Stratum-level intervention impact is NOT extracted here (it has its own
#   dedicated script TB_STRATUM_IMPACT.R for both scenarios).
#
# EXPECTED RUNTIME: ~2-4 minutes total
#
# DATE: May 2026
#
################################################################################

suppressPackageStartupMessages({
  library(deSolve)
})

cat("================================================================\n")
cat("  TB MODEL — STRATEGY EXTRACTION (BOTH SCENARIOS)\n")
cat("================================================================\n\n")

# --- Source the base model (runs calibration automatically) ---
cat("Loading base model and running calibration...\n")
invisible(capture.output(suppressMessages(
  source("TB_VACCINATION_MODEL_COMPLETE.R")
)))
cat("Base model loaded.\n")
cat(sprintf("  Baseline beta = %.4f\n\n", calibration$beta))

# ---- Strategies ----------------------------------------------------------
# Plausible: heterogeneous psi, VE = 0.50
# Optimistic: uniform psi = 0.50, VE = 0.70 (theta vectors same)
strategies_plaus <- list(
  "Baseline"          = list(theta = c(0, 0, 0, 0),       psi = 0,    VE = 0.50),
  "All High-Risk"     = list(theta = c(1, 1, 1, 0),       psi = 0.05, VE = 0.50),
  "All Foreign-Born"  = list(theta = c(0.19, 0.18, 1, 0), psi = 0.05, VE = 0.50),
  "FB Stratum Only"   = list(theta = c(0, 0, 1, 0),       psi = 0.05, VE = 0.50),
  "Universal"         = list(theta = c(1, 1, 1, 1),       psi = 0.02, VE = 0.50),
  "Medical Only"      = list(theta = c(0, 1, 0, 0),       psi = 0.05, VE = 0.50),
  "HIV Only"          = list(theta = c(1, 0, 0, 0),       psi = 0.10, VE = 0.50),
  "HIV + Medical"     = list(theta = c(1, 1, 0, 0),       psi = 0.05, VE = 0.50)
)

OPT_PSI <- 0.50
OPT_VE  <- 0.70
strategies_opt <- list(
  "Baseline"          = list(theta = c(0, 0, 0, 0),       psi = 0,       VE = OPT_VE),
  "All High-Risk"     = list(theta = c(1, 1, 1, 0),       psi = OPT_PSI, VE = OPT_VE),
  "All Foreign-Born"  = list(theta = c(0.19, 0.18, 1, 0), psi = OPT_PSI, VE = OPT_VE),
  "FB Stratum Only"   = list(theta = c(0, 0, 1, 0),       psi = OPT_PSI, VE = OPT_VE),
  "Universal"         = list(theta = c(1, 1, 1, 1),       psi = OPT_PSI, VE = OPT_VE),
  "Medical Only"      = list(theta = c(0, 1, 0, 0),       psi = OPT_PSI, VE = OPT_VE),
  "HIV Only"          = list(theta = c(1, 0, 0, 0),       psi = OPT_PSI, VE = OPT_VE),
  "HIV + Medical"     = list(theta = c(1, 1, 0, 0),       psi = OPT_PSI, VE = OPT_VE)
)

# ---- Define scenarios ----------------------------------------------------
scenarios <- list(
  plausible = list(
    label       = "PLAUSIBLE",
    strategies  = strategies_plaus,
    csv_comp    = "TB_COMPARTMENT_SIZES.csv",
    csv_summary = "TB_STRATEGY_SUMMARY.csv",
    csv_sources = "TB_CASE_SOURCES.csv",
    csv_ltbi    = "TB_LTBI_POOL_ANALYSIS.csv",
    csv_flows   = "TB_ANNUAL_FLOWS.csv"
  ),
  optimistic = list(
    label       = "OPTIMISTIC",
    strategies  = strategies_opt,
    csv_comp    = "TB_COMPARTMENT_SIZES_OPTIMISTIC.csv",
    csv_summary = "TB_STRATEGY_SUMMARY_OPTIMISTIC.csv",
    csv_sources = "TB_CASE_SOURCES_OPTIMISTIC.csv",
    csv_ltbi    = "TB_LTBI_POOL_ANALYSIS_OPTIMISTIC.csv",
    csv_flows   = "TB_ANNUAL_FLOWS_OPTIMISTIC.csv"
  )
)

# ---- Helper functions ----------------------------------------------------
format_compartments <- function(state, scenario_name) {
  comp <- extract_compartments(state)
  data.frame(
    Scenario = scenario_name,
    Stratum  = stratum_names,
    N_pop    = N_vec,
    S = comp$S, L_f = comp$L_f, L_s = comp$L_s,
    L_fv = comp$L_fv, L_sv = comp$L_sv,
    I = comp$I, I_v = comp$I_v, R = comp$R,
    Total_LTBI       = comp$L_f + comp$L_s + comp$L_fv + comp$L_sv,
    Total_LTBI_unvax = comp$L_f + comp$L_s,
    Total_LTBI_vax   = comp$L_fv + comp$L_sv,
    Total_Infectious = comp$I + comp$I_v,
    Total_N          = comp$S + comp$L_f + comp$L_s + comp$L_fv +
                       comp$L_sv + comp$I + comp$I_v + comp$R,
    LTBI_prevalence_pct = (comp$L_f + comp$L_s + comp$L_fv + comp$L_sv) /
                          (comp$S + comp$L_f + comp$L_s + comp$L_fv +
                           comp$L_sv + comp$I + comp$I_v + comp$R) * 100,
    Lf_fraction_of_LTBI = comp$L_f /
                          (comp$L_f + comp$L_s + comp$L_fv + comp$L_sv + 0.001) * 100
  )
}

calc_annual_flows <- function(state, params) {
  comp <- extract_compartments(state)
  VE <- params$VE
  psi <- params$psi
  theta <- params$theta

  flows <- data.frame(
    Stratum         = stratum_names,
    Flow_Lf_to_I    = eps_f * comp$L_f,
    Flow_Ls_to_I    = eps_s * comp$L_s,
    Flow_Lfv_to_I   = (1 - VE) * eps_f * comp$L_fv,
    Flow_Lsv_to_I   = (1 - VE) * eps_s * comp$L_sv,
    Flow_R_to_I     = rho * comp$R,
    Flow_Lf_to_Lfv  = psi * theta * comp$L_f,
    Flow_Ls_to_Lsv  = psi * theta * comp$L_s,
    Flow_Lfv_to_Lf  = omega * comp$L_fv,
    Flow_Lsv_to_Ls  = omega * comp$L_sv,
    Flow_Lf_to_Ls   = nu_f_trans * comp$L_f,
    Flow_Lfv_to_Lsv = nu_f_trans * comp$L_fv,
    Flow_Lf_clear   = nu_f_clear * comp$L_f,
    Flow_Ls_clear   = nu_s_clear * comp$L_s,
    Flow_I_to_R     = gamma * comp$I,
    Flow_Iv_to_R    = gamma * comp$I_v,
    Flow_R_to_Ls    = alpha_stab * comp$R
  )

  flows$Total_new_cases     <- flows$Flow_Lf_to_I + flows$Flow_Ls_to_I +
                               flows$Flow_Lfv_to_I + flows$Flow_Lsv_to_I +
                               flows$Flow_R_to_I
  flows$Total_vaccinations  <- flows$Flow_Lf_to_Lfv + flows$Flow_Ls_to_Lsv

  flows
}

################################################################################
# RUN EACH SCENARIO
################################################################################

for (scen_name in names(scenarios)) {

  scen <- scenarios[[scen_name]]

  cat("\n\n")
  cat("################################################################\n")
  cat(sprintf("#  SCENARIO: %s\n", scen$label))
  cat("################################################################\n\n")

  # ---- Initialize collectors ---------------------------------------------
  all_compartments     <- data.frame()
  all_strategy_summary <- data.frame()
  all_case_sources     <- data.frame()
  all_flows            <- data.frame()

  # ---- Get baseline first (we need baseline_total_cases for prevented) ----
  baseline_strat <- scen$strategies[["Baseline"]]
  params_base <- list(
    beta = calibration$beta, kappa = c(1, 1, 1, 1),
    psi  = baseline_strat$psi, theta = baseline_strat$theta,
    mtb_prev = calibration$mtb_prev,
    M = M_default, N = N_vec,
    VE = baseline_strat$VE
  )
  baseline_result    <- run_model(params_base, years = 50)
  baseline_state     <- baseline_result$final_state
  baseline_cases_obj <- calc_annual_cases(baseline_state, params_base)
  baseline_total_cases <- baseline_cases_obj$total

  # ---- Process each strategy ---------------------------------------------
  for (strat_name in names(scen$strategies)) {
    strat <- scen$strategies[[strat_name]]
    cat(sprintf("  %-20s  (psi = %.2f, VE = %.2f) ... ",
                strat_name, strat$psi, strat$VE))

    params <- list(
      beta = calibration$beta, kappa = c(1, 1, 1, 1),
      psi  = strat$psi, theta = strat$theta,
      mtb_prev = calibration$mtb_prev,
      M = M_default, N = N_vec,
      VE = strat$VE
    )

    result      <- run_model(params, years = 50)
    final_state <- result$final_state

    # --- Compartment sizes ---
    comp_df <- format_compartments(final_state, strat_name)
    comp_df$psi          <- strat$psi
    comp_df$theta_HIV    <- strat$theta[1]
    comp_df$theta_Medical <- strat$theta[2]
    comp_df$theta_FB     <- strat$theta[3]
    comp_df$theta_USB    <- strat$theta[4]
    comp_df$VE           <- strat$VE
    all_compartments <- rbind(all_compartments, comp_df)

    # --- Annual cases and sources ---
    cases_obj   <- calc_annual_cases(final_state, params)
    total_cases <- cases_obj$total

    case_source_df <- data.frame(
      Scenario         = strat_name,
      Total_cases      = total_cases,
      From_Lf          = cases_obj$from_Lf,
      From_Ls          = cases_obj$from_Ls,
      From_relapse     = cases_obj$from_relapse,
      Pct_from_Lf      = cases_obj$pct_Lf,
      Pct_from_Ls      = cases_obj$pct_Ls,
      Pct_from_relapse = cases_obj$pct_relapse,
      Pct_reactivation = cases_obj$pct_reactivation
    )
    all_case_sources <- rbind(all_case_sources, case_source_df)

    # --- Strategy summary ---
    vax_obj         <- calc_vaccinations(final_state, strat$psi, strat$theta)
    annual_vax      <- vax_obj$total
    cases_prevented <- baseline_total_cases - total_cases
    pct_reduction   <- cases_prevented / baseline_total_cases * 100
    nnv             <- ifelse(cases_prevented > 0, annual_vax / cases_prevented, NA)

    summary_df <- data.frame(
      Scenario           = strat_name,
      psi                = strat$psi,
      VE                 = strat$VE,
      Baseline_cases     = baseline_total_cases,
      Annual_cases       = total_cases,
      Cases_prevented    = cases_prevented,
      Pct_reduction      = pct_reduction,
      Annual_vaccinations = annual_vax,
      NNV                = nnv
    )
    all_strategy_summary <- rbind(all_strategy_summary, summary_df)

    # --- Annual flows ---
    flows_df <- calc_annual_flows(final_state, params)
    flows_df$Scenario <- strat_name
    all_flows <- rbind(all_flows, flows_df)

    cat("done\n")
  }

  # ---- LTBI pool analysis (aggregated across strata per scenario) --------
  ltbi_summary <- aggregate(
    cbind(L_f, L_s, L_fv, L_sv, Total_LTBI, Total_LTBI_unvax, Total_LTBI_vax,
          Total_Infectious, S, R) ~ Scenario,
    data = all_compartments,
    FUN  = sum
  )

  baseline_ltbi <- ltbi_summary$Total_LTBI[ltbi_summary$Scenario == "Baseline"]
  baseline_Lf   <- ltbi_summary$L_f[ltbi_summary$Scenario == "Baseline"]
  baseline_Ls   <- ltbi_summary$L_s[ltbi_summary$Scenario == "Baseline"]
  baseline_inf  <- ltbi_summary$Total_Infectious[ltbi_summary$Scenario == "Baseline"]

  ltbi_summary$LTBI_pct_of_baseline       <- ltbi_summary$Total_LTBI / baseline_ltbi * 100
  ltbi_summary$Lf_pct_of_baseline         <- ltbi_summary$L_f       / baseline_Lf   * 100
  ltbi_summary$Ls_pct_of_baseline         <- ltbi_summary$L_s       / baseline_Ls   * 100
  ltbi_summary$Pct_LTBI_vaccinated        <- ltbi_summary$Total_LTBI_vax /
                                              ltbi_summary$Total_LTBI * 100
  ltbi_summary$Lf_fraction                <- ltbi_summary$L_f / ltbi_summary$Total_LTBI * 100
  ltbi_summary$Infectious_pct_of_baseline <- ltbi_summary$Total_Infectious /
                                              baseline_inf * 100

  # ---- Write CSVs --------------------------------------------------------
  write.csv(all_compartments,     scen$csv_comp,    row.names = FALSE)
  write.csv(all_strategy_summary, scen$csv_summary, row.names = FALSE)
  write.csv(all_case_sources,     scen$csv_sources, row.names = FALSE)
  write.csv(ltbi_summary,         scen$csv_ltbi,    row.names = FALSE)
  write.csv(all_flows,            scen$csv_flows,   row.names = FALSE)

  cat(sprintf("\n  Outputs (%s):\n", scen$label))
  cat(sprintf("    %s\n", scen$csv_comp))
  cat(sprintf("    %s\n", scen$csv_summary))
  cat(sprintf("    %s\n", scen$csv_sources))
  cat(sprintf("    %s\n", scen$csv_ltbi))
  cat(sprintf("    %s\n", scen$csv_flows))
}

cat("\n\n================================================================\n")
cat("  ALL SCENARIOS COMPLETE\n")
cat("================================================================\n")
