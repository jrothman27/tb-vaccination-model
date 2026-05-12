################################################################################
#
#  TB VACCINATION MODEL — LHS UNCERTAINTY (PLAUSIBLE + OPTIMISTIC)
#
################################################################################
#
# PURPOSE:
#   Latin Hypercube Sampling (LHS) uncertainty intervals for both the plausible
#   and optimistic scenarios in a single run. Each scenario:
#     - Varies 12 parameters simultaneously across plausible ranges
#     - Recalibrates (beta, mtb_prev) for each draw to maintain incidence fit
#     - Runs all 6 strategies per draw -> 95% uncertainty intervals
#     - Reuses each draw's calibration across all 6 strategies
#
#   Scenarios differ only in the strategy-level psi (vaccination rate). The
#   LHS sweep itself uses the same 12-parameter ranges in both scenarios so
#   uncertainty bands are directly comparable on every parameter except psi.
#
#   Plausible scenario:   strategy-specific psi (0.02 to 0.10/yr)
#   Optimistic scenario:  psi = 0.50/yr for all strategies
#
# REQUIRES:
#   - deSolve, lhs packages
#   - TB_VACCINATION_MODEL_COMPLETE.R in working directory
#
# USAGE:
#   1. Set working directory to folder containing the model script
#   2. Source this file: source("TB_LHS.R")
#   3. Results saved as CSVs in working directory
#
#   To run only one scenario, comment out the unwanted entry in the
#   `scenarios` list below.
#
# OUTPUT FILES:
#   Plausible scenario:
#     - TB_LHS_UNCERTAINTY_INTERVALS.csv
#     - TB_PRCC_RESULTS.csv
#     - TB_LHS_RESULTS.RData
#   Optimistic scenario:
#     - TB_LHS_UNCERTAINTY_INTERVALS_OPTIMISTIC.csv
#     - TB_LHS_OPTIMISTIC_RESULTS.RData
#
# NOTE ON STRATEGY NAMES:
#   This script writes strategy names in their manuscript labels (e.g.,
#   "All Mtb-Infected", "PLWH", "Medical comorbidities") rather than raw
#   R-code labels (e.g., "Universal", "HIV Only", "Medical Only"). Any
#   downstream script that joins on the old raw labels needs to be updated
#   accordingly. The figures script's clean_strategy() already handles both.
#
# EXPECTED RUNTIME: ~60-120 min for 500 draws x 2 scenarios
#
# DATE: May 2026
#
################################################################################

library(deSolve)
library(lhs)

cat("================================================================\n")
cat("  TB VACCINATION MODEL — LHS UNCERTAINTY (BOTH SCENARIOS)\n")
cat("================================================================\n\n")

# --- Source the base model (runs calibration automatically) ---
cat("Loading base model and running calibration...\n")
invisible(capture.output(suppressMessages(
  source("TB_VACCINATION_MODEL_COMPLETE.R")
)))
cat("Base model loaded.\n")
cat(sprintf("  Baseline beta = %.4f\n", calibration$beta))
cat(sprintf("  Baseline mtb_prev = [%.4f, %.4f, %.4f, %.4f]\n",
            calibration$mtb_prev[1], calibration$mtb_prev[2],
            calibration$mtb_prev[3], calibration$mtb_prev[4]))
cat(sprintf("  Actual 2024 cases = %d\n\n", actual_2024_cases))

# Store baseline for reference
baseline_cal <- calibration


# ===========================================================================
# 1A. DEFINE PARAMETER DISTRIBUTIONS (shared across scenarios)
# ===========================================================================
#
# 12 uncertain parameters varied simultaneously, with ranges drawn from
# Ekramnia 2024 Table 3/4 95% UIs where available, or +/-50-60% of point
# estimates for parameters without published CIs.
#
# NOT varied (fixed or minimal impact from one-way analysis):
#   - mu, mu_TB, rho: demographic/mortality (<=1% impact in one-way)
#   - p: fraction to fast latent (absorbed into eps_f structure)
#   - N_vec: population sizes (fixed from Census)
#   - iota: immigration (<=3% impact, correlated with N_FB)
#   - nu_f_trans, nu_f_clear: fast compartment dynamics (<=1% impact)
# ===========================================================================

param_defs <- data.frame(
  name = c(
    "eps_s_HIV",    # 1. Slow reactivation: HIV+
    "eps_s_Med",    # 2. Slow reactivation: Medical
    "eps_s_FB",     # 3. Slow reactivation: Foreign-born
    "eps_s_USB",    # 4. Slow reactivation: US-born
    "eps_f_base",   # 5. Fast progression base rate
    "nu_s_clear",   # 6. Immune clearance rate
    "sigma",        # 7. Reinfection protection
    "epsilon_mix",  # 8. Assortative mixing (nativity)
    "gamma_tx",     # 9. Treatment rate
    "alpha_stab",   # 10. Post-treatment stabilization
    "VE",           # 11. Vaccine efficacy
    "duration_yr"   # 12. Vaccine duration (years) - converted to omega below
  ),
  lower = c(
    0.0041,   # HIV: Ekramnia Table 4, 95% UI lower (0.41/100 PY)
    0.00042,  # Med: ~60% below point estimate (composite stratum)
    0.00058,  # FB:  Ekramnia Table 3, non-USB 95% UI lower (0.058/100 PY)
    0.00029,  # USB: Ekramnia Table 3, USB 95% UI lower (0.029/100 PY)
    0.12,     # eps_f: ~2.5yr mean duration in L_f
    0.0005,   # nu_s_clear: near-zero clearance (half-life ~1400yr)
    0.50,     # sigma: 50% reinfection protection (lower bound)
    0.40,     # epsilon_mix: moderate assortative mixing
    1.50,     # gamma: 8-month average treatment
    0.10,     # alpha_stab: ~10yr elevated relapse risk
    0.30,     # VE: pessimistic vaccine efficacy
    5         # duration: 5-year protection (pessimistic)
  ),
  upper = c(
    0.013,    # HIV: Ekramnia Table 4, 95% UI upper (1.3/100 PY)
    0.00166,  # Med: ~60% above point estimate
    0.00097,  # FB:  Ekramnia Table 3, non-USB 95% UI upper (0.097/100 PY)
    0.0020,   # USB: conservative upper bound (truncated from Ekramnia's 0.87)
    0.28,     # eps_f: ~1.5yr mean duration in L_f
    0.005,    # nu_s_clear: meaningful clearance (half-life ~140yr)
    1.00,     # sigma: complete reinfection protection
    0.90,     # epsilon_mix: high assortative mixing
    2.50,     # gamma: 5-month average treatment
    0.40,     # alpha_stab: ~2.5yr elevated relapse risk
    0.70,     # VE: optimistic vaccine efficacy
    15        # duration: 15-year protection (optimistic)
  ),
  baseline = c(
    0.0069,   # HIV
    0.00104,  # Med
    0.00074,  # FB
    0.00068,  # USB
    0.19,     # eps_f_base
    0.002,    # nu_s_clear
    0.79,     # sigma
    0.70,     # epsilon_mix
    2.0,      # gamma
    0.20,     # alpha_stab
    0.50,     # VE
    10        # duration: 10-year protection (baseline)
  ),
  stringsAsFactors = FALSE
)

cat("LHS Parameter Ranges (shared across scenarios):\n")
cat(sprintf("  %-15s  %12s  %12s  %12s\n",
            "Parameter", "Lower", "Upper", "Baseline"))
cat(paste(rep("-", 56), collapse = ""), "\n")
for (i in 1:nrow(param_defs)) {
  cat(sprintf("  %-15s  %12.6f  %12.6f  %12.6f\n",
              param_defs$name[i], param_defs$lower[i],
              param_defs$upper[i], param_defs$baseline[i]))
}
cat("\n")

# ===========================================================================
# 1B. GENERATE LHS SAMPLE (shared across scenarios)
# ===========================================================================

N_samples <- 500   # Increase for more precision (1000+), decrease for speed
n_params  <- nrow(param_defs)

set.seed(42)  # Reproducibility — same seed used for both scenarios so draws
              # sample the same points in parameter space; only psi differs
lhs_raw <- randomLHS(N_samples, n_params)

# Scale uniform [0,1] samples to parameter ranges
lhs_params <- matrix(NA, nrow = N_samples, ncol = n_params)
colnames(lhs_params) <- param_defs$name
for (j in 1:n_params) {
  lhs_params[, j] <- param_defs$lower[j] +
    lhs_raw[, j] * (param_defs$upper[j] - param_defs$lower[j])
}

# Convert duration_yr to omega (waning rate = 1/duration)
duration_col <- which(colnames(lhs_params) == "duration_yr")
lhs_params <- cbind(lhs_params, omega = 1 / lhs_params[, duration_col])
colnames(lhs_params)[ncol(lhs_params)] <- "omega"

cat(sprintf("Generated %d LHS samples across %d parameters.\n", N_samples, n_params))
cat("  Note: duration_yr sampled uniformly [5, 15] years, converted to omega = 1/duration\n")
cat(sprintf("  Median sampled duration = %.1f years (baseline = 10 years)\n",
            median(lhs_params[, "duration_yr"])))
cat(sprintf("  Median sampled omega = %.4f (baseline = 0.10)\n\n",
            median(lhs_params[, "omega"])))

# ===========================================================================
# 1C. DEFINE STRATEGIES (theta vectors shared across scenarios)
# ===========================================================================
#
# Each strategy has:
#   - name: manuscript label (used in CSV output)
#   - theta: 4-vector of stratum coverage proportions
#   - psi_plausible: vaccination rate under the plausible scenario
#
# Under the optimistic scenario, all psi values are overridden to 0.50/yr.
# ===========================================================================

strategies <- list(
  all_high_risk = list(
    name = "All High-Risk",
    theta = c(1, 1, 1, 0),
    psi_plausible = 0.05
  ),
  plwh_nusb = list(
    name = "PLWH + Non-U.S.-Born",
    theta = c(1, 0.18, 1, 0),
    psi_plausible = 0.05
  ),
  all_mtb = list(
    name = "All Mtb-Infected",
    theta = c(1, 1, 1, 1),
    psi_plausible = 0.02
  ),
  plwh_medical = list(
    name = "PLWH + Medical",
    theta = c(1, 1, 0, 0),
    psi_plausible = 0.05
  ),
  medical = list(
    name = "Medical comorbidities",
    theta = c(0, 1, 0, 0),
    psi_plausible = 0.05
  ),
  plwh = list(
    name = "PLWH",
    theta = c(1, 0, 0, 0),
    psi_plausible = 0.10
  )
)

n_strategies <- length(strategies)

# ===========================================================================
# 1D. DEFINE SCENARIOS
# ===========================================================================
#
# Each scenario specifies how strategy psi values are determined and where
# results are written. Comment out an entry to skip that scenario.
# ===========================================================================

scenarios <- list(
  plausible = list(
    label           = "PLAUSIBLE",
    psi_override    = NULL,   # NULL = use each strategy's psi_plausible
    eps_mix_range   = NULL,   # NULL = use param_defs range as-is
    csv_ui          = "TB_LHS_UNCERTAINTY_INTERVALS.csv",
    rdata           = "TB_LHS_RESULTS.RData",
    run_prcc        = TRUE,
    csv_prcc        = "TB_PRCC_RESULTS.csv",
    prcc_outcomes   = c("cases_prevented", "NNV")  # both outcomes
  ),
  optimistic = list(
    label           = "OPTIMISTIC",
    psi_override    = 0.50,   # all strategies use this psi
    eps_mix_range   = c(0.50, 0.90),  # symmetric +/- 0.20 around 0.70 baseline.
                                      # Wider plausible range under-weights the
                                      # 0.70 baseline and shifts cases-prevented
                                      # UIs below the deterministic point estimates.
    csv_ui          = "TB_LHS_UNCERTAINTY_INTERVALS_OPTIMISTIC.csv",
    rdata           = "TB_LHS_OPTIMISTIC_RESULTS.RData",
    run_prcc        = TRUE,
    csv_prcc        = "TB_PRCC_RESULTS_OPTIMISTIC.csv",
    prcc_outcomes   = c("cases_prevented")          # cases_prevented only
  )
)

# Save ALL original globals that will be modified (shared across scenarios)
orig_eps_s       <- eps_s
orig_eps_f       <- eps_f
orig_nu_s_clear  <- nu_s_clear
orig_sigma       <- sigma
orig_gamma       <- gamma
orig_alpha_stab  <- alpha_stab
orig_omega       <- omega
orig_M_default   <- M_default
orig_epsilon_mix <- epsilon_mix


################################################################################
#
# RUN EACH SCENARIO
#
################################################################################

for (scen_name in names(scenarios)) {

  scen <- scenarios[[scen_name]]

  cat("\n\n")
  cat("################################################################\n")
  cat(sprintf("#  SCENARIO: %s\n", scen$label))
  cat("################################################################\n\n")

  # --- Resolve strategy psi values for this scenario ---
  scen_strategies <- lapply(strategies, function(s) {
    s$psi <- if (is.null(scen$psi_override)) s$psi_plausible else scen$psi_override
    s
  })

  cat("Strategies for this scenario:\n")
  for (s in names(scen_strategies)) {
    cat(sprintf("  %-30s  theta = [%g, %g, %g, %g]  psi = %.2f/yr\n",
                scen_strategies[[s]]$name,
                scen_strategies[[s]]$theta[1], scen_strategies[[s]]$theta[2],
                scen_strategies[[s]]$theta[3], scen_strategies[[s]]$theta[4],
                scen_strategies[[s]]$psi))
  }
  cat("\n")

  # --- Apply per-scenario epsilon_mix range override (if specified) ---
  # Rescales the epsilon_mix column from the global param_defs range to the
  # scenario-specific range, preserving the same [0,1] LHS draws (so all
  # other parameters sample identical points across scenarios).
  lhs_params_scen <- lhs_params
  if (!is.null(scen$eps_mix_range)) {
    eps_col <- which(colnames(lhs_params_scen) == "epsilon_mix")
    eps_idx <- which(param_defs$name == "epsilon_mix")
    raw_uniform <- (lhs_params_scen[, eps_col] - param_defs$lower[eps_idx]) /
                   (param_defs$upper[eps_idx] - param_defs$lower[eps_idx])
    new_lo <- scen$eps_mix_range[1]
    new_hi <- scen$eps_mix_range[2]
    lhs_params_scen[, eps_col] <- new_lo + raw_uniform * (new_hi - new_lo)
    cat(sprintf("Override: epsilon_mix range = [%.2f, %.2f] for this scenario\n",
                new_lo, new_hi))
    cat(sprintf("  (Default param_defs range: [%.2f, %.2f])\n\n",
                param_defs$lower[eps_idx], param_defs$upper[eps_idx]))
  }

  # =========================================================================
  # RUN LHS LOOP
  # =========================================================================

  cat("================================================================\n")
  cat(sprintf("  RUNNING LHS UNCERTAINTY ANALYSIS (%s)\n", scen$label))
  cat(sprintf("  %d samples x %d strategies = %d model runs + %d baselines\n",
              N_samples, n_strategies,
              N_samples * n_strategies, N_samples))
  cat("================================================================\n\n")

  # Pre-allocate storage
  lhs_results <- list()
  for (s in names(scen_strategies)) {
    lhs_results[[s]] <- data.frame(
      cases_prevented    = rep(NA_real_, N_samples),
      pct_reduction      = rep(NA_real_, N_samples),
      NNV                = rep(NA_real_, N_samples),
      baseline_cases     = rep(NA_real_, N_samples),
      intervention_cases = rep(NA_real_, N_samples)
    )
  }

  cal_success <- logical(N_samples)
  cal_beta    <- rep(NA_real_, N_samples)

  t_start <- Sys.time()
  progress_at <- seq(50, N_samples, by = 50)

  for (i in 1:N_samples) {

    # --- Extract this draw ---
    draw <- lhs_params_scen[i, ]

    # --- Modify globals ---
    eps_s_draw <- c(draw["eps_s_HIV"], draw["eps_s_Med"],
                    draw["eps_s_FB"],  draw["eps_s_USB"])
    names(eps_s_draw) <- NULL

    ef_base <- draw["eps_f_base"]
    eps_f_draw <- ef_base * c(2.5, 1.25, 1.0, 1.0)   # RR multipliers preserved

    assign("eps_s",      eps_s_draw,             envir = .GlobalEnv)
    assign("eps_f",      eps_f_draw,             envir = .GlobalEnv)
    assign("nu_s_clear", draw["nu_s_clear"],     envir = .GlobalEnv)
    assign("sigma",      draw["sigma"],          envir = .GlobalEnv)
    assign("gamma",      draw["gamma_tx"],       envir = .GlobalEnv)
    assign("alpha_stab", draw["alpha_stab"],     envir = .GlobalEnv)
    assign("omega",      draw["omega"],          envir = .GlobalEnv)

    # Rebuild mixing matrix with new epsilon_mix for nativity strata
    eps_vec_draw <- c(0, 0, draw["epsilon_mix"], draw["epsilon_mix"])
    assign("epsilon_mix", eps_vec_draw, envir = .GlobalEnv)
    assign("M_default", calc_mixing_matrix(N_vec, eps_vec_draw), envir = .GlobalEnv)

    # --- Recalibrate beta and mtb_prev ---
    cal_i <- tryCatch({
      calibrate_model(
        initial_beta     = baseline_cal$beta,
        initial_mtb_prev = baseline_cal$mtb_prev,
        verbose = FALSE
      )
    }, error = function(e) NULL)

    if (is.null(cal_i) || cal_i$mape > 5) {
      cal_success[i] <- FALSE
      next
    }
    cal_success[i] <- TRUE
    cal_beta[i]    <- cal_i$beta

    # --- Run baseline (no vaccination) ---
    params_base <- list(
      beta = cal_i$beta, kappa = c(1, 1, 1, 1),
      psi = 0, theta = c(0, 0, 0, 0),
      mtb_prev = cal_i$mtb_prev,
      M = M_default, N = N_vec,
      VE = draw["VE"]
    )
    base_run <- tryCatch(run_model(params_base, years = 50), error = function(e) NULL)

    if (is.null(base_run)) { cal_success[i] <- FALSE; next }

    base_cases_raw <- calc_annual_cases(base_run$final_state, params_base)$total
    sf <- actual_2024_cases / base_cases_raw   # scale factor for this draw
    base_cases_scaled <- actual_2024_cases

    # --- Run each strategy ---
    for (s in names(scen_strategies)) {
      strat <- scen_strategies[[s]]

      params_int <- list(
        beta = cal_i$beta, kappa = c(1, 1, 1, 1),
        psi = strat$psi, theta = strat$theta,
        mtb_prev = cal_i$mtb_prev,
        M = M_default, N = N_vec,
        VE = draw["VE"]
      )

      int_run <- tryCatch(run_model(params_int, years = 30), error = function(e) NULL)
      if (is.null(int_run)) next

      int_cases_raw    <- calc_annual_cases(int_run$final_state, params_int)$total
      int_cases_scaled <- int_cases_raw * sf

      prevented <- base_cases_scaled - int_cases_scaled
      pct_red   <- prevented / base_cases_scaled * 100

      vax_annual <- calc_vaccinations(int_run$final_state,
                                      strat$psi, strat$theta)$total
      nnv <- if (prevented > 0) round(vax_annual / prevented) else NA

      lhs_results[[s]]$cases_prevented[i]    <- round(prevented)
      lhs_results[[s]]$pct_reduction[i]      <- pct_red
      lhs_results[[s]]$NNV[i]                <- nnv
      lhs_results[[s]]$baseline_cases[i]     <- base_cases_scaled
      lhs_results[[s]]$intervention_cases[i] <- round(int_cases_scaled)
    }

    # --- Progress report ---
    if (i %in% progress_at) {
      elapsed   <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      rate      <- i / elapsed
      remaining <- (N_samples - i) / rate
      n_ok      <- sum(cal_success[1:i])
      cat(sprintf("  [%d/%d] %.1f min elapsed, ~%.1f min remaining | ",
                  i, N_samples, elapsed, remaining))
      cat(sprintf("beta range: [%.2f, %.2f] | Success: %d/%d\n",
                  min(cal_beta[1:i], na.rm = TRUE),
                  max(cal_beta[1:i], na.rm = TRUE),
                  n_ok, i))
    }
  }

  # --- Restore all globals ---
  assign("eps_s",       orig_eps_s,       envir = .GlobalEnv)
  assign("eps_f",       orig_eps_f,       envir = .GlobalEnv)
  assign("nu_s_clear",  orig_nu_s_clear,  envir = .GlobalEnv)
  assign("sigma",       orig_sigma,       envir = .GlobalEnv)
  assign("gamma",       orig_gamma,       envir = .GlobalEnv)
  assign("alpha_stab",  orig_alpha_stab,  envir = .GlobalEnv)
  assign("omega",       orig_omega,       envir = .GlobalEnv)
  assign("M_default",   orig_M_default,   envir = .GlobalEnv)
  assign("epsilon_mix", orig_epsilon_mix, envir = .GlobalEnv)

  elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  n_success <- sum(cal_success)
  cat(sprintf("\nLHS complete (%s): %d/%d successful (%.1f%%) in %.1f minutes\n\n",
              scen$label, n_success, N_samples,
              n_success / N_samples * 100, elapsed_total))

  # =========================================================================
  # COMPUTE 95% UNCERTAINTY INTERVALS
  # =========================================================================

  cat("================================================================\n")
  cat(sprintf("  LHS RESULTS: 95%% UNCERTAINTY INTERVALS (%s)\n", scen$label))
  cat("================================================================\n\n")

  ui_summary <- data.frame(
    Strategy         = character(),
    Prevented_median = numeric(),
    Prevented_lo     = numeric(),
    Prevented_hi     = numeric(),
    PctRed_median    = numeric(),
    PctRed_lo        = numeric(),
    PctRed_hi        = numeric(),
    NNV_median       = numeric(),
    NNV_lo           = numeric(),
    NNV_hi           = numeric(),
    N_valid          = integer(),
    stringsAsFactors = FALSE
  )

  for (s in names(scen_strategies)) {
    df <- lhs_results[[s]]
    valid <- !is.na(df$cases_prevented) & !is.na(df$NNV)
    df_v  <- df[valid, ]

    if (nrow(df_v) < 50) {
      cat(sprintf("WARNING: %s -- only %d valid runs (skipping)\n",
                  scen_strategies[[s]]$name, nrow(df_v)))
      next
    }

    prev_q <- quantile(df_v$cases_prevented, c(0.025, 0.50, 0.975))
    pct_q  <- quantile(df_v$pct_reduction,   c(0.025, 0.50, 0.975))
    nnv_q  <- quantile(df_v$NNV,             c(0.025, 0.50, 0.975))

    ui_summary <- rbind(ui_summary, data.frame(
      Strategy         = scen_strategies[[s]]$name,
      Prevented_median = round(prev_q[2]),
      Prevented_lo     = round(prev_q[1]),
      Prevented_hi     = round(prev_q[3]),
      PctRed_median    = round(pct_q[2], 1),
      PctRed_lo        = round(pct_q[1], 1),
      PctRed_hi        = round(pct_q[3], 1),
      NNV_median       = round(nnv_q[2]),
      NNV_lo           = round(nnv_q[1]),
      NNV_hi           = round(nnv_q[3]),
      N_valid          = nrow(df_v),
      stringsAsFactors = FALSE
    ))
  }

  # --- Print formatted table ---
  cat(sprintf("%-30s  %25s  %20s  %20s\n",
              "Strategy", "Cases Prevented/yr", "% Reduction", "NNV"))
  cat(sprintf("%-30s  %25s  %20s  %20s\n",
              "", "Median (95% UI)", "Median (95% UI)", "Median (95% UI)"))
  cat(paste(rep("-", 100), collapse = ""), "\n")

  for (i in 1:nrow(ui_summary)) {
    r <- ui_summary[i, ]
    cat(sprintf("%-30s  %5d (%5d - %5d)      %4.1f%% (%4.1f - %4.1f)    %4d (%4d - %4d)\n",
                r$Strategy,
                r$Prevented_median, r$Prevented_lo, r$Prevented_hi,
                r$PctRed_median, r$PctRed_lo, r$PctRed_hi,
                r$NNV_median, r$NNV_lo, r$NNV_hi))
  }

  write.csv(ui_summary, scen$csv_ui, row.names = FALSE)
  cat(sprintf("\nSaved: %s\n", scen$csv_ui))

  # =========================================================================
  # PARTIAL RANK CORRELATION COEFFICIENTS (PRCC) — plausible only
  # =========================================================================

  prcc_prev <- prcc_nnv <- NULL  # default for scenarios that skip PRCC

  if (isTRUE(scen$run_prcc)) {

    cat("\n================================================================\n")
    cat(sprintf("  PARTIAL RANK CORRELATION COEFFICIENTS (PRCC) — %s\n", scen$label))
    cat("================================================================\n\n")

    ref_strat <- "all_high_risk"
    # Validity depends on which outcomes this scenario reports
    needs_nnv <- "NNV" %in% scen$prcc_outcomes
    ref_valid <- !is.na(lhs_results[[ref_strat]]$cases_prevented)
    if (needs_nnv) {
      ref_valid <- ref_valid & !is.na(lhs_results[[ref_strat]]$NNV)
    }

    if (sum(ref_valid) > 100) {

      param_cols_for_prcc <- setdiff(colnames(lhs_params_scen), "duration_yr")
      X_valid <- lhs_params_scen[ref_valid, param_cols_for_prcc]

      compute_prcc <- function(X, y) {
        X_r <- apply(X, 2, rank)
        y_r <- rank(y)
        n_p <- ncol(X_r)
        prcc <- pval <- numeric(n_p)

        for (j in 1:n_p) {
          other <- setdiff(1:n_p, j)
          resid_y <- residuals(lm(y_r ~ X_r[, other]))
          resid_x <- residuals(lm(X_r[, j] ~ X_r[, other]))
          ct <- cor.test(resid_x, resid_y, method = "pearson")
          prcc[j] <- ct$estimate
          pval[j] <- ct$p.value
        }
        data.frame(parameter = colnames(X), PRCC = prcc, p_value = pval)
      }

      print_prcc <- function(df, label) {
        cat(sprintf("\nPRCC for %s (All High-Risk strategy):\n\n", label))
        cat(sprintf("  %-15s  %8s  %12s  %s\n",
                    "Parameter", "PRCC", "p-value", "Significance"))
        cat(paste(rep("-", 50), collapse = ""), "\n")
        for (k in 1:nrow(df)) {
          r <- df[k, ]
          sig <- ifelse(r$p_value < 0.001, "***",
                 ifelse(r$p_value < 0.01, "**",
                 ifelse(r$p_value < 0.05, "*", "")))
          cat(sprintf("  %-15s  %8.3f  %12.2e  %s\n",
                      r$parameter, r$PRCC, r$p_value, sig))
        }
      }

      prcc_tables <- list()

      if ("cases_prevented" %in% scen$prcc_outcomes) {
        y_prev <- lhs_results[[ref_strat]]$cases_prevented[ref_valid]
        prcc_prev <- compute_prcc(X_valid, y_prev)
        prcc_prev <- prcc_prev[order(-abs(prcc_prev$PRCC)), ]
        prcc_prev$outcome <- "cases_prevented"
        print_prcc(prcc_prev, "Cases Prevented")
        prcc_tables[["cases_prevented"]] <- prcc_prev
      } else {
        prcc_prev <- NULL
      }

      if ("NNV" %in% scen$prcc_outcomes) {
        y_nnv <- lhs_results[[ref_strat]]$NNV[ref_valid]
        prcc_nnv <- compute_prcc(X_valid[!is.na(y_nnv), ],
                                  y_nnv[!is.na(y_nnv)])
        prcc_nnv <- prcc_nnv[order(-abs(prcc_nnv$PRCC)), ]
        prcc_nnv$outcome <- "NNV"
        print_prcc(prcc_nnv, "NNV")
        prcc_tables[["NNV"]] <- prcc_nnv
      } else {
        prcc_nnv <- NULL
      }

      prcc_all <- do.call(rbind, prcc_tables)
      write.csv(prcc_all, scen$csv_prcc, row.names = FALSE)
      cat(sprintf("\nSaved: %s\n", scen$csv_prcc))

    } else {
      cat(sprintf("Too few valid runs (%d) for PRCC analysis.\n", sum(ref_valid)))
      prcc_prev <- prcc_nnv <- NULL
    }
  }

  # =========================================================================
  # SAVE SCENARIO RESULTS
  # =========================================================================

  save(lhs_params, lhs_params_scen, lhs_results, ui_summary, cal_success, cal_beta,
       param_defs, prcc_prev, prcc_nnv,
       scen_strategies, baseline_cal,
       file = scen$rdata)

  cat(sprintf("\nSaved: %s\n", scen$rdata))

  cat(sprintf("\n--- Scenario %s complete ---\n", scen$label))
  cat(sprintf("  Total time: %.1f minutes\n", elapsed_total))
  cat(sprintf("  Successful calibrations: %d / %d (%.1f%%)\n",
              n_success, N_samples, n_success / N_samples * 100))
}


cat("\n\n================================================================\n")
cat("  ALL SCENARIOS COMPLETE\n")
cat("================================================================\n")
cat("\n  Output files:\n")
for (scen_name in names(scenarios)) {
  scen <- scenarios[[scen_name]]
  cat(sprintf("    [%s]\n", scen$label))
  cat(sprintf("      %s\n", scen$csv_ui))
  if (isTRUE(scen$run_prcc)) cat(sprintf("      %s\n", scen$csv_prcc))
  cat(sprintf("      %s\n", scen$rdata))
}
cat("================================================================\n")
