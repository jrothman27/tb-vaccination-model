# ============================================================================
#
#        CLEAN ALTERNATIVE CALIBRATION: Fix mtb_prev, Calibrate epsilon_s
#
# ============================================================================
#
# PURPOSE:
#   This script performs a "clean" alternative calibration that:
#     1. Fixes mtb_prev at Aim 1 ecological estimates for each stratum
#     2. Calibrates epsilon_s (and beta) to match incidence targets
#
#   This is conceptually cleaner than the multi-target optimization because:
#     - Each stratum gets exactly its Aim 1 prevalence (no residual allocation)
#     - The rate-pool tradeoff is tested directly without optimization artifacts
#     - USB gets its correct Aim 1 value (0.79%), not the squeezed value (0.01%)
#
# COMPARISON:
#   Primary calibration:    Fix eps_s (Ekramnia), calibrate mtb_prev -> larger pools
#   Clean alternative:      Fix mtb_prev (Aim 1),  calibrate eps_s   -> higher rates
#
# EXPECTED RESULTS:
#   Stratum   Pool Ratio (Primary/Aim1)   Expected eps_s increase
#   HIV       2.41x                        ~2.4x
#   Medical   1.76x                        ~1.8x
#   FB        1.18x                        ~1.2x
#   USB       1.19x                        ~1.2x (NOT 2.94x as in multi-target!)
#
# USAGE:
#   1. Source the main model file first: source("TB_VACCINATION_MODEL_COMPLETE.R")
#   2. Source this file: source("TB_CLEAN_ALTERNATIVE_CALIBRATION.R")
#   3. Results will be saved to CSV files
#
# ============================================================================

cat("\n")
cat("======================================================================\n")
cat("  CLEAN ALTERNATIVE CALIBRATION\n")
cat("  Fix mtb_prev at Aim 1 estimates, calibrate epsilon_s to incidence\n")
cat("======================================================================\n\n")

# ============================================================================
# SECTION 1: AIM 1 ECOLOGICAL ESTIMATES
# ============================================================================
#
# These are the Mtb infection prevalence estimates from Aim 1, calculated as:
#   LTBI_i = (Cases_i Ã— p_reactivation) / eps_s_i
#
# Source: Model guide v4.4, derived from Ekramnia et al. rates and NTSS cases
# ============================================================================

# Aim 1 prevalence estimates (as fractions)
aim1_mtb_prev <- c(
  0.0431,   # HIV:     4.31% (51,707 / 1.2M)
  0.0988,   # Medical: 9.88% (3.95M / 40M)
  0.2110,   # FB:      21.1% (using gross FB estimate)
  0.0079    # USB:     0.79% (2.27M / 287.4M, using gross USB)
)

# Primary calibration mtb_prev (for comparison)
primary_mtb_prev <- c(0.1117, 0.1807, 0.3014, 0.0086)

# Primary eps_s (Ekramnia values)
primary_eps_s <- c(0.0069, 0.00104, 0.00074, 0.00068)

cat("Aim 1 mtb_prev targets:\n")
cat(sprintf("  HIV:     %.2f%% (vs primary %.2f%%, ratio %.2fx)\n", 
            aim1_mtb_prev[1]*100, primary_mtb_prev[1]*100, 
            primary_mtb_prev[1]/aim1_mtb_prev[1]))
cat(sprintf("  Medical: %.2f%% (vs primary %.2f%%, ratio %.2fx)\n", 
            aim1_mtb_prev[2]*100, primary_mtb_prev[2]*100,
            primary_mtb_prev[2]/aim1_mtb_prev[2]))
cat(sprintf("  FB:      %.2f%% (vs primary %.2f%%, ratio %.2fx)\n", 
            aim1_mtb_prev[3]*100, primary_mtb_prev[3]*100,
            primary_mtb_prev[3]/aim1_mtb_prev[3]))
cat(sprintf("  USB:     %.2f%% (vs primary %.2f%%, ratio %.2fx)\n\n", 
            aim1_mtb_prev[4]*100, primary_mtb_prev[4]*100,
            primary_mtb_prev[4]/aim1_mtb_prev[4]))

# ============================================================================
# SECTION 2: CALIBRATION FUNCTION - FIX mtb_prev, CALIBRATE eps_s
# ============================================================================

calibrate_alternative <- function(fixed_mtb_prev = aim1_mtb_prev,
                                   initial_beta = 5.0,
                                   initial_eps_s = primary_eps_s,
                                   verbose = TRUE) {
  #' Calibrate eps_s and beta while fixing mtb_prev at Aim 1 values
  #'
  #' @param fixed_mtb_prev Vector of 4 LTBI prevalences (fixed, not calibrated)
  #' @param initial_beta Starting value for transmission rate
  #' @param initial_eps_s Starting values for slow reactivation rates
  #' @param verbose Print progress messages
  #' @return List with calibrated parameters and fit statistics
  
  if (verbose) {
    cat("=======================================================================\n")
    cat("ALTERNATIVE CALIBRATION: Fix mtb_prev, calibrate eps_s\n")
    cat("=======================================================================\n\n")
    cat("Fixed parameters:   mtb_prev (at Aim 1 ecological estimates)\n")
    cat("Free parameters:    beta, eps_s[1:4] (5 total)\n")
    cat("Targets:            Stratum-specific incidence (4 targets)\n\n")
  }
  
  # Store original eps_s to restore later
  orig_eps_s <- eps_s
  orig_eps_f <- eps_f
  
  # --- Objective function ---
  # Calibrate beta and eps_s to match incidence, with mtb_prev FIXED
  objective <- function(theta) {
    beta <- exp(theta[1])
    new_eps_s <- exp(theta[2:5])
    
    # Bounds check for beta
    if (beta < 0.1 || beta > 100) return(1e10)
    
    # Bounds check for eps_s (allow up to 10x Ekramnia values)
    if (any(new_eps_s < 0.0001) || any(new_eps_s > 0.10)) return(1e10)
    
    # Update global eps_s for model run
    # (This is a bit hacky but necessary given model structure)
    assign("eps_s", new_eps_s, envir = .GlobalEnv)
    
    # Also update eps_f proportionally to maintain fast/slow ratio
    # (Optional: could also calibrate eps_f, but keep it simple)
    eps_f_ratio <- new_eps_s / orig_eps_s
    new_eps_f <- orig_eps_f * sqrt(eps_f_ratio)  # Scale eps_f more gently
    assign("eps_f", new_eps_f, envir = .GlobalEnv)
    
    # Create parameter list with FIXED mtb_prev
    params <- list(
      beta = beta,
      kappa = c(1, 1, 1, 1),
      psi = 0,
      theta = c(0, 0, 0, 0),
      mtb_prev = fixed_mtb_prev,  # FIXED at Aim 1 values
      M = M_default, 
      N = N_vec,
      VE = VE_d
    )
    
    # Run model
    tryCatch({
      result <- run_model(params, years = 50)
      model_inc <- calc_incidence(result$final_state)$stratum
      
      if (any(is.na(model_inc)) || any(model_inc <= 0)) return(1e10)
      
      # Sum of squared relative errors
      rel_errors <- (model_inc - target_inc) / target_inc
      sse <- sum(rel_errors^2)
      
      return(sse)
    }, error = function(e) return(1e10))
  }
  
  # --- Run optimization ---
  theta0 <- log(c(initial_beta, initial_eps_s))
  
  if (verbose) cat("Running optimization...\n")
  
  result <- optim(
    par = theta0,
    fn = objective,
    method = "L-BFGS-B",
    lower = log(c(0.1, rep(0.0001, 4))),
    upper = log(c(100, rep(0.10, 4))),
    control = list(maxit = 1000, factr = 1e7)
  )
  
  # --- Extract calibrated parameters ---
  beta_cal <- exp(result$par[1])
  eps_s_cal <- exp(result$par[2:5])
  
  # Update global eps_s with calibrated values for validation
  assign("eps_s", eps_s_cal, envir = .GlobalEnv)
  eps_f_ratio <- eps_s_cal / orig_eps_s
  eps_f_cal <- orig_eps_f * sqrt(eps_f_ratio)
  assign("eps_f", eps_f_cal, envir = .GlobalEnv)
  
  # --- Validate calibration ---
  params_cal <- list(
    beta = beta_cal,
    kappa = c(1, 1, 1, 1),
    psi = 0,
    theta = c(0, 0, 0, 0),
    mtb_prev = fixed_mtb_prev,
    M = M_default, 
    N = N_vec,
    VE = VE_d
  )
  
  final_result <- run_model(params_cal, years = 50)
  final_inc <- calc_incidence(final_result$final_state)
  final_cases <- calc_annual_cases(final_result$final_state, params_cal)
  
  # Extract compartment sizes
  comp <- extract_compartments(final_result$final_state)
  total_Ls <- sum(comp$L_s)
  total_Lf <- sum(comp$L_f)
  total_LTBI <- total_Ls + total_Lf
  
  # Calculate fit statistics
  abs_errors <- abs(final_inc$stratum - target_inc)
  rel_errors_pct <- abs_errors / target_inc * 100
  mape <- mean(rel_errors_pct)
  
  if (verbose) {
    cat(sprintf("\nConvergence: %s\n", 
                ifelse(result$convergence == 0, "Successful", "Check results")))
    cat(sprintf("Final objective value: %.2e\n", result$value))
    cat(sprintf("Function evaluations: %d\n\n", result$counts[1]))
    
    cat("CALIBRATED PARAMETERS:\n")
    cat(sprintf("  beta = %.4f (primary: %.4f)\n", beta_cal, 5.699))
    cat("\n  Reactivation rates (eps_s):\n")
    cat(sprintf("  %-10s  %12s  %12s  %8s\n", 
                "Stratum", "Calibrated", "Primary", "Ratio"))
    cat("  ------------------------------------------------\n")
    for (i in 1:4) {
      ratio <- eps_s_cal[i] / orig_eps_s[i]
      cat(sprintf("  %-10s  %12.5f  %12.5f  %7.2fx\n",
                  stratum_names[i], eps_s_cal[i], orig_eps_s[i], ratio))
    }
    
    cat("\n  Fast progression rates (eps_f):\n")
    cat(sprintf("  %-10s  %12s  %12s  %8s\n", 
                "Stratum", "Calibrated", "Primary", "Ratio"))
    cat("  ------------------------------------------------\n")
    for (i in 1:4) {
      ratio <- eps_f_cal[i] / orig_eps_f[i]
      cat(sprintf("  %-10s  %12.5f  %12.5f  %7.2fx\n",
                  stratum_names[i], eps_f_cal[i], orig_eps_f[i], ratio))
    }
    
    cat("\nINCIDENCE FIT:\n")
    cat(sprintf("  %-10s  %8s  %8s  %10s\n", "Stratum", "Model", "Target", "Error"))
    cat("  ------------------------------------------\n")
    for (i in 1:4) {
      cat(sprintf("  %-10s  %8.2f  %8.2f  %9.2f%%\n",
                  stratum_names[i], final_inc$stratum[i], target_inc[i], rel_errors_pct[i]))
    }
    cat(sprintf("  MAPE: %.4f%%\n", mape))
    
    cat("\nLTBI POOL SIZES:\n")
    cat(sprintf("  %-10s  %14s  %14s\n", "Stratum", "L_s", "LTBI Prev"))
    cat("  ------------------------------------------\n")
    for (i in 1:4) {
      prev <- (comp$L_s[i] + comp$L_f[i]) / N_vec[i] * 100
      cat(sprintf("  %-10s  %14.0f  %13.2f%%\n",
                  stratum_names[i], comp$L_s[i], prev))
    }
    cat(sprintf("  Total LTBI: %.2f million\n", total_LTBI / 1e6))
    
    cat("\nCASE SOURCE DISTRIBUTION:\n")
    cat(sprintf("  From L_f (recent):     %.1f%%\n", final_cases$pct_Lf))
    cat(sprintf("  From L_s (reactivation): %.1f%%\n", final_cases$pct_Ls))
    cat(sprintf("  From R (relapse):      %.1f%%\n", final_cases$pct_relapse))
  }
  
  # Restore original eps_s and eps_f
  assign("eps_s", orig_eps_s, envir = .GlobalEnv)
  assign("eps_f", orig_eps_f, envir = .GlobalEnv)
  
  return(list(
    beta = beta_cal,
    eps_s = eps_s_cal,
    eps_f = eps_f_cal,
    mtb_prev = fixed_mtb_prev,
    params = params_cal,
    mape = mape,
    model_inc = final_inc$stratum,
    compartments = comp,
    total_LTBI = total_LTBI,
    pct_fast = final_cases$pct_Lf,
    pct_slow = final_cases$pct_Ls,
    pct_relapse = final_cases$pct_relapse,
    baseline_cases = final_cases$total,
    convergence = result$convergence
  ))
}


# ============================================================================
# SECTION 3: RUN ALTERNATIVE CALIBRATION
# ============================================================================

cat("\n")
cat("======================================================================\n")
cat("  RUNNING CLEAN ALTERNATIVE CALIBRATION\n")
cat("======================================================================\n\n")

# Run the alternative calibration
alt_cal <- calibrate_alternative(
  fixed_mtb_prev = aim1_mtb_prev,
  initial_beta = 7.0,  # Start higher since we expect higher beta
  initial_eps_s = primary_eps_s * 1.5,  # Start at 1.5x Ekramnia
  verbose = TRUE
)


# ============================================================================
# SECTION 4: RUN PRIMARY CALIBRATION FOR COMPARISON
# ============================================================================

cat("\n")
cat("======================================================================\n")
cat("  RUNNING PRIMARY CALIBRATION (for comparison)\n")
cat("======================================================================\n\n")

primary_cal <- calibrate_model(verbose = TRUE)


# ============================================================================
# SECTION 5: RUN INTERVENTION SCENARIOS WITH BOTH CALIBRATIONS
# ============================================================================

cat("\n")
cat("======================================================================\n")
cat("  COMPARING INTERVENTION RESULTS\n")
cat("======================================================================\n\n")

# Define key scenarios to test
test_scenarios <- list(
  list(name = "All High-Risk", psi = 0.05, theta = c(1, 1, 1, 0)),
  list(name = "All Foreign-Born", psi = 0.05, theta = c(0.19, 0.18, 1, 0)),
  list(name = "FB Stratum Only", psi = 0.05, theta = c(0, 0, 1, 0)),
  list(name = "Universal", psi = 0.02, theta = c(1, 1, 1, 1)),
  list(name = "Medical Only", psi = 0.05, theta = c(0, 1, 0, 0)),
  list(name = "HIV Only", psi = 0.10, theta = c(1, 0, 0, 0))
)

# Function to run a scenario and return results
run_scenario <- function(cal_result, scenario, cal_type = "primary") {
  
  # Set eps_s based on calibration type
  if (cal_type == "alternative") {
    assign("eps_s", alt_cal$eps_s, envir = .GlobalEnv)
    assign("eps_f", alt_cal$eps_f, envir = .GlobalEnv)
    mtb_prev_use <- alt_cal$mtb_prev
    beta_use <- alt_cal$beta
  } else {
    # Primary uses global defaults
    mtb_prev_use <- primary_cal$mtb_prev
    beta_use <- primary_cal$beta
  }
  
  # Baseline (no vaccination)
  params_base <- list(
    beta = beta_use,
    kappa = c(1, 1, 1, 1),
    psi = 0,
    theta = c(0, 0, 0, 0),
    mtb_prev = mtb_prev_use,
    M = M_default, N = N_vec,
    VE = 0.50
  )
  
  result_base <- run_model(params_base, years = 50)
  cases_base <- calc_annual_cases(result_base$final_state, params_base)
  
  # Intervention
  params_int <- list(
    beta = beta_use,
    kappa = c(1, 1, 1, 1),
    psi = scenario$psi,
    theta = scenario$theta,
    mtb_prev = mtb_prev_use,
    M = M_default, N = N_vec,
    VE = 0.50
  )
  
  result_int <- run_model(params_int, years = 50)
  cases_int <- calc_annual_cases(result_int$final_state, params_int)
  vax <- calc_vaccinations(result_int$final_state, scenario$psi, scenario$theta)
  
  # Calculate impact
  prevented <- cases_base$total - cases_int$total
  pct_reduction <- prevented / cases_base$total * 100
  NNV <- if (prevented > 0) vax$total / prevented else NA
  
  # Restore original eps_s
  assign("eps_s", primary_eps_s, envir = .GlobalEnv)
  assign("eps_f", c(0.475, 0.2375, 0.19, 0.19), envir = .GlobalEnv)
  
  return(list(
    strategy = scenario$name,
    calibration = cal_type,
    baseline_cases = cases_base$total,
    intervention_cases = cases_int$total,
    cases_prevented = prevented,
    pct_reduction = pct_reduction,
    annual_vaccinations = vax$total,
    NNV = NNV
  ))
}

# Run all scenarios with both calibrations
results <- data.frame()

for (scenario in test_scenarios) {
  # Primary calibration
  res_primary <- run_scenario(primary_cal, scenario, "primary")
  results <- rbind(results, as.data.frame(res_primary))
  
  # Alternative calibration
  res_alt <- run_scenario(alt_cal, scenario, "alternative")
  results <- rbind(results, as.data.frame(res_alt))
}

# ============================================================================
# SECTION 6: DISPLAY AND SAVE RESULTS
# ============================================================================

cat("\nINTERVENTION COMPARISON: Primary vs Alternative Calibration\n")
cat("(VE = 50%, 10-year protection)\n")
cat("======================================================================\n\n")

# Reshape for side-by-side comparison
primary_results <- results[results$calibration == "primary", ]
alt_results <- results[results$calibration == "alternative", ]

cat(sprintf("%-20s | %12s %8s %6s | %12s %8s %6s | %6s\n",
            "Strategy", "Prim Prev", "Prim NNV", "Prim %", 
            "Alt Prev", "Alt NNV", "Alt %", "Ratio"))
cat(paste(rep("-", 95), collapse = ""), "\n")

for (i in 1:nrow(primary_results)) {
  ratio <- alt_results$cases_prevented[i] / primary_results$cases_prevented[i]
  cat(sprintf("%-20s | %12.0f %8.0f %5.1f%% | %12.0f %8.0f %5.1f%% | %5.2fx\n",
              primary_results$strategy[i],
              primary_results$cases_prevented[i],
              primary_results$NNV[i],
              primary_results$pct_reduction[i],
              alt_results$cases_prevented[i],
              alt_results$NNV[i],
              alt_results$pct_reduction[i],
              ratio))
}

# ============================================================================
# SECTION 7: SAVE RESULTS TO CSV
# ============================================================================

# Save calibration comparison
calibration_comparison <- data.frame(
  Parameter = c("beta", 
                "eps_s_HIV", "eps_s_Medical", "eps_s_FB", "eps_s_USB",
                "mtb_prev_HIV", "mtb_prev_Medical", "mtb_prev_FB", "mtb_prev_USB",
                "Total_LTBI_millions", "MAPE_pct"),
  Primary = c(primary_cal$beta,
              primary_eps_s[1], primary_eps_s[2], primary_eps_s[3], primary_eps_s[4],
              primary_cal$mtb_prev[1], primary_cal$mtb_prev[2], 
              primary_cal$mtb_prev[3], primary_cal$mtb_prev[4],
              21.7, primary_cal$mape),
  Alternative = c(alt_cal$beta,
                  alt_cal$eps_s[1], alt_cal$eps_s[2], alt_cal$eps_s[3], alt_cal$eps_s[4],
                  alt_cal$mtb_prev[1], alt_cal$mtb_prev[2],
                  alt_cal$mtb_prev[3], alt_cal$mtb_prev[4],
                  alt_cal$total_LTBI / 1e6, alt_cal$mape),
  Ratio = c(alt_cal$beta / primary_cal$beta,
            alt_cal$eps_s[1] / primary_eps_s[1],
            alt_cal$eps_s[2] / primary_eps_s[2],
            alt_cal$eps_s[3] / primary_eps_s[3],
            alt_cal$eps_s[4] / primary_eps_s[4],
            alt_cal$mtb_prev[1] / primary_cal$mtb_prev[1],
            alt_cal$mtb_prev[2] / primary_cal$mtb_prev[2],
            alt_cal$mtb_prev[3] / primary_cal$mtb_prev[3],
            alt_cal$mtb_prev[4] / primary_cal$mtb_prev[4],
            (alt_cal$total_LTBI / 1e6) / 21.7,
            NA)
)

write.csv(calibration_comparison, "TB_CLEAN_ALTERNATIVE_CALIBRATION_PARAMS.csv", 
          row.names = FALSE)
cat("\nSaved calibration comparison to: TB_CLEAN_ALTERNATIVE_CALIBRATION_PARAMS.csv\n")

# Save intervention comparison
write.csv(results, "TB_CLEAN_ALTERNATIVE_INTERVENTION_RESULTS.csv", 
          row.names = FALSE)
cat("Saved intervention results to: TB_CLEAN_ALTERNATIVE_INTERVENTION_RESULTS.csv\n")

# Save summary comparison table
summary_table <- data.frame(
  Strategy = primary_results$strategy,
  Primary_Prevented = round(primary_results$cases_prevented),
  Primary_NNV = round(primary_results$NNV),
  Primary_PctRed = round(primary_results$pct_reduction, 1),
  Alt_Prevented = round(alt_results$cases_prevented),
  Alt_NNV = round(alt_results$NNV),
  Alt_PctRed = round(alt_results$pct_reduction, 1),
  Prevented_Ratio = round(alt_results$cases_prevented / primary_results$cases_prevented, 2),
  NNV_Ratio = round(primary_results$NNV / alt_results$NNV, 2)
)

write.csv(summary_table, "TB_CLEAN_ALTERNATIVE_SUMMARY.csv", row.names = FALSE)
cat("Saved summary table to: TB_CLEAN_ALTERNATIVE_SUMMARY.csv\n")

cat("\n")
cat("======================================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("======================================================================\n")
cat("\nKey findings:\n")
cat(sprintf("  - Primary calibration total LTBI: %.1f million\n", 21.7))
cat(sprintf("  - Alternative calibration total LTBI: %.1f million\n", alt_cal$total_LTBI / 1e6))
cat(sprintf("  - Pool ratio: %.2fx\n", 21.7 / (alt_cal$total_LTBI / 1e6)))
cat("\neps_s ratios (Alternative / Primary):\n")
for (i in 1:4) {
  cat(sprintf("  %s: %.2fx\n", stratum_names[i], alt_cal$eps_s[i] / primary_eps_s[i]))
}
cat("\nStrategy rankings should be IDENTICAL between calibrations.\n")
cat("NNV should be LOWER (more favorable) under alternative calibration.\n")
