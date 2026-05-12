################################################################################
#
#  TB VACCINATION MODEL — CALIBRATION PARAMETER EXTRACTION
#
################################################################################
#
# PURPOSE:
#   Dumps the calibrated parameter set and target incidence values to CSV.
#   Scenario-independent (calibration depends only on the natural-history
#   parameters and 2024 surveillance targets, not on vaccination strategy).
#
# REQUIRES:
#   - TB_VACCINATION_MODEL_COMPLETE.R in working directory
#
# USAGE:
#   1. Set working directory to folder containing the model script
#   2. Source this file: source("TB_EXTRACT_CALIBRATION.R")
#   3. Result saved as CSV in working directory
#
# OUTPUT FILES:
#   - TB_CALIBRATION_RESULTS.csv
#
# EXPECTED RUNTIME: < 30 seconds (just dumps the calibration object)
#
# DATE: May 2026
#
################################################################################

cat("================================================================\n")
cat("  TB VACCINATION MODEL — CALIBRATION EXTRACTION\n")
cat("================================================================\n\n")

# --- Source the base model (runs calibration automatically) ---
cat("Loading base model and running calibration...\n")
invisible(capture.output(suppressMessages(
  source("TB_VACCINATION_MODEL_COMPLETE.R")
)))
cat("Base model loaded.\n")
cat(sprintf("  Baseline beta = %.4f\n", calibration$beta))
cat(sprintf("  MAPE = %.2f%%\n\n", calibration$mape))

# ---- Build calibration table ---------------------------------------------
calibration_df <- data.frame(
  Parameter = c(
    "beta",
    "mtb_prev_HIV", "mtb_prev_Medical", "mtb_prev_FB", "mtb_prev_USB",
    "target_inc_HIV", "target_inc_Medical", "target_inc_FB", "target_inc_USB",
    "N_HIV", "N_Medical", "N_FB", "N_USB",
    "eps_s_HIV", "eps_s_Medical", "eps_s_FB", "eps_s_USB",
    "eps_f_HIV", "eps_f_Medical", "eps_f_FB", "eps_f_USB",
    "gamma", "omega", "VE_reference", "sigma", "alpha_stab",
    "nu_f_trans", "nu_s_clear"
  ),
  Value = c(
    calibration$beta,
    calibration$mtb_prev[1], calibration$mtb_prev[2],
    calibration$mtb_prev[3], calibration$mtb_prev[4],
    target_inc[1], target_inc[2], target_inc[3], target_inc[4],
    N_vec[1], N_vec[2], N_vec[3], N_vec[4],
    eps_s[1], eps_s[2], eps_s[3], eps_s[4],
    eps_f[1], eps_f[2], eps_f[3], eps_f[4],
    gamma, omega, VE_d, sigma, alpha_stab,
    nu_f_trans, nu_s_clear
  ),
  Description = c(
    "Calibrated transmission rate",
    "Calibrated Mtb infection prevalence - HIV",
    "Calibrated Mtb infection prevalence - Medical",
    "Calibrated Mtb infection prevalence - Non-U.S.-Born",
    "Calibrated Mtb infection prevalence - U.S.-Born",
    "Target incidence per 100k - HIV",
    "Target incidence per 100k - Medical",
    "Target incidence per 100k - Non-U.S.-Born",
    "Target incidence per 100k - U.S.-Born",
    "Population size - HIV",
    "Population size - Medical",
    "Population size - Non-U.S.-Born",
    "Population size - U.S.-Born",
    "Slow reactivation rate - HIV",
    "Slow reactivation rate - Medical",
    "Slow reactivation rate - Non-U.S.-Born",
    "Slow reactivation rate - U.S.-Born",
    "Fast progression rate - HIV",
    "Fast progression rate - Medical",
    "Fast progression rate - Non-U.S.-Born",
    "Fast progression rate - U.S.-Born",
    "Treatment rate (/yr)",
    "Vaccine waning rate (/yr)",
    "Reference vaccine efficacy",
    "Reinfection protection",
    "Post-treatment stabilization rate (/yr)",
    "Lf to Ls transition rate (/yr)",
    "LTBI clearance rate (/yr)"
  ),
  stringsAsFactors = FALSE
)

write.csv(calibration_df, "TB_CALIBRATION_RESULTS.csv", row.names = FALSE)

cat("================================================================\n")
cat("  CALIBRATION EXTRACTION COMPLETE\n")
cat("================================================================\n")
cat(sprintf("  Wrote %d parameter rows\n", nrow(calibration_df)))
cat("  Output: TB_CALIBRATION_RESULTS.csv\n")
cat("================================================================\n")
