################################################################################
#
#  TB VACCINATION MODEL — DIRECT/INDIRECT + PSI SENSITIVITY + THRESHOLD
#                         (PLAUSIBLE + OPTIMISTIC)
#
################################################################################
#
# PURPOSE:
#   Combined script bundling three sub-analyses across both scenarios:
#
#     A. Direct/indirect effect decomposition (static-FOI counterfactual)
#     B. Vaccination rate (psi) dose-response (psi 1%-20%/yr at fixed VE)
#     C. Threshold analysis (VE x duration grid + minimum VE table)
#
#   Plausible scenario:
#     - VE = 0.50 reference (sub-analyses A and B)
#     - Strategy-specific psi (each strategy's plausible default)
#
#   Optimistic scenario:
#     - VE = 0.70 reference (sub-analyses A and B)
#     - psi = 0.50/yr for all strategies (sub-analyses A and C)
#     - For sub-analysis B, psi is the swept variable; VE held at 0.70
#
# REQUIRES:
#   - deSolve package
#   - TB_VACCINATION_MODEL_COMPLETE.R in working directory
#
# USAGE:
#   1. Set working directory to folder containing the model script
#   2. Source this file: source("TB_DIRECT_INDIRECT_PSI_THRESHOLD.R")
#   3. Results saved as CSVs in working directory
#
#   To run only one scenario, comment out the unwanted entry in the
#   `scenarios` list below.
#
# OUTPUT FILES:
#   Plausible:
#     - TB_DIRECT_INDIRECT_EFFECTS.csv
#     - TB_PSI_SENSITIVITY.csv
#     - TB_THRESHOLD_ANALYSIS.csv
#     - TB_THRESHOLD_MINVE.csv
#   Optimistic:
#     - TB_DIRECT_INDIRECT_EFFECTS_OPTIMISTIC.csv
#     - TB_PSI_SENSITIVITY_OPTIMISTIC.csv
#     - TB_THRESHOLD_ANALYSIS_OPTIMISTIC.csv
#     - TB_THRESHOLD_MINVE_OPTIMISTIC.csv
#   Combined:
#     - TB_DIPT.RData
#
# EXPECTED RUNTIME: ~20-30 minutes total
#
# DATE: May 2026
#
################################################################################

suppressPackageStartupMessages({
  library(deSolve)
})

cat("================================================================\n")
cat("  TB MODEL — DIRECT/INDIRECT + PSI + THRESHOLD (BOTH SCENARIOS)\n")
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

# ---- Strategy definitions (theta vectors are scenario-independent) -------
# Per-scenario psi values are set inside the scenarios list below.
strategy_thetas <- list(
  all_high_risk = list(name = "All High-Risk (HIV+Med+FB)", theta = c(1, 1, 1, 0)),
  plwh_all_nusb = list(name = "PLWH + All NUSB",            theta = c(1, 0.18, 1, 0)),
  all_fb        = list(name = "All Foreign-Born",           theta = c(0.19, 0.18, 1, 0)),
  fb_stratum    = list(name = "FB Stratum Only",            theta = c(0, 0, 1, 0)),
  universal     = list(name = "Universal",                  theta = c(1, 1, 1, 1)),
  medical       = list(name = "Medical Only",               theta = c(0, 1, 0, 0)),
  hiv           = list(name = "HIV Only",                   theta = c(1, 0, 0, 0)),
  hiv_medical   = list(name = "HIV + Medical",              theta = c(1, 1, 0, 0))
)

# Per-scenario psi values for sub-analyses A and C
plaus_psi <- list(
  all_high_risk = 0.05, plwh_all_nusb = 0.05, all_fb = 0.05,
  fb_stratum = 0.05, universal = 0.02, medical = 0.05,
  hiv = 0.10, hiv_medical = 0.05
)

opt_psi <- list(
  all_high_risk = 0.50, plwh_all_nusb = 0.50, all_fb = 0.50,
  fb_stratum = 0.50, universal = 0.50, medical = 0.50,
  hiv = 0.50, hiv_medical = 0.50
)

# ---- Define scenarios ----------------------------------------------------
scenarios <- list(
  plausible = list(
    label    = "PLAUSIBLE",
    VE       = 0.50,
    psi      = plaus_psi,
    csv_di   = "TB_DIRECT_INDIRECT_EFFECTS.csv",
    csv_psi  = "TB_PSI_SENSITIVITY.csv",
    csv_thr  = "TB_THRESHOLD_ANALYSIS.csv",
    csv_minve = "TB_THRESHOLD_MINVE.csv"
  ),
  optimistic = list(
    label    = "OPTIMISTIC",
    VE       = 0.70,
    psi      = opt_psi,
    csv_di   = "TB_DIRECT_INDIRECT_EFFECTS_OPTIMISTIC.csv",
    csv_psi  = "TB_PSI_SENSITIVITY_OPTIMISTIC.csv",
    csv_thr  = "TB_THRESHOLD_ANALYSIS_OPTIMISTIC.csv",
    csv_minve = "TB_THRESHOLD_MINVE_OPTIMISTIC.csv"
  )
)

# ---- Compute baseline equilibrium FOI (lambda) -- shared infrastructure --
cat("Computing baseline equilibrium force of infection...\n")

params_base <- list(
  beta = baseline_cal$beta, kappa = c(1, 1, 1, 1),
  psi = 0, theta = c(0, 0, 0, 0),
  mtb_prev = baseline_cal$mtb_prev,
  M = M_default, N = N_vec, VE = 0.50
)

base_run <- run_model(params_base, years = 50)
base_state <- base_run$final_state
base_comp <- extract_compartments(base_state)

N_base <- base_comp$S + base_comp$L_f + base_comp$L_s + base_comp$L_fv +
          base_comp$L_sv + base_comp$I + base_comp$I_v + base_comp$R

lambda_baseline <- numeric(4)
for (i in 1:4) {
  foi_sum <- 0
  for (j in 1:4) {
    if (N_base[j] > 0) {
      foi_sum <- foi_sum + M_default[i, j] *
        (base_comp$I[j] + base_comp$I_v[j]) / N_base[j]
    }
  }
  lambda_baseline[i] <- baseline_cal$beta * 1.0 * foi_sum
}

cat(sprintf("  lambda_baseline = [%.6f, %.6f, %.6f, %.6f]\n\n",
            lambda_baseline[1], lambda_baseline[2],
            lambda_baseline[3], lambda_baseline[4]))

# Scaling factor
base_cases_raw <- calc_annual_cases(base_state, params_base)$total
sf <- actual_2024_cases / base_cases_raw
base_cases_scaled <- actual_2024_cases

# ---- Define static-FOI ODE system (used by sub-analysis A) ---------------
tb_model_fixed_lambda <- function(t, y, params) {
  beta   <- params$beta
  kappa  <- params$kappa
  psi    <- params$psi
  theta  <- params$theta
  mtb_prev <- params$mtb_prev
  M      <- params$M
  VE     <- params$VE

  S    <- y[c(1, 9, 17, 25)]
  L_f  <- y[c(2, 10, 18, 26)]
  L_s  <- y[c(3, 11, 19, 27)]
  L_fv <- y[c(4, 12, 20, 28)]
  L_sv <- y[c(5, 13, 21, 29)]
  I    <- y[c(6, 14, 22, 30)]
  I_v  <- y[c(7, 15, 23, 31)]
  R    <- y[c(8, 16, 24, 32)]

  N <- S + L_f + L_s + L_fv + L_sv + I + I_v + R
  lambda <- params$fixed_lambda
  N_fixed <- params$N
  Lambda <- mu * N_fixed + iota * N_fixed
  Gamma_S  <- Lambda * (1 - mtb_prev)
  Gamma_Ls <- Lambda * mtb_prev

  dS <- dL_f <- dL_s <- dL_fv <- dL_sv <- dI <- dI_v <- dR <- numeric(4)

  for (i in 1:4) {
    mu_out <- mu[i] + iota[i]
    dS[i] <- Gamma_S[i] - lambda[i] * S[i] - mu_out * S[i] +
             nu_f_clear * L_f[i] + nu_s_clear * L_s[i] +
             nu_f_clear * L_fv[i] + nu_s_clear * L_sv[i]
    dL_f[i] <- p[i] * lambda[i] * S[i] +
               p[i] * lambda[i] * R[i] * (1 - sigma) -
               eps_f[i] * L_f[i] - nu_f_trans * L_f[i] -
               nu_f_clear * L_f[i] - psi * theta[i] * L_f[i] +
               omega * L_fv[i] - mu_out * L_f[i]
    dL_s[i] <- Gamma_Ls[i] + nu_f_trans * L_f[i] + alpha_stab * R[i] +
               (1 - p[i]) * lambda[i] * S[i] +
               (1 - p[i]) * lambda[i] * R[i] * (1 - sigma) -
               eps_s[i] * L_s[i] - nu_s_clear * L_s[i] -
               psi * theta[i] * L_s[i] + omega * L_sv[i] - mu_out * L_s[i]
    dL_fv[i] <- psi * theta[i] * L_f[i] -
                (1 - VE) * eps_f[i] * L_fv[i] -
                nu_f_trans * L_fv[i] - nu_f_clear * L_fv[i] -
                omega * L_fv[i] - mu_out * L_fv[i]
    dL_sv[i] <- psi * theta[i] * L_s[i] + nu_f_trans * L_fv[i] -
                (1 - VE) * eps_s[i] * L_sv[i] -
                nu_s_clear * L_sv[i] - omega * L_sv[i] - mu_out * L_sv[i]
    dI[i] <- eps_f[i] * L_f[i] + eps_s[i] * L_s[i] + rho[i] * R[i] -
             gamma * I[i] - (mu_out + mu_TB[i]) * I[i]
    dI_v[i] <- (1 - VE) * eps_f[i] * L_fv[i] +
               (1 - VE) * eps_s[i] * L_sv[i] -
               gamma * I_v[i] - (mu_out + mu_TB[i]) * I_v[i]
    dR[i] <- gamma * I[i] + gamma * I_v[i] - alpha_stab * R[i] -
             lambda[i] * R[i] * (1 - sigma) - rho[i] * R[i] - mu_out * R[i]
  }

  dy <- numeric(32)
  for (i in 1:4) {
    idx <- (i - 1) * 8 + 1
    dy[idx]     <- dS[i]
    dy[idx + 1] <- dL_f[i]
    dy[idx + 2] <- dL_s[i]
    dy[idx + 3] <- dL_fv[i]
    dy[idx + 4] <- dL_sv[i]
    dy[idx + 5] <- dI[i]
    dy[idx + 6] <- dI_v[i]
    dy[idx + 7] <- dR[i]
  }
  list(dy)
}

run_model_fixed_lambda <- function(params, years = 30) {
  y0 <- calc_initial_conditions(N_vec, params$mtb_prev, target_inc)
  times <- seq(0, years, by = 1)
  out <- ode(y = y0, times = times, func = tb_model_fixed_lambda,
             parms = params, method = "lsoda")
  final_state <- as.numeric(out[nrow(out), -1])
  list(trajectory = out, final_state = final_state, params = params)
}

# Common grids
psi_values    <- c(0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20)
ve_grid       <- seq(0.10, 0.95, by = 0.05)
duration_grid <- c(5, 10, 15, 20, 30)
omega_grid    <- 1 / duration_grid

orig_omega <- omega

# Storage for cross-scenario results
results_all <- list()

################################################################################
# RUN EACH SCENARIO
################################################################################

for (scen_name in names(scenarios)) {

  scen <- scenarios[[scen_name]]

  cat("\n\n")
  cat("################################################################\n")
  cat(sprintf("#  SCENARIO: %s  (VE = %.0f%%)\n",
              scen$label, scen$VE * 100))
  cat("################################################################\n\n")

  # ============================================================================
  # SUB-ANALYSIS A: DIRECT VS. INDIRECT
  # ============================================================================

  cat("================================================================\n")
  cat(sprintf("  SUB-ANALYSIS A: DIRECT VS. INDIRECT (%s)\n", scen$label))
  cat("================================================================\n\n")

  decomp_results <- data.frame(
    Strategy            = character(),
    Total_prevented     = numeric(),
    Direct_prevented    = numeric(),
    Indirect_prevented  = numeric(),
    Pct_direct          = numeric(),
    Pct_indirect        = numeric(),
    Total_pct_reduction = numeric(),
    stringsAsFactors    = FALSE
  )

  for (s in names(strategy_thetas)) {
    strat <- strategy_thetas[[s]]
    psi_s <- scen$psi[[s]]
    cat(sprintf("  %-30s ... ", strat$name))

    params_dyn <- list(
      beta = baseline_cal$beta, kappa = c(1, 1, 1, 1),
      psi = psi_s, theta = strat$theta,
      mtb_prev = baseline_cal$mtb_prev,
      M = M_default, N = N_vec, VE = scen$VE
    )
    dyn_run <- run_model(params_dyn, years = 30)
    dyn_cases_scaled <- calc_annual_cases(dyn_run$final_state, params_dyn)$total * sf

    params_static <- list(
      beta = baseline_cal$beta, kappa = c(1, 1, 1, 1),
      psi = psi_s, theta = strat$theta,
      mtb_prev = baseline_cal$mtb_prev,
      M = M_default, N = N_vec, VE = scen$VE,
      fixed_lambda = lambda_baseline
    )
    static_run <- run_model_fixed_lambda(params_static, years = 30)
    static_cases_scaled <- calc_annual_cases(static_run$final_state, params_static)$total * sf

    total_prevented    <- actual_2024_cases - dyn_cases_scaled
    direct_prevented   <- actual_2024_cases - static_cases_scaled
    indirect_prevented <- static_cases_scaled - dyn_cases_scaled

    pct_direct   <- if (total_prevented > 0) direct_prevented / total_prevented * 100 else NA
    pct_indirect <- if (total_prevented > 0) indirect_prevented / total_prevented * 100 else NA
    total_pct_red <- total_prevented / actual_2024_cases * 100

    decomp_results <- rbind(decomp_results, data.frame(
      Strategy            = strat$name,
      Total_prevented     = round(total_prevented),
      Direct_prevented    = round(direct_prevented),
      Indirect_prevented  = round(indirect_prevented),
      Pct_direct          = round(pct_direct, 1),
      Pct_indirect        = round(pct_indirect, 1),
      Total_pct_reduction = round(total_pct_red, 1),
      stringsAsFactors    = FALSE
    ))
    cat("done\n")
  }

  write.csv(decomp_results, scen$csv_di, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", scen$csv_di))

  # ============================================================================
  # SUB-ANALYSIS B: PSI DOSE-RESPONSE
  # ============================================================================

  cat(sprintf("\n================================================================\n"))
  cat(sprintf("  SUB-ANALYSIS B: PSI DOSE-RESPONSE (VE = %.0f%%)\n", scen$VE * 100))
  cat(sprintf("================================================================\n\n"))

  psi_all_results <- data.frame(
    strategy        = character(),
    psi_pct         = numeric(),
    annual_cases    = numeric(),
    cases_prevented = numeric(),
    pct_reduction   = numeric(),
    annual_vax      = numeric(),
    NNV             = numeric(),
    marginal_NNV    = numeric(),
    stringsAsFactors = FALSE
  )

  for (s in names(strategy_thetas)) {
    strat <- strategy_thetas[[s]]
    cat(sprintf("  psi sweep: %-30s ... ", strat$name))

    prev_cases <- base_cases_raw
    prev_vax   <- 0

    for (psi_val in psi_values) {
      params <- list(
        beta = baseline_cal$beta, kappa = c(1, 1, 1, 1),
        psi = psi_val, theta = strat$theta,
        mtb_prev = baseline_cal$mtb_prev,
        M = M_default, N = N_vec, VE = scen$VE
      )

      result <- tryCatch(run_model(params, years = 30), error = function(e) NULL)
      if (is.null(result)) next

      cases_raw <- calc_annual_cases(result$final_state, params)$total
      cases_scaled <- cases_raw * sf
      vax <- calc_vaccinations(result$final_state, psi_val, strat$theta)$total

      prevented <- base_cases_scaled - cases_scaled
      pct_red   <- prevented / base_cases_scaled * 100
      nnv       <- if (prevented > 0) round(vax / prevented) else NA

      marginal_prev <- (prev_cases - cases_raw) * sf
      marginal_vax  <- vax - prev_vax
      marg_nnv <- if (!is.na(marginal_prev) && marginal_prev > 0) {
        round(marginal_vax / marginal_prev)
      } else { NA }

      psi_all_results <- rbind(psi_all_results, data.frame(
        strategy        = strat$name,
        psi_pct         = psi_val * 100,
        annual_cases    = round(cases_scaled),
        cases_prevented = round(prevented),
        pct_reduction   = round(pct_red, 2),
        annual_vax      = round(vax),
        NNV             = nnv,
        marginal_NNV    = marg_nnv,
        stringsAsFactors = FALSE
      ))

      prev_cases <- cases_raw
      prev_vax   <- vax
    }
    cat("done\n")
  }

  write.csv(psi_all_results, scen$csv_psi, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", scen$csv_psi))

  # ============================================================================
  # SUB-ANALYSIS C: THRESHOLD (VE x DURATION GRID)
  # ============================================================================

  cat(sprintf("\n================================================================\n"))
  cat(sprintf("  SUB-ANALYSIS C: THRESHOLD ANALYSIS (%s)\n", scen$label))
  cat(sprintf("================================================================\n\n"))

  grid_results <- data.frame(
    strategy        = character(),
    VE              = numeric(),
    duration_yr     = numeric(),
    cases_prevented = numeric(),
    pct_reduction   = numeric(),
    NNV             = numeric(),
    stringsAsFactors = FALSE
  )

  total_runs <- length(strategy_thetas) * length(ve_grid) * length(duration_grid)
  cat(sprintf("Running %d scenarios (%d strategies x %d VE x %d durations)...\n\n",
              total_runs, length(strategy_thetas),
              length(ve_grid), length(duration_grid)))

  for (s in names(strategy_thetas)) {
    strat <- strategy_thetas[[s]]
    psi_s <- scen$psi[[s]]
    cat(sprintf("  %-30s ... ", strat$name))

    for (d_idx in seq_along(duration_grid)) {
      dur <- duration_grid[d_idx]
      om  <- omega_grid[d_idx]
      assign("omega", om, envir = .GlobalEnv)

      for (ve in ve_grid) {
        params <- list(
          beta = baseline_cal$beta, kappa = c(1, 1, 1, 1),
          psi = psi_s, theta = strat$theta,
          mtb_prev = baseline_cal$mtb_prev,
          M = M_default, N = N_vec, VE = ve
        )

        result <- tryCatch(run_model(params, years = 30), error = function(e) NULL)
        if (is.null(result)) next

        cases_raw <- calc_annual_cases(result$final_state, params)$total
        cases_scaled <- cases_raw * sf
        vax <- calc_vaccinations(result$final_state, psi_s, strat$theta)$total

        prevented <- base_cases_scaled - cases_scaled
        pct_red   <- prevented / base_cases_scaled * 100
        nnv       <- if (prevented > 0) round(vax / prevented) else NA

        grid_results <- rbind(grid_results, data.frame(
          strategy        = strat$name,
          VE              = ve * 100,
          duration_yr     = dur,
          cases_prevented = round(prevented),
          pct_reduction   = round(pct_red, 2),
          NNV             = nnv,
          stringsAsFactors = FALSE
        ))
      }
    }
    cat("done\n")
  }

  # Restore omega after this scenario's threshold runs
  assign("omega", orig_omega, envir = .GlobalEnv)

  # Extract minimum-VE table
  thresholds <- data.frame(
    strategy    = character(),
    duration_yr = numeric(),
    target      = character(),
    min_VE_pct  = numeric(),
    stringsAsFactors = FALSE
  )

  targets <- list(
    ">=5% reduction"  = function(df) df$pct_reduction >= 5,
    ">=10% reduction" = function(df) df$pct_reduction >= 10,
    "NNV < 500"       = function(df) !is.na(df$NNV) & df$NNV < 500,
    "NNV < 1000"      = function(df) !is.na(df$NNV) & df$NNV < 1000
  )

  for (s in names(strategy_thetas)) {
    sname <- strategy_thetas[[s]]$name

    for (dur in duration_grid) {
      sub <- grid_results[grid_results$strategy == sname &
                          grid_results$duration_yr == dur, ]
      if (nrow(sub) == 0) next

      for (tname in names(targets)) {
        meets <- targets[[tname]](sub)
        if (any(meets)) {
          min_ve <- min(sub$VE[meets])
        } else {
          min_ve <- NA
        }
        thresholds <- rbind(thresholds, data.frame(
          strategy    = sname,
          duration_yr = dur,
          target      = tname,
          min_VE_pct  = min_ve,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  write.csv(grid_results, scen$csv_thr,   row.names = FALSE)
  write.csv(thresholds,   scen$csv_minve, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", scen$csv_thr))
  cat(sprintf("  Saved: %s\n", scen$csv_minve))

  # Store for RData
  results_all[[scen_name]] <- list(
    decomp_results  = decomp_results,
    psi_all_results = psi_all_results,
    grid_results    = grid_results,
    thresholds      = thresholds
  )
}

# ---- Save R objects ------------------------------------------------------
save(results_all, lambda_baseline, strategy_thetas,
     plaus_psi, opt_psi,
     file = "TB_DIPT.RData")

cat("\n\n================================================================\n")
cat("  ALL SUB-ANALYSES COMPLETE (BOTH SCENARIOS)\n")
cat("================================================================\n")
cat("  Saved: TB_DIPT.RData\n")
cat("================================================================\n")
