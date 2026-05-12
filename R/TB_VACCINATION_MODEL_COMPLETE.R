# ============================================================================
#
#                   TB VACCINATION MODEL - COMPLETE IMPLEMENTATION
#                        WITH INTERVENTION SCENARIOS
#
# ============================================================================
#
# DESCRIPTION:
#   This script implements a stratified compartmental model for TB transmission
#   and vaccination in the United States. It includes:
#     - Complete model specification with 4 strata and 8 compartments each
#     - Calibration of beta and mtb_prev to stratum-specific incidence
#     - Validation of case source proportions against US TB epidemiology
#     - Multiple intervention scenarios with scaled case counts
#     - Calculation of cases prevented and NNV (number needed to vaccinate)
#
# MODEL STRUCTURE:
#   - 4 mutually exclusive strata: HIV+, Medical conditions, Foreign-born, US-born
#   - 8 compartments per stratum: S, L_f, L_s, L_fv, L_sv, I, I_v, R
#   - 32 total ODEs
#
# CALIBRATION APPROACH (consistent with published TB modeling literature):
#   - Natural history parameters (eps_f, eps_s, gamma, etc.) FIXED from literature
#   - Transmission rate (beta) and LTBI prevalence (mtb_prev) CALIBRATED to incidence
#   - Post-treatment stabilization (alpha_stab): R -> L_s transition limits
#     relapse window to ~5 years, matching clinical data (Millet 2013, Dobler 2009)
#   - Proportion from fast progression reported as VALIDATION metric
#   - Relapse proportion reported as VALIDATION metric (~3-6%, Kim 2013)
#   - Case counts SCALED to 2024 NTSS surveillance data (10,388 cases)
#
#   Key references:
#     Menzies et al. 2018 Lancet ID (systematic review of TB latency modeling)
#     Trauer et al. 2017 Epidemics (optimal latency structures)
#     Menzies et al. 2021 Epidemiology (US TB transmission dynamics)
#     Menzies et al. 2024 Epidemiol Infect (US reactivation TB rates)
#     Ekramnia et al. 2024 (meta-analysis of progression rates)
#
# AUTHOR: [Your name]
# DATE: January 2026
# VERSION: 2.3 (Stratum-specific mixing, demographic balance, post-treatment stabilization)
#
# USAGE:
#   1. Source this file: source("TB_VACCINATION_MODEL_COMPLETE.R")
#   2. Run calibration: cal <- calibrate_model()
#   3. Run scenarios: results <- run_all_scenarios(cal)
#   4. View results: print_scenario_summary(results)
#
# ============================================================================


# ============================================================================
# SECTION 1: SETUP AND DEPENDENCIES
# ============================================================================

# Clear workspace (optional - comment out if you want to preserve objects)
# rm(list = ls())

# Check for deSolve package
# deSolve provides efficient ODE solvers; we fall back to RK4 if unavailable
use_desolve <- requireNamespace("deSolve", quietly = TRUE)
if (use_desolve) {
  library(deSolve)
  cat("Using deSolve package (lsoda solver - recommended)\n")
} else {
  cat("WARNING: deSolve not available - using built-in RK4 solver\n")
  cat("For better performance, install deSolve: install.packages('deSolve')\n")
}
cat("\n")


# ============================================================================
# SECTION 2: POPULATION PARAMETERS
# ============================================================================
# 
# The model stratifies the US population into 4 mutually exclusive groups.
# IMPORTANT: A person can only be in ONE stratum. For example, a foreign-born
# person with HIV is counted in the HIV stratum, NOT the FB stratum.
#
# Sources:
#   - HIV: CDC HIV Surveillance Report 2023
#   - Medical: Literature estimates for diabetes, ESRD, immunosuppression
#   - FB: US Census Bureau, adjusted for overlap with HIV/Medical
#   - USB: Residual US population
# ============================================================================

# Stratum names (for labeling output)
stratum_names <- c("HIV", "Medical", "FB", "USB")
n_strata <- 4

# ACTUAL 2024 TB CASE COUNT (for scaling model output)
# Source: NTSS 2024 provisional data (10,388 reported cases)
# Model equilibrium = 10,340 â†’ scale factor â‰ˆ 1.005 (negligible adjustment)
actual_2024_cases <- 10388

# Population sizes (in persons)
# HIV+:     1.2 million people living with HIV in the US
# Medical:  40.0 million with TB risk factors (diabetes, ESRD, TNF-Î± inhibitors, etc.)
# FB:       40.6 million foreign-born WITHOUT HIV or medical conditions
# USB:      253.4 million US-born WITHOUT HIV or medical conditions
N_vec <- c(
  1.2e6,    # HIV+ (stratum 1)
  40.0e6,   # Medical conditions (stratum 2)
  40.6e6,   # Foreign-born, pure (stratum 3)
  253.4e6   # US-born, general (stratum 4)
)


# ============================================================================
# SECTION 3: DISEASE PROGRESSION PARAMETERS
# ============================================================================
#
# These parameters govern how individuals progress through infection states.
# Key distinction:
#   - Fast (recent) LTBI: High risk in first 2-3 years after infection
#   - Slow (remote) LTBI: Lower risk, can persist for decades
#
# Sources:
#   - eps_s: Ekramnia et al. 2024 (meta-analysis for your dissertation)
#   - eps_f: Adjusted to achieve correct reactivation/recent split
#   - nu_f_trans: Calibrated to produce ~2 year duration in L_f
# ============================================================================

# --- Slow (remote) progression rates [per year] ---
# These are the reactivation rates for people with established LTBI (>2-3 years)
# From Ekramnia et al. 2024 meta-analysis
eps_s <- c(
  0.0069,   # HIV+: 0.69% per year (highest due to immunosuppression)
  0.00104,  # Medical: 0.104% per year
  0.00074,  # FB: 0.074% per year
  0.00068   # USB: 0.068% per year (lowest, reference group)
)

# --- Fast (recent) progression rates [per year] ---
# FIXED from biological argument, consistent with published TB models.
#
# Standard approach in TB modeling literature:
#   - Menzies et al. 2018 (Lancet ID): ~5% of healthy adults develop
#     active TB within 2 years of infection ("fast progressors")
#   - Trauer et al. 2017 (Epidemics): Two-compartment latency structures
#     require L_f -> L_s "stabilisation" transition; first latent phase
#     should be shorter than traditionally assumed
#   - Most published models fix progression rates from literature, then
#     calibrate only beta and/or LTBI prevalence to incidence data
#
# Derivation:
#   Mean time in L_f: T = 1/(eps_f + nu_f_trans + mu)
#   For ~2 year mean duration in immunocompetent strata:
#     eps_f = 1/T - nu_f_trans - mu â‰ˆ 1/2 - 0.30 - 0.01 = 0.19/yr
#   Stratum-specific RR multipliers reflect differential immune status:
#     HIV (RR=2.5): severe immunosuppression -> faster progression
#     Medical (RR=1.25): moderate risk elevation
#     FB, USB (RR=1.0): immunocompetent baseline
#
# These values are NOT calibrated. The proportion of cases from fast
# progression (~8-12%) is reported as a VALIDATION metric.
eps_f_base <- 0.19    # Derived: 1/2yr - 0.30 - ~0.01
eps_f_RR <- c(2.50, 1.25, 1.00, 1.00)  # Relative risk multipliers
eps_f <- eps_f_base * eps_f_RR
# Result: eps_f = c(0.475, 0.238, 0.19, 0.19)

# --- L_f to L_s transition rate [per year] ---
# Individuals move from "fast" to "slow" latent pool over time
# This represents immune system establishing long-term containment
# Rate of 0.30/year gives mean duration in L_f of ~2 years
nu_f_trans <- 0.30

# --- Clearance rates [per year] ---
# Rate at which LTBI spontaneously clears (very slow)
# nu_s = 0.002 gives half-life of ~350 years (effectively lifelong)
# Source: Vynnycky & Fine 1997, Horsburgh 2010
nu_f_clear <- 1.40e-6  # Fast clearance (negligible)
nu_s_clear <- 0.002    # Slow clearance

# --- Fraction of new infections entering fast latent (L_f) ---
# Remainder (1-p) enters slow latent (L_s) directly
# Higher p for HIV+ reflects higher risk of rapid progression
p <- c(
  0.25,   # HIV+: 25% to L_f
  0.20,   # Medical: 20% to L_f
  0.15,   # FB: 15% to L_f
  0.15    # USB: 15% to L_f
)


# ============================================================================
# SECTION 4: TREATMENT AND MORTALITY PARAMETERS
# ============================================================================

# --- Treatment rate [per year] ---
# gamma = 2.0 means average treatment duration of 6 months
# Source: Standard TB treatment duration
gamma <- 2.0

# --- Relapse rate [per year] ---
# Rate at which recovered individuals develop TB again
# Higher for HIV+ due to immunosuppression
# Source: Romanowski 2019, WHO 2023
rho <- c(
  0.015,  # HIV+: 1.5% per year
  0.007,  # Medical: 0.7% per year
  0.005,  # FB: 0.5% per year
  0.005   # USB: 0.5% per year
)

# --- Post-treatment stabilization rate [per year] ---
# Rate at which recovered individuals transition from R -> L_s
# Represents residual bacilli stabilizing into dormancy after treatment
# This is NOT reinfection (no new exposure required) - it reflects the
# biological process by which treated individuals' relapse risk declines
# over time as persister bacilli enter a quiescent state
# Mean time in R = 1/alpha; most relapses occur within 2-5 years
# Source: Barcelona cohort (Millet 2013): 88% of relapses within 3yr
#         Taiwan DOT era (Chen 2023): 1.54% recurrence within 2yr
#         Australia (Dobler 2009): 0.4% recurrence over 5.7yr follow-up
alpha_stab <- 0.20  # Mean ~5 years at elevated relapse risk

# --- Background mortality rate [per year] ---
# All-cause mortality (excluding TB)
# Source: CDC life tables, condition-specific mortality
mu <- c(
  0.012,  # HIV+: 1.2% per year (higher than general population)
  0.025,  # Medical: 2.5% per year (comorbidities increase mortality)
  0.009,  # FB: 0.9% per year
  0.010   # USB: 1.0% per year
)

# --- TB-specific mortality rate [per year] ---
# Additional mortality for individuals with active TB
# Source: NTSS case fatality ratios
mu_TB <- c(
  0.25,   # HIV+: 25% per year (very high)
  0.10,   # Medical: 10% per year
  0.05,   # FB: 5% per year
  0.05    # USB: 5% per year
)


# ============================================================================
# SECTION 5: DEMOGRAPHIC AND TRANSMISSION PARAMETERS
# ============================================================================

# --- Immigration rate [per year] ---
# Only applies to FB stratum
# Represents annual immigration as fraction of FB population
# --- Immigration rate [per year] ---
# Rate at which new individuals enter each stratum from outside the model
# Immigration is balanced by emigration (applied as outflow from all
# compartments at rate iota), maintaining approximately constant population
# Source: Census Bureau immigration statistics
iota <- c(
  0,      # HIV+: no direct immigration to this stratum
  0,      # Medical: no direct immigration
  0.023,  # FB: 2.3% annual immigration rate (balanced by emigration)
  0       # USB: no immigration
)

# --- Reinfection protection ---
# Individuals with prior infection (in R compartment) have partial immunity
# sigma = 0.79 means 79% reduction in susceptibility upon re-exposure
# Source: Andrews et al. 2012
sigma <- 0.79

# --- Mixing parameter ---
# epsilon_mix = within-group mixing preference for each stratum
# Key insight: Assortative mixing only makes sense for SOCIAL/GEOGRAPHIC
# categories (FB vs USB), not CLINICAL categories (HIV, Medical).
# A person with HIV or diabetes contacts the general population, not
# predominantly other PLHIV or diabetics. Their elevated TB risk comes
# from faster progression, not differential exposure.
# Following Menzies et al. (2018, 2021), preferential mixing is applied
# only along the nativity axis (foreign-born vs US-born).
epsilon_mix <- c(
  0.00,   # HIV: proportionate mixing (clinical category, not social)
  0.00,   # Medical: proportionate mixing (clinical category)
  0.70,   # FB: assortative mixing (geographic/social clustering)
  0.70    # USB: assortative mixing (geographic/social clustering)
)


# ============================================================================
# SECTION 6: VACCINATION PARAMETERS
# ============================================================================
#
# These parameters define the vaccination intervention.
# At baseline (no vaccination), psi = 0.
# For intervention scenarios, psi > 0 and theta defines targeting.
# ============================================================================

# --- Vaccine waning rate [per year] ---
# omega = 0.10 means average protection duration of 10 years
omega <- 0.10

# --- Vaccine efficacy against disease progression ---
# VE_d = 0.50 means 50% reduction in progression from L to I
# Based on M72/AS01E trial results (Van Der Meeren 2018)
# Can be varied in scenarios (e.g., 0.30, 0.50, 0.70, 0.90)
VE_d <- 0.50

# --- Baseline vaccination rate [per year] ---
# psi = 0 means no vaccination at baseline
# Intervention scenarios set psi > 0 (e.g., 0.05 = 5% per year)
psi_baseline <- 0

# --- Vaccination targeting by stratum ---
# theta[i] = fraction of stratum i eligible for vaccination
# For pure targeting: theta[i] = 1 (all) or 0 (none)
# For FB-focused targeting: theta reflects FB proportion within each stratum
# Default: target only FB stratum
theta_default <- c(0, 0, 1, 0)  # Only FB stratum targeted

# --- Proportion foreign-born within each stratum ---
# This is CRITICAL for FB-targeted interventions!
# When we want to vaccinate "all foreign-born persons," we need to:
#   - Vaccinate 100% of the FB stratum (theta[3] = 1)
#   - Vaccinate the FB fraction of HIV stratum (theta[1] = prop_fb_hiv)
#   - Vaccinate the FB fraction of Medical stratum (theta[2] = prop_fb_med)
#   - Vaccinate 0% of USB stratum (theta[4] = 0)
#
# Sources:
#   - HIV: Shaji 2020 estimated ~17-19% of PLHIV are foreign-born;
#     CDC surveillance 2010-2017 (Espinoza et al. 2020) reports 18.9% of
#     HIV diagnoses among non-US-born after imputation. We use 19%.
#   - Medical: NHANES 2009-2018 (Han 2022): FB diabetes prevalence 17.1%
#     vs USB 13.2%, yielding ~17.5% FB share of the diabetes population.
#     For the broader Medical stratum (96% diabetes), we use 18%.
prop_fb_in_hiv <- 0.19      # 19% of HIV+ stratum is foreign-born
prop_fb_in_medical <- 0.18  # 18% of Medical stratum is foreign-born

# Theta vector for targeting ALL foreign-born (including FB within HIV/Medical)
theta_all_fb <- c(prop_fb_in_hiv, prop_fb_in_medical, 1, 0)


# ============================================================================
# SECTION 7: CALIBRATION TARGETS
# ============================================================================
#
# These are the incidence rates the model must match.
# Source: National TB Surveillance System (NTSS) 2024
# ============================================================================

# Target incidence [per 100,000 population per year]
target_inc <- c(
  34.83,  # HIV+: 34.83 per 100,000
  9.35,   # Medical: 9.35 per 100,000
  11.67,  # FB: 11.67 per 100,000
  0.57    # USB: 0.57 per 100,000
)

# Overall US incidence (for reference)
target_overall <- 3.08  # per 100,000


# ============================================================================
# SECTION 8: MIXING MATRIX CALCULATION
# ============================================================================
#
# The mixing matrix M[i,j] gives the probability that a contact made by
# someone in stratum i is with someone in stratum j.
#
# With preferential mixing (epsilon_mix > 0):
#   - M[i,i] = epsilon + (1-epsilon) * N[i]/N_total  (within-group)
#   - M[i,j] = (1-epsilon) * N[j]/N_total            (between-group)
#
# Rows sum to 1 (contacts must go somewhere).
# ============================================================================

calc_mixing_matrix <- function(N_vec, epsilon) {
  #' Calculate the mixing matrix for heterogeneous contact patterns
  #'
  #' @param N_vec Vector of population sizes by stratum
  #' @param epsilon Within-group mixing preference: scalar (applied to all)
  #'        or vector of length 4 (stratum-specific)
  #' @return 4x4 mixing matrix M where M[i,j] = P(contact with j | in stratum i)
  
  # Allow scalar or vector epsilon
  if (length(epsilon) == 1) epsilon <- rep(epsilon, 4)
  
  total_N <- sum(N_vec)
  M <- matrix(0, nrow = 4, ncol = 4)
  
  for (i in 1:4) {
    for (j in 1:4) {
      if (i == j) {
        # Within-group contact probability
        M[i, j] <- epsilon[i] + (1 - epsilon[i]) * (N_vec[j] / total_N)
      } else {
        # Between-group contact probability
        M[i, j] <- (1 - epsilon[i]) * (N_vec[j] / total_N)
      }
    }
  }
  
  return(M)
}

# Calculate mixing matrix with default epsilon
M_default <- calc_mixing_matrix(N_vec, epsilon_mix)


# ============================================================================
# SECTION 9: INITIAL CONDITIONS
# ============================================================================
#
# The model needs starting values for all 32 state variables.
# We initialize based on:
#   - mtb_prev: fraction of each stratum with LTBI (goes to L_s)
#   - target_inc: used to estimate initial I and R compartments
#   - Remainder goes to S (susceptible)
# ============================================================================

calc_initial_conditions <- function(N_vec, mtb_prev, target_inc) {
  #' Calculate initial conditions for the ODE system
  #'
  #' @param N_vec Vector of population sizes
  #' @param mtb_prev Vector of LTBI prevalence by stratum (calibrated parameter)
  #' @param target_inc Vector of target incidence rates (per 100,000)
  #' @return Vector of 32 initial values (8 compartments x 4 strata)
  
  y0 <- numeric(32)
  
  for (i in 1:4) {
    # Index for this stratum (compartments are stored sequentially)
    idx <- (i - 1) * 8 + 1
    
    N_i <- N_vec[i]
    
    # --- L_s (slow latent): Initialize from mtb_prev ---
    # All LTBI starts in L_s because immigrants have established infection
    L_f_init <- 0
    L_s_init <- N_i * mtb_prev[i]
    
    # --- I (infectious): Estimate from incidence ---
    # At equilibrium: I â‰ˆ annual_cases / (treatment_rate + death_rate)
    annual_cases <- N_i * (target_inc[i] / 100000)
    I_init <- annual_cases / (gamma + mu[i] + mu_TB[i])
    
    # --- R (recovered): Estimate from treatment flow ---
    # At equilibrium: R â‰ˆ gamma * I / (stabilization + relapse + death)
    R_init <- gamma * I_init / (alpha_stab + rho[i] + mu[i])
    
    # --- S (susceptible): Remainder of population ---
    S_init <- N_i - L_f_init - L_s_init - I_init - R_init
    S_init <- max(S_init, 0)  # Ensure non-negative
    
    # Assign to state vector
    # Order: S, L_f, L_s, L_fv, L_sv, I, I_v, R
    y0[idx]     <- S_init     # S
    y0[idx + 1] <- L_f_init   # L_f
    y0[idx + 2] <- L_s_init   # L_s
    y0[idx + 3] <- 0          # L_fv (vaccinated - starts at 0)
    y0[idx + 4] <- 0          # L_sv (vaccinated - starts at 0)
    y0[idx + 5] <- I_init     # I
    y0[idx + 6] <- 0          # I_v (vaccinated - starts at 0)
    y0[idx + 7] <- R_init     # R
  }
  
  return(y0)
}


# ============================================================================
# SECTION 10: MODEL DIFFERENTIAL EQUATIONS
# ============================================================================
#
# This is the core of the model: the system of 32 ODEs.
#
# For each stratum i, the compartments evolve according to:
#
#   dS/dt   = Recruitment - Infection - Death + Clearance
#   dL_f/dt = Fast infection - Progression - Transition - Vaccination - Death
#   dL_s/dt = Slow infection + Transition - Progression - Vaccination - Death
#   dL_fv/dt = Vaccination of L_f - Progression (reduced) - Transition - Waning - Death
#   dL_sv/dt = Vaccination of L_s + Transition - Progression (reduced) - Waning - Death
#   dI/dt   = Progression from L_f, L_s + Relapse - Treatment - Death
#   dI_v/dt = Progression from L_fv, L_sv - Treatment - Death
#   dR/dt   = Treatment - Reinfection - Relapse - Death
#
# ============================================================================

tb_model <- function(t, y, params) {
  #' TB model differential equations
  #'
  #' @param t Current time (required by ODE solver, not used directly)
  #' @param y Current state vector (32 values)
  #' @param params List of parameters including beta, kappa, psi, theta, mtb_prev, M
  #' @return List containing vector of derivatives (dy/dt)
  
  # --- Unpack parameters ---
  beta <- params$beta           # Transmission rate
  kappa <- params$kappa         # Relative susceptibility by stratum
  psi <- params$psi             # Vaccination rate
  theta <- params$theta         # Vaccination targeting
  mtb_prev <- params$mtb_prev   # LTBI prevalence for recruitment
  M <- params$M                 # Mixing matrix
  VE <- params$VE               # Vaccine efficacy
  
  # --- Extract state variables by stratum ---
  # Each stratum has 8 compartments stored sequentially
  S    <- y[c(1, 9, 17, 25)]    # Susceptible
  L_f  <- y[c(2, 10, 18, 26)]   # Fast latent (unvaccinated)
  L_s  <- y[c(3, 11, 19, 27)]   # Slow latent (unvaccinated)
  L_fv <- y[c(4, 12, 20, 28)]   # Fast latent (vaccinated)
  L_sv <- y[c(5, 13, 21, 29)]   # Slow latent (vaccinated)
  I    <- y[c(6, 14, 22, 30)]   # Infectious (unvaccinated)
  I_v  <- y[c(7, 15, 23, 31)]   # Infectious (vaccinated)
  R    <- y[c(8, 16, 24, 32)]   # Recovered
  
  # --- Calculate total population by stratum ---
  N <- S + L_f + L_s + L_fv + L_sv + I + I_v + R
  
  # --- Calculate force of infection (lambda) ---
  # lambda[i] = rate at which susceptibles in stratum i become infected
  # Depends on: transmission rate, susceptibility, mixing, prevalence of I
  lambda <- numeric(4)
  for (i in 1:4) {
    foi_sum <- 0
    for (j in 1:4) {
      if (N[j] > 0) {
        # Probability of contact with j, times prevalence of infection in j
        foi_sum <- foi_sum + M[i, j] * (I[j] + I_v[j]) / N[j]
      }
    }
    lambda[i] <- beta * kappa[i] * foi_sum
  }
  
  # --- Calculate recruitment (births + immigration) ---
  # Use FIXED population sizes (N_vec from params) for recruitment rates.
  # This ensures constant population: births replace deaths (mu*N_fixed),
  # and immigration adds new entrants (iota*N_fixed) balanced by implicit
  # emigration/aging-out at the same rate, maintaining N â‰ˆ N_fixed.
  N_fixed <- params$N
  Lambda <- mu * N_fixed + iota * N_fixed  # Constant recruitment rate
  
  # Recruits are split between S and L_s based on mtb_prev
  # (immigrants from endemic countries may arrive with LTBI)
  # All LTBI-positive recruits enter L_s (remote latency), consistent with
  # published models and the assumption that immigrant LTBI was acquired
  # years/decades before arrival.
  Gamma_S  <- Lambda * (1 - mtb_prev)  # Recruits to susceptible
  Gamma_Ls <- Lambda * mtb_prev        # Recruits with LTBI (to L_s)
  
  # --- Initialize derivative vectors ---
  dS <- dL_f <- dL_s <- dL_fv <- dL_sv <- dI <- dI_v <- dR <- numeric(4)
  
  # --- Calculate derivatives for each stratum ---
  for (i in 1:4) {
    
    # Combined outflow rate: background death + emigration
    # For strata with immigration (iota > 0), emigration balances inflow
    # to maintain approximately constant population size
    mu_out <- mu[i] + iota[i]
    
    # ----- Susceptible (S) -----
    # Inflow: Recruitment without LTBI, clearance from L_f and L_s
    # Outflow: Infection, death/emigration
    dS[i] <- Gamma_S[i] -                           # Recruitment (susceptible)
             lambda[i] * S[i] -                     # New infections
             mu_out * S[i] +                        # Death + emigration
             nu_f_clear * L_f[i] +                  # Clearance from L_f
             nu_s_clear * L_s[i] +                  # Clearance from L_s
             nu_f_clear * L_fv[i] +                 # Clearance from L_fv
             nu_s_clear * L_sv[i]                   # Clearance from L_sv
    
    # ----- Fast Latent, Unvaccinated (L_f) -----
    # Inflow: Fraction p of new infections from S and R
    # Outflow: Progression to I, transition to L_s, clearance, vaccination, death/emigration
    dL_f[i] <- p[i] * lambda[i] * S[i] +               # New infections (fast, from S)
               p[i] * lambda[i] * R[i] * (1 - sigma) - # Reinfection (fast, from R)
               eps_f[i] * L_f[i] -                     # Progression to active TB
               nu_f_trans * L_f[i] -                   # Transition to slow latent
               nu_f_clear * L_f[i] -                   # Clearance
               psi * theta[i] * L_f[i] +               # Vaccination (outflow)
               omega * L_fv[i] -                       # Waning from vaccinated (inflow)
               mu_out * L_f[i]                         # Death + emigration
    
    # ----- Slow Latent, Unvaccinated (L_s) -----
    # Inflow: Recruitment with LTBI, transition from L_f, fraction (1-p) of new infections,
    #         stabilization from R (post-treatment bacilli enter dormancy)
    # Outflow: Progression to I, clearance, vaccination, death/emigration
    dL_s[i] <- Gamma_Ls[i] +                          # Recruitment with LTBI
               nu_f_trans * L_f[i] +                # Transition from fast latent
               alpha_stab * R[i] +                  # Post-treatment stabilization (R -> L_s)
               (1 - p[i]) * lambda[i] * S[i] +      # New infections (slow, from S)
               (1 - p[i]) * lambda[i] * R[i] * (1 - sigma) - # Reinfection (slow)
               eps_s[i] * L_s[i] -                  # Progression to active TB
               nu_s_clear * L_s[i] -                # Clearance
               psi * theta[i] * L_s[i] +            # Vaccination (outflow)
               omega * L_sv[i] -                    # Waning from vaccinated (inflow)
               mu_out * L_s[i]                      # Death + emigration
    
    # ----- Fast Latent, Vaccinated (L_fv) -----
    # Inflow: Vaccination of L_f
    # Outflow: Progression (reduced by VE), transition to L_sv, clearance, waning, death/emigration
    dL_fv[i] <- psi * theta[i] * L_f[i] -           # Vaccination (inflow)
                (1 - VE) * eps_f[i] * L_fv[i] -     # Progression (reduced)
                nu_f_trans * L_fv[i] -              # Transition to slow
                nu_f_clear * L_fv[i] -              # Clearance
                omega * L_fv[i] -                   # Waning (returns to L_f)
                mu_out * L_fv[i]                    # Death + emigration
    
    # ----- Slow Latent, Vaccinated (L_sv) -----
    # Inflow: Vaccination of L_s, transition from L_fv
    # Outflow: Progression (reduced by VE), clearance, waning, death/emigration
    dL_sv[i] <- psi * theta[i] * L_s[i] +           # Vaccination (inflow)
                nu_f_trans * L_fv[i] -              # Transition from fast
                (1 - VE) * eps_s[i] * L_sv[i] -     # Progression (reduced)
                nu_s_clear * L_sv[i] -              # Clearance
                omega * L_sv[i] -                   # Waning (returns to L_s)
                mu_out * L_sv[i]                    # Death + emigration
    
    # ----- Infectious, Unvaccinated (I) -----
    # Inflow: Progression from L_f and L_s, relapse from R
    # Outflow: Treatment, death (background + TB) + emigration
    dI[i] <- eps_f[i] * L_f[i] +                    # Progression from L_f
             eps_s[i] * L_s[i] +                    # Progression from L_s
             rho[i] * R[i] -                        # Relapse
             gamma * I[i] -                         # Treatment
             (mu_out + mu_TB[i]) * I[i]             # Death + emigration + TB death
    
    # ----- Infectious, Vaccinated (I_v) -----
    # Inflow: Progression from L_fv and L_sv (reduced by VE)
    # Outflow: Treatment, death + emigration
    # Note: No relapse to I_v (relapse goes to I)
    dI_v[i] <- (1 - VE) * eps_f[i] * L_fv[i] +      # Progression from L_fv
               (1 - VE) * eps_s[i] * L_sv[i] -      # Progression from L_sv
               gamma * I_v[i] -                     # Treatment
               (mu_out + mu_TB[i]) * I_v[i]         # Death + emigration + TB death
    
    # ----- Recovered (R) -----
    # Inflow: Treatment of I and I_v
    # Outflow: Stabilization to L_s, reinfection, relapse, death/emigration
    dR[i] <- gamma * I[i] +                         # Treatment (from I)
             gamma * I_v[i] -                       # Treatment (from I_v)
             alpha_stab * R[i] -                    # Stabilization (R -> L_s)
             lambda[i] * R[i] * (1 - sigma) -       # Reinfection
             rho[i] * R[i] -                        # Relapse
             mu_out * R[i]                          # Death + emigration
  }
  
  # --- Pack derivatives into output vector ---
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
  
  # Return as list (required by deSolve)
  return(list(dy))
}


# ============================================================================
# SECTION 11: MODEL RUNNER FUNCTION
# ============================================================================

run_model <- function(params, years = 50, dt = 0.1) {
  #' Run the TB model to equilibrium
  #'
  #' @param params List containing: beta, kappa, psi, theta, mtb_prev, M, VE
  #' @param years Number of years to simulate
  #' @param dt Time step for RK4 solver (ignored if using deSolve)
  #' @return List with trajectory and final_state
  
  # Calculate initial conditions
  y0 <- calc_initial_conditions(N_vec, params$mtb_prev, target_inc)
  
  if (use_desolve) {
    # Use deSolve (preferred - adaptive step size)
    times <- seq(0, years, by = 1)
    out <- ode(y = y0, times = times, func = tb_model, parms = params,
               method = "lsoda")
    final_state <- as.numeric(out[nrow(out), -1])
    trajectory <- out
  } else {
    # Fallback: 4th-order Runge-Kutta
    n_steps <- as.integer(years / dt) + 1
    y <- y0
    
    for (step in 2:n_steps) {
      t <- (step - 1) * dt
      k1 <- unlist(tb_model(t, y, params))
      k2 <- unlist(tb_model(t + dt/2, y + dt/2 * k1, params))
      k3 <- unlist(tb_model(t + dt/2, y + dt/2 * k2, params))
      k4 <- unlist(tb_model(t + dt, y + dt * k3, params))
      y <- y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
      y <- pmax(y, 0)  # Ensure non-negative
    }
    
    final_state <- y
    trajectory <- NULL
  }
  
  return(list(
    trajectory = trajectory,
    final_state = final_state,
    params = params
  ))
}


# ============================================================================
# SECTION 12: OUTPUT CALCULATION FUNCTIONS
# ============================================================================

extract_compartments <- function(state) {
  #' Extract compartment vectors from state vector
  #'
  #' @param state Vector of 32 state values
  #' @return List with vectors for each compartment type
  
  list(
    S    = state[c(1, 9, 17, 25)],
    L_f  = state[c(2, 10, 18, 26)],
    L_s  = state[c(3, 11, 19, 27)],
    L_fv = state[c(4, 12, 20, 28)],
    L_sv = state[c(5, 13, 21, 29)],
    I    = state[c(6, 14, 22, 30)],
    I_v  = state[c(7, 15, 23, 31)],
    R    = state[c(8, 16, 24, 32)]
  )
}


calc_incidence <- function(state) {
  #' Calculate TB incidence from model state
  #'
  #' @param state Vector of 32 state values
  #' @return List with stratum-specific and overall incidence (per 100,000)
  
  comp <- extract_compartments(state)
  
  # Total population by stratum
  N <- comp$S + comp$L_f + comp$L_s + comp$L_fv + comp$L_sv + 
       comp$I + comp$I_v + comp$R
  
  # Incidence = prevalent cases / population * 100,000
  # (This is point prevalence used as proxy for incidence at equilibrium)
  inc_stratum <- ((comp$I + comp$I_v) / N) * 100000
  
  # Overall incidence
  total_I <- sum(comp$I) + sum(comp$I_v)
  total_N <- sum(N)
  inc_overall <- (total_I / total_N) * 100000
  
  return(list(
    stratum = inc_stratum,
    overall = inc_overall,
    names = stratum_names
  ))
}


calc_annual_cases <- function(state, params) {
  #' Calculate annual TB cases
  #'
  #' @param state Vector of 32 state values
  #' @param params Parameter list (needs VE)
  #' @return List with case counts by source and stratum
  #'
  #' Annual cases are calculated from the model's incidence rate (I/N * 100,000)
  #' multiplied by population. This ensures consistency with how we calibrate
  #' and report incidence.
  
  comp <- extract_compartments(state)
  VE <- params$VE
  
  # ----- POPULATION BY STRATUM -----
  N <- comp$S + comp$L_f + comp$L_s + comp$L_fv + comp$L_sv + 
       comp$I + comp$I_v + comp$R
  
  # ----- INCIDENCE RATE (per 100,000) -----
  inc_rate <- ((comp$I + comp$I_v) / N) * 100000
  
  # ----- ANNUAL CASES FROM INCIDENCE -----
  # Annual cases = incidence_rate * N / 100,000
  actual_cases_by_stratum <- inc_rate * N / 100000
  actual_total <- sum(actual_cases_by_stratum)
  
  # ----- SOURCE BREAKDOWN -----
  # Calculate the FLOW into I compartment from each source
  
  # Flow from fast latent (recent infection)
  flow_Lf <- eps_f * comp$L_f + (1 - VE) * eps_f * comp$L_fv
  
  # Flow from slow latent (reactivation)
  flow_Ls <- eps_s * comp$L_s + (1 - VE) * eps_s * comp$L_sv
  
  # Flow from relapse
  flow_relapse <- rho * comp$R
  
  # Total flow
  total_flow <- sum(flow_Lf) + sum(flow_Ls) + sum(flow_relapse)
  
  # Calculate proportions
  pct_Lf <- sum(flow_Lf) / total_flow * 100
  pct_Ls <- sum(flow_Ls) / total_flow * 100
  pct_relapse <- sum(flow_relapse) / total_flow * 100
  
  # Scale source-specific cases to match actual total
  from_Lf <- actual_total * (sum(flow_Lf) / total_flow)
  from_Ls <- actual_total * (sum(flow_Ls) / total_flow)
  from_relapse <- actual_total * (sum(flow_relapse) / total_flow)
  
  return(list(
    by_stratum = actual_cases_by_stratum,
    from_Lf = from_Lf,
    from_Ls = from_Ls,
    from_relapse = from_relapse,
    total = actual_total,
    pct_Lf = pct_Lf,
    pct_Ls = pct_Ls,
    pct_relapse = pct_relapse,
    pct_reactivation = pct_Ls + pct_relapse
  ))
}


calc_vaccinations <- function(state, psi, theta) {
  #' Calculate annual vaccinations
  #'
  #' @param state Vector of 32 state values
  #' @param psi Vaccination rate
  #' @param theta Targeting vector
  #' @return Annual vaccinations by stratum
  

  comp <- extract_compartments(state)
  
  # Vaccinations = psi * theta * (L_f + L_s) for each stratum
  annual_vax <- psi * theta * (comp$L_f + comp$L_s)
  
  return(list(
    by_stratum = annual_vax,
    total = sum(annual_vax)
  ))
}


# ============================================================================
# SECTION 13: CALIBRATION FUNCTION
# ============================================================================
#
# Calibration finds the values of beta and mtb_prev that minimize the
# difference between model-predicted and observed incidence.
#
# Method: L-BFGS-B optimization (box-constrained)
# Objective: Sum of squared relative errors
# ============================================================================

calibrate_model <- function(initial_beta = 5.0,
                            initial_mtb_prev = c(0.04, 0.11, 0.21, 0.008),
                            verbose = TRUE) {
  #' Calibrate the model to match target incidence
  #'
  #' Calibrates beta and mtb_prev[1:4] to match stratum-specific incidence.
  #' This follows the standard approach in TB modeling:
  #'   - Fix natural history parameters from literature (eps_f, eps_s, etc.)
  #'   - Calibrate transmission rate and LTBI prevalence to incidence data
  #'   - Report proportion from fast progression as validation metric
  #'
  #' References:
  #'   - Menzies et al. 2018 Lancet ID (systematic review of TB models)
  #'   - Menzies et al. 2021 Epidemiology (US TB transmission model)
  #'   - Hill et al. 2021 PLoS ONE (US LTBI prevalence estimation)
  #'
  #' @param initial_beta Starting value for transmission rate
  #' @param initial_mtb_prev Starting values for LTBI prevalence

  #' @param verbose Print progress messages
  #' @return List with calibrated parameters and fit statistics
  
  if (verbose) {
    cat("=======================================================================\n")
    cat("MODEL CALIBRATION\n")
    cat("=======================================================================\n\n")
    cat("Free parameters: beta, mtb_prev[1:4] (5 total)\n")
    cat("Targets: Stratum-specific incidence (4 targets)\n")
    cat("Fixed: eps_f (from literature), eps_s (Ekramnia 2024)\n")
    cat("Note: More parameters than targets => exact fit is guaranteed.\n")
    cat("      Model validity assessed via VALIDATION metrics (case sources).\n\n")
  }
  
  # --- Objective function ---
  # Parameters are on log scale to ensure positivity
  objective <- function(theta) {
    beta <- exp(theta[1])
    mtb_prev <- exp(theta[2:5])
    
    # Bounds check
    if (beta < 0.1 || beta > 50) return(1e10)
    if (any(mtb_prev < 0.0001) || any(mtb_prev > 0.6)) return(1e10)
    
    # Create parameter list
    params <- list(
      beta = beta,
      kappa = c(1, 1, 1, 1),
      psi = 0,
      theta = c(0, 0, 0, 0),
      mtb_prev = mtb_prev,
      M = M_default, N = N_vec,
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
  theta0 <- log(c(initial_beta, initial_mtb_prev))
  
  if (verbose) cat("Running optimization...\n")
  
  result <- optim(
    par = theta0,
    fn = objective,
    method = "L-BFGS-B",
    lower = log(c(0.1, 0.0001, 0.0001, 0.0001, 0.0001)),
    upper = log(c(50, 0.6, 0.6, 0.6, 0.6)),
    control = list(maxit = 500)
  )
  
  # --- Extract calibrated parameters ---
  beta_cal <- exp(result$par[1])
  mtb_prev_cal <- exp(result$par[2:5])
  
  # --- Validate calibration ---
  params_cal <- list(
    beta = beta_cal,
    kappa = c(1, 1, 1, 1),
    psi = 0,
    theta = c(0, 0, 0, 0),
    mtb_prev = mtb_prev_cal,
    M = M_default, N = N_vec,
    VE = VE_d
  )
  
  final_result <- run_model(params_cal, years = 50)
  final_inc <- calc_incidence(final_result$final_state)
  final_cases <- calc_annual_cases(final_result$final_state, params_cal)
  
  # Calculate fit statistics
  abs_errors <- abs(final_inc$stratum - target_inc)
  rel_errors <- abs_errors / target_inc * 100
  mape <- mean(rel_errors)
  
  # R-squared
  ss_res <- sum((final_inc$stratum - target_inc)^2)
  ss_tot <- sum((target_inc - mean(target_inc))^2)
  r_squared <- 1 - ss_res / ss_tot
  
  if (verbose) {
    cat(sprintf("\nConvergence: %s\n", 
                ifelse(result$convergence == 0, "Successful", "Check results")))
    cat(sprintf("Function evaluations: %d\n\n", result$counts[1]))
    
    cat("CALIBRATED PARAMETERS:\n")
    cat(sprintf("  beta = %.4f\n", beta_cal))
    cat(sprintf("  mtb_prev = [%.4f, %.4f, %.4f, %.4f]\n",
                mtb_prev_cal[1], mtb_prev_cal[2], mtb_prev_cal[3], mtb_prev_cal[4]))
    cat(sprintf("           = [%.2f%%, %.2f%%, %.2f%%, %.3f%%]\n\n",
                mtb_prev_cal[1]*100, mtb_prev_cal[2]*100, 
                mtb_prev_cal[3]*100, mtb_prev_cal[4]*100))
    
    cat("INCIDENCE TARGETS (matched by construction):\n")
    cat("  Note: With 5 free parameters and 4 targets, exact fit is guaranteed.\n")
    cat("  Model validity depends on VALIDATION metrics below, not fit statistics.\n")
    cat(sprintf("  %-10s  %8s  %8s\n", "Stratum", "Model", "Target"))
    cat("  ------------------------------\n")
    for (i in 1:4) {
      cat(sprintf("  %-10s  %8.2f  %8.2f\n",
                  stratum_names[i], final_inc$stratum[i], target_inc[i]))
    }
    
    cat("\nVALIDATION METRICS (not calibrated to):\n")
    cat(sprintf("  From L_f (recent infection): %.1f%%\n", final_cases$pct_Lf))
    cat(sprintf("  From L_s (remote reactivation): %.1f%%\n", final_cases$pct_Ls))
    cat(sprintf("  From R (relapse): %.1f%%\n", final_cases$pct_relapse))
    cat(sprintf("  Combined reactivation (L_s + R): %.1f%%\n", final_cases$pct_reactivation))
    cat("  Expected: >80%% reactivation in US (Menzies 2024)\n")
    cat("  Expected: relapse â‰¤6%% of cases (Kim 2013, CDC NTSS)\n")
    if (final_cases$pct_reactivation >= 80 && final_cases$pct_reactivation <= 95) {
      cat("  >> VALIDATION PASSED: Reactivation proportion in expected range\n\n")
    } else {
      cat("  >> NOTE: Reactivation proportion outside typical range\n\n")
    }
  }
  
  return(list(
    beta = beta_cal,
    mtb_prev = mtb_prev_cal,
    params = params_cal,
    mape = mape,
    r_squared = r_squared,
    model_inc = final_inc$stratum,
    pct_fast = final_cases$pct_Lf,
    pct_relapse = final_cases$pct_relapse,
    pct_reactivation = final_cases$pct_reactivation,
    baseline_cases = final_cases$total,
    convergence = result$convergence
  ))
}


# ============================================================================
# SECTION 14: INTERVENTION SCENARIO DEFINITIONS
# ============================================================================
#
# Each scenario is defined by:
#   - psi: Vaccination rate (fraction of eligible LTBI vaccinated per year)
#   - theta: Targeting vector (which strata are eligible)
#   - VE: Vaccine efficacy
#   - description: Human-readable description
# ============================================================================

define_scenarios <- function() {
  #' Define intervention scenarios
  #'
  #' @return List of scenario definitions
  #'
  #' NOTE ON THETA VALUES:
  #' theta[i] represents the proportion of stratum i eligible for vaccination.
  #' 
  #' For "FB stratum only" scenarios: theta = c(0, 0, 1, 0)
  #'   - Only vaccinates the pure FB stratum (FB without HIV/Medical)
  #'
  #' For "All foreign-born" scenarios: theta = c(0.19, 0.18, 1, 0)
  #'   - Vaccinates 100% of FB stratum
  #'   - Vaccinates 19% of HIV stratum (the FB within HIV)
  #'   - Vaccinates 18% of Medical stratum (the FB within Medical)
  #'   - This captures ALL foreign-born persons across all strata
  
  scenarios <- list(
    
    # =========================================================================
    # BASELINE (no vaccination)
    # =========================================================================
    baseline = list(
      name = "Baseline (no vaccination)",
      psi = 0,
      theta = c(0, 0, 0, 0),
      VE = 0.50,
      description = "No vaccination intervention"
    ),
    
    # =========================================================================
    # FB STRATUM ONLY - Targets only the "pure" FB stratum
    # (Does NOT include FB persons in HIV or Medical strata)
    # =========================================================================
    fb_stratum_50ve = list(
      name = "FB stratum only, 50% VE",
      psi = 0.05,
      theta = c(0, 0, 1, 0),
      VE = 0.50,
      description = "5%/year vaccination of FB stratum only (excludes FB in HIV/Medical)"
    ),
    
    fb_stratum_70ve = list(
      name = "FB stratum only, 70% VE",
      psi = 0.05,
      theta = c(0, 0, 1, 0),
      VE = 0.70,
      description = "5%/year vaccination of FB stratum only (excludes FB in HIV/Medical)"
    ),
    
    fb_stratum_90ve = list(
      name = "FB stratum only, 90% VE",
      psi = 0.05,
      theta = c(0, 0, 1, 0),
      VE = 0.90,
      description = "5%/year vaccination of FB stratum only (excludes FB in HIV/Medical)"
    ),
    
    # =========================================================================
    # ALL FOREIGN-BORN - Targets FB persons in ALL strata
    # theta = c(prop_fb_hiv, prop_fb_med, 1, 0) = c(0.19, 0.18, 1, 0)
    # This captures the FB within HIV and Medical strata too
    # =========================================================================
    all_fb_50ve = list(
      name = "All FB (incl. HIV/Med), 50% VE",
      psi = 0.05,
      theta = c(0.19, 0.18, 1, 0),  # FB proportion in each stratum
      VE = 0.50,
      description = "5%/year vaccination of ALL foreign-born (including FB in HIV/Medical)"
    ),
    
    all_fb_70ve = list(
      name = "All FB (incl. HIV/Med), 70% VE",
      psi = 0.05,
      theta = c(0.19, 0.18, 1, 0),  # FB proportion in each stratum
      VE = 0.70,
      description = "5%/year vaccination of ALL foreign-born (including FB in HIV/Medical)"
    ),
    
    all_fb_90ve = list(
      name = "All FB (incl. HIV/Med), 90% VE",
      psi = 0.05,
      theta = c(0.19, 0.18, 1, 0),  # FB proportion in each stratum
      VE = 0.90,
      description = "5%/year vaccination of ALL foreign-born (including FB in HIV/Medical)"
    ),
    
    all_fb_intensive = list(
      name = "All FB intensive (10%/yr), 70% VE",
      psi = 0.10,
      theta = c(0.19, 0.18, 1, 0),  # FB proportion in each stratum
      VE = 0.70,
      description = "10%/year vaccination of ALL foreign-born, 70% vaccine efficacy"
    ),
    
    # =========================================================================
    # HIV STRATUM TARGETED (all HIV+, regardless of nativity)
    # =========================================================================
    hiv_all = list(
      name = "All HIV+",
      psi = 0.10,
      theta = c(1, 0, 0, 0),
      VE = 0.70,
      description = "10%/year vaccination of all HIV+ persons, 70% VE"
    ),
    
    # =========================================================================
    # MEDICAL STRATUM TARGETED (all with medical conditions)
    # =========================================================================
    med_all = list(
      name = "All Medical",
      psi = 0.05,
      theta = c(0, 1, 0, 0),
      VE = 0.70,
      description = "5%/year vaccination of all with medical conditions, 70% VE"
    ),
    
    # =========================================================================
    # HIGH-RISK COMBINATIONS
    # =========================================================================
    hiv_and_med = list(
      name = "HIV + Medical",
      psi = 0.05,
      theta = c(1, 1, 0, 0),
      VE = 0.50,
      description = "5%/year vaccination of HIV and Medical strata, 50% VE"
    ),
    
    # All high-risk: HIV + Medical + FB (but NOT FB within HIV/Med, that's covered)
    all_high_risk = list(
      name = "All high-risk (HIV + Med + FB)",
      psi = 0.05,
      theta = c(1, 1, 1, 0),
      VE = 0.70,
      description = "5%/year vaccination of HIV, Medical, and FB strata, 70% VE"
    ),
    
    # =========================================================================
    # UNIVERSAL VACCINATION
    # =========================================================================
    universal_50ve = list(
      name = "Universal, 50% VE",
      psi = 0.02,
      theta = c(1, 1, 1, 1),
      VE = 0.50,
      description = "2%/year vaccination of all strata, 50% vaccine efficacy"
    ),
    
    universal_70ve = list(
      name = "Universal, 70% VE",
      psi = 0.02,
      theta = c(1, 1, 1, 1),
      VE = 0.70,
      description = "2%/year vaccination of all strata, 70% vaccine efficacy"
    )
  )
  
  return(scenarios)
}


# ============================================================================
# SECTION 15: RUN INTERVENTION SCENARIOS
# ============================================================================

run_scenario <- function(scenario, calibration, years = 30) {
  #' Run a single intervention scenario
  #'
  #' @param scenario Scenario definition (from define_scenarios)
  #' @param calibration Calibration results (from calibrate_model)
  #' @param years Simulation horizon
  #' @return List with scenario results
  
  # Build parameter list
  params <- list(
    beta = calibration$beta,
    kappa = c(1, 1, 1, 1),
    psi = scenario$psi,
    theta = scenario$theta,
    mtb_prev = calibration$mtb_prev,
    M = M_default, N = N_vec,
    VE = scenario$VE
  )
  
  # Run model
  result <- run_model(params, years = years)
  
  # Calculate outcomes
  inc <- calc_incidence(result$final_state)
  cases <- calc_annual_cases(result$final_state, params)
  vax <- calc_vaccinations(result$final_state, scenario$psi, scenario$theta)
  
  # Calculate cumulative values over time horizon
  # (Simplified: assumes equilibrium reached, multiply annual by years)
  cumulative_cases <- cases$total * years
  cumulative_vax <- vax$total * years
  
  return(list(
    name = scenario$name,
    description = scenario$description,
    psi = scenario$psi,
    theta = scenario$theta,
    VE = scenario$VE,
    inc_stratum = inc$stratum,
    inc_overall = inc$overall,
    annual_cases = cases$total,
    cases_by_stratum = cases$by_stratum,
    annual_vaccinations = vax$total,
    vax_by_stratum = vax$by_stratum,
    cumulative_cases = cumulative_cases,
    cumulative_vaccinations = cumulative_vax,
    pct_reactivation = cases$pct_reactivation
  ))
}


run_all_scenarios <- function(calibration, years = 30, verbose = TRUE) {
  #' Run all intervention scenarios
  #'
  #' @param calibration Calibration results
  #' @param years Simulation horizon
  #' @param verbose Print progress
  #' @return List with all scenario results and comparison metrics
  
  scenarios <- define_scenarios()
  results <- list()
  
  if (verbose) {
    cat("=======================================================================\n")
    cat("RUNNING INTERVENTION SCENARIOS\n")
    cat("=======================================================================\n\n")
    cat(sprintf("Time horizon: %d years\n", years))
    cat(sprintf("Number of scenarios: %d\n\n", length(scenarios)))
  }
  
  # Run baseline first
  if (verbose) cat("Running: Baseline... ")
  baseline <- run_scenario(scenarios$baseline, calibration, years)
  results$baseline <- baseline
  if (verbose) cat("done\n")
  
  # Run intervention scenarios
  for (name in names(scenarios)) {
    if (name == "baseline") next
    
    if (verbose) cat(sprintf("Running: %s... ", scenarios[[name]]$name))
    results[[name]] <- run_scenario(scenarios[[name]], calibration, years)
    if (verbose) cat("done\n")
  }
  
  # Calculate comparison metrics
  if (verbose) cat("\nCalculating comparison metrics...\n")
  
  # Scaling factor to convert model cases to actual 2024 cases
  baseline_cases <- baseline$annual_cases
  scale_factor <- actual_2024_cases / baseline_cases
  
  comparison <- data.frame(
    scenario = character(),
    annual_cases_model = numeric(),
    annual_cases_scaled = numeric(),
    cases_prevented = numeric(),
    pct_reduction = numeric(),
    cumulative_vax = numeric(),
    NNV = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(results)) {
    r <- results[[name]]
    
    cases_prevented_model <- baseline_cases - r$annual_cases
    cases_prevented_scaled <- cases_prevented_model * scale_factor
    pct_reduction <- cases_prevented_model / baseline_cases * 100
    
    # NNV = cumulative vaccinations / cumulative cases prevented
    cumulative_prevented <- cases_prevented_scaled * years
    if (cumulative_prevented > 0) {
      nnv <- r$cumulative_vaccinations / cumulative_prevented
    } else {
      nnv <- NA
    }
    
    comparison <- rbind(comparison, data.frame(
      scenario = r$name,
      annual_cases_model = round(r$annual_cases),
      annual_cases_scaled = round(r$annual_cases * scale_factor),
      cases_prevented = round(cases_prevented_scaled),
      pct_reduction = round(pct_reduction, 1),
      cumulative_vax = round(r$cumulative_vaccinations / 1e6, 2),
      NNV = round(nnv),
      stringsAsFactors = FALSE
    ))
  }
  
  results$comparison <- comparison
  results$baseline_cases <- baseline_cases
  results$years <- years
  
  if (verbose) cat("done\n")
  
  return(results)
}


# ============================================================================
# SECTION 16: OUTPUT AND REPORTING FUNCTIONS
# ============================================================================

print_scenario_summary <- function(results) {
  #' Print summary table of scenario results
  #'
  #' @param results Output from run_all_scenarios
  
  cat("\n")
  cat("=======================================================================\n")
  cat("INTERVENTION SCENARIO RESULTS\n")
  cat("=======================================================================\n\n")
  
  # Calculate scale factor
  scale_factor <- actual_2024_cases / results$baseline_cases
  
  cat(sprintf("Baseline (model equilibrium): %s cases/year\n", 
              format(round(results$baseline_cases), big.mark = ",")))
  cat(sprintf("Baseline (scaled to 2024):    %s cases/year\n",
              format(actual_2024_cases, big.mark = ",")))
  cat(sprintf("Time horizon: %d years\n\n", results$years))
  
  comp <- results$comparison
  
  cat(sprintf("%-30s  %10s  %10s  %8s  %10s  %6s\n",
              "Scenario", "Cases/yr", "Prevented", "% Red.", "Cum. Vax", "NNV"))
  cat(paste(rep("-", 90), collapse = ""), "\n")
  
  for (i in 1:nrow(comp)) {
    r <- comp[i,]
    nnv_str <- ifelse(is.na(r$NNV), "â€”", as.character(r$NNV))
    cat(sprintf("%-30s  %10s  %10s  %7.1f%%  %9sM  %6s\n",
                r$scenario,
                format(r$annual_cases_scaled, big.mark = ","),
                format(r$cases_prevented, big.mark = ","),
                r$pct_reduction,
                r$cumulative_vax,
                nnv_str))
  }
  
  cat("\n")
  cat("Cases/yr and Prevented are scaled to match 2024 US TB case count\n")
  cat("NNV = Number Needed to Vaccinate to prevent one TB case\n")
  cat("Cum. Vax = Cumulative vaccinations over time horizon (millions)\n")
}


print_stratum_impact <- function(results, scenario_name) {
  #' Print stratum-level impact for a specific scenario
  #'
  #' @param results Output from run_all_scenarios
  #' @param scenario_name Name of scenario to display
  
  if (!(scenario_name %in% names(results))) {
    cat(sprintf("Scenario '%s' not found\n", scenario_name))
    return(invisible(NULL))
  }
  
  baseline <- results$baseline
  scenario <- results[[scenario_name]]
  
  cat(sprintf("\nSTRATUM-LEVEL IMPACT: %s\n", scenario$name))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat(sprintf("%-12s  %12s  %12s  %12s  %8s\n",
              "Stratum", "Baseline", "Intervention", "Prevented", "% Red."))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  for (i in 1:4) {
    base_cases <- baseline$cases_by_stratum[i]
    int_cases <- scenario$cases_by_stratum[i]
    prevented <- base_cases - int_cases
    pct_red <- prevented / base_cases * 100
    
    cat(sprintf("%-12s  %12.0f  %12.0f  %12.0f  %7.1f%%\n",
                stratum_names[i], base_cases, int_cases, prevented, pct_red))
  }
  
  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat(sprintf("%-12s  %12.0f  %12.0f  %12.0f  %7.1f%%\n",
              "TOTAL",
              baseline$annual_cases,
              scenario$annual_cases,
              baseline$annual_cases - scenario$annual_cases,
              (baseline$annual_cases - scenario$annual_cases) / baseline$annual_cases * 100))
}


# ============================================================================
# SECTION 17: SENSITIVITY ANALYSIS FUNCTIONS
# ============================================================================

run_ve_sensitivity <- function(calibration, psi = 0.05, theta = c(0, 0, 1, 0),
                                ve_values = c(0.30, 0.50, 0.70, 0.90),
                                years = 30) {
  #' Run sensitivity analysis on vaccine efficacy
  #'
  #' @param calibration Calibration results
  #' @param psi Vaccination rate
  #' @param theta Targeting vector
  #' @param ve_values Vector of VE values to test
  #' @param years Simulation horizon
  #' @return Data frame with results
  
  # Get baseline
  params_base <- list(
    beta = calibration$beta, kappa = c(1,1,1,1), psi = 0,
    theta = c(0,0,0,0), mtb_prev = calibration$mtb_prev,
    M = M_default, N = N_vec, VE = 0.50
  )
  baseline <- run_model(params_base, years = years)
  baseline_cases <- calc_annual_cases(baseline$final_state, params_base)$total
  
  results <- data.frame(
    VE = numeric(),
    annual_cases = numeric(),
    cases_prevented = numeric(),
    pct_reduction = numeric(),
    NNV = numeric()
  )
  
  for (ve in ve_values) {
    params <- list(
      beta = calibration$beta, kappa = c(1,1,1,1), psi = psi,
      theta = theta, mtb_prev = calibration$mtb_prev,
      M = M_default, N = N_vec, VE = ve
    )
    
    result <- run_model(params, years = years)
    cases <- calc_annual_cases(result$final_state, params)$total
    vax <- calc_vaccinations(result$final_state, psi, theta)$total
    
    prevented <- baseline_cases - cases
    cumulative_prevented <- prevented * years
    cumulative_vax <- vax * years
    nnv <- ifelse(cumulative_prevented > 0, cumulative_vax / cumulative_prevented, NA)
    
    results <- rbind(results, data.frame(
      VE = ve * 100,
      annual_cases = round(cases),
      cases_prevented = round(prevented),
      pct_reduction = round(prevented / baseline_cases * 100, 1),
      NNV = round(nnv)
    ))
  }
  
  return(results)
}


run_psi_sensitivity <- function(calibration, VE = 0.70, theta = c(0, 0, 1, 0),
                                 psi_values = c(0.01, 0.02, 0.05, 0.10, 0.15),
                                 years = 30) {
  #' Run sensitivity analysis on vaccination rate
  #'
  #' @param calibration Calibration results
  #' @param VE Vaccine efficacy
  #' @param theta Targeting vector
  #' @param psi_values Vector of psi values to test
  #' @param years Simulation horizon
  #' @return Data frame with results
  
  # Get baseline
  params_base <- list(
    beta = calibration$beta, kappa = c(1,1,1,1), psi = 0,
    theta = c(0,0,0,0), mtb_prev = calibration$mtb_prev,
    M = M_default, N = N_vec, VE = VE
  )
  baseline <- run_model(params_base, years = years)
  baseline_cases <- calc_annual_cases(baseline$final_state, params_base)$total
  
  results <- data.frame(
    psi_pct = numeric(),
    annual_cases = numeric(),
    cases_prevented = numeric(),
    pct_reduction = numeric(),
    annual_vax = numeric(),
    NNV = numeric()
  )
  
  for (psi in psi_values) {
    params <- list(
      beta = calibration$beta, kappa = c(1,1,1,1), psi = psi,
      theta = theta, mtb_prev = calibration$mtb_prev,
      M = M_default, N = N_vec, VE = VE
    )
    
    result <- run_model(params, years = years)
    cases <- calc_annual_cases(result$final_state, params)$total
    vax <- calc_vaccinations(result$final_state, psi, theta)$total
    
    prevented <- baseline_cases - cases
    cumulative_prevented <- prevented * years
    cumulative_vax <- vax * years
    nnv <- ifelse(cumulative_prevented > 0, cumulative_vax / cumulative_prevented, NA)
    
    results <- rbind(results, data.frame(
      psi_pct = psi * 100,
      annual_cases = round(cases),
      cases_prevented = round(prevented),
      pct_reduction = round(prevented / baseline_cases * 100, 1),
      annual_vax = round(vax),
      NNV = round(nnv)
    ))
  }
  
  return(results)
}


# ============================================================================
# SECTION 18: MAIN EXECUTION
# ============================================================================
#
# This section runs when the script is executed directly.
# It performs calibration, runs scenarios, and displays results.
# ============================================================================

cat("\n")
cat("########################################################################\n")
cat("#                                                                      #\n")
cat("#           TB VACCINATION MODEL - COMPLETE ANALYSIS                   #\n")
cat("#                                                                      #\n")
cat("########################################################################\n\n")

# --- Step 1: Calibrate the model ---
cat("STEP 1: MODEL CALIBRATION\n")
cat("=======================================================================\n\n")

calibration <- calibrate_model(verbose = TRUE)

# --- Step 2: Run intervention scenarios ---
cat("\nSTEP 2: INTERVENTION SCENARIOS\n")

scenario_results <- run_all_scenarios(calibration, years = 30, verbose = TRUE)

# --- Step 3: Display results ---
cat("\nSTEP 3: RESULTS\n")

print_scenario_summary(scenario_results)

# --- Step 4: Show stratum-level impact for key scenarios ---
print_stratum_impact(scenario_results, "fb_stratum_70ve")
print_stratum_impact(scenario_results, "all_fb_70ve")

# --- Step 5: Sensitivity analyses ---
cat("\n=======================================================================\n")
cat("SENSITIVITY ANALYSIS: VACCINE EFFICACY\n")
cat("=======================================================================\n\n")

ve_sens <- run_ve_sensitivity(calibration)
cat("FB-targeted vaccination (5%/year), varying VE:\n\n")
cat(sprintf("%-8s  %12s  %12s  %8s  %6s\n", 
            "VE", "Annual Cases", "Prevented", "% Red.", "NNV"))
cat(paste(rep("-", 55), collapse = ""), "\n")
for (i in 1:nrow(ve_sens)) {
  cat(sprintf("%-8s  %12s  %12s  %7.1f%%  %6d\n",
              paste0(ve_sens$VE[i], "%"),
              format(ve_sens$annual_cases[i], big.mark = ","),
              format(ve_sens$cases_prevented[i], big.mark = ","),
              ve_sens$pct_reduction[i],
              ve_sens$NNV[i]))
}

cat("\n=======================================================================\n")
cat("SENSITIVITY ANALYSIS: VACCINATION RATE\n")
cat("=======================================================================\n\n")

psi_sens <- run_psi_sensitivity(calibration)
cat("FB-targeted vaccination (70% VE), varying Ïˆ:\n\n")
cat(sprintf("%-10s  %12s  %12s  %8s  %12s  %6s\n",
            "Ïˆ (%/yr)", "Annual Cases", "Prevented", "% Red.", "Annual Vax", "NNV"))
cat(paste(rep("-", 75), collapse = ""), "\n")
for (i in 1:nrow(psi_sens)) {
  cat(sprintf("%-10s  %12s  %12s  %7.1f%%  %12s  %6d\n",
              paste0(psi_sens$psi_pct[i], "%"),
              format(psi_sens$annual_cases[i], big.mark = ","),
              format(psi_sens$cases_prevented[i], big.mark = ","),
              psi_sens$pct_reduction[i],
              format(psi_sens$annual_vax[i], big.mark = ","),
              psi_sens$NNV[i]))
}

# --- Final summary ---
cat("\n=======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=======================================================================\n\n")

cat("Key objects created:\n")
cat("  - calibration: Calibrated model parameters\n")
cat("  - scenario_results: All intervention scenario outcomes\n")
cat("  - ve_sens: Vaccine efficacy sensitivity results\n")
cat("  - psi_sens: Vaccination rate sensitivity results\n\n")

cat("To explore further:\n")
cat("  - print_stratum_impact(scenario_results, 'scenario_name')\n")
cat("  - run_ve_sensitivity(calibration, psi = X, theta = c(...))\n")
cat("  - run_psi_sensitivity(calibration, VE = X, theta = c(...))\n")
