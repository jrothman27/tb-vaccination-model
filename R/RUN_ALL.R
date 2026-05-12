################################################################################
#
#  RUN_ALL.R — TB VACCINATION MODEL: COMPLETE ANALYSIS PIPELINE
#
################################################################################
#
# Runs every analysis in the modular pipeline in the correct order. Uses the
# unified scripts (which run both plausible and optimistic scenarios in one
# pass) where they exist.
#
# REQUIREMENTS:
#   - All R scripts in the working directory
#   - Packages: deSolve, lhs, ggplot2, scales, patchwork, RColorBrewer,
#     viridis, dplyr, tidyr, ggtext
#
# ESTIMATED RUNTIME:
#   - ~60-180 minutes total
#   - LHS analysis dominates (45-90 min for both scenarios at 500 draws each)
#   - One-way sensitivity is second-longest (45-90 min)
#   - All other analyses combined: ~30 min
#
# To run only a subset of analyses, comment out the corresponding source()
# calls below.
#
# DATE: May 2026
#
################################################################################

# Optional: clear workspace and set wd
# rm(list = ls())
# setwd("~/your/working/directory")

# Load packages once at the top
suppressPackageStartupMessages({
  pkgs <- c("deSolve", "lhs", "ggplot2", "scales", "patchwork",
            "RColorBrewer", "viridis", "dplyr", "tidyr", "ggtext")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }
})

# Track total time
pipeline_start <- Sys.time()

cat("\n")
cat("########################################################################\n")
cat("#                                                                      #\n")
cat("#         TB VACCINATION MODEL — COMPLETE ANALYSIS PIPELINE            #\n")
cat("#                                                                      #\n")
cat("########################################################################\n\n")
cat(sprintf("Start: %s\n\n", format(pipeline_start, "%Y-%m-%d %H:%M:%S")))

# ---- Helper to run one step with timing ----------------------------------
run_step <- function(step_num, total_steps, script_path, description) {
  cat(sprintf("\n========================================================================\n"))
  cat(sprintf("STEP %d/%d: %s\n", step_num, total_steps, description))
  cat(sprintf("  source(\"%s\")\n", script_path))
  cat(sprintf("========================================================================\n\n"))

  t0 <- Sys.time()
  source(script_path)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  cat(sprintf("\n>>> Step %d completed in %.1f minutes\n", step_num, dt))
  invisible(dt)
}

n_total_steps <- 11

################################################################################
# STEP 1: Calibration extraction (scenario-independent)
################################################################################
run_step(1, n_total_steps, "TB_EXTRACT_CALIBRATION.R",
         "Extract calibrated parameter values")

################################################################################
# STEP 2: Strategy extraction (compartments, summary, case sources, LTBI, flows)
################################################################################
run_step(2, n_total_steps, "TB_EXTRACT_STRATEGIES.R",
         "Extract per-strategy outputs (both scenarios)")

################################################################################
# STEP 3: Stratum-level intervention impact
################################################################################
run_step(3, n_total_steps, "TB_STRATUM_IMPACT.R",
         "Stratum-level intervention impact (both scenarios)")

################################################################################
# STEP 4: Time-to-impact / transient dynamics
################################################################################
run_step(4, n_total_steps, "TB_TIME_TO_IMPACT.R",
         "Time-to-impact analysis (both scenarios)")

################################################################################
# STEP 5: Direct/indirect + psi sensitivity + threshold (3-in-1)
################################################################################
run_step(5, n_total_steps, "TB_DIRECT_INDIRECT_PSI_THRESHOLD.R",
         "Direct/indirect + psi + threshold (both scenarios)")

################################################################################
# STEP 6: Effective coverage
################################################################################
run_step(6, n_total_steps, "TB_EFFECTIVE_COVERAGE.R",
         "Effective coverage analysis (both scenarios)")

################################################################################
# STEP 7: NNV benchmarking (optimistic only - plausible NNV comes from earlier)
################################################################################
run_step(7, n_total_steps, "TB_NNV_BENCHMARK_OPTIMISTIC.R",
         "NNV benchmark equilibrium estimates (optimistic)")

################################################################################
# STEP 8: One-way sensitivity (LONG - 45-90 min)
################################################################################
run_step(8, n_total_steps, "TB_ONEWAY.R",
         "One-way sensitivity sweeps + published NNV (both scenarios)")

################################################################################
# STEP 9: LHS uncertainty + PRCC (LONG - 45-90 min)
################################################################################
run_step(9, n_total_steps, "TB_LHS.R",
         "LHS uncertainty intervals + PRCC (both scenarios)")

################################################################################
# STEP 10: Master output merge
################################################################################
run_step(10, n_total_steps, "TB_EXTRACT_MASTER.R",
         "Master output merge (both scenarios)")

################################################################################
# STEP 11: Alternative calibration (project-level diagnostic)
################################################################################
run_step(11, n_total_steps, "TB_CLEAN_ALTERNATIVE_CALIBRATION.R",
         "Alternative calibration sensitivity check")

################################################################################
# WRAP UP
################################################################################

pipeline_end <- Sys.time()
total_min <- as.numeric(difftime(pipeline_end, pipeline_start, units = "mins"))

cat("\n\n")
cat("########################################################################\n")
cat("#                          PIPELINE COMPLETE                           #\n")
cat("########################################################################\n\n")
cat(sprintf("End:   %s\n", format(pipeline_end, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Total: %.1f minutes (%.2f hours)\n\n", total_min, total_min / 60))
cat("All CSV outputs are in the working directory.\n")
cat("To regenerate figures, run the dissertation figures script next:\n")
cat("  source(\"TB_DISSERTATION_FIGURES.R\")\n\n")
