################################################################################
#
#  TB VACCINATION MODEL — MASTER OUTPUT MERGE (PLAUSIBLE + OPTIMISTIC)
#
################################################################################
#
# PURPOSE:
#   Merges the per-scenario strategy outputs (strategy summary + LTBI pool +
#   case sources) into single consolidated master files. Pure file merge, no
#   ODE runs. Used as a convenience aggregation for downstream consumers that
#   want all key per-strategy metrics in one place.
#
# REQUIRES:
#   The following input CSVs in working directory (produced by
#   TB_EXTRACT_STRATEGIES.R):
#     Plausible:
#       - TB_STRATEGY_SUMMARY.csv
#       - TB_LTBI_POOL_ANALYSIS.csv
#       - TB_CASE_SOURCES.csv
#     Optimistic:
#       - TB_STRATEGY_SUMMARY_OPTIMISTIC.csv
#       - TB_LTBI_POOL_ANALYSIS_OPTIMISTIC.csv
#       - TB_CASE_SOURCES_OPTIMISTIC.csv
#
# USAGE:
#   1. Set working directory to folder containing the input CSVs
#   2. Source this file: source("TB_EXTRACT_MASTER.R")
#   3. Master CSVs saved to working directory
#
#   To produce only one scenario's master file, comment out the unwanted
#   entry in the `scenarios` list below.
#
# OUTPUT FILES:
#   - TB_MODEL_OUTPUTS_COMPLETE.csv             (plausible)
#   - TB_MODEL_OUTPUTS_COMPLETE_OPTIMISTIC.csv  (optimistic)
#
# EXPECTED RUNTIME: < 5 seconds
#
# DATE: May 2026
#
################################################################################

cat("================================================================\n")
cat("  TB MODEL — MASTER OUTPUT MERGE\n")
cat("================================================================\n\n")

# ---- Define scenarios ----------------------------------------------------
scenarios <- list(
  plausible = list(
    label       = "PLAUSIBLE",
    csv_summary = "TB_STRATEGY_SUMMARY.csv",
    csv_ltbi    = "TB_LTBI_POOL_ANALYSIS.csv",
    csv_sources = "TB_CASE_SOURCES.csv",
    csv_master  = "TB_MODEL_OUTPUTS_COMPLETE.csv"
  ),
  optimistic = list(
    label       = "OPTIMISTIC",
    csv_summary = "TB_STRATEGY_SUMMARY_OPTIMISTIC.csv",
    csv_ltbi    = "TB_LTBI_POOL_ANALYSIS_OPTIMISTIC.csv",
    csv_sources = "TB_CASE_SOURCES_OPTIMISTIC.csv",
    csv_master  = "TB_MODEL_OUTPUTS_COMPLETE_OPTIMISTIC.csv"
  )
)

# ---- Helper: try to read a CSV, return NULL on miss ----------------------
read_or_null <- function(path) {
  if (!file.exists(path)) {
    cat(sprintf("  WARNING: %s not found - skipping merge\n", path))
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

# ---- Process each scenario -----------------------------------------------
for (scen_name in names(scenarios)) {

  scen <- scenarios[[scen_name]]

  cat(sprintf("\nProcessing scenario: %s\n", scen$label))
  cat("----------------------------------------\n")

  strategy_summary <- read_or_null(scen$csv_summary)
  ltbi_summary     <- read_or_null(scen$csv_ltbi)
  case_sources     <- read_or_null(scen$csv_sources)

  # If any required input is missing, skip this scenario
  if (is.null(strategy_summary) || is.null(ltbi_summary) || is.null(case_sources)) {
    cat(sprintf("  Skipping %s - one or more inputs missing\n", scen$label))
    next
  }

  # Build the master output by merging on Scenario
  master_output <- merge(strategy_summary, ltbi_summary,
                         by = "Scenario", all = TRUE)

  master_output <- merge(
    master_output,
    case_sources[, c("Scenario", "Pct_from_Lf", "Pct_from_Ls", "Pct_reactivation")],
    by = "Scenario", all = TRUE
  )

  write.csv(master_output, scen$csv_master, row.names = FALSE)
  cat(sprintf("  Wrote: %s (%d strategies)\n",
              scen$csv_master, nrow(master_output)))
}

cat("\n================================================================\n")
cat("  MERGE COMPLETE\n")
cat("================================================================\n")
