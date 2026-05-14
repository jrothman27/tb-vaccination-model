################################################################################
#
#  TB VACCINATION MODEL -- MANUSCRIPT AND APPENDIX FIGURES
#
################################################################################
#
# Generates publication-quality figures for Aim 3 dissertation chapter.
# All figures saved as PNG (300 dpi) for embedding in Word/LaTeX.
#
# 
# FIGURE LIST:
#   Fig 1:  Model structure diagram (TB_Transmission_Model.png -- pre-made)
#   Fig 2:  Case source validation (model emergent vs. literature ranges)
#   Fig S:  Recent transmission validation (appendix, NTSS genotyping)
#   Fig 3:  Strategy comparison [manuscript Figure 2]
#           2 rows (cases prevented top / NNV bottom),
#           each with faceted panels: Optimistic (left) | Plausible (right)
#   Fig 4:  Stratum-level impact breakdown (grouped bar by strategy x stratum)
#   Fig 5:  Incidence trajectory over time (baseline vs. vaccination scenarios)
#   Fig 6:  Cumulative cases prevented over time
#   Fig 7:  Time-to-impact: annual cases prevented trajectory [was Fig 5]
#   Fig 8:  Tornado diagram (one-way sensitivity) [was Fig 2]
#   Fig 9:  LHS forest plot (95% uncertainty intervals) [was Fig 3]
#   Fig 10: PRCC bar plot [was Fig 4]
#   Fig 11: Direct vs. indirect effects [was Fig 6]
#   Fig 12: VE x Duration threshold heatmap [manuscript Figure 3]
#           Two-panel: optimistic (left, 50%/yr), plausible (right, 5%/yr)
#   Fig 13: NNV comparison with established vaccines
#   Fig 14: Psi dose-response [was Fig 8]
#
# REQUIRES:
#   - ggplot2, scales, patchwork, RColorBrewer, viridis, dplyr, tidyr
#   - RData/CSV files from previous analyses (in working directory)
#
# NOTE: All labels use plain ASCII to avoid UTF-8 encoding issues.
#
################################################################################

# Install packages if needed
for (pkg in c("ggplot2", "scales", "patchwork", "RColorBrewer", "viridis",
              "dplyr", "tidyr", "stringr", "ggtext")) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

cat("================================================================\n")
cat("  GENERATING DISSERTATION FIGURES (v4, 4/21/26)\n")
cat("================================================================\n\n")

# --- Load data from previous analyses (assumes WD is set) ---

# Load tornado data from CSV (was previously in RData)
tornado <- read.csv("TB_ONEWAY_TORNADO.csv")
cat("Loaded: TB_ONEWAY_TORNADO.csv\n")

# Load CSVs
lhs_ui    <- read.csv("TB_LHS_UNCERTAINTY_INTERVALS.csv")
prcc      <- read.csv("TB_PRCC_RESULTS.csv")
tti       <- read.csv("TB_TIME_TO_IMPACT_RESULTS.csv")
di        <- read.csv("TB_DIRECT_INDIRECT_EFFECTS.csv")
psi_s     <- read.csv("TB_PSI_SENSITIVITY.csv")
thresh    <- read.csv("TB_THRESHOLD_ANALYSIS.csv")
thresh_mv <- read.csv("TB_THRESHOLD_MINVE.csv")

# Optimistic threshold grid (psi = 50%/yr) -- optional; Fig 12 falls back to
# single-panel if this file is not present
thresh_opt <- tryCatch(
  read.csv("TB_THRESHOLD_ANALYSIS_OPTIMISTIC.csv"),
  error = function(e) {
    cat("  Note: TB_THRESHOLD_ANALYSIS_OPTIMISTIC.csv not found -- Fig 12 will be single-panel.\n")
    cat("        Run TB_OPTIMISTIC_THRESHOLD_ANALYSIS.R to generate it.\n")
    NULL
  }
)

# Optimistic scenario strategy outcomes (VE = 70%, psi = 50%/yr) -- required for
# the two-panel Fig 3. Contains all strategies at psi sweep; we filter to
# VE = 70, psi = 50 below. Falls back to single-panel plausible-only if absent.
opt_strat <- tryCatch(
  read.csv("TB_OPTIMISTIC_SCENARIOS.csv"),
  error = function(e) {
    cat("  Note: TB_OPTIMISTIC_SCENARIOS.csv not found -- Fig 3 will be plausible-only.\n")
    NULL
  }
)

cal       <- read.csv("TB_CALIBRATION_RESULTS.csv")
strat_imp <- read.csv("TB_STRATUM_LEVEL_IMPACT.csv")
strat_sum <- read.csv("TB_STRATEGY_SUMMARY.csv")
nnv_pub   <- read.csv("TB_PUBLISHED_NNV_COMPARISON.csv")

# Optimistic-scenario stratum-level impact (VE = 70%, psi = 50%/yr). Required
# for the two-panel Figure S2 (manuscript Fig 4). Falls back to plausible-only
# if absent. Generate via TB_OPTIMISTIC_STRATUM_IMPACT.R.
strat_imp_opt <- tryCatch(
  read.csv("TB_STRATUM_LEVEL_IMPACT_OPTIMISTIC.csv"),
  error = function(e) {
    cat("  Note: TB_STRATUM_LEVEL_IMPACT_OPTIMISTIC.csv not found -- Fig 4 will be plausible-only.\n")
    cat("        Run TB_OPTIMISTIC_STRATUM_IMPACT.R to generate it.\n")
    NULL
  }
)

cat("Loaded: All CSV data files\n\n")


###############################################################################
# ENCODING FIX: Replace garbled tornado labels with clean ASCII
###############################################################################

clean_tornado_labels <- c(
  # Tornado / one-way params
  "eps_s_FB"    = "Non-U.S.-born reactivation rate (eps_s, /yr)",
  "nu_s_clear"  = "Mtb infection clearance rate (nu_s_clear, /yr)",
  "eps_s_scale" = "Slow reactivation (eps_s, x multiplier)",
  "eps_s_USB"   = "U.S.-born reactivation rate (eps_s, /yr)",
  "gamma"       = "Treatment rate (gamma, /yr)",
  "epsilon_mix" = "Assortative mixing (eps_mix)",
  "p_scale"     = "Fraction to fast Mtb infection (p, x multiplier)",
  "eps_f_base"  = "Fast progression base (eps_f_base, /yr)",
  "nu_f_trans"  = "L_f to L_s transition (nu_f_trans, /yr)",
  "alpha_stab"  = "Post-treatment stabilization (alpha_stab, /yr)",
  "sigma"       = "Reinfection protection (sigma)",
  # PRCC-specific params
  "gamma_tx"           = "Treatment rate (gamma, /yr)",
  "eps_s_HIV"          = "PLWH reactivation rate (eps_s, /yr)",
  "eps_s_Med"          = "Medical reactivation rate (eps_s, /yr)",
  "eps_s_Non-U.S.-born"= "Non-U.S.-born reactivation rate (eps_s, /yr)",
  "omega"              = "Vaccine waning rate (omega, /yr)",
  "VE"                 = "Vaccine efficacy (VE)"
)

if (exists("tornado")) {
  for (i in seq_len(nrow(tornado))) {
    param <- as.character(tornado$parameter[i])
    if (param %in% names(clean_tornado_labels)) {
      tornado$label[i] <- clean_tornado_labels[param]
    }
  }
  cat("Fixed: tornado labels cleaned (ASCII)\n\n")
}


###############################################################################
# COMMON THEME AND PALETTES
###############################################################################

theme_diss <- theme_bw(base_size = 11) +
  theme(
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10),
    plot.title       = element_text(size = 13, face = "bold", hjust = 0,
                                     margin = margin(b = 6)),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 10),
    strip.text       = element_text(size = 10, face = "bold"),
    plot.margin      = margin(10, 15, 10, 10)
  )

# Strategy colors -- Revised Nature Medicine palette (colorblind-accessible)
# 3 swaps from original NM pastels: darker blue, deeper rose, richer purple
strat_colors <- c(
  "All High-Risk"          = "#2B6A99",
  "PLWH + Non-U.S.-Born"   = "#5DA88C",
  "All Non-U.S.-Born"      = "#7EC4CF",
  "All Mtb-Infected"       = "#F7D08A",
  "PLWH + Medical"         = "#F2A586",
  "Medical"                = "#C44E6C",
  "PLWH"                   = "#7B68AE"
)

# Strategy display order -- ranked by cases prevented (most impact first)
strat_order <- c("All High-Risk", "PLWH + Non-U.S.-Born", "All Mtb-Infected",
                 "PLWH + Medical", "Medical", "PLWH")

# Extended order including exploratory All Non-U.S.-Born
strat_order_full <- c("All High-Risk", "PLWH + Non-U.S.-Born",
                      "All Non-U.S.-Born", "All Mtb-Infected",
                      "PLWH + Medical", "Medical", "PLWH")

# Strategy point shapes (for colorblind accessibility on line charts)
# circle=16, triangle=17, square=15, diamond=18, star=8, cross=4, plus=3
strat_shapes <- c(
  "All High-Risk"          = 16,
  "PLWH + Non-U.S.-Born"   = 15,
  "All Non-U.S.-Born"      = 17,
  "All Mtb-Infected"       = 18,
  "PLWH + Medical"         = 8,
  "Medical"                = 4,
  "PLWH"                   = 3
)

# Stratum colors -- matched to revised NM palette
stratum_colors <- c(
  "HIV"           = "#7B68AE",
  "Medical"       = "#F2A586",
  "Non-U.S.-born" = "#2B6A99",
  "U.S.-born"     = "#91D0A4"
)

outdir <- "."

# Helper: clean strategy names consistently
clean_strategy <- function(x) {
  x <- gsub(" \\(HIV\\+Med\\+FB\\)", "", x)
  x <- gsub(" \\(2%/yr\\)", "", x)
  # Longest-to-shortest replacement order to avoid partial matches
  x <- gsub("PLWH \\+ All NUSB", "PLWH + Non-U.S.-Born", x)
  x <- gsub("All Foreign-Born", "All Non-U.S.-Born", x)
  x <- gsub("FB Stratum Only", "NUSB Stratum", x)
  x <- gsub("Medical Only", "Medical", x)
  x <- gsub("^Medical comorbidities$", "Medical", x)  # v9: optimistic UI uses full name
  x <- gsub("HIV Only", "PLWH", x)
  x <- gsub("HIV \\+ Medical", "PLWH + Medical", x)
  x <- gsub("^Universal$", "All Mtb-Infected", x)
  x <- gsub("All High-Risk", "All High-Risk", x)
  return(x)
}

# Drop Non-U.S.-Born-only exploratory strategies from primary figures.
# These are not primary policy options and should only appear in exploratory tables.
# Patterns:
#   - "All NUSB" / "All Non-U.S.-Born": vaccinate only the NUSB stratum
#   - "NUSB Stratum Only" / "NUSB Only": same, alternative naming
#   - "FB Stratum Only": legacy raw name for the same strategy
drop_nusb_only <- function(df, col = "strategy") {
  df[!grepl("^All NUSB$|^All Non-U\\.S\\.-Born$|NUSB Stratum|NUSB Only|FB Stratum Only", df[[col]]), ]
}


########################################################################
# FIGURE 1: Model Structure Diagram
########################################################################
# NOTE: This is a pre-made PNG (TB_Transmission_Model.png).
# It is included in the LaTeX document directly via \includegraphics.
# No R code needed -- just copy to output directory with Fig numbering.

cat("Figure 1: Model structure diagram (pre-made PNG)...\n")
if (file.exists("TB_Transmission_Model.png")) {
  file.copy("TB_Transmission_Model.png", file.path(outdir, "Fig1_model_diagram.png"),
            overwrite = TRUE)
  cat("  Copied to Fig1_model_diagram.png\n")
} else {
  cat("  WARNING: TB_Transmission_Model.png not found. Skipping.\n")
}



########################################################################
# FIGURE 2: Case Source Validation (main text)
#   Model-emergent case source proportions vs. literature ranges
########################################################################

cat("Figure 2: Case source validation...\n")

# Load case sources for validation
case_src  <- read.csv("TB_CASE_SOURCES.csv")
ann_flows <- read.csv("TB_ANNUAL_FLOWS.csv")

cs_base <- case_src[case_src$Scenario == "Baseline", ]

case_df <- data.frame(
  Source_type = factor(c("Recent\ntransmission\n(Lf \u2192 I)",
                         "Remote\nreactivation\n(Ls \u2192 I)",
                         "Post-treatment\nrelapse\n(R \u2192 I)"),
                       levels = c("Recent\ntransmission\n(Lf \u2192 I)",
                                  "Remote\nreactivation\n(Ls \u2192 I)",
                                  "Post-treatment\nrelapse\n(R \u2192 I)")),
  Model_pct = c(cs_base$Pct_from_Lf, cs_base$Pct_from_Ls, cs_base$Pct_from_relapse),
  Lit_lo    = c(10, 85, 3),
  Lit_hi    = c(15, 90, 6)
)

fig2 <- ggplot(case_df, aes(x = Source_type)) +
  geom_rect(aes(xmin = as.numeric(Source_type) - 0.35,
                xmax = as.numeric(Source_type) + 0.35,
                ymin = Lit_lo, ymax = Lit_hi),
            fill = "#D9D9D9", alpha = 0.7) +
  geom_point(aes(y = Model_pct), shape = 18, size = 6, color = "#C44E6C") +
  geom_text(aes(y = Model_pct, label = sprintf("%.1f%%", Model_pct)),
            vjust = -1.2, size = 4, color = "#C44E6C", fontface = "bold") +
  geom_text(aes(y = (Lit_lo + Lit_hi) / 2,
                label = paste0(Lit_lo, "\u2013", Lit_hi, "%")),
            hjust = 1.4, size = 3.5, color = "grey50", fontface = "italic") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Proportion of TB Cases (%)") +
  theme_diss +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 12)
  ) +
  annotate("point", x = 2.6, y = 97, shape = 18, size = 5, color = "#C44E6C") +
  annotate("text", x = 2.72, y = 97, label = "Model (emergent)",
           size = 3.5, hjust = 0, color = "#C44E6C") +
  annotate("rect", xmin = 2.47, xmax = 2.6, ymin = 90.5, ymax = 93.5,
           fill = "#D9D9D9") +
  annotate("text", x = 2.72, y = 92, label = "Literature range",
           size = 3.5, hjust = 0, color = "grey50")

ggsave(file.path(outdir, "Fig2_case_source_validation.pdf"), fig2,
       width = 7, height = 5.5, device = cairo_pdf, bg = "white")


########################################################################
# SUPPLEMENTARY FIGURE: Recent Transmission Validation (appendix)
#   Stratum-level Lf->I fraction vs. NTSS genotyping estimates
########################################################################

cat("Figure S: Recent transmission validation...\n")

base_flows <- ann_flows[ann_flows$Scenario == "Baseline", ]
base_flows$total_new <- base_flows$Flow_Lf_to_I + base_flows$Flow_Ls_to_I +
                         base_flows$Flow_Lfv_to_I + base_flows$Flow_Lsv_to_I +
                         base_flows$Flow_R_to_I
base_flows$pct_Lf <- 100 * base_flows$Flow_Lf_to_I / base_flows$total_new

rt_df <- data.frame(
  Group = rep(c("Non-U.S.-\nborn", "U.S.-\nborn", "HIV+", "Overall"), 2),
  Source = rep(c("NTSS Genotyping\n(2016-2019)", "Model (Lf \u2192 I)"), each = 4),
  Pct = c(
    7.8, 25.5, 16.5, 13.0,
    base_flows$pct_Lf[base_flows$Stratum == "FB"],
    base_flows$pct_Lf[base_flows$Stratum == "USB"],
    base_flows$pct_Lf[base_flows$Stratum == "HIV"],
    cs_base$Pct_from_Lf
  )
)

rt_df$Group <- factor(rt_df$Group,
                      levels = c("Non-U.S.-\nborn", "U.S.-\nborn", "HIV+", "Overall"))
rt_df$Source <- factor(rt_df$Source,
                       levels = c("NTSS Genotyping\n(2016-2019)", "Model (Lf \u2192 I)"))

fig_s_rt <- ggplot(rt_df, aes(x = Group, y = Pct, fill = Source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Pct)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("NTSS Genotyping\n(2016-2019)" = "#B0B0B0",
                                "Model (Lf \u2192 I)" = "#C44E6C")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  annotate("segment", x = 0.65, xend = 1.35, y = 41, yend = 41,
           color = "#2B6A99", linewidth = 1.5) +
  annotate("text", x = 1, y = 43, label = "7.8% vs 8.2%",
           color = "#2B6A99", size = 3.5, fontface = "bold") +
  labs(x = NULL, y = "% Recent Transmission",
       fill = NULL) +
  theme_diss +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 11)
  )

ggsave(file.path(outdir, "Fig_S_recent_transmission.pdf"), fig_s_rt,
       width = 7, height = 5.5, device = cairo_pdf, bg = "white")

########################################################################
# FIGURE 3: Strategy Comparison [manuscript Figure 2]
#   Two-row faceted layout (facet_wrap on Scenario within each row):
#     a. Annual cases prevented       (facets: Optimistic | Plausible)
#     b. Number needed to vaccinate   (facets: Optimistic | Plausible)
#   One bold row-title per row; gray facet strips carry the scenario text.
#   Strategies ordered top-to-bottom by cases prevented under OPTIMISTIC.
#   Strategy colours consistent across all four panels.
#   Falls back to plausible-only layout (no facets) if optimistic CSV absent.
########################################################################

cat("Figure 3: Strategy comparison (2-row faceted)...\n")

# Scenario strip labels (gray facet headers) -- two-line text to match
# the tight header block in the target layout
scen_opt_label   <- "Optimistic (VE = 70%,\nvaccination rate = 50%/yr)"
scen_plaus_label <- "Plausible (VE = 50%,\nstrategy-specific vaccination rate)"

# ---- Plausible panel data (VE = 50%, duration = 10 yr, strategy-specific psi) ----
fig3_plaus_raw <- thresh[thresh$VE == 50 & thresh$duration_yr == 10, ]
fig3_plaus_raw$strategy_short <- clean_strategy(fig3_plaus_raw$strategy)
fig3_plaus_raw <- drop_nusb_only(fig3_plaus_raw, "strategy_short")

# Italicize "Mtb" in the display label (rendered via ggtext::element_markdown on
# the y-axis). Strategy name in `Strategy` column is kept plain-text so it can
# still be used as a factor and matched against strat_colors.
italicize_mtb <- function(x) gsub("Mtb", "*Mtb*", x, fixed = TRUE)

fig3_plaus <- data.frame(
  Strategy  = fig3_plaus_raw$strategy_short,
  Prevented = fig3_plaus_raw$cases_prevented,
  NNV       = fig3_plaus_raw$NNV,
  Pct       = round(fig3_plaus_raw$pct_reduction, 1),
  Scenario  = scen_plaus_label
)

# ---- Optimistic panel data (VE = 70%, psi = 50%/yr) ----
if (!is.null(opt_strat)) {
  fig3_opt_raw <- opt_strat[opt_strat$VE_pct == 70 & opt_strat$psi_pct == 50, ]
  fig3_opt_raw$strategy_short <- clean_strategy(fig3_opt_raw$strategy)
  fig3_opt_raw <- drop_nusb_only(fig3_opt_raw, "strategy_short")

  fig3_opt <- data.frame(
    Strategy  = fig3_opt_raw$strategy_short,
    Prevented = fig3_opt_raw$cases_prevented,
    NNV       = fig3_opt_raw$NNV,
    Pct       = round(fig3_opt_raw$pct_reduction, 1),
    Scenario  = scen_opt_label
  )

  # Ordering: top-to-bottom by cases prevented under OPTIMISTIC
  opt_order <- fig3_opt$Strategy[order(-fig3_opt$Prevented)]
  have_opt <- TRUE
} else {
  fig3_opt  <- fig3_plaus[0, ]
  opt_order <- fig3_plaus$Strategy[order(-fig3_plaus$Prevented)]
  have_opt  <- FALSE
}

# Combine into a single long-format frame. factor(Scenario) sets panel
# order: Optimistic (left), Plausible (right).
if (have_opt) {
  fig3_long <- rbind(fig3_opt, fig3_plaus)
  scenario_levels <- c(scen_opt_label, scen_plaus_label)
} else {
  fig3_long <- fig3_plaus
  scenario_levels <- scen_plaus_label
}
fig3_long$Scenario <- factor(fig3_long$Scenario, levels = scenario_levels)
# coord_flip() inverts the y-axis, so reverse factor order
fig3_long$Strategy <- factor(fig3_long$Strategy, levels = rev(opt_order))
fig3_long <- fig3_long[!is.na(fig3_long$Strategy), ]

# Shared x-axis maxima within each metric (read directly across facets).
# Small headroom + right-hand label offset prevents the "(X%)" text from
# clipping past the axis without requiring an extra-wide figure.
prev_max <- max(fig3_long$Prevented, na.rm = TRUE)
nnv_max  <- max(fig3_long$NNV,       na.rm = TRUE)

prev_label_offset <- prev_max * 0.015
nnv_label_offset  <- nnv_max  * 0.015

# Pct labels: integer-like values drop trailing ".0" (e.g. "7%", not "7.0%"),
# matching the target image. Non-integer values keep one decimal ("57.8%").
# as.character() avoids format()'s right-align padding that would produce
# "( 7%)" with a leading space instead of "(7%)".
fmt_pct <- function(p) {
  ifelse(abs(p - round(p)) < 0.05,
         paste0(as.character(round(p)), "%"),
         paste0(formatC(round(p, 1), format = "f", digits = 1), "%"))
}

# ---- Shared row-theme tweaks ------------------------------------
# Cleaner than the previous 2x2 grid: no y-axis line/ticks (strategy names
# speak for themselves), no horizontal gridlines (bars are self-anchoring),
# tighter strip/margin spacing, slightly larger y-text for readability.
# axis.text.y uses element_markdown so "*Mtb*" -> italic Mtb in strategy labels.
row_theme <- theme_diss +
  theme(
    plot.title         = element_text(size = 14, face = "bold", hjust = 0,
                                      margin = margin(b = 8)),
    strip.background   = element_rect(fill = "grey94", colour = NA),
    strip.text         = element_text(size = 10.5, face = "bold",
                                      colour = "grey15",
                                      margin = margin(t = 6, b = 6)),
    panel.spacing.x    = unit(1.2, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    axis.title.x       = element_text(size = 10.5, colour = "grey25",
                                      margin = margin(t = 8), hjust = 0.5),
    axis.text.y        = ggtext::element_markdown(size = 10.5, colour = "grey15"),
    axis.text.x        = element_text(size = 10,   colour = "grey25"),
    axis.ticks.y       = element_blank(),
    axis.line.y        = element_blank()
  )

# ---- Row a: Annual cases prevented -------------------------------
# Label format: "5,385 (51.8%)" -- the number first, then percent reduction in
# parens. Integer-like percents drop trailing ".0" via fmt_pct().
p_cases <- ggplot(fig3_long,
                  aes(x = Strategy, y = Prevented, fill = Strategy)) +
  geom_col(width = 0.72, colour = NA) +
  geom_text(aes(y = Prevented + prev_label_offset,
                label = paste0(scales::comma(Prevented), " (", fmt_pct(Pct), ")")),
            hjust = 0, size = 3.5, colour = "grey15") +
  coord_flip(ylim = c(0, prev_max * 1.38), clip = "off") +
  scale_x_discrete(labels = italicize_mtb) +
  scale_fill_manual(values = strat_colors, guide = "none") +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL,
       y = "Annual Cases Prevented\n(% reduction in overall TB burden)",
       title = "a. Annual cases prevented") +
  row_theme +
  theme(plot.margin = margin(5, 14, 10, 5))

# Facet only when we have both scenarios; otherwise show a single panel
# without a strip header.
if (have_opt) {
  p_cases <- p_cases + facet_wrap(~ Scenario, nrow = 1)
}

# ---- Row b: Number needed to vaccinate ---------------------------
p_nnv <- ggplot(fig3_long,
                aes(x = Strategy, y = NNV, fill = Strategy)) +
  geom_col(width = 0.72, colour = NA) +
  geom_text(aes(y = NNV + nnv_label_offset,
                label = scales::comma(NNV)),
            hjust = 0, size = 3.7, colour = "grey15") +
  coord_flip(ylim = c(0, nnv_max * 1.18), clip = "off") +
  scale_x_discrete(labels = italicize_mtb) +
  scale_fill_manual(values = strat_colors, guide = "none") +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Number Needed to Vaccinate",
       title = "b. Number needed to vaccinate") +
  row_theme +
  theme(plot.margin = margin(10, 14, 5, 5))

if (have_opt) {
  p_nnv <- p_nnv + facet_wrap(~ Scenario, nrow = 1)
}

# ---- Assemble two rows --------------------------------------------
fig3 <- p_cases / p_nnv + plot_layout(heights = c(1, 1))

# Wider when we have two facets; narrower for plausible-only fallback.
fig3_width  <- if (have_opt) 9   else 6.5
fig3_height <- 7.4

ggsave(file.path(outdir, "Fig3_strategy_comparison.pdf"), fig3,
       width = fig3_width, height = fig3_height,
       device = cairo_pdf, bg = "white")
ggsave(file.path(outdir, "Fig3_strategy_comparison.png"), fig3,
       width = fig3_width, height = fig3_height, dpi = 300, bg = "white")

cat(sprintf("  Fig 3 saved (2-row faceted, %s)\n",
            if (have_opt) "optimistic + plausible" else "plausible only"))


########################################################################
# FIGURE 4: Stratum-Level Impact Breakdown  [manuscript Figure S2]
#   Two-column facet layout (plausible left, optimistic right):
#     Row a (top):    Annual cases prevented by stratum (stacked)
#     Row b (bottom): Stratum-level NNV (dodged bars)
#   Falls back to single-column (plausible-only) if optimistic CSV is absent.
########################################################################

cat("Figure 4: Stratum-level impact breakdown (two-scenario)...\n")

# Build combined long-format data frame with a Scenario column
si_plaus <- strat_imp[strat_imp$Scenario != "Baseline", ]
si_plaus$Strategy <- clean_strategy(si_plaus$Scenario)
si_plaus <- drop_nusb_only(si_plaus, "Strategy")
si_plaus$Scenario_Type <- "Plausible"

if (!is.null(strat_imp_opt)) {
  si_opt <- strat_imp_opt[strat_imp_opt$Scenario != "Baseline", ]
  si_opt$Strategy <- clean_strategy(si_opt$Scenario)
  si_opt <- drop_nusb_only(si_opt, "Strategy")
  si_opt$Scenario_Type <- "Optimistic"

  # Retain a consistent column set before rbind (defensive: drop cols that
  # differ, keep only the union minimum)
  common_cols <- intersect(colnames(si_plaus), colnames(si_opt))
  si <- rbind(si_plaus[, common_cols], si_opt[, common_cols])
} else {
  si <- si_plaus
}

si$Strategy      <- factor(si$Strategy, levels = strat_order)
si               <- si[!is.na(si$Strategy), ]
si$Stratum       <- factor(si$Stratum,
                           levels = c("HIV", "Medical", "FB", "USB"),
                           labels = c("HIV", "Medical", "Non-U.S.-born", "U.S.-born"))
si$Scenario_Type <- factor(si$Scenario_Type, levels = c("Optimistic", "Plausible"))

# Facet labels with inline scenario details
scenario_labels_s2 <- c(
  "Optimistic" = "Optimistic: VE = 70%, \u03c8 = 50%/yr",
  "Plausible"  = "Plausible: VE = 50%, strategy-specific \u03c8"
)

# Panel A: cases prevented by stratum (stacked)
fig4a <- ggplot(si, aes(x = Strategy, y = Cases_prevented, fill = Stratum)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = stratum_colors) +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~ Scenario_Type, nrow = 1,
             labeller = labeller(Scenario_Type = scenario_labels_s2),
             scales = "free_x") +
  labs(x = NULL, y = "Annual Cases Prevented",
       fill = "Stratum", title = "a. Annual cases prevented by stratum") +
  theme_diss +
  theme(legend.position  = "right",
        strip.text       = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "grey95", color = NA))

# Panel B: stratum-level NNV (exclude strata with no vaccinations)
si_nnv <- si[si$Annual_vaccinations > 0 & !is.na(si$NNV) & si$NNV > 0, ]

fig4b <- ggplot(si_nnv, aes(x = Strategy, y = NNV, fill = Stratum)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = stratum_colors) +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~ Scenario_Type, nrow = 1,
             labeller = labeller(Scenario_Type = scenario_labels_s2),
             scales = "free_x") +
  labs(x = NULL, y = "Stratum-Level NNV",
       fill = "Stratum", title = "b. Stratum-level number needed to vaccinate") +
  theme_diss +
  theme(legend.position  = "right",
        strip.text       = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "grey95", color = NA))

fig4 <- fig4a / fig4b

# Wider canvas when both scenarios are rendered
fig4_width  <- if (!is.null(strat_imp_opt)) 11 else 8
fig4_height <- if (!is.null(strat_imp_opt)) 9  else 8

ggsave(file.path(outdir, "Fig4_stratum_impact.pdf"), fig4,
       width = fig4_width, height = fig4_height,
       device = cairo_pdf, bg = "white")
ggsave(file.path(outdir, "Fig4_stratum_impact.png"), fig4,
       width = fig4_width, height = fig4_height, dpi = 300, bg = "white")

cat(sprintf("  Fig 4 saved (%s panel%s)\n",
            if (!is.null(strat_imp_opt)) "two" else "one",
            if (!is.null(strat_imp_opt)) "s" else ""))


########################################################################
# FIGURE 5: Incidence Trajectory Over Time (REVISED)
#   Panel A: Total annual TB cases by strategy
#   Panel B: Indirect effects -- cases prevented in U.S.-born (untargeted)
#   Legend ordered by impact (most cases prevented first)
########################################################################

cat("Figure 5: Incidence trajectory (two-panel)...\n")

tti_plot <- tti
tti_plot$strategy_short <- clean_strategy(tti_plot$strategy)
tti_plot <- drop_nusb_only(tti_plot, "strategy_short")
tti_plot$strategy_short <- factor(tti_plot$strategy_short, levels = strat_order)
tti_plot <- tti_plot[!is.na(tti_plot$strategy_short), ]

# --- Panel A: Total incidence trajectories ---

# Baseline (same for all strategies)
inc_base <- data.frame(
  year     = tti_plot$year[tti_plot$strategy_short == "All High-Risk"],
  cases    = tti_plot$baseline_cases[tti_plot$strategy_short == "All High-Risk"],
  Scenario = "No Vaccination (Baseline)"
)

# Intervention trajectories
inc_vax <- data.frame(
  year     = tti_plot$year,
  cases    = tti_plot$intervention_cases,
  Scenario = as.character(tti_plot$strategy_short)
)

inc_all <- rbind(inc_base, inc_vax)
inc_all$Scenario <- factor(inc_all$Scenario,
                           levels = c("No Vaccination (Baseline)", strat_order))

# Color, shape, and linetype palettes including baseline
inc_colors <- c("No Vaccination (Baseline)" = "black", strat_colors[strat_order])
inc_shapes <- c("No Vaccination (Baseline)" = 1, strat_shapes[strat_order])
inc_linetypes <- setNames(
  c("dashed", rep("solid", length(strat_order))),
  c("No Vaccination (Baseline)", strat_order)
)

panelA <- ggplot(inc_all, aes(x = year, y = cases, color = Scenario,
                               linetype = Scenario, shape = Scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  scale_color_manual(values = inc_colors) +
  scale_shape_manual(values = inc_shapes) +
  scale_linetype_manual(values = inc_linetypes) +
  scale_x_continuous(breaks = c(1, 2, 3, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Years After Vaccination Begins",
       y = "Annual TB Cases",
       color = "Scenario", linetype = "Scenario", shape = "Scenario",
       title = "a") +
  theme_diss +
  theme(legend.position = "right",
        legend.key.width = unit(1.5, "cm"))

# --- Panel B: Indirect effects on U.S.-born (untargeted population) ---

# prev_USB = cases prevented in the U.S.-born stratum at each year
# For strategies that don't target USB, these are entirely indirect effects
usb_indirect <- tti_plot[, c("year", "strategy_short", "prev_USB")]
names(usb_indirect) <- c("year", "Strategy", "cases_prevented_usb")

# All strategies except All Mtb-Infected do not target USB,
# so prev_USB is fully indirect for those. All Mtb-Infected targets USB
# directly, so we exclude it from the indirect effects panel.
usb_indirect <- usb_indirect[usb_indirect$Strategy != "All Mtb-Infected", ]
usb_indirect$Strategy <- factor(usb_indirect$Strategy,
                                 levels = strat_order[strat_order != "All Mtb-Infected"])

panelB <- ggplot(usb_indirect, aes(x = year, y = cases_prevented_usb,
                                    color = Strategy, shape = Strategy)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_manual(values = strat_colors) +
  scale_shape_manual(values = strat_shapes) +
  scale_x_continuous(breaks = c(1, 2, 3, 5, 10, 15, 20, 25, 30)) +
  labs(x = "Years After Vaccination Begins",
       y = "Cases Prevented in U.S.-Born",
       title = "b") +
  theme_diss +
  theme(legend.position = "right",
        legend.key.width = unit(1.5, "cm"))

fig5 <- panelA / panelB

ggsave(file.path(outdir, "Fig5_incidence_trajectory.pdf"), fig5,
       width = 10, height = 10, device = cairo_pdf, bg = "white")


########################################################################
# FIGURE 6: Cumulative Cases Prevented Over Time (NEW)
########################################################################

cat("Figure 6: Cumulative cases prevented...\n")

tti_cumul <- tti_plot

fig6 <- ggplot(tti_cumul, aes(x = year, y = cumul_prevented,
                               color = strategy_short, shape = strategy_short)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  scale_color_manual(values = strat_colors, name = "Strategy") +
  scale_shape_manual(values = strat_shapes, name = "Strategy") +
  scale_x_continuous(breaks = c(1, 2, 3, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Years After Vaccination Begins",
       y = "Cumulative Cases Prevented") +
  theme_diss +
  theme(legend.position = "right")

ggsave(file.path(outdir, "Fig6_cumulative_prevented.pdf"), fig6,
       width = 9, height = 5, device = cairo_pdf, bg = "white")


########################################################################
# FIGURE 7: Time-to-Impact Curves [was Fig 5]
#   Annual cases prevented per year (approach to equilibrium)
########################################################################

cat("Figure 7: Time-to-impact...\n")

fig7 <- ggplot(tti_cumul, aes(x = year, y = cases_prevented,
                               color = strategy_short, shape = strategy_short)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  scale_color_manual(values = strat_colors, name = "Strategy") +
  scale_shape_manual(values = strat_shapes, name = "Strategy") +
  scale_x_continuous(breaks = c(1, 2, 3, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Years After Vaccination Begins",
       y = "Cases Prevented per Year") +
  theme_diss +
  theme(legend.position = "right")

ggsave(file.path(outdir, "Fig7_time_to_impact.pdf"), fig7,
       width = 9, height = 5, device = cairo_pdf, bg = "white")


########################################################################
# FIGURE 8: Tornado Diagram (One-Way Sensitivity) [was Fig 2]
########################################################################

cat("Figure 8: Tornado diagram...\n")

torn <- tornado[order(tornado$range_prevented), ]
torn$label <- factor(torn$label, levels = torn$label)

torn$low_delta  <- torn$min_prevented - torn$baseline_prev
torn$high_delta <- torn$max_prevented - torn$baseline_prev

fig8_data <- data.frame(
  Parameter = rep(torn$label, 2),
  Direction = rep(c("Low", "High"), each = nrow(torn)),
  Delta     = c(torn$low_delta, torn$high_delta)
)

fig8 <- ggplot(fig8_data, aes(x = Parameter, y = Delta, fill = Direction)) +
  geom_col(position = "identity", width = 0.65) +
  coord_flip() +
  scale_fill_manual(values = c("Low" = "#F2A586", "High" = "#2B6A99"),
                    labels = c("Low" = "Lower bound", "High" = "Upper bound")) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  labs(x = NULL,
       y = paste0("Change in Annual Cases Prevented (from baseline = ",
                  format(round(torn$baseline_prev[1]), big.mark = ","), ")"),
       fill = "Parameter\nValue") +
  theme_diss +
  theme(legend.position = if (packageVersion("ggplot2") >= "3.5.0") "inside" else c(0.85, 0.15),
        legend.background = element_rect(fill = "white", color = "grey80"))

# Set legend.position.inside if ggplot2 >= 3.5.0
if (packageVersion("ggplot2") >= "3.5.0") {
  fig8 <- fig8 + theme(legend.position.inside = c(0.85, 0.15))
}

ggsave(file.path(outdir, "Fig8_tornado_oneway.pdf"), fig8,
       width = 9, height = 5.5, device = cairo_pdf, bg = "white")


########################################################################
# FIGURE 9: LHS Forest Plot [was Fig 3]
########################################################################

cat("Figure 9: LHS forest plot...\n")

lhs_ui$strategy_short <- clean_strategy(lhs_ui$Strategy)
lhs_ui$strategy_short <- gsub("All Mtb-Infected \\(2%/yr\\)", "All Mtb-Infected", lhs_ui$strategy_short)
lhs_ui <- drop_nusb_only(lhs_ui, "strategy_short")
lhs_ui$strategy_short <- factor(lhs_ui$strategy_short, levels = rev(strat_order))

fig9 <- ggplot(lhs_ui, aes(x = strategy_short, y = Prevented_median,
                             color = strategy_short)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_errorbar(aes(ymin = Prevented_lo, ymax = Prevented_hi),
                width = 0.3, linewidth = 0.8, show.legend = FALSE) +
  coord_flip() +
  scale_color_manual(values = rev(strat_colors)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL,
       y = "Annual Cases Prevented (Median and 95% UI)") +
  theme_diss

ggsave(file.path(outdir, "Fig9_lhs_forest_plot.pdf"), fig9,
       width = 8, height = 4.5, device = cairo_pdf, bg = "white")


########################################################################
# FIGURE 10: PRCC Bar Plot (dual-scenario) [was Fig 4]
#   Two-panel: optimistic (LEFT, psi = 50%/yr) + plausible (RIGHT, psi = 5%/yr)
#   Matches optimistic-first reorder convention used throughout manuscript.
#   Falls back to single-panel (plausible) if TB_PRCC_RESULTS_OPTIMISTIC.csv
#   is not present in the working directory.
########################################################################

cat("Figure 10: PRCC bar plot (dual-scenario)...\n")

# ---- Load optimistic PRCC CSV if available -------------------------------
prcc_opt <- tryCatch(
  read.csv("TB_PRCC_RESULTS_OPTIMISTIC.csv", stringsAsFactors = FALSE),
  error = function(e) {
    cat("  Note: TB_PRCC_RESULTS_OPTIMISTIC.csv not found -- Fig 10 plausible-only.\n")
    NULL
  }
)

# ---- Panel builder (shared across scenarios) -----------------------------
build_prcc_panel <- function(prcc_df, panel_title) {
  pp <- prcc_df
  pp$parameter <- gsub("eps_s_FB", "eps_s_Non-U.S.-born", pp$parameter)
  pp$display <- ifelse(
    pp$parameter %in% names(clean_tornado_labels),
    clean_tornado_labels[pp$parameter],
    pp$parameter
  )
  pp <- pp[order(abs(pp$PRCC), decreasing = TRUE), ]
  pp$display <- factor(pp$display, levels = rev(pp$display))
  pp$sig <- ifelse(pp$p_value < 0.001, "***",
            ifelse(pp$p_value < 0.01, "**",
            ifelse(pp$p_value < 0.05, "*", "ns")))
  pp$fill_color <- ifelse(pp$PRCC > 0, "Positive", "Negative")

  ggplot(pp, aes(x = display, y = PRCC, fill = fill_color)) +
    geom_col(width = 0.65, show.legend = FALSE) +
    geom_text(aes(label = sig, y = PRCC + sign(PRCC) * 0.03),
              size = 3, color = "grey30") +
    coord_flip(ylim = c(-1, 1)) +
    scale_fill_manual(values = c("Positive" = "#2B6A99",
                                 "Negative" = "#F2A586")) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(title = panel_title,
         x = NULL, y = "Partial Rank Correlation Coefficient (PRCC)") +
    theme_diss +
    theme(plot.title = element_text(size = 11, face = "bold"))
}

# ---- Build panels --------------------------------------------------------
prcc_plaus_df <- prcc[prcc$outcome == "cases_prevented", ]
fig10_plaus <- build_prcc_panel(
  prcc_plaus_df,
  "Plausible (VE varies 30-70%, psi = 5%/yr)")

if (!is.null(prcc_opt)) {
  prcc_opt_df <- prcc_opt[prcc_opt$outcome == "cases_prevented", ]
  fig10_opt <- build_prcc_panel(
    prcc_opt_df,
    "Optimistic (VE varies 50-90%, psi = 50%/yr)")

  # Panel order: OPTIMISTIC LEFT, PLAUSIBLE RIGHT
  fig10 <- fig10_opt + fig10_plaus + plot_layout(ncol = 2)
  fig_width <- 14
} else {
  fig10 <- fig10_plaus
  fig_width <- 8
}

ggsave(file.path(outdir, "Fig10_prcc_barplot.pdf"), fig10,
       width = fig_width, height = 5, device = cairo_pdf, bg = "white")
ggsave(file.path(outdir, "Fig10_prcc_barplot.png"), fig10,
       width = fig_width, height = 5, dpi = 300, bg = "white")

cat(sprintf("  Saved: Fig10_prcc_barplot.pdf/png (%s)\n",
            if (!is.null(prcc_opt)) "dual-scenario, OPT LEFT / PLAUS RIGHT"
            else "plausible-only fallback"))


########################################################################
# FIGURE 11: Direct vs. Indirect Effects [was Fig 6]
########################################################################

cat("Figure 11: Direct vs. indirect effects...\n")

di$strategy_short <- clean_strategy(di$Strategy)
di <- drop_nusb_only(di, "strategy_short")
di$strategy_short <- factor(di$strategy_short, levels = rev(strat_order))

fig11_data <- data.frame(
  Strategy = rep(di$strategy_short, 2),
  Effect   = factor(rep(c("Direct", "Indirect"), each = nrow(di)),
                    levels = c("Indirect", "Direct")),
  Cases    = c(di$Direct_prevented, di$Indirect_prevented)
)

fig11 <- ggplot(fig11_data, aes(x = Strategy, y = Cases, fill = Effect)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Direct" = "#2B6A99", "Indirect" = "#F2A586")) +
  geom_text(data = di, aes(x = strategy_short, y = Total_prevented,
                            label = paste0(round(Pct_indirect, 1), "% indirect")),
            hjust = -0.05, size = 3, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::comma, limits = c(0, 1500)) +
  labs(x = NULL, y = "Annual Cases Prevented",
       fill = "Effect Type") +
  theme_diss +
  theme(legend.position = if (packageVersion("ggplot2") >= "3.5.0") "inside" else c(0.85, 0.15),
        legend.background = element_rect(fill = "white", color = "grey80"))

if (packageVersion("ggplot2") >= "3.5.0") {
  fig11 <- fig11 + theme(legend.position.inside = c(0.85, 0.15))
}

ggsave(file.path(outdir, "Fig11_direct_indirect.pdf"), fig11,
       width = 9, height = 4.5, device = cairo_pdf, bg = "white")


########################################################################
# FIGURE 12: Threshold Heatmap (VE x Duration) -- side-by-side two-panel:
#            optimistic (vaccination rate = 50%/yr, LEFT) vs.
#            plausible  (vaccination rate = 5%/yr,  RIGHT)
#            Both panels: All High-Risk strategy
########################################################################

cat("Figure 12: Threshold heatmap (optimistic + plausible)...\n")

# ---- Plausible panel data ----
thr_hr_plaus <- thresh[thresh$strategy == "All High-Risk", ]
# Restrict display to VE >= 30% (plausible licensing range); underlying
# analysis in TB_THRESHOLD_ANALYSIS.csv and Table 6 retain full 10-95% range.
thr_hr_plaus <- thr_hr_plaus[thr_hr_plaus$VE >= 30, ]
thr_hr_plaus$scenario <- "Plausible (vaccination rate = 5%/yr)"

# ---- Optimistic panel data (if available) ----
if (!is.null(thresh_opt)) {
  thr_hr_opt <- thresh_opt[thresh_opt$strategy == "All High-Risk", ]
  thr_hr_opt <- thr_hr_opt[thr_hr_opt$VE >= 30, ]
  thr_hr_opt$scenario <- "Optimistic (vaccination rate = 50%/yr)"
  thr_hr_all <- rbind(thr_hr_plaus, thr_hr_opt)
} else {
  thr_hr_all <- thr_hr_plaus
}

if (nrow(thr_hr_all) > 0) {

  thr_hr_all$dur_label <- paste0(thr_hr_all$duration_yr, " yr")
  thr_hr_all$dur_label <- factor(thr_hr_all$dur_label,
                                  levels = c("5 yr", "10 yr", "15 yr",
                                             "20 yr", "30 yr"))

  # Panel order: optimistic (left), plausible (right)
  thr_hr_all$scenario <- factor(thr_hr_all$scenario,
                                 levels = c("Optimistic (vaccination rate = 50%/yr)",
                                            "Plausible (vaccination rate = 5%/yr)"))

  # Shared color scale: max across both panels so colors are directly
  # comparable between plausible and optimistic
  shared_max <- max(thr_hr_all$pct_reduction, na.rm = TRUE)

  # Text color: white when cell is darker (pct_reduction above a threshold
  # scaled to the shared max), black otherwise
  text_cutoff <- 0.35 * shared_max

  fig12 <- ggplot(thr_hr_all,
                  aes(x = factor(VE), y = dur_label, fill = pct_reduction)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", pct_reduction)),
              size = 2.6,
              color = ifelse(thr_hr_all$pct_reduction > text_cutoff,
                              "white", "black")) +
    scale_fill_gradientn(
      colors = c("#FFFFFF", "#F2A586", "#C44E6C", "#2B6A99", "#1A4A6E"),
      values = scales::rescale(c(0, 0.30, 0.55, 0.85, 1.0)),
      name = "% Reduction",
      limits = c(0, shared_max),
      # Tall, narrow colorbar so the gradient is easy to read at docx page width.
      guide = guide_colorbar(barheight = unit(14, "cm"),
                             barwidth  = unit(0.5, "cm"),
                             ticks.colour = "grey40",
                             frame.colour = "grey40")
    ) +
    # Panels stacked vertically (optimistic top, plausible below) so cells are
    # large enough to render cell labels legibly at column width in the docx.
    facet_wrap(~ scenario, ncol = 1) +
    labs(x = "Vaccine Efficacy (%)", y = "Duration of Protection") +
    theme_diss +
    theme(panel.grid  = element_blank(),
          axis.ticks  = element_blank(),
          strip.text  = element_text(face = "bold", size = 12),
          strip.background = element_rect(fill = "grey95", color = NA),
          legend.position  = "right",
          legend.title     = element_text(size = 11),
          legend.text      = element_text(size = 10))

  # Vertical stacking: taller rather than wider. If only plausible is available,
  # fall back to the single-panel wide layout.
  if (!is.null(thresh_opt)) {
    fig_width  <- 10
    fig_height <- 8.5
  } else {
    fig_width  <- 10
    fig_height <- 4.5
  }
  ggsave(file.path(outdir, "Fig12_threshold_heatmap.pdf"), fig12,
         width = fig_width, height = fig_height, device = cairo_pdf, bg = "white")
  ggsave(file.path(outdir, "Fig12_threshold_heatmap.png"), fig12,
         width = fig_width, height = fig_height, dpi = 300, bg = "white")

  cat(sprintf("  Fig 12 saved (%s panel%s)\n",
              if (!is.null(thresh_opt)) "two" else "one",
              if (!is.null(thresh_opt)) "s" else ""))
} else {
  cat("  WARNING: No All High-Risk data in threshold. Skipping.\n")
}




########################################################################
# FIGURE 13: NNV Comparison with Established Vaccines + Harris 2019
#
# v9: two-scenario TB entries (6 strategies x 2 scenarios = 12 TB rows).
# Optimistic entries: dark-blue triangles, UI from LHS_OPTIMISTIC.
# Plausible  entries: light-blue circles,  UI from LHS (plausible).
#
# Fixes v8 missing-whisker bugs:
#   (1) "Medical comorbidities" in the optimistic UI file now normalizes to
#       "Medical" via clean_strategy() (see update above), so Medical
#       (optimistic) merges correctly and gets whiskers.
#   (2) Optimistic UI is loaded from its own file.
#
# KNOWN DATA GAP: the plausible UI file (TB_LHS_UNCERTAINTY_INTERVALS.csv)
# does not contain a "PLWH + Non-U.S.-Born" row -- it has "All Non-U.S.-Born"
# instead. That plausible TB entry therefore has no LHS UI whiskers until the
# plausible LHS analysis is rerun with the current strategy list. Every other
# TB entry will render with full UI.
########################################################################

cat("Figure 13 (v9): NNV comparison, two-scenario...\n")

if (!require("ggtext", quietly = TRUE)) {
  install.packages("ggtext", repos = "https://cloud.r-project.org", quiet = TRUE)
  library(ggtext)
}

# Load optimistic UI (required for optimistic whiskers)
lhs_ui_opt <- tryCatch(
  read.csv("TB_LHS_UNCERTAINTY_INTERVALS_OPTIMISTIC.csv"),
  error = function(e) {
    cat("  WARNING: TB_LHS_UNCERTAINTY_INTERVALS_OPTIMISTIC.csv not found;",
        "optimistic TB whiskers will be absent.\n")
    NULL
  }
)

# --- Plausible TB points (from strat_sum) + plausible UI ---
plaus_pts <- strat_sum[strat_sum$Scenario != "Baseline", ]
plaus_pts$Strategy <- clean_strategy(plaus_pts$Scenario)
plaus_pts <- drop_nusb_only(plaus_pts, "Strategy")

lhs_plaus <- lhs_ui
lhs_plaus$Strategy <- clean_strategy(lhs_plaus$Strategy)
lhs_plaus <- drop_nusb_only(lhs_plaus, "Strategy")

plaus_merge <- merge(
  plaus_pts[, c("Strategy", "NNV")],
  lhs_plaus[, c("Strategy", "NNV_lo", "NNV_hi")],
  by = "Strategy", all.x = TRUE
)
plaus_merge$NNV_lo[is.na(plaus_merge$NNV_lo)] <- plaus_merge$NNV[is.na(plaus_merge$NNV_lo)]
plaus_merge$NNV_hi[is.na(plaus_merge$NNV_hi)] <- plaus_merge$NNV[is.na(plaus_merge$NNV_hi)]

plaus_df <- data.frame(
  Label    = paste0("**TB**: ", plaus_merge$Strategy, " (plausible)"),
  NNV_mid  = plaus_merge$NNV,
  NNV_lo   = plaus_merge$NNV_lo,
  NNV_hi   = plaus_merge$NNV_hi,
  Category = "TB Vaccination (This Study, Plausible)",
  stringsAsFactors = FALSE
)

# --- Optimistic TB points (from opt_strat at VE=70, psi=50) + optimistic UI ---
if (!is.null(opt_strat)) {
  opt_pts <- opt_strat[opt_strat$VE_pct == 70 & opt_strat$psi_pct == 50, ]
  opt_pts$Strategy <- clean_strategy(opt_pts$strategy)
  opt_pts <- drop_nusb_only(opt_pts, "Strategy")

  if (!is.null(lhs_ui_opt)) {
    lhs_opt_m <- lhs_ui_opt
    lhs_opt_m$Strategy <- clean_strategy(lhs_opt_m$Strategy)
    lhs_opt_m <- drop_nusb_only(lhs_opt_m, "Strategy")
    opt_merge <- merge(
      opt_pts[, c("Strategy", "NNV")],
      lhs_opt_m[, c("Strategy", "NNV_lo", "NNV_hi")],
      by = "Strategy", all.x = TRUE
    )
  } else {
    opt_merge <- opt_pts[, c("Strategy", "NNV")]
    opt_merge$NNV_lo <- NA_real_
    opt_merge$NNV_hi <- NA_real_
  }
  opt_merge$NNV_lo[is.na(opt_merge$NNV_lo)] <- opt_merge$NNV[is.na(opt_merge$NNV_lo)]
  opt_merge$NNV_hi[is.na(opt_merge$NNV_hi)] <- opt_merge$NNV[is.na(opt_merge$NNV_hi)]

  opt_df <- data.frame(
    Label    = paste0("**TB**: ", opt_merge$Strategy, " (optimistic)"),
    NNV_mid  = opt_merge$NNV,
    NNV_lo   = opt_merge$NNV_lo,
    NNV_hi   = opt_merge$NNV_hi,
    Category = "TB Vaccination (This Study, Optimistic)",
    stringsAsFactors = FALSE
  )
} else {
  opt_df <- plaus_df[0, ]
}

our_df <- rbind(opt_df, plaus_df)

# Harris 2019 (Lancet Glob Health)
harris_df <- data.frame(
  Label = c(
    "**TB**: Pre-infection (China, older adult)",
    "**TB**: Post-infection latency (China, older adult)",
    "**TB**: Post-infection L&R (China, older adult)",
    "**TB**: Pre- & post-infection (China, older adult)"
  ),
  NNV_mid = c(1022, 574, 292, 230),
  NNV_lo  = c(752,  350, 257, 199),
  NNV_hi  = c(1318, 2886, 365, 269),
  Category = "TB Modeling (Published)",
  stringsAsFactors = FALSE
)

# Established vaccine comparators
comp_df <- data.frame(
  Label = c(
    "**Shingrix**: HZ infection (Adults \u226550 yr)",
    "**Influenza**: Illness (Adults \u226565 yr, Australia)",
    "**Influenza**: Hospitalization (Adults \u226565 yr, Australia)",
    "**HPV**: Cervical cancer (Girls age 12, Canada)",
    "**HPV**: Cervical cancer death (Girls age 12, Canada)",
    "**PCV-13**: All CAP, 5 yr (Adults \u226565 yr, US)",
    "**PCV-13**: Hospital CAP, 5 yr (Adults \u226565 yr, US)",
    "**PCV-13**: All CAP, 1 yr (Adults \u226565 yr, US)",
    "**PCV-13**: Hospital CAP, 1 yr (Adults \u226565 yr, US)",
    "**Influenza**: Death (Adults \u226565 yr, Australia)"
  ),
  NNV_mid = c(11, 43, 777, 324, 729, 234, 576, 656, 1620, 5338),
  NNV_lo  = c(8,  16, 470, 195, 411, NA,  NA,  454, 1110, 1429),
  NNV_hi  = c(15, 192, 1684, 757, 1921, NA, NA, 2110, 5130, 19231),
  Category = "Other Recommended Vaccines",
  stringsAsFactors = FALSE
)
comp_df$NNV_lo[is.na(comp_df$NNV_lo)] <- comp_df$NNV_mid[is.na(comp_df$NNV_lo)]
comp_df$NNV_hi[is.na(comp_df$NNV_hi)] <- comp_df$NNV_mid[is.na(comp_df$NNV_hi)]

bench_df <- rbind(our_df, harris_df, comp_df)

# Order: optimistic TB at top, then plausible TB, then Harris, then comparators
# (each group internally sorted by NNV). coord_flip() inverts axis, so reverse.
opt_order    <- opt_df$Label[order(opt_df$NNV_mid)]
plaus_order  <- plaus_df$Label[order(plaus_df$NNV_mid)]
harris_order <- harris_df$Label[order(harris_df$NNV_mid)]
comp_order   <- comp_df$Label[order(comp_df$NNV_mid)]

bench_df$Label <- factor(bench_df$Label,
                         levels = c(comp_order, harris_order,
                                    plaus_order, opt_order))

cat_colors <- c(
  "TB Vaccination (This Study, Optimistic)" = "#1F4E79",
  "TB Vaccination (This Study, Plausible)"  = "#6FA8DC",
  "TB Modeling (Published)"                 = "#D4A84B",
  "Other Recommended Vaccines"              = "#888888"
)
cat_shapes <- c(
  "TB Vaccination (This Study, Optimistic)" = 17,
  "TB Vaccination (This Study, Plausible)"  = 16,
  "TB Modeling (Published)"                 = 16,
  "Other Recommended Vaccines"              = 16
)

fig13 <- ggplot(bench_df, aes(x = Label, y = NNV_mid,
                              color = Category, shape = Category)) +
  geom_errorbar(aes(ymin = NNV_lo, ymax = NNV_hi),
                width = 0.3, linewidth = 0.7) +
  geom_point(size = 3) +
  coord_flip() +
  scale_color_manual(values = cat_colors, name = NULL) +
  scale_shape_manual(values = cat_shapes, name = NULL) +
  scale_y_log10(
    labels = scales::comma,
    breaks = c(10, 50, 100, 200, 500, 1000, 2000, 5000, 20000)
  ) +
  labs(x = NULL,
       y = "Number Needed to Vaccinate (NNV, log scale)") +
  theme_diss +
  theme(legend.position = "top",
        axis.text.y     = ggtext::element_markdown(size = 8)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE))

ggsave(file.path(outdir, "Fig13_nnv_benchmarking.pdf"), fig13,
       width = 11, height = 9, device = cairo_pdf, bg = "white")
ggsave(file.path(outdir, "Fig13_nnv_benchmarking.png"), fig13,
       width = 11, height = 9, dpi = 300, bg = "white")

cat("  Saved: Fig13_nnv_benchmarking.pdf/png\n")


########################################################################
# FIGURE 14: Psi Dose-Response [was Fig 8]
########################################################################

cat("Figure 14: Psi dose-response...\n")

psi_plot <- psi_s
psi_plot$strategy_short <- clean_strategy(psi_plot$strategy)
psi_plot <- drop_nusb_only(psi_plot, "strategy_short")
psi_plot$strategy_short <- factor(psi_plot$strategy_short, levels = strat_order)
psi_plot <- psi_plot[!is.na(psi_plot$strategy_short), ]

fig14 <- ggplot(psi_plot, aes(x = psi_pct, y = NNV, color = strategy_short,
                               shape = strategy_short)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  scale_color_manual(values = strat_colors, name = "Strategy") +
  scale_shape_manual(values = strat_shapes, name = "Strategy") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Vaccination Rate (psi, %/year)",
       y = "Number Needed to Vaccinate (NNV)") +
  theme_diss +
  theme(legend.position = "right")

ggsave(file.path(outdir, "Fig14_psi_doseresponse.pdf"), fig14,
       width = 9, height = 5, device = cairo_pdf, bg = "white")


########################################################################
# SUMMARY
########################################################################

cat("\n================================================================\n")
cat("  ALL FIGURES GENERATED (v2)\n")
cat("================================================================\n")
cat("  Output directory:", outdir, "\n\n")
cat("  Fig 1:  Fig1_model_diagram.png         (pre-made, copied)\n")
cat("  Fig 2:  Fig2_case_source_validation.pdf  (case source validation)\n")
cat("  Fig S:  Fig_S_recent_transmission.pdf    (appendix: recent transmission)\n")
cat("  Fig 3:  Fig3_strategy_comparison_grouped.pdf  (grouped bars, optimistic + plausible)\n")
cat("  Fig 4:  Fig4_stratum_impact.pdf         ** NEW **\n")
cat("  Fig 5:  Fig5_incidence_trajectory.pdf   ** NEW **\n")
cat("  Fig 6:  Fig6_cumulative_prevented.pdf   ** NEW **\n")
cat("  Fig 7:  Fig7_time_to_impact.pdf         (updated from old Fig 5)\n")
cat("  Fig 8:  Fig8_tornado_oneway.pdf         (updated from old Fig 2)\n")
cat("  Fig 9:  Fig9_lhs_forest_plot.pdf        (updated from old Fig 3)\n")
cat("  Fig 10: Fig10_prcc_barplot.pdf          (updated from old Fig 4)\n")
cat("  Fig 11: Fig11_direct_indirect.pdf       (updated from old Fig 6)\n")
cat("  Fig 12: Fig12_threshold_heatmap.pdf     (updated from old Fig 7)\n")
cat("  Fig 13: Fig13_nnv_benchmarking.pdf      ** UPDATED: established vaccine comparison **\n")
cat("  Fig 14: Fig14_psi_doseresponse.pdf      (updated from old Fig 8)\n")
cat("================================================================\n")
