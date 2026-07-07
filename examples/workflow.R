# examples/workflow.R
# usage examples extracted verbatim from the original script tail.
# NOT sourced by run.R: it contains executable example calls that would
# otherwise run the full analysis on load. run these interactively.

# =============================================================================
# USAGE
# =============================================================================

# # ============================================================================
# # RECOMMENDED WORKFLOW (Quick Start)
# # ============================================================================
#
# # full pipeline with optimal settings
# results <- run_riclpm_analysis(
#   config,
#   waves = 5,
#   ses_method = "binary", # balanced groups (~50/50) - CRITICAL for H3
#   substance_method = "sum", # count variable (0-6) for adequate variance
#   estimator = "MLR", # robust standard errors
#   check_coverage = TRUE, # run coverage diagnostics
#   use_weasel = TRUE, # principled subset selection (requires library(weasel))
#   weasel_measured = "extr" # align selection with personality measurement waves
# )
#
# the weasel methods paragraph and audit tables live in:
# results$weasel_selection$justification
# results$weasel_selection$comparison
#
# create publication-ready table
# pub_table <- create_publication_table(results)
# print_publication_table(pub_table)
#
# ============================================================================
# IMPORTANT METHODOLOGICAL NOTES
# ============================================================================
#
# 1. SES METHOD MATTERS:
#    - "binary" = median split (~50/50) - RECOMMENDED for reliable moderation tests
#    - "composite" = 3 levels (low/middle/high) - often creates severe imbalance
#    - If using 3-level SES, expect multi-group χ² tests to be unreliable
#
# 2. SUBSTANCE USE CODING MATTERS:
#    - "sum" = count 0-6 - RECOMMENDED for variance
#    - "any" = binary 0/1 - creates zero variance in subgroups (93% zeros)
#    - If using "any", expect model failures in multi-group analysis
#
# 3. INTERPRETING H3 (SES MODERATION):
#    - Check "Reliable" column in publication table
#    - If "No*", results are unreliable due to group imbalance
#    - Imbalance > 5:1 = interpret with extreme caution
#    - Consider interaction-based moderation as alternative
#
# ============================================================================
# DETAILED STEP-BY-STEP WORKFLOW
# ============================================================================
#
# 1. Load data
# load_liss_data() reads the frozen per-module merges in config$merged_dir
# (data_merged/{cp,ch,ci}_merged.sav) plus config$avars_file; no credentials or
# 2fa at runtime. to rebuild the merges from raw per-wave files (off the analysis
# path) run fetch_liss_raw(config) then merge_liss_to_disk(config) once.
liss <- load_liss_data(config)
#
# 2. Clean income with outlier detection
income_clean <- clean_income(liss$i, liss$b,
  income_cap = 150000,
  dataset_outliers = TRUE
)
summary <- summarize_outliers(income_clean)
print_outlier_summary(summary)
#
# 3. Construct SES - BINARY RECOMMENDED for multi-group analysis
income_clean <- construct_ses(income_clean, method = "binary")
table(income_clean$ses) # should be ~50/50 split
#
# 4. Big Five scores
b5_scores <- compute_b5_scores(liss$p)
#
# 5. Substance use - SUM (count) RECOMMENDED
# NOTE: "any" (binary) creates zero variance in subgroups due to low base rate (~7%)
substance <- process_substance_use(liss$h, method = "sum")
table(substance$substance_use)
#
# 6. Prepare analysis data
liss$i_clean <- income_clean
analysis_data <- prepare_analysis_data(liss, substance_method = "sum")
#
# 6b. weasel subset selection (showcase) - principled, reproducible, citable.
# restricts to waves where extr is measured, re-indexes to a dense per-person
# axis, scores scenarios, and emits a methods paragraph (sel$justification).
# requires library(weasel) or a sourced weasel_all_in_one.R on disk.
sel <- select_weasel_subset(analysis_data, measured = "extr",
  span = "core", core_len = 5L,
  author = "van den Bosch", year = "2026"
)
cat(sel$justification, "\n") # paste-ready methods text
#
# 7. Reshape for RI-CLPM (weasel's window length drives the wave count)
data_wide <- reshape_for_riclpm(sel$data, max_waves = sel$L)
#
# 8. CHECK VARIANCE BEFORE FITTING - CRITICAL!
sapply(paste0("substance_use_", 1:5), function(v) var(data_wide[[v]], na.rm = TRUE))
# Values should be > 0.01; if any are 0 or near-0, use method = "sum"
#
# 9. Check coverage
coverage <- check_coverage(data_wide,
  vars = c("extr", "open", "cons", "substance_use"),
  threshold = 0.10
)
print_coverage_summary(coverage)
#
# 10. Check variance by SES group (for multi-group models)
var_check <- check_variance_by_group(data_wide, "ses", "substance_use")
if (var_check$has_zero_variance) {
  message("WARNING: Zero variance in some SES groups - multi-group model will fail")
  message("Consider using interaction model instead")
}
#
# 11. Fit single RI-CLPM with robust estimator
riclpm_extr <- fit_riclpm(data_wide, "extr", "substance_use",
  waves = 5, estimator = "MLR"
)
lavaan::summary(riclpm_extr$fit, fit.measures = TRUE, standardized = TRUE)
#
# 12. Test all hypotheses
results <- run_hypothesis_tests(data_wide, waves = 5, estimator = "MLR")
print_hypothesis_summary(results)
#
# ============================================================================
# CFA OPTIONS FOR POOR-FITTING MEASUREMENT MODELS
# ============================================================================
#
# Standard CFA (may have poor fit for Big Five)
cfa_standard <- fit_b5_cfa(liss$p, "open", model_type = "standard")
#
# Bifactor CFA (often better fit - allows cross-loadings)
cfa_bifactor <- fit_b5_cfa(liss$p, "open", model_type = "bifactor")
#
# Parcel-based CFA (reduces complexity)
cfa_parcel <- fit_b5_cfa(liss$p, "open", model_type = "parcel")
#
# Compare all models for a dimension
cfa_comparison <- compare_cfa_models(liss$p, "open")
print(cfa_comparison$comparison)
#
# ============================================================================
# SENSITIVITY ANALYSIS
# ============================================================================
#
# Compare different specifications
sensitivity <- run_sensitivity_analysis(config, waves = 5)
#
# ============================================================================
# INTERPRETING NULL FINDINGS
# ============================================================================
#
# If all cross-lagged paths are non-significant, check:
# 1. Random Intercept loadings - high values (>0.80) indicate trait-like stability
#    dominates, leaving little within-person variance for cross-lagged effects
# 2. This is a VALID finding, not a failure - it means personality-substance use
#    associations are due to stable between-person differences, not dynamic
#    reciprocal causation within persons over time
# 3. Report the RI loadings alongside the null cross-lagged effects