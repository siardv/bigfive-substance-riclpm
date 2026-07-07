# full_record.R
# one-shot capture of every output the manuscript patch slots consume
source("run.R")

# main run: reproduction check, ses stability (P6), attrition (P10), ci bounds (P4)
results <- run_riclpm_analysis(config, waves = 5)

# substance battery availability by year (P3)
describe_substance_items(results$liss$h)

# longitudinal invariance on the selected occasions (P8)
inv <- run_longitudinal_invariance(
  results$liss$p, results$weasel_selection,
  dims = c("extr", "open", "cons")
)

# sensitivity battery: illicit_five (P5), lenient (P2), annual_interval (P9),
# ses_equivalized (P6), parcels (P1), plus the binary and tertile robustness
sens <- run_sensitivity_analysis(config, waves = 5)