# run.R
# entry point: sources the R/ modules in numeric order and exposes the pipeline.
# defines objects only; runs no analysis on load.
#
# packages used via :: are haven, dplyr, purrr, tibble, lavaan, tidyr,
# sjlabelled, imputeTS, magrittr and lissr. weasel must be attached
# (library(weasel)) for module 06.

.riclpm_modules <- sort(
  list.files("R", pattern = "^[0-9].*\\.R$", full.names = TRUE)
)
invisible(lapply(.riclpm_modules, source))
rm(.riclpm_modules)

message("riclpm modules loaded. entry points:")
message("  run_riclpm_analysis(config, ...)   # full pipeline (weasel on by default)")
message("  run_sensitivity_analysis(config)   # specification comparison")
message("see examples/workflow.R for an interactive step-by-step walkthrough,")
message("or run 'Rscript run_all.R' to render the full analysis transcript (run_all.md).")
