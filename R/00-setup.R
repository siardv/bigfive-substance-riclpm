# 00-setup.R
# configuration and income category bounds. sourced first; defines objects only.

config <- list(
  merged_dir   = "data_merged",     # frozen per-module merges read by load_liss_data
  avars_file   = "data_init/b.sav", # background variables (not fetched by lissr)
  raw_dir      = "data",            # per-wave source files, used only to rebuild merges
  liss_modules = c("cp", "ch", "ci")
)

# income category bounds (LISS coding 1-7)
income_bounds <- list(
  lower = c(0, 8000, 16000, 24000, 36000, 48000, 60000),
  upper = c(8000, 16000, 24000, 36000, 48000, 60000, 120000)
)
