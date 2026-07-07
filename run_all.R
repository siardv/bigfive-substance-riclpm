#!/usr/bin/env Rscript
# render the RI-CLPM analysis (Big Five, substance use, SES) as a
# github-renderable transcript: run_all.Rmd -> run_all.md.
# run from the project root, next to run.R.

# ---- 0. packages ----
# cran packages used by the R/ modules plus the render toolchain.
pkgs <- c(
  "dplyr", "haven", "imputeTS", "lavaan", "magrittr", "purrr",
  "sjlabelled", "tibble", "tidyr", "knitr", "rmarkdown"
)
missing <- setdiff(pkgs, rownames(utils::installed.packages()))
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

# weasel powers module 06 (subset selection) and is not on cran: it must be
# attached (library(weasel)). when it is absent the transcript still renders,
# with the subset-selection step marked skipped.
weasel_ready <- requireNamespace("weasel", quietly = TRUE)
if (!weasel_ready) {
  message(
    "note: weasel not found. the subset-selection step will render as skipped.\n",
    "install weasel to include it."
  )
}

# lissr is only needed to rebuild the frozen merges (merge_liss_to_disk). this
# transcript reads data_merged/*.sav directly via haven and does not require it.

# ---- 1. render the transcript ----
if (!rmarkdown::pandoc_available()) {
  stop("pandoc not found: install pandoc or render run_all.Rmd inside RStudio")
}
rmarkdown::render("run_all.Rmd", envir = new.env(), quiet = TRUE)
cat("wrote run_all.md (figures under run_all_files/, tables under tables/)\n")
