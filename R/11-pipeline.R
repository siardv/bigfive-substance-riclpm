# 11-pipeline.R
# top-level orchestration: run_riclpm_analysis() and run_sensitivity_analysis().
# sourced last; ties the modules together.

#' run complete RI-CLPM analysis
#' @param config configuration list
#' @param waves number of waves for RI-CLPM (also the weasel target window length)
#' @param ses_method SES coding: "binary" (recommended), "composite", "continuous"
#' @param substance_method substance use coding: "sum" (count, RECOMMENDED) or "any" (binary)
#' @param estimator "MLR" (robust, recommended) or "ML"
#' @param check_coverage whether to run coverage diagnostics
#' @param use_weasel run weasel subset selection between prepare and reshape
#' @param weasel_measured variable(s) that must be observed for a wave to count
#'   as present (default "extr"); aligns selection with personality measurement
#' @param weasel_span weasel span rule: "core" or "full"
#' @param weasel_scenario weasel scenario to apply; NULL uses the recommended
#'   one; "lenient_info_max" gives the interior-gap-tolerant selection
#' @param traits personality traits to model (extend with "neur" and "agre"
#'   for the full-domain sensitivity)
#' @param substance_items substance item variables (default: full six-item
#'   battery; pass illicit_five_vars for the sedatives-excluded sensitivity)
#' @param substance_min_valid minimum answered items for a non-NA substance score
#' @param ses_equivalize equivalize household income for household composition
#'   before the binary split (OECD-modified scale when counts are available)
#' @param ses_split_within_year classify each person-wave against its own
#'   calendar year's income median rather than the pooled median
#' @param interval_filter "none" (default) or "annual": restrict to respondents
#'   whose selected window is spaced at exactly one year throughout
#' @param measurement "composite" (single-indicator, default) or "parcels"
#'   (multiple-indicator RI-CLPM correcting within-person estimates for
#'   parcel unreliability; H3 stays on the composite specification)
#' @param report_attrition compare selected respondents with the excluded
#'   remainder of the eligible pool (requires weasel)
#' @param run_invariance fit the item-level longitudinal invariance hierarchy
#'   on the selected occasions (requires weasel)
#' @note Use substance_method="sum" (not "any") to ensure adequate variance
run_riclpm_analysis <- function(config = config, waves = 5,
                                ses_method = "binary",
                                substance_method = "sum",
                                estimator = "MLR",
                                check_coverage = TRUE,
                                use_weasel = TRUE,
                                weasel_measured = "extr",
                                weasel_span = "core",
                                weasel_scenario = NULL,
                                traits = c("extr", "open", "cons"),
                                substance_items = drugs_use_vars,
                                substance_min_valid = 1,
                                ses_equivalize = FALSE,
                                ses_split_within_year = FALSE,
                                interval_filter = c("none", "annual"),
                                measurement = c("composite", "parcels"),
                                report_attrition = TRUE,
                                run_invariance = FALSE) {
  interval_filter <- match.arg(interval_filter)
  measurement <- match.arg(measurement)

  message("=", strrep("=", 59))
  message("RI-CLPM ANALYSIS: Big Five, Substance Use, and SES")
  message(strrep("=", 60), "\n")
  message(
    "Settings: waves=", waves, ", SES=", ses_method,
    ", substance=", substance_method, ", estimator=", estimator,
    ", weasel=", use_weasel,
    ", traits=", paste(traits, collapse = "/"),
    ", items=", length(substance_items),
    ", measurement=", measurement,
    ", interval=", interval_filter, "\n"
  )

  # load and prepare data
  message("Step 1: Loading data...")
  liss <- load_liss_data(config)

  message("\nStep 2: Cleaning income and constructing SES...")
  liss$i_clean <- clean_income(liss$i, liss$b, verbose = FALSE)
  liss$i_clean <- construct_ses(liss$i_clean,
    method = ses_method,
    equivalize = ses_equivalize,
    split_within_year = ses_split_within_year
  )

  # verify SES levels
  ses_tab <- table(liss$i_clean$ses, useNA = "ifany")
  message("  SES after construct_ses: ", paste(names(ses_tab), "=", ses_tab, collapse = ", "))
  if (ses_method == "binary") {
    describe_ses_stability(liss$i_clean)
  }

  message("\nStep 3: Preparing analysis dataset...")
  # pass substance settings and the parcel flag to prepare_analysis_data
  analysis_data <- prepare_analysis_data(liss,
    substance_method = substance_method,
    substance_items = substance_items,
    substance_min_valid = substance_min_valid,
    include_parcels = (measurement == "parcels")
  )

  # verify SES still correct after merge
  ses_check <- table(analysis_data$ses, useNA = "ifany")
  message("  SES after prepare_analysis_data: ", paste(names(ses_check), "=", ses_check, collapse = ", "))

  # keep the pre-selection frame for the attrition comparison
  analysis_data_pre <- analysis_data

  # Step 3b: principled subset selection with weasel (showcase). operates on the
  # long analysis data, restricts to waves where weasel_measured is observed,
  # re-indexes to a dense per-person axis, and selects a scenario. the chosen
  # window length then drives the reshape, so the RI-CLPM wave count is rule-based
  weasel_selection <- NULL
  if (isTRUE(use_weasel)) {
    message("\nStep 3b: Selecting analysis subset with weasel...")
    weasel_selection <- select_weasel_subset(
      analysis_data,
      measured = weasel_measured,
      span = weasel_span,
      core_len = waves,
      scenario = weasel_scenario
    )
    analysis_data <- weasel_selection$data
    waves <- weasel_selection$L
  }

  # Step 3c: optional uniform-interval restriction; keeps respondents whose
  # consecutive selected waves are all one calendar year apart, so the lag-1
  # coefficients carry a single, well-defined interval
  uniform_interval <- NULL
  if (interval_filter == "annual" && !is.null(weasel_selection)) {
    message("\nStep 3c: Restricting to uniform one-year intervals...")
    uniform_interval <- filter_uniform_interval(weasel_selection)
    if (!is.null(uniform_interval)) {
      analysis_data <- dplyr::filter(
        analysis_data,
        nomem_encr %in% uniform_interval$ids
      )
    }
  }

  # attrition: selected respondents against the excluded remainder of the
  # eligible pool (anyone with an observed anchoring variable)
  attrition <- NULL
  if (isTRUE(report_attrition) && !is.null(weasel_selection)) {
    message("\nStep 3d: Selection attrition comparison...")
    attrition <- describe_selection_attrition(
      analysis_data_pre,
      selected_ids = unique(analysis_data$nomem_encr),
      anchor = weasel_measured[1],
      traits = traits
    )
  }

  message("\nStep 4: Reshaping for RI-CLPM...")
  reshape_vars <- unique(c(traits, "substance_use"))
  if (measurement == "parcels") {
    reshape_vars <- unique(c(reshape_vars, parcel_var_names(traits, 3)))
  }
  data_wide <- reshape_for_riclpm(analysis_data,
    vars = reshape_vars,
    max_waves = waves
  )

  # descriptive statistics
  message("\nDescriptive statistics:")
  message("  N persons: ", nrow(data_wide))
  message("  SES distribution:")
  ses_wide <- table(data_wide$ses, useNA = "ifany")
  print(ses_wide)

  # verify binary SES has only 2 levels
  if (ses_method == "binary" && length(names(ses_wide)[!is.na(names(ses_wide))]) != 2) {
    warning(
      "SES should have 2 levels for binary method, but has: ",
      paste(names(ses_wide), collapse = ", ")
    )
  }

  # check substance use variance - critical for model convergence
  sub_cols <- grep("^substance_use_", names(data_wide), value = TRUE)
  if (length(sub_cols) > 0) {
    sub_means <- sapply(sub_cols, function(col) mean(data_wide[[col]], na.rm = TRUE))
    sub_vars <- sapply(sub_cols, function(col) stats::var(data_wide[[col]], na.rm = TRUE))
    message(
      "  Substance use mean per wave: ",
      paste(round(sub_means, 3), collapse = ", ")
    )
    message(
      "  Substance use var per wave: ",
      paste(round(sub_vars, 3), collapse = ", ")
    )

    # warn if variance is too low
    if (any(sub_vars < 0.01, na.rm = TRUE)) {
      message("\n  WARNING: Low variance in substance_use. Consider using method='sum' instead of 'any'")
    }
  }

  # coverage diagnostics
  coverage <- NULL
  if (check_coverage) {
    message("\nStep 4b: Checking pairwise coverage...")
    coverage <- check_coverage(data_wide,
      vars = unique(c(traits, "substance_use")),
      threshold = 0.10
    )
    print_coverage_summary(coverage)
  }

  message("\nStep 5: Running CFA for Big Five dimensions...")
  cfa_results <- purrr::map(traits, function(dim) {
    fit_b5_cfa(liss$p, dim)
  })
  cfa_fit_table <- extract_cfa_fit(cfa_results)
  print(cfa_fit_table)

  # warn if CFA fit is poor and compare alternative models
  poor_cfa <- cfa_fit_table$cfi < 0.90
  cfa_comparisons <- NULL

  if (any(poor_cfa, na.rm = TRUE)) {
    poor_dims <- cfa_fit_table$dimension[poor_cfa]
    message(
      "\n⚠️  WARNING: Poor CFA fit (CFI < 0.90) for: ",
      paste(poor_dims, collapse = ", ")
    )
    message("\nComparing alternative measurement models...")

    # try alternative models for poor-fitting dimensions
    cfa_comparisons <- list()
    for (dim in poor_dims) {
      message("  Comparing models for ", dim, "...")
      comparison <- tryCatch(
        {
          compare_cfa_models(liss$p, dim)
        },
        error = function(e) {
          message("    comparison failed: ", e$message)
          NULL
        }
      )

      if (!is.null(comparison)) {
        cfa_comparisons[[dim]] <- comparison

        # find best model
        best_idx <- which.max(comparison$comparison$cfi)
        best_model <- comparison$comparison$model_type[best_idx]
        best_cfi <- comparison$comparison$cfi[best_idx]

        message("    Best model: ", best_model, " (CFI=", round(best_cfi, 3), ")")

        if (best_model == "parcel" && best_cfi >= 0.95) {
          # NOTE: the three-parcel single-factor CFA is just-identified, so its
          # fit indices are 1.0 by construction and carry no evidential weight;
          # the substantive parcel check is measurement="parcels", which
          # re-estimates the structural model with parcel indicators
          message(
            "    NOTE: parcel CFA is just-identified (fit indices are ",
            "uninformative); use measurement='parcels' for the structural check"
          )
        }
      }
    }

    message("\nIMPLICATIONS OF POOR CFA FIT:")
    message("  - Cross-lagged effects may be ATTENUATED (biased toward zero)")
    message("  - True effects might exist but not be detectable")
    message("  - This is a common issue with Big Five measurement in panel data")
    message("  - The parcel measurement model (measurement='parcels') tests this directly")
  }

  # optional item-level longitudinal invariance on the selected occasions
  invariance <- NULL
  if (isTRUE(run_invariance) && !is.null(weasel_selection)) {
    message("\nStep 5b: Longitudinal measurement invariance...")
    invariance <- run_longitudinal_invariance(
      liss$p, weasel_selection,
      dims = traits, estimator = estimator
    )
  }

  message("\nStep 6: Testing hypotheses...")
  hypothesis_results <- run_hypothesis_tests(data_wide, waves,
    estimator = estimator,
    traits = traits,
    measurement = measurement
  )

  message("\nStep 7: Generating summary and interpretation...")
  print_hypothesis_summary(hypothesis_results)

  # add interpretation of null findings
  interpret_riclpm_results(hypothesis_results, cfa_fit_table, coverage)

  invisible(list(
    liss = liss,
    analysis_data = analysis_data,
    data_wide = data_wide,
    weasel_selection = weasel_selection,
    uniform_interval = uniform_interval,
    attrition = attrition,
    cfa_results = cfa_results,
    cfa_fit = cfa_fit_table,
    cfa_comparisons = cfa_comparisons,
    invariance = invariance,
    hypothesis_results = hypothesis_results,
    coverage = coverage,
    settings = list(
      waves = waves, ses_method = ses_method,
      substance_method = substance_method, estimator = estimator,
      use_weasel = use_weasel,
      weasel_scenario = if (!is.null(weasel_selection)) weasel_selection$scenario else NA_character_,
      traits = traits,
      substance_items = substance_items,
      substance_min_valid = substance_min_valid,
      ses_equivalize = ses_equivalize,
      ses_split_within_year = ses_split_within_year,
      interval_filter = interval_filter,
      measurement = measurement
    )
  ))
}

#' run sensitivity analysis with different specifications
#' each specification changes exactly one thing against the main analysis:
#'   main            binary SES + count outcome (reference)
#'   three_ses       3-level SES grouping
#'   binary_sub      binary any-use outcome
#'   illicit_five    sedatives item excluded from the count
#'   lenient         interior-gap-tolerant weasel selection (larger sample)
#'   annual_interval respondents with uniform one-year wave spacing only
#'   ses_equivalized equivalized, within-year income median split
#'   parcels         multiple-indicator (parcel) measurement model
#' @param config configuration list
#' @param waves number of waves
#' @param specs which specifications to run (character subset of the above)
#' @param traits personality traits threaded into every specification
run_sensitivity_analysis <- function(config, waves = 5,
                                     specs = c(
                                       "main", "three_ses", "binary_sub",
                                       "illicit_five", "lenient",
                                       "annual_interval", "ses_equivalized",
                                       "parcels"
                                     ),
                                     traits = c("extr", "open", "cons")) {
  message("\n", strrep("=", 60))
  message("SENSITIVITY ANALYSIS")
  message(strrep("=", 60), "\n")

  # specification registry: label plus the arguments that deviate from main
  registry <- list(
    main = list(
      label = "Binary SES + Count sub (MAIN)",
      args = list(ses_method = "binary", substance_method = "sum")
    ),
    three_ses = list(
      label = "3-level SES + Count sub",
      args = list(ses_method = "composite", substance_method = "sum")
    ),
    binary_sub = list(
      label = "Binary SES + Binary sub",
      args = list(ses_method = "binary", substance_method = "any")
    ),
    illicit_five = list(
      label = "Sedatives excluded (five-item count)",
      args = list(
        ses_method = "binary", substance_method = "sum",
        substance_items = illicit_five_vars
      )
    ),
    lenient = list(
      label = "Lenient selection (interior gaps allowed)",
      args = list(
        ses_method = "binary", substance_method = "sum",
        weasel_scenario = "lenient_info_max"
      )
    ),
    annual_interval = list(
      label = "Uniform one-year intervals only",
      args = list(
        ses_method = "binary", substance_method = "sum",
        interval_filter = "annual"
      )
    ),
    ses_equivalized = list(
      label = "Equivalized within-year SES split",
      args = list(
        ses_method = "binary", substance_method = "sum",
        ses_equivalize = TRUE, ses_split_within_year = TRUE
      )
    ),
    parcels = list(
      label = "Parcel measurement model",
      args = list(
        ses_method = "binary", substance_method = "sum",
        measurement = "parcels"
      )
    )
  )

  specs <- intersect(specs, names(registry))
  results <- list()

  for (i in seq_along(specs)) {
    spec <- specs[i]
    message(
      "\n--- Analysis ", i, ": ", registry[[spec]]$label, " ---"
    )
    results[[spec]] <- tryCatch(
      do.call(
        run_riclpm_analysis,
        c(list(config = config, waves = waves, traits = traits),
          registry[[spec]]$args
        )
      ),
      error = function(e) {
        message("Analysis ", i, " (", spec, ") failed: ", e$message)
        NULL
      }
    )
  }

  message("\n", strrep("=", 60))
  message("SENSITIVITY ANALYSIS COMPLETE")
  message(strrep("=", 60))

  # compare results across specifications
  message("\nCross-specification comparison:")

  for (trait in traits) {
    message("\n", toupper(trait), ":")
    for (spec in specs) {
      res <- results[[spec]]
      if (!is.null(res) &&
        !is.null(res$hypothesis_results$by_trait[[trait]]$params)) {
        h1 <- res$hypothesis_results$by_trait[[trait]]$params$cl_summary %>%
          dplyr::filter(path_type == "H1_pers_to_sub")
        h2 <- res$hypothesis_results$by_trait[[trait]]$params$cl_summary %>%
          dplyr::filter(path_type == "H2_sub_to_pers")
        bounds <- res$hypothesis_results$by_trait[[trait]]$cl_bounds

        if (nrow(h1) > 0) {
          sig_h1 <- if (h1$pvalue_mean < 0.05) "*" else ""
          sig_h2 <- if (nrow(h2) > 0 && h2$pvalue_mean < 0.05) "*" else ""
          message("  ", registry[[spec]]$label, " (n=", nrow(res$data_wide), "):")
          b1 <- if (!is.null(bounds)) {
            bb <- bounds[bounds$path_type == "H1_pers_to_sub", ]
            if (nrow(bb) > 0) paste0(", |β| < ", formatC(bb$excl_bound, format = "f", digits = 3)) else ""
          } else {
            ""
          }
          message(
            "    H1 (pers->sub): β=", round(h1$est_std, 3),
            ", p=", round(h1$pvalue_mean, 3), sig_h1, b1
          )
          if (nrow(h2) > 0) {
            b2 <- if (!is.null(bounds)) {
              bb <- bounds[bounds$path_type == "H2_sub_to_pers", ]
              if (nrow(bb) > 0) paste0(", |β| < ", formatC(bb$excl_bound, format = "f", digits = 3)) else ""
            } else {
              ""
            }
            message(
              "    H2 (sub->pers): β=", round(h2$est_std, 3),
              ", p=", round(h2$pvalue_mean, 3), sig_h2, b2
            )
          }
        }
      } else {
        message("  ", registry[[spec]]$label, ": FAILED or no results")
      }
    }
  }

  invisible(results)
}
