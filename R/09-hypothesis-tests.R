# 09-hypothesis-tests.R
# parameter extraction and H1/H2/H3 testing across traits.

# =============================================================================
# HYPOTHESIS TESTING
# =============================================================================

#' extract RI-CLPM parameters for hypothesis testing
#' @param fit lavaan fit object
#' @param personality personality trait name
extract_riclpm_params <- function(fit, personality) {
  if (is.null(fit)) {
    return(NULL)
  }

  params <- lavaan::parameterEstimates(fit, standardized = TRUE)

  # identify key paths
  # c = personality -> outcome (H1)
  # d = outcome -> personality (H2)
  key_paths <- params %>%
    dplyr::filter(
      grepl("^wp_", lhs) & op == "~" & grepl("^wp_", rhs)
    ) %>%
    dplyr::mutate(
      path_type = dplyr::case_when(
        grepl("substance_use", lhs) & grepl(personality, rhs) ~ "H1_pers_to_sub",
        grepl(personality, lhs) & grepl("substance_use", rhs) ~ "H2_sub_to_pers",
        grepl(personality, lhs) & grepl(personality, rhs) ~ "AR_personality",
        grepl("substance_use", lhs) & grepl("substance_use", rhs) ~ "AR_substance",
        TRUE ~ "other"
      )
    )

  # summarize cross-lagged effects
  # NOTE: with constrain_ar/constrain_cl = TRUE each path_type collapses one
  # labeled parameter (identical est and pvalue across rows; std.all varies
  # slightly per wave through the wave-specific latent variances); with
  # unconstrained paths these means average genuinely distinct parameters
  cl_effects <- key_paths %>%
    dplyr::filter(path_type %in% c("H1_pers_to_sub", "H2_sub_to_pers")) %>%
    dplyr::group_by(path_type) %>%
    dplyr::summarize(
      est_mean = mean(est),
      est_std = mean(std.all, na.rm = TRUE),
      pvalue_mean = mean(pvalue),
      significant = sum(pvalue < 0.05),
      n_paths = dplyr::n(),
      .groups = "drop"
    )

  list(
    all_params = params,
    key_paths = key_paths,
    cl_summary = cl_effects
  )
}

#' test all hypotheses for a given personality trait
#' @param data wide format data
#' @param personality personality trait
#' @param waves number of waves
#' @param estimator "ML" or "MLR" (robust)
test_hypotheses <- function(data, personality, waves = 5, estimator = "MLR",
                            measurement = c("composite", "parcels")) {
  measurement <- match.arg(measurement)
  message(
    "\n=== Testing hypotheses for ", toupper(personality),
    if (measurement == "parcels") " (parcel measurement model)" else "",
    " ===\n"
  )

  # check variance first
  var_check <- check_variance_by_group(data, "ses", "substance_use")
  if (var_check$has_zero_variance) {
    message("WARNING: Zero variance detected in substance_use for some SES groups")
    message("  Multi-group model may fail. Will attempt interaction model as fallback.")
  }

  # fit basic RI-CLPM (single-indicator composite, or the multiple-indicator
  # parcel specification correcting within-person estimates for unreliability)
  riclpm <- if (measurement == "parcels") {
    fit_riclpm_mi(data, personality, "substance_use", waves,
      estimator = estimator
    )
  } else {
    fit_riclpm(data, personality, "substance_use", waves,
      estimator = estimator
    )
  }

  # extract parameters
  params <- NULL
  ri_loadings <- NULL
  cl_bounds <- NULL
  if (!is.null(riclpm$fit)) {
    params <- extract_riclpm_params(riclpm$fit, personality)

    # extract random intercept loadings (key for interpretation)
    all_params <- lavaan::parameterEstimates(riclpm$fit, standardized = TRUE)
    ri_loadings <- all_params %>%
      dplyr::filter(op == "=~", grepl("^RI_", lhs)) %>%
      dplyr::group_by(lhs) %>%
      dplyr::summarize(
        mean_loading = mean(std.all, na.rm = TRUE),
        .groups = "drop"
      )

    # report RI loadings (important for understanding null findings)
    message("Random Intercept (trait-like stability):")
    for (i in seq_len(nrow(ri_loadings))) {
      pct_between <- round(ri_loadings$mean_loading[i]^2 * 100, 1)
      message(
        "  ", ri_loadings$lhs[i], ": λ = ",
        round(ri_loadings$mean_loading[i], 3),
        " (", pct_between, "% between-person variance)"
      )
    }

    message("\nH1 (", personality, " -> substance use):")
    h1 <- params$cl_summary %>%
      dplyr::filter(path_type == "H1_pers_to_sub")
    if (nrow(h1) > 0) {
      message(
        "  β = ", round(h1$est_std, 3),
        ", p = ", round(h1$pvalue_mean, 4),
        " (", h1$significant, "/", h1$n_paths, " significant)"
      )
    }

    message("H2 (substance use -> ", personality, "):")
    h2 <- params$cl_summary %>%
      dplyr::filter(path_type == "H2_sub_to_pers")
    if (nrow(h2) > 0) {
      message(
        "  β = ", round(h2$est_std, 3),
        ", p = ", round(h2$pvalue_mean, 4),
        " (", h2$significant, "/", h2$n_paths, " significant)"
      )
    }

    # model fit
    fit_idx <- lavaan::fitmeasures(
      riclpm$fit,
      c("cfi", "tli", "rmsea", "srmr")
    )
    message(
      "\nModel fit: CFI=", round(fit_idx["cfi"], 3),
      ", TLI=", round(fit_idx["tli"], 3),
      ", RMSEA=", round(fit_idx["rmsea"], 3),
      ", SRMR=", round(fit_idx["srmr"], 3)
    )

    # 95% CI exclusion bounds for the constrained cross-lagged paths: reframes
    # the null as "effects beyond |bound| are inconsistent with the data"
    cl_bounds <- extract_cl_bounds(riclpm$fit, personality)
    print_cl_bounds(cl_bounds)
  }

  # fit multi-group model for H3 (composite specification only; the parcel
  # model tests measurement robustness of H1/H2, moderation stays on the main
  # composite models)
  mg <- NULL
  mg_success <- FALSE
  interaction_result <- NULL

  if (measurement == "parcels") {
    message("\nH3 (SES moderation): skipped under the parcel measurement model")
  } else {
  message("\nH3 (SES moderation):")
  mg <- fit_riclpm_multigroup(data, personality, "substance_use", waves,
    estimator = estimator
  )

  # check if multi-group succeeded
  mg_success <- !is.null(mg$comparison) && !is.null(mg$fit_configural)

  if (mg_success) {
    # test 1: all regressions (original, includes autoregressive)
    chi_diff <- mg$comparison$`Chisq diff`[2]
    df_diff <- mg$comparison$`Df diff`[2]
    p_mod <- mg$comparison$`Pr(>Chisq)`[2]

    message("  Multi-group test (all regressions):")
    message(
      "    χ² diff = ", round(chi_diff, 2),
      ", df = ", df_diff,
      ", p = ", format.pval(p_mod, digits = 3)
    )

    # test 2: cross-lagged paths only (targeted H3 test)
    if (!is.null(mg$comparison_cl)) {
      chi_diff_cl <- mg$comparison_cl$`Chisq diff`[2]
      df_diff_cl <- mg$comparison_cl$`Df diff`[2]
      p_mod_cl <- mg$comparison_cl$`Pr(>Chisq)`[2]

      message("  Targeted test (cross-lagged paths only):")
      message(
        "    χ² diff = ", round(chi_diff_cl, 2),
        ", df = ", df_diff_cl,
        ", p = ", format.pval(p_mod_cl, digits = 3)
      )

      # interpret difference between tests
      if (!is.na(p_mod) && p_mod < 0.05 &&
        !is.na(p_mod_cl) && p_mod_cl >= 0.05) {
        message("\n    ⚠️  IMPORTANT: Overall test significant but cross-lagged test not.")
        message("       SES may moderate autoregressive stability rather than")
        message("       the hypothesized cross-lagged effects (variance and")
        message("       covariance components are held equal across groups here).")
      }
    }

    # report reliability status
    if (!isTRUE(mg$results_reliable)) {
      message("\n    ⚠️  RELIABILITY WARNING:")
      if (!is.null(mg$reliability_issues$severe_imbalance)) {
        message(
          "       Severe group imbalance (",
          round(mg$imbalance_ratio, 1), ":1 ratio)"
        )
      }
      if (!is.null(mg$reliability_issues$moderate_imbalance)) {
        message(
          "       Moderate group imbalance (",
          round(mg$imbalance_ratio, 1), ":1 ratio)"
        )
      }
      if (!is.null(mg$reliability_issues$small_groups)) {
        message("       Some groups have <100 observations")
      }
      if (!is.null(mg$reliability_issues$majority_dominates)) {
        message("       Majority group dominates (>90%)")
      }
      message("       Interpret χ² test with caution")
    }

    # primary conclusion based on targeted test if available
    p_primary <- if (!is.null(mg$comparison_cl)) {
      mg$comparison_cl$`Pr(>Chisq)`[2]
    } else {
      p_mod
    }

    if (!is.na(p_primary) && p_primary < 0.05) {
      if (isTRUE(mg$results_reliable)) {
        message("    -> SES MODERATES the cross-lagged relationships")
      } else {
        message("    -> SES appears to moderate, but result may be UNRELIABLE")
      }
    } else if (!is.na(p_mod) && p_mod < 0.05) {
      message("    -> SES moderates overall model, but NOT specifically cross-lagged paths")
    } else {
      message("    -> No significant SES moderation")
    }

    # report group-specific effects
    if (!is.null(mg$group_params)) {
      report_group_effects(mg$group_params, personality, data)
    }
  } else {
    message("  Multi-group model failed (likely due to zero variance in subgroups)")
    message("  Attempting interaction-based moderation test...")

    # fallback to interaction model
    interaction_result <- fit_riclpm_interaction(
      data, personality,
      "substance_use", waves,
      estimator
    )

    if (!is.null(interaction_result$fit)) {
      message("\n  Interaction model results:")
      if (!is.null(interaction_result$mod_params) &&
        nrow(interaction_result$mod_params) > 0) {
        mod_po <- interaction_result$mod_params[
          interaction_result$mod_params$label == "mod_po",
        ]
        mod_op <- interaction_result$mod_params[
          interaction_result$mod_params$label == "mod_op",
        ]

        if (nrow(mod_po) > 0) {
          sig_po <- if (mod_po$pvalue[1] < 0.05) "*" else ""
          message(
            "    SES moderates ", personality, "->substance: β=",
            round(mod_po$std.all[1], 3), ", p=",
            round(mod_po$pvalue[1], 3), sig_po
          )
        }
        if (nrow(mod_op) > 0) {
          sig_op <- if (mod_op$pvalue[1] < 0.05) "*" else ""
          message(
            "    SES moderates substance->", personality, ": β=",
            round(mod_op$std.all[1], 3), ", p=",
            round(mod_op$pvalue[1], 3), sig_op
          )
        }
      }
    } else {
      message("  Interaction model also failed")
    }
  }
  }

  list(
    personality = personality,
    measurement = measurement,
    riclpm = riclpm,
    multigroup = mg,
    interaction = interaction_result,
    params = params,
    ri_loadings = ri_loadings,
    cl_bounds = cl_bounds,
    mg_success = mg_success
  )
}

#' helper to report group-specific effects
report_group_effects <- function(group_params, personality, data) {
  message("\nGroup-specific cross-lagged effects:")
  cl_by_group <- group_params %>%
    dplyr::filter(grepl("substance_use", lhs) & grepl(personality, rhs) |
      grepl(personality, lhs) & grepl("substance_use", rhs)) %>%
    dplyr::select(lhs, op, rhs, est, std.all, pvalue, group)

  if (nrow(cl_by_group) > 0) {
    # H1 by group
    h1_grp <- cl_by_group %>%
      dplyr::filter(grepl("substance_use", lhs)) %>%
      dplyr::group_by(group) %>%
      dplyr::summarize(
        path = "H1 (pers->sub)",
        beta = mean(std.all, na.rm = TRUE),
        p = mean(pvalue, na.rm = TRUE),
        .groups = "drop"
      )

    # H2 by group
    h2_grp <- cl_by_group %>%
      dplyr::filter(grepl(personality, lhs)) %>%
      dplyr::group_by(group) %>%
      dplyr::summarize(
        path = "H2 (sub->pers)",
        beta = mean(std.all, na.rm = TRUE),
        p = mean(pvalue, na.rm = TRUE),
        .groups = "drop"
      )

    grp_summary <- dplyr::bind_rows(h1_grp, h2_grp)
    if (nrow(grp_summary) > 0) {
      ses_levels <- levels(factor(data$ses[!is.na(data$ses)]))
      for (g in unique(grp_summary$group)) {
        grp_name <- if (g <= length(ses_levels)) ses_levels[g] else paste0("Group ", g)
        message("  ", grp_name, " SES:")
        grp_data <- grp_summary[grp_summary$group == g, ]
        for (i in seq_len(nrow(grp_data))) {
          sig_star <- if (grp_data$p[i] < 0.05) "*" else ""
          message(
            "    ", grp_data$path[i], ": β=",
            round(grp_data$beta[i], 3),
            ", p=", round(grp_data$p[i], 3), sig_star
          )
        }
      }

      # check for moderation paradox
      any_sig <- any(grp_summary$p < 0.05)
      if (!any_sig) {
        message("\n  ⚠️  MODERATION INTERPRETATION NOTE:")
        message("     The individual group cross-lagged paths are all")
        message("     non-significant. If the overall χ² test is significant,")
        message("     then, because the variance and covariance components are")
        message("     held equal across groups here, this points to SES")
        message("     moderating the autoregressive stability paths rather")
        message("     than the cross-lagged effects specifically.")
        message("     Report this nuance when interpreting H3.")
      }
    }
  }
}

#' run complete hypothesis testing for all relevant traits
#' H1a/H2a: extraversion, H1b/H2b: openness, H1c/H2c: conscientiousness
#' @param data wide format data
#' @param waves number of waves
#' @param estimator "ML" or "MLR" (robust)
#' @param traits personality traits to test (extend with "neur" and "agre" for
#'   the full-domain sensitivity)
#' @param measurement "composite" (single-indicator, default) or "parcels"
#'   (multiple-indicator RI-CLPM; requires parcel columns in data)
run_hypothesis_tests <- function(data, waves = 5, estimator = "MLR",
                                 traits = c("extr", "open", "cons"),
                                 measurement = "composite") {
  results <- purrr::map(
    traits,
    ~ test_hypotheses(data, .x, waves, estimator, measurement = measurement)
  )
  names(results) <- traits

  # summary table
  summary_tbl <- purrr::map_dfr(results, function(r) {
    if (is.null(r$params)) {
      return(NULL)
    }

    r$params$cl_summary %>%
      dplyr::mutate(personality = r$personality)
  })

  list(
    by_trait = results,
    summary = summary_tbl
  )
}



# =============================================================================
# CROSS-LAGGED CONFIDENCE BOUNDS
# =============================================================================

#' extract standardized cross-lagged estimates with confidence-interval bounds
#' with the equality-constrained specification the H1 and H2 paths are single
#' parameters whose standardized value varies slightly by wave; this collapses
#' each path to its mean standardized estimate, the union of the wave-wise CIs,
#' and the exclusion bound (the largest standardized effect the data leave
#' compatible in either direction). delta-method CIs via
#' lavaan::standardizedSolution(), so they inherit the robust (MLR) covariance.
#' @param fit lavaan fit object
#' @param personality personality trait name
#' @param outcome outcome variable name
#' @param level confidence level (default 0.95)
extract_cl_bounds <- function(fit, personality, outcome = "substance_use",
                              level = 0.95) {
  if (is.null(fit)) {
    return(NULL)
  }

  std <- tryCatch(
    lavaan::standardizedSolution(fit, type = "std.all", level = level),
    error = function(e) NULL
  )
  if (is.null(std)) {
    return(NULL)
  }

  reg <- std[std$op == "~" &
    grepl("^wp_", std$lhs) & grepl("^wp_", std$rhs), , drop = FALSE]
  if (nrow(reg) == 0) {
    return(NULL)
  }

  path_type <- ifelse(
    grepl(outcome, reg$lhs) & grepl(personality, reg$rhs), "H1_pers_to_sub",
    ifelse(
      grepl(personality, reg$lhs) & grepl(outcome, reg$rhs), "H2_sub_to_pers",
      "other"
    )
  )
  keep <- path_type != "other"
  reg <- reg[keep, , drop = FALSE]
  path_type <- path_type[keep]
  if (nrow(reg) == 0) {
    return(NULL)
  }

  out <- do.call(rbind, lapply(split(seq_len(nrow(reg)), path_type), function(ix) {
    r <- reg[ix, , drop = FALSE]
    data.frame(
      path_type = path_type[ix][1],
      est_std = round(mean(r$est.std), 4),
      ci_lower = round(min(r$ci.lower), 4),
      ci_upper = round(max(r$ci.upper), 4),
      excl_bound = round(max(abs(c(r$ci.lower, r$ci.upper))), 4),
      n_rows = nrow(r),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

#' print cross-lagged confidence bounds
#' @param bounds return value of extract_cl_bounds()
#' @param level confidence level used (display only)
print_cl_bounds <- function(bounds, level = 0.95) {
  if (is.null(bounds) || nrow(bounds) == 0) {
    return(invisible(NULL))
  }
  message(
    "\n", round(level * 100),
    "% CIs on standardized cross-lagged paths (exclusion bounds):"
  )
  for (i in seq_len(nrow(bounds))) {
    message(
      "  ", bounds$path_type[i],
      ": β = ", formatC(bounds$est_std[i], format = "f", digits = 3),
      ", CI [", formatC(bounds$ci_lower[i], format = "f", digits = 3),
      ", ", formatC(bounds$ci_upper[i], format = "f", digits = 3),
      "]; effects beyond |",
      formatC(bounds$excl_bound[i], format = "f", digits = 3), "| excluded"
    )
  }
  invisible(bounds)
}
