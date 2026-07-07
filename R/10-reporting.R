# 10-reporting.R
# results formatting, interpretation, and publication-ready tables.

# =============================================================================
# RESULTS REPORTING
# =============================================================================

#' format RI-CLPM results as table
format_riclpm_results <- function(results) {
  purrr::map_dfr(names(results$by_trait), function(trait) {
    r <- results$by_trait[[trait]]
    if (is.null(r$params)) {
      return(NULL)
    }

    r$params$cl_summary %>%
      dplyr::mutate(
        trait = trait,
        hypothesis = dplyr::case_when(
          path_type == "H1_pers_to_sub" ~ paste0("H1", substr(trait, 1, 1)),
          path_type == "H2_sub_to_pers" ~ paste0("H2", substr(trait, 1, 1)),
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::select(hypothesis, trait, est_std, pvalue_mean, significant, n_paths)
  })
}

#' print hypothesis test summary
print_hypothesis_summary <- function(results) {
  cat("\n", strrep("=", 60), "\n")
  cat("HYPOTHESIS TEST SUMMARY\n")
  cat(strrep("=", 60), "\n\n")

  cat("PREDICTIONS (based on theoretical framework):\n")
  cat("-", strrep("-", 58), "\n")
  cat("H1: Personality → Substance Use (within-person)\n")
  cat("  H1a: Extraversion  → INCREASES substance use (+)\n")
  cat("  H1b: Openness      → INCREASES substance use (+)\n")
  cat("  H1c: Conscientiousness → DECREASES substance use (-)\n\n")

  cat("H2: Substance Use → Personality (within-person)\n")
  cat("  H2a-b: Substance use changes extraversion/openness (±)\n")
  cat("  H2c: Substance use DECREASES conscientiousness (-)\n\n")

  cat("H3: SES Moderation\n")
  cat("  Effects STRONGER for LOW SES individuals\n\n")

  if (!is.null(results$summary) && nrow(results$summary) > 0) {
    cat("OBSERVED RESULTS:\n")
    cat("-", strrep("-", 58), "\n")

    formatted <- format_riclpm_results(results)
    if (!is.null(formatted) && nrow(formatted) > 0) {
      for (i in seq_len(nrow(formatted))) {
        row <- formatted[i, ]
        sig_star <- if (row$pvalue_mean < 0.05) "*" else ""
        supported <- if (row$pvalue_mean < 0.05) "SIGNIFICANT" else "not supported"

        cat(sprintf(
          "  %s (%s): β = %6.3f, p = %.3f %s [%s]\n",
          row$hypothesis, row$trait,
          row$est_std, row$pvalue_mean, sig_star, supported
        ))
      }
    }

    # overall summary
    n_sig <- sum(formatted$pvalue_mean < 0.05, na.rm = TRUE)
    cat(
      "\n  Summary: ", n_sig, " of ", nrow(formatted),
      " cross-lagged paths significant (α = .05)\n"
    )
  }

  # H3 summary
  cat("\nH3 (SES Moderation):\n")
  cat("-", strrep("-", 58), "\n")
  for (trait in names(results$by_trait)) {
    r <- results$by_trait[[trait]]
    mg <- r$multigroup

    if (!is.null(mg)) {
      # prefer targeted cross-lagged test
      comp_cl <- mg$comparison_cl
      comp_all <- mg$comparison

      # use targeted test if available, otherwise overall
      if (!is.null(comp_cl) && nrow(comp_cl) >= 2) {
        p_val <- comp_cl$`Pr(>Chisq)`[2]
        test_type <- "cross-lagged"
      } else if (!is.null(comp_all) && nrow(comp_all) >= 2) {
        p_val <- comp_all$`Pr(>Chisq)`[2]
        test_type <- "all regressions"
      } else {
        next
      }

      sig <- if (!is.na(p_val) && p_val < 0.05) "SIGNIFICANT" else "ns"

      # check reliability
      reliable <- isTRUE(mg$results_reliable)
      reliability_note <- if (!reliable) " [⚠️ unreliable]" else ""

      cat(sprintf(
        "  %s: χ² diff p = %.4f [%s] (%s test)%s\n",
        toupper(trait), p_val, sig, test_type, reliability_note
      ))

      # report group sizes and imbalance
      if (!is.null(mg$group_ns)) {
        ns <- mg$group_ns
        cat(sprintf(
          "       Groups: %s\n",
          paste(names(ns), "=", ns, collapse = ", ")
        ))
      }
      if (!is.null(mg$imbalance_ratio) &&
        mg$imbalance_ratio > 5) {
        cat(sprintf(
          "       ⚠️  Imbalance ratio: %.1f:1 (interpret with caution)\n",
          mg$imbalance_ratio
        ))
      }

      # warn if overall significant but targeted not
      if (!is.null(comp_cl) && !is.null(comp_all) &&
        nrow(comp_cl) >= 2 && nrow(comp_all) >= 2) {
        p_overall <- comp_all$`Pr(>Chisq)`[2]
        p_target <- comp_cl$`Pr(>Chisq)`[2]
        if (!is.na(p_overall) && !is.na(p_target) &&
          p_overall < 0.05 && p_target >= 0.05) {
          cat("       Note: Overall test sig but cross-lagged test not\n")
          cat("             (SES moderates AR/variances, not cross-lagged)\n")
        }
      }
    } else if (!is.null(r$interaction) && !is.null(r$interaction$mod_params)) {
      # report interaction results if multi-group failed
      mp <- r$interaction$mod_params
      if (nrow(mp) > 0) {
        cat(sprintf("  %s (interaction model): ", toupper(trait)))
        for (j in seq_len(nrow(mp))) {
          sig <- if (mp$pvalue[j] < 0.05) "*" else ""
          cat(sprintf("%s β=%.3f%s ", mp$label[j], mp$std.all[j], sig))
        }
        cat("\n")
      }
    }
  }

  cat("\n", strrep("=", 60), "\n")
}


#' interpret RI-CLPM results, especially null findings
#' @param results hypothesis test results
#' @param cfa_fit CFA fit table
#' @param coverage coverage diagnostics
interpret_riclpm_results <- function(results, cfa_fit = NULL, coverage = NULL) {
  message("\n", strrep("=", 60))
  message("INTERPRETATION OF RESULTS")
  message(strrep("=", 60))

  # check for null main effects
  if (!is.null(results$summary)) {
    n_sig <- sum(results$summary$significant > 0, na.rm = TRUE)
    n_total <- nrow(results$summary)

    if (n_sig == 0) {
      message("\n*** NULL FINDINGS FOR CROSS-LAGGED EFFECTS (H1/H2) ***")
      message("\nAll cross-lagged paths are non-significant. This can indicate:")
      message("")
      message("1. TRAIT-LIKE STABILITY DOMINATES")
      message("   The Random Intercept captures stable between-person differences.")
      message("   High RI loadings (>0.80) suggest ~64%+ of variance is trait-like.")
      message("   After accounting for this stability, there's minimal within-person")
      message("   reciprocal causation between personality and substance use.")
      message("")
      message("2. THEORETICAL INTERPRETATION")
      message("   People who use more substances tend to have certain personality")
      message("   profiles (between-person association), BUT changes in personality")
      message("   do NOT predict changes in substance use within the same person")
      message("   (and vice versa). This is a valid and informative null finding.")
      message("")
      message("3. METHODOLOGICAL CONSIDERATIONS")
    }
  }

  # discuss CFA fit issues
  if (!is.null(cfa_fit)) {
    poor_cfa <- any(cfa_fit$cfi < 0.90, na.rm = TRUE)
    if (poor_cfa) {
      message("   - Poor CFA fit may attenuate true effects")
      message("     (measurement error -> biased estimates toward zero)")
    }
  }

  # discuss coverage issues
  if (!is.null(coverage) && coverage$min_coverage < 0.10) {
    message("   - Low pairwise coverage (<10%) may affect estimation")
    message("     Consider sensitivity analysis with fewer waves")
  }

  # H3 interpretation
  message("\n*** H3 (SES MODERATION) INTERPRETATION ***")

  # check multigroup reliability across traits
  has_unreliable <- FALSE
  for (trait in names(results$by_trait)) {
    r <- results$by_trait[[trait]]
    if (!is.null(r$multigroup) && !isTRUE(r$multigroup$results_reliable)) {
      has_unreliable <- TRUE
      break
    }
  }

  if (has_unreliable) {
    message("\n⚠️  RELIABILITY CONCERNS WITH SES MODERATION TESTS")
    message("   The multi-group χ² difference tests may be unreliable due to:")
    message("   - Severe group imbalance (one group dominates)")
    message("   - Small sample sizes in some SES groups")
    message("")
    message("   RECOMMENDATIONS:")
    message("   1. Use BINARY SES (median split) for more balanced groups")
    message("   2. Use interaction-based moderation instead of multi-group")
    message("   3. If reporting these results, clearly note the limitation")
    message("   4. Consider that 'significant' χ² may be driven by")
    message("      the large middle group, not true SES differences")
  } else {
    message("\n   Multi-group tests appear reliable (balanced groups)")
  }

  # discuss RI loadings if available
  message("\n*** KEY STATISTICS TO REPORT ***")
  message("   - Random Intercept loadings (proportion of between-person variance)")
  message("   - Cross-lagged standardized effects and CIs")
  message("   - Model fit indices (CFI, RMSEA, SRMR)")
  message("   - SES moderation test (chi-square difference or interaction term)")
  message("   - Group sizes and imbalance ratios for multi-group analyses")

  message("\n", strrep("=", 60))
}

#' create publication-ready summary table
#' @param results output from run_riclpm_analysis
#' @return data.frame suitable for publication
create_publication_table <- function(results) {
  hyp_results <- results$hypothesis_results

  # H1/H2 results
  h1h2_table <- purrr::map_dfr(names(hyp_results$by_trait), function(trait) {
    r <- hyp_results$by_trait[[trait]]
    if (is.null(r$params)) {
      return(NULL)
    }

    cl <- r$params$cl_summary
    ri <- r$ri_loadings

    # get RI loadings for interpretation
    ri_pers <- ri[grepl("RI_", ri$lhs) & !grepl("substance", ri$lhs), ]
    ri_sub <- ri[grepl("substance", ri$lhs), ]

    pers_between <- if (nrow(ri_pers) > 0) ri_pers$mean_loading[1]^2 else NA
    sub_between <- if (nrow(ri_sub) > 0) ri_sub$mean_loading[1]^2 else NA

    h1 <- cl[cl$path_type == "H1_pers_to_sub", ]
    h2 <- cl[cl$path_type == "H2_sub_to_pers", ]

    tibble::tibble(
      Trait = toupper(trait),
      `H1 (Pers→Sub) β` = if (nrow(h1) > 0) round(h1$est_std, 3) else NA,
      `H1 p` = if (nrow(h1) > 0) format.pval(h1$pvalue_mean, digits = 3) else NA,
      `H2 (Sub→Pers) β` = if (nrow(h2) > 0) round(h2$est_std, 3) else NA,
      `H2 p` = if (nrow(h2) > 0) format.pval(h2$pvalue_mean, digits = 3) else NA,
      `RI Pers (% between)` = round(pers_between * 100, 1),
      `RI Sub (% between)` = round(sub_between * 100, 1)
    )
  })

  # H3 results - use TARGETED cross-lagged test (comparison_cl), not overall test
  h3_table <- purrr::map_dfr(names(hyp_results$by_trait), function(trait) {
    r <- hyp_results$by_trait[[trait]]
    mg <- r$multigroup

    if (is.null(mg)) {
      return(tibble::tibble(
        Trait = toupper(trait),
        `H3 χ² diff (cross-lagged)` = NA,
        `H3 p` = NA,
        `H3 Overall χ²` = NA,
        `Overall p` = NA,
        `Group Ns` = NA,
        `Imbalance` = NA,
        Reliable = NA
      ))
    }

    # prefer targeted cross-lagged test (comparison_cl) over overall test
    comp_cl <- mg$comparison_cl
    comp_all <- mg$comparison

    # targeted test values
    chi_diff_cl <- if (!is.null(comp_cl) && nrow(comp_cl) >= 2) comp_cl$`Chisq diff`[2] else NA
    p_val_cl <- if (!is.null(comp_cl) && nrow(comp_cl) >= 2) comp_cl$`Pr(>Chisq)`[2] else NA

    # overall test values (for comparison)
    chi_diff_all <- if (!is.null(comp_all) && nrow(comp_all) >= 2) comp_all$`Chisq diff`[2] else NA
    p_val_all <- if (!is.null(comp_all) && nrow(comp_all) >= 2) comp_all$`Pr(>Chisq)`[2] else NA

    group_str <- paste(names(mg$group_ns), "=", mg$group_ns, collapse = "; ")

    tibble::tibble(
      Trait = toupper(trait),
      `H3 χ² diff (cross-lagged)` = round(chi_diff_cl, 2),
      `H3 p` = format.pval(p_val_cl, digits = 3),
      `H3 Overall χ²` = round(chi_diff_all, 2),
      `Overall p` = format.pval(p_val_all, digits = 3),
      `Group Ns` = group_str,
      `Imbalance` = paste0(round(mg$imbalance_ratio, 1), ":1"),
      Reliable = ifelse(isTRUE(mg$results_reliable), "Yes", "No*")
    )
  })

  list(
    h1_h2 = h1h2_table,
    h3 = h3_table,
    cfa_fit = results$cfa_fit,
    settings = results$settings
  )
}

#' print publication table
#' @param pub_table output from create_publication_table
print_publication_table <- function(pub_table) {
  cat("\n", strrep("=", 80), "\n")
  cat("PUBLICATION-READY RESULTS SUMMARY\n")
  cat(strrep("=", 80), "\n\n")

  cat(
    "Settings: Waves=", pub_table$settings$waves,
    ", SES=", pub_table$settings$ses_method,
    ", Substance=", pub_table$settings$substance_method,
    ", Estimator=", pub_table$settings$estimator, "\n\n"
  )

  cat("H1/H2: Cross-Lagged Effects (Within-Person)\n")
  cat(strrep("-", 80), "\n")
  print(pub_table$h1_h2, n = Inf)
  cat("\n")

  cat("H3: SES Moderation (Targeted Cross-Lagged Test)\n")
  cat(strrep("-", 80), "\n")
  cat("NOTE: 'H3 χ² diff (cross-lagged)' tests ONLY cross-lagged paths (the hypothesis).\n")
  cat("      'H3 Overall χ²' tests ALL regressions including autoregressive paths.\n")
  cat("      Use the targeted test (cross-lagged) for H3 conclusions.\n\n")
  print(pub_table$h3, n = Inf)

  # check for unreliable results
  if (any(pub_table$h3$Reliable == "No*", na.rm = TRUE)) {
    cat("\n* Results marked as unreliable due to severe group imbalance.\n")
    cat("  Interpret with caution or use interaction-based moderation.\n")
  }

  # check for moderation paradox (overall sig but targeted not)
  h3_data <- pub_table$h3
  paradox_traits <- c()
  for (i in seq_len(nrow(h3_data))) {
    overall_p <- suppressWarnings(as.numeric(h3_data$`Overall p`[i]))
    target_p <- suppressWarnings(as.numeric(h3_data$`H3 p`[i]))
    if (!is.na(overall_p) && !is.na(target_p) &&
      overall_p < 0.05 && target_p >= 0.05) {
      paradox_traits <- c(paradox_traits, h3_data$Trait[i])
    }
  }
  if (length(paradox_traits) > 0) {
    cat(
      "\n⚠️  MODERATION INTERPRETATION NOTE for ",
      paste(paradox_traits, collapse = ", "), ":\n"
    )
    cat("   Overall test significant but targeted test not.\n")
    cat("   SES moderates autoregressive stability or variances,\n")
    cat("   NOT the hypothesized cross-lagged effects.\n")
  }

  cat("\nCFA Fit:\n")
  cat(strrep("-", 80), "\n")
  print(pub_table$cfa_fit)

  # interpretation notes
  cat("\n", strrep("=", 80), "\n")
  cat("INTERPRETATION NOTES:\n")
  cat("- RI loadings show % variance due to stable between-person differences\n")
  cat("- High RI (>80%) suggests trait-like stability dominates\n")
  cat("- Null cross-lagged effects with high RI = valid null finding\n")
  cat("- CFA CFI < 0.90 may attenuate true effects\n")
  cat(strrep("=", 80), "\n")
}

