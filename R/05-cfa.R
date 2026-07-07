# 05-cfa.R
# confirmatory factor analysis, measurement invariance, fit extraction.

# CFA AND MEASUREMENT INVARIANCE
# =============================================================================

#' build CFA model syntax
build_cfa_model <- function(factor_name, indicators) {
  paste0(factor_name, " =~ ", paste(indicators, collapse = " + "))
}

#' build bifactor CFA model (general + specific factors)
#' improves fit when items have cross-loadings
build_bifactor_model <- function(factor_name, indicators, n_specific = 2) {
  n_items <- length(indicators)
  items_per_specific <- ceiling(n_items / n_specific)

  # general factor loads on all items
  general <- paste0(
    "g_", factor_name, " =~ ",
    paste(indicators, collapse = " + ")
  )

  # specific factors load on subsets
  specific <- sapply(1:n_specific, function(i) {
    start_idx <- (i - 1) * items_per_specific + 1
    end_idx <- min(i * items_per_specific, n_items)
    spec_items <- indicators[start_idx:end_idx]
    paste0(
      "s", i, "_", factor_name, " =~ ",
      paste(spec_items, collapse = " + ")
    )
  })

  # orthogonality constraints
  ortho <- c(
    paste0("g_", factor_name, " ~~ 0*s1_", factor_name),
    paste0("g_", factor_name, " ~~ 0*s2_", factor_name),
    paste0("s1_", factor_name, " ~~ 0*s2_", factor_name)
  )

  paste(c(general, specific, ortho), collapse = "\n")
}

#' build parcel-based CFA model
#' reduces model complexity by averaging items into parcels
build_parcel_model <- function(factor_name, parcels) {
  paste0(factor_name, " =~ ", paste(parcels, collapse = " + "))
}

#' fit Big Five CFA with multiple model options
#' @param data personality data
#' @param dimension Big Five dimension
#' @param model_type "standard", "bifactor", or "parcel"
fit_b5_cfa <- function(data, dimension, model_type = "standard", ...) {
  items <- intersect(liss_cols("cp", b5_items[[dimension]]), names(data))

  if (model_type == "bifactor") {
    model <- build_bifactor_model(dimension, items)
  } else if (model_type == "parcel") {
    # create parcels first
    data <- create_b5_parcels(data, dimension, n_parcels = 3)
    parcel_names <- paste0(dimension, "_p", 1:3)
    model <- build_parcel_model(dimension, parcel_names)
  } else {
    model <- build_cfa_model(dimension, items)
  }

  # bifactor and parcel models are exploratory and known to be unstable for some
  # big five scales; suppress their non-convergence warnings and flag failure
  # explicitly rather than letting raw lavaan warnings reach the console
  fit <- tryCatch(
    suppressWarnings(lavaan::cfa(model, data = data, ...)),
    error = function(e) {
      warning("CFA failed for ", dimension, " (", model_type, "): ", e$message)
      NULL
    }
  )

  # fit measures are only meaningful for a converged solution; fitMeasures()
  # errors on a non-converged fit, so guard it
  converged <- !is.null(fit) && isTRUE(lavaan::lavInspect(fit, "converged"))
  if (!is.null(fit) && !converged) {
    message(
      "  CFA did not converge for ", dimension, " (", model_type,
      "); reporting NA fit measures"
    )
  }

  na_measures <- c(
    cfi = NA, tli = NA, rmsea = NA, srmr = NA,
    chisq = NA, df = NA, pvalue = NA
  )
  fit_measures <- if (converged) {
    tryCatch(
      lavaan::fitMeasures(fit, names(na_measures)),
      error = function(e) na_measures
    )
  } else {
    na_measures
  }

  list(
    dimension = dimension,
    model_type = model_type,
    model = model,
    fit = fit,
    fit_measures = fit_measures
  )
}

#' compare CFA model types for a dimension
compare_cfa_models <- function(data, dimension) {
  message("comparing CFA models for ", dimension, "...")

  results <- list(
    standard = fit_b5_cfa(data, dimension, "standard"),
    bifactor = fit_b5_cfa(data, dimension, "bifactor"),
    parcel = fit_b5_cfa(data, dimension, "parcel")
  )

  # compare fit
  comparison <- purrr::map_dfr(names(results), function(type) {
    r <- results[[type]]
    tibble::tibble(
      dimension = dimension,
      model_type = type,
      cfi = r$fit_measures["cfi"],
      tli = r$fit_measures["tli"],
      rmsea = r$fit_measures["rmsea"],
      srmr = r$fit_measures["srmr"]
    )
  })

  list(
    models = results,
    comparison = comparison
  )
}


# =============================================================================
# CFA FIT EXTRACTION HELPER
# =============================================================================

#' extract CFA fit indices from list of CFA results
#' @param cfa_list list with dimension, model, fit components
#' @param measures fit measures to extract
extract_cfa_fit <- function(cfa_list,
                            measures = c(
                              "cfi", "tli", "rmsea", "srmr",
                              "chisq", "df", "pvalue"
                            )) {
  purrr::map_dfr(cfa_list, function(x) {
    # handle new structure where fit_measures may be pre-computed
    if (!is.null(x$fit_measures) && !all(is.na(x$fit_measures))) {
      fm <- x$fit_measures
      return(tibble::tibble(
        dimension = x$dimension,
        model_type = x$model_type %||% "standard",
        cfi = fm["cfi"],
        tli = fm["tli"],
        rmsea = fm["rmsea"],
        srmr = fm["srmr"],
        chisq = fm["chisq"],
        df = fm["df"],
        pvalue = fm["pvalue"]
      ))
    }

    # fallback to extracting from fit object
    fit_vals <- tryCatch(
      {
        if (is.null(x$fit)) {
          stats::setNames(as.list(rep(NA_real_, length(measures))), measures)
        } else {
          fm <- lavaan::fitmeasures(x$fit, measures)
          as.list(fm)
        }
      },
      error = function(e) {
        stats::setNames(as.list(rep(NA_real_, length(measures))), measures)
      }
    )

    tibble::tibble(
      dimension = x$dimension,
      model_type = x$model_type %||% "standard",
      cfi = fit_vals$cfi,
      tli = fit_vals$tli,
      rmsea = fit_vals$rmsea,
      srmr = fit_vals$srmr,
      chisq = fit_vals$chisq,
      df = fit_vals$df,
      pvalue = fit_vals$pvalue
    )
  })
}

