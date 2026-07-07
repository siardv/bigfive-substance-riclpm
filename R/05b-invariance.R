# 05b-invariance.R
# item-level longitudinal measurement invariance for the Big Five scales,
# run on the weasel-selected occasions. moves the invariance evidence into the
# reproducible pipeline so the manuscript claim traces to the rendered
# transcript rather than an interactive session.

# =============================================================================
# DATA PREPARATION
# =============================================================================

#' reshape reverse-scored items of one dimension to person-by-wave-index wide
#' rows are the weasel-selected person-waves; columns are cp items suffixed
#' with the dense wave index (_w1.._wL), one row per person.
#' @param p personality module data (long, item columns present)
#' @param wave_map weasel wave_map (nomem_encr, wave_idx, wavenr_orig)
#' @param dimension Big Five dimension name
prepare_invariance_items <- function(p, wave_map, dimension) {
  if (is.null(wave_map)) {
    warning("prepare_invariance_items: wave_map is NULL")
    return(NULL)
  }

  items <- liss_cols("cp", b5_items[[dimension]])
  items_present <- intersect(items, names(p))
  if (length(items_present) < 3) {
    warning("prepare_invariance_items: too few items found for ", dimension)
    return(NULL)
  }

  d <- p
  if (!"wavenr" %in% names(d) && "wave" %in% names(d)) {
    d$wavenr <- wave_to_wavenr(d$wave)
  }

  # reverse-score in place so the invariance models see the analysis metric
  orig_positions <- match(items_present, items)
  reverse_here <- which(orig_positions %in% b5_reverse[[dimension]])
  for (i in reverse_here) {
    d[[items_present[i]]] <- 5 + 1 - d[[items_present[i]]]
  }

  d <- d %>%
    dplyr::select(dplyr::all_of(c("nomem_encr", "wavenr", items_present))) %>%
    dplyr::distinct(nomem_encr, wavenr, .keep_all = TRUE)

  wm <- wave_map
  names(wm)[names(wm) == "wavenr_orig"] <- "wavenr"

  long_sel <- dplyr::inner_join(wm, d, by = c("nomem_encr", "wavenr"))

  # manual wide pivot: one merge per dense wave index
  waves <- sort(unique(long_sel$wave_idx))
  wide <- NULL
  for (t in waves) {
    dt <- long_sel[long_sel$wave_idx == t, c("nomem_encr", items_present), drop = FALSE]
    names(dt)[-1] <- paste0(items_present, "_w", t)
    wide <- if (is.null(wide)) dt else merge(wide, dt, by = "nomem_encr", all = TRUE)
  }

  attr(wide, "items") <- items_present
  attr(wide, "waves") <- waves
  wide
}

# =============================================================================
# MODEL SYNTAX
# =============================================================================

#' lavaan syntax for one level of the longitudinal invariance hierarchy
#' fixed-factor identification: the wave-1 factor variance is 1 and its mean 0
#' at every level; equality constraints are added stepwise. same-item residual
#' covariances across all wave pairs are included at every level, the standard
#' specification for repeated indicators.
#' @param items item column stems (without the _wT suffix)
#' @param waves dense wave indices
#' @param level "configural", "metric", "scalar", or "strict"
build_longitudinal_invariance_model <- function(items, waves,
                                                level = c(
                                                  "configural", "metric",
                                                  "scalar", "strict"
                                                )) {
  level <- match.arg(level)
  n_items <- length(items)
  ind <- function(k, t) paste0(items[k], "_w", t)
  fac <- function(t) paste0("f", t)

  metric_up <- level %in% c("metric", "scalar", "strict")
  scalar_up <- level %in% c("scalar", "strict")
  strict_up <- level == "strict"

  # loadings: all freed from the marker convention (NA* on the first), scale
  # set through the factor variance instead; labels shared across waves from
  # the metric level upward
  meas <- sapply(waves, function(t) {
    rhs <- sapply(seq_len(n_items), function(k) {
      lab <- if (metric_up) paste0("l", k, "*") else ""
      pre <- if (k == 1) paste0("NA*", lab) else lab
      paste0(pre, ind(k, t))
    })
    paste0(fac(t), " =~ ", paste(rhs, collapse = " + "))
  })

  # factor variances: wave 1 fixed to 1 always; later waves fixed to 1 in the
  # configural model (per-wave standardized metric) and freed once loadings
  # carry the metric
  fvar <- sapply(waves, function(t) {
    if (t == waves[1] || !metric_up) {
      paste0(fac(t), " ~~ 1*", fac(t))
    } else {
      paste0(fac(t), " ~~ ", fac(t))
    }
  })

  # factor means: wave 1 fixed to 0 always; later waves fixed to 0 until
  # intercept invariance identifies them
  fmean <- sapply(waves, function(t) {
    if (t == waves[1] || !scalar_up) {
      paste0(fac(t), " ~ 0*1")
    } else {
      paste0(fac(t), " ~ 1")
    }
  })

  # intercepts: free until the scalar level, then shared labels per item
  ints <- unlist(lapply(waves, function(t) {
    sapply(seq_len(n_items), function(k) {
      if (scalar_up) {
        paste0(ind(k, t), " ~ i", k, "*1")
      } else {
        paste0(ind(k, t), " ~ 1")
      }
    })
  }))

  # residual variances: free until the strict level, then shared per item
  resv <- unlist(lapply(waves, function(t) {
    sapply(seq_len(n_items), function(k) {
      if (strict_up) {
        paste0(ind(k, t), " ~~ r", k, "*", ind(k, t))
      } else {
        paste0(ind(k, t), " ~~ ", ind(k, t))
      }
    })
  }))

  # correlated uniqueness: same item across every wave pair
  pairs <- utils::combn(waves, 2)
  cu <- unlist(lapply(seq_len(n_items), function(k) {
    apply(pairs, 2, function(st) paste0(ind(k, st[1]), " ~~ ", ind(k, st[2])))
  }))

  # factor covariances across waves are freely estimated (lavaan default for
  # exogenous latents), so no explicit lines are needed
  paste(c(
    paste0("# longitudinal invariance: ", level),
    meas, fvar, fmean, ints, resv,
    "# correlated uniqueness (same item across waves)",
    cu
  ), collapse = "\n")
}

# =============================================================================
# FITTING AND COMPARISON
# =============================================================================

#' fit the longitudinal invariance hierarchy for one dimension
#' @param items_wide return value of prepare_invariance_items()
#' @param estimator "ML" or "MLR" (robust)
#' @param levels invariance levels to fit, in order
fit_longitudinal_invariance <- function(items_wide,
                                        estimator = "MLR",
                                        levels = c(
                                          "configural", "metric",
                                          "scalar", "strict"
                                        )) {
  if (is.null(items_wide)) {
    return(NULL)
  }
  items <- attr(items_wide, "items")
  waves <- attr(items_wide, "waves")

  fits <- list()
  for (lv in levels) {
    model <- build_longitudinal_invariance_model(items, waves, level = lv)
    fits[[lv]] <- tryCatch(
      lavaan::cfa(model,
        data = items_wide,
        missing = "fiml",
        estimator = estimator
      ),
      error = function(e) {
        warning("invariance fit failed at ", lv, ": ", e$message)
        NULL
      }
    )
  }

  ok <- !vapply(fits, is.null, logical(1))
  fit_tab <- do.call(rbind, lapply(names(fits)[ok], function(lv) {
    fm <- lavaan::fitmeasures(
      fits[[lv]],
      c("chisq.scaled", "df", "cfi.robust", "rmsea.robust", "srmr")
    )
    data.frame(
      level = lv,
      chisq = round(unname(fm["chisq.scaled"]), 1),
      df = unname(fm["df"]),
      cfi = round(unname(fm["cfi.robust"]), 4),
      rmsea = round(unname(fm["rmsea.robust"]), 4),
      srmr = round(unname(fm["srmr"]), 4),
      stringsAsFactors = FALSE
    )
  }))
  if (!is.null(fit_tab)) {
    fit_tab$delta_cfi <- c(NA, round(diff(fit_tab$cfi), 4))
  }

  lrt <- if (sum(ok) >= 2) {
    tryCatch(
      do.call(lavaan::lavTestLRT, unname(fits[ok])),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  # label rows by invariance level; the default deparses the model objects
  # into the row names and inflates the printed table to megabytes per line
  if (!is.null(lrt)) rownames(lrt) <- names(fits)[ok]

  list(fits = fits, fit_table = fit_tab, lrt = lrt)
}

#' run the invariance hierarchy for several dimensions on the selected sample
#' pipeline-facing entry point; prints fit tables and scaled LRTs per dimension
#' @param p personality module data
#' @param weasel_selection return value of select_weasel_subset()
#' @param dims dimensions to test
#' @param estimator "ML" or "MLR"
run_longitudinal_invariance <- function(p, weasel_selection,
                                        dims = c("extr", "open", "cons"),
                                        estimator = "MLR") {
  results <- list()
  for (dim in dims) {
    message("\nLongitudinal invariance: ", toupper(dim))
    items_wide <- prepare_invariance_items(p, weasel_selection$wave_map, dim)
    res <- fit_longitudinal_invariance(items_wide, estimator = estimator)
    if (!is.null(res$fit_table)) {
      print(res$fit_table)
    }
    if (!is.null(res$lrt)) {
      message("  scaled LRTs (each level against the previous):")
      print(res$lrt)
    }
    results[[dim]] <- res
  }
  invisible(results)
}
