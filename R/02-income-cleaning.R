# 02-income-cleaning.R
# household-income outlier detection and constrained imputation.

# =============================================================================
# HELPER FUNCTIONS: INCOME CLEANING
# =============================================================================

#' round to nearest power of 10
#' @param x numeric vector
#' @param type "nearest", "trunc", or "ceiling"
#' @param cutoff threshold for rounding up (default 0.75)
get_power10 <- function(x, type = "nearest", cutoff = 0.75) {
  x <- rna(x)
  if (length(x) == 0) {
    return(numeric(0))
  }

  # handle zeros and negatives
  x[x <= 0] <- NA
  x <- rna(x)
  if (length(x) == 0) {
    return(numeric(0))
  }

  trunc_x <- 10^trunc(log10(x))
  ceiling_x <- 10^ceiling(log10(x))

  if (type == "trunc") {
    trunc_x
  } else if (type == "ceiling") {
    ceiling_x
  } else {
    # nearest
    sapply(seq_along(x), function(i) {
      if (x[i] / trunc_x[i] > (cutoff * 10)) {
        ceiling_x[i]
      } else {
        trunc_x[i]
      }
    })
  }
}

#' calculate lag difference (log change)
#' @param df data frame with nethh column
#' @param power10 also compute power of 10
lag_diff <- function(df, power10 = TRUE) {
  x <- rna(df$nethh)
  if (length(x) < 2) {
    df$diff[!is.na(df$nethh)] <- 0
    if (power10) df$power10[!is.na(df$nethh)] <- get_power10(x)
    return(df)
  }

  # log difference relative to lagged value
  change <- function(k) {
    k <- rna(k)
    if (length(k) < 2) {
      return(0)
    }
    delt <- log(k / dplyr::lag(k, 1))
    signif(mean(abs(rna(delt))), 2)
  }

  d <- sapply(seq_along(x), function(i) {
    start <- max(1, i - 1)
    end <- min(length(x), i + 1)
    change(x[start:end])
  })

  df$diff[!is.na(df$nethh)] <- d

  if (power10) {
    df$power10[!is.na(df$nethh)] <- get_power10(x)
  }

  df
}

#' find similar cases for imputation
#' @param target target row data
#' @param data full data frame
#' @param value_col column to impute
#' @param key_cols columns to match on
#' @param method matching method for key columns
#' @param agg aggregation function for matched values
#' @param min_n minimum matches required
similar_cases <- function(target, data, value_col = "nethh",
                          key_cols = c(
                            "aantalhh", "positiehh", "belbezig",
                            "leeftijd", "oplmet"
                          ),
                          method = max, agg = median, min_n = 1) {
  # start with rows that have non-NA values
  base <- which(!is.na(data[[value_col]]))
  if (length(base) < min_n) {
    return(NA)
  }

  # filter by matching keys
  key_cols <- intersect(key_cols, names(data))
  for (key in key_cols) {
    if (is.na(target[[key]][1])) next
    k <- which(data[[key]] == method(target[[key]]))
    if (length(intersect(base, k)) >= min_n) {
      base <- intersect(base, k)
    }
  }

  if (length(base) == 0) {
    return(NA)
  }
  agg(data[[value_col]][base], na.rm = TRUE)
}

# =============================================================================
# ROBUST OUTLIER DETECTION METHODS
# =============================================================================
# following Tukey (1977), Leys et al. (2013), Breunig et al. (2000)

#' detect univariate outliers using IQR, z-score, or MAD
#' @param x numeric vector
#' @param method "iqr" (Tukey), "zscore", or "mad" (Leys et al.)
#' @param threshold multiplier for bounds (default: 1.5 for IQR, 2.5 for MAD/z)
#' @return logical vector indicating outliers
detect_univariate_outliers <- function(x, method = c("iqr", "mad", "zscore"),
                                       threshold = NULL) {
  method <- match.arg(method)
  x_clean <- x[!is.na(x) & is.finite(x)]

  if (length(x_clean) < 4) {
    return(rep(FALSE, length(x)))
  }

  # default thresholds per method
  if (is.null(threshold)) {
    threshold <- if (method == "iqr") 1.5 else 2.5
  }

  # calculate bounds
  bounds <- switch(method,
    "iqr" = {
      q <- unname(stats::quantile(x_clean, c(0.25, 0.75), na.rm = TRUE))
      iqr <- diff(q)
      c(lower = q[1] - threshold * iqr, upper = q[2] + threshold * iqr)
    },
    "zscore" = {
      m <- mean(x_clean)
      s <- stats::sd(x_clean)
      c(lower = m - threshold * s, upper = m + threshold * s)
    },
    "mad" = {
      med <- stats::median(x_clean)
      mad_val <- stats::mad(x_clean, constant = 1.4826) # scale to match SD
      c(lower = med - threshold * mad_val, upper = med + threshold * mad_val)
    }
  )

  # identify outliers
  is_outlier <- !is.na(x) & is.finite(x) & (x < bounds["lower"] | x > bounds["upper"])

  attr(is_outlier, "bounds") <- bounds
  attr(is_outlier, "method") <- method
  attr(is_outlier, "threshold") <- threshold

  is_outlier
}

#' detect outliers using multiple methods and return consensus
#' @param x numeric vector
#' @param methods character vector of methods to use
#' @param consensus minimum methods that must agree (default: 2)
#' @return list with outlier flags, bounds, and method details
detect_outliers_consensus <- function(x, methods = c("iqr", "mad"),
                                      consensus = 1) {
  results <- lapply(methods, function(m) {
    detect_univariate_outliers(x, method = m)
  })
  names(results) <- methods

  # count how many methods flag each observation
  outlier_counts <- Reduce(`+`, lapply(results, as.integer))

  # consensus outliers
  is_outlier <- outlier_counts >= consensus

  # extract bounds from each method
  bounds_list <- lapply(results, function(r) attr(r, "bounds"))

  list(
    outliers = is_outlier,
    counts = outlier_counts,
    by_method = results,
    bounds = bounds_list
  )
}

#' calculate robust z-score using MAD
#' @param x numeric vector
#' @return vector of robust z-scores (modified z-scores)
robust_zscore <- function(x) {
  med <- stats::median(x, na.rm = TRUE)
  mad_val <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  if (mad_val == 0) mad_val <- stats::sd(x, na.rm = TRUE)
  if (mad_val == 0) {
    return(rep(0, length(x)))
  }
  (x - med) / mad_val
}

# =============================================================================
# INCOME DATA CLEANING
# =============================================================================

#' prepare income data for cleaning
#' @param i income data frame (from i.sav)
#' @param b background data frame (from b.sav)
prepare_income <- function(i, b = NULL) {
  # merge with background if provided
  if (!is.null(b)) {
    # align the background wave key to the module wave key. the analysis modules
    # are adapted to an annual wavenr (wave_year - 2007 via lissr_year_to_wavenr),
    # whereas the LISS background file keys months as YYYYMM. routing the month
    # through wave_to_wavenr() yields a monthly index on a different scale, so a
    # join on wavenr attaches an arbitrary background month to each annual wave
    # (e.g. an income wave of 2008 receives the 2007-11 record). derive the
    # background wavenr on the same annual scale instead, then keep one background
    # row per person-year (latest month within the year) so the join stays
    # one-to-one. a single-month background file (no `wave` column) is left as-is
    # and joined on nomem_encr only, matching the lissr-recommended attachment.
    if ("wave" %in% names(b) && !"wavenr" %in% names(b)) {
      b$wavenr <- as.integer(b$wave %/% 100L - 2007L)
      b <- b[order(b$nomem_encr, b$wavenr, b$wave), , drop = FALSE]
      b <- b[!duplicated(b[c("nomem_encr", "wavenr")], fromLast = TRUE), ,
             drop = FALSE]
    }

    # lissr: join demographics on nomem_encr (and the aligned annual wavenr when
    # the background is wave-stamped), never on nohouse_encr, because the household
    # id is not a stable person-level key and changes with household composition.
    possible_keys <- c("nomem_encr", "wavenr")
    join_keys <- intersect(possible_keys, intersect(names(i), names(b)))

    if (length(join_keys) >= 1) {
      # add only the background columns not already present in i
      b_cols <- c(join_keys, setdiff(names(b), names(i)))
      i <- dplyr::left_join(i, b[b_cols], by = join_keys)
    }
  }

  # ensure required columns exist
  required <- c("nomem_encr", "nohouse_encr", "wavenr", "nethh")
  if (!all(required %in% names(i))) {
    # try to identify columns by pattern
    if (!"nethh" %in% names(i)) {
      # ci00a339 = total net household income
      nethh_col <- grep("ci00a339|nethh$", names(i), value = TRUE)[1]
      if (!is.na(nethh_col)) names(i)[names(i) == nethh_col] <- "nethh"
    }
  }

  # convert to numeric and take absolute values
  income_cols <- grep("net|brut|nethh", names(i), value = TRUE)
  for (col in income_cols) {
    i[[col]] <- abs(as.numeric(sjlabelled::remove_all_labels(i[[col]])))
  }

  # initialize cleaning columns if not present
  if (!"is_na" %in% names(i)) i$is_na <- NA
  if (!"user_na" %in% names(i)) i$user_na <- NA
  if (!"diff" %in% names(i)) i$diff <- 0
  if (!"outlier" %in% names(i)) i$outlier <- 0
  if (!"power10" %in% names(i)) i$power10 <- 0
  if (!"valid_hh" %in% names(i)) i$valid_hh <- 0

  # map income categories to bounds
  if ("nethh_min" %in% names(i)) {
    # if nethh_min contains category codes (1-7), convert to bounds
    if (max(rna(i$nethh_min)) <= 7) {
      i$nethh_max <- income_bounds$upper[match(i$nethh_min, 1:7)]
      i$nethh_min <- income_bounds$lower[match(i$nethh_min, 1:7)]
    }
  }

  # bound columns may be absent (continuous income such as ci00a339 has no
  # category brackets); create as NA so the downstream bound checks no-op
  if (!"nethh_min" %in% names(i)) i$nethh_min <- NA_real_
  if (!"nethh_max" %in% names(i)) i$nethh_max <- NA_real_

  # remove labels for faster computation
  i <- as.data.frame(sjlabelled::remove_all_labels(i))

  # birth year is time-invariant; carry it within person and backfill age for any
  # wave-year the (time-bounded) background file does not cover, so the income
  # imputation's age matching does not silently drop late waves. age = wave_year -
  # birth_year, with wave_year = wavenr + 2007 (inverse of lissr_year_to_wavenr).
  if ("gebjaar" %in% names(i) && "wavenr" %in% names(i)) {
    gp <- stats::ave(
      i$gebjaar, i$nomem_encr,
      FUN = function(g) {
        v <- g[!is.na(g)]
        if (length(v)) v[1] else NA_real_
      }
    )
    if ("leeftijd" %in% names(i)) {
      yr <- i$wavenr + 2007L
      fill <- is.na(i$leeftijd) & !is.na(gp)
      i$leeftijd[fill] <- yr[fill] - gp[fill]
    }
  }

  i
}

#' detect potential outliers in income data
#' @param income prepared income data frame
detect_income_outliers <- function(income) {
  # mark invalid income category codes
  if ("nethh_min" %in% names(income)) {
    invalid_cat <- which(income$nethh_min < 0 |
      (!is.na(income$nethh_min) & income$nethh_min > 120000))
    income <- dummy_na(income, invalid_cat, "nethh_min", "user_na")
    income <- dummy_na(income, invalid_cat, "nethh_max", "user_na")
  }

  # find potential outliers based on multiple criteria
  outlier_idx <- which(with(
    income,
    # very low values
    (!is.na(nethh) & nethh < 10) |
      # low values when individual income is present
      (!is.na(nethh) & !is.na(brutoink) & nethh <= 100) |
      # household income suspiciously close to individual income
      (!is.na(nethh) & nethh < 10000 & (
        (!is.na(nettoink) & abs(nethh - nettoink) < 100) |
          (!is.na(brutoink) & abs(nethh - brutoink) < 100)
      ))
  ))

  # store original values and mark as NA
  income$is_na[outlier_idx] <- income$nethh[outlier_idx]
  income$nethh[outlier_idx] <- NA

  income
}

#' clean household income with outlier detection and imputation
#'
#' constraints:
#' - absolute cap on income (default 150k EUR for Netherlands)
#' - corrections must respect nethh_min/nethh_max category bounds when available
#' - scaling corrections only applied if result is plausible
#'
#' @param income prepared income data frame
#' @param verbose print progress
#' @param income_cap maximum plausible annual household income
#' @param min_income minimum plausible annual household income
clean_household_income <- function(income, verbose = TRUE,
                                   income_cap = 150000, min_income = 8000) {
  households <- unique(income$nohouse_encr)
  n_hh <- length(households)

  if (verbose) message("cleaning income for ", n_hh, " households...")

  for (h in seq_along(households)) {
    hh_id <- households[h]
    hh_idx <- which(income$nohouse_encr == hh_id)
    df <- income[hh_idx, ]

    # DECISION POINT 1: skip if insufficient data
    if (length(rna(df$nethh)) < 2) next

    # calculate lag differences and power of 10
    df <- lag_diff(df)

    # DECISION POINT 2: determine household's typical income range
    # use category bounds if available, otherwise use observed values
    hh_min_bound <- min(rna(c(df$nethh_min, df$nethh)), na.rm = TRUE)
    hh_max_bound <- max(rna(c(df$nethh_max, df$nethh)), na.rm = TRUE)

    # apply absolute cap
    hh_max_bound <- min(hh_max_bound, income_cap, na.rm = TRUE)
    if (is.infinite(hh_min_bound)) hh_min_bound <- min_income
    if (is.infinite(hh_max_bound)) hh_max_bound <- income_cap

    # DECISION POINT 3: scale correction for systematically low values
    # ONLY if: (a) all values < 10k, (b) scaling would fall within bounds
    mode_power10 <- stat_mode(df$power10[df$power10 > 0])

    if (all(rna(df$power10) < 10000) && length(rna(df$power10)) > 0) {
      low_idx <- which(df$power10 > 0 & !is.na(df$nethh))

      for (idx in low_idx) {
        original_val <- df$nethh[idx]
        row_max <- df$nethh_max[idx]
        row_min <- df$nethh_min[idx]

        # try scaling by 10
        scaled_val <- original_val * 10

        # CONSTRAINT: only scale if result is plausible
        scale_ok <- FALSE
        if (!is.na(row_max) && !is.na(row_min)) {
          # category bounds available: scaled value must fall within
          scale_ok <- scaled_val >= row_min && scaled_val <= row_max
        } else if (scaled_val <= income_cap && scaled_val >= min_income) {
          # no bounds: check against absolute limits and consistency
          # also check if scaled value is more consistent with other waves
          other_vals <- rna(df$nethh[-idx])
          if (length(other_vals) > 0) {
            med_other <- median(other_vals)
            # scale if it brings value closer to other observations
            scale_ok <- abs(scaled_val - med_other) < abs(original_val - med_other)
          }
        }

        if (scale_ok) {
          df$outlier[idx] <- original_val
          df$nethh[idx] <- scaled_val
        }
      }
    }

    # recalculate after potential scaling
    df <- lag_diff(df)

    # DECISION POINT 4: iterative outlier detection and imputation
    # uses robust methods: IQR (Tukey 1977), MAD (Leys et al. 2013)
    max_iter <- nrow(df)
    iter <- 0
    corrected_rows <- c() # track to avoid infinite loops

    while (iter < max_iter) {
      iter <- iter + 1
      df <- lag_diff(df)

      mode_power10 <- stat_mode(df$power10[df$power10 > 0])
      if (is.na(mode_power10)) mode_power10 <- 10000

      # OUTLIER DETECTION: hierarchical criteria (most to least specific)
      err_idx <- integer(0)
      detection_method <- NA_character_

      # CRITERION 1: scale error (high lag diff + different power of 10)
      # most specific - indicates decimal point error
      crit1 <- which(df$diff >= 0.6 &
        df$power10 > 0 &
        df$power10 != mode_power10 &
        !(seq_len(nrow(df)) %in% corrected_rows))
      if (length(crit1) > 0) {
        err_idx <- crit1
        detection_method <- "scale_error"
      }

      # CRITERION 2: exceeds category bounds (data entry error)
      if (length(err_idx) == 0) {
        for (row in seq_len(nrow(df))) {
          if (row %in% corrected_rows) next
          val <- df$nethh[row]
          max_b <- df$nethh_max[row]
          min_b <- df$nethh_min[row]
          # flag if >50% above upper bound or <50% below lower bound
          if (!is.na(val) && !is.na(max_b) && val > max_b * 1.5) {
            err_idx <- c(err_idx, row)
            detection_method <- "exceeds_category"
          } else if (!is.na(val) && !is.na(min_b) && min_b > 0 && val < min_b * 0.5) {
            err_idx <- c(err_idx, row)
            detection_method <- "below_category"
          }
        }
      }

      # CRITERION 3: exceeds absolute cap (implausible value)
      if (length(err_idx) == 0) {
        over_cap <- which(df$nethh > income_cap &
          !(seq_len(nrow(df)) %in% corrected_rows))
        if (length(over_cap) > 0) {
          err_idx <- over_cap
          detection_method <- "exceeds_cap"
        }
      }

      # CRITERION 4: robust univariate detection (IQR + MAD consensus)
      # only applied if enough data points and no specific errors found
      if (length(err_idx) == 0 && length(rna(df$nethh)) >= 4) {
        # use both IQR and MAD for robust detection
        outlier_check <- detect_outliers_consensus(
          df$nethh,
          methods = c("iqr", "mad"),
          consensus = 1 # flagged by at least one method
        )

        # require high lag difference as additional confirmation
        robust_out <- which(outlier_check$outliers &
          df$diff >= 0.4 &
          !(seq_len(nrow(df)) %in% corrected_rows))

        if (length(robust_out) > 0) {
          err_idx <- robust_out
          detection_method <- "robust_univariate"
        }
      }

      # CRITERION 5: extreme robust z-score (MAD-based)
      # catches outliers that passed IQR but have extreme deviation
      if (length(err_idx) == 0 && length(rna(df$nethh)) >= 3) {
        z_robust <- robust_zscore(df$nethh)
        extreme_z <- which(abs(z_robust) > 3 &
          !(seq_len(nrow(df)) %in% corrected_rows))
        if (length(extreme_z) > 0) {
          err_idx <- extreme_z
          detection_method <- "extreme_zscore"
        }
      }

      if (length(err_idx) == 0) break

      # select worst outlier based on detection method
      if (detection_method %in% c("scale_error", "robust_univariate")) {
        # for scale/statistical errors: pick highest lag difference
        i <- err_idx[which.max(df$diff[err_idx])]
      } else if (detection_method == "extreme_zscore") {
        # for z-score: pick most extreme
        z_vals <- abs(robust_zscore(df$nethh)[err_idx])
        i <- err_idx[which.max(z_vals)]
      } else {
        # for bound violations: pick largest deviation
        i <- err_idx[1]
      }
      corrected_rows <- c(corrected_rows, i)

      # IMPUTATION: generate candidates constrained by bounds
      original_val <- df$nethh[i]
      row_min <- df$nethh_min[i]
      row_max <- df$nethh_max[i]

      # determine valid range for this row
      valid_min <- if (!is.na(row_min)) row_min else min_income
      valid_max <- if (!is.na(row_max)) row_max else income_cap
      valid_max <- min(valid_max, income_cap)

      # generate candidate values
      candidates <- c()

      # candidate 1: other household values (most reliable)
      other_vals <- rna(df$nethh[-i])
      if (length(other_vals) > 0) {
        candidates <- c(candidates, median(other_vals), mean(other_vals))
      }

      # candidate 2: category midpoint
      if (!is.na(row_min) && !is.na(row_max)) {
        candidates <- c(candidates, (row_min + row_max) / 2)
      }

      # candidate 3: scale correction (if original was misscaled)
      if (!is.na(original_val) && original_val > 0) {
        # try dividing by 10, 100
        for (divisor in c(10, 100)) {
          scaled <- original_val / divisor
          if (scaled >= valid_min && scaled <= valid_max) {
            candidates <- c(candidates, scaled)
          }
        }
        # try multiplying by 10 (only if very low)
        if (original_val < 1000) {
          scaled <- original_val * 10
          if (scaled >= valid_min && scaled <= valid_max) {
            candidates <- c(candidates, scaled)
          }
        }
      }

      # candidate 4: Kalman smoothing (if enough data)
      imp <- tryCatch(
        {
          df_tmp <- df$nethh
          df_tmp[i] <- NA
          if (length(rna(df_tmp)) >= 2) {
            imputeTS::na_ma(df_tmp, k = 2, weighting = "linear")[i]
          } else {
            NA
          }
        },
        error = function(e) NA
      )
      if (!is.na(imp)) candidates <- c(candidates, imp)

      # candidate 5: similar cases
      sim_val <- tryCatch(
        similar_cases(df[i, ], income, "nethh"),
        error = function(e) NA
      )
      if (!is.na(sim_val)) candidates <- c(candidates, sim_val)

      # CONSTRAINT: filter candidates to valid range
      candidates <- unique(candidates)
      candidates <- candidates[!is.na(candidates) &
        candidates >= valid_min &
        candidates <= valid_max]

      if (length(candidates) == 0) {
        # fallback: use midpoint of valid range
        candidates <- (valid_min + valid_max) / 2
      }

      # select best candidate (closest to household median)
      hh_median <- median(rna(df$nethh[-i]), na.rm = TRUE)
      if (is.na(hh_median)) hh_median <- (valid_min + valid_max) / 2

      best_val <- candidates[which.min(abs(candidates - hh_median))]

      # APPLY CORRECTION
      if (df$outlier[i] == 0 || is.na(df$outlier[i])) {
        df$outlier[i] <- original_val
      }
      df$nethh[i] <- best_val
    }

    # WRITE BACK: only if corrections were made
    if (any(df$outlier > 0, na.rm = TRUE)) {
      income[hh_idx, c("nethh", "diff", "outlier", "power10")] <-
        df[, c("nethh", "diff", "outlier", "power10")]
    }

    if (verbose && h %% 1000 == 0) {
      message("  processed ", h, "/", n_hh, " households")
    }
  }

  if (verbose) {
    n_corrected <- sum(income$outlier > 0, na.rm = TRUE)
    n_with_nethh <- sum(!is.na(income$nethh))
    message("income cleaning complete")
    message("  ", n_corrected, " values corrected out of ", n_with_nethh, " non-NA nethh values")

    # distribution check
    if (n_with_nethh > 0) {
      over_cap <- sum(income$nethh > income_cap, na.rm = TRUE)
      under_min <- sum(income$nethh < min_income & income$nethh > 0, na.rm = TRUE)
      if (over_cap > 0) message("  warning: ", over_cap, " values still exceed cap")
      if (under_min > 0) message("  note: ", under_min, " values below ", min_income)
    }
  }
  income
}

#' detect dataset-level outliers after household cleaning
#' uses robust methods on the full distribution
#' @param income cleaned income data frame
#' @param methods univariate methods to use
#' @param verbose print summary
detect_dataset_outliers <- function(income,
                                    methods = c("iqr", "mad"),
                                    verbose = TRUE) {
  if (!"dataset_outlier" %in% names(income)) {
    income$dataset_outlier <- NA_character_
  }

  # UNIVARIATE: detect on nethh using multiple methods
  if (verbose) message("detecting dataset-level univariate outliers...")

  uni_results <- detect_outliers_consensus(income$nethh,
    methods = methods,
    consensus = 1
  )

  uni_outlier_idx <- which(uni_results$outliers & is.na(income$dataset_outlier))
  if (length(uni_outlier_idx) > 0) {
    income$dataset_outlier[uni_outlier_idx] <- paste0(
      "univariate_",
      sapply(uni_outlier_idx, function(i) {
        flagged_by <- names(which(sapply(uni_results$by_method, `[`, i)))
        paste(flagged_by, collapse = "+")
      })
    )
  }

  if (verbose) {
    message("  univariate outliers detected: ", length(uni_outlier_idx))
    for (m in names(uni_results$bounds)) {
      b <- uni_results$bounds[[m]]
      message(
        "    ", toupper(m), " bounds: [",
        round(b["lower"]), ", ", round(b["upper"]), "]"
      )
    }
  }

  # SUMMARY
  if (verbose) {
    n_flagged <- sum(!is.na(income$dataset_outlier))
    message("total dataset-level outliers flagged: ", n_flagged)
  }

  income
}

#' full income cleaning pipeline
#' @param i income data
#' @param b background data (optional)
#' @param verbose print progress
#' @param income_cap maximum plausible annual household income (default 150k)
#' @param min_income minimum plausible annual household income (default 8k)
#' @param dataset_outliers also run dataset-level outlier detection
clean_income <- function(i, b = NULL, verbose = TRUE,
                         income_cap = 150000, min_income = 8000,
                         dataset_outliers = TRUE) {
  income <- prepare_income(i, b)
  income <- detect_income_outliers(income)
  income <- clean_household_income(income, verbose, income_cap, min_income)

  # hard plausibility ceiling: values that survive household correction above
  # income_cap (or non-positive) are unrecoverable and set to NA
  bad <- !is.na(income$nethh) & (income$nethh <= 0 | income$nethh > income_cap)
  if (verbose && sum(bad) > 0) message("  hard cap: ", sum(bad), " implausible nethh values set to NA")
  income$nethh[bad] <- NA

  if (dataset_outliers) {
    income <- detect_dataset_outliers(income, verbose = verbose)
  }

  income
}

#' summarize outlier detection results
#' @param income cleaned income data frame
#' @return list with summary statistics
summarize_outliers <- function(income) {
  n_total <- nrow(income)
  n_nethh <- sum(!is.na(income$nethh))

  # household-level corrections
  n_hh_corrected <- sum(income$outlier > 0, na.rm = TRUE)

  # dataset-level outliers
  n_dataset_outliers <- if ("dataset_outlier" %in% names(income)) {
    sum(!is.na(income$dataset_outlier))
  } else {
    0
  }

  # distribution summary
  nethh_vals <- rna(income$nethh)
  dist_summary <- if (length(nethh_vals) > 0) {
    list(
      min = min(nethh_vals),
      q1 = stats::quantile(nethh_vals, 0.25),
      median = stats::median(nethh_vals),
      mean = mean(nethh_vals),
      q3 = stats::quantile(nethh_vals, 0.75),
      max = max(nethh_vals),
      sd = stats::sd(nethh_vals),
      mad = stats::mad(nethh_vals)
    )
  } else {
    NULL
  }

  list(
    n_total = n_total,
    n_with_income = n_nethh,
    n_household_corrections = n_hh_corrected,
    n_dataset_outliers = n_dataset_outliers,
    pct_corrected = round(100 * n_hh_corrected / n_nethh, 2),
    distribution = dist_summary
  )
}

#' print outlier summary
#' @param summary output from summarize_outliers()
print_outlier_summary <- function(summary) {
  cat("\n=== Income Outlier Detection Summary ===\n")
  cat("Total observations:", summary$n_total, "\n")
  cat("Observations with income:", summary$n_with_income, "\n")
  cat(
    "Household-level corrections:", summary$n_household_corrections,
    "(", summary$pct_corrected, "%)\n"
  )
  cat("Dataset-level outliers flagged:", summary$n_dataset_outliers, "\n")

  if (!is.null(summary$distribution)) {
    cat("\nIncome Distribution (after cleaning):\n")
    d <- summary$distribution
    cat("  Min:", round(d$min), "\n")
    cat("  Q1:", round(d$q1), "\n")
    cat("  Median:", round(d$median), "\n")
    cat("  Mean:", round(d$mean), "\n")
    cat("  Q3:", round(d$q3), "\n")
    cat("  Max:", round(d$max), "\n")
    cat("  SD:", round(d$sd), "\n")
    cat("  MAD:", round(d$mad), "\n")
  }
  cat("========================================\n")
}

