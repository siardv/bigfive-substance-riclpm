# 07-assembly-reshape.R
# merge sources into the long analysis dataset, reshape to wide for
# RI-CLPM, and pairwise-coverage / variance diagnostics.

# =============================================================================
# DATA MERGING FOR RI-CLPM
# =============================================================================

#' prepare analysis dataset by merging all data sources
#' NOTE: LISS modules are collected at different calendar times within a panel wave
#' - Personality (cp): typically May/June
#' - Health (ch): typically November
#' - Income (ci): monthly
#' We merge on nomem_encr + wavenr (panel wave number)
#' @param liss list with p (personality), h (health), i_clean (income), b (background)
#' @param waves waves to include (NULL = all)
#' @param substance_method scoring method for substance use ("sum" recommended, "any", "mean")
#' @param substance_items substance item variables to score over (default: the
#'   full six-item battery; pass illicit_five_vars for the sedatives-excluded run)
#' @param substance_min_valid minimum answered items for a non-NA substance score
#' @param include_parcels also carry three balanced item parcels per Big Five
#'   dimension ({dim}_p1..{dim}_p3), required by the multiple-indicator RI-CLPM
prepare_analysis_data <- function(liss, waves = NULL, substance_method = "sum",
                                  substance_items = drugs_use_vars,
                                  substance_min_valid = 1,
                                  include_parcels = FALSE) {
  message("preparing analysis dataset...")

  # 1. Big Five scores
  b5 <- compute_b5_scores(liss$p, min_valid = 8)
  message("  B5 scores: ", nrow(b5), " observations")

  # ensure wavenr exists
  if (!"wavenr" %in% names(b5) && "wave" %in% names(b5)) {
    # derive wavenr from wave (YYYYMM -> sequential number)
    b5$wavenr <- wave_to_wavenr(b5$wave)
  }

  # 2. Substance use - use specified method
  substance <- process_substance_use(liss$h,
    method = substance_method,
    items = substance_items, min_valid = substance_min_valid
  )
  message("  Substance use: ", nrow(substance), " observations (method=", substance_method, ")")

  # ensure wavenr exists in substance data
  if (!is.null(substance) && nrow(substance) > 0) {
    if (!"wavenr" %in% names(substance) && "wave" %in% names(substance)) {
      substance$wavenr <- wave_to_wavenr(substance$wave)
    }

    # keep only relevant columns for merge
    substance_cols <- c(
      "nomem_encr", "wavenr", "substance_use",
      intersect(substance_items, names(substance))
    )
    substance <- dplyr::select(substance, dplyr::any_of(substance_cols))
    substance <- dplyr::distinct(substance, nomem_encr, wavenr, .keep_all = TRUE)
    message("  Substance distinct by person-wave: ", nrow(substance))
  }

  # 4. SES (from cleaned income data)
  ses_data <- liss$i_clean %>%
    dplyr::select(dplyr::any_of(c(
      "nomem_encr", "wavenr", "wave",
      "ses", "edu_ses", "inc_ses",
      "nethh", "oplmet"
    )))

  # ensure wavenr exists
  if (!"wavenr" %in% names(ses_data) && "wave" %in% names(ses_data)) {
    ses_data$wavenr <- wave_to_wavenr(ses_data$wave)
  }

  # get person-level SES (most common across waves)
  # preserve original factor levels
  orig_ses_levels <- levels(factor(ses_data$ses[!is.na(ses_data$ses)]))

  ses_person <- ses_data %>%
    dplyr::filter(!is.na(ses)) %>%
    dplyr::group_by(nomem_encr) %>%
    dplyr::summarize(
      ses = stat_mode(ses),
      edu_ses = stat_mode(edu_ses),
      inc_ses = stat_mode(inc_ses),
      .groups = "drop"
    )

  # ensure factor levels are preserved
  if (length(orig_ses_levels) > 0) {
    ses_person$ses <- factor(ses_person$ses, levels = orig_ses_levels)
  }

  message("  Persons with SES: ", nrow(ses_person))

  # also get wave-level income for merging
  ses_wave <- ses_data %>%
    dplyr::select(nomem_encr, wavenr, nethh, oplmet) %>%
    dplyr::distinct(nomem_encr, wavenr, .keep_all = TRUE)

  # keep only essential B5 columns
  b5_slim <- b5 %>%
    dplyr::select(dplyr::any_of(c(
      "nomem_encr", "wave", "wavenr",
      "extr", "agre", "cons", "neur", "open"
    )))

  # MERGE: join on nomem_encr + wavenr (NOT wave!)
  analysis_data <- b5_slim

  # optional item parcels for the multiple-indicator RI-CLPM; the b5 base
  # stays the row spine and parcels attach by person-wave
  if (isTRUE(include_parcels)) {
    parcels <- add_b5_parcels(liss$p, n_parcels = 3)
    if (!is.null(parcels)) {
      analysis_data <- dplyr::left_join(analysis_data, parcels,
        by = c("nomem_encr", "wavenr")
      )
      message(
        "  Parcels merged: ",
        sum(grepl("_p[0-9]+$", names(analysis_data))), " parcel columns"
      )
    }
  }

  if (!is.null(substance) && nrow(substance) > 0) {
    analysis_data <- dplyr::left_join(
      analysis_data,
      substance,
      by = c("nomem_encr", "wavenr")
    )
    n_matched <- sum(!is.na(analysis_data$substance_use))
    message("  Merged substance use: ", n_matched, " matched")
  }

  # add wave-level income
  analysis_data <- dplyr::left_join(
    analysis_data,
    ses_wave,
    by = c("nomem_encr", "wavenr")
  )

  # add person-level SES
  analysis_data <- dplyr::left_join(
    analysis_data,
    ses_person,
    by = "nomem_encr"
  )

  # filter waves if specified
  if (!is.null(waves)) {
    analysis_data <- dplyr::filter(analysis_data, wavenr %in% waves)
  }

  # add wave indicator as factor
  analysis_data$wave_f <- factor(analysis_data$wavenr)

  message("\nFinal dataset:")
  message("  observations: ", nrow(analysis_data))
  message("  unique persons: ", length(unique(analysis_data$nomem_encr)))
  message("  waves: ", paste(sort(unique(analysis_data$wavenr)), collapse = ", "))
  message("  substance_use non-NA: ", sum(!is.na(analysis_data$substance_use)))
  message("  SES non-NA: ", sum(!is.na(analysis_data$ses)))

  analysis_data
}

#' reshape data to wide format for RI-CLPM
#' @param data long format analysis data
#' @param vars variables to reshape (in addition to IDs)
#' @param max_waves maximum number of waves to include
#' @param require_substance require non-NA substance use to include person
reshape_for_riclpm <- function(data,
                               vars = c(
                                 "extr", "agre", "cons", "neur", "open",
                                 "substance_use"
                               ),
                               max_waves = 5,
                               require_substance = TRUE) {
  # select relevant columns
  keep_cols <- c("nomem_encr", "wavenr", intersect(vars, names(data)), "ses")
  data_long <- data %>%
    dplyr::select(dplyr::any_of(keep_cols)) %>%
    dplyr::filter(!is.na(wavenr))

  # optionally filter to persons with substance use data
  if (require_substance && "substance_use" %in% names(data_long)) {
    # keep persons who have at least one substance use observation
    persons_with_sub <- data_long %>%
      dplyr::filter(!is.na(substance_use)) %>%
      dplyr::pull(nomem_encr) %>%
      unique()

    data_long <- dplyr::filter(data_long, nomem_encr %in% persons_with_sub)
    message(
      "  filtered to ", length(persons_with_sub),
      " persons with substance use data"
    )
  }

  # create sequential wave number within person (1, 2, 3...)
  data_long <- data_long %>%
    dplyr::group_by(nomem_encr) %>%
    dplyr::arrange(wavenr) %>%
    dplyr::mutate(wave_seq = dplyr::row_number()) %>%
    dplyr::filter(wave_seq <= max_waves) %>%
    dplyr::ungroup()

  # get person-level SES (mode across waves)
  # preserve original factor levels
  orig_ses_levels <- levels(factor(data_long$ses[!is.na(data_long$ses)]))

  ses_person <- data_long %>%
    dplyr::filter(!is.na(ses)) %>%
    dplyr::group_by(nomem_encr) %>%
    dplyr::summarize(
      ses = stat_mode(ses),
      .groups = "drop"
    )

  # ensure factor levels are preserved
  if (length(orig_ses_levels) > 0) {
    ses_person$ses <- factor(ses_person$ses, levels = orig_ses_levels)
  }

  # identify which vars actually have data
  vars_with_data <- vars[vars %in% names(data_long)]

  # reshape to wide
  data_wide <- data_long %>%
    dplyr::select(-dplyr::any_of(c("ses", "wavenr"))) %>%
    tidyr::pivot_wider(
      id_cols = nomem_encr,
      names_from = wave_seq,
      values_from = dplyr::all_of(vars_with_data),
      names_sep = "_"
    ) %>%
    dplyr::left_join(ses_person, by = "nomem_encr")

  # report on data availability
  message(
    "reshaped to wide format: ", nrow(data_wide), " persons, ",
    ncol(data_wide) - 2, " variables"
  )

  # check substance use availability
  sub_cols <- grep("^substance_use_", names(data_wide), value = TRUE)
  if (length(sub_cols) > 0) {
    sub_n <- sapply(sub_cols, function(col) sum(!is.na(data_wide[[col]])))
    message(
      "  substance_use observations per wave: ",
      paste(sub_n, collapse = ", ")
    )
  }

  # check personality availability
  extr_cols <- grep("^extr_", names(data_wide), value = TRUE)
  if (length(extr_cols) > 0) {
    extr_n <- sapply(extr_cols, function(col) sum(!is.na(data_wide[[col]])))
    message(
      "  extr observations per wave: ",
      paste(extr_n, collapse = ", ")
    )
  }

  # check SES
  message(
    "  SES distribution: ",
    paste(names(table(data_wide$ses)), "=", table(data_wide$ses), collapse = ", ")
  )

  data_wide
}

# =============================================================================
# COVERAGE DIAGNOSTICS
# =============================================================================

#' check pairwise coverage in wide format data
#' @param data wide format data
#' @param vars variables to check (patterns like "extr_", "substance_use_")
#' @param threshold minimum acceptable coverage (default 0.10)
check_coverage <- function(data, vars = c("extr", "substance_use"), threshold = 0.10) {
  # identify relevant columns
  all_cols <- c()
  for (v in vars) {
    cols <- grep(paste0("^", v, "_\\d+$"), names(data), value = TRUE)
    all_cols <- c(all_cols, cols)
  }

  if (length(all_cols) == 0) {
    warning("no matching columns found")
    return(NULL)
  }

  # compute pairwise coverage
  n <- nrow(data)
  coverage_matrix <- matrix(NA, length(all_cols), length(all_cols),
    dimnames = list(all_cols, all_cols)
  )

  for (i in seq_along(all_cols)) {
    for (j in seq_along(all_cols)) {
      both_obs <- sum(!is.na(data[[all_cols[i]]]) & !is.na(data[[all_cols[j]]]))
      coverage_matrix[i, j] <- both_obs / n
    }
  }

  # identify low coverage pairs
  low_coverage <- which(coverage_matrix < threshold &
    lower.tri(coverage_matrix), arr.ind = TRUE)

  low_pairs <- if (nrow(low_coverage) > 0) {
    data.frame(
      var1 = all_cols[low_coverage[, 1]],
      var2 = all_cols[low_coverage[, 2]],
      coverage = coverage_matrix[low_coverage]
    )
  } else {
    data.frame(var1 = character(), var2 = character(), coverage = numeric())
  }

  # summary statistics
  diag_coverage <- diag(coverage_matrix)
  off_diag <- coverage_matrix[lower.tri(coverage_matrix)]

  list(
    matrix = coverage_matrix,
    low_pairs = low_pairs,
    n_low = nrow(low_pairs),
    min_coverage = min(off_diag, na.rm = TRUE),
    mean_coverage = mean(off_diag, na.rm = TRUE),
    var_coverage = setNames(diag_coverage, all_cols),
    threshold = threshold
  )
}

#' print coverage summary
print_coverage_summary <- function(coverage) {
  cat("\n=== Pairwise Coverage Diagnostics ===\n")
  cat("Threshold: ", coverage$threshold * 100, "%\n")
  cat("Minimum coverage: ", round(coverage$min_coverage * 100, 1), "%\n")
  cat("Mean coverage: ", round(coverage$mean_coverage * 100, 1), "%\n")
  cat("Pairs below threshold: ", coverage$n_low, "\n")

  if (coverage$n_low > 0) {
    cat("\nLow coverage pairs:\n")
    print(coverage$low_pairs[order(coverage$low_pairs$coverage), ])
    cat("\nWARNING: Low coverage may bias estimates. Consider:\n")
    cat("  - Reducing number of waves\n")
    cat("  - Using more complete data subsets\n")
    cat("  - Using multiple imputation\n")
  } else {
    cat("\nAll pairs meet coverage threshold.\n")
  }
  cat("=====================================\n")
}

#' check variance by group for multi-group models
#' @param data wide format data
#' @param group_var grouping variable
#' @param vars variables to check
check_variance_by_group <- function(data, group_var = "ses", vars = "substance_use") {
  if (!group_var %in% names(data)) {
    warning("group variable '", group_var, "' not found")
    return(NULL)
  }

  data <- data[!is.na(data[[group_var]]), ]
  groups <- unique(data[[group_var]])

  # find relevant columns
  all_cols <- c()
  for (v in vars) {
    cols <- grep(paste0("^", v, "_\\d+$"), names(data), value = TRUE)
    all_cols <- c(all_cols, cols)
  }

  # compute variance by group
  results <- list()
  for (g in groups) {
    grp_data <- data[data[[group_var]] == g, ]
    grp_vars <- sapply(all_cols, function(col) {
      stats::var(grp_data[[col]], na.rm = TRUE)
    })
    results[[as.character(g)]] <- grp_vars
  }

  # identify zero variance
  zero_var <- list()
  for (g in names(results)) {
    zero_cols <- names(results[[g]])[results[[g]] == 0]
    if (length(zero_cols) > 0) {
      zero_var[[g]] <- zero_cols
    }
  }

  list(
    variance_by_group = results,
    zero_variance = zero_var,
    has_zero_variance = length(zero_var) > 0
  )
}



# =============================================================================
# PARCELS AND SELECTION DIAGNOSTICS
# =============================================================================

#' parcel column names for a set of dimensions
#' @param dims Big Five dimension names
#' @param n_parcels parcels per dimension
parcel_var_names <- function(dims = c("extr", "open", "cons"), n_parcels = 3) {
  as.vector(sapply(dims, function(d) paste0(d, "_p", seq_len(n_parcels))))
}

#' build person-wave item parcels for all Big Five dimensions
#' wraps create_b5_parcels() over the long personality frame and returns one
#' distinct row per person-wave with the parcel columns only
#' @param p personality module data (long, item columns present)
#' @param n_parcels parcels per dimension (3 recommended for the
#'   multiple-indicator RI-CLPM)
add_b5_parcels <- function(p, n_parcels = 3) {
  if (is.null(p) || nrow(p) == 0) {
    return(NULL)
  }
  d <- p
  if (!"wavenr" %in% names(d) && "wave" %in% names(d)) {
    d$wavenr <- wave_to_wavenr(d$wave)
  }
  if (!all(c("nomem_encr", "wavenr") %in% names(d))) {
    warning("add_b5_parcels: person-wave keys missing")
    return(NULL)
  }

  for (dim in names(b5_items)) {
    d <- create_b5_parcels(d, dim, n_parcels = n_parcels)
  }

  keep <- c("nomem_encr", "wavenr", parcel_var_names(names(b5_items), n_parcels))
  d %>%
    dplyr::select(dplyr::any_of(keep)) %>%
    dplyr::distinct(nomem_encr, wavenr, .keep_all = TRUE)
}

#' compare selected respondents with the excluded remainder of the eligible pool
#' the eligible pool is every person in the pre-selection long data with at
#' least one observed value on the anchoring variable; comparisons cover the
#' first observed trait scores, ever-any-substance use, and low-SES share.
#' this backs the attrition paragraph accompanying the completers-only design.
#' @param pre_data long analysis data before weasel selection
#' @param selected_ids person ids retained by the selection rule
#' @param anchor variable defining eligibility (default "extr")
#' @param traits trait columns to compare at first observation
describe_selection_attrition <- function(pre_data, selected_ids,
                                         anchor = "extr",
                                         traits = c("extr", "open", "cons")) {
  if (!anchor %in% names(pre_data)) {
    message("describe_selection_attrition: anchor variable not found")
    return(invisible(NULL))
  }

  eligible <- pre_data %>%
    dplyr::filter(!is.na(.data[[anchor]])) %>%
    dplyr::arrange(nomem_encr, wavenr)

  first_obs <- eligible %>%
    dplyr::group_by(nomem_encr) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  person_level <- pre_data %>%
    dplyr::group_by(nomem_encr) %>%
    dplyr::summarize(
      ever_substance = as.integer(any(substance_use > 0, na.rm = TRUE)),
      has_substance = any(!is.na(substance_use)),
      ses_low = as.integer(any(ses == "low", na.rm = TRUE) & !any(ses == "high", na.rm = TRUE)),
      has_ses = any(!is.na(ses)),
      .groups = "drop"
    )

  comp <- first_obs %>%
    dplyr::left_join(person_level, by = "nomem_encr") %>%
    dplyr::mutate(selected = nomem_encr %in% selected_ids)

  summarize_group <- function(dd) {
    out <- c(n = nrow(dd))
    for (tr in intersect(traits, names(dd))) {
      out[paste0(tr, "_mean")] <- round(mean(dd[[tr]], na.rm = TRUE), 3)
    }
    out["ever_substance_pct"] <- round(
      100 * mean(dd$ever_substance[dd$has_substance], na.rm = TRUE), 1
    )
    out["ses_low_pct"] <- round(
      100 * mean(dd$ses_low[dd$has_ses], na.rm = TRUE), 1
    )
    out
  }

  tab <- rbind(
    selected = summarize_group(comp[comp$selected, , drop = FALSE]),
    excluded = summarize_group(comp[!comp$selected, , drop = FALSE])
  )

  message("selection attrition (eligible pool = >=1 observed ", anchor, "):")
  print(tab)

  # standardized mean differences on first-observed traits
  smd <- sapply(intersect(traits, names(comp)), function(tr) {
    a <- comp[[tr]][comp$selected]
    b <- comp[[tr]][!comp$selected]
    sp <- sqrt((stats::var(a, na.rm = TRUE) + stats::var(b, na.rm = TRUE)) / 2)
    round((mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sp, 3)
  })
  message("  standardized mean differences (selected - excluded):")
  print(smd)

  invisible(list(table = tab, smd = smd))
}
