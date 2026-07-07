# 03-ses.R
# socio-economic status construction (education/income recodes).

# =============================================================================
# SES CONSTRUCTION
# =============================================================================

#' recode education to 3-level SES
#' @param oplmet numeric education codes
recode_education <- function(oplmet) {
  dplyr::case_when(
    oplmet %in% c(1, 2, 8, 9) ~ 1L, # low
    oplmet %in% c(3, 4) ~ 2L, # middle
    oplmet %in% c(5, 6) ~ 3L, # high
    TRUE ~ NA_integer_
  )
}

#' recode income to 3-level SES
#' @param nethh net household income
#' @param low_cutoff upper bound for low SES
#' @param high_cutoff lower bound for high SES
recode_income <- function(nethh, low_cutoff = 24000, high_cutoff = 48000) {
  dplyr::case_when(
    nethh < low_cutoff ~ 1L,
    nethh >= high_cutoff ~ 3L,
    TRUE ~ 2L
  )
}

#' equivalize household income for household size and composition
#' returns nethh divided by an equivalence factor. with both household size and
#' children counts available, the OECD-modified scale is used
#' (1 + 0.5 * (adults - 1) + 0.3 * children); with household size only, the
#' square-root scale is the fallback; with neither, income is returned
#' unchanged with a warning so the caller can see equivalization did not apply.
#' @param nethh net household income
#' @param hh_size number of household members (LISS background: aantalhh)
#' @param n_children number of children living at home (LISS background: aantalki)
equivalize_income <- function(nethh, hh_size = NULL, n_children = NULL) {
  if (is.null(hh_size) || all(is.na(hh_size))) {
    warning("no household size available; income not equivalized")
    return(nethh)
  }
  hh <- pmax(as.numeric(hh_size), 1)
  if (!is.null(n_children) && !all(is.na(n_children))) {
    kids <- pmin(pmax(as.numeric(n_children), 0), hh - 1)
    kids[is.na(kids)] <- 0
    adults <- pmax(hh - kids, 1)
    factor_eq <- 1 + 0.5 * (adults - 1) + 0.3 * kids
  } else {
    factor_eq <- sqrt(hh)
  }
  # rows without a usable size keep a factor of 1 (unequivalized) rather than NA
  factor_eq[is.na(factor_eq)] <- 1
  nethh / factor_eq
}

#' construct SES variable with multiple coding options
#' @param data data frame with oplmet and nethh columns
#' @param method "composite", "education", "income", "binary", or "continuous"
#' @param binary_var for binary method: "income", "education", or "composite"
#' @param equivalize for the income-based binary split: divide nethh by an
#'   equivalence factor (OECD-modified when household size and children counts
#'   are present, square-root of household size otherwise) before splitting
#' @param split_within_year for the income-based binary split: classify each
#'   person-wave against the median of its own calendar year rather than the
#'   pooled median, removing the nominal income growth confound on a
#'   calendar-staggered panel
#' @param hh_size_col column holding household size (default "aantalhh")
#' @param children_col column holding children living at home (default "aantalki")
construct_ses <- function(data, method = c(
                            "binary", "composite", "education",
                            "income", "continuous"
                          ),
                          binary_var = "income",
                          equivalize = FALSE,
                          split_within_year = FALSE,
                          hh_size_col = "aantalhh",
                          children_col = "aantalki") {
  # NOTE: binary is now the DEFAULT - recommended for balanced multi-group analyses
  method <- match.arg(method)

  # compute component scores
  if ("oplmet" %in% names(data)) {
    data$edu_ses <- recode_education(data$oplmet)
  }

  if ("nethh" %in% names(data)) {
    data$inc_ses <- recode_income(data$nethh)
    data$inc_log <- log(pmax(data$nethh, 1))
  }

  if (method == "education") {
    data$ses <- factor(data$edu_ses,
      levels = 1:3,
      labels = c("low", "middle", "high")
    )
  } else if (method == "income") {
    data$ses <- factor(data$inc_ses,
      levels = 1:3,
      labels = c("low", "middle", "high")
    )
  } else if (method == "composite") {
    # composite: average of education and income, rounded
    data$ses <- round((data$edu_ses + data$inc_ses) / 2)
    data$ses <- factor(data$ses,
      levels = 1:3,
      labels = c("low", "middle", "high")
    )
  } else if (method == "binary") {
    # binary split for better power in multi-group models
    if (binary_var == "income" && "nethh" %in% names(data)) {
      # income basis: raw nethh by default; equivalized when requested
      inc_basis <- data$nethh
      if (isTRUE(equivalize)) {
        hh_size <- if (hh_size_col %in% names(data)) data[[hh_size_col]] else NULL
        n_children <- if (children_col %in% names(data)) data[[children_col]] else NULL
        inc_basis <- equivalize_income(data$nethh, hh_size, n_children)
        data$nethh_eq <- inc_basis
        message("  income equivalized (", 
          if (!is.null(n_children)) "OECD-modified scale" else "square-root scale", ")")
      }

      if (isTRUE(split_within_year) && "wavenr" %in% names(data)) {
        # per-year median: classify each person-wave against the median of its
        # own calendar year (wave_year = wavenr + 2007), so classification does
        # not encode nominal income growth across the panel period
        year <- data$wavenr + 2007L
        med_by_year <- stats::ave(inc_basis, year,
          FUN = function(x) stats::median(x, na.rm = TRUE)
        )
        data$ses <- factor(ifelse(inc_basis < med_by_year, "low", "high"),
          levels = c("low", "high")
        )
        yr_meds <- tapply(inc_basis, year, stats::median, na.rm = TRUE)
        yr_meds <- yr_meds[!is.na(yr_meds)]
        message(
          "  binary SES (within-year income median split): medians ",
          round(min(yr_meds)), "-", round(max(yr_meds)),
          " across ", length(yr_meds), " years"
        )
      } else {
        med_inc <- stats::median(inc_basis, na.rm = TRUE)
        data$ses <- factor(ifelse(inc_basis < med_inc, "low", "high"),
          levels = c("low", "high")
        )
        message("  binary SES (income median split): ", round(med_inc))
      }
    } else if (binary_var == "education" && "edu_ses" %in% names(data)) {
      # low/middle vs high education
      data$ses <- factor(ifelse(data$edu_ses < 3, "low", "high"),
        levels = c("low", "high")
      )
      message("  binary SES (education: low/middle vs high)")
    } else if (binary_var == "composite") {
      # standardized composite score, then median split
      if ("edu_ses" %in% names(data) && "inc_ses" %in% names(data)) {
        edu_z <- scale(data$edu_ses)[, 1]
        inc_z <- scale(data$inc_ses)[, 1]
        ses_score <- rowMeans(cbind(edu_z, inc_z), na.rm = TRUE)
        data$ses <- factor(ifelse(ses_score < 0, "low", "high"),
          levels = c("low", "high")
        )
        message("  binary SES (composite median split)")
      }
    }
  } else if (method == "continuous") {
    # continuous SES score for use as covariate
    if ("edu_ses" %in% names(data) && "nethh" %in% names(data)) {
      edu_z <- scale(data$edu_ses)[, 1]
      inc_z <- scale(log(pmax(data$nethh, 1)))[, 1]
      data$ses_cont <- rowMeans(cbind(edu_z, inc_z), na.rm = TRUE)
      # also create categorical for reference
      data$ses <- cut(data$ses_cont,
        breaks = c(-Inf, -0.5, 0.5, Inf),
        labels = c("low", "middle", "high")
      )
      message("  continuous SES created (ses_cont)")
    }
  }

  data
}

#' report classification stability of the binary income split
#' shows how many persons ever cross the median across their observed waves;
#' modal person-level assignment (used downstream) discards this mobility, so
#' the share of movers belongs in the paper as a transparency statistic.
#' @param data person-wave data with nomem_encr and ses columns
describe_ses_stability <- function(data) {
  if (!all(c("nomem_encr", "ses") %in% names(data))) {
    message("describe_ses_stability: needs nomem_encr and ses columns")
    return(invisible(NULL))
  }
  d <- data[!is.na(data$ses), c("nomem_encr", "ses")]
  n_levels <- tapply(d$ses, d$nomem_encr, function(x) length(unique(x)))
  n_persons <- length(n_levels)
  n_movers <- sum(n_levels > 1)
  message(
    "  SES classification stability: ", n_persons - n_movers, " of ",
    n_persons, " classified persons stable (",
    round(100 * (n_persons - n_movers) / n_persons, 1), "%); ",
    n_movers, " ever cross the split (",
    round(100 * n_movers / n_persons, 1), "%)"
  )
  invisible(list(n_persons = n_persons, n_movers = n_movers))
}
