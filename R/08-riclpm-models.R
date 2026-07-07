# 08-riclpm-models.R
# RI-CLPM lavaan syntax builders and single/multigroup/interaction fits.

# =============================================================================
# RI-CLPM MODEL SPECIFICATION
# =============================================================================

#' build RI-CLPM model syntax for personality and substance use
#' based on Hamaker et al. (2015) specification
#'
#' @param personality personality trait name (e.g., "extr")
#' @param outcome outcome variable name (e.g., "substance_use")
#' @param waves number of waves
#' @param constrain_ar constrain autoregressive paths equal across time
#' @param constrain_cl constrain cross-lagged paths equal across time
build_riclpm_model <- function(personality = "extr",
                               outcome = "substance_use",
                               waves = 5,
                               constrain_ar = TRUE,
                               constrain_cl = TRUE) {
  p <- personality
  o <- outcome

  # variable names at each wave
  p_vars <- paste0(p, "_", 1:waves)
  o_vars <- paste0(o, "_", 1:waves)

  # --- RANDOM INTERCEPTS (between-person) ---
  ri_p <- paste0("RI_", p, " =~ ", paste0("1*", p_vars, collapse = " + "))
  ri_o <- paste0("RI_", o, " =~ ", paste0("1*", o_vars, collapse = " + "))

  # --- WITHIN-PERSON CENTERED VARIABLES ---
  wp_p <- sapply(1:waves, function(t) paste0("wp_", p, t, " =~ 1*", p_vars[t]))
  wp_o <- sapply(1:waves, function(t) paste0("wp_", o, t, " =~ 1*", o_vars[t]))

  # --- AUTOREGRESSIVE PATHS ---
  if (constrain_ar) {
    ar_p <- sapply(2:waves, function(t) {
      paste0("wp_", p, t, " ~ a*wp_", p, t - 1)
    })
    ar_o <- sapply(2:waves, function(t) {
      paste0("wp_", o, t, " ~ b*wp_", o, t - 1)
    })
  } else {
    ar_p <- sapply(2:waves, function(t) {
      paste0("wp_", p, t, " ~ wp_", p, t - 1)
    })
    ar_o <- sapply(2:waves, function(t) {
      paste0("wp_", o, t, " ~ wp_", o, t - 1)
    })
  }

  # --- CROSS-LAGGED PATHS ---
  if (constrain_cl) {
    # personality -> outcome (H1)
    cl_p_to_o <- sapply(2:waves, function(t) {
      paste0("wp_", o, t, " ~ c*wp_", p, t - 1)
    })
    # outcome -> personality (H2)
    cl_o_to_p <- sapply(2:waves, function(t) {
      paste0("wp_", p, t, " ~ d*wp_", o, t - 1)
    })
  } else {
    cl_p_to_o <- sapply(2:waves, function(t) {
      paste0("wp_", o, t, " ~ wp_", p, t - 1)
    })
    cl_o_to_p <- sapply(2:waves, function(t) {
      paste0("wp_", p, t, " ~ wp_", o, t - 1)
    })
  }

  # --- COVARIANCES ---
  # random intercept covariance
  ri_cov <- paste0("RI_", p, " ~~ RI_", o)

  # within-person residual covariances at wave 1
  wp_cov1 <- paste0("wp_", p, "1 ~~ wp_", o, "1")

  # residual covariances at subsequent waves
  wp_cov <- sapply(2:waves, function(t) {
    paste0("wp_", p, t, " ~~ wp_", o, t)
  })

  # --- VARIANCE CONSTRAINTS ---
  # fix observed variable residual variances to 0
  obs_var_p <- paste0(p_vars, " ~~ 0*", p_vars)
  obs_var_o <- paste0(o_vars, " ~~ 0*", o_vars)

  # --- ASSEMBLE MODEL ---
  model <- paste(c(
    "# Random intercepts (between-person)",
    ri_p, ri_o,
    "",
    "# Within-person centered variables",
    wp_p, wp_o,
    "",
    "# Autoregressive paths",
    ar_p, ar_o,
    "",
    "# Cross-lagged paths",
    "# Personality -> Outcome (H1)",
    cl_p_to_o,
    "# Outcome -> Personality (H2)",
    cl_o_to_p,
    "",
    "# Covariances",
    ri_cov, wp_cov1, wp_cov,
    "",
    "# Fix observed residual variances to 0",
    obs_var_p, obs_var_o
  ), collapse = "\n")

  model
}

#' fit RI-CLPM model
#' @param data wide format data
#' @param personality personality trait
#' @param outcome outcome variable
#' @param waves number of waves
#' @param estimator "ML" or "MLR" (robust)
#' @param em_iter maximum EM iterations for missing data
#' @param ... additional arguments to lavaan
fit_riclpm <- function(data, personality = "extr", outcome = "substance_use",
                       waves = 5, estimator = "MLR", em_iter = 10000, ...) {
  model <- build_riclpm_model(personality, outcome, waves)

  fit <- tryCatch(
    {
      lavaan::sem(model,
        data = data,
        missing = "fiml",
        estimator = estimator,
        em.h1.iter.max = em_iter,
        ...
      )
    },
    error = function(e) {
      warning("RI-CLPM fitting failed: ", e$message)
      NULL
    }
  )

  list(
    personality = personality,
    outcome = outcome,
    waves = waves,
    model = model,
    fit = fit
  )
}

#' build multi-group RI-CLPM syntax with explicit per-group labels
#' companion to build_riclpm_model(); used only for the SES moderation (H3) fits.
#' single labels are never used (those force cross-group equality AND emit a
#' lavaan note); instead every parameter carries a c(...) label vector with one
#' entry per group, so cross-group equality is imposed deliberately and visibly.
#'
#' variance and covariance components are pooled across groups by default
#' (one shared label per component). pooling estimates these nuisance parameters
#' from both groups jointly, which is what keeps the sparse per-group cells (the
#' thin wave-5 substance data) from driving within-person residual variances
#' negative. the autoregressive and cross-lagged paths are time-invariant within
#' group and are switched between free and equal across groups to build the
#' nested models for the chi-square difference tests.
#'
#' @param personality personality trait name (e.g., "extr")
#' @param outcome outcome variable name (e.g., "substance_use")
#' @param waves number of waves
#' @param groups character vector of group levels (e.g., levels(data$ses))
#' @param ar_free autoregressive paths free across groups (TRUE) or equal (FALSE)
#' @param cl_free cross-lagged paths free across groups (TRUE) or equal (FALSE)
#' @param pool_variances pool RI and within-person variances across groups
#' @param pool_covariances pool RI and within-person covariances across groups
build_riclpm_model_mg <- function(personality = "extr",
                                  outcome = "substance_use",
                                  waves = 5,
                                  groups,
                                  ar_free = TRUE,
                                  cl_free = TRUE,
                                  pool_variances = TRUE,
                                  pool_covariances = TRUE) {
  p <- personality
  o <- outcome
  g <- length(groups)
  gn <- gsub("[^A-Za-z0-9]", "_", as.character(groups)) # safe label suffixes

  p_vars <- paste0(p, "_", 1:waves)
  o_vars <- paste0(o, "_", 1:waves)

  # per-group label vector: free -> one label per group; equal -> one shared label
  glab <- function(base, free) {
    labs <- if (free) paste0(base, "_", gn) else rep(paste0(base, "_eq"), g)
    paste0("c(", paste(labs, collapse = ","), ")*")
  }

  # --- RANDOM INTERCEPTS (between-person), unit loadings ---
  ri_p <- paste0("RI_", p, " =~ ", paste0("1*", p_vars, collapse = " + "))
  ri_o <- paste0("RI_", o, " =~ ", paste0("1*", o_vars, collapse = " + "))

  # --- WITHIN-PERSON CENTERED VARIABLES, unit loadings ---
  wp_p <- sapply(1:waves, function(t) paste0("wp_", p, t, " =~ 1*", p_vars[t]))
  wp_o <- sapply(1:waves, function(t) paste0("wp_", o, t, " =~ 1*", o_vars[t]))

  # --- AUTOREGRESSIVE PATHS (time-invariant within group) ---
  ar_p <- sapply(2:waves, function(t) {
    paste0("wp_", p, t, " ~ ", glab("arp", ar_free), "wp_", p, t - 1)
  })
  ar_o <- sapply(2:waves, function(t) {
    paste0("wp_", o, t, " ~ ", glab("aro", ar_free), "wp_", o, t - 1)
  })

  # --- CROSS-LAGGED PATHS (time-invariant within group) ---
  # personality -> outcome (H1), outcome -> personality (H2)
  cl_p_to_o <- sapply(2:waves, function(t) {
    paste0("wp_", o, t, " ~ ", glab("clpo", cl_free), "wp_", p, t - 1)
  })
  cl_o_to_p <- sapply(2:waves, function(t) {
    paste0("wp_", p, t, " ~ ", glab("clop", cl_free), "wp_", o, t - 1)
  })

  # --- VARIANCE COMPONENTS (pooled across groups by default) ---
  vf <- !pool_variances # free across groups?
  ri_var_p <- paste0("RI_", p, " ~~ ", glab("vrip", vf), "RI_", p)
  ri_var_o <- paste0("RI_", o, " ~~ ", glab("vrio", vf), "RI_", o)
  wp_var_p <- sapply(1:waves, function(t) {
    paste0("wp_", p, t, " ~~ ", glab(paste0("vwp", t), vf), "wp_", p, t)
  })
  wp_var_o <- sapply(1:waves, function(t) {
    paste0("wp_", o, t, " ~~ ", glab(paste0("vwo", t), vf), "wp_", o, t)
  })

  # --- COVARIANCE COMPONENTS (pooled across groups by default) ---
  cf <- !pool_covariances
  ri_cov <- paste0("RI_", p, " ~~ ", glab("cri", cf), "RI_", o)
  wp_cov1 <- paste0("wp_", p, "1 ~~ ", glab("cw1", cf), "wp_", o, "1")
  wp_cov <- sapply(2:waves, function(t) {
    paste0("wp_", p, t, " ~~ ", glab(paste0("cw", t), cf), "wp_", o, t)
  })

  # --- OBSERVED RESIDUALS FIXED TO 0 (identical to single-group spec) ---
  obs_var_p <- paste0(p_vars, " ~~ 0*", p_vars)
  obs_var_o <- paste0(o_vars, " ~~ 0*", o_vars)

  paste(c(
    "# Random intercepts (between-person)",
    ri_p, ri_o,
    "",
    "# Within-person centered variables",
    wp_p, wp_o,
    "",
    "# Autoregressive paths",
    ar_p, ar_o,
    "",
    "# Cross-lagged paths",
    "# Personality -> Outcome (H1)",
    cl_p_to_o,
    "# Outcome -> Personality (H2)",
    cl_o_to_p,
    "",
    "# Variance components (pooled across groups by default)",
    ri_var_p, ri_var_o, wp_var_p, wp_var_o,
    "",
    "# Covariance components (pooled across groups by default)",
    ri_cov, wp_cov1, wp_cov,
    "",
    "# Fix observed residual variances to 0",
    obs_var_p, obs_var_o
  ), collapse = "\n")
}

#' fit multi-group RI-CLPM by SES
#' @param data wide format data with SES variable
#' @param personality personality trait
#' @param outcome outcome variable
#' @param waves number of waves
#' @param estimator "ML" or "MLR" (robust)
#' @param em_iter maximum EM iterations
fit_riclpm_multigroup <- function(data, personality = "extr",
                                  outcome = "substance_use",
                                  waves = 5,
                                  estimator = "MLR",
                                  em_iter = 10000) {
  # ensure SES is factor and remove NA
  data <- data[!is.na(data$ses), ]
  data$ses <- factor(data$ses)

  n_per_group <- table(data$ses)
  n_total <- sum(n_per_group)
  n_groups <- length(n_per_group)

  message("  SES group sizes: ", paste(names(n_per_group), "=", n_per_group, collapse = ", "))

  # compute imbalance metrics
  props <- n_per_group / n_total
  min_prop <- min(props)
  max_prop <- max(props)
  imbalance_ratio <- max_prop / min_prop

  # determine reliability of multi-group test
  reliability_issues <- list()

  # CRITICAL: check for severe imbalance
  if (imbalance_ratio > 10) {
    reliability_issues$severe_imbalance <- TRUE
    warning(
      "\n  SEVERE GROUP IMBALANCE DETECTED (ratio: ",
      round(imbalance_ratio, 1), ":1)",
      "\n  The largest group is ", round(imbalance_ratio, 1),
      "x larger than the smallest.",
      "\n  Multi-group χ² test results are UNRELIABLE.",
      "\n  Recommend: Use binary SES (method='binary') or interaction model."
    )
  } else if (imbalance_ratio > 5) {
    reliability_issues$moderate_imbalance <- TRUE
    warning(
      "  Moderate group imbalance (ratio: ", round(imbalance_ratio, 1),
      ":1). Interpret with caution."
    )
  }

  # check for minimum group size
  if (any(n_per_group < 100)) {
    reliability_issues$small_groups <- TRUE
    small_grps <- names(n_per_group)[n_per_group < 100]
    warning(
      "  Groups with <100 observations: ", paste(small_grps, collapse = ", "),
      "\n  Parameter estimates in small groups may be unstable."
    )
  }

  # check if majority group dominates
  if (max_prop > 0.90) {
    reliability_issues$majority_dominates <- TRUE
    warning(
      "  One group contains >90% of observations (",
      round(max_prop * 100, 1), "%).",
      "\n  Multi-group comparison is essentially testing that group alone."
    )
  }

  # check for zero variance in outcome by group - critical issue
  outcome_vars <- paste0(outcome, "_", 1:waves)
  outcome_vars <- intersect(outcome_vars, names(data))

  zero_var_groups <- c()
  for (g in levels(data$ses)) {
    grp_data <- data[data$ses == g, outcome_vars, drop = FALSE]
    grp_vars <- sapply(grp_data, function(x) stats::var(x, na.rm = TRUE))
    if (any(grp_vars == 0, na.rm = TRUE)) {
      zero_var_groups <- c(zero_var_groups, g)
    }
  }

  if (length(zero_var_groups) > 0) {
    warning(
      "Zero variance in ", outcome, " for SES groups: ",
      paste(zero_var_groups, collapse = ", "),
      "\n  Multi-group model cannot be estimated.",
      "\n  Consider using method='sum' instead of method='any' for more variance."
    )
    return(list(
      personality = personality,
      outcome = outcome,
      fit_configural = NULL,
      fit_constrained = NULL,
      comparison = NULL,
      group_params = NULL,
      group_ns = n_per_group,
      imbalance_ratio = imbalance_ratio,
      reliability_issues = c(reliability_issues, list(zero_variance = TRUE)),
      results_reliable = FALSE,
      error = "zero_variance_in_groups"
    ))
  }

  # three nested multi-group models give MLR-scaled LRTs for H3. variance and
  # covariance components are pooled across groups (one shared label each) so the
  # sparse wave-5 substance cell cannot drive a within-person residual negative;
  # only the regression paths differ between models, via per-group label vectors:
  #   m_free  : autoregressive + cross-lagged free across groups (configural)
  #   m_cleq  : cross-lagged equal across groups, autoregressive free
  #   m_alleq : autoregressive + cross-lagged equal across groups
  groups <- levels(data$ses)

  message(
    "fitting multi-group RI-CLPM for ", personality, " ~ ", outcome,
    " (", length(groups), " SES groups; variance components pooled across groups)"
  )

  m_free <- build_riclpm_model_mg(personality, outcome, waves, groups,
    ar_free = TRUE, cl_free = TRUE
  )
  m_cleq <- build_riclpm_model_mg(personality, outcome, waves, groups,
    ar_free = TRUE, cl_free = FALSE
  )
  m_alleq <- build_riclpm_model_mg(personality, outcome, waves, groups,
    ar_free = FALSE, cl_free = FALSE
  )

  fit_one <- function(model, what) {
    tryCatch(
      lavaan::sem(model,
        data = data, group = "ses",
        missing = "fiml", estimator = estimator,
        em.h1.iter.max = em_iter
      ),
      error = function(e) {
        warning(what, " model failed: ", e$message)
        NULL
      }
    )
  }

  fit_config <- fit_one(m_free, "configural (free)")
  fit_constrain_cl <- fit_one(m_cleq, "cross-lagged equal")
  fit_constrain_all <- fit_one(m_alleq, "all-regressions equal")

  # admissibility check: convergence and non-negative variance estimates
  adm <- function(fit) {
    if (is.null(fit)) {
      return(list(ok = FALSE, converged = FALSE, n_neg = NA_integer_))
    }
    conv <- isTRUE(lavaan::lavInspect(fit, "converged"))
    pe <- lavaan::parameterEstimates(fit)
    nneg <- sum(pe$op == "~~" & pe$lhs == pe$rhs & pe$est < 0, na.rm = TRUE)
    list(ok = conv && nneg == 0, converged = conv, n_neg = nneg)
  }
  a_free <- adm(fit_config)
  a_cl <- adm(fit_constrain_cl)
  a_all <- adm(fit_constrain_all)

  if (!a_free$converged || !a_cl$converged || !a_all$converged) {
    warning("  a multi-group model did not converge")
  }
  nneg_total <- sum(c(a_free$n_neg, a_cl$n_neg, a_all$n_neg), na.rm = TRUE)
  if (nneg_total > 0) {
    warning(
      "  Heywood case detected: ", nneg_total,
      " negative variance estimate(s) across multi-group models"
    )
    reliability_issues$heywood_case <- TRUE
  }

  # targeted H3 test (cross-lagged group-invariance), MLR-scaled via lavaan::anova
  comparison_cl <- if (!is.null(fit_config) && !is.null(fit_constrain_cl)) {
    tryCatch(
      lavaan::anova(fit_constrain_cl, fit_config),
      error = function(e) {
        warning("cross-lagged comparison failed: ", e$message)
        NULL
      }
    )
  } else {
    NULL
  }

  # broad test (all structural regressions group-invariant: autoregressive + cross-lagged)
  comparison <- if (!is.null(fit_config) && !is.null(fit_constrain_all)) {
    tryCatch(
      lavaan::anova(fit_constrain_all, fit_config),
      error = function(e) {
        warning("model comparison failed: ", e$message)
        NULL
      }
    )
  } else {
    NULL
  }

  # group-specific cross-lagged estimates from the configural (free) model
  group_params <- if (!is.null(fit_config)) {
    tryCatch(
      {
        params <- lavaan::parameterEstimates(fit_config, standardized = TRUE)
        params[grepl("^wp_", params$lhs) &
          params$op == "~" &
          grepl("^wp_", params$rhs), ]
      },
      error = function(e) NULL
    )
  } else {
    NULL
  }

  list(
    personality = personality,
    outcome = outcome,
    fit_configural = fit_config,
    fit_constrained = fit_constrain_all,
    fit_constrained_cl = fit_constrain_cl,
    comparison = comparison,
    comparison_cl = comparison_cl,
    group_params = group_params,
    group_ns = n_per_group,
    imbalance_ratio = imbalance_ratio,
    reliability_issues = reliability_issues,
    results_reliable = length(reliability_issues) == 0
  )
}

#' fit RI-CLPM with SES as continuous moderator (interaction approach)
#' alternative to multi-group when groups have insufficient variance
#' @param data wide format data with ses_cont (continuous SES)
#' @param personality personality trait
#' @param outcome outcome variable
#' @param waves number of waves
#' @param estimator "ML" or "MLR"
fit_riclpm_interaction <- function(data, personality = "extr",
                                   outcome = "substance_use",
                                   waves = 5,
                                   estimator = "MLR") {
  # need continuous SES
  if (!"ses_cont" %in% names(data)) {
    # create from factor if needed
    if ("ses" %in% names(data)) {
      data$ses_cont <- as.numeric(factor(data$ses))
      data$ses_cont <- scale(data$ses_cont)[, 1]
    } else {
      warning("no SES variable found")
      return(NULL)
    }
  }

  p <- personality
  o <- outcome

  # create interaction terms (centered)
  for (t in 1:waves) {
    p_var <- paste0(p, "_", t)
    o_var <- paste0(o, "_", t)

    if (p_var %in% names(data)) {
      # center personality
      data[[paste0(p_var, "_c")]] <- scale(data[[p_var]], scale = FALSE)[, 1]
      # interaction: personality x SES
      data[[paste0(p_var, "_x_ses")]] <- data[[paste0(p_var, "_c")]] * data$ses_cont
    }
    if (o_var %in% names(data)) {
      # center outcome
      data[[paste0(o_var, "_c")]] <- scale(data[[o_var]], scale = FALSE)[, 1]
      # interaction: outcome x SES
      data[[paste0(o_var, "_x_ses")]] <- data[[paste0(o_var, "_c")]] * data$ses_cont
    }
  }

  # build model with interaction terms as predictors
  p_vars <- paste0(p, "_", 1:waves)
  o_vars <- paste0(o, "_", 1:waves)

  # base RI-CLPM
  ri_p <- paste0("RI_", p, " =~ ", paste0("1*", p_vars, collapse = " + "))
  ri_o <- paste0("RI_", o, " =~ ", paste0("1*", o_vars, collapse = " + "))

  wp_p <- sapply(1:waves, function(t) paste0("wp_", p, t, " =~ 1*", p_vars[t]))
  wp_o <- sapply(1:waves, function(t) paste0("wp_", o, t, " =~ 1*", o_vars[t]))

  # autoregressive with interaction moderation
  ar_p <- sapply(2:waves, function(t) {
    paste0("wp_", p, t, " ~ wp_", p, t - 1, " + ", p_vars[t - 1], "_x_ses")
  })
  ar_o <- sapply(2:waves, function(t) {
    paste0("wp_", o, t, " ~ wp_", o, t - 1, " + ", o_vars[t - 1], "_x_ses")
  })

  # cross-lagged with interaction moderation (key tests)
  cl_p_to_o <- sapply(2:waves, function(t) {
    paste0("wp_", o, t, " ~ cl_po*wp_", p, t - 1, " + mod_po*", p_vars[t - 1], "_x_ses")
  })
  cl_o_to_p <- sapply(2:waves, function(t) {
    paste0("wp_", p, t, " ~ cl_op*wp_", o, t - 1, " + mod_op*", o_vars[t - 1], "_x_ses")
  })

  # covariances
  ri_cov <- paste0("RI_", p, " ~~ RI_", o)
  wp_cov1 <- paste0("wp_", p, "1 ~~ wp_", o, "1")
  wp_cov <- sapply(2:waves, function(t) paste0("wp_", p, t, " ~~ wp_", o, t))

  # fix observed residual variances
  obs_var_p <- paste0(p_vars, " ~~ 0*", p_vars)
  obs_var_o <- paste0(o_vars, " ~~ 0*", o_vars)

  model <- paste(c(
    "# Random intercepts",
    ri_p, ri_o,
    "# Within-person",
    wp_p, wp_o,
    "# Autoregressive with SES moderation",
    ar_p, ar_o,
    "# Cross-lagged with SES moderation",
    "# cl_po = personality->outcome, mod_po = SES moderation of this path",
    cl_p_to_o,
    "# cl_op = outcome->personality, mod_op = SES moderation of this path",
    cl_o_to_p,
    "# Covariances",
    ri_cov, wp_cov1, wp_cov,
    "# Fix observed variances",
    obs_var_p, obs_var_o
  ), collapse = "\n")

  fit <- tryCatch(
    {
      lavaan::sem(model, data = data, missing = "fiml", estimator = estimator)
    },
    error = function(e) {
      warning("interaction model failed: ", e$message)
      NULL
    }
  )

  # extract moderation parameters
  mod_params <- if (!is.null(fit)) {
    params <- lavaan::parameterEstimates(fit, standardized = TRUE)
    params[params$label %in% c("cl_po", "cl_op", "mod_po", "mod_op"), ]
  } else {
    NULL
  }

  list(
    personality = personality,
    outcome = outcome,
    model = model,
    fit = fit,
    mod_params = mod_params
  )
}


# =============================================================================
# MULTIPLE-INDICATOR RI-CLPM (PARCELS)
# =============================================================================

#' build a multiple-indicator RI-CLPM with parcel indicators for personality
#' each occasion carries a latent personality factor measured by n_parcels
#' balanced item parcels ({p}_p1_{t} .. {p}_pK_{t}), with loadings and
#' intercepts constrained equal over occasions (metric and scalar longitudinal
#' invariance); the random-intercept/within decomposition is then placed on the
#' latent occasion factors, so within-person estimates are corrected for parcel
#' unreliability. the outcome stays single-indicator, identical to
#' build_riclpm_model(). same-parcel residual covariances across occasions are
#' included by default (correlated uniqueness), the standard specification for
#' repeated indicators.
#' identification: parcel 1 anchors each occasion factor (loading 1,
#' intercept 0); occasion factor means are free; exogenous latent means are 0.
#' @param personality personality trait stem (parcel columns {p}_pK_{t})
#' @param outcome outcome variable
#' @param waves number of waves
#' @param n_parcels parcels per occasion (default 3)
#' @param constrain_ar constrain autoregressive paths equal across time
#' @param constrain_cl constrain cross-lagged paths equal across time
#' @param correlated_uniqueness same-parcel residual covariances across occasions
build_riclpm_model_mi <- function(personality = "extr",
                                  outcome = "substance_use",
                                  waves = 5,
                                  n_parcels = 3,
                                  constrain_ar = TRUE,
                                  constrain_cl = TRUE,
                                  correlated_uniqueness = TRUE) {
  p <- personality
  o <- outcome

  parcel <- function(k, t) paste0(p, "_p", k, "_", t)
  eta <- function(t) paste0("eta_", p, t)
  o_vars <- paste0(o, "_", 1:waves)

  # --- MEASUREMENT: latent occasion factors with invariant loadings/intercepts ---
  meas <- sapply(1:waves, function(t) {
    rhs <- c(
      paste0("1*", parcel(1, t)),
      sapply(2:n_parcels, function(k) paste0("lp", k, "*", parcel(k, t)))
    )
    paste0(eta(t), " =~ ", paste(rhs, collapse = " + "))
  })

  intercepts <- unlist(lapply(1:waves, function(t) {
    c(
      paste0(parcel(1, t), " ~ 0*1"),
      sapply(2:n_parcels, function(k) paste0(parcel(k, t), " ~ ip", k, "*1"))
    )
  }))

  # occasion factor means carry the occasion-specific grand means
  eta_means <- sapply(1:waves, function(t) paste0(eta(t), " ~ 1"))

  # --- CORRELATED UNIQUENESS: same parcel across occasion pairs ---
  cu <- character(0)
  if (correlated_uniqueness && waves > 1) {
    pairs <- utils::combn(1:waves, 2)
    cu <- unlist(lapply(1:n_parcels, function(k) {
      apply(pairs, 2, function(st) {
        paste0(parcel(k, st[1]), " ~~ ", parcel(k, st[2]))
      })
    }))
  }

  # --- RANDOM INTERCEPTS (between-person) ---
  ri_p <- paste0(
    "RI_", p, " =~ ",
    paste0("1*", sapply(1:waves, eta), collapse = " + ")
  )
  ri_o <- paste0("RI_", o, " =~ ", paste0("1*", o_vars, collapse = " + "))

  # --- WITHIN-PERSON CENTERED VARIABLES ---
  wp_p <- sapply(1:waves, function(t) paste0("wp_", p, t, " =~ 1*", eta(t)))
  wp_o <- sapply(1:waves, function(t) paste0("wp_", o, t, " =~ 1*", o_vars[t]))

  # --- AUTOREGRESSIVE PATHS ---
  ar_lab_p <- if (constrain_ar) "a*" else ""
  ar_lab_o <- if (constrain_ar) "b*" else ""
  ar_p <- sapply(2:waves, function(t) {
    paste0("wp_", p, t, " ~ ", ar_lab_p, "wp_", p, t - 1)
  })
  ar_o <- sapply(2:waves, function(t) {
    paste0("wp_", o, t, " ~ ", ar_lab_o, "wp_", o, t - 1)
  })

  # --- CROSS-LAGGED PATHS ---
  cl_lab_c <- if (constrain_cl) "c*" else ""
  cl_lab_d <- if (constrain_cl) "d*" else ""
  cl_p_to_o <- sapply(2:waves, function(t) {
    paste0("wp_", o, t, " ~ ", cl_lab_c, "wp_", p, t - 1)
  })
  cl_o_to_p <- sapply(2:waves, function(t) {
    paste0("wp_", p, t, " ~ ", cl_lab_d, "wp_", o, t - 1)
  })

  # --- COVARIANCES ---
  ri_cov <- paste0("RI_", p, " ~~ RI_", o)
  wp_cov1 <- paste0("wp_", p, "1 ~~ wp_", o, "1")
  wp_cov <- sapply(2:waves, function(t) {
    paste0("wp_", p, t, " ~~ wp_", o, t)
  })

  # --- VARIANCE CONSTRAINTS ---
  # occasion factors and observed outcome decompose exactly into RI + wp
  eta_var0 <- sapply(1:waves, function(t) paste0(eta(t), " ~~ 0*", eta(t)))
  obs_var_o <- paste0(o_vars, " ~~ 0*", o_vars)

  paste(c(
    "# Measurement: latent occasion factors (metric + scalar invariance)",
    meas, intercepts, eta_means,
    "",
    "# Correlated uniqueness (same parcel across occasions)",
    cu,
    "",
    "# Random intercepts (between-person)",
    ri_p, ri_o,
    "",
    "# Within-person centered variables",
    wp_p, wp_o,
    "",
    "# Autoregressive paths",
    ar_p, ar_o,
    "",
    "# Cross-lagged paths",
    "# Personality -> Outcome (H1)",
    cl_p_to_o,
    "# Outcome -> Personality (H2)",
    cl_o_to_p,
    "",
    "# Covariances",
    ri_cov, wp_cov1, wp_cov,
    "",
    "# Fix decomposition residuals to 0",
    eta_var0, obs_var_o
  ), collapse = "\n")
}

#' fit the multiple-indicator RI-CLPM
#' @param data wide format data containing parcel columns {p}_pK_{t}
#' @param personality personality trait stem
#' @param outcome outcome variable
#' @param waves number of waves
#' @param n_parcels parcels per occasion
#' @param estimator "ML" or "MLR" (robust)
#' @param em_iter maximum EM iterations for missing data
#' @param ... additional arguments to lavaan
fit_riclpm_mi <- function(data, personality = "extr", outcome = "substance_use",
                          waves = 5, n_parcels = 3,
                          estimator = "MLR", em_iter = 10000, ...) {
  needed <- as.vector(sapply(1:n_parcels, function(k) {
    paste0(personality, "_p", k, "_", 1:waves)
  }))
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    warning(
      "fit_riclpm_mi: parcel columns missing (",
      paste(utils::head(missing_cols, 3), collapse = ", "),
      if (length(missing_cols) > 3) ", ..." else "",
      "); run prepare_analysis_data(include_parcels = TRUE)"
    )
    return(list(
      personality = personality, outcome = outcome, waves = waves,
      model = NULL, fit = NULL
    ))
  }

  model <- build_riclpm_model_mi(personality, outcome, waves, n_parcels)

  fit <- tryCatch(
    {
      lavaan::sem(model,
        data = data,
        missing = "fiml",
        estimator = estimator,
        em.h1.iter.max = em_iter,
        ...
      )
    },
    error = function(e) {
      warning("multiple-indicator RI-CLPM fitting failed: ", e$message)
      NULL
    }
  )

  list(
    personality = personality,
    outcome = outcome,
    waves = waves,
    model = model,
    fit = fit
  )
}
