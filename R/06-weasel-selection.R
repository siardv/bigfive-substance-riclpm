# 06-weasel-selection.R
# rule-based respondent/wave selection with weasel, inserted between
# prepare_analysis_data() and reshape_for_riclpm().

# =============================================================================
# WEASEL SUBSET SELECTION
# =============================================================================
# weasel operates on long-format panel data keyed by an id and an integer wave
# column, the contract of analysis_data (nomem_encr, wavenr). it scores
# competing inclusion scenarios, applies one, and emits a methods paragraph
# plus audit tables; the chosen window length L then drives the wide reshape.
# all calls are namespace-qualified, so the package needs installing only.

.load_weasel <- function() {
  if (!requireNamespace("weasel", quietly = TRUE)) {
    stop(
      "the weasel package is required but not installed. install it with:\n",
      "  install.packages(\"weasel\")   # or remotes::install_github(\"siardv/weasel\")",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' translate the selected waves to calendar years via the weasel wave_map and
#' report whether the analysis waves are annual or biennial. wave_map carries the
#' dense wave_idx and the original wavenr; calendar year is wavenr_orig + 2007.
describe_selected_waves <- function(wave_map) {
  if (is.null(wave_map)) {
    cat("wave_map is NULL: weasel did not preserve the original-wave column, so ",
        "selected-wave calendar years cannot be recovered from the selection ",
        "step alone (only module-level coverage from describe_wave_calendar()).\n",
        sep = "")
    return(invisible(NULL))
  }
  wm <- wave_map
  wm$year <- wm$wavenr_orig + 2007L

  # (a) which calendar years sit behind each dense wave index
  per_idx <- wm %>%
    dplyr::group_by(wave_idx) %>%
    dplyr::summarize(
      years    = paste(sort(unique(year)), collapse = "/"),
      n_years  = dplyr::n_distinct(year),
      n_people = dplyr::n_distinct(nomem_encr),
      .groups  = "drop"
    ) %>%
    dplyr::arrange(wave_idx)
  cat("\ndense wave_idx -> calendar year(s):\n")
  print(as.data.frame(per_idx), row.names = FALSE)

  # (b) within-person spacing between consecutive selected waves (calendar years)
  gap_tab <- wm %>%
    dplyr::arrange(nomem_encr, wave_idx) %>%
    dplyr::group_by(nomem_encr) %>%
    dplyr::summarize(g = list(diff(sort(unique(year)))), .groups = "drop") %>%
    tidyr::unnest(g) %>%
    dplyr::count(g, name = "n_intervals") %>%
    dplyr::arrange(g)
  cat("\nwithin-person gaps between consecutive selected waves (years):\n")
  print(as.data.frame(gap_tab), row.names = FALSE)

  # (c) verdict
  staggered <- any(per_idx$n_years > 1)
  annual    <- nrow(gap_tab) == 1 && gap_tab$g[1] == 1
  cat("\n  grid: ", if (staggered)
        "staggered (a wave index spans multiple calendar years across respondents)"
      else "fixed (each wave index is one calendar year for everyone)", "\n", sep = "")
  cat("  cadence: ", if (annual) "annual"
      else paste0("not annual (modal gap = ",
                  gap_tab$g[which.max(gap_tab$n_intervals)], " years)"),
      "\n", sep = "")
  invisible(list(per_idx = per_idx, gaps = gap_tab))
}

#' select a principled longitudinal subset with weasel
#'
#' @param analysis_data long-format output of prepare_analysis_data()
#' @param id id column name (LISS: "nomem_encr")
#' @param wave integer wave column name (LISS: "wavenr")
#' @param measured character vector of variable(s) that must be observed for a
#'   person-wave row to count as present; defaults to "extr". this aligns
#'   weasel's wave-presence logic with actual measurement, since LISS modules
#'   are fielded on different schedules. set to character(0) to disable.
#' @param span "core" (highest-coverage window of core_len) or "full"
#' @param core_len desired window length when span = "core"
#' @param scenario scenario to apply; NULL uses weasel's recommended one
#' @param author,year citation fields woven into the methods paragraph
#' @param print whether to print comparison, audit, and justification output
#' @return list: data (selected long subset), L (chosen window length),
#'   scenario, lower, upper, wave_map (dense index to original wavenr mapping),
#'   plan, comparison, summary, justification
select_weasel_subset <- function(analysis_data,
                                  id = "nomem_encr",
                                  wave = "wavenr",
                                  measured = "extr",
                                  span = "core",
                                  core_len = 5L,
                                  scenario = NULL,
                                  author = "van den Bosch",
                                  year = "2026",
                                  print = TRUE) {
  .load_weasel()

  # weasel needs the bare id/wave panel; drop rows missing either key
  keys_ok <- !is.na(analysis_data[[id]]) & !is.na(analysis_data[[wave]])
  panel <- analysis_data[keys_ok, , drop = FALSE]

  # modules are fielded on different schedules, so a row can exist at a wavenr
  # while the modeled variable was never measured there. keep only rows where
  # `measured` is observed and re-index each person's waves to a dense order
  # (1, 2, ...); weasel then plans on that axis, where "present at waves 1:L"
  # means L consecutive measurements, the RI-CLPM unit that
  # reshape_for_riclpm() expects. set measured = character(0) to opt out.
  measured <- intersect(measured, names(panel))
  if (length(measured) > 0) {
    obs <- rowSums(!is.na(panel[, measured, drop = FALSE])) == length(measured)
    n_before <- nrow(panel)
    panel <- panel[obs, , drop = FALSE]
    if (isTRUE(print)) {
      cat("weasel: restricted to rows with ",
          paste(measured, collapse = " + "), " observed (",
          nrow(panel), " of ", n_before, " person-wave rows)\n", sep = "")
    }
  }
  if (nrow(panel) == 0) {
    stop("no rows remain after requiring measured variable(s); ",
         "check the `measured` argument.", call. = FALSE)
  }

  # dense per-person measurement-order index; keep the original wave for mapping
  orig_wave <- paste0(".", wave, "_orig")
  panel[[orig_wave]] <- panel[[wave]]
  ord <- order(panel[[id]], panel[[wave]])
  panel <- panel[ord, , drop = FALSE]
  panel[[wave]] <- stats::ave(panel[[wave]], panel[[id]],
                              FUN = function(w) seq_along(w))
  panel[[wave]] <- as.integer(panel[[wave]])

  # 1) build and score scenarios on the dense axis
  plan_obj <- weasel::weasel_plan(panel, id = id, wave = wave,
                                  span = span, core_len = as.integer(core_len))
  comparison <- weasel::weasel_compare_scenarios(plan_obj)

  # guard: a degenerate panel can leave every scenario empty (no recommendation)
  if (!any(comparison$recommended) || all(comparison$n_ids == 0)) {
    stop("no scenario retained any respondents on the dense measurement axis; ",
         "try a smaller core_len or check measurement coverage.", call. = FALSE)
  }

  # 2) resolve the scenario (recommended unless overridden)
  if (is.null(scenario)) {
    scenario <- comparison$scenario[comparison$recommended][1]
  } else {
    scenario <- weasel::weasel_match_scenario(scenario, comparison$scenario)
  }
  row <- comparison[comparison$scenario == scenario, , drop = FALSE]
  L <- as.integer(row$L[[1]])
  lower <- as.integer(row$lower[[1]])
  upper <- as.integer(row$upper[[1]])

  # 3) audit the chosen scenario and generate the methods paragraph
  summ <- weasel::weasel_summarize_subset(plan_obj, scenario,
                                          data = panel, id = id, wave = wave)
  selectivity <- tryCatch(
    weasel::weasel_selectivity(plan_obj, scenario),
    error = function(e) NULL
  )
  sensitivity <- tryCatch(
    weasel::weasel_sensitivity(plan_obj),
    error = function(e) NULL
  )
  justification <- weasel::weasel_justify_subset(
    plan_obj, scenario,
    style = "methods", author = author, year = year
  )
  # normalise en/em dashes to hyphen for paper-bound text
  justification <- gsub("[\u2012\u2013\u2014\u2015]", "-", justification)

  # 4) extract the selected long subset (carried on the dense axis so it feeds
  #    reshape_for_riclpm() directly). keep the wave mapping separately so the
  #    returned data is clean to pivot; an extra column would become a spurious
  #    wide variable in reshape_for_riclpm(). read wave_map to confirm the dense
  #    window maps onto a clean annual run of original waves.
  subset_long <- weasel::weasel_apply(plan_obj, scenario)
  wave_map <- NULL
  if (orig_wave %in% names(subset_long)) {
    wave_map <- unique(subset_long[, c(id, wave, orig_wave), drop = FALSE])
    names(wave_map) <- c(id, "wave_idx", "wavenr_orig")
    subset_long[[orig_wave]] <- NULL
  }

  if (isTRUE(print)) {
    cat("\n", strrep("=", 80), "\n", sep = "")
    cat("WEASEL SUBSET SELECTION\n")
    cat(strrep("=", 80), "\n\n", sep = "")

    weasel::weasel_print_table(
      comparison[, c("scenario", "n_ids", "mean_prop_present",
                     "endpoint_rate", "max_missing", "n_gap_max",
                     "max_gap_max", "score", "recommended")],
      title = "Scenario comparison (scored)", digits = 3
    )
    cat(weasel::weasel_compare_to_sentence(comparison), "\n")

    cat("\nChosen scenario: ", scenario,
        " (waves ", lower, ":", upper, ", L = ", L, ")\n", sep = "")
    weasel::weasel_print_table(summ$headline, title = "Chosen subset headline")
    weasel::weasel_print_table(summ$per_wave_coverage, title = "Per-wave coverage")
    weasel::weasel_print_table(summ$missing_distribution,
                               title = "Missingness distribution")
    if (!is.null(selectivity)) {
      weasel::weasel_print_table(
        selectivity,
        title = "Selectivity: retained vs excluded (standardized mean differences)",
        digits = 3
      )
    }
    if (!is.null(sensitivity)) {
      weasel::weasel_print_table(
        utils::head(sensitivity[order(-sensitivity$n_ids), ], 8),
        title = "Sensitivity sweep: selection tolerances vs retained sample",
        digits = 3
      )
    }

    # translate the dense selection back to calendar years and report cadence
    describe_selected_waves(wave_map)

    cat("\nMethods-section justification (paste-ready):\n")
    cat(strrep("-", 80), "\n")
    cat(justification, "\n")
    cat(strrep("=", 80), "\n\n")
  }

  list(
    data          = subset_long,
    L             = L,
    scenario      = scenario,
    lower         = lower,
    upper         = upper,
    wave_map      = wave_map,
    plan          = plan_obj,
    comparison    = comparison,
    summary       = summ,
    selectivity   = selectivity,
    sensitivity   = sensitivity,
    justification = justification
  )
}

# =============================================================================
# UNIFORM-INTERVAL SUBSAMPLE
# =============================================================================

#' identify respondents whose selected window is spaced at exactly one year
#' the calendar-staggered grid mixes one-year with multi-year transitions
#' inside the same lag-1 coefficient; this helper reads the weasel wave_map
#' and returns the ids whose consecutive selected waves are all one calendar
#' year apart, for the uniform-interval sensitivity refit.
#' @param weasel_selection return value of select_weasel_subset()
filter_uniform_interval <- function(weasel_selection) {
  wm <- weasel_selection$wave_map
  if (is.null(wm) || !all(c("nomem_encr", "wavenr_orig") %in% names(wm))) {
    message("filter_uniform_interval: wave_map unavailable; no filter applied")
    return(NULL)
  }

  gaps_by_person <- tapply(wm$wavenr_orig, wm$nomem_encr, function(w) {
    diff(sort(unique(w)))
  })
  all_annual <- vapply(gaps_by_person, function(g) {
    length(g) > 0 && all(g == 1L)
  }, logical(1))
  ids <- as.numeric(names(all_annual)[all_annual])

  gap_pool <- unlist(gaps_by_person, use.names = FALSE)
  message(
    "  uniform-interval filter: ", length(ids), " of ",
    length(all_annual), " respondents have all one-year transitions (",
    round(100 * length(ids) / length(all_annual), 1), "%)"
  )
  message(
    "  transition-gap distribution before filtering: ",
    paste(names(table(gap_pool)), "y=", as.vector(table(gap_pool)),
      collapse = ", ", sep = ""
    )
  )

  list(
    ids = ids,
    n_uniform = length(ids),
    n_total = length(all_annual),
    gap_table = table(gap_pool)
  )
}
