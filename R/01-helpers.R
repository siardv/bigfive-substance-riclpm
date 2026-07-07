# 01-helpers.R
# data loading (load_liss_data) and general-purpose utilities.

`%>%`  <- magrittr::`%>%`
`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# DATA LOADING
# =============================================================================

#' load LISS data from the frozen per-module merges in config$merged_dir, adapt to
#' the project schema, and read background variables from disk. returns a
#' per-module list ($p, $h, $i, $b, $s, $bck) consumed by clean_income() and
#' prepare_analysis_data(). refresh the merged files with merge_liss_to_disk();
#' $s and $bck are unused downstream and kept NULL.
#' @param config list with merged_dir, avars_file, liss_modules
load_liss_data <- function(config) {
  modules <- config$liss_modules %||% c("cp", "ch", "ci")
  if (!all(c("cp", "ch", "ci") %in% modules)) {
    stop("liss_modules must include 'cp', 'ch', and 'ci'", call. = FALSE)
  }
  merged <- lapply(modules, function(m) {
    path <- file.path(config$merged_dir, paste0(m, "_merged.sav"))
    if (!file.exists(path)) {
      stop("merged file not found: ", path,
           "\nrebuild it with merge_liss_to_disk(config).", call. = FALSE)
    }
    adapt_lissr_module(haven::read_sav(path), m)
  })
  names(merged) <- modules

  list(
    p   = merged[["cp"]],
    h   = merged[["ch"]],
    i   = merged[["ci"]],
    b   = haven::read_sav(config$avars_file),
    s   = NULL,
    bck = NULL
  )
}

#' map a LISS wave_year to the project wavenr. module-independent so cp and ch of
#' one panel year share a wavenr, which pairs extr with substance_use in
#' prepare_analysis_data(). the 2007 offset yields the biennial personality
#' structure that select_weasel_subset()'s dense re-indexing expects.
lissr_year_to_wavenr <- function(year) as.integer(year - 2007L)

#' adapt one merged module to the project schema: rename harmonized items
#' s<NNN> to <module>00a<NNN> (the names 04-measures.R scores on) and add wavenr.
#' @param d data frame read from a *_merged.sav file
#' @param module two-letter module code ("cp", "ch", "ci")
adapt_lissr_module <- function(d, module,
                               year_to_wavenr = lissr_year_to_wavenr) {
  nm <- names(d)
  is_item <- grepl("^s[0-9]+$", nm)
  nm[is_item] <- paste0(module, "00a", sub("^s", "", nm[is_item]))
  names(d) <- nm
  if (!"wave_year" %in% names(d)) {
    stop("merged file for module '", module, "' lacks wave_year", call. = FALSE)
  }
  d$wavenr <- year_to_wavenr(d$wave_year)
  d
}

# =============================================================================
# DATA BUILD (run only to refresh the frozen merges; not on the analysis path)
# =============================================================================

#' fetch raw LISS per-wave files into config$raw_dir (run once; interactive login
#' + 2FA). store credentials first with lissr::liss_store_credentials(username).
#' @param config list with raw_dir and liss_modules
#' @param waves optional wave codes passed to lissr::liss_download(.waves = )
fetch_liss_raw <- function(config, waves = NULL) {
  if (!dir.exists(config$raw_dir)) dir.create(config$raw_dir, recursive = TRUE)
  lissr::liss_login()
  lissr::liss_download(
    .dir     = config$raw_dir,
    .modules = config$liss_modules %||% c("cp", "ch", "ci"),
    .waves   = waves
  )
  invisible(config$raw_dir)
}

#' rebuild the frozen per-module merges in config$merged_dir from the raw per-wave
#' files in config$raw_dir via lissr. run after fetch_liss_raw() when source data
#' changes; load_liss_data() then reads the refreshed files.
#' @param config list with raw_dir, merged_dir, liss_modules
merge_liss_to_disk <- function(config) {
  modules <- config$liss_modules %||% c("cp", "ch", "ci")
  if (!dir.exists(config$merged_dir)) dir.create(config$merged_dir, recursive = TRUE)
  for (m in modules) {
    lissr::merge_liss_module(
      lissr::liss_recipe(m),
      data_dir   = config$raw_dir,
      output_dir = config$merged_dir
    )
  }
  invisible(config$merged_dir)
}

# =============================================================================
# HELPER FUNCTIONS: GENERAL
# =============================================================================

#' most frequent value
stat_mode <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' drop NA and infinite values
rna <- function(x, flat = FALSE) {
  if (flat) x <- unlist(x)
  if (is.list(x)) stats::na.omit(x) else x[!is.na(x) & !is.infinite(x)]
}

#' move values at `ind` from `col` into `dummy`, then set `col` to NA
dummy_na <- function(df, ind, col, dummy, copy_value = TRUE) {
  if (length(ind) == 0) return(df)
  df[[dummy]][ind] <- if (copy_value) df[[col]][ind] else 1
  df[[col]][ind] <- NA
  df
}

#' convert wave (YYYYMM) to wavenr; reference 200711 = 1
wave_to_wavenr <- function(wave) {
  year  <- as.integer(wave %/% 100)
  month <- as.integer(wave %% 100)
  (year - 2007) * 12 + (month - 11) + 1L
}

# =============================================================================
# WAVE CALENDAR PROVENANCE (diagnostic)
# =============================================================================

#' print the calendar provenance of each module's waves. wavenr is defined as
#' wave_year - 2007 (see lissr_year_to_wavenr); month is not carried on the
#' analysis path (only documented: cp ~ May/June, ch ~ November, ci monthly).
#' reads wave_year, which adapt_lissr_module() leaves in place.
#' @param liss per-module list from load_liss_data() ($p, $h, $i)
#' @param modules named vector mapping list slots to printable labels
describe_wave_calendar <- function(liss,
                                   modules = c(p = "personality (cp)",
                                               h = "health (ch)",
                                               i = "income (ci)")) {
  for (key in names(modules)) {
    d <- liss[[key]]
    if (is.null(d) || !all(c("wave_year", "wavenr") %in% names(d))) next
    tab <- d %>%
      dplyr::distinct(wave_year, wavenr) %>%
      dplyr::arrange(wavenr)
    cat("\n", modules[[key]], ", wave_year -> wavenr (= wave_year - 2007):\n",
        sep = "")
    print(as.data.frame(tab), row.names = FALSE)
    yrs  <- sort(unique(tab$wavenr)) + 2007L
    gaps <- diff(yrs)
    cat("  calendar years present: ", paste(yrs, collapse = ", "), "\n", sep = "")
    cat("  year-to-year spacing: ", paste(gaps, collapse = ", "),
        if (length(gaps) > 0 && all(gaps == 1)) "  (annual)"
        else "  (not a clean annual run)", "\n", sep = "")
  }
  invisible(NULL)
}
