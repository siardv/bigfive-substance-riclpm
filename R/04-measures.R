# 04-measures.R
# measurement: Big Five items/scoring and substance use.

# =============================================================================
# BIG FIVE DEFINITIONS
# =============================================================================

# IPIP 50-item Big Five items
b5_items <- list(
  extr = c("020", "025", "030", "035", "040", "045", "050", "055", "060", "065"),
  agre = c("021", "026", "031", "036", "041", "046", "051", "056", "061", "066"),
  cons = c("022", "027", "032", "037", "042", "047", "052", "057", "062", "067"),
  neur = c("023", "028", "033", "038", "043", "048", "053", "058", "063", "068"),
  open = c("024", "029", "034", "039", "044", "049", "054", "059", "064", "069")
)

# reverse-scored item positions (1-based)
b5_reverse <- list(
  extr = c(2, 4, 6, 8, 10),
  agre = c(1, 3, 5, 7),
  cons = c(2, 4, 6, 8),
  neur = c(2, 4),
  open = c(2, 4, 6)
)

# substance use items (ch module)
drugs_use_items <- c("159", "160", "161", "162", "163", "270")

# sedatives-excluded variant for the pre-registered-style sensitivity run:
# drops item 159 (sedatives/tranquillizers, which can reflect prescribed use)
illicit_five_items <- setdiff(drugs_use_items, "159")

# =============================================================================
# BIG FIVE EXTRACTION AND SCORING
# =============================================================================

#' construct LISS column names
liss_cols <- function(module, items) {
  paste0(module, "00a", items)
}

#' extract columns from data frame
liss_select <- function(data, module, items = NULL, ids = TRUE) {
  col_names <- names(data)
  id_cols <- intersect(c("nomem_encr", "nohouse_encr", "wave", "wavenr"), col_names)

  if (is.null(items)) {
    mod_cols <- grep(paste0("^", module, "00a"), col_names, value = TRUE)
  } else {
    mod_cols <- intersect(liss_cols(module, items), col_names)
  }

  if (length(mod_cols) == 0) {
    return(tibble::tibble())
  }

  cols <- if (ids) c(id_cols, mod_cols) else mod_cols
  out <- dplyr::select(data, dplyr::all_of(cols))
  out[rowSums(!is.na(out[mod_cols])) > 0, ]
}

#' reverse score
reverse_score <- function(x, max_val) {
  max_val + 1 - x
}

#' compute scale score
compute_scale <- function(data, items, reverse_idx = NULL, max_val = 5,
                          method = c("mean", "sum"), min_valid = NULL) {
  method <- match.arg(method)
  item_data <- as.data.frame(data[items])

  if (!is.null(reverse_idx)) {
    for (i in reverse_idx) {
      item_data[[i]] <- reverse_score(item_data[[i]], max_val)
    }
  }

  n_valid <- rowSums(!is.na(item_data))
  min_req <- min_valid %||% length(items)

  score <- if (method == "mean") {
    rowMeans(item_data, na.rm = TRUE)
  } else {
    rowSums(item_data, na.rm = TRUE)
  }

  score[n_valid < min_req] <- NA
  score
}

#' create Big Five parcels for improved CFA fit
#' parceling reduces model complexity and improves fit when items are numerous
#' @param data data with B5 items
#' @param dimension personality dimension
#' @param n_parcels number of parcels (2 or 3)
create_b5_parcels <- function(data, dimension, n_parcels = 2) {
  items <- liss_cols("cp", b5_items[[dimension]])
  items_present <- intersect(items, names(data))

  if (length(items_present) < n_parcels) {
    warning("not enough items for parceling")
    return(data)
  }

  # assign items to parcels (balanced)
  parcel_assign <- rep(1:n_parcels, length.out = length(items_present))

  for (p in 1:n_parcels) {
    parcel_items <- items_present[parcel_assign == p]
    parcel_name <- paste0(dimension, "_p", p)

    # need to reverse score first
    orig_positions <- match(parcel_items, items)
    reverse_idx <- which(orig_positions %in% b5_reverse[[dimension]])

    item_data <- as.data.frame(data[parcel_items])
    if (length(reverse_idx) > 0) {
      for (i in reverse_idx) {
        item_data[[i]] <- 5 + 1 - item_data[[i]]
      }
    }

    data[[parcel_name]] <- rowMeans(item_data, na.rm = TRUE)
    data[[parcel_name]][rowSums(!is.na(item_data)) == 0] <- NA
  }

  data
}

#' compute Big Five scores with optional parceling
#' @param p personality data
#' @param max_val maximum item value
#' @param min_valid minimum valid items for score
#' @param use_parcels create parcel variables for CFA
compute_b5_scores <- function(p, max_val = 5, min_valid = 8, use_parcels = FALSE) {
  scores <- tibble::tibble(
    nomem_encr = p$nomem_encr,
    wave = if ("wave" %in% names(p)) p$wave else NA,
    wavenr = if ("wavenr" %in% names(p)) p$wavenr else NA
  )

  for (dim in names(b5_items)) {
    items <- liss_cols("cp", b5_items[[dim]])
    items_present <- intersect(items, names(p))
    orig_positions <- match(items_present, items)
    reverse_adj <- which(orig_positions %in% b5_reverse[[dim]])

    scores[[dim]] <- compute_scale(p, items_present, reverse_adj,
      max_val,
      min_valid = min_valid
    )
  }

  # optionally add parcels for CFA
  if (use_parcels) {
    for (dim in names(b5_items)) {
      scores <- create_b5_parcels(
        cbind(scores, p[intersect(liss_cols("cp", b5_items[[dim]]), names(p))]),
        dim,
        n_parcels = 2
      )
    }
  }

  scores
}

# =============================================================================
# SUBSTANCE USE EXTRACTION
# =============================================================================

#' drug use variable names (ch module)
#' ch00a159 = sedatives, ch00a160 = soft drugs (cannabis), ch00a161 = XTC,
#' ch00a162 = hallucinogens, ch00a163 = hard drugs, ch00a270 = laughing gas (added 2020, mostly missing pre-2020)
drugs_use_vars <- c("ch00a159", "ch00a160", "ch00a161", "ch00a162", "ch00a163", "ch00a270")

#' sedatives-excluded variable set (companion to illicit_five_items)
illicit_five_vars <- setdiff(drugs_use_vars, "ch00a159")

#' extract substance use data from health module
#' @param h health data (ch module)
#' @param items substance item variables to extract (default: full battery)
extract_substance_use <- function(h, items = drugs_use_vars) {
  # identify available columns
  id_cols <- c("nomem_encr", "wave", "wavenr")
  id_cols <- intersect(id_cols, names(h))

  use_cols <- intersect(items, names(h))

  if (length(use_cols) == 0) {
    warning("no substance use variables found")
    return(NULL)
  }

  # select and clean
  out <- dplyr::select(h, dplyr::all_of(c(id_cols, use_cols)))

  # convert labelled to numeric
  for (col in use_cols) {
    if (col %in% names(out)) {
      out[[col]] <- sjlabelled::as_numeric(out[[col]], use.labels = FALSE)
    }
  }

  # filter to rows with at least one non-NA substance use variable
  out[rowSums(!is.na(out[use_cols])) > 0, ]
}

#' recode substance use to binary (0 = no use, 1 = any use)
#' LISS scale is three-point (1 = never, 2 = sometimes, 3 = regularly);
#' collapse to any past-year use: 1 -> 0, 2 or higher -> 1
#' @param use_data substance use data from extract_substance_use()
#' @param items substance item variables to recode (default: full battery)
recode_substance_binary <- function(use_data, items = drugs_use_vars) {
  use_cols <- intersect(items, names(use_data))

  for (col in use_cols) {
    use_data[[col]] <- dplyr::case_when(
      use_data[[col]] == 1 ~ 0L,
      use_data[[col]] == 2 ~ 1L,
      use_data[[col]] >= 3 ~ 1L, # any other positive response = use
      TRUE ~ NA_integer_
    )
  }
  use_data
}

#' compute substance use composite score
#' options: sum (count of substances used), any (binary any use)
#' @param use_data substance use data (binary coded)
#' @param method "sum" for count, "any" for binary, "mean" for proportion
#' @param items substance item variables to score over (default: full battery)
#' @param min_valid minimum answered items for a non-NA score (default 1,
#'   the historical behavior); note that the sixth item (ch00a270) exists only
#'   from 2020, so most person-waves can answer at most five items
compute_substance_score <- function(use_data, method = c("sum", "any", "mean"),
                                    items = drugs_use_vars, min_valid = 1) {
  method <- match.arg(method)
  use_cols <- intersect(items, names(use_data))

  use_matrix <- as.matrix(use_data[use_cols])

  score <- switch(method,
    "sum" = rowSums(use_matrix, na.rm = TRUE),
    "any" = as.integer(rowSums(use_matrix, na.rm = TRUE) > 0),
    "mean" = rowMeans(use_matrix, na.rm = TRUE)
  )

  # set to NA below the answered-items threshold
  n_valid <- rowSums(!is.na(use_matrix))
  score[n_valid < min_valid] <- NA

  use_data$substance_use <- score
  use_data
}

#' full substance use processing pipeline
#' @param h health data
#' @param method scoring method ("sum", "any", "mean")
#' @param items substance item variables (default: full six-item battery;
#'   pass illicit_five_vars for the sedatives-excluded sensitivity)
#' @param min_valid minimum answered items for a non-NA score
process_substance_use <- function(h, method = "sum",
                                  items = drugs_use_vars, min_valid = 1) {
  h %>%
    extract_substance_use(items = items) %>%
    recode_substance_binary(items = items) %>%
    compute_substance_score(method = method, items = items, min_valid = min_valid)
}

#' describe item availability of the substance battery by calendar year
#' reports, per item, the person-wave counts with an observed response by
#' wave year, plus the distribution of answered-items-per-row. this backs the
#' manuscript footnote that ch00a270 (laughing gas) enters the battery in 2020,
#' so the effective range is 0-5 for earlier person-waves.
#' @param h health data (ch module) with a wavenr or wave column
#' @param items substance item variables to describe
describe_substance_items <- function(h, items = drugs_use_vars) {
  use_cols <- intersect(items, names(h))
  if (length(use_cols) == 0) {
    message("describe_substance_items: no substance items found")
    return(invisible(NULL))
  }

  d <- h
  if (!"wavenr" %in% names(d) && "wave" %in% names(d)) {
    d$wavenr <- wave_to_wavenr(d$wave)
  }
  if (!"wavenr" %in% names(d)) {
    message("describe_substance_items: no wave information found")
    return(invisible(NULL))
  }
  d$year <- d$wavenr + 2007L

  for (col in use_cols) {
    d[[col]] <- sjlabelled::as_numeric(d[[col]], use.labels = FALSE)
  }

  # restrict to rows that touched the battery at all
  answered <- rowSums(!is.na(d[use_cols]))
  d <- d[answered > 0, , drop = FALSE]
  answered <- answered[answered > 0]

  message("substance battery availability (person-waves with an observed response):")
  avail <- sapply(use_cols, function(col) {
    tapply(!is.na(d[[col]]), d$year, sum)
  })
  print(avail)

  message("\nanswered items per person-wave:")
  print(table(answered))

  if ("ch00a270" %in% use_cols) {
    n270 <- sum(!is.na(d$ch00a270))
    message(
      "\nch00a270 (laughing gas, added 2020): observed on ", n270, " of ",
      nrow(d), " battery rows (",
      round(100 * n270 / nrow(d), 1), "%)"
    )
  }

  invisible(list(availability_by_year = avail, answered_per_row = table(answered)))
}
