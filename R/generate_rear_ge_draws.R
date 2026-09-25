#' Generate rear-specific daily GE draws from a joint rear-type fit
#'
#' Extracts one rear group's weekly guidance-efficiency posterior from
#' [fit_ge_rear_model2()] and returns coherent posterior columns for SCRAPI2.
#' When daily spill is supplied, each draw retains its fitted weekly residual
#' while replacing weekly spill with observed daily spill.
#'
#' @param ge_fit Result from [fit_ge_rear_model2()].
#' @param rear_type Rear group to extract. Defaults to the fit's target rear.
#' @param strat_assign Date-to-model-row mapping. Defaults to
#'   `ge_fit$strat_assign`. Accepts `Week`/`Collapse` or
#'   `date`/`stratum_idx` columns.
#' @param pass_dates Passage-data dates in output row order.
#' @param B Number of posterior columns returned.
#' @param daily_spill Optional data frame containing `Date` and `spill.per`.
#'   When supplied, daily spill is applied using the fitted spill slope.
#' @param clip_to_ci Restrict sampling to complete MCMC rows where all weekly
#'   psi values for `rear_type` fall inside their marginal intervals. Default
#'   `FALSE`; full-posterior sampling is recommended for uncertainty propagation.
#' @param ci_lower,ci_upper Marginal bounds used only when `clip_to_ci = TRUE`.
#' @param strict_dates Stop when a passage date has no model-row assignment.
#'   Default `TRUE`.
#' @param strict_daily_spill When daily spill is supplied, stop if any modeled
#'   passage date lacks a daily value. If `FALSE`, use that week's fitted mean.
#' @param seed Optional sampling seed.
#'
#' @return Data frame with `SampleEndDate` followed by `boot_1` through
#'   `boot_B`. Every column is one coherent posterior realization.
#'
#' @examples
#' \dontrun{
#' ge_W <- generate_rear_ge_draws(
#'   rear_fit, rear_type = "W",
#'   pass_dates = passageData$SampleEndDate,
#'   B = 5000, daily_spill = lgr_daily_spill,
#'   clip_to_ci = FALSE, seed = 11
#' )
#' }
#' @importFrom stats plogis qlogis quantile
#' @export
generate_rear_ge_draws <- function(
    ge_fit,
    rear_type = NULL,
    strat_assign = ge_fit$strat_assign,
    pass_dates,
    B = 2000,
    daily_spill = NULL,
    clip_to_ci = FALSE,
    ci_lower = 0.025,
    ci_upper = 0.975,
    strict_dates = TRUE,
    strict_daily_spill = TRUE,
    seed = NULL) {

  parse_dates <- function(x, label) {
    if (inherits(x, "Date")) return(x)
    raw <- trimws(as.character(x))
    out <- as.Date(rep(NA_character_, length(raw)))
    for (fmt in c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y")) {
      todo <- which(is.na(out) & !is.na(raw) & nzchar(raw))
      if (!length(todo)) break
      out[todo] <- as.Date(raw[todo], format = fmt)
    }
    if (anyNA(out)) stop("Could not parse every ", label, ".", call. = FALSE)
    out
  }
  whole <- function(x, name, minimum = 1L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < minimum || x != floor(x)) {
      stop(name, " must be a finite integer >= ", minimum, ".", call. = FALSE)
    }
    as.integer(x)
  }

  B <- whole(B, "B")
  for (nm in c("clip_to_ci", "strict_dates", "strict_daily_spill")) {
    value <- get(nm)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(nm, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (!is.numeric(ci_lower) || !is.numeric(ci_upper) ||
      length(ci_lower) != 1L || length(ci_upper) != 1L ||
      !is.finite(ci_lower) || !is.finite(ci_upper) ||
      ci_lower < 0 || ci_upper > 1 || ci_lower >= ci_upper) {
    stop("ci_lower and ci_upper must satisfy 0 <= lower < upper <= 1.",
         call. = FALSE)
  }
  needed <- c("samples", "rear_levels", "weeks", "spill_mean", "spill_sd",
              "lgr_spill_std")
  missing <- setdiff(needed, names(ge_fit))
  if (length(missing)) {
    stop("ge_fit is missing: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }
  if (is.null(rear_type)) rear_type <- ge_fit$target_rear
  r <- match(rear_type, ge_fit$rear_levels)
  if (length(r) != 1L || is.na(r)) {
    stop("rear_type must be one of: ",
         paste(ge_fit$rear_levels, collapse = ", "), ".", call. = FALSE)
  }

  pass_dates_d <- parse_dates(pass_dates, "passage date")
  if (!is.data.frame(strat_assign)) {
    stop("strat_assign must be a data frame.", call. = FALSE)
  }
  if (all(c("Week", "Collapse") %in% names(strat_assign)) &&
      !"date" %in% names(strat_assign)) {
    map <- strat_assign[, c("Week", "Collapse"), drop = FALSE]
    if (anyDuplicated(map$Week)) stop("strat_assign has duplicate Week rows.",
                                     call. = FALSE)
    s <- map$Collapse[match(as.integer(format(pass_dates_d, "%V")), map$Week)]
  } else if (all(c("date", "stratum_idx") %in% names(strat_assign))) {
    map_dates <- parse_dates(strat_assign$date, "strat_assign date")
    if (anyDuplicated(map_dates)) stop("strat_assign has duplicate dates.",
                                      call. = FALSE)
    s <- strat_assign$stratum_idx[match(pass_dates_d, map_dates)]
  } else {
    stop("strat_assign needs Week/Collapse or date/stratum_idx columns.",
         call. = FALSE)
  }
  valid_s <- !is.na(s) & s == floor(s) & s >= 1L & s <= length(ge_fit$weeks)
  if (strict_dates && any(!valid_s)) {
    stop(sum(!valid_s), " passage date(s) have no fitted model row.",
         call. = FALSE)
  }

  mat <- do.call(rbind, lapply(ge_fit$samples, as.matrix))
  psi_cols <- paste0("psi[", r, ",", seq_along(ge_fit$weeks), "]")
  if (!all(c(psi_cols, "beta") %in% colnames(mat))) {
    stop("The requested rear-specific psi or beta draws are absent.",
         call. = FALSE)
  }
  psi_all <- mat[, psi_cols, drop = FALSE]
  valid_rows <- seq_len(nrow(mat))
  if (clip_to_ci) {
    lo <- apply(psi_all, 2L, stats::quantile, probs = ci_lower)
    hi <- apply(psi_all, 2L, stats::quantile, probs = ci_upper)
    inside <- sweep(psi_all, 2L, lo, `>=`) & sweep(psi_all, 2L, hi, `<=`)
    valid_rows <- which(rowSums(inside) == ncol(inside))
    if (!length(valid_rows)) {
      stop("No complete posterior rows satisfy the requested CI restriction.",
           call. = FALSE)
    }
  }
  if (!is.null(seed)) set.seed(seed)
  draw_id <- sample(valid_rows, B, replace = length(valid_rows) < B)
  selected <- mat[draw_id, , drop = FALSE]
  psi <- pmin(pmax(selected[, psi_cols, drop = FALSE], 1e-9), 1 - 1e-9)
  beta <- selected[, "beta"]

  # Weighted season fallback is used only when strict_dates = FALSE.
  target_counts <- ge_fit$N_seen[r, ]
  if (sum(target_counts) > 0) {
    season <- as.numeric(psi %*% (target_counts / sum(target_counts)))
  } else {
    season <- rowMeans(psi)
  }
  ge <- matrix(season, nrow = length(pass_dates_d), ncol = B, byrow = TRUE)

  if (is.null(daily_spill)) {
    for (d in which(valid_s)) ge[d, ] <- psi[, s[d]]
  } else {
    if (!is.data.frame(daily_spill) ||
        !all(c("Date", "spill.per") %in% names(daily_spill))) {
      stop("daily_spill must contain Date and spill.per columns.", call. = FALSE)
    }
    spill_dates <- parse_dates(daily_spill$Date, "daily spill date")
    if (anyDuplicated(spill_dates)) {
      stop("daily_spill contains duplicate dates.", call. = FALSE)
    }
    spill <- suppressWarnings(as.numeric(daily_spill$spill.per))
    spill_day <- spill[match(pass_dates_d, spill_dates)]
    missing_spill <- valid_s & !is.finite(spill_day)
    if (strict_daily_spill && any(missing_spill)) {
      stop(sum(missing_spill),
           " modeled passage date(s) lack finite daily spill values.",
           call. = FALSE)
    }
    for (d in which(valid_s)) {
      ss <- s[d]
      x <- spill_day[d]
      if (!is.finite(x)) {
        x_std <- ge_fit$lgr_spill_std[ss]
      } else {
        x_std <- (x - ge_fit$spill_mean) / ge_fit$spill_sd
      }
      ge[d, ] <- stats::plogis(
        stats::qlogis(psi[, ss]) +
          beta * (x_std - ge_fit$lgr_spill_std[ss]))
    }
  }
  ge <- pmin(pmax(ge, 0), 1)
  colnames(ge) <- paste0("boot_", seq_len(B))
  data.frame(SampleEndDate = pass_dates, ge, check.names = FALSE)
}
