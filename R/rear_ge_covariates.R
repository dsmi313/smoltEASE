#' Daily PIT passage weights by rear group
#'
#' Counts each eligible tag once, on its first main GRJ/GRS passage date.
#' GRJ takes route priority for a tag detected at both sites. Transported fish
#' must already have their bypass passage records in source = "main"; loading
#' records are not counted. The caller must first exclude dam-tagged fish,
#' adult releases, reverse movements, and unusable transport histories.
#' Weights are G + S / p_mean, using a fixed posterior mean weekly GRS
#' detection probability. Uncertainty in these plug-in weights is not propagated.
#'
#' @param events Cleaned data frame with tag, rear_type, source, site, det_time.
#'   det_time must be POSIXct; these are passage times, not loading times.
#' @param p_week Data frame with rear_type, week, p_mean; unique rear-week rows.
#' @param weeks Optional ISO weeks to keep, in model order.
#' @param tz Time zone used to turn passage times into calendar dates. Default
#'   "UTC" matches PTAGIS local clock strings parsed with a UTC label.
#' @return Data frame with rear_type, date, week, G, S, p_mean, weight.
#' @export
calc_daily_passage_weights <- function(events, p_week, weeks = NULL, tz = "UTC") {
  .ge_require_columns(events, c("tag", "rear_type", "source", "site", "det_time"), "events")
  .ge_require_columns(p_week, c("rear_type", "week", "p_mean"), "p_week")
  if (!inherits(events$det_time, "POSIXct")) stop("events$det_time must be POSIXct.", call. = FALSE)
  if (!is.character(tz) || length(tz) != 1L || is.na(tz) || !nzchar(tz) ||
      !tz %in% c("UTC", "GMT", OlsonNames())) stop("tz must be a recognized time zone.", call. = FALSE)
  .ge_check_weeks(p_week$week)
  if (!is.null(weeks)) .ge_check_weeks(weeks, unique = TRUE)
  p_week$rear_type <- as.character(p_week$rear_type)
  if (anyNA(p_week$rear_type) || any(!nzchar(p_week$rear_type)) ||
      anyDuplicated(p_week[c("rear_type", "week")])) stop("p_week needs unique, nonmissing rear-week rows.", call. = FALSE)
  if (!is.numeric(p_week$p_mean) || any(!is.finite(p_week$p_mean) |
      p_week$p_mean <= 0 | p_week$p_mean > 1)) stop("p_mean must be finite and in (0, 1].", call. = FALSE)
  use <- !is.na(events$source) & events$source == "main" &
    !is.na(events$site) & events$site %in% c("GRJ", "GRS")
  x <- events[use, , drop = FALSE]
  x$tag <- as.character(x$tag); x$rear_type <- as.character(x$rear_type)
  if (anyNA(x[c("tag", "rear_type", "det_time")]) ||
      any(!nzchar(x$tag) | !nzchar(x$rear_type))) stop("Main passage records need tag, rear_type and det_time.", call. = FALSE)
  if (anyDuplicated(unique(x[c("tag", "rear_type")])$tag)) stop("A tag belongs to more than one rear group.", call. = FALSE)
  empty <- data.frame(rear_type = character(), date = as.Date(character()),
                      week = integer(), G = integer(), S = integer(),
                      p_mean = numeric(), weight = numeric())
  if (!nrow(x)) return(empty)
  fish <- dplyr::summarise(dplyr::group_by(x, tag, rear_type),
    date = as.Date(min(det_time), tz = tz),
    route = if (any(site == "GRJ")) "G" else "S", .groups = "drop")
  fish$week <- as.integer(format(fish$date, "%V"))
  if (!is.null(weeks)) fish <- fish[fish$week %in% weeks, , drop = FALSE]
  if (!nrow(fish)) return(empty)
  .ge_single_year(fish$date)
  out <- dplyr::summarise(dplyr::group_by(fish, rear_type, date, week),
                          G = sum(route == "G"), S = sum(route == "S"), .groups = "drop")
  out <- dplyr::left_join(out, p_week[c("rear_type", "week", "p_mean")], by = c("rear_type", "week"))
  if (anyNA(out$p_mean)) stop("Some passage rear-weeks lack p_mean.", call. = FALSE)
  out$weight <- out$G + out$S / out$p_mean
  if (any(!is.finite(out$weight))) stop("Passage weight overflowed; inspect p_mean.", call. = FALSE)
  as.data.frame(dplyr::arrange(out, rear_type, date))
}

#' Apply a configurable spill-decay window
#'
#' For date d, computes sum(decay^k * spill[d-k]) / sum(decay^k), with
#' k = 0,...,lag_days. At the beginning of the supplied series only available
#' preceding days are used. Provide earlier dates for a full initial window.
#' decay = 0 gives same-day spill; decay = 1 gives an equally weighted window.
#' Calendar dates must be consecutive. Outflow is always left at its daily value.
#'
#' @param daily_spill Data frame with Date, spill.per (percent, 0 to 100),
#'   outflow (nonnegative daily mean flow rate). One row per consecutive date.
#' @param decay Spill-decay factor between zero and one. Default 0.5.
#' @param lag_days Number of preceding calendar days; default 7.
#' @return Sorted daily_spill with spill.per replaced by the decayed spill.
#'   Other columns are retained. Attribute spill_decay records the settings.
#' @export
apply_spill_decay <- function(daily_spill, decay = 0.5, lag_days = 7L) {
  x <- .ge_daily_covariates(daily_spill)
  if (!is.numeric(decay) || length(decay) != 1L || !is.finite(decay) ||
      decay < 0 || decay > 1) stop("decay must be one number in [0, 1].", call. = FALSE)
  if (!is.numeric(lag_days) || length(lag_days) != 1L || !is.finite(lag_days) ||
      lag_days < 0 || lag_days != floor(lag_days) || lag_days > .Machine$integer.max) {
    stop("lag_days must be one nonnegative integer.", call. = FALSE)
  }
  if (any(diff(as.numeric(x$Date)) != 1)) stop("daily_spill must have consecutive calendar dates; fill operations gaps first.", call. = FALSE)
  original <- x$spill.per
  x$spill.per <- vapply(seq_along(original), function(i) {
    k <- seq.int(0L, min(lag_days, i - 1L))
    w <- decay^k
    sum(w * original[i - k]) / sum(w)
  }, numeric(1))
  attr(x, "spill_decay") <- list(decay = decay, lag_days = as.integer(lag_days))
  x
}

#' Build rear-weighted weekly GE covariate matrices
#'
#' Uses daily PIT passage weights to average the supplied daily spill and
#' outflow separately by rear and ISO week. Does not apply spill decay itself.
#' Rear-weeks with no positive weight keep their shared weekly reference;
#' the audit identifies every such fallback. Matrices are in exactly the
#' supplied rear_levels and weeks order, ready for psi_spill and psi_outflow.
#'
#' @param daily_weights Result of [calc_daily_passage_weights()], or a data
#'   frame with rear_type, date, week, weight; unique rear-date rows.
#' @param daily_spill Daily covariates, possibly from [apply_spill_decay()].
#' @param weeks Unique fitted ISO weeks, in model order.
#' @param rear_levels Unique rear groups, in model order.
#' @param shared_spill,shared_outflow Finite shared weekly reference vectors,
#'   one per week. Spill must be in [0,100] and outflow nonnegative.
#' @return List with psi_spill, psi_outflow, and audit (rear, week,
#'   total_weight, n_days, fallback, spill.per, outflow).
#' @export
build_weighted_weekly_covariates <- function(daily_weights, daily_spill, weeks,
                                            rear_levels, shared_spill, shared_outflow) {
  .ge_require_columns(daily_weights, c("rear_type", "date", "week", "weight"), "daily_weights")
  .ge_check_weeks(weeks, unique = TRUE)
  if (!is.character(rear_levels) || !length(rear_levels) || anyNA(rear_levels) ||
      any(!nzchar(rear_levels)) || anyDuplicated(rear_levels)) stop("rear_levels must be unique nonempty rear names.", call. = FALSE)
  if (!is.numeric(shared_spill) || length(shared_spill) != length(weeks) ||
      any(!is.finite(shared_spill) | shared_spill < 0 | shared_spill > 100) ||
      !is.numeric(shared_outflow) || length(shared_outflow) != length(weeks) ||
      any(!is.finite(shared_outflow) | shared_outflow < 0)) stop("Shared covariates need one valid value per fitted week.", call. = FALSE)
  ds <- .ge_daily_covariates(daily_spill)
  w <- as.data.frame(daily_weights)
  w$date <- .ge_dates(w$date, "daily_weights$date")
  w$rear_type <- as.character(w$rear_type)
  if (anyNA(w$rear_type) || any(!w$rear_type %in% rear_levels) ||
      anyDuplicated(w[c("rear_type", "date")])) stop("Weights need unique rear-date rows and known rear groups.", call. = FALSE)
  .ge_check_weeks(w$week)
  if (any(w$week != as.integer(format(w$date, "%V")))) stop("Weight weeks do not match their ISO passage dates.", call. = FALSE)
  if (!is.numeric(w$weight) || any(!is.finite(w$weight) | w$weight < 0)) stop("Weights must be finite and nonnegative.", call. = FALSE)
  w <- w[w$week %in% weeks, , drop = FALSE]
  .ge_single_year(c(w$date, ds$Date))
  ii <- match(w$date, ds$Date)
  if (anyNA(ii)) stop("Some passage dates lack daily spill or outflow.", call. = FALSE)
  nr <- length(rear_levels); ns <- length(weeks)
  spill <- matrix(shared_spill, nr, ns, byrow = TRUE, dimnames = list(rear_levels, as.character(weeks)))
  flow <- matrix(shared_outflow, nr, ns, byrow = TRUE, dimnames = dimnames(spill))
  audit <- vector("list", nr * ns)
  for (r in seq_len(nr)) for (ss in seq_len(ns)) {
    j <- which(w$rear_type == rear_levels[r] & w$week == weeks[ss] & w$weight > 0)
    total <- sum(w$weight[j])
    if (!is.finite(total)) stop("Total passage weight overflowed.", call. = FALSE)
    if (total > 0) {
      spill[r, ss] <- sum(w$weight[j] / total * ds$spill.per[ii[j]])
      flow[r, ss] <- sum(w$weight[j] / total * ds$outflow[ii[j]])
    }
    audit[[(r - 1L) * ns + ss]] <- data.frame(rear = rear_levels[r], week = weeks[ss],
      total_weight = total, n_days = length(j), fallback = total == 0,
      spill.per = spill[r, ss], outflow = flow[r, ss])
  }
  list(psi_spill = spill, psi_outflow = flow, audit = do.call(rbind, audit))
}

#' Prepare passage-weighted GE covariates with configurable spill decay
#'
#' Combines [calc_daily_passage_weights()], [apply_spill_decay()] and
#' [build_weighted_weekly_covariates()]. Extracts posterior mean weekly p
#' from a saved reference rear fit, normally the production fit. The reference
#' fit supplies rear/week order and fallback covariates, not new fitted models.
#' Use the returned daily_spill both for weighting and for daily GE generation.
#' Changing decay requires a corresponding refit; changing daily covariates
#' alone on an old fit does not perform a new decay calibration.
#'
#' @param events Cleaned PIT passage records; see [calc_daily_passage_weights()].
#' @param daily_spill Unlagged daily Date, spill.per, outflow.
#' @param ge_fit Saved reference rear fit with samples, rear_levels, weeks,
#'   lgr_spill_pct and the shared weekly outflow vector.
#' @param decay Spill-decay factor in [0,1]. Default 0.5.
#' @param lag_days Number of preceding calendar days. Default 7.
#' @param tz Clock-time zone for passage dates. Default "UTC".
#' @return List with psi_spill, psi_outflow, audit, daily_weights, p_week,
#'   daily_spill, and settings. The matrices can be passed directly to
#'   [fit_ge_rear_daynight()] or [fit_ge_rear_model2()].
#' @export
prep_rear_ge_covariates <- function(events, daily_spill, ge_fit, decay = 0.5,
                                    lag_days = 7L, tz = "UTC") {
  if (!is.list(ge_fit) || !all(c("samples", "rear_levels", "weeks", "lgr_spill_pct") %in% names(ge_fit))) {
    stop("ge_fit must contain samples, rear_levels, weeks, lgr_spill_pct.", call. = FALSE)
  }
  rears <- ge_fit$rear_levels; weeks <- ge_fit$weeks
  .ge_check_weeks(weeks, unique = TRUE)
  samples <- ge_fit$samples
  if (inherits(samples, "mcmc") || is.matrix(samples)) samples <- list(samples)
  if (!is.list(samples) || !length(samples)) stop("ge_fit has no posterior chains.", call. = FALSE)
  mats <- lapply(samples, as.matrix)
  if (any(vapply(mats, nrow, integer(1)) == 0L) ||
      !all(vapply(mats, function(x) identical(colnames(x), colnames(mats[[1L]])), logical(1)))) {
    stop("Reference chains are empty or parameter columns differ.", call. = FALSE)
  }
  post <- do.call(rbind, mats)
  p_week <- do.call(rbind, lapply(seq_along(rears), function(r) {
    cols <- sprintf("p[%d,%d]", r, seq_along(weeks))
    if (!all(cols %in% colnames(post))) stop("Reference fit lacks rear-week p draws.", call. = FALSE)
    pp <- post[, cols, drop = FALSE]
    if (any(!is.finite(pp) | pp <= 0 | pp > 1)) stop("Reference p draws must be in (0,1].", call. = FALSE)
    data.frame(rear_type = rears[r], week = weeks, p_mean = colMeans(pp))
  }))
  shared_flow <- ge_fit$outflow
  if (is.null(shared_flow) && is.matrix(ge_fit$psi_outflow)) {
    # Day/night fits store the shared series through its standardization.
    # A weighted rear series is not an appropriate shared fallback.
    stop("Reference fit lacks shared outflow; use the production rear fit.", call. = FALSE)
  }
  weights <- calc_daily_passage_weights(events, p_week, weeks = weeks, tz = tz)
  daily <- apply_spill_decay(daily_spill, decay = decay, lag_days = lag_days)
  weekly <- build_weighted_weekly_covariates(weights, daily, weeks, rears,
                                             ge_fit$lgr_spill_pct, shared_flow)
  c(weekly, list(daily_weights = weights, p_week = p_week, daily_spill = daily,
                settings = list(decay = decay, lag_days = as.integer(lag_days),
                  tz = tz, weight_method = "G + S / posterior mean weekly p", outflow_lagged = FALSE)))
}

.ge_require_columns <- function(x, columns, label) {
  if (!is.data.frame(x) || !all(columns %in% names(x))) stop(label, " needs columns: ", paste(columns, collapse = ", "), ".", call. = FALSE)
}
.ge_check_weeks <- function(x, unique = FALSE) {
  if (!is.numeric(x) || any(!is.finite(x) | x < 1 | x > 53 | x != floor(x)) ||
      (unique && (!length(x) || anyDuplicated(x)))) stop("weeks must be integer ISO weeks in 1:53, unique when defining model order.", call. = FALSE)
}
.ge_dates <- function(x, label) {
  if (!inherits(x, "Date")) {
    if (!is.character(x)) stop(label, " must be Date or YYYY-MM-DD strings.", call. = FALSE)
    if (anyNA(x) || any(!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x))) stop(label, " contains invalid dates.", call. = FALSE)
    raw <- x
    x <- as.Date(raw, format = "%Y-%m-%d")
    if (anyNA(x) || any(format(x, "%Y-%m-%d") != raw)) stop(label, " contains invalid dates.", call. = FALSE)
  }
  if (anyNA(x)) stop(label, " contains missing or invalid dates.", call. = FALSE)
  x
}
.ge_single_year <- function(dates) {
  if (length(unique(format(dates, "%Y"))) > 1L) stop("Use one calendar year per call to avoid mixing ISO weeks across years.", call. = FALSE)
}
.ge_daily_covariates <- function(x) {
  .ge_require_columns(x, c("Date", "spill.per", "outflow"), "daily_spill")
  if (!nrow(x)) stop("daily_spill is empty.", call. = FALSE)
  x <- as.data.frame(x); x$Date <- .ge_dates(x$Date, "daily_spill$Date")
  if (anyDuplicated(x$Date)) stop("daily_spill contains duplicate dates.", call. = FALSE)
  if (!is.numeric(x$spill.per) || any(!is.finite(x$spill.per) | x$spill.per < 0 | x$spill.per > 100) ||
      !is.numeric(x$outflow) || any(!is.finite(x$outflow) | x$outflow < 0)) stop("Daily spill must be finite in [0,100]; outflow finite and nonnegative.", call. = FALSE)
  x <- x[order(x$Date), , drop = FALSE]; rownames(x) <- NULL
  x
}

utils::globalVariables(c("tag", "rear_type", "det_time", "site", "route", "date", "week"))
