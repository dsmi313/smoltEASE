#' Daily GE draws from a day/night fit
#'
#' Adjusts day and night GE to daily covariates, centered on the fit's rear
#' weekly covariates, then mixes by the inferred night share of fish.
#' The night share is inferred from spilled GRS detections only; guided
#' counts never enter it. Its detected-spill share is shrunk toward the
#' weekly model value with shrink_k pseudo-fish. Days with no detections
#' use that weekly detected-spill share. At reference covariates this
#' recovers the week's night share of fish exactly.
#'
#' Weekly psi_day and psi_night are passage-weighted means of daily GE, not GE
#' at the weekly mean covariate. Because GE is nonlinear in its covariates,
#' shifting the weekly logit by each day's covariate deviation makes the
#' passage-weighted mean of the daily values drift from the weekly estimate
#' (upward when GE is small and varies within the week). When daily_weights is
#' supplied, one logit offset per week and draw, separately for day and night,
#' is added so the weighted mean of daily day (night) GE equals that draw's
#' weekly psi_day (psi_night). Without daily_weights the uncalibrated
#' behaviour of earlier versions is kept.
#' @param ge_fit Result of [fit_ge_rear_daynight()].
#' @param rear_type Rear group; defaults to the target rear.
#' @param pass_dates Dates in output row order.
#' @param B Number of posterior columns.
#' @param daily_spill Data frame with Date, spill.per, outflow, using the
#'   same covariate definition (including any lag) as the fit.
#' @param spill_counts Data frame with Date, S_day, S_night for this rear.
#' @param shrink_k Positive pseudo-count shrinking daily night share; default 10.
#' @param strict_daily_spill Stop if a date lacks spill or outflow.
#' @param seed Posterior sampling seed.
#' @param daily_weights Optional daily passage weights, normally
#'   `prep_rear_ge_covariates(...)$daily_weights`: rear_type, date, week, weight.
#'   Supply it to calibrate daily GE to the weekly estimates (recommended).
#' @return Data frame with SampleEndDate and boot_1 through boot_B. Attributes
#'   draw_id, night_fraction, and p_effective preserve draw alignment and
#'   effective GRS detection probabilities for route checks. With
#'   daily_weights, attribute calibration_offset holds the day and night
#'   week-by-draw logit offsets.
#' @export
generate_rear_ge_draws_daynight <- function(
    ge_fit, rear_type = NULL, pass_dates, B = 2000, daily_spill, spill_counts,
    shrink_k = 10, strict_daily_spill = TRUE, seed = NULL, daily_weights = NULL) {
  if (!inherits(ge_fit, "ge_rear_daynight")) {
    stop("ge_fit must come from fit_ge_rear_daynight().", call. = FALSE)
  }
  if (is.null(rear_type)) rear_type <- ge_fit$target_rear
  r <- match(rear_type, ge_fit$rear_levels)
  if (is.na(r)) stop("Unknown rear_type.", call. = FALSE)
  if (!all(c("Date", "spill.per", "outflow") %in% names(daily_spill))) {
    stop("daily_spill needs Date, spill.per, outflow.", call. = FALSE)
  }
  if (!all(c("Date", "S_day", "S_night") %in% names(spill_counts))) {
    stop("spill_counts needs Date, S_day, S_night.", call. = FALSE)
  }
  if (!is.numeric(shrink_k) || length(shrink_k) != 1L ||
      !is.finite(shrink_k) || shrink_k <= 0) {
    stop("shrink_k must be one positive finite number.", call. = FALSE)
  }
  dates <- as.Date(pass_dates)
  s <- match(as.integer(format(dates, "%V")), ge_fit$weeks)
  if (anyNA(s)) stop(sum(is.na(s)), " date(s) fall outside the fitted weeks.", call. = FALSE)
  mat <- do.call(rbind, lapply(ge_fit$samples, as.matrix))
  if (!is.null(seed)) set.seed(seed)
  draw_id <- sample(seq_len(nrow(mat)), B, replace = nrow(mat) < B)
  sel <- mat[draw_id, , drop = FALSE]
  col <- function(nm, ss) sel[, sprintf("%s[%d,%d]", nm, r, ss)]
  clamp <- function(x) pmin(pmax(x, 1e-9), 1 - 1e-9)
  beta <- sel[, "beta"]
  beta_outflow <- sel[, "beta_outflow"]
  ds_dates <- as.Date(daily_spill$Date)
  spill_day <- as.numeric(daily_spill$spill.per)[match(dates, ds_dates)]
  flow_day <- as.numeric(daily_spill$outflow)[match(dates, ds_dates)]
  if (strict_daily_spill && any(!is.finite(spill_day) | !is.finite(flow_day))) {
    stop("Some dates lack daily spill or outflow.", call. = FALSE)
  }
  sc_dates <- as.Date(spill_counts$Date)
  S_n <- as.numeric(spill_counts$S_night)[match(dates, sc_dates)]
  S_d <- as.numeric(spill_counts$S_day)[match(dates, sc_dates)]
  S_n[is.na(S_n)] <- 0
  S_d[is.na(S_d)] <- 0
  p_off <- if (is.null(ge_fit$p_night_offset)) 0 else ge_fit$p_night_offset
  day_adj <- function(dd, ss) {
    i <- match(dd, ds_dates)
    sp <- as.numeric(daily_spill$spill.per)[i]; fl <- as.numeric(daily_spill$outflow)[i]
    zs <- (sp - ge_fit$spill_mean) / ge_fit$spill_sd - ge_fit$psi_spill_std[r, ss]
    zf <- (fl - ge_fit$outflow_mean) / ge_fit$outflow_sd - ge_fit$psi_outflow_std[r, ss]
    outer(zs, beta) + outer(zf, beta_outflow)
  }
  nw <- length(ge_fit$weeks)
  off_d <- off_n <- matrix(0, nw, B)
  if (!is.null(daily_weights)) {
    if (!is.data.frame(daily_weights) ||
        !all(c("rear_type", "date", "week", "weight") %in% names(daily_weights))) {
      stop("daily_weights needs rear_type, date, week, weight.", call. = FALSE)
    }
    w <- daily_weights[as.character(daily_weights$rear_type) == rear_type &
                         is.finite(daily_weights$weight) & daily_weights$weight > 0, , drop = FALSE]
    w$date <- as.Date(w$date)
    if (anyDuplicated(w$date)) stop("daily_weights has duplicate dates for this rear.", call. = FALSE)
    for (ss in seq_len(nw)) {
      k <- which(w$week == ge_fit$weeks[ss])
      if (!length(k)) next
      i <- match(w$date[k], ds_dates)
      ok <- !is.na(i) & is.finite(as.numeric(daily_spill$spill.per)[i]) &
        is.finite(as.numeric(daily_spill$outflow)[i])
      if (strict_daily_spill && !all(ok)) {
        stop("Some weighted passage dates lack daily spill or outflow.", call. = FALSE)
      }
      k <- k[ok]
      if (!length(k)) next
      A <- day_adj(w$date[k], ss)
      pd <- clamp(col("psi_day", ss)); pn <- clamp(col("psi_night", ss))
      off_d[ss, ] <- .ge_calibrate_offset(sweep(A, 2, stats::qlogis(pd), "+"), w$weight[k], pd)
      off_n[ss, ] <- .ge_calibrate_offset(sweep(A, 2, stats::qlogis(pn), "+"), w$weight[k], pn)
    }
  }
  ge <- matrix(NA_real_, length(dates), B)
  f_mat <- matrix(NA_real_, length(dates), B)
  p_eff <- matrix(NA_real_, length(dates), B)
  for (d in seq_along(dates)) {
    ss <- s[d]
    pd <- clamp(col("psi_day", ss))
    pn <- clamp(col("psi_night", ss))
    a <- col("a", ss)
    p_d <- clamp(col("p", ss))
    p_n <- if (isTRUE(ge_fit$joint_grs_detection)) clamp(col("p_night", ss)) else
      stats::plogis(stats::qlogis(p_d) + p_off)
    adj <- if (is.finite(spill_day[d]) && is.finite(flow_day[d])) {
      beta * ((spill_day[d] - ge_fit$spill_mean) / ge_fit$spill_sd - ge_fit$psi_spill_std[r, ss]) +
        beta_outflow * ((flow_day[d] - ge_fit$outflow_mean) / ge_fit$outflow_sd -
                          ge_fit$psi_outflow_std[r, ss])
    } else 0
    gd <- stats::plogis(stats::qlogis(pd) + adj + off_d[ss, ])
    gn <- stats::plogis(stats::qlogis(pn) + adj + off_n[ss, ])
    s_week <- a * (1 - pn) * p_n / (a * (1 - pn) * p_n + (1 - a) * (1 - pd) * p_d)
    s_d <- (S_n[d] + shrink_k * s_week) / (S_n[d] + S_d[d] + shrink_k)
    f <- s_d * (1 - gd) * p_d / (s_d * (1 - gd) * p_d + (1 - s_d) * (1 - gn) * p_n)
    f_mat[d, ] <- f
    ge[d, ] <- f * gn + (1 - f) * gd
    sp_n <- f * (1 - gn)
    sp_d <- (1 - f) * (1 - gd)
    p_eff[d, ] <- (sp_n * p_n + sp_d * p_d) / (sp_n + sp_d)
  }
  ge <- pmin(pmax(ge, 0), 1)
  colnames(ge) <- paste0("boot_", seq_len(B))
  out <- data.frame(SampleEndDate = pass_dates, ge, check.names = FALSE)
  attr(out, "draw_id") <- draw_id
  attr(out, "night_fraction") <- f_mat
  attr(out, "p_effective") <- p_eff
  if (!is.null(daily_weights)) attr(out, "calibration_offset") <- list(day = off_d, night = off_n)
  out
}

# One logit offset per draw so the weighted mean of plogis(L + offset) over a
# week's days equals target. L is days by draws; Newton steps, capped at 2.
.ge_calibrate_offset <- function(L, w, target) {
  w <- w / sum(w)
  cc <- rep(0, ncol(L))
  for (i in seq_len(100L)) {
    P <- stats::plogis(sweep(L, 2, cc, "+"))
    step <- (colSums(w * P) - target) / pmax(colSums(w * P * (1 - P)), 1e-12)
    step <- pmin(pmax(step, -2), 2)
    cc <- cc - step
    if (max(abs(step)) < 1e-10) break
  }
  cc
}
