#' @title Generate per-day GE posterior draws for use in SCRAPI2
#'
#' @description Post-processes MCMC output from \code{\link{fit_ge_model}} into a
#'   data frame of guidance efficiency (GE) posterior draws for the
#'   \code{geDraws} argument of \code{\link{SCRAPI2}}. Each \code{boot_*} column
#'   is one posterior draw of daily GE; \code{SCRAPI2} consumes one column per
#'   bootstrap iteration, exactly mirroring how GSI draws are propagated.
#'
#'   GE is the route-selection probability \eqn{\psi} from the multistate model
#'   (P(JBS route)). Strata estimated in the JAGS likelihood use their \eqn{\psi}
#'   posterior directly. Strata with no PIT data (\code{n_pool == 0}) borrow
#'   strength from the spill regression \emph{and} the hierarchical process
#'   variance: their \eqn{\psi} is drawn as
#'   \code{plogis(alpha + beta * spill_std + eps)} with
#'   \code{eps ~ N(0, sigma_psi)} per draw (partial pooling, after Hance et al.).
#'   Spill is standardised with the \code{spill_mean}/\code{spill_sd} stored in
#'   \code{ge_fit}, matching the scale on which \code{alpha}/\code{beta} were fit.
#'
#'   When \code{daily_spill} is supplied, GE is resolved to a \emph{daily} value
#'   rather than held constant within a stratum. For each day the stratum's
#'   fitted process deviation \code{eps} (recovered per draw as
#'   \code{logit(psi) - (alpha + beta * stratum_spill_std)}) is held fixed while
#'   the spill covariate is replaced by that day's spill:
#'   \deqn{GE_d = logit^{-1}(\alpha + \beta \cdot spillstd_d + eps_{stratum(d)}).}
#'   This lets GE track day-to-day spill within a stratum while staying anchored
#'   to the stratum's mark-recapture estimate. All draws for a given day share
#'   the same \code{(alpha, beta, eps)}, so daily GE values are correlated
#'   through the coefficients rather than carrying independent daily noise (which
#'   would understate uncertainty).
#'
#'   When \code{clip_to_ci = TRUE}, sampling is restricted to complete MCMC
#'   rows in which \emph{every} observed stratum's \eqn{\psi} falls within that
#'   stratum's own marginal \code{ci_lower}/\code{ci_upper} empirical quantile.
#'   Because this is an intersection across strata, the retained set is not
#'   "95% of the posterior" -- with multiple strata the retained fraction can
#'   be well below 95%. Selecting whole rows (rather than sampling each
#'   stratum's \eqn{\psi} independently) preserves the fitted joint dependence
#'   among strata and shared parameters (\code{alpha}, \code{beta},
#'   \code{sigma_psi}); the same \code{(alpha, beta, eps)} used for a given
#'   draw's daily projection came from one real posterior sample, not a
#'   composite of unrelated ones.
#'
#' @param ge_fit list returned by \code{\link{fit_ge_model}}.
#' @param strat_assign data frame mapping each calendar date to a trap stratum.
#'   Same object passed to \code{\link{prep_ge_data}}. Accepts the
#'   \code{Week}/\code{Collapse} format or the \code{date}/\code{stratum}/
#'   \code{stratum_idx} format. To estimate GE \strong{weekly}, pass a mapping in
#'   which \code{Collapse} equals \code{Week} (one stratum per week) -- the same
#'   mapping must be passed to \code{\link{prep_ge_data}}.
#' @param pass_dates vector of dates corresponding to rows of \code{passageData}
#'   (i.e. \code{passageData$SampleEndDate}). Accepts \code{Date} objects or
#'   character strings in any of the formats \code{"\%m/\%d/\%Y"},
#'   \code{"\%Y-\%m-\%d"}, or \code{"\%d/\%m/\%Y"}. No \code{as.Date()} call
#'   is needed by the caller.
#' @param B number of GE draws to return. Default 2000.
#' @param daily_spill optional data frame of daily LGR spill with columns
#'   \code{Date} (Date or character) and \code{spill.per} (numeric percentage) --
#'   the same spill series passed to \code{\link{prep_ge_data}} as
#'   \code{spill_data}. When supplied, GE is resolved to daily resolution within
#'   each stratum using the fitted spill regression (see Details). When
#'   \code{NULL} (default) GE is held constant within each stratum (the original
#'   behaviour). Days with no matching spill value fall back to the stratum mean.
#' @param clip_to_ci logical. When \code{TRUE}, restricts sampling to complete
#'   MCMC rows where every observed stratum's \eqn{\psi} falls within that
#'   stratum's own \code{ci_lower}/\code{ci_upper} marginal empirical quantile
#'   (the intersection of the per-stratum central intervals, not "95% of the
#'   posterior" -- see Details). Default \code{FALSE} (draws from the full
#'   posterior).
#' @param ci_lower,ci_upper quantiles defining each stratum's marginal
#'   interval when \code{clip_to_ci = TRUE}. Defaults 0.025 and 0.975.
#' @param seed optional integer. When supplied, \code{set.seed(seed)} is
#'   called immediately before row sampling, making the returned draws (and
#'   the downstream \code{rnorm()} draws for zero-data strata) reproducible
#'   across repeated calls with identical inputs. Default \code{NULL} (no
#'   seeding, original behaviour).
#'
#' @return A data frame with \code{length(pass_dates)} rows and \code{B + 1}
#'   columns: \code{SampleEndDate}, \code{boot_1}, \ldots, \code{boot_B}.
#'   All GE values are in \eqn{[0, 1]}.
#'
#' @importFrom stats plogis qlogis rnorm quantile
#' @export
generate_ge_draws <- function(ge_fit,
                              strat_assign,
                              pass_dates,
                              B           = 2000,
                              daily_spill = NULL,
                              clip_to_ci  = FALSE,
                              ci_lower    = 0.025,
                              ci_upper    = 0.975,
                              seed        = NULL) {

  # Internal date parser: handles Date objects and the three character formats
  # that appear in SCRAPI/smoltEASE pipelines without requiring the caller to
  # run as.Date() themselves.
  parse_dates <- function(x) {
    if (inherits(x, "Date")) return(x)
    as.Date(x, tryFormats = c("%m/%d/%Y", "%Y-%m-%d", "%d/%m/%Y"))
  }

  ge_data    <- ge_fit$ge_data
  obs_strata <- ge_fit$obs_strata
  spill_mean <- ge_fit$spill_mean
  spill_sd   <- ge_fit$spill_sd

  # --- Normalise strat_assign to date/stratum/stratum_idx format ---
  if (all(c("Week", "Collapse") %in% names(strat_assign)) &&
      !("date" %in% names(strat_assign))) {

    all_dates <- sort(unique(parse_dates(pass_dates)))
    wk_num    <- as.integer(format(all_dates, "%V"))

    week_to_strat <- strat_assign
    names(week_to_strat)[names(week_to_strat) == "Collapse"] <- "stratum"
    week_to_strat$stratum_idx <- as.integer(
      factor(week_to_strat$stratum, levels = sort(unique(week_to_strat$stratum)))
    )

    strat_assign <- merge(
      data.frame(date = all_dates, Week = wk_num),
      week_to_strat,
      by = "Week", all.x = FALSE
    )[, c("date", "stratum", "stratum_idx")]
  }

  # obs_strata must preserve ge_data row order so JAGS psi[1..k] columns align
  stopifnot(all(diff(obs_strata$stratum_idx) > 0))

  # --- Extract psi posterior draws + regression / process-error params ---
  samp_mat  <- do.call(rbind, lapply(ge_fit$samples, as.matrix))
  psi_cols  <- grep("^psi\\[", colnames(samp_mat))
  psi_mat   <- samp_mat[, psi_cols, drop = FALSE]   # n_mcmc x n_obs_strata
  alpha_col <- grep("^alpha$",     colnames(samp_mat))
  beta_col  <- grep("^beta$",      colnames(samp_mat))
  sigma_col <- grep("^sigma_psi$", colnames(samp_mat))

  n_all    <- nrow(ge_data)
  obs_idx  <- match(obs_strata$stratum_idx, ge_data$stratum_idx)
  zero_idx <- setdiff(seq_len(n_all), obs_idx)

  # Standardise stratum spill exactly as fit_ge_model did
  strat_spill_std <- (ge_data$spill_val - spill_mean) / spill_sd

  # --- Row selection ---
  # clip_to_ci = TRUE restricts to complete MCMC rows where every observed
  # stratum's psi falls within that stratum's own marginal central interval
  # (ci_lower/ci_upper). This is the intersection of the stratum-specific
  # intervals, not "95% of the joint posterior" -- with multiple strata the
  # retained fraction can be well below 95%. Preserving whole rows keeps the
  # fitted correlation among strata and shared parameters (alpha, beta,
  # sigma_psi) intact; sampling each stratum independently would break it.
  psi_lo <- apply(psi_mat, 2, quantile, probs = ci_lower, na.rm = TRUE)
  psi_hi <- apply(psi_mat, 2, quantile, probs = ci_upper, na.rm = TRUE)

  if (clip_to_ci) {
    inside_ci <- sweep(psi_mat, 2, psi_lo, `>=`) &
                 sweep(psi_mat, 2, psi_hi, `<=`)
    valid_rows <- which(rowSums(inside_ci) == ncol(inside_ci))
    if (length(valid_rows) == 0)
      stop("No complete posterior rows satisfy the requested CI restriction.")
  } else {
    valid_rows <- seq_len(nrow(samp_mat))
  }

  if (!is.null(seed)) set.seed(seed)
  idx <- sample(valid_rows, size = B, replace = length(valid_rows) < B)

  selected  <- samp_mat[idx, , drop = FALSE]
  psi_draws <- selected[, psi_cols, drop = FALSE]
  alpha_v   <- selected[, alpha_col]
  beta_v    <- selected[, beta_col]
  sigma_v   <- selected[, sigma_col]

  # --- Per-stratum process deviation (eps) on the logit-psi scale ---
  eps_full    <- matrix(NA_real_, nrow = B, ncol = n_all)
  psi_clamped <- pmin(pmax(psi_draws, 1e-6), 1 - 1e-6)
  for (s in seq_along(obs_idx))
    eps_full[, obs_idx[s]] <- qlogis(psi_clamped[, s]) -
      (alpha_v + beta_v * strat_spill_std[obs_idx[s]])
  for (zi in zero_idx)
    eps_full[, zi] <- rnorm(B, 0, sigma_v)

  # --- Stratum-constant GE draws for fallback / non-daily mode ---
  ge_strat <- matrix(NA_real_, nrow = B, ncol = n_all)
  ge_strat[, obs_idx] <- psi_draws
  for (zi in zero_idx)
    ge_strat[, zi] <- plogis(alpha_v + beta_v * strat_spill_std[zi] + eps_full[, zi])
  ge_strat <- pmin(pmax(ge_strat, 0), 1)

  # Season-mean fallback: n_pool-weighted mean across strata
  pool_wts   <- ifelse(ge_data$n_pool > 0, ge_data$n_pool, 0)
  wt_mat     <- matrix(pool_wts, nrow = B, ncol = n_all, byrow = TRUE)
  wt_rowsums <- rowSums(wt_mat)
  season_ge  <- ifelse(wt_rowsums > 0,
                       rowSums(ge_strat * wt_mat) / wt_rowsums,
                       rowMeans(ge_strat))

  # --- Map strata to daily passageData rows ---
  pass_dates_d  <- parse_dates(pass_dates)
  strat_dates_d <- as.Date(strat_assign$date)
  day_match     <- match(pass_dates_d, strat_dates_d)
  strat_idx_v   <- strat_assign$stratum_idx[day_match]

  valid   <- !is.na(strat_idx_v)
  ge_days <- matrix(season_ge, nrow = B, ncol = length(pass_dates))

  if (is.null(daily_spill)) {
    ge_days[, valid] <- ge_strat[, strat_idx_v[valid]]
  } else {
    ds        <- daily_spill
    ds$Date   <- as.Date(ds$Date)
    spill_day <- ds$spill.per[match(pass_dates_d, ds$Date)]
    for (d in which(valid)) {
      s_col <- strat_idx_v[d]
      x_d   <- spill_day[d]
      if (is.na(x_d)) x_d <- ge_data$spill_val[s_col]
      std_d <- (x_d - spill_mean) / spill_sd
      ge_days[, d] <- plogis(alpha_v + beta_v * std_d + eps_full[, s_col])
    }
    ge_days <- pmin(pmax(ge_days, 0), 1)
  }

  result   <- t(ge_days)
  draws_df <- as.data.frame(result)
  colnames(draws_df) <- paste0("boot_", seq_len(B))
  cbind(data.frame(SampleEndDate = pass_dates), draws_df)
}
