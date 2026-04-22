#' @title Generate per-day GE posterior draws for use in SCRAPI2
#'
#' @description Post-processes MCMC output from \code{\link{fit_ge_model}} into a
#'   data frame of guidance efficiency (GE) draws for the \code{geDraws} argument
#'   of \code{\link{SCRAPI2}}.
#'
#'   GE is computed outside JAGS to avoid numerical instability:
#'   \deqn{GE_s = \frac{n_{GRJ,s}}{n_{GRJ,s} + n_{GRS,s} / \psi_s}}
#'
#'   Strata excluded from the JAGS likelihood (\code{n_pool == 0}) receive psi
#'   values predicted by the regression: \code{plogis(alpha + beta * spill_val)}.
#'   Days not matched to any stratum receive a n_pool-weighted season-mean GE.
#'
#' @param ge_fit list returned by \code{\link{fit_ge_model}}.
#' @param strat_assign data frame mapping each calendar date to a trap stratum.
#'   Same object passed to \code{\link{prep_ge_data}}. Required columns:
#'   \code{date} (Date), \code{stratum}, \code{stratum_idx} (integer).
#' @param pass_dates vector of dates (character or Date) corresponding to rows of
#'   \code{passageData} (i.e. \code{passageData$SampleEndDate}).
#' @param B number of GE draws to return. Default 2000.
#'
#' @return A data frame with \code{length(pass_dates)} rows and \code{B + 1}
#'   columns: \code{SampleEndDate}, \code{boot_1}, \ldots, \code{boot_B}.
#'   All GE values are in \eqn{[0, 1]}.
#'
#' @importFrom stats plogis
#' @export
generate_ge_draws <- function(ge_fit,
                               strat_assign,
                               pass_dates,
                               B = 2000) {

  ge_data    <- ge_fit$ge_data
  obs_strata <- ge_fit$obs_strata

  # --- Normalise strat_assign to date/stratum/stratum_idx format ---
  if (all(c("Week", "Collapse") %in% names(strat_assign)) &&
      !("date" %in% names(strat_assign))) {

    all_dates <- sort(unique(as.Date(pass_dates)))
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

  # --- Extract psi posterior draws ---
  samp_mat  <- do.call(rbind, lapply(ge_fit$samples, as.matrix))
  psi_cols  <- grep("^psi\\[", colnames(samp_mat))
  psi_mat   <- samp_mat[, psi_cols, drop = FALSE]   # n_mcmc x n_obs_strata

  idx       <- sample(nrow(psi_mat), B, replace = nrow(psi_mat) < B)
  psi_draws <- psi_mat[idx, , drop = FALSE]          # B x n_obs_strata
  alpha_v   <- samp_mat[idx, grep("^alpha$", colnames(samp_mat))]
  beta_v    <- samp_mat[idx, grep("^beta$",  colnames(samp_mat))]

  # --- Build full psi matrix for ALL strata ---
  n_all    <- nrow(ge_data)
  psi_full <- matrix(NA_real_, nrow = B, ncol = n_all)

  obs_idx  <- match(obs_strata$stratum_idx, ge_data$stratum_idx)
  psi_full[, obs_idx] <- psi_draws

  # Zero-pool strata: use regression-predicted psi from alpha/beta draws
  zero_idx <- setdiff(seq_len(n_all), obs_idx)
  for (zi in zero_idx)
    psi_full[, zi] <- plogis(alpha_v + beta_v * ge_data$spill_val[zi])

  # --- Post-hoc GE: vectorised ---
  n_GRJ <- ge_data$n_GRJ_obs
  n_GRS <- ge_data$n_GRS_obs

  true_grs <- sweep(1 / psi_full, 2, n_GRS, "*")           # B x n_all
  denom    <- sweep(true_grs,     2, n_GRJ, "+")           # B x n_all
  ge_mat   <- matrix(n_GRJ, B, n_all, byrow = TRUE) / denom
  ge_mat[denom <= 0] <- NA_real_
  ge_mat   <- pmin(pmax(ge_mat, 0), 1)                     # clip to [0, 1]

  # Season-mean fallback: n_pool-weighted mean across strata
  pool_wts   <- ifelse(ge_data$n_pool > 0, ge_data$n_pool, 0)
  wt_mat     <- matrix(pool_wts, nrow = B, ncol = n_all, byrow = TRUE)
  wt_mat[is.na(ge_mat)] <- 0
  ge_nona    <- ge_mat; ge_nona[is.na(ge_nona)] <- 0
  wt_rowsums <- rowSums(wt_mat)
  season_ge  <- ifelse(wt_rowsums > 0,
                       rowSums(ge_nona * wt_mat) / wt_rowsums,
                       rowMeans(ge_mat, na.rm = TRUE))   # B-length vector

  # --- Map strata to daily passageData rows ---
  # Explicit Date coercion prevents silent NA from character format mismatches
  pass_dates_d  <- as.Date(pass_dates)
  strat_dates_d <- as.Date(strat_assign$date)
  day_match     <- match(pass_dates_d, strat_dates_d)
  strat_idx_v   <- strat_assign$stratum_idx[day_match]   # NA = unmatched date

  valid   <- !is.na(strat_idx_v)
  ge_days <- matrix(season_ge, nrow = B, ncol = length(pass_dates))  # fallback
  ge_days[, valid] <- ge_mat[, strat_idx_v[valid]]
  result  <- t(ge_days)   # n_days x B

  draws_df <- as.data.frame(result)
  colnames(draws_df) <- paste0("boot_", seq_len(B))
  cbind(data.frame(SampleEndDate = pass_dates), draws_df)
}
