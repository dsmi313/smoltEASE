#' @title Fit Bayesian multistate model for guidance efficiency at LGR
#'
#' @description Estimates guidance efficiency (psi) per trap stratum using a
#'   multistate mark-recapture model. Detection histories from GRJ (bypass),
#'   GRS (spillway bay 1), and GOJ (Little Goose bypass) are used to jointly
#'   estimate psi (GE), p (GRS detection probability), and phi (GOJ detection
#'   probability). This replaces the original binomial psi-pool likelihood with
#'   a more complete observation model that accounts for fish missed by both
#'   LGR arrays but detected downstream.
#'
#'   The logit of each parameter is modelled as a linear function of spill
#'   (LGR spill for psi and p; LGS spill for phi) plus stratum-level process
#'   error. Informative priors on phi are derived from McCann et al. (2023)
#'   Table 8.3 LGS steelhead coefficients.
#'
#'   GE is estimated directly as psi. \code{\link{generate_ge_draws}} then
#'   uses the psi posterior draws and raw GRJ/GRS counts to compute per-draw
#'   GE values for SCRAPI2.
#'
#'   Weekly counts are capped at \code{max_n} per stratum before fitting to
#'   prevent the multinomial likelihood from collapsing the posterior in
#'   data-rich strata (following data-thinning approach of McCann et al. 2023).
#'
#'   The model is hierarchical: every stratum's \code{logit(psi)} is drawn from
#'   the shared spill regression plus process error, so strata may be defined as
#'   finely as individual weeks. Weeks with sparse PIT data are partially pooled
#'   toward the regression, and weeks with no detection-history data contribute
#'   no likelihood term but still receive a regression-based psi (their draws are
#'   completed in \code{\link{generate_ge_draws}}). To estimate GE weekly, pass a
#'   \code{strat_assign} to \code{\link{prep_ge_data}} in which \code{Collapse}
#'   equals \code{Week} (one stratum per week).
#'
#' @param ge_data data frame returned by \code{\link{prep_ge_data}}. Must
#'   contain columns \code{stratum_idx}, \code{n_GRJ_obs}, \code{n_GRS_obs},
#'   \code{n_pool}, \code{spill_val}, and \code{lgs_spill_val}.
#' @param max_n integer. Maximum effective sample size per stratum for the
#'   multinomial likelihood. Counts are scaled proportionally if N_seen > max_n.
#'   Default 300.
#' @param n_iter number of MCMC iterations after burn-in per chain. Default 10000.
#' @param n_burnin number of adaptation/burn-in iterations. Default 5000.
#' @param n_chains number of MCMC chains. Default 3.
#' @param n_thin thinning interval. Default 5.
#' @param seed random seed. Default 42.
#'
#' @return A list with elements:
#'   \item{samples}{an \code{mcmc.list} with posterior draws of
#'     \code{psi[1..S]}, \code{p[1..S]}, \code{phi[1..S]},
#'     \code{alpha} (= alpha_psi), \code{beta} (= beta_psi),
#'     \code{sigma_psi}, \code{alpha_p}, \code{beta_p}, \code{sigma_p},
#'     \code{alpha_phi}, \code{beta_phi}, \code{sigma_phi}}
#'   \item{ge_data}{the input \code{ge_data} (passed through)}
#'   \item{obs_strata}{subset of \code{ge_data} with \code{n_pool > 0}}
#'   \item{spill_mean}{mean of \code{spill_val} used for standardisation}
#'   \item{spill_sd}{SD of \code{spill_val} used for standardisation}
#'
#' @details The five detection history cell probabilities (conditioned on
#'   being observed at least once) are:
#'   \enumerate{
#'     \item (GRJ, no GOJ): psi * (1 - phi)
#'     \item (GRJ, GOJ):    psi * phi
#'     \item (GRS, no GOJ): (1 - psi) * p * (1 - phi)
#'     \item (GRS, GOJ):    (1 - psi) * p * phi
#'     \item (no LGR, GOJ): (1 - psi) * (1 - p) * phi
#'   }
#'
#' @export
fit_ge_model <- function(ge_data,
                         max_n    = 300,
                         n_iter   = 10000,
                         n_burnin = 5000,
                         n_chains = 3,
                         n_thin   = 5,
                         seed     = 42) {

  if (!requireNamespace("rjags", quietly = TRUE))
    stop("Package 'rjags' is required. Install with install.packages('rjags') ",
         "and ensure JAGS is installed on your system (apt/brew install jags).")

  required_cols <- c("stratum_idx", "n_GRJ_obs", "n_GRS_obs", "n_pool",
                     "spill_val", "lgs_spill_val",
                     "n_h1", "n_h2", "n_h3", "n_h4", "n_h5")
  missing <- setdiff(required_cols, names(ge_data))
  if (length(missing) > 0)
    stop("ge_data is missing required columns: ", paste(missing, collapse = ", "),
         "\nEnsure prep_ge_data() was called with GOJ detection data.")

  obs_strata <- ge_data[ge_data$n_pool > 0, ]
  if (nrow(obs_strata) == 0)
    stop("No strata have psi-pool observations (n_pool > 0).")

  S <- nrow(obs_strata)

  # Standardise LGR spill on observed strata
  spill_mean <- mean(obs_strata$spill_val, na.rm = TRUE)
  spill_sd   <- sd(obs_strata$spill_val, na.rm = TRUE)
  if (!is.finite(spill_sd) || spill_sd == 0)
    stop("spill_val has zero or undefined SD among observed strata.")

  lgr_spill_std <- (obs_strata$spill_val     - spill_mean) / spill_sd
  lgs_spill_std <- (obs_strata$lgs_spill_val - mean(obs_strata$lgs_spill_val, na.rm = TRUE)) /
                    sd(obs_strata$lgs_spill_val, na.rm = TRUE)

  # Build history count matrix and cap N_seen
  n_mat <- as.matrix(obs_strata[, c("n_h1","n_h2","n_h3","n_h4","n_h5")])
  N_seen <- rowSums(n_mat)
  scale  <- ifelse(N_seen > 0, pmin(1, max_n / N_seen), 1)
  n_mat  <- round(sweep(n_mat, 1, scale, "*"))
  N_seen <- rowSums(n_mat)

  # Strata with detection-history data enter the multinomial likelihood. Strata
  # with none (possible when estimating weekly with sparse PIT data) still get a
  # hierarchical psi from the spill regression + process error, but contribute
  # no likelihood term -- this keeps low-data weeks estimable via partial
  # pooling rather than dropping or breaking them.
  lik_idx <- which(N_seen > 0)
  S_lik   <- length(lik_idx)
  if (S_lik == 0)
    stop("No strata have detection-history data (all N_seen == 0).")

  jags_data <- list(
    S             = S,
    S_lik         = S_lik,
    lik_idx       = as.array(lik_idx),
    n             = n_mat,
    N_seen        = N_seen,
    lgr_spill_std = as.numeric(lgr_spill_std),
    lgs_spill_std = as.numeric(lgs_spill_std)
  )

  model_string <- "
  model {

    eps     <- 1.0E-9
    tau_psi <- pow(sigma_psi, -2)
    tau_p   <- pow(sigma_p,   -2)
    tau_phi <- pow(sigma_phi, -2)

    for (s in 1:S) {

      logit_psi[s] ~ dnorm(alpha + beta * lgr_spill_std[s], tau_psi)
      psi[s] <- ilogit(logit_psi[s])

      logit_p[s] ~ dnorm(alpha_p + beta_p * lgr_spill_std[s], tau_p)
      p[s] <- ilogit(logit_p[s])

      logit_phi[s] ~ dnorm(alpha_phi + beta_phi * lgs_spill_std[s], tau_phi)
      phi[s] <- ilogit(logit_phi[s])

      pi_raw[s, 1] <- psi[s] * (1 - phi[s])               + eps
      pi_raw[s, 2] <- psi[s] * phi[s]                     + eps
      pi_raw[s, 3] <- (1 - psi[s]) * p[s] * (1 - phi[s])  + eps
      pi_raw[s, 4] <- (1 - psi[s]) * p[s] * phi[s]        + eps
      pi_raw[s, 5] <- (1 - psi[s]) * (1 - p[s]) * phi[s]  + eps

      p_seen[s] <- pi_raw[s,1] + pi_raw[s,2] + pi_raw[s,3] +
                   pi_raw[s,4] + pi_raw[s,5]

      for (k in 1:5) {
        pi_obs[s, k] <- pi_raw[s, k] / p_seen[s]
      }
    }

    # multinomial likelihood only for strata with detection-history data
    for (j in 1:S_lik) {
      n[lik_idx[j], 1:5] ~ dmulti(pi_obs[lik_idx[j], 1:5], N_seen[lik_idx[j]])
    }

    # psi (GE): named alpha/beta to match generate_ge_draws expectations
    alpha     ~ dt(0, pow(2, -2), 7)
    beta      ~ dt(0, pow(1, -2), 7) T(, 0)
    sigma_psi ~ dunif(0.05, 3)

    alpha_p  ~ dt(0, pow(2, -2), 7)
    beta_p   ~ dt(0, pow(1, -2), 7)
    sigma_p  ~ dunif(0.05, 3)

    # Informative priors from McCann et al. (2023) Table 8.3 LGS steelhead
    alpha_phi ~ dnorm(-2.60, pow(0.30, -2))
    beta_phi  ~ dnorm(-1.78, pow(0.25, -2)) T(, 0)
    sigma_phi ~ dunif(0.05, 1)
  }"

  set.seed(seed)
  jags_fit <- rjags::jags.model(
    textConnection(model_string),
    data     = jags_data,
    n.chains = n_chains,
    n.adapt  = n_burnin,
    quiet    = TRUE
  )

  samples <- rjags::coda.samples(
    jags_fit,
    variable.names = c("psi", "p", "phi",
                       "alpha", "beta", "sigma_psi",
                       "alpha_p", "beta_p", "sigma_p",
                       "alpha_phi", "beta_phi", "sigma_phi"),
    n.iter = n_iter,
    thin   = n_thin
  )

  list(
    samples    = samples,
    ge_data    = ge_data,
    obs_strata = obs_strata,
    spill_mean = spill_mean,
    spill_sd   = spill_sd
  )
}
