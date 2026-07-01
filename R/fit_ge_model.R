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
#'   (LGR spill for psi; LGR spill + bay-1 spill for p; LGS spill for phi)
#'   plus stratum-level process error. When \code{ge_data} contains a
#'   \code{bay1_spill_val} column (added by \code{\link{prep_ge_data}} when
#'   \code{spill_data} carries \code{bay1.per}), \code{p} is modelled as a
#'   function of both total LGR spill and bay-1 spill percent, following the
#'   approach of Hance et al. (2024) who treat bay-1 operation as a direct
#'   modifier of GRS detection probability. When absent, only total LGR spill
#'   is used for \code{p}. The phi prior is derived from McCann et al. (2023)
#'   Table 8.3 LGS steelhead coefficients, reprojected onto the standardized
#'   LGS spill scale computed from this run's data (see Details).
#'
#'   GE is estimated directly as psi. \code{\link{generate_ge_draws}} then
#'   uses the psi posterior draws and raw GRJ/GRS counts to compute per-draw
#'   GE values for SCRAPI2.
#'
#'   Weekly counts are capped at \code{max_n} per stratum before fitting to
#'   prevent the multinomial likelihood from collapsing the posterior in
#'   data-rich strata.
#'
#'   \strong{Nested (two-level) shrinkage.} When \code{ge_data} carries a
#'   \code{parent} column (added by \code{\link{prep_ge_data}} via its
#'   \code{parent_strata} argument), the model becomes nested: each stratum's
#'   \code{logit(psi)} shrinks toward its \emph{parent} group's mean
#'   (\code{mu_strat}), and the parent means shrink toward the global mean
#'   \code{alpha}. Without a \code{parent} column the original single-level
#'   model is used unchanged.
#'
#' @param ge_data data frame returned by \code{\link{prep_ge_data}}. Must
#'   contain columns \code{stratum_idx}, \code{n_GRJ_obs}, \code{n_GRS_obs},
#'   \code{n_pool}, \code{spill_val}, and \code{lgs_spill_val}. Optional
#'   \code{bay1_spill_val} column triggers the two-covariate p regression.
#'   Optional \code{parent} column triggers the nested model.
#' @param max_n integer. Maximum effective sample size per stratum for the
#'   multinomial likelihood. Counts are scaled proportionally if N_seen > max_n.
#'   Default 1000. Sensitivity testing across max_n = 200, 500, 1000, and 2000
#'   on MY2025 steelhead showed psi/p/phi point estimates stable to within
#'   typical posterior uncertainty across that range, with CI width decreasing
#'   smoothly; 1000 was chosen as a value that retains meaningful CI precision
#'   without showing signs of overconfidence.
#' @param n_iter number of MCMC iterations after burn-in per chain. Default 10000.
#' @param n_burnin number of adaptation/burn-in iterations. Default 5000.
#' @param n_chains number of MCMC chains. Default 3.
#' @param n_thin thinning interval. Default 5.
#' @param seed random seed. Default 42.
#'
#' @return A list with elements:
#'   \item{samples}{an \code{mcmc.list} of posterior draws}
#'   \item{ge_data}{the input \code{ge_data} (passed through)}
#'   \item{obs_strata}{subset of \code{ge_data} with \code{n_pool > 0}}
#'   \item{spill_mean}{mean of \code{spill_val} used for standardisation}
#'   \item{spill_sd}{SD of \code{spill_val} used for standardisation}
#'   \item{lgs_spill_mean}{mean of \code{lgs_spill_val} used for standardisation}
#'   \item{lgs_spill_sd}{SD of \code{lgs_spill_val} used for standardisation}
#'   \item{phi_prior}{list with McCann-derived phi prior means and raw coefficients}
#'   \item{has_bay1}{logical: whether the two-covariate p regression was used}
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
#'   \strong{Bay-1 spill covariate for p.} GRS is the bay-1 spillway antenna.
#'   The probability that a spillway fish is detected at GRS depends on both
#'   total spill magnitude (which determines how many fish take the spillway)
#'   and the fraction of spillway flow through bay 1 (which determines how
#'   many of those fish are attracted to the monitored bay). Following
#'   Hance et al. (2024), when \code{bay1_spill_val} is present in
#'   \code{ge_data} the p regression includes both covariates:
#'   \code{logit(p) = alpha_p + beta_p * lgr_spill_std + beta2_p * bay1_spill_std}.
#'   Both slopes are assigned the same weakly-informative Student-t prior as
#'   the other regression coefficients (Hance et al. 2024).
#'
#'   \strong{phi prior derivation.} McCann et al. (2023) Table 8.3 (LGS,
#'   steelhead): FGE = 0.81, beta_PropSpill = -0.07718, beta_Weir = -0.42932,
#'   PropSpill on 0-100 percent scale (confirmed by worked example). stFlow
#'   dropped (held at reference); Weir = 1. Raw-scale intercept and slope
#'   reprojected onto standardized LGS spill each run:
#'   \code{alpha_phi_mean = alpha_raw + beta_raw * mean(lgs_spill_val)},
#'   \code{beta_phi_mean = beta_raw * sd(lgs_spill_val)}.
#'
#' @export
fit_ge_model <- function(ge_data,
                         max_n    = 1000,
                         n_iter   = 10000,
                         n_burnin = 5000,
                         n_chains = 3,
                         n_thin   = 5,
                         seed     = 42) {

  if (!requireNamespace("rjags", quietly = TRUE))
    stop("Package 'rjags' is required. Install with install.packages('rjags') ",
         "and ensure JAGS is installed on your system.")

  required_cols <- c("stratum_idx", "n_GRJ_obs", "n_GRS_obs", "n_pool",
                     "spill_val", "lgs_spill_val",
                     "n_h1", "n_h2", "n_h3", "n_h4", "n_h5")
  missing <- setdiff(required_cols, names(ge_data))
  if (length(missing) > 0)
    stop("ge_data is missing required columns: ", paste(missing, collapse = ", "))

  has_bay1 <- "bay1_spill_val" %in% names(ge_data)
  if (has_bay1) {
    message("bay1_spill_val found: using two-covariate p regression ",
            "(total LGR spill + bay-1 spill percent), following Hance et al. (2024).")
  } else {
    message("bay1_spill_val not found: using single-covariate p regression ",
            "(total LGR spill only). Add bay1.per to spill_data to enable the ",
            "two-covariate p regression.")
  }

  obs_strata <- ge_data[ge_data$n_pool > 0, ]
  if (nrow(obs_strata) == 0)
    stop("No strata have psi-pool observations (n_pool > 0).")

  S <- nrow(obs_strata)

  # Standardise LGR spill
  spill_mean <- mean(obs_strata$spill_val, na.rm = TRUE)
  spill_sd   <- sd(obs_strata$spill_val, na.rm = TRUE)
  if (!is.finite(spill_sd) || spill_sd == 0)
    stop("spill_val has zero or undefined SD among observed strata.")
  lgr_spill_std <- (obs_strata$spill_val - spill_mean) / spill_sd

  # Standardise LGS spill
  lgs_spill_mean <- mean(obs_strata$lgs_spill_val, na.rm = TRUE)
  lgs_spill_sd   <- sd(obs_strata$lgs_spill_val,   na.rm = TRUE)
  if (!is.finite(lgs_spill_sd) || lgs_spill_sd == 0)
    stop("lgs_spill_val has zero or undefined SD among observed strata.")
  lgs_spill_std <- (obs_strata$lgs_spill_val - lgs_spill_mean) / lgs_spill_sd

  # Standardise bay-1 spill (when present)
  if (has_bay1) {
    bay1_spill_mean <- mean(obs_strata$bay1_spill_val, na.rm = TRUE)
    bay1_spill_sd   <- sd(obs_strata$bay1_spill_val,   na.rm = TRUE)
    if (!is.finite(bay1_spill_sd) || bay1_spill_sd == 0)
      stop("bay1_spill_val has zero or undefined SD among observed strata.")
    bay1_spill_std <- (obs_strata$bay1_spill_val - bay1_spill_mean) / bay1_spill_sd
  }

  # McCann et al. (2023) Table 8.3 LGS steelhead phi prior, reprojected
  #   alpha_raw = logit(0.81) + (-0.42932)*1 = 1.0207
  #   beta_raw  = -0.07718  (per % spill)
  #   alpha_phi_mean = alpha_raw + beta_raw * mean(lgs_spill_val)
  #   beta_phi_mean  = beta_raw  * sd(lgs_spill_val)
  mccann_FGE_LGS_sthd       <- 0.81
  mccann_weir_LGS_sthd      <- -0.42932
  mccann_propspill_LGS_sthd <- -0.07718
  mccann_weir_on            <- 1

  alpha_phi_raw        <- qlogis(mccann_FGE_LGS_sthd) + mccann_weir_LGS_sthd * mccann_weir_on
  beta_phi_raw         <- mccann_propspill_LGS_sthd
  alpha_phi_prior_mean <- alpha_phi_raw + beta_phi_raw * lgs_spill_mean
  beta_phi_prior_mean  <- beta_phi_raw  * lgs_spill_sd

  phi_prior <- list(
    alpha_phi_mean   = alpha_phi_prior_mean,
    beta_phi_mean    = beta_phi_prior_mean,
    alpha_phi_sd     = 0.30,
    beta_phi_sd      = 0.25,
    alpha_raw        = alpha_phi_raw,
    beta_raw         = beta_phi_raw,
    mccann_FGE       = mccann_FGE_LGS_sthd,
    mccann_weir      = mccann_weir_LGS_sthd,
    mccann_propspill = mccann_propspill_LGS_sthd,
    weir_on          = mccann_weir_on
  )

  message(sprintf(
    "phi prior (McCann 2023 Table 8.3, LGS steelhead, reprojected): alpha_phi ~ N(%.3f, %.2f), beta_phi ~ N(%.3f, %.2f)",
    alpha_phi_prior_mean, phi_prior$alpha_phi_sd, beta_phi_prior_mean, phi_prior$beta_phi_sd
  ))

  # Build history count matrix and cap N_seen
  n_mat  <- as.matrix(obs_strata[, c("n_h1","n_h2","n_h3","n_h4","n_h5")])
  N_seen <- rowSums(n_mat)
  scale  <- ifelse(N_seen > 0, pmin(1, max_n / N_seen), 1)
  n_mat  <- round(sweep(n_mat, 1, scale, "*"))
  N_seen <- rowSums(n_mat)

  lik_idx <- which(N_seen > 0)
  S_lik   <- length(lik_idx)
  if (S_lik == 0)
    stop("No strata have detection-history data (all N_seen == 0).")

  nested <- "parent" %in% names(obs_strata) &&
            length(unique(obs_strata$parent)) > 1
  if (nested) {
    parent_obs <- as.integer(factor(obs_strata$parent,
                                    levels = sort(unique(obs_strata$parent))))
    n_strat    <- max(parent_obs)
  }

  jags_data <- list(
    S              = S,
    S_lik          = S_lik,
    lik_idx        = as.array(lik_idx),
    n              = n_mat,
    N_seen         = N_seen,
    lgr_spill_std  = as.numeric(lgr_spill_std),
    lgs_spill_std  = as.numeric(lgs_spill_std),
    alpha_phi_mean = alpha_phi_prior_mean,
    alpha_phi_sd   = phi_prior$alpha_phi_sd,
    beta_phi_mean  = beta_phi_prior_mean,
    beta_phi_sd    = phi_prior$beta_phi_sd
  )
  if (has_bay1)  jags_data$bay1_spill_std <- as.numeric(bay1_spill_std)
  if (nested)  { jags_data$parent  <- as.array(parent_obs)
                 jags_data$n_strat <- n_strat }

  # p regression block (shared across flat/nested): single or two-covariate
  model_flat_single <- "
  model {

    eps       <- 1.0E-9
    tau_psi   <- pow(sigma_psi, -2)
    tau_p     <- pow(sigma_p,   -2)
    tau_phi   <- pow(sigma_phi, -2)

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
      for (k in 1:5) { pi_obs[s, k] <- pi_raw[s, k] / p_seen[s] }
    }

    for (j in 1:S_lik) {
      n[lik_idx[j], 1:5] ~ dmulti(pi_obs[lik_idx[j], 1:5], N_seen[lik_idx[j]])
    }

    # Weakly-informative Student-t priors following Hance et al. (2024):
    # intercepts ~ t(df=7, scale=2); slopes ~ t(df=7, scale=1)
    alpha     ~ dt(0, pow(2, -2), 7)
    beta      ~ dt(0, pow(1, -2), 7) T(, 0)
    sigma_psi ~ dunif(0.05, 3)

    alpha_p  ~ dt(0, pow(2, -2), 7)
    beta_p   ~ dt(0, pow(1, -2), 7)
    sigma_p  ~ dunif(0.05, 3)

    alpha_phi ~ dnorm(alpha_phi_mean, pow(alpha_phi_sd, -2))
    beta_phi  ~ dnorm(beta_phi_mean,  pow(beta_phi_sd,  -2)) T(, 0)
    sigma_phi ~ dunif(0.05, 1)
  }"

  model_flat_bay1 <- "
  model {

    eps       <- 1.0E-9
    tau_psi   <- pow(sigma_psi, -2)
    tau_p     <- pow(sigma_p,   -2)
    tau_phi   <- pow(sigma_phi, -2)

    for (s in 1:S) {

      logit_psi[s] ~ dnorm(alpha + beta * lgr_spill_std[s], tau_psi)
      psi[s] <- ilogit(logit_psi[s])

      logit_p[s] ~ dnorm(alpha_p + beta_p * lgr_spill_std[s] + beta2_p * bay1_spill_std[s], tau_p)
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
      for (k in 1:5) { pi_obs[s, k] <- pi_raw[s, k] / p_seen[s] }
    }

    for (j in 1:S_lik) {
      n[lik_idx[j], 1:5] ~ dmulti(pi_obs[lik_idx[j], 1:5], N_seen[lik_idx[j]])
    }

    # Weakly-informative Student-t priors following Hance et al. (2024):
    # intercepts ~ t(df=7, scale=2); slopes ~ t(df=7, scale=1)
    alpha     ~ dt(0, pow(2, -2), 7)
    beta      ~ dt(0, pow(1, -2), 7) T(, 0)
    sigma_psi ~ dunif(0.05, 3)

    alpha_p  ~ dt(0, pow(2, -2), 7)
    beta_p   ~ dt(0, pow(1, -2), 7)
    beta2_p  ~ dt(0, pow(1, -2), 7)
    sigma_p  ~ dunif(0.05, 3)

    alpha_phi ~ dnorm(alpha_phi_mean, pow(alpha_phi_sd, -2))
    beta_phi  ~ dnorm(beta_phi_mean,  pow(beta_phi_sd,  -2)) T(, 0)
    sigma_phi ~ dunif(0.05, 1)
  }"

  model_nested_single <- "
  model {

    eps       <- 1.0E-9
    tau_psi   <- pow(sigma_psi,   -2)
    tau_strat <- pow(sigma_strat, -2)
    tau_p     <- pow(sigma_p,     -2)
    tau_phi   <- pow(sigma_phi,   -2)

    for (s in 1:S) {

      logit_psi[s] ~ dnorm(mu_strat[parent[s]] + beta * lgr_spill_std[s], tau_psi)
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
      for (k in 1:5) { pi_obs[s, k] <- pi_raw[s, k] / p_seen[s] }
    }

    for (j in 1:S_lik) {
      n[lik_idx[j], 1:5] ~ dmulti(pi_obs[lik_idx[j], 1:5], N_seen[lik_idx[j]])
    }

    for (g in 1:n_strat) { mu_strat[g] ~ dnorm(alpha, tau_strat) }

    # Weakly-informative Student-t priors following Hance et al. (2024):
    # intercepts ~ t(df=7, scale=2); slopes ~ t(df=7, scale=1)
    alpha       ~ dt(0, pow(2, -2), 7)
    beta        ~ dt(0, pow(1, -2), 7) T(, 0)
    sigma_psi   ~ dunif(0.05, 3)
    sigma_strat ~ dunif(0.05, 3)

    alpha_p  ~ dt(0, pow(2, -2), 7)
    beta_p   ~ dt(0, pow(1, -2), 7)
    sigma_p  ~ dunif(0.05, 3)

    alpha_phi ~ dnorm(alpha_phi_mean, pow(alpha_phi_sd, -2))
    beta_phi  ~ dnorm(beta_phi_mean,  pow(beta_phi_sd,  -2)) T(, 0)
    sigma_phi ~ dunif(0.05, 1)
  }"

  model_nested_bay1 <- "
  model {

    eps       <- 1.0E-9
    tau_psi   <- pow(sigma_psi,   -2)
    tau_strat <- pow(sigma_strat, -2)
    tau_p     <- pow(sigma_p,     -2)
    tau_phi   <- pow(sigma_phi,   -2)

    for (s in 1:S) {

      logit_psi[s] ~ dnorm(mu_strat[parent[s]] + beta * lgr_spill_std[s], tau_psi)
      psi[s] <- ilogit(logit_psi[s])

      logit_p[s] ~ dnorm(alpha_p + beta_p * lgr_spill_std[s] + beta2_p * bay1_spill_std[s], tau_p)
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
      for (k in 1:5) { pi_obs[s, k] <- pi_raw[s, k] / p_seen[s] }
    }

    for (j in 1:S_lik) {
      n[lik_idx[j], 1:5] ~ dmulti(pi_obs[lik_idx[j], 1:5], N_seen[lik_idx[j]])
    }

    for (g in 1:n_strat) { mu_strat[g] ~ dnorm(alpha, tau_strat) }

    # Weakly-informative Student-t priors following Hance et al. (2024):
    # intercepts ~ t(df=7, scale=2); slopes ~ t(df=7, scale=1)
    alpha       ~ dt(0, pow(2, -2), 7)
    beta        ~ dt(0, pow(1, -2), 7) T(, 0)
    sigma_psi   ~ dunif(0.05, 3)
    sigma_strat ~ dunif(0.05, 3)

    alpha_p  ~ dt(0, pow(2, -2), 7)
    beta_p   ~ dt(0, pow(1, -2), 7)
    beta2_p  ~ dt(0, pow(1, -2), 7)
    sigma_p  ~ dunif(0.05, 3)

    alpha_phi ~ dnorm(alpha_phi_mean, pow(alpha_phi_sd, -2))
    beta_phi  ~ dnorm(beta_phi_mean,  pow(beta_phi_sd,  -2)) T(, 0)
    sigma_phi ~ dunif(0.05, 1)
  }"

  model_string <- if (nested && has_bay1)  model_nested_bay1   else
                  if (nested)              model_nested_single  else
                  if (has_bay1)            model_flat_bay1      else
                                           model_flat_single

  monitor <- c("psi", "p", "phi",
               "alpha", "beta", "sigma_psi",
               "alpha_p", "beta_p", "sigma_p",
               "alpha_phi", "beta_phi", "sigma_phi")
  if (has_bay1) monitor <- c(monitor, "beta2_p")
  if (nested)   monitor <- c(monitor, "mu_strat", "sigma_strat")

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
    variable.names = monitor,
    n.iter = n_iter,
    thin   = n_thin
  )

  list(
    samples        = samples,
    ge_data        = ge_data,
    obs_strata     = obs_strata,
    spill_mean     = spill_mean,
    spill_sd       = spill_sd,
    lgs_spill_mean = lgs_spill_mean,
    lgs_spill_sd   = lgs_spill_sd,
    phi_prior      = phi_prior,
    has_bay1       = has_bay1
  )
}
