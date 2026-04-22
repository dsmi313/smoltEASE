#' @title Fit Bayesian hierarchical model for spillway detection probability
#'
#' @description Estimates spillway detection probability (psi) per trap stratum
#'   using a hierarchical JAGS model. The logit of psi is modelled as a linear
#'   function of spill (regression covariate) plus stratum-level random effects
#'   that share a common precision. Only strata with at least one psi-pool
#'   observation are included in the binomial likelihood; zero-data strata receive
#'   regression-predicted psi in \code{\link{generate_ge_draws}}.
#'
#' @param ge_data data frame returned by \code{\link{prep_ge_data}}.
#' @param n_iter number of MCMC iterations after burn-in per chain. Default 10000.
#' @param n_burnin number of adaptation/burn-in iterations. Default 3000.
#' @param n_chains number of MCMC chains. Default 3.
#' @param n_thin thinning interval. Default 2.
#' @param seed random seed passed to \code{set.seed()} for reproducibility. Default 42.
#'
#' @return A list with elements:
#'   \item{samples}{an \code{mcmc.list} object with posterior draws of
#'     \code{psi[1..k]}, \code{alpha}, \code{beta}, \code{sigma_psi}}
#'   \item{ge_data}{the input \code{ge_data} (passed through for downstream use)}
#'   \item{obs_strata}{the subset of \code{ge_data} with \code{n_pool > 0}}
#'
#' @details Requires JAGS to be installed on the system. Install with
#'   \code{sudo apt install jags} (Linux) or \code{brew install jags} (macOS),
#'   then \code{install.packages("rjags")} in R.
#'
#'   Check convergence with \code{coda::gelman.diag(ge_fit$samples)} —
#'   R-hat values should be below 1.1 for all parameters.
#'
#' @export
fit_ge_model <- function(ge_data,
                         n_iter   = 10000,
                         n_burnin = 3000,
                         n_chains = 3,
                         n_thin   = 2,
                         seed     = 42) {

  if (!requireNamespace("rjags", quietly = TRUE))
    stop("Package 'rjags' is required. Install with install.packages('rjags') ",
         "and ensure JAGS is installed on your system (apt/brew install jags).")

  obs_strata <- ge_data[ge_data$n_pool > 0, ]
  if (nrow(obs_strata) == 0)
    stop("No strata have psi-pool observations (n_pool > 0). ",
         "Check that dat_up and strat_assign are correctly aligned.")

  spill_mean <- mean(obs_strata$spill_val, na.rm = TRUE)
  spill_sd   <- sd(obs_strata$spill_val, na.rm = TRUE)

  if (!is.finite(spill_sd) || spill_sd == 0)
    stop("spill_val has zero or undefined SD among observed strata; cannot scale.")

  spill_sc <- (obs_strata$spill_val - spill_mean) / spill_sd

  jags_data <- list(
    n_strata   = nrow(obs_strata),
    n_GRS_pool = obs_strata$n_GRS_pool,
    n_pool     = obs_strata$n_pool,
    spill_sc   = spill_sc
  )

  model_string <- "
  model {
    for (s in 1:n_strata) {
      n_GRS_pool[s] ~ dbin(psi[s], n_pool[s])
      logit(psi[s]) <- logit_psi[s]
      logit_psi[s]  ~ dnorm(mu_psi[s], tau_psi)
      mu_psi[s]     <- alpha + beta * spill_sc[s]
    }

    alpha     ~ dnorm(0, 0.001)
    beta      ~ dnorm(0, 0.001)
    tau_psi   <- pow(sigma_psi, -2)
    sigma_psi ~ dnorm(0, 1) T(0,)
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
    variable.names = c("psi", "alpha", "beta", "sigma_psi"),
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
