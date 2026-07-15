#' @title Fit Bayesian multistate model for guidance efficiency at LGR
#'
#' @description Estimates guidance efficiency (psi) per trap stratum using a
#'   multistate mark-recapture model. Detection histories from GRJ (bypass),
#'   GRS (spillway bay 1), and GOJ (Little Goose bypass) are used to jointly
#'   estimate psi (GE), p (GRS detection probability), and phi (GOJ survival
#'   and bypass probability).
#'
#'   The logit of each parameter is modelled as a linear function of spill
#'   (LGR spill for psi and p; LGS spill for phi) plus stratum-level i.i.d.
#'   process error. A random walk extension on logit(psi) was evaluated but
#'   found to be non-identifiable given the nested hierarchical structure: the
#'   global intercept alpha and the walk's initial level delta[1] are collinear
#'   and both float freely in the posterior. The i.i.d. process error model
#'   with nested shrinkage (mu_strat) was retained.
#'
#'   The phi prior is derived from McCann et al. (2023) Table 8.3 LGS steelhead
#'   coefficients, reprojected onto the standardized LGS spill scale computed
#'   from this run's data (see Details). McCann et al. (2023) report only point
#'   estimates (medians across 1,000 resampling iterations) with no standard
#'   errors; prior SDs of 0.50 for both alpha_phi and beta_phi were chosen to
#'   allow meaningful departure from the McCann estimates while still anchoring
#'   phi in data-sparse strata.
#'
#'   \strong{Nested (two-level) shrinkage.} When \code{ge_data} carries a
#'   \code{parent} column (added by \code{\link{prep_ge_data}} via its
#'   \code{parent_strata} argument), each stratum's logit(psi) shrinks toward
#'   its parent group's mean (\code{mu_strat}), which in turn shrinks toward
#'   the global mean \code{alpha}. This part is unchanged and auto-detected,
#'   exactly as before. Set \code{nested_p = TRUE} to apply the same nested
#'   shrinkage to \code{p} via a parallel \code{mu_strat_p}, using the same
#'   \code{parent} grouping. This is a new, separate, opt-in argument (default
#'   \code{FALSE}) so existing callers relying on the original psi-only
#'   nesting see no change in behavior.
#'
#'   \strong{Route-specific phi.} Set \code{route_effect = TRUE} to let phi
#'   differ by LGR passage route. The bypass-route value \code{phi_B} follows
#'   the usual spill regression; the spillway-route value \code{phi_S} is
#'   offset from it by a per-stratum random effect,
#'   \code{delta_route[s] ~ N(delta_route0, sigma_route)}. A deterministic
#'   linear-in-spill form for \code{delta_route} was evaluated first and
#'   rejected: the slope was unstable across refits (driven almost entirely by
#'   a single high-leverage, atypically-low-spill stratum) even though a
#'   descriptive correlation between the random-effect deviations and spill
#'   was later found to be real and robust to that same point (Pearson r
#'   dropped from 0.93 to 0.82 excluding it, still p < 0.001; Spearman rho
#'   0.83). That correlation is a documented, unresolved lead for future
#'   refinement, not yet a reliable functional form -- the unstructured
#'   random-effect version is what's implemented here. Simulate-then-recover
#'   testing on MY2025 steelhead showed it matches the spill-parameterised
#'   version on phi_S recovery while being far more stable run to run, and
#'   clearly beats assuming phi_S = phi_B (coverage 0.74 vs 0.97, wins in all
#'   14 strata, not concentrated in the one high-leverage week). Default
#'   \code{FALSE}, reproducing the original shared-phi behavior.
#'
#' @param ge_data data frame returned by \code{\link{prep_ge_data}}. Must
#'   contain columns \code{stratum_idx}, \code{n_GRJ_obs}, \code{n_GRS_obs},
#'   \code{n_pool}, \code{spill_val}, and \code{lgs_spill_val}. Optional
#'   \code{parent} column (with more than one distinct value among
#'   observed strata) enables psi nesting (auto-detected) and is required if
#'   \code{nested_p = TRUE}.
#' @param nested_p logical. If \code{TRUE}, apply nested shrinkage to \code{p}
#'   (mirrors the existing psi nesting, same \code{parent} grouping, separate
#'   \code{mu_strat_p} hierarchy). Requires a usable \code{parent} column;
#'   errors otherwise. Default \code{FALSE}.
#' @param route_effect logical. If \code{TRUE}, split phi into \code{phi_B}
#'   (bypass route) and \code{phi_S} (spillway route), differing by a
#'   per-stratum random effect (see Details). Default \code{FALSE}, which
#'   reproduces the original single shared-\code{phi} behavior.
#' @param max_n integer. Maximum effective sample size per stratum for the
#'   multinomial likelihood. Default 1000. Sensitivity testing across
#'   max_n = 200, 500, 1000, and 2000 on MY2025 steelhead showed psi/p/phi
#'   point estimates stable across that range.
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
#'   \item{nested}{logical: whether psi used the nested (mu_strat) structure}
#'   \item{nested_p}{logical: whether p used the nested (mu_strat_p) structure}
#'   \item{route_effect}{logical: whether phi was split into phi_B/phi_S}
#'
#' @details The five detection history cell probabilities (conditioned on
#'   being observed at least once) are, when \code{route_effect = FALSE}:
#'   \enumerate{
#'     \item (GRJ, no GOJ): psi * (1 - phi)
#'     \item (GRJ, GOJ):    psi * phi
#'     \item (GRS, no GOJ): (1 - psi) * p * (1 - phi)
#'     \item (GRS, GOJ):    (1 - psi) * p * phi
#'     \item (no LGR, GOJ): (1 - psi) * (1 - p) * phi
#'   }
#'   and with \code{phi_B} in place of \code{phi} for cells 1-2 and
#'   \code{phi_S} in place of \code{phi} for cells 3-5 when
#'   \code{route_effect = TRUE}.
#'
#'   \strong{phi prior derivation.} McCann et al. (2023) Table 8.3 (LGS,
#'   steelhead): FGE = 0.81, beta_PropSpill = -0.07718, beta_Weir = -0.42932,
#'   PropSpill on 0-100 percent scale. stFlow dropped; Weir = 1.
#'   alpha_raw = logit(0.81) + (-0.42932)*1 = 1.0207; beta_raw = -0.07718.
#'   Reprojected onto standardized LGS spill each run:
#'   alpha_phi_mean = alpha_raw + beta_raw * mean(lgs_spill_val),
#'   beta_phi_mean  = beta_raw * sd(lgs_spill_val).
#'   Prior SDs of 0.50 for both alpha_phi and beta_phi reflect that McCann
#'   et al. (2023) report only point estimates with no standard errors. This
#'   prior anchors \code{phi_B} (or the shared \code{phi} when
#'   \code{route_effect = FALSE}) unchanged; \code{delta_route} is a
#'   separate, weakly-informative offset layered on top of it, not a
#'   reparameterisation of it.
#'
#' @export
fit_ge_model <- function(ge_data,
                         nested_p     = FALSE,
                         route_effect = FALSE,
                         max_n    = 1000,
                         n_iter   = 30000,
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

  # McCann et al. (2023) Table 8.3 LGS steelhead phi prior, reprojected.
  # Table 8.3 reports medians only (no SEs); prior SDs of 0.50 chosen to allow
  # meaningful departure from the McCann point estimates while anchoring phi.
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
    alpha_phi_sd     = 0.50,
    beta_phi_sd      = 0.50,
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

  if (nested_p && !nested)
    stop("nested_p = TRUE requires a usable 'parent' column (more than one ",
         "distinct value among observed strata), same requirement as psi nesting.")

  if (nested || nested_p) {
    parent_obs <- as.integer(factor(obs_strata$parent,
                                    levels = sort(unique(obs_strata$parent))))
    n_strat    <- max(parent_obs)
  }

  jags_data <- list(
    S             = S,
    S_lik         = S_lik,
    lik_idx       = as.array(lik_idx),
    n             = n_mat,
    N_seen        = N_seen,
    lgr_spill_std = as.numeric(lgr_spill_std),
    lgs_spill_std = as.numeric(lgs_spill_std),
    alpha_phi_mean = alpha_phi_prior_mean,
    alpha_phi_sd   = phi_prior$alpha_phi_sd,
    beta_phi_mean  = beta_phi_prior_mean,
    beta_phi_sd    = phi_prior$beta_phi_sd
  )
  if (nested || nested_p) {
    jags_data$parent  <- as.array(parent_obs)
    jags_data$n_strat <- n_strat
  }

  # ---- assemble the model string from conditional blocks --------------------
  # Two independent toggles (nested_p, route_effect; nested is auto-detected)
  # would mean several near-duplicate hardcoded model strings under the
  # original one-string-per-combination pattern. Assembling from small blocks
  # instead keeps each toggle's effect in exactly one place.

  psi_line <- if (nested)
    "      logit_psi[s] ~ dnorm(mu_strat[parent[s]] + beta * lgr_spill_std[s], tau_psi)\n"
  else
    "      logit_psi[s] ~ dnorm(alpha + beta * lgr_spill_std[s], tau_psi)\n"

  p_mean <- if (nested_p) "mu_strat_p[parent[s]]" else "alpha_p"
  p_line <- sprintf(
    "      logit_p[s] ~ dnorm(%s + beta_p * lgr_spill_std[s], tau_p)\n",
    p_mean
  )

  phi_block <- if (route_effect) "
      logit_phi[s] ~ dnorm(alpha_phi + beta_phi * lgs_spill_std[s], tau_phi)
      phi_B[s] <- ilogit(logit_phi[s])
      delta_route[s] ~ dnorm(delta_route0, tau_route)
      phi_S[s] <- ilogit(logit_phi[s] + delta_route[s])
      phi[s]   <- phi_B[s]
" else "
      logit_phi[s] ~ dnorm(alpha_phi + beta_phi * lgs_spill_std[s], tau_phi)
      phi[s] <- ilogit(logit_phi[s])
"

  pi_raw_block <- if (route_effect) "
      pi_raw[s, 1] <- psi[s] * (1 - phi_B[s])               + eps
      pi_raw[s, 2] <- psi[s] * phi_B[s]                     + eps
      pi_raw[s, 3] <- (1 - psi[s]) * p[s] * (1 - phi_S[s])  + eps
      pi_raw[s, 4] <- (1 - psi[s]) * p[s] * phi_S[s]        + eps
      pi_raw[s, 5] <- (1 - psi[s]) * (1 - p[s]) * phi_S[s]  + eps
" else "
      pi_raw[s, 1] <- psi[s] * (1 - phi[s])               + eps
      pi_raw[s, 2] <- psi[s] * phi[s]                     + eps
      pi_raw[s, 3] <- (1 - psi[s]) * p[s] * (1 - phi[s])  + eps
      pi_raw[s, 4] <- (1 - psi[s]) * p[s] * phi[s]        + eps
      pi_raw[s, 5] <- (1 - psi[s]) * (1 - p[s]) * phi[s]  + eps
"

  tau_extra <- paste0(
    if (nested)       "    tau_strat   <- pow(sigma_strat,   -2)\n" else "",
    if (nested_p)     "    tau_strat_p <- pow(sigma_strat_p, -2)\n" else "",
    if (route_effect) "    tau_route   <- pow(sigma_route,   -2)\n" else ""
  )

  hierarchy_block <- paste0(
    if (nested)   "    for (g in 1:n_strat) { mu_strat[g]   ~ dnorm(alpha,   tau_strat)   }\n" else "",
    if (nested_p) "    for (g in 1:n_strat) { mu_strat_p[g] ~ dnorm(alpha_p, tau_strat_p) }\n" else ""
  )

  priors_extra <- paste0(
    if (nested)       "    sigma_strat   ~ dunif(0.05, 3)\n" else "",
    if (nested_p)     "    sigma_strat_p ~ dunif(0.05, 3)\n" else "",
    if (route_effect) "    delta_route0 ~ dnorm(0, pow(0.5, -2))\n    sigma_route  ~ dunif(0.05, 3)\n" else ""
  )

  model_string <- paste0("
  model {

    eps       <- 1.0E-9
    tau_psi   <- pow(sigma_psi, -2)
    tau_p     <- pow(sigma_p,   -2)
    tau_phi   <- pow(sigma_phi, -2)
", tau_extra, "
    for (s in 1:S) {

", psi_line, "      psi[s] <- ilogit(logit_psi[s])

", p_line, "      p[s] <- ilogit(logit_p[s])
", phi_block, "
", pi_raw_block, "
      p_seen[s] <- pi_raw[s,1] + pi_raw[s,2] + pi_raw[s,3] +
                   pi_raw[s,4] + pi_raw[s,5]
      for (k in 1:5) { pi_obs[s, k] <- pi_raw[s, k] / p_seen[s] }
    }

    for (j in 1:S_lik) {
      n[lik_idx[j], 1:5] ~ dmulti(pi_obs[lik_idx[j], 1:5], N_seen[lik_idx[j]])
    }
", hierarchy_block, "
    # Weakly-informative Student-t priors following Hance et al. (2024):
    # intercepts ~ t(df=7, scale=2); slopes ~ t(df=7, scale=1)
    alpha     ~ dt(0, pow(2, -2), 7)
    beta      ~ dt(0, pow(1, -2), 7) T(, 0)
    sigma_psi ~ dunif(0.05, 3)

    alpha_p  ~ dt(0, pow(2, -2), 7)
    beta_p   ~ dt(0, pow(1, -2), 7)
    sigma_p  ~ dunif(0.05, 3)

    # phi prior: McCann et al. (2023) Table 8.3 LGS steelhead, reprojected.
    #   alpha_raw = logit(0.81) + (-0.42932)*1 = 1.0207
    #   beta_raw  = -0.07718  (per % spill)
    alpha_phi ~ dnorm(alpha_phi_mean, pow(alpha_phi_sd, -2))
    beta_phi  ~ dnorm(beta_phi_mean,  pow(beta_phi_sd,  -2)) T(, 0)
    sigma_phi ~ dunif(0.05, 3)
", priors_extra, "
  }")

  monitor <- c("psi", "p", "phi",
               "alpha", "beta", "sigma_psi",
               "alpha_p", "beta_p", "sigma_p",
               "alpha_phi", "beta_phi", "sigma_phi")
  if (nested)       monitor <- c(monitor, "mu_strat", "sigma_strat")
  if (nested_p)     monitor <- c(monitor, "mu_strat_p", "sigma_strat_p")
  if (route_effect) monitor <- c(monitor, "phi_B", "phi_S", "delta_route",
                                 "delta_route0", "sigma_route")

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
    nested         = nested,
    nested_p       = nested_p,
    route_effect   = route_effect
  )
}
