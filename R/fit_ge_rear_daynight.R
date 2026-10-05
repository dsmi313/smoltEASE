#' Fit the rear-type six-cell GE model with a day/night split
#'
#' Extends [fit_ge_rear_model2()] (full rear structure, rear-flexible phi,
#' free route offset) by splitting each rear-week into fish that passed LGR by
#' day and by night. GE has a day value and a night value:
#'   logit psi_day[r,s]   = weekly GE regression (as in fit_ge_rear_model2)
#'   logit psi_night[r,s] = logit psi_day[r,s] + night_eff[r,s]
#'   night_eff[r,s] ~ Normal(night_mu[r], sigma_night)
#' a[r,s] is the share of that rear's fish passing LGR at night. The weekly GE
#' for all fish is psi[r,s] = (1 - a) psi_day + a psi_night.
#'
#' Downstream recovery, the route offset, and transport are the same for day
#' and night fish within a rear-week, and so is GRS detection (p) unless
#' `p_night_offset` is set. Because p is common by default,
#' GOJ-only fish (cell 5), which have no LGR time, enter as one pooled cell
#' whose probability sums over day and night. This common-p assumption cannot
#' be tested with these data; check sensitivity if it matters.
#'
#' Likelihood per rear-week: multinomial over 11 cells, day c1 c2 c3 c4 c6,
#' night c1 c2 c3 c4 c6, and pooled c5, conditional on being seen.
#'
#' @param ge_dn Named list (one per rear) of [prep_ge_data_daynight()] results.
#' @param night_sd Prior SD of each rear's mean night effect (logit scale).
#' @param night_week_sd_scale Half-normal prior scale for the weekly SD of the
#'   night effect around its rear mean.
#' @param day_ge "rear" (default): each rear-week has its own day GE residual.
#'   "shared": the weekly day GE residual is shared by all rears, so rears
#'   differ only by a constant offset (plus their own covariate values):
#'   logit psi_day[r,s] = eta[s] + rear_psi[r] + covariates[r,s],
#'   eta[s] ~ Normal(mu_strat[parent[s]], sigma_psi). Use when one rear has
#'   too few daytime guided fish to estimate its own weekly day GE.
#' @param p_night_offset Fixed logit-scale offset for GRS detection of
#'   night-passing spilled fish: logit p_night = logit p + p_night_offset.
#'   0 (default) is the common-p model; nonzero values are for sensitivity.
#' @inheritParams fit_ge_rear_model2
#' @return A list like [fit_ge_rear_model2()], with `psi_summary` for the
#'   weekly mixture, `period_summary` (psi_day, psi_night, a), and the fields
#'   [generate_rear_ge_draws_daynight()] needs.
#' @export
fit_ge_rear_daynight <- function(
    ge_dn, weeks = NULL, target_rear = "W", parent = NULL,
    psi_spill = NULL, psi_outflow = NULL,
    delta_sd = 0.5, rear_sd_scale = 1, phi_slope_sd_scale = 0.5,
    night_sd = 2.5, night_week_sd_scale = 1,
    day_ge = c("rear", "shared"), p_night_offset = 0,
    trans_on = NULL, phi_prior = NULL,
    n_iter = 30000, n_adapt = 2000, n_burnin = 10000, n_chains = 3,
    n_thin = 10, seed = 11, rhat_threshold = 1.01, verbose = TRUE) {

  if (!is.list(ge_dn) || is.null(names(ge_dn)) || length(ge_dn) < 2L ||
      anyDuplicated(names(ge_dn))) {
    stop("ge_dn must be a uniquely named list with at least two rear groups.", call. = FALSE)
  }
  rear_levels <- names(ge_dn)
  if (!target_rear %in% rear_levels) stop("target_rear must name one element of ge_dn.", call. = FALSE)
  for (nm in rear_levels) {
    if (!all(c("base", "n_day", "n_night", "n_c5") %in% names(ge_dn[[nm]]))) {
      stop("ge_dn[['", nm, "']] is not a prep_ge_data_daynight() result.", call. = FALSE)
    }
  }
  day_ge <- match.arg(day_ge)
  if (!is.numeric(p_night_offset) || length(p_night_offset) != 1L || !is.finite(p_night_offset)) {
    stop("p_night_offset must be one finite number.", call. = FALSE)
  }
  if (n_iter %% n_thin != 0L) stop("n_iter must be divisible by n_thin.", call. = FALSE)
  target <- ge_dn[[target_rear]]$base
  S <- nrow(as.matrix(target$n))
  R <- length(rear_levels)
  if (is.null(weeks)) weeks <- target$weeks
  weeks <- as.integer(weeks)
  if (length(weeks) != S) stop("weeks must have one value per six-cell row.", call. = FALSE)

  n <- array(0L, dim = c(R, S, 11L),
             dimnames = list(rear_levels, paste0("wk", weeks),
                             c(paste0("day_c", c(1:4, 6)), paste0("night_c", c(1:4, 6)), "c5")))
  for (r in seq_len(R)) {
    x <- ge_dn[[r]]
    if (!isTRUE(all.equal(as.numeric(x$base$lgr_spill_pct), as.numeric(target$lgr_spill_pct))) ||
        !isTRUE(all.equal(as.numeric(x$base$lgs_spill_pct), as.numeric(target$lgs_spill_pct)))) {
      stop("All rear groups must use identical LGR and LGS spill values.", call. = FALSE)
    }
    n[r, , 1:5] <- as.matrix(x$n_day)[, c(1:4, 6)]
    n[r, , 6:10] <- as.matrix(x$n_night)[, c(1:4, 6)]
    n[r, , 11] <- x$n_c5
  }
  storage.mode(n) <- "integer"
  N_seen <- apply(n, c(1L, 2L), sum)
  lik <- which(N_seen > 0, arr.ind = TRUE)

  if (!is.null(trans_on)) {
    if (!is.matrix(trans_on) || !identical(dim(trans_on), c(R, S)) ||
        anyNA(trans_on) || !all(trans_on %in% 0:1)) {
      stop("trans_on must be an R x S matrix of 0/1.", call. = FALSE)
    }
    trans_on <- matrix(as.integer(trans_on), R, S)
    if (any((n[, , 5] + n[, , 10]) > 0L & trans_on == 0L)) {
      stop("Transported fish occur in a rear-week marked transport-off.", call. = FALSE)
    }
  }

  lgr <- as.numeric(target$lgr_spill_pct)
  lgs <- as.numeric(target$lgs_spill_pct)
  lgr_mean <- mean(lgr); lgr_sd <- stats::sd(lgr)
  lgs_mean <- mean(lgs); lgs_sd <- stats::sd(lgs)
  lgr_std <- (lgr - lgr_mean) / lgr_sd
  lgs_std <- (lgs - lgs_mean) / lgs_sd
  outflow <- as.numeric(target$lgr_outflow)
  outflow_mean <- mean(outflow); outflow_sd <- stats::sd(outflow)
  outflow_std <- (outflow - outflow_mean) / outflow_sd
  rear_matrix <- function(x, name, shared, lower = -Inf, upper = Inf) {
    dn <- list(rear_levels, as.character(weeks))
    if (is.null(x)) return(matrix(shared, R, S, byrow = TRUE, dimnames = dn))
    if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(R, S)) ||
        any(!is.finite(x)) || any(x < lower | x > upper)) {
      stop(name, " must be an R x S matrix of finite values, ordered by rear and week.", call. = FALSE)
    }
    if (!is.null(rownames(x)) && !identical(rownames(x), rear_levels)) {
      stop(name, " row names must match the ge_dn rear order.", call. = FALSE)
    }
    matrix(as.numeric(x), R, S, dimnames = dn)
  }
  psi_spill <- rear_matrix(psi_spill, "psi_spill", lgr, 0, 100)
  psi_outflow <- rear_matrix(psi_outflow, "psi_outflow", outflow)
  psi_spill_std <- (psi_spill - lgr_mean) / lgr_sd
  psi_outflow_std <- (psi_outflow - outflow_mean) / outflow_sd

  if (is.null(parent)) parent <- target$parent
  parent <- match(parent, unique(parent))
  n_strat <- max(parent)
  prior_names <- c("alpha_phi_mean", "alpha_phi_sd", "beta_phi_mean", "beta_phi_sd")
  if (is.null(phi_prior)) phi_prior <- target[prior_names]
  phi_prior <- phi_prior[prior_names]

  jd <- c(list(
    R = R, S = S, n = n, N_seen = N_seen, N_lik = nrow(lik),
    lik_r = as.integer(lik[, "row"]), lik_s = as.integer(lik[, "col"]),
    parent = as.integer(parent), n_strat = n_strat,
    lgr_spill_std = lgr_std, lgs_spill_std = lgs_std,
    psi_spill_std = unname(psi_spill_std), psi_outflow_std = unname(psi_outflow_std),
    rear_sd_scale = rear_sd_scale, phi_slope_sd_scale = phi_slope_sd_scale,
    delta_sd = delta_sd, night_sd = night_sd, night_week_sd_scale = night_week_sd_scale,
    p_night_offset = p_night_offset
  ), phi_prior)
  if (!is.null(trans_on)) jd$trans_on <- unname(trans_on)
  transport_block <- if (is.null(trans_on)) {
    "trans[r,s] ~ dbeta(1, 1)"
  } else {
    "trans_free[r,s] ~ dbeta(1, 1)
        trans[r,s] <- trans_on[r,s] * trans_free[r,s]"
  }

  day_block <- if (day_ge == "rear") {
    "logit_psi[r,s] ~ dnorm(
          mu_strat[parent[s]] + rear_psi[r] +
          beta * psi_spill_std[r,s] +
          beta_outflow * psi_outflow_std[r,s],
          tau_psi)"
  } else {
    "logit_psi[r,s] <- eta[s] + rear_psi[r] +
          beta * psi_spill_std[r,s] +
          beta_outflow * psi_outflow_std[r,s]"
  }
  eta_block <- if (day_ge == "shared") {
    "for (s in 1:S) { eta[s] ~ dnorm(mu_strat[parent[s]], tau_psi) }"
  } else {
    ""
  }

  model_string <- paste0("model {
    eps <- 1.0E-9

    # Rear effects (as in fit_ge_rear_model2, full structure)
    tau_rear_psi <- pow(sigma_rear_psi + 1.0E-6, -2)
    tau_rear_p <- pow(sigma_rear_p + 1.0E-6, -2)
    tau_rear_phi <- pow(sigma_rear_phi + 1.0E-6, -2)
    for (r in 1:R) {
      rear_psi_raw[r] ~ dnorm(0, tau_rear_psi)
      rear_p_raw[r] ~ dnorm(0, tau_rear_p)
      rear_phi_raw[r] ~ dnorm(0, tau_rear_phi)
    }
    rear_psi_mean <- mean(rear_psi_raw[])
    rear_p_mean <- mean(rear_p_raw[])
    rear_phi_mean <- mean(rear_phi_raw[])
    for (r in 1:R) {
      rear_psi[r] <- rear_psi_raw[r] - rear_psi_mean
      rear_p[r] <- rear_p_raw[r] - rear_p_mean
      rear_phi[r] <- rear_phi_raw[r] - rear_phi_mean
    }
    sigma_rear_psi ~ dnorm(0, pow(rear_sd_scale, -2)) T(0,)
    sigma_rear_p ~ dnorm(0, pow(rear_sd_scale, -2)) T(0,)
    sigma_rear_phi ~ dnorm(0, pow(rear_sd_scale, -2)) T(0,)

    # Day GE: weekly regression on rear-specific covariates (rear-specific or
    # shared weekly residual). Night GE: day GE plus a weekly night effect.
    tau_psi <- pow(sigma_psi, -2)
    tau_strat <- pow(sigma_strat, -2)
    tau_night <- pow(sigma_night + 1.0E-6, -2)
    ", eta_block, "
    for (r in 1:R) {
      night_mu[r] ~ dnorm(0, pow(night_sd, -2))
      for (s in 1:S) {
        ", day_block, "
        night_eff[r,s] ~ dnorm(night_mu[r], tau_night)
        psi_day[r,s] <- ilogit(logit_psi[r,s])
        psi_night[r,s] <- ilogit(logit_psi[r,s] + night_eff[r,s])
        a[r,s] ~ dbeta(1, 1)
        psi[r,s] <- (1 - a[r,s]) * psi_day[r,s] + a[r,s] * psi_night[r,s]
      }
    }
    sigma_night ~ dnorm(0, pow(night_week_sd_scale, -2)) T(0,)
    for (g in 1:n_strat) { mu_strat[g] ~ dnorm(alpha, tau_strat) }
    alpha ~ dt(0, pow(2, -2), 7)
    beta ~ dt(0, pow(1, -2), 7) T(, 0)
    beta_outflow ~ dt(0, pow(1, -2), 7)
    sigma_psi ~ dnorm(0, pow(0.5, -2)) T(0,)
    sigma_strat ~ dunif(0.05, 3)

    # GRS detection: common to day and night fish within a rear-week
    tau_p <- pow(sigma_p, -2)
    tau_strat_p <- pow(sigma_strat_p, -2)
    for (r in 1:R) {
      for (s in 1:S) {
        logit_p[r,s] ~ dnorm(mu_strat_p[parent[s]] + rear_p[r] + beta_p * lgr_spill_std[s], tau_p)
        p[r,s] <- ilogit(logit_p[r,s])
        p_night[r,s] <- ilogit(logit_p[r,s] + p_night_offset)
      }
    }
    for (g in 1:n_strat) { mu_strat_p[g] ~ dnorm(alpha_p, tau_strat_p) }
    alpha_p ~ dt(0, pow(2, -2), 7)
    beta_p ~ dt(0, pow(1, -2), 7)
    sigma_p ~ dunif(0.05, 3)
    sigma_strat_p ~ dunif(0.05, 3)

    # Recovery: rear-flexible LGS slopes, as in fit_ge_rear_model2
    tau_beta_phi_rear <- pow(sigma_beta_phi_rear + 1.0E-6, -2)
    for (r in 1:R) {
      beta_phi_rear_raw[r] ~ dnorm(0, tau_beta_phi_rear)
      sigma_phi_rear[r] ~ dunif(0.05, 3)
    }
    beta_phi_rear_mean <- mean(beta_phi_rear_raw[])
    for (r in 1:R) {
      beta_phi_rear[r] <- beta_phi + beta_phi_rear_raw[r] - beta_phi_rear_mean
      tau_phi_rear[r] <- pow(sigma_phi_rear[r], -2)
      for (s in 1:S) {
        logit_phi_S[r,s] ~ dnorm(alpha_phi + rear_phi[r] + beta_phi_rear[r] * lgs_spill_std[s],
                                 tau_phi_rear[r])
        phi_S[r,s] <- ilogit(logit_phi_S[r,s])
        phi_B[r,s] <- ilogit(logit_phi_S[r,s] - delta[r,s])
      }
    }
    alpha_phi ~ dnorm(alpha_phi_mean, pow(alpha_phi_sd, -2))
    beta_phi ~ dnorm(beta_phi_mean, pow(beta_phi_sd, -2)) T(, 0)
    sigma_beta_phi_rear ~ dnorm(0, pow(phi_slope_sd_scale, -2)) T(0,)

    # Route offset, free, rear-specific
    tau_rear_delta <- pow(sigma_rear_delta + 1.0E-6, -2)
    for (r in 1:R) { rear_delta_raw[r] ~ dnorm(0, tau_rear_delta) }
    rear_delta_mean <- mean(rear_delta_raw[])
    for (r in 1:R) { rear_delta[r] <- rear_delta_raw[r] - rear_delta_mean }
    sigma_rear_delta ~ dnorm(0, pow(rear_sd_scale, -2)) T(0,)
    delta0 ~ dnorm(0, pow(delta_sd, -2))
    sigma_delta ~ dunif(0, 3)
    tau_delta <- pow(sigma_delta + 1.0E-6, -2)
    for (r in 1:R) {
      for (s in 1:S) { delta[r,s] ~ dnorm(delta0 + rear_delta[r], tau_delta) }
    }

    # Eleven-cell likelihood: day c1 c2 c3 c4 c6, night c1 c2 c3 c4 c6, pooled c5.
    # Night spilled fish are detected with p_night (equal to p when the offset is 0).
    for (r in 1:R) {
      for (s in 1:S) {
        ", transport_block, "
        pi_raw[r,s,1] <- (1-a[r,s]) * psi_day[r,s] * (1-trans[r,s]) * (1-phi_B[r,s]) + eps
        pi_raw[r,s,2] <- (1-a[r,s]) * psi_day[r,s] * (1-trans[r,s]) * phi_B[r,s] + eps
        pi_raw[r,s,3] <- (1-a[r,s]) * (1-psi_day[r,s]) * p[r,s] * (1-phi_S[r,s]) + eps
        pi_raw[r,s,4] <- (1-a[r,s]) * (1-psi_day[r,s]) * p[r,s] * phi_S[r,s] + eps
        pi_raw[r,s,5] <- (1-a[r,s]) * psi_day[r,s] * trans[r,s] + eps
        pi_raw[r,s,6] <- a[r,s] * psi_night[r,s] * (1-trans[r,s]) * (1-phi_B[r,s]) + eps
        pi_raw[r,s,7] <- a[r,s] * psi_night[r,s] * (1-trans[r,s]) * phi_B[r,s] + eps
        pi_raw[r,s,8] <- a[r,s] * (1-psi_night[r,s]) * p_night[r,s] * (1-phi_S[r,s]) + eps
        pi_raw[r,s,9] <- a[r,s] * (1-psi_night[r,s]) * p_night[r,s] * phi_S[r,s] + eps
        pi_raw[r,s,10] <- a[r,s] * psi_night[r,s] * trans[r,s] + eps
        pi_raw[r,s,11] <- ((1-a[r,s]) * (1-psi_day[r,s]) * (1-p[r,s]) +
                           a[r,s] * (1-psi_night[r,s]) * (1-p_night[r,s])) * phi_S[r,s] + eps
        p_seen[r,s] <- sum(pi_raw[r,s,1:11])
        for (cc in 1:11) { pi_obs[r,s,cc] <- pi_raw[r,s,cc] / p_seen[r,s] }
      }
    }
    for (j in 1:N_lik) {
      n[lik_r[j],lik_s[j],1:11] ~ dmulti(pi_obs[lik_r[j],lik_s[j],1:11], N_seen[lik_r[j],lik_s[j]])
    }
  }")

  for (pkg in c("rjags", "coda")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install R package '", pkg, "' and JAGS.", call. = FALSE)
  }
  monitor <- c("psi", "psi_day", "psi_night", "a", "night_eff", "night_mu", "sigma_night",
               "p", "p_night", if (day_ge == "shared") "eta", "phi_S", "phi_B", "trans", "alpha", "beta", "beta_outflow",
               "sigma_psi", "mu_strat", "sigma_strat", "alpha_p", "beta_p", "sigma_p",
               "mu_strat_p", "sigma_strat_p", "alpha_phi", "beta_phi", "beta_phi_rear",
               "sigma_beta_phi_rear", "sigma_phi_rear", "rear_psi", "rear_p",
               "sigma_rear_psi", "sigma_rear_p", "rear_phi", "sigma_rear_phi",
               "delta", "delta0", "sigma_delta", "rear_delta", "sigma_rear_delta")
  if (verbose) {
    message("Day/night six-cell fit (day GE ", day_ge, ", p night offset ", p_night_offset,
            "): ", R, " rears x ", S, " weeks; ", nrow(lik), " nonempty rear-weeks; ", n_chains, " chains x ", n_iter,
            " post-burn iterations, thin ", n_thin, ".")
  }
  con <- textConnection(model_string)
  on.exit(close(con), add = TRUE)
  inits <- lapply(seq_len(n_chains), function(i) {
    list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = as.integer(seed + i))
  })
  jm <- rjags::jags.model(con, data = jd, inits = inits, n.chains = n_chains,
                          n.adapt = n_adapt, quiet = !verbose)
  if (n_burnin > 0L) stats::update(jm, n.iter = n_burnin, progress.bar = if (verbose) "text" else "none")
  samples <- rjags::coda.samples(jm, variable.names = monitor, n.iter = n_iter, thin = n_thin,
                                 progress.bar = if (verbose) "text" else "none")

  mat <- do.call(rbind, lapply(samples, as.matrix))
  psrf <- tryCatch(
    coda::gelman.diag(samples, autoburnin = FALSE, multivariate = FALSE)$psrf,
    error = function(e) matrix(NA_real_, nrow = ncol(mat), ncol = 2L,
                               dimnames = list(colnames(mat), c("Point est.", "Upper C.I."))))
  ess <- tryCatch(coda::effectiveSize(samples),
                  error = function(e) stats::setNames(rep(NA_real_, ncol(mat)), colnames(mat)))
  q <- t(apply(mat, 2L, stats::quantile, probs = c(0.025, 0.5, 0.975)))
  summary <- data.frame(parameter = colnames(mat), mean = colMeans(mat),
                        sd = apply(mat, 2L, stats::sd), median = q[, 2L],
                        lci = q[, 1L], uci = q[, 3L],
                        Rhat = psrf[colnames(mat), 1L], n_eff = ess[colnames(mat)],
                        row.names = NULL, check.names = FALSE)
  bad <- !is.finite(summary$Rhat) | summary$Rhat > rhat_threshold
  core <- grepl(paste0("^psi\\[|^psi_day\\[|^psi_night\\[|^a\\[|^p\\[|^phi_[SB]\\[|",
                       "^delta\\[|^night_mu\\[|^rear_(psi|p|phi|delta)\\[|",
                       "^(alpha|beta|beta_outflow|sigma_psi|sigma_night|sigma_rear_psi|sigma_rear_p)$"),
                summary$parameter)
  diagnostics <- list(threshold = rhat_threshold, all_rhat_pass = !any(bad),
                      core_ge_rhat_pass = !any(bad & core),
                      flagged_core_ge_parameters = summary[bad & core, , drop = FALSE])
  if (any(bad & core)) {
    warning(sum(bad & core), " core parameter(s) have R-hat > ", rhat_threshold, ".", call. = FALSE)
  }

  col_summary <- function(nm, r) {
    cols <- paste0(nm, "[", r, ",", seq_len(S), "]")
    x <- mat[, cols, drop = FALSE]
    qq <- t(apply(x, 2L, stats::quantile, probs = c(0.025, 0.5, 0.975)))
    list(mean = colMeans(x), median = qq[, 2L], lci = qq[, 1L], uci = qq[, 3L],
         Rhat = summary$Rhat[match(cols, summary$parameter)],
         n_eff = summary$n_eff[match(cols, summary$parameter)])
  }
  psi_summary <- do.call(rbind, lapply(seq_len(R), function(r) {
    x <- col_summary("psi", r)
    data.frame(rear = rear_levels[r], week = weeks, parent = parent,
               N_seen = as.numeric(N_seen[r, ]), mean = x$mean, median = x$median,
               lci = x$lci, uci = x$uci, Rhat = x$Rhat, n_eff = x$n_eff, row.names = NULL)
  }))
  period_summary <- do.call(rbind, lapply(seq_len(R), function(r) {
    d <- col_summary("psi_day", r); nn <- col_summary("psi_night", r); aa <- col_summary("a", r)
    data.frame(rear = rear_levels[r], week = weeks,
               n_day = rowSums(n[r, , 1:5]), n_night = rowSums(n[r, , 6:10]), n_c5 = n[r, , 11],
               psi_day = d$mean, psi_day_lci = d$lci, psi_day_uci = d$uci,
               psi_night = nn$mean, psi_night_lci = nn$lci, psi_night_uci = nn$uci,
               night_share = aa$mean, night_share_lci = aa$lci, night_share_uci = aa$uci,
               row.names = NULL)
  }))

  result <- list(
    samples = samples, summary = summary, diagnostics = diagnostics,
    psi_summary = psi_summary, period_summary = period_summary,
    rear_levels = rear_levels, target_rear = target_rear, weeks = weeks,
    day_ge = day_ge, p_night_offset = p_night_offset,
    parent = parent, n = n, N_seen = N_seen,
    spill_mean = lgr_mean, spill_sd = lgr_sd, lgr_spill_pct = lgr, lgr_spill_std = lgr_std,
    outflow_mean = outflow_mean, outflow_sd = outflow_sd, outflow = outflow,
    outflow_std = outflow_std,
    psi_spill_pct = psi_spill, psi_spill_std = psi_spill_std,
    psi_outflow = psi_outflow, psi_outflow_std = psi_outflow_std,
    lgs_spill_mean = lgs_mean, lgs_spill_sd = lgs_sd,
    strat_assign = data.frame(Week = weeks, Collapse = seq_len(S)),
    model_string = model_string, jags_data = jd, trans_on = trans_on,
    settings = list(model = "six-cell rear-type, day/night split",
                    day_ge = day_ge, p_night_offset = p_night_offset,
                    night_sd = night_sd, night_week_sd_scale = night_week_sd_scale,
                    n_iter = n_iter, n_adapt = n_adapt, n_burnin = n_burnin,
                    n_chains = n_chains, n_thin = n_thin, seed = seed))
  class(result) <- c("ge_rear_daynight", "list")
  if (verbose) print(period_summary, row.names = FALSE)
  result
}
