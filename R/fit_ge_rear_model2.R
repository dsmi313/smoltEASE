#' Fit a rear-type hierarchical six-cell GE model
#'
#' Fits the six-cell guidance-efficiency model jointly to two or more rear-type
#' groups. Weekly guidance efficiency (psi) and transport probability are
#' rear-specific. The GE regression contains standardized percent spill,
#' standardized outflow, and their interaction. The seasonal hierarchy,
#' spillway detection,
#' downstream recovery, and route offset are shared across rear groups. Rear
#' effects on logit(psi) are centered and partially pooled.
#'
#' @param ge_data Named list of six-cell data lists returned by
#'   [prep_ge_data2()], one per rear group. For example, `list(H = ge_H,
#'   W = ge_W, U = ge_U)`. Every group must use the same weeks and spill series.
#' @param weeks Increasing ISO week numbers. If `NULL`, they are read from the
#'   target-rear object or from count-matrix row names such as `wk13`.
#' @param outflow Optional numeric vector of mean LGR outflow, one value per
#'   modeled week and in the same row order as `ge_data`. When `NULL`, use
#'   `ge_data[[target_rear]]$lgr_outflow`. Values are standardized internally.
#'   Use a rate (for example, mean hourly flow), not a sum of rates.
#' @param target_rear Rear code whose parent hierarchy supplies the default
#'   grouping and whose GE will usually be passed to SCRAPI2. Default `"W"`.
#' @param parent Optional common parent-group vector. When `NULL`, the parent
#'   vector from `ge_data[[target_rear]]` is used.
#' @param delta_mode Estimate (`"free"`) or fix (`"fixed"`) the route offset.
#' @param delta_sd Prior SD for the mean route offset when estimated.
#' @param rear_sd_scale Half-normal prior scale for the among-rear SD on the
#'   logit-GE scale. Default 1.
#' @param phi_prior Optional list containing `alpha_phi_mean`, `alpha_phi_sd`,
#'   `beta_phi_mean`, and `beta_phi_sd`. Otherwise these are read from the
#'   target-rear data or use the provisional six-cell defaults.
#' @param n_iter,n_adapt,n_burnin,n_chains,n_thin MCMC controls. `n_iter` is
#'   post-burn sampling iterations per chain before thinning.
#' @param seed Nonnegative integer used to seed the JAGS chains.
#' @param rhat_threshold R-hat warning threshold. Default 1.01.
#' @param verbose Print model and diagnostic summaries.
#'
#' @details This is an exploratory partial-pooling model. It treats `U` as an
#' observed unknown-rear classification, not as a distinct biological
#' population or a latent H/W mixture. The shared nuisance-process assumptions
#' should be checked with posterior predictive diagnostics and recovery tests.
#'
#' @return A list containing posterior `samples`, weekly `psi_summary`, full
#' parameter `summary`, convergence `diagnostics`, model inputs and scaling,
#' and a `strat_assign` object compatible with [generate_rear_ge_draws()].
#'
#' @examples
#' \dontrun{
#' fit <- fit_ge_rear_model2(
#'   list(H = ge_H, W = ge_W, U = ge_U),
#'   weeks = 13:26,
#'   outflow = weekly_mean_outflow,
#'   target_rear = "W",
#'   n_iter = 50000,
#'   n_burnin = 20000
#' )
#' subset(fit$psi_summary, rear == "W")
#' }
#' @export
fit_ge_rear_model2 <- function(
    ge_data,
    weeks = NULL,
    outflow = NULL,
    target_rear = "W",
    parent = NULL,
    delta_mode = c("free", "fixed"),
    delta_sd = 0.5,
    rear_sd_scale = 1,
    phi_prior = NULL,
    n_iter = 30000,
    n_adapt = 2000,
    n_burnin = 10000,
    n_chains = 3,
    n_thin = 10,
    seed = 11,
    rhat_threshold = 1.01,
    verbose = TRUE) {

  delta_mode <- match.arg(delta_mode)
  whole <- function(x, name, minimum = 0L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < minimum || x != floor(x) || x > .Machine$integer.max) {
      stop(name, " must be a finite integer >= ", minimum, ".", call. = FALSE)
    }
    as.integer(x)
  }
  positive <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
      stop(name, " must be a finite positive number.", call. = FALSE)
    }
    x
  }

  n_iter <- whole(n_iter, "n_iter", 1L)
  n_adapt <- whole(n_adapt, "n_adapt", 1L)
  n_burnin <- whole(n_burnin, "n_burnin")
  n_chains <- whole(n_chains, "n_chains", 2L)
  n_thin <- whole(n_thin, "n_thin", 1L)
  seed <- whole(seed, "seed")
  delta_sd <- positive(delta_sd, "delta_sd")
  rear_sd_scale <- positive(rear_sd_scale, "rear_sd_scale")
  rhat_threshold <- positive(rhat_threshold, "rhat_threshold")
  if (rhat_threshold <= 1) {
    stop("rhat_threshold must exceed 1.", call. = FALSE)
  }
  if (n_iter %% n_thin != 0L || n_iter / n_thin < 4L) {
    stop("n_iter must be divisible by n_thin and retain at least four draws per chain.",
         call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.list(ge_data) || length(ge_data) < 1L ||
      is.null(names(ge_data)) || any(!nzchar(names(ge_data))) ||
      anyDuplicated(names(ge_data))) {
    stop("ge_data must be a uniquely named list containing at least one rear group.",
         call. = FALSE)
  }
  rear_levels <- names(ge_data)
  if (!is.character(target_rear) || length(target_rear) != 1L ||
      !target_rear %in% rear_levels) {
    stop("target_rear must name one element of ge_data.", call. = FALSE)
  }

  target <- ge_data[[target_rear]]
  required <- c("n", "lgr_spill_pct", "lgs_spill_pct")
  missing_target <- setdiff(required, names(target))
  if (length(missing_target)) {
    stop("Target-rear data are missing: ", paste(missing_target, collapse = ", "),
         ".", call. = FALSE)
  }
  target_n <- as.matrix(target$n)
  if (ncol(target_n) != 6L || nrow(target_n) < 2L) {
    stop("Each rear group must contain an S x 6 count matrix (S >= 2).",
         call. = FALSE)
  }
  S <- nrow(target_n)
  R <- length(rear_levels)

  if (is.null(weeks)) {
    for (nm in c("weeks", "Week", "week")) {
      if (is.null(weeks) && !is.null(target[[nm]])) weeks <- target[[nm]]
    }
    rn <- rownames(target_n)
    if (is.null(weeks) && !is.null(rn) && all(grepl("^wk[0-9]+$", rn))) {
      weeks <- as.integer(sub("^wk", "", rn))
    }
  }
  if (is.null(weeks) || !is.numeric(weeks) || length(weeks) != S ||
      any(!is.finite(weeks)) || any(weeks != floor(weeks)) ||
      any(weeks < 1 | weeks > 53) || any(diff(weeks) <= 0)) {
    stop("weeks must be S unique increasing ISO week numbers.", call. = FALSE)
  }
  weeks <- as.integer(weeks)

  n <- array(0L, dim = c(R, S, 6L),
             dimnames = list(rear_levels, paste0("wk", weeks), paste0("c", 1:6)))
  for (r in seq_len(R)) {
    x <- ge_data[[r]]
    missing <- setdiff(required, names(x))
    if (length(missing)) {
      stop("Rear group ", rear_levels[r], " is missing: ",
           paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    nr <- as.matrix(x$n)
    if (!identical(dim(nr), c(S, 6L)) || !is.numeric(nr) ||
        any(!is.finite(nr)) || any(nr < 0) || any(nr != floor(nr))) {
      stop("Rear group ", rear_levels[r],
           " must have an S x 6 matrix of nonnegative integer counts.",
           call. = FALSE)
    }
    if (!isTRUE(all.equal(as.numeric(x$lgr_spill_pct),
                          as.numeric(target$lgr_spill_pct))) ||
        !isTRUE(all.equal(as.numeric(x$lgs_spill_pct),
                          as.numeric(target$lgs_spill_pct)))) {
      stop("All rear groups must use identical LGR and LGS spill values.",
           call. = FALSE)
    }
    n[r, , ] <- nr
  }
  storage.mode(n) <- "integer"
  N_seen <- apply(n, c(1L, 2L), sum)
  if (!any(N_seen > 0)) stop("No six-cell histories were supplied.", call. = FALSE)
  lik <- which(N_seen > 0, arr.ind = TRUE)

  lgr <- as.numeric(target$lgr_spill_pct)
  lgs <- as.numeric(target$lgs_spill_pct)
  for (item in list(LGR = lgr, LGS = lgs)) {
    if (length(item) != S || any(!is.finite(item)) || any(item < 0 | item > 100) ||
        !is.finite(stats::sd(item)) || stats::sd(item) == 0) {
      stop("LGR and LGS spill must each contain S finite percentages with nonzero SD.",
           call. = FALSE)
    }
  }
  lgr_mean <- mean(lgr)
  lgr_sd <- stats::sd(lgr)
  lgs_mean <- mean(lgs)
  lgs_sd <- stats::sd(lgs)
  lgr_std <- (lgr - lgr_mean) / lgr_sd
  lgs_std <- (lgs - lgs_mean) / lgs_sd
  if (is.null(outflow)) outflow <- target$lgr_outflow
  if (!is.numeric(outflow) || length(outflow) != S ||
      any(!is.finite(outflow)) || !is.finite(stats::sd(outflow)) ||
      stats::sd(outflow) == 0) {
    stop("outflow must contain S finite values with nonzero SD.", call. = FALSE)
  }
  outflow_mean <- mean(outflow)
  outflow_sd <- stats::sd(outflow)
  outflow_std <- (outflow - outflow_mean) / outflow_sd
  interaction_std <- lgr_std * outflow_std

  if (is.null(parent)) parent <- target$parent
  if (is.null(parent)) parent <- seq_len(S)
  if (!is.numeric(parent) || length(parent) != S || any(!is.finite(parent)) ||
      any(parent < 1 | parent != floor(parent)) || anyDuplicated(rle(parent)$values)) {
    stop("parent must be S positive integers defining contiguous groups.",
         call. = FALSE)
  }
  parent <- match(parent, unique(parent))
  n_strat <- max(parent)

  prior_names <- c("alpha_phi_mean", "alpha_phi_sd",
                   "beta_phi_mean", "beta_phi_sd")
  if (is.null(phi_prior)) {
    if (all(prior_names %in% names(target))) {
      phi_prior <- target[prior_names]
      prior_source <- paste0("ge_data[['", target_rear, "']]")
    } else {
      phi_prior <- list(alpha_phi_mean = 0, alpha_phi_sd = 2,
                        beta_phi_mean = 0, beta_phi_sd = 1)
      prior_source <- "provisional six-cell defaults"
    }
  } else {
    prior_source <- "explicit phi_prior argument"
  }
  if (!is.list(phi_prior) || !all(prior_names %in% names(phi_prior))) {
    stop("phi_prior must contain ", paste(prior_names, collapse = ", "), ".",
         call. = FALSE)
  }
  phi_prior <- phi_prior[prior_names]
  for (nm in prior_names) {
    value <- phi_prior[[nm]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        (grepl("_sd$", nm) && value <= 0)) {
      stop("Invalid phi prior value: ", nm, ".", call. = FALSE)
    }
  }

  jd <- c(list(
    R = R, S = S, n = n, N_seen = N_seen,
    N_lik = nrow(lik), lik_r = as.integer(lik[, "row"]),
    lik_s = as.integer(lik[, "col"]), parent = as.integer(parent),
    n_strat = n_strat, lgr_spill_std = lgr_std,
    outflow_std = outflow_std, interaction_std = interaction_std,
    lgs_spill_std = lgs_std, rear_sd_scale = rear_sd_scale
  ), phi_prior)

  delta_block <- if (delta_mode == "free") {
    jd$delta_sd <- delta_sd
    "
    delta0 ~ dnorm(0, pow(delta_sd, -2))
    sigma_delta ~ dunif(0, 3)
    tau_delta <- pow(sigma_delta + 1.0E-6, -2)
    for (s in 1:S) { delta[s] ~ dnorm(delta0, tau_delta) }
    "
  } else {
    "for (s in 1:S) { delta[s] <- 0 }"
  }

  rear_block <- if (R > 1L) {
    "
    tau_rear <- pow(sigma_rear_psi + 1.0E-6, -2)
    for (r in 1:R) { rear_psi_raw[r] ~ dnorm(0, tau_rear) }
    rear_mean <- mean(rear_psi_raw[])
    for (r in 1:R) { rear_psi[r] <- rear_psi_raw[r] - rear_mean }
    sigma_rear_psi ~ dnorm(0, pow(rear_sd_scale, -2)) T(0,)
    "
  } else {
    "rear_psi[1] <- 0"
  }

  model_string <- paste0("model {
    eps <- 1.0E-9

    # Rear-specific GE with centered, partially pooled rear effects.
    tau_psi <- pow(sigma_psi, -2)
    tau_strat <- pow(sigma_strat, -2)
    ", rear_block, "
    for (r in 1:R) {
      for (s in 1:S) {
        logit_psi[r,s] ~ dnorm(
          mu_strat[parent[s]] + rear_psi[r] +
          beta * lgr_spill_std[s] +
          beta_outflow * outflow_std[s] +
          beta_interaction * interaction_std[s],
          tau_psi)
        psi[r,s] <- ilogit(logit_psi[r,s])
      }
    }
    for (g in 1:n_strat) { mu_strat[g] ~ dnorm(alpha, tau_strat) }
    alpha ~ dt(0, pow(2, -2), 7)
    beta ~ dt(0, pow(1, -2), 7) T(, 0)
    beta_outflow ~ dt(0, pow(1, -2), 7)
    beta_interaction ~ dt(0, pow(1, -2), 7)
    sigma_psi ~ dnorm(0, pow(0.5, -2)) T(0,)
    sigma_strat ~ dunif(0.05, 3)

    # Shared spillway-array detection.
    tau_p <- pow(sigma_p, -2)
    tau_strat_p <- pow(sigma_strat_p, -2)
    for (s in 1:S) {
      logit_p[s] ~ dnorm(mu_strat_p[parent[s]] + beta_p * lgr_spill_std[s], tau_p)
      p[s] <- ilogit(logit_p[s])
    }
    for (g in 1:n_strat) { mu_strat_p[g] ~ dnorm(alpha_p, tau_strat_p) }
    alpha_p ~ dt(0, pow(2, -2), 7)
    beta_p ~ dt(0, pow(1, -2), 7)
    sigma_p ~ dunif(0.05, 3)
    sigma_strat_p ~ dunif(0.05, 3)

    # Shared downstream recovery and route offset.
    tau_phi <- pow(sigma_phi, -2)
    for (s in 1:S) {
      logit_phi_S[s] ~ dnorm(alpha_phi + beta_phi * lgs_spill_std[s], tau_phi)
      phi_S[s] <- ilogit(logit_phi_S[s])
      phi_B[s] <- ilogit(logit_phi_S[s] - delta[s])
    }
    alpha_phi ~ dnorm(alpha_phi_mean, pow(alpha_phi_sd, -2))
    beta_phi ~ dnorm(beta_phi_mean, pow(beta_phi_sd, -2)) T(, 0)
    sigma_phi ~ dunif(0.05, 3)
    ", delta_block, "

    # Transport and the six-cell likelihood are rear-specific.
    for (r in 1:R) {
      for (s in 1:S) {
        trans[r,s] ~ dbeta(1, 1)
        pi_raw[r,s,1] <- psi[r,s] * (1-trans[r,s]) * (1-phi_B[s]) + eps
        pi_raw[r,s,2] <- psi[r,s] * (1-trans[r,s]) * phi_B[s] + eps
        pi_raw[r,s,3] <- (1-psi[r,s]) * p[s] * (1-phi_S[s]) + eps
        pi_raw[r,s,4] <- (1-psi[r,s]) * p[s] * phi_S[s] + eps
        pi_raw[r,s,5] <- (1-psi[r,s]) * (1-p[s]) * phi_S[s] + eps
        pi_raw[r,s,6] <- psi[r,s] * trans[r,s] + eps
        p_seen[r,s] <- sum(pi_raw[r,s,1:6])
        for (cc in 1:6) { pi_obs[r,s,cc] <- pi_raw[r,s,cc] / p_seen[r,s] }
      }
    }
    for (j in 1:N_lik) {
      n[lik_r[j],lik_s[j],1:6] ~ dmulti(
        pi_obs[lik_r[j],lik_s[j],1:6], N_seen[lik_r[j],lik_s[j]])
    }
  }")

  for (pkg in c("rjags", "coda")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Install R package '", pkg,
           "' and ensure the separate JAGS application is installed.",
           call. = FALSE)
    }
  }
  monitor <- c(
    "psi", "p", "phi_S", "phi_B",
    "trans", "pi_obs", "alpha", "beta", "sigma_psi", "mu_strat",
    "beta_outflow", "beta_interaction",
    "sigma_strat", "alpha_p", "beta_p", "sigma_p", "mu_strat_p",
    "sigma_strat_p", "alpha_phi", "beta_phi", "sigma_phi"
  )
  if (R > 1L) monitor <- c(monitor, "rear_psi", "sigma_rear_psi")
  if (delta_mode == "free") {
    monitor <- c(monitor, "delta", "delta0", "sigma_delta")
  }
  if (verbose) {
    message("Rear-type six-cell fit: ", R, " rear groups x ", S,
            " weeks; ", nrow(lik), " nonempty rear-week rows; ",
            n_chains, " chains x ", n_iter,
            " post-burn iterations, thin ", n_thin, ".")
    message("Target rear: ", target_rear,
            "; phi-prior source: ", prior_source, ".")
  }

  con <- textConnection(model_string)
  on.exit(close(con), add = TRUE)
  inits <- lapply(seq_len(n_chains), function(i) {
    list(.RNG.name = "base::Mersenne-Twister",
         .RNG.seed = as.integer(seed + i))
  })
  jm <- rjags::jags.model(con, data = jd, inits = inits,
                          n.chains = n_chains, n.adapt = n_adapt,
                          quiet = !verbose)
  if (n_burnin > 0L) {
    stats::update(jm, n.iter = n_burnin,
                  progress.bar = if (verbose) "text" else "none")
  }
  samples <- rjags::coda.samples(
    jm, variable.names = monitor, n.iter = n_iter, thin = n_thin,
    progress.bar = if (verbose) "text" else "none"
  )

  mat <- do.call(rbind, lapply(samples, as.matrix))
  psrf <- tryCatch(
    coda::gelman.diag(samples, autoburnin = FALSE,
                      multivariate = FALSE)$psrf,
    error = function(e) matrix(
      NA_real_, nrow = ncol(mat), ncol = 2L,
      dimnames = list(colnames(mat), c("Point est.", "Upper C.I."))))
  ess <- tryCatch(coda::effectiveSize(samples), error = function(e) {
    stats::setNames(rep(NA_real_, ncol(mat)), colnames(mat))
  })
  q <- t(apply(mat, 2L, stats::quantile, probs = c(0.025, 0.5, 0.975)))
  summary <- data.frame(
    parameter = colnames(mat), mean = colMeans(mat),
    sd = apply(mat, 2L, stats::sd), median = q[, 2L],
    lci = q[, 1L], uci = q[, 3L],
    Rhat = psrf[colnames(mat), 1L],
    Rhat_upper = psrf[colnames(mat), 2L],
    n_eff = ess[colnames(mat)], row.names = NULL,
    check.names = FALSE
  )
  bad <- !is.finite(summary$Rhat) | summary$Rhat > rhat_threshold
  core <- grepl(
    "^psi\\[|^rear_psi\\[|^(alpha|beta|beta_outflow|beta_interaction|sigma_psi|sigma_rear_psi)$",
    summary$parameter)
  diagnostics <- list(
    method = "classical coda R-hat, autoburnin=FALSE; coda ESS",
    threshold = rhat_threshold,
    all_rhat_pass = !any(bad),
    core_ge_rhat_pass = !any(bad & core),
    flagged_parameters = summary[bad, , drop = FALSE],
    flagged_core_ge_parameters = summary[bad & core, , drop = FALSE]
  )
  if (any(bad)) {
    warning(sum(bad), " parameter(s) have R-hat > ", rhat_threshold,
            " or unavailable R-hat. Inspect chains before using GE.",
            call. = FALSE)
  }

  psi_summary <- do.call(rbind, lapply(seq_len(R), function(r) {
    cols <- paste0("psi[", r, ",", seq_len(S), "]")
    qq <- t(apply(mat[, cols, drop = FALSE], 2L, stats::quantile,
                  probs = c(0.025, 0.5, 0.975)))
    data.frame(
      rear = rear_levels[r], week = weeks, parent = parent,
      N_seen = as.numeric(N_seen[r, ]),
      mean = colMeans(mat[, cols, drop = FALSE]),
      median = qq[, 2L], lci = qq[, 1L], uci = qq[, 3L],
      Rhat = summary$Rhat[match(cols, summary$parameter)],
      n_eff = summary$n_eff[match(cols, summary$parameter)],
      row.names = NULL
    )
  }))

  rear_summary <- summary[grepl("^rear_psi\\[", summary$parameter), , drop = FALSE]
  if (R > 1L) {
    rear_summary$rear <- rear_levels
    rear_summary <- rear_summary[, c("rear", setdiff(names(rear_summary), "rear"))]
  }

  result <- list(
    samples = samples, summary = summary, diagnostics = diagnostics,
    psi_summary = psi_summary, rear_summary = rear_summary,
    rear_levels = rear_levels, target_rear = target_rear,
    weeks = weeks, parent = parent, n = n, N_seen = N_seen,
    spill_mean = lgr_mean, spill_sd = lgr_sd,
    lgr_spill_pct = lgr, lgr_spill_std = lgr_std,
    outflow_mean = outflow_mean, outflow_sd = outflow_sd,
    outflow = as.numeric(outflow), outflow_std = outflow_std,
    interaction_std = interaction_std,
    lgs_spill_mean = lgs_mean, lgs_spill_sd = lgs_sd,
    strat_assign = data.frame(Week = weeks, Collapse = seq_len(S)),
    model_string = model_string, jags_data = jd,
    settings = list(
      model = "six-cell rear-type partial pooling",
      delta_mode = delta_mode, delta_sd = delta_sd,
      rear_sd_scale = rear_sd_scale, phi_prior_source = prior_source,
      psi_covariates = c("percent spill", "outflow", "spill x outflow"),
      shared_processes = c("seasonal pattern", "covariate slopes", "p",
                           "phi_S", "phi_B", "delta"),
      rear_specific_processes = c("psi", "trans"),
      n_iter = n_iter, n_adapt = n_adapt, n_burnin = n_burnin,
      n_chains = n_chains, n_thin = n_thin, seed = seed
    )
  )
  class(result) <- c("ge_rear_model2", "list")
  if (verbose) print(psi_summary, row.names = FALSE)
  result
}
