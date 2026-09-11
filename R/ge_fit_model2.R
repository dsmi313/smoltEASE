# Six-cell extension for smoltEASE. The five-cell fit_ge_model() is unchanged.
# This file only defines a function; sourcing it does not fit or read any files.
#
# Standalone use, with the `dat` object from your six-cell build_dat():
#   source("fit_ge_model2.R")
#   ge_fit <- fit_ge_model2(dat, n_iter = 30000)
#   ge_fit$psi_summary
#
# To integrate into the package, put this file in R/, export fit_ge_model2
# (the roxygen tag below does this), and declare rjags and coda in Suggests.
# JAGS must also be installed separately on the machine running the fit.

#' Fit the wild six-cell guidance-efficiency model
#'
#' Fits the supplied six-cell model, including an explicit transported-fish
#' category. Both psi and p have nested parent-group/weekly random effects.
#' Downstream recovery is anchored on the spillway route:
#' phi_B = plogis(logit_phi_S - delta). The default estimates the route offset;
#' delta_mode = "fixed" sets it to zero. No likelihood-count cap is applied.
#'
#' @param ge_data Either the list returned by the supplied six-cell build_dat()
#'   (with n, lgr_spill_pct, lgs_spill_pct, and optionally parent and the four
#'   phi-prior settings), or a data frame with six count columns c1,...,c6
#'   (alternatively n_h1,...,n_h6), spill_val (LGR percent spill), and
#'   lgs_spill_val (LGS percent spill). A supplied stratum_idx must be 1,...,S
#'   in row order; a supplied stratum must contain unique increasing integers.
#'   Optional parent assigns consecutive rows to coarser groups. Data must
#'   already be wild-only, with mutually exclusive histories and transported
#'   fish removed from cells 1-5. The older prep_ge_data() output is NOT a
#'   six-cell input: adding a zero sixth column does not reconstruct transport.
#' @param weeks Optional increasing ISO week numbers, one per input row. If
#'   omitted, read from ge_data$weeks, ge_data$Week, ge_data$week, or matrix
#'   row names such as wk13. No global WEEKS variable is consulted. Needed to
#'   build the convenient Week/Collapse mapping in the returned strat_assign.
#' @param delta_mode "free" (default) or "fixed".
#' @param delta_sd Prior SD of delta0 on the logit scale; default 0.5.
#' @param parent_floor If parent is absent, merge adjacent rows forward until
#'   BOTH accumulated c4+c5 and c3+c4 exceed this floor (strictly >, matching
#'   the supplied script). Attach a trailing remainder to the last group.
#'   Default 10. This GE grouping is not the SCRAPI composition Collapse map.
#' @param phi_prior Optional list with alpha_phi_mean, alpha_phi_sd,
#'   beta_phi_mean, and beta_phi_sd. Otherwise use those four values in
#'   ge_data when all are present; otherwise use 0, 2, 0, 1, respectively.
#'   Those defaults reproduce the PROVISIONAL priors in the supplied six-cell
#'   script. They are not externally calibrated priors. beta_phi is truncated
#'   above at zero, as in that script. Review these settings for production.
#' @param n_iter Post-burn sampling iterations PER CHAIN, BEFORE thinning.
#'   Default 30000; with n_thin = 10 and three chains this yields 9000 draws.
#' @param n_adapt Adaptation iterations, separate from burn-in. Default 2000.
#' @param n_burnin Discarded burn-in iterations. Default 10000.
#' @param n_chains Number of independently seeded chains; at least 2. Default 3.
#' @param n_thin Thinning interval. Default 10.
#' @param seed Nonnegative integer. Chain i uses JAGS seed seed + i. Default 11.
#' @param rhat_threshold Warn when the classical coda R-hat exceeds this value
#'   or is unavailable. Default 1.01. More iterations do not ensure convergence.
#' @param verbose Print fit settings and the weekly psi summary. Default TRUE.
#'
#' @details Cell order and unnormalised probabilities (before the original
#'   1e-9 numerical offset, which is retained in every cell):
#'   1. GRJ only, not transported: psi * (1-trans) * (1-phi_B).
#'   2. GRJ and GOJ, not transported: psi * (1-trans) * phi_B.
#'   3. GRS only: (1-psi) * p * (1-phi_S).
#'   4. GRS and GOJ: (1-psi) * p * phi_S.
#'   5. GOJ only: (1-psi) * (1-p) * phi_S.
#'   6. Transported: psi * trans.
#'
#'   The likelihood conditions on detection in one of these six categories.
#'   phi_B and trans cancel from the normalising constant. This cancellation
#'   alone does NOT prove posterior independence of psi and delta: phi_S is
#'   shared with the route-offset model. The function makes no such guarantee.
#'
#'   Rows with zero histories are retained in the full hierarchy and have
#'   posterior psi draws; only their multinomial likelihood is omitted.
#'   Scaling uses every supplied row, exactly as in the six-cell build_dat().
#'   If the list contains standardised spill vectors, these are checked against
#'   the raw percent-spill vectors rather than silently using another scale.
#'
#'   Compatibility with smoltEASE::generate_ge_draws(): samples contains psi
#'   columns in NUMERIC index order. obs_strata is a compatibility view of ALL
#'   modelled rows, including zero-history rows, so their fitted hierarchical
#'   draws are used rather than regenerated independently. likelihood_strata
#'   identifies rows actually entering the multinomial. n_pool in the returned
#'   ge_data is N_seen, used only for the existing generator's fallback weights;
#'   it is not the old five-cell upstream-pool count. Match every passage date
#'   to a fitted stratum to avoid the generator's season-mean fallback.
#'
#'   For the original weekly-constant handoff use daily_spill = NULL and
#'   clip_to_ci = FALSE in generate_ge_draws(). Daily downscaling is optional
#'   and is an additional modelling choice, not part of this six-cell fit.
#'
#' @return A list with samples (coda mcmc.list), ge_data, obs_strata,
#'   likelihood_strata, scaling metadata, strat_assign, psi_summary, summary,
#'   diagnostics, model_string, jags_data, and settings. Also contains
#'   sims.list, Rhat, and n.eff for the supplied psi_table() convention. Rhat
#'   is the classical univariate Gelman-Rubin point estimate from coda, with
#'   autoburnin = FALSE, not rank-normalised split R-hat. ESS is coda's
#'   effectiveSize, not bulk/tail ESS. Inspect chain traces as well.
#'
#' @examples
#' \dontrun{
#' # dat already contains the SIX mutually exclusive history counts.
#' ge_fit <- fit_ge_model2(dat, weeks = 13:26, n_iter = 30000)
#' print(ge_fit$psi_summary)
#'
#' # Only proceed after reviewing diagnostics and the phi priors.
#' pass_dates <- as.Date(STHD.FPC$SampleEndDate, format = "%m/%d/%Y")
#' stopifnot(all(as.integer(format(pass_dates, "%V")) %in% ge_fit$weeks))
#' ge_draws <- smoltEASE::generate_ge_draws(
#'   ge_fit, strat_assign = ge_fit$strat_assign, pass_dates = pass_dates,
#'   B = 5000, daily_spill = NULL, clip_to_ci = FALSE, seed = 11)
#' # Pass ge_draws to smoltEASE::SCRAPI2(..., geDraws = ge_draws).
#' }
#' @export
fit_ge_model2 <- function(ge_data, weeks = NULL,
                          delta_mode = c("free", "fixed"), delta_sd = 0.5,
                          parent_floor = 10, phi_prior = NULL,
                          n_iter = 30000, n_adapt = 2000, n_burnin = 10000,
                          n_chains = 3, n_thin = 10, seed = 11,
                          rhat_threshold = 1.01, verbose = TRUE) {
  delta_mode <- match.arg(delta_mode)
  whole <- function(x, name, minimum = 0) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < minimum || x != floor(x) || x > .Machine$integer.max)
      stop(name, " must be a finite integer >= ", minimum, ".", call. = FALSE)
    as.integer(x)
  }
  finite_scalar <- function(x, name, positive = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        (positive && x <= 0))
      stop(name, " must be a finite", if (positive) " positive" else "",
           " number.", call. = FALSE)
    x
  }
  n_iter <- whole(n_iter, "n_iter", 1)
  n_adapt <- whole(n_adapt, "n_adapt", 1)
  n_burnin <- whole(n_burnin, "n_burnin")
  n_chains <- whole(n_chains, "n_chains", 2)
  n_thin <- whole(n_thin, "n_thin", 1)
  seed <- whole(seed, "seed")
  parent_floor <- whole(parent_floor, "parent_floor")
  delta_sd <- finite_scalar(delta_sd, "delta_sd", TRUE)
  rhat_threshold <- finite_scalar(rhat_threshold, "rhat_threshold", TRUE)
  if (rhat_threshold <= 1) stop("rhat_threshold must exceed 1.", call. = FALSE)
  if (as.double(seed) + n_chains > .Machine$integer.max)
    stop("seed + n_chains exceeds the JAGS integer seed range.", call. = FALSE)
  if (n_iter %% n_thin != 0L || n_iter / n_thin < 4)
    stop("n_iter must be divisible by n_thin and retain at least 4 draws per chain.",
         call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("verbose must be TRUE or FALSE.", call. = FALSE)

  # Normalise only already-built SIX-CELL data. Never reconstruct c6 from c1.
  cell_names <- paste0("c", 1:6)
  hist_names <- paste0("n_h", 1:6)
  if (!is.list(ge_data))
    stop("ge_data must be a six-cell data frame or build_dat() list.", call. = FALSE)
  if (is.data.frame(ge_data)) {
    columns <- if (all(cell_names %in% names(ge_data))) cell_names else hist_names
    missing <- setdiff(c(columns, "spill_val", "lgs_spill_val"), names(ge_data))
    if (length(missing))
      stop("Six-cell input is missing: ", paste(missing, collapse = ", "),
           ". The five-cell prep_ge_data() output cannot be used unchanged.",
           call. = FALSE)
    n <- as.matrix(ge_data[, columns, drop = FALSE])
    if (all(c(cell_names, hist_names) %in% names(ge_data)) &&
        !isTRUE(all.equal(unname(n),
                         unname(as.matrix(ge_data[, hist_names, drop = FALSE])))))
      stop("c1:c6 and n_h1:n_h6 disagree; supply one unambiguous count table.",
           call. = FALSE)
    lgr <- ge_data$spill_val
    lgs <- ge_data$lgs_spill_val
  } else {
    required <- c("n", "lgr_spill_pct", "lgs_spill_pct")
    missing <- setdiff(required, names(ge_data))
    if (length(missing))
      stop("Six-cell list is missing: ", paste(missing, collapse = ", "),
           ". Raw percent-spill vectors are needed for the GE-draw handoff.",
           call. = FALSE)
    n <- as.matrix(ge_data$n)
    if (ncol(n) != 6L) stop("n must have exactly SIX columns.", call. = FALSE)
    if (!is.null(colnames(n))) {
      if (setequal(colnames(n), cell_names)) n <- n[, cell_names, drop = FALSE]
      else if (setequal(colnames(n), hist_names)) n <- n[, hist_names, drop = FALSE]
      else stop("Name n's columns c1,...,c6 or n_h1,...,n_h6, or leave them unnamed",
                " in the documented cell order.", call. = FALSE)
    }
    lgr <- ge_data$lgr_spill_pct
    lgs <- ge_data$lgs_spill_pct
  }
  if (!is.numeric(n) || nrow(n) < 2L || ncol(n) != 6L ||
      any(!is.finite(n)) || any(n < 0) || any(n != floor(n)))
    stop("Counts must be a numeric S x 6 matrix of nonnegative integers (S >= 2).",
         call. = FALSE)
  S <- nrow(n)
  N_seen <- rowSums(n)
  if (!any(N_seen > 0) || any(N_seen > .Machine$integer.max))
    stop("At least one row must have histories; row totals must fit in an integer.",
         call. = FALSE)
  for (name in c("S", "N_seen")) {
    expected <- if (name == "S") S else N_seen
    if (!is.null(ge_data[[name]]) &&
        !isTRUE(all.equal(as.numeric(ge_data[[name]]), as.numeric(expected))))
      stop("ge_data$", name, " disagrees with the six-cell count matrix.", call. = FALSE)
  }
  storage.mode(n) <- "integer"
  colnames(n) <- cell_names
  for (name in c("lgr", "lgs")) {
    x <- if (name == "lgr") lgr else lgs
    if (!is.numeric(x) || length(x) != S || any(!is.finite(x)) ||
        any(x < 0 | x > 100) || !is.finite(stats::sd(x)) || stats::sd(x) == 0)
      stop(name, " spill must have S finite percent values in [0,100] and nonzero SD.",
           call. = FALSE)
  }
  lgr_mean <- mean(lgr); lgr_sd <- stats::sd(lgr)
  lgs_mean <- mean(lgs); lgs_sd <- stats::sd(lgs)
  lgr_std <- (lgr - lgr_mean) / lgr_sd
  lgs_std <- (lgs - lgs_mean) / lgs_sd
  for (name in c("lgr_spill_std", "lgs_spill_std")) {
    supplied <- ge_data[[name]]
    expected <- if (name == "lgr_spill_std") lgr_std else lgs_std
    if (!is.null(supplied) &&
        (!is.numeric(supplied) || length(supplied) != S ||
         any(!is.finite(supplied)) || max(abs(supplied - expected)) > 1e-8))
      stop(name, " disagrees with standardisation of the raw spill across all rows.",
           call. = FALSE)
  }

  # Model-row identity is distinct from the coarser parent grouping.
  if (is.null(weeks)) {
    for (name in c("weeks", "Week", "week"))
      if (is.null(weeks) && !is.null(ge_data[[name]])) weeks <- ge_data[[name]]
    rn <- rownames(n)
    if (is.null(weeks) && !is.null(rn) && all(grepl("^wk[0-9]+$", rn)))
      weeks <- as.integer(sub("^wk", "", rn))
  }
  if (!is.null(weeks)) {
    if (!is.numeric(weeks) || length(weeks) != S || any(!is.finite(weeks)) ||
        any(weeks != floor(weeks)) || any(weeks < 1 | weeks > 53) ||
        any(diff(weeks) <= 0))
      stop("weeks must be S unique increasing ISO week numbers (one migration year).",
           call. = FALSE)
    weeks <- as.integer(weeks)
  }
  if (!is.null(ge_data$stratum_idx) &&
      !identical(as.numeric(ge_data$stratum_idx), as.numeric(seq_len(S))))
    stop("stratum_idx must be 1,...,S in input row order. Do not reuse a collapsed index.",
         call. = FALSE)
  labels <- ge_data$stratum
  if (is.null(labels)) labels <- if (!is.null(weeks)) weeks else seq_len(S)
  if (!is.numeric(labels) || length(labels) != S || any(!is.finite(labels)) ||
      any(labels != floor(labels)) || any(diff(labels) <= 0))
    stop("stratum must contain S unique increasing integer labels.", call. = FALSE)

  parent <- ge_data$parent
  if (is.null(parent)) {
    parent <- integer(S)
    group <- 1L; id_p <- 0; id_phi <- 0; open <- integer()
    for (s in seq_len(S)) {
      id_p <- id_p + n[s, 4] + n[s, 5]
      id_phi <- id_phi + n[s, 3] + n[s, 4]
      open <- c(open, s)
      if (id_p > parent_floor && id_phi > parent_floor) {
        parent[open] <- group
        group <- group + 1L; id_p <- 0; id_phi <- 0; open <- integer()
      }
    }
    if (length(open)) parent[open] <- max(1L, group - 1L)
  }
  if (!is.numeric(parent) || length(parent) != S || any(!is.finite(parent)) ||
      any(parent < 1 | parent != floor(parent)) || anyDuplicated(rle(parent)$values))
    stop("parent must be S positive integers defining contiguous groups.", call. = FALSE)
  parent <- match(parent, unique(parent))
  n_strat <- max(parent)

  # Retain the attached model's priors, including its provisional phi settings.
  prior_names <- c("alpha_phi_mean", "alpha_phi_sd", "beta_phi_mean", "beta_phi_sd")
  prior_source <- "explicit phi_prior argument"
  if (is.null(phi_prior)) {
    present <- prior_names %in% names(ge_data)
    if (any(present) && !all(present))
      stop("Supply all four phi-prior settings or none of them.", call. = FALSE)
    if (all(present)) {
      phi_prior <- ge_data[prior_names]
      prior_source <- "ge_data"
    } else {
      phi_prior <- list(alpha_phi_mean = 0, alpha_phi_sd = 2,
                        beta_phi_mean = 0, beta_phi_sd = 1)
      prior_source <- "provisional defaults from the supplied six-cell script"
    }
  }
  if (!is.list(phi_prior) || !all(prior_names %in% names(phi_prior)))
    stop("phi_prior must be a list containing ", paste(prior_names, collapse = ", "),
         ".", call. = FALSE)
  phi_prior <- phi_prior[prior_names]
  for (name in prior_names)
    phi_prior[[name]] <- finite_scalar(phi_prior[[name]], name, grepl("_sd$", name))

  gd <- data.frame(stratum_idx = seq_len(S), stratum = labels,
                   spill_val = as.numeric(lgr), lgs_spill_val = as.numeric(lgs),
                   parent = parent, N_seen = as.integer(N_seen),
                   has_history = N_seen > 0, n_pool = as.integer(N_seen))
  gd[hist_names] <- as.data.frame(n)
  if (!is.null(weeks)) gd$Week <- weeks
  lik_idx <- which(N_seen > 0)
  jd <- c(list(S = S, S_lik = length(lik_idx), lik_idx = as.array(lik_idx),
               n = n, N_seen = as.integer(N_seen), parent = as.array(parent),
               n_strat = n_strat, lgr_spill_std = as.numeric(lgr_std),
               lgs_spill_std = as.numeric(lgs_std)), phi_prior)

  delta_block <- if (delta_mode == "free") "
  delta0 ~ dnorm(0, pow(delta_sd, -2))
  sigma_delta ~ dunif(0, 3)
  tau_delta <- pow(sigma_delta + 1.0E-6, -2)
  for (s in 1:S) { delta[s] ~ dnorm(delta0, tau_delta) }
" else "
  for (s in 1:S) { delta[s] <- 0 }
"
  if (delta_mode == "free") jd$delta_sd <- delta_sd

  model_string <- paste0("model {
  eps <- 1.0E-9

  # Guidance efficiency: nested parent-group and weekly variation.
  tau_psi <- pow(sigma_psi, -2)
  tau_strat <- pow(sigma_strat, -2)
  for (s in 1:S) {
    logit_psi[s] ~ dnorm(mu_strat[parent[s]] + beta * lgr_spill_std[s], tau_psi)
    psi[s] <- ilogit(logit_psi[s])
  }
  for (g in 1:n_strat) { mu_strat[g] ~ dnorm(alpha, tau_strat) }
  alpha ~ dt(0, pow(2, -2), 7)
  beta ~ dt(0, pow(1, -2), 7) T(, 0)
  sigma_psi ~ dnorm(0, pow(0.5, -2)) T(0,)
  sigma_strat ~ dunif(0.05, 3)

  # Spillway-array detection: the same nesting, separate parameters.
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

  # Downstream recovery anchored on the SPILLWAY route.
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
  for (s in 1:S) {
    trans[s] ~ dbeta(1, 1)
    pi_raw[s, 1] <- psi[s] * (1 - trans[s]) * (1 - phi_B[s]) + eps
    pi_raw[s, 2] <- psi[s] * (1 - trans[s]) * phi_B[s] + eps
    pi_raw[s, 3] <- (1 - psi[s]) * p[s] * (1 - phi_S[s]) + eps
    pi_raw[s, 4] <- (1 - psi[s]) * p[s] * phi_S[s] + eps
    pi_raw[s, 5] <- (1 - psi[s]) * (1 - p[s]) * phi_S[s] + eps
    pi_raw[s, 6] <- psi[s] * trans[s] + eps
    p_seen[s] <- sum(pi_raw[s, 1:6])
    for (c in 1:6) { pi_obs[s, c] <- pi_raw[s, c] / p_seen[s] }
  }
  for (j in 1:S_lik) {
    n[lik_idx[j], 1:6] ~ dmulti(pi_obs[lik_idx[j], 1:6], N_seen[lik_idx[j]])
  }
}")

  for (pkg in c("rjags", "coda"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Install R package '", pkg,
           "' and ensure the separate JAGS application is installed.", call. = FALSE)
  monitor <- c("psi", "p", "phi_S", "phi_B", "trans",
                "alpha", "beta", "sigma_psi", "mu_strat", "sigma_strat",
                "alpha_p", "beta_p", "sigma_p", "mu_strat_p", "sigma_strat_p",
                "alpha_phi", "beta_phi", "sigma_phi")
  if (delta_mode == "free") monitor <- c(monitor, "delta", "delta0", "sigma_delta")
  settings <- list(model = "six-cell spillway-anchored", delta_mode = delta_mode,
                    delta_sd = delta_sd, parent_floor = parent_floor,
                    n_iter = n_iter, n_adapt = n_adapt, n_burnin = n_burnin,
                    n_chains = n_chains, n_thin = n_thin, seed = seed,
                    rhat_threshold = rhat_threshold, phi_prior_source = prior_source,
                    likelihood_count_cap = Inf)
  if (verbose) {
    message("Six-cell fit: ", S, " rows, ", n_strat, " parents, delta ", delta_mode,
            "; ", n_chains, " chains x ", n_iter,
            " post-burn iterations, thin ", n_thin, ".")
    message("Phi-prior source: ", prior_source,
            ". Review these settings before production use.")
  }
  con <- textConnection(model_string)
  on.exit(close(con), add = TRUE)
  inits <- lapply(seq_len(n_chains), function(i)
    list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = as.integer(seed + i)))
  jm <- rjags::jags.model(con, data = jd, inits = inits, n.chains = n_chains,
                          n.adapt = n_adapt, quiet = !verbose)
  if (n_burnin > 0L)
    stats::update(jm, n.iter = n_burnin, progress.bar = if (verbose) "text" else "none")
  samples <- rjags::coda.samples(jm, variable.names = monitor,
                                n.iter = n_iter, thin = n_thin,
                                progress.bar = if (verbose) "text" else "none")

  # The existing generator takes psi columns in their stored order. JAGS names
  # must map to row 1,...,S, not a lexical order such as 1,10,11,...,2.
  psi_names <- paste0("psi[", seq_len(S), "]")
  sample_names <- colnames(as.matrix(samples[[1L]]))
  if (!all(psi_names %in% sample_names))
    stop("JAGS did not return all expected indexed psi columns.", call. = FALSE)
  column_order <- c(psi_names, setdiff(sample_names, psi_names))
  samples <- coda::mcmc.list(lapply(samples, function(ch) ch[, column_order, drop = FALSE]))
  mat <- do.call(rbind, lapply(samples, as.matrix))
  psrf <- tryCatch(coda::gelman.diag(samples, autoburnin = FALSE,
                                    multivariate = FALSE)$psrf,
                    error = function(e) {
                      warning("R-hat unavailable: ", conditionMessage(e), call. = FALSE)
                      matrix(NA_real_, nrow = ncol(mat), ncol = 2,
                             dimnames = list(colnames(mat), NULL))
                    })
  ess <- tryCatch(coda::effectiveSize(samples), error = function(e) {
    warning("ESS unavailable: ", conditionMessage(e), call. = FALSE)
    stats::setNames(rep(NA_real_, ncol(mat)), colnames(mat))
  })
  quant <- t(apply(mat, 2, stats::quantile, probs = c(0.025, 0.5, 0.975)))
  sm <- data.frame(parameter = colnames(mat), mean = colMeans(mat),
                    sd = apply(mat, 2, stats::sd), median = quant[, 2],
                    lci = quant[, 1], uci = quant[, 3],
                    Rhat = psrf[colnames(mat), 1], Rhat_upper = psrf[colnames(mat), 2],
                    n_eff = ess[colnames(mat)], row.names = NULL)
  psi_summary <- cbind(gd[, c("stratum_idx", "stratum", "parent", "N_seen")],
                        sm[match(psi_names, sm$parameter),
                           c("mean", "median", "lci", "uci", "Rhat", "n_eff")])
  psi_summary$width <- psi_summary$uci - psi_summary$lci
  if (!is.null(weeks)) psi_summary$week <- weeks
  bad <- !is.finite(sm$Rhat) | sm$Rhat > rhat_threshold
  bad_psi <- !is.finite(psi_summary$Rhat) | psi_summary$Rhat > rhat_threshold
  diagnostics <- list(method = "classical coda R-hat, autoburnin=FALSE; coda ESS",
                       threshold = rhat_threshold, all_rhat_pass = !any(bad),
                       psi_rhat_pass = !any(bad_psi),
                       flagged_parameters = sm[bad, , drop = FALSE])
  if (any(bad))
    warning(sum(bad), " monitored parameter(s), including ", sum(bad_psi),
            " psi value(s), have R-hat > ", rhat_threshold, " or unavailable R-hat.",
            " Inspect diagnostics and traces before SCRAPI2; 30k is not a guarantee.",
            call. = FALSE)
  if (verbose) print(psi_summary, row.names = FALSE, digits = 4)

  # Convenience aliases for the original psi_table()/simulation code.
  sims_list <- rhat_list <- neff_list <- stats::setNames(vector("list", length(monitor)), monitor)
  for (v in monitor) {
    cols <- grep(paste0("^", v, "(\\[[0-9]+\\])?$"), colnames(mat), value = TRUE)
    if (length(cols) > 1L) {
      index <- as.integer(sub(".*\\[([0-9]+)\\]$", "\\1", cols))
      cols <- cols[order(index)]
    }
    sims_list[[v]] <- if (length(cols) == 1L && cols == v) mat[, cols] else mat[, cols, drop = FALSE]
    rhat_list[[v]] <- unname(sm$Rhat[match(cols, sm$parameter)])
    neff_list[[v]] <- unname(sm$n_eff[match(cols, sm$parameter)])
  }
  list(samples = samples, ge_data = gd, obs_strata = gd,
        likelihood_strata = gd[gd$has_history, , drop = FALSE],
        spill_mean = lgr_mean, spill_sd = lgr_sd,
        lgs_spill_mean = lgs_mean, lgs_spill_sd = lgs_sd,
        phi_prior = phi_prior, nested = TRUE, nested_p = TRUE,
        route_effect = delta_mode == "free", delta_mode = delta_mode,
        weeks = weeks,
        strat_assign = if (is.null(weeks)) NULL else data.frame(Week = weeks, Collapse = labels),
        psi_summary = psi_summary, summary = sm, diagnostics = diagnostics,
        sims.list = sims_list, Rhat = rhat_list, n.eff = neff_list,
        model_string = model_string, jags_data = jd, settings = settings,
        package_versions = c(rjags = as.character(utils::packageVersion("rjags")),
                             coda = as.character(utils::packageVersion("coda"))))
}
