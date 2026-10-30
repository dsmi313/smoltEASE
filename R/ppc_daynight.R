#' Posterior predictive check for a day/night GE fit
#'
#' Checks a [fit_ge_rear_daynight()] result against its own data, using the
#' saved posterior draws (no refit). For each draw it rebuilds the eleven cell
#' probabilities exactly as the JAGS model does, simulates a new set of
#' detection histories for every nonempty rear-week, and scores observed and
#' simulated counts against their expected values with the Freeman-Tukey
#' statistic (Kery and Schaub 2012, section 7.10.2). The Bayesian p-value is
#' the share of draws where the simulated data fit worse than the observed
#' data: near 0.5 is fine, near 0 means the observed counts sit further from
#' the model than its own simulated data usually do.
#'
#' The model has close to one free parameter per cell in each rear-week, so
#' it can match almost any counts and the overall check has little power. What
#' the model does restrict is that transport, bypass recovery, and (unless
#' `phi_night` was set) spill recovery are the same day and night within a
#' rear-week. `day_night` tests those: it compares the observed day minus
#' night difference in each proportion (season totals by rear) with the same
#' difference in the simulated data. A two-sided p near 0 means the data break
#' that assumption. Season cell totals (`by_cell`) show which history type is
#' over- or under-predicted.
#'
#' Covers the eleven-cell histories only. If the fit used `grs_detection`,
#' the row-history likelihood is not checked.
#'
#' @param fit A [fit_ge_rear_daynight()] result.
#' @param n_draws Number of posterior draws to use, sampled without
#'   replacement. `Inf` uses all of them.
#' @param seed Random seed for the draw subsample and the simulated data.
#' @return An object of class `ppc_daynight`: a list with `overall` (Bayesian
#'   p for all rear-weeks and by rear), `day_night`, `by_week` (Bayesian p per
#'   rear-week), `by_cell` (season totals by cell with 95% predictive
#'   intervals), the per-draw discrepancies `ft_obs` and `ft_new`, and
#'   `draws_used`. It has print and plot methods.
#' @examples
#' \dontrun{
#' fit <- readRDS("fit_MY2023_daynight_rear_decay0.50_raw.rds")
#' ppc <- ppc_daynight(fit)
#' ppc
#' plot(ppc)
#' }
#' @export
ppc_daynight <- function(fit, n_draws = 2000, seed = 1) {
  if (!inherits(fit, "ge_rear_daynight")) stop("fit must be a fit_ge_rear_daynight() result.", call. = FALSE)
  mat <- do.call(rbind, lapply(fit$samples, as.matrix))
  set.seed(seed)
  if (is.finite(n_draws) && nrow(mat) > n_draws) mat <- mat[sort(sample(nrow(mat), n_draws)), , drop = FALSE]
  D <- nrow(mat)
  R <- length(fit$rear_levels)
  S <- length(fit$weeks)
  eps <- as.numeric(regmatches(fit$model_string, regexpr("(?<=eps <- )[0-9.eE+-]+", fit$model_string, perl = TRUE)))
  if (length(eps) != 1L || !is.finite(eps)) stop("Could not read eps from fit$model_string.", call. = FALSE)
  if (isTRUE(fit$joint_grs_detection)) message("Joint GRS detection fit: only the eleven-cell histories are checked.")

  # draws x rear x week array for one saved parameter
  draws_of <- function(nm) {
    cols <- paste0(nm, "[", rep(seq_len(R), S), ",", rep(seq_len(S), each = R), "]")
    miss <- setdiff(cols, colnames(mat))
    if (length(miss)) stop("Parameter not saved in fit$samples: ", miss[1], call. = FALSE)
    array(mat[, cols], c(D, R, S))
  }
  a <- draws_of("a"); pd <- draws_of("psi_day"); pn <- draws_of("psi_night")
  tr <- draws_of("trans"); pB <- draws_of("phi_B"); pS <- draws_of("phi_S")
  p <- draws_of("p"); pnt <- draws_of("p_night")
  pSn <- if (!is.null(fit$phi_night) && fit$phi_night != "none") draws_of("phi_S_night") else pS

  cells <- c(paste0(c("c1", "c2", "c3", "c4", "c6"), "_day"),
             paste0(c("c1", "c2", "c3", "c4", "c6"), "_night"), "c5")
  lik <- which(fit$N_seen > 0, arr.ind = TRUE)
  J <- nrow(lik)
  ft_obs <- ft_new <- matrix(NA_real_, D, J)
  rep_tot <- array(0, c(D, R, 11), dimnames = list(NULL, fit$rear_levels, cells))
  obs_tot <- matrix(0, R, 11, dimnames = list(fit$rear_levels, cells))

  for (j in seq_len(J)) {
    r <- lik[j, 1]; s <- lik[j, 2]
    pr <- daynight_cell_probs(a[, r, s], pd[, r, s], pn[, r, s], tr[, r, s], pB[, r, s],
                              pS[, r, s], p[, r, s], pnt[, r, s], eps, pSn[, r, s])
    N <- fit$N_seen[r, s]
    n_obs <- fit$n[r, s, ]
    E <- N * pr
    n_new <- t(apply(pr, 1, function(q) stats::rmultinom(1, N, q)))
    ft_obs[, j] <- rowSums((sqrt(matrix(n_obs, D, 11, byrow = TRUE)) - sqrt(E))^2)
    ft_new[, j] <- rowSums((sqrt(n_new) - sqrt(E))^2)
    rep_tot[, r, ] <- rep_tot[, r, ] + n_new
    obs_tot[r, ] <- obs_tot[r, ] + n_obs
  }

  bayes_p <- function(cols) mean(rowSums(ft_new[, cols, drop = FALSE]) > rowSums(ft_obs[, cols, drop = FALSE]))
  overall <- data.frame(rear = c("all", fit$rear_levels), rear_weeks = c(J, tabulate(lik[, 1], R)),
                        bayes_p = c(bayes_p(seq_len(J)),
                                    vapply(seq_len(R), function(r) bayes_p(which(lik[, 1] == r)), 1)))
  by_week <- data.frame(rear = fit$rear_levels[lik[, 1]], week = fit$weeks[lik[, 2]],
                        N_seen = fit$N_seen[lik], bayes_p = vapply(seq_len(J), bayes_p, 1))
  by_week <- by_week[order(by_week$rear, by_week$week), ]
  by_cell <- do.call(rbind, lapply(seq_len(R), function(r) {
    q <- apply(rep_tot[, r, ], 2, stats::quantile, c(0.025, 0.5, 0.975))
    obs <- matrix(obs_tot[r, ], D, 11, byrow = TRUE)
    data.frame(rear = fit$rear_levels[r], cell = cells, observed = obs_tot[r, ],
               predicted_median = q[2, ], lower_95 = q[1, ], upper_95 = q[3, ],
               p_obs_above = colMeans(rep_tot[, r, ] < obs) + 0.5 * colMeans(rep_tot[, r, ] == obs),
               row.names = NULL)
  }))

  # Day vs night: proportions the model holds equal across periods
  props <- list(transport = list(num = "c6", den = c("c1", "c2", "c6")),
                bypass_recovery = list(num = "c2", den = c("c1", "c2")),
                spill_recovery = list(num = "c4", den = c("c3", "c4")))
  day_night <- do.call(rbind, lapply(names(props), function(nm) {
    do.call(rbind, lapply(seq_len(R), function(r) {
      prop_of <- function(tot, period) {
        num <- tot[, paste0(props[[nm]]$num, "_", period), drop = FALSE]
        den <- rowSums(tot[, paste0(props[[nm]]$den, "_", period), drop = FALSE])
        ifelse(den > 0, rowSums(num) / den, NA_real_)
      }
      obs <- matrix(obs_tot[r, ], 1, dimnames = list(NULL, cells))
      rep <- matrix(rep_tot[, r, ], D, dimnames = list(NULL, cells))
      d_obs <- prop_of(obs, "day") - prop_of(obs, "night")
      d_rep <- prop_of(rep, "day") - prop_of(rep, "night")
      ok <- is.finite(d_rep)
      p2 <- if (is.finite(d_obs) && any(ok)) {
        min(1, 2 * min(mean(d_rep[ok] >= d_obs), mean(d_rep[ok] <= d_obs)))
      } else NA_real_
      data.frame(rear = fit$rear_levels[r], proportion = nm,
                 day = prop_of(obs, "day"), night = prop_of(obs, "night"), observed_diff = d_obs,
                 simulated_diff_lower_95 = if (any(ok)) unname(stats::quantile(d_rep[ok], 0.025)) else NA_real_,
                 simulated_diff_upper_95 = if (any(ok)) unname(stats::quantile(d_rep[ok], 0.975)) else NA_real_,
                 p_two_sided = p2, row.names = NULL)
    }))
  }))
  structure(list(overall = overall, day_night = day_night, by_week = by_week, by_cell = by_cell,
                 ft_obs = rowSums(ft_obs), ft_new = rowSums(ft_new), draws_used = D),
            class = "ppc_daynight")
}

# Eleven cell probabilities for one rear-week, one row per draw, in the JAGS
# model's order: day c1 c2 c3 c4 c6, night c1 c2 c3 c4 c6, pooled c5. PSn is
# night spill recovery (equal to PS unless the fit set phi_night).
daynight_cell_probs <- function(A, Pd, Pn, Tr, PB, PS, P, Pnt, eps, PSn = PS) {
  pr <- cbind((1 - A) * Pd * (1 - Tr) * (1 - PB),
              (1 - A) * Pd * (1 - Tr) * PB,
              (1 - A) * (1 - Pd) * P * (1 - PS),
              (1 - A) * (1 - Pd) * P * PS,
              (1 - A) * Pd * Tr,
              A * Pn * (1 - Tr) * (1 - PB),
              A * Pn * (1 - Tr) * PB,
              A * (1 - Pn) * Pnt * (1 - PSn),
              A * (1 - Pn) * Pnt * PSn,
              A * Pn * Tr,
              (1 - A) * (1 - Pd) * (1 - P) * PS + A * (1 - Pn) * (1 - Pnt) * PSn) + eps
  pr / rowSums(pr)
}

#' @export
print.ppc_daynight <- function(x, ...) {
  cat("Posterior predictive check,", x$draws_used, "draws\n\n")
  print(format(x$overall, digits = 3), row.names = FALSE)
  cat("\nDay vs night proportions the model holds equal (season totals by rear)\n")
  print(format(x$day_night, digits = 3), row.names = FALSE)
  flagged <- x$by_week[x$by_week$bayes_p < 0.05 | x$by_week$bayes_p > 0.95, ]
  cat("\nRear-weeks with p < 0.05 or > 0.95:", nrow(flagged), "of", nrow(x$by_week),
      "(about", round(0.1 * nrow(x$by_week), 1), "expected by chance)\n")
  if (nrow(flagged)) print(format(flagged, digits = 3), row.names = FALSE)
  off <- x$by_cell[x$by_cell$observed < x$by_cell$lower_95 | x$by_cell$observed > x$by_cell$upper_95, ]
  cat("\nSeason cell totals outside the 95% predictive interval:", nrow(off), "of", nrow(x$by_cell), "\n")
  if (nrow(off)) print(format(off, digits = 3), row.names = FALSE)
  invisible(x)
}

#' @export
plot.ppc_daynight <- function(x, main = "Posterior predictive check", ...) {
  lim <- range(c(x$ft_obs, x$ft_new))
  graphics::plot(x$ft_obs, x$ft_new, pch = 16, cex = 0.4, col = grDevices::adjustcolor("steelblue", 0.4),
                 xlim = lim, ylim = lim, xlab = "Discrepancy, observed counts",
                 ylab = "Discrepancy, simulated counts", main = main, ...)
  graphics::abline(0, 1, lty = 2)
  graphics::mtext(sprintf("Bayesian p = %.3f (points below the line: data fit worse than the model's own)",
                          x$overall$bayes_p[1]), side = 3, line = 0.2, cex = 0.8)
  invisible(x)
}
