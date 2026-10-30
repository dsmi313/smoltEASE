make_group <- function(night_swap_c3c4 = FALSE) {
  day <- rbind(c(8, 2, 12, 2, 0, 0), c(7, 2, 10, 2, 0, 3),
               c(6, 1, 8, 1, 0, 4), c(7, 1, 10, 1, 0, 0))
  night <- day
  night[, 1] <- night[, 1] + 2L
  # Swapping c3 and c4 at night makes night spill recovery far above day's
  if (night_swap_c3c4) night[, 3:4] <- day[, 4:3]
  c5 <- rep(2L, 4)
  n <- day + night
  n[, 5] <- c5
  list(base = list(n = n, weeks = 13:16, parent = 1:4,
                   lgr_spill_pct = c(25, 45, 65, 35),
                   lgs_spill_pct = c(20, 40, 60, 30),
                   lgr_outflow = c(70, 85, 100, 80),
                   alpha_phi_mean = 0, alpha_phi_sd = 2,
                   beta_phi_mean = 0, beta_phi_sd = 1),
       n_day = day, n_night = night, n_c5 = c5)
}

test_that("phi_night = 'none' leaves the model text as before", {
  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")
  dn <- list(H = make_group(), W = make_group())
  for (st in c("original", "consistent")) {
    fit <- suppressWarnings(fit_ge_rear_daynight(
      dn, structure = st, n_adapt = 10, n_burnin = 0, n_iter = 10,
      n_chains = 1, n_thin = 1, verbose = FALSE))
    expect_identical(fit$phi_night, "none")
    expect_false(grepl("phi_S_night|phi_night", fit$model_string))
    expect_true(grepl("(1-p_night[r,s])) * phi_S[r,s] + eps", fit$model_string, fixed = TRUE))
    expect_false("phi_night_sd" %in% names(fit$jags_data))
  }
})

test_that("phi_night options compile, sample, and save the night recovery nodes", {
  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")
  dn <- list(H = make_group(), W = make_group())
  for (st in c("original", "consistent")) for (pn in c("shared", "rear")) {
    fit <- suppressWarnings(fit_ge_rear_daynight(
      dn, structure = st, phi_night = pn, day_ge = "shared",
      n_adapt = 50, n_burnin = 50, n_iter = 100, n_chains = 2, n_thin = 1, verbose = FALSE))
    pars <- fit$summary$parameter
    expect_identical(fit$phi_night, pn)
    expect_identical(fit$settings$phi_night, pn)
    expect_true(all(c("phi_S_night[1,1]", "phi_S_night[2,4]", "phi_night_eff[1]", "phi_night_eff[2]") %in% pars))
    expect_true(all(is.finite(fit$summary$mean[grepl("^phi_S_night\\[", pars)])))
    if (pn == "shared") {
      expect_equal(fit$summary$mean[pars == "phi_night_eff[1]"], fit$summary$mean[pars == "phi_night_eff[2]"])
    }
  }
  expect_error(fit_ge_rear_daynight(dn, phi_night_sd = 0), "phi_night_sd")
})

test_that("ppc cell probabilities match the JAGS model's pi_obs", {
  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")
  dn <- list(H = make_group(), W = make_group(night_swap_c3c4 = TRUE))
  for (pn in c("none", "rear")) {
    fit <- suppressWarnings(fit_ge_rear_daynight(
      dn, structure = "consistent", delta_prior = "loose", phi_night = pn, day_ge = "shared",
      p_night_offset = 0.5, n_adapt = 100, n_burnin = 100, n_iter = 100,
      n_chains = 1, n_thin = 1, verbose = FALSE))
    jm <- rjags::jags.model(textConnection(fit$model_string), data = fit$jags_data,
                            n.chains = 1, quiet = TRUE)
    stats::update(jm, 100, progress.bar = "none")
    nodes <- c("pi_obs", "a", "psi_day", "psi_night", "trans", "phi_B", "phi_S", "p", "p_night",
               if (pn != "none") "phi_S_night")
    sm <- as.matrix(rjags::coda.samples(jm, nodes, n.iter = 50, progress.bar = "none")[[1]])
    for (r in 1:2) for (s in 1:4) {
      g <- function(nm) sm[, sprintf("%s[%d,%d]", nm, r, s)]
      mine <- smoltEASE:::daynight_cell_probs(
        g("a"), g("psi_day"), g("psi_night"), g("trans"), g("phi_B"), g("phi_S"), g("p"), g("p_night"),
        1e-9, if (pn != "none") g("phi_S_night") else g("phi_S"))
      expect_lt(max(abs(mine - sm[, sprintf("pi_obs[%d,%d,%d]", r, s, 1:11)])), 1e-12)
    }
  }
})

test_that("ppc_daynight flags unequal day/night spill recovery and phi_night absorbs it", {
  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")
  dn <- list(H = make_group(), W = make_group(night_swap_c3c4 = TRUE))
  fit_it <- function(pn) suppressWarnings(fit_ge_rear_daynight(
    dn, structure = "consistent", delta_prior = "loose", phi_night = pn, day_ge = "shared",
    n_adapt = 300, n_burnin = 500, n_iter = 2000, n_chains = 2, n_thin = 2, verbose = FALSE))
  ppc0 <- ppc_daynight(fit_it("none"), n_draws = 1000)
  expect_s3_class(ppc0, "ppc_daynight")
  expect_named(ppc0, c("overall", "day_night", "by_week", "by_cell", "ft_obs", "ft_new", "draws_used"))
  expect_equal(nrow(ppc0$by_cell), 22L)
  spill0 <- ppc0$day_night[ppc0$day_night$proportion == "spill_recovery", ]
  expect_lt(spill0$p_two_sided[spill0$rear == "W"], 0.05)
  expect_output(print(ppc0), "Day vs night")
  fit_r <- fit_it("rear")
  eff <- fit_r$summary$mean[fit_r$summary$parameter == "phi_night_eff[2]"]
  expect_gt(eff, 0)
  ppc_r <- ppc_daynight(fit_r, n_draws = 1000)
  spill_r <- ppc_r$day_night[ppc_r$day_night$proportion == "spill_recovery", ]
  expect_gt(spill_r$p_two_sided[spill_r$rear == "W"], 0.05)
  expect_error(ppc_daynight(list()), "fit_ge_rear_daynight")
})
