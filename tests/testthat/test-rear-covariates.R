test_that("daily projection centers on each rear's fitted covariates", {
  skip_if_not_installed("coda")
  m <- cbind(`psi[1,1]` = rep(.2, 4), beta = rep(-1, 4),
             beta_outflow = rep(.5, 4))
  fit <- list(samples = coda::mcmc.list(coda::mcmc(m)),
              rear_levels = "W", target_rear = "W", weeks = 16L,
              spill_mean = 70, spill_sd = 20, lgr_spill_std = 0,
              outflow_mean = 80, outflow_sd = 15, outflow_std = 0,
              psi_spill_std = matrix(1, 1, 1), psi_outflow_std = matrix(1, 1, 1),
              N_seen = matrix(10, 1, 1),
              strat_assign = data.frame(Week = 16L, Collapse = 1L))
  date <- as.Date("2025-04-14")
  daily <- data.frame(Date = date, spill.per = 90, outflow = 95)
  out <- generate_rear_ge_draws(fit, pass_dates = date, B = 4,
                                daily_spill = daily, seed = 1)
  expect_equal(as.numeric(out[1, -1]), rep(.2, 4), tolerance = 1e-12)
  # Older fits still use their shared reference.
  fit$psi_spill_std <- fit$psi_outflow_std <- NULL
  old <- generate_rear_ge_draws(fit, pass_dates = date, B = 4,
                                daily_spill = daily, seed = 1)
  expect_equal(as.numeric(old[1, -1]), rep(plogis(qlogis(.2) - .5), 4), tolerance = 1e-12)
  centered <- generate_rear_ge_draws(fit, pass_dates = date, B = 4,
                                     daily_spill = daily, seed = 1, centering = TRUE,
                                     center_weights = data.frame(Date = date, weight = 10))
  expect_equal(as.numeric(centered[1, -1]), rep(.2, 4), tolerance = 1e-12)
  expect_equal(attr(centered, "center_reference")$total_weight, 10)
})
