test_that("rear-specific GE draws preserve the requested rear", {
  skip_if_not_installed("coda")
  m <- cbind(
    `psi[1,1]` = c(.7, .71, .72, .73),
    `psi[1,2]` = c(.6, .61, .62, .63),
    `psi[2,1]` = c(.2, .21, .22, .23),
    `psi[2,2]` = c(.1, .11, .12, .13),
    beta = rep(-.5, 4), beta_outflow = rep(.2, 4),
    beta_interaction = rep(.1, 4)
  )
  fit <- list(
    samples = coda::mcmc.list(coda::mcmc(m), coda::mcmc(m)),
    rear_levels = c("H", "W"), target_rear = "W", weeks = 13:14,
    spill_mean = 80, spill_sd = 10, lgr_spill_std = c(-1, 1),
    outflow_mean = 90, outflow_sd = 10, outflow_std = c(-1, 1),
    interaction_std = c(1, 1),
    N_seen = matrix(c(10, 10, 5, 5), nrow = 2, byrow = TRUE),
    strat_assign = data.frame(Week = 13:14, Collapse = 1:2)
  )
  dates <- as.Date(c("2025-03-26", "2025-04-02"))
  out <- generate_rear_ge_draws(fit, rear_type = "W", pass_dates = dates,
                                B = 20, seed = 4)
  expect_equal(nrow(out), 2)
  expect_true(all(unlist(out[1, -1]) %in% m[, "psi[2,1]"]))
  expect_true(all(unlist(out[2, -1]) %in% m[, "psi[2,2]"]))
})

test_that("daily spill equal to weekly spill leaves GE unchanged", {
  skip_if_not_installed("coda")
  m <- cbind(`psi[1,1]` = c(.2, .25, .3, .35), beta = rep(-.7, 4),
             beta_outflow = rep(.3, 4), beta_interaction = rep(.2, 4))
  fit <- list(
    samples = coda::mcmc.list(coda::mcmc(m), coda::mcmc(m)),
    rear_levels = "W", target_rear = "W", weeks = 13,
    spill_mean = 80, spill_sd = 10, lgr_spill_std = 0,
    outflow_mean = 90, outflow_sd = 10, outflow_std = 0,
    interaction_std = 0,
    N_seen = matrix(10, nrow = 1),
    strat_assign = data.frame(Week = 13, Collapse = 1)
  )
  date <- as.Date("2025-03-26")
  a <- generate_rear_ge_draws(fit, pass_dates = date, B = 20, seed = 9)
  b <- generate_rear_ge_draws(
    fit, pass_dates = date, B = 20,
    daily_spill = data.frame(Date = date, spill.per = 80, outflow = 90),
    seed = 9)
  expect_equal(a, b, tolerance = 1e-12)
})

test_that("missing daily spill stops by default", {
  skip_if_not_installed("coda")
  m <- cbind(`psi[1,1]` = rep(.2, 4), beta = rep(-.7, 4),
             beta_outflow = rep(.3, 4), beta_interaction = rep(.2, 4))
  fit <- list(
    samples = coda::mcmc.list(coda::mcmc(m), coda::mcmc(m)),
    rear_levels = "W", target_rear = "W", weeks = 13,
    spill_mean = 80, spill_sd = 10, lgr_spill_std = 0,
    outflow_mean = 90, outflow_sd = 10, outflow_std = 0,
    interaction_std = 0,
    N_seen = matrix(10, nrow = 1),
    strat_assign = data.frame(Week = 13, Collapse = 1)
  )
  expect_error(
    generate_rear_ge_draws(
      fit, pass_dates = as.Date("2025-03-26"), B = 4,
      daily_spill = data.frame(Date = as.Date("2025-03-27"),
                               spill.per = 80, outflow = 90)),
    "lack finite daily spill or outflow"
  )
})

test_that("daily projection carries the spill-by-outflow interaction", {
  skip_if_not_installed("coda")
  m <- cbind(`psi[1,1]` = rep(.2, 4), beta = 0,
             beta_outflow = 0, beta_interaction = 1)
  fit <- list(
    samples = coda::mcmc.list(coda::mcmc(m), coda::mcmc(m)),
    rear_levels = "W", target_rear = "W", weeks = 13,
    spill_mean = 80, spill_sd = 10, lgr_spill_std = 0,
    outflow_mean = 90, outflow_sd = 10, outflow_std = 0,
    interaction_std = 0, N_seen = matrix(10, nrow = 1),
    strat_assign = data.frame(Week = 13, Collapse = 1)
  )
  date <- as.Date("2025-03-26")
  out <- generate_rear_ge_draws(
    fit, pass_dates = date, B = 4,
    daily_spill = data.frame(Date = date, spill.per = 90, outflow = 100),
    seed = 1)
  expect_equal(as.numeric(out[1, -1]), rep(plogis(qlogis(.2) + 1), 4))
})
