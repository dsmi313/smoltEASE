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


test_that("hierarchical fit supports full and reduced rear structures", {
  fn <- paste(deparse(body(fit_ge_rear_model2)), collapse = "\n")
  expect_match(fn, "rear_structure", fixed = TRUE)
  expect_match(fn, "full_structure", fixed = TRUE)
  expect_match(fn, "p\\[r,s\\]", fixed = FALSE)
  expect_match(fn, "phi_S\\[r,s\\]", fixed = FALSE)
  expect_match(fn, "phi_S\\[s\\]", fixed = FALSE)
  expect_match(fn, "rear_p\\[r\\]", fixed = FALSE)
  expect_match(fn, "rear_phi\\[r\\]", fixed = FALSE)
  expect_false(grepl("beta_interaction ~", fn, fixed = TRUE))
})

test_that("rear-flexible recovery fits and retains rear-specific slopes", {
  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")
  n_H <- rbind(
    c(22, 3, 30, 2, 1, 0), c(18, 3, 28, 2, 1, 8),
    c(12, 2, 22, 2, 1, 12), c(15, 2, 25, 2, 1, 0))
  n_W <- rbind(
    c(8, 1, 18, 1, 1, 0), c(7, 1, 16, 1, 1, 4),
    c(5, 1, 13, 1, 1, 6), c(6, 1, 15, 1, 1, 0))
  make_group <- function(n) {
    list(n = n, lgr_spill_pct = c(25, 45, 65, 35),
         lgs_spill_pct = c(20, 40, 60, 30),
         lgr_outflow = c(70, 85, 100, 80),
         parent = 1:4)
  }
  ge_data <- list(H = make_group(n_H), W = make_group(n_W))
  trans_on <- matrix(c(0, 1, 1, 0, 0, 1, 1, 0),
                     nrow = 2, byrow = TRUE,
                     dimnames = list(c("H", "W"), as.character(13:16)))
  expect_error(
    fit_ge_rear_model2(ge_data, weeks = 13:16,
                       rear_structure = "reduced",
                       phi_structure = "rear_flexible"),
    "requires rear_structure"
  )
  bad_schedule <- trans_on
  bad_schedule["H", "14"] <- 0
  expect_error(
    fit_ge_rear_model2(ge_data, weeks = 13:16, trans_on = bad_schedule),
    "Transported \\(c6\\)"
  )
  fit <- suppressWarnings(fit_ge_rear_model2(
    ge_data, weeks = 13:16, rear_structure = "full",
    phi_structure = "rear_flexible", delta_mode = "fixed",
    trans_on = trans_on,
    n_adapt = 100, n_burnin = 100, n_iter = 200,
    n_chains = 2, n_thin = 2, verbose = FALSE
  ))
  expect_s3_class(fit, "ge_rear_model2")
  expect_equal(fit$settings$phi_structure, "rear_flexible")
  expect_true(all(c("beta_phi_rear[1]", "beta_phi_rear[2]",
                    "sigma_phi_rear[1]", "sigma_phi_rear[2]") %in%
                  fit$summary$parameter))
  expect_true(all(is.finite(fit$psi_summary$mean)))
  expect_equal(fit$settings$scheduled_transport, TRUE)
  off_nodes <- c("trans[1,1]", "trans[1,4]", "trans[2,1]", "trans[2,4]")
  expect_true(all(off_nodes %in% fit$summary$parameter))
  expect_true(all(fit$summary$mean[
    fit$summary$parameter %in% off_nodes] == 0))
})

test_that("new additive fits need no interaction column", {
  skip_if_not_installed("coda")
  m <- cbind(`psi[1,1]` = rep(.2, 4), beta = rep(-.5, 4),
             beta_outflow = rep(.25, 4))
  fit <- list(
    samples = coda::mcmc.list(coda::mcmc(m), coda::mcmc(m)),
    rear_levels = "W", target_rear = "W", weeks = 13,
    spill_mean = 80, spill_sd = 10, lgr_spill_std = 0,
    outflow_mean = 90, outflow_sd = 10, outflow_std = 0,
    N_seen = matrix(10, nrow = 1),
    strat_assign = data.frame(Week = 13, Collapse = 1)
  )
  date <- as.Date("2025-03-26")
  out <- generate_rear_ge_draws(
    fit, pass_dates = date, B = 4,
    daily_spill = data.frame(Date = date, spill.per = 90, outflow = 100),
    seed = 1)
  expected <- plogis(qlogis(.2) - .5 + .25)
  expect_equal(as.numeric(out[1, -1]), rep(expected, 4))
})
