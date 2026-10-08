test_that("classify_daynight splits at sunrise and sunset", {
  skip_if_not_installed("suncalc")
  t <- as.POSIXct(c("2025-04-15 12:00:00", "2025-04-15 02:00:00",
                   "2025-04-15 23:00:00", NA), tz = "UTC")
  expect_equal(classify_daynight(t), c("day", "night", "night", NA))
  t2 <- as.POSIXct(c("2025-04-15 19:30:00", "2025-04-15 19:45:00"), tz = "UTC")
  expect_equal(classify_daynight(t2), c("day", "night"))
  sun <- suncalc::getSunlightTimes(as.Date("2025-04-15"),
                                  lat = 46.660394, lon = -117.436261,
                                  keep = c("sunrise", "sunset"),
                                  tz = "America/Los_Angeles")
  boundaries <- c(sun$sunrise - 1, sun$sunrise, sun$sunset - 1, sun$sunset)
  clocks <- as.POSIXct(format(boundaries, "%Y-%m-%d %H:%M:%OS6"),
                      format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
  expect_equal(classify_daynight(clocks), c("night", "day", "day", "night"))
})

make_dn_fit <- function(pd, pn, a, p = 0.5, p_off = 0) {
  m <- cbind(`psi_day[1,1]` = pd, `psi_night[1,1]` = pn, `a[1,1]` = a,
             `p[1,1]` = rep(p, length(pd)),
             beta = rep(-1, length(pd)), beta_outflow = rep(0.5, length(pd)))
  fit <- list(samples = coda::mcmc.list(coda::mcmc(m)),
              rear_levels = "W", target_rear = "W", weeks = 16L,
              spill_mean = 70, spill_sd = 20, outflow_mean = 80, outflow_sd = 15,
              psi_spill_std = matrix(0, 1, 1), psi_outflow_std = matrix(0, 1, 1),
              p_night_offset = p_off)
  class(fit) <- c("ge_rear_daynight", "list")
  fit
}
spill_on_ref <- data.frame(Date = as.Date(c("2025-04-14", "2025-04-15")),
                           spill.per = 70, outflow = 80)

test_that("equal day and night GE give the same GE whatever the night share", {
  skip_if_not_installed("coda")
  fit <- make_dn_fit(pd = c(.1, .2), pn = c(.1, .2), a = c(.3, .3))
  counts <- data.frame(Date = spill_on_ref$Date, S_day = c(100, 0), S_night = c(0, 100))
  out <- generate_rear_ge_draws_daynight(fit, pass_dates = spill_on_ref$Date, B = 10,
                                        daily_spill = spill_on_ref, spill_counts = counts, seed = 1)
  g <- as.matrix(out[, -1])
  expect_true(all(g >= 0.1 - 1e-12 & g <= 0.2 + 1e-12))
  expect_equal(g[1, ], g[2, ])
})

test_that("a more nocturnal day gets higher GE, and attributes are returned", {
  skip_if_not_installed("coda")
  fit <- make_dn_fit(pd = rep(.01, 4), pn = rep(.3, 4), a = rep(.2, 4))
  counts <- data.frame(Date = spill_on_ref$Date, S_day = c(90, 10), S_night = c(10, 90))
  out <- generate_rear_ge_draws_daynight(fit, pass_dates = spill_on_ref$Date, B = 8,
                                        daily_spill = spill_on_ref, spill_counts = counts,
                                        shrink_k = 10, seed = 2)
  g <- as.matrix(out[, -1])
  expect_true(all(g[2, ] > g[1, ]))
  expect_equal(dim(attr(out, "night_fraction")), c(2L, 8L))
  expect_equal(dim(attr(out, "p_effective")), c(2L, 8L))
  expect_true(all(abs(attr(out, "p_effective") - 0.5) < 1e-12))
})

test_that("no spilled detections recover the weekly fish share at reference covariates", {
  skip_if_not_installed("coda")
  pd <- .02; pn <- .3; a <- .25
  fit <- make_dn_fit(pd = rep(pd, 2), pn = rep(pn, 2), a = rep(a, 2))
  counts <- data.frame(Date = as.Date("2025-04-14"), S_day = 0, S_night = 0)
  out <- generate_rear_ge_draws_daynight(fit, pass_dates = counts$Date, B = 2,
                                        daily_spill = spill_on_ref, spill_counts = counts, seed = 3)
  expect_equal(unname(unlist(out[1, -1])), rep(a * pn + (1 - a) * pd, 2), tolerance = 1e-10)
  expect_equal(unname(attr(out, "night_fraction")), matrix(a, 1, 2), tolerance = 1e-10)
})

test_that("no detections retain the weekly detected-spill share off reference", {
  skip_if_not_installed("coda")
  fit <- make_dn_fit(pd = rep(.02, 2), pn = rep(.3, 2), a = rep(.25, 2), p_off = .5)
  counts <- data.frame(Date = as.Date("2025-04-14"), S_day = 0, S_night = 0)
  cov <- spill_on_ref
  cov$spill.per <- 90
  out <- generate_rear_ge_draws_daynight(fit, pass_dates = counts$Date, B = 2,
                                        daily_spill = cov, spill_counts = counts, seed = 3)
  pd <- .5; pn <- plogis(qlogis(pd) + .5)
  sw <- .25 * (1 - .3) * pn / (.25 * (1 - .3) * pn + .75 * (1 - .02) * pd)
  gd <- plogis(qlogis(.02) - 1); gn <- plogis(qlogis(.3) - 1)
  f <- sw * (1 - gd) * pd / (sw * (1 - gd) * pd + (1 - sw) * (1 - gn) * pn)
  expect_equal(unname(attr(out, "night_fraction")), matrix(f, 1, 2), tolerance = 1e-10)
  expect_equal(unname(unlist(out[1, -1])), rep(f * gn + (1 - f) * gd, 2), tolerance = 1e-10)
})

test_that("both day GE structures compile and sample with scheduled transport", {
  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")
  make_group <- function() {
    day <- rbind(c(8, 2, 12, 2, 0, 0), c(7, 2, 10, 2, 0, 3),
                 c(6, 1, 8, 1, 0, 4), c(7, 1, 10, 1, 0, 0))
    night <- day
    night[, 1] <- night[, 1] + 2L
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
  dn <- list(H = make_group(), W = make_group())
  schedule <- matrix(rep(c(0, 1, 1, 0), 2), 2, 4, byrow = TRUE)
  for (mode in c("rear", "shared")) {
    fit <- suppressWarnings(fit_ge_rear_daynight(
      dn, day_ge = mode, p_night_offset = .5, trans_on = schedule,
      n_adapt = 100, n_burnin = 100, n_iter = 200,
      n_chains = 2, n_thin = 2, verbose = FALSE))
    expect_s3_class(fit, "ge_rear_daynight")
    expect_true(all(is.finite(fit$period_summary$psi_day)))
    expect_true(all(is.finite(fit$period_summary$psi_night)))
    expect_true(all(fit$summary$mean[fit$summary$parameter %in%
                    c("trans[1,1]", "trans[1,4]", "trans[2,1]", "trans[2,4]")] == 0))
  }
})

test_that("day/night preparation preserves six-cell counts and LGR clock times", {
  skip_if_not_installed("suncalc")
  events <- data.frame(
    tag = c("b", "br", "br", "s", "sr", "sr", "u", "t", "t"),
    site = c("GRJ", "GRJ", "GOJ", "GRS", "GRS", "GOJ", "GOJ", "GRJ", "GRJ"),
    det_time = as.POSIXct(c("2025-04-14 12:00:00", "2025-04-14 23:00:00",
                           "2025-04-15 23:00:00", "2025-04-14 12:00:00",
                           "2025-04-14 23:00:00", "2025-04-15 23:00:00",
                           "2025-04-15 12:00:00", "2025-04-14 12:00:00",
                           "2025-04-15 12:00:00"), tz = "UTC"),
    antenna = c(rep("01", 8), "61"), release_site = "UPSTREAM",
    source = c(rep("main", 8), "transport"))
  spill <- data.frame(Site = rep(c("LGR", "LGS"), each = 2),
                      DateTime = rep(c("2025-04-14 12:00:00", "2025-04-21 12:00:00"), 2),
                      HourlySpill = c(30, 50, 20, 40), HourlyFlow = c(80, 100, 70, 90))
  dn <- prep_ge_data_daynight(events, spill, weeks = 16:17)
  rebuilt <- dn$n_day + dn$n_night
  rebuilt[, 5] <- dn$n_c5
  expect_equal(rebuilt, dn$base$n)
  expect_s3_class(dn$tags$lgr_time, "POSIXct")
  expect_equal(nrow(dn$tags), 6L)
  expect_equal(dn$tags$period[match("br", dn$tags$tag)], "night")
  expect_true(is.na(dn$tags$period[match("u", dn$tags$tag)]))
  expect_equal(dn$tags$lgr_time[match("t", dn$tags$tag)],
               as.POSIXct("2025-04-14 12:00:00", tz = "UTC"))
})

test_that("daily_weights calibrates daily GE to the weekly estimate", {
  skip_if_not_installed("coda")
  pd <- c(.02, .03); pn <- c(.1, .12)
  fit <- make_dn_fit(pd = pd, pn = pn, a = c(.3, .3))
  days <- as.Date("2025-04-14") + 0:6
  cov <- data.frame(Date = days, spill.per = c(95, 90, 85, 70, 50, 30, 10), outflow = 80)
  counts <- data.frame(Date = days, S_day = 0, S_night = 0)
  wts <- data.frame(rear_type = "W", date = days, week = 16L, weight = c(5, 8, 3, 6, 2, 9, 4))
  out <- generate_rear_ge_draws_daynight(fit, pass_dates = days, B = 2, daily_spill = cov,
                                        spill_counts = counts, seed = 1, daily_weights = wts)
  off <- attr(out, "calibration_offset")
  id <- attr(out, "draw_id")
  zs <- (cov$spill.per - 70) / 20
  for (b in 1:2) {
    gd <- plogis(qlogis(pd[id[b]]) - zs + off$day[1, b])
    gn <- plogis(qlogis(pn[id[b]]) - zs + off$night[1, b])
    expect_equal(sum(wts$weight * gd) / sum(wts$weight), pd[id[b]], tolerance = 1e-9)
    expect_equal(sum(wts$weight * gn) / sum(wts$weight), pn[id[b]], tolerance = 1e-9)
  }
  expect_true(all(off$day < 0))
  plain <- generate_rear_ge_draws_daynight(fit, pass_dates = days, B = 2, daily_spill = cov,
                                          spill_counts = counts, seed = 1)
  expect_true(all(as.matrix(out[, -1]) < as.matrix(plain[, -1])))
  expect_null(attr(plain, "calibration_offset"))
})

test_that("calibration leaves GE unchanged at constant covariates", {
  skip_if_not_installed("coda")
  fit <- make_dn_fit(pd = rep(.02, 2), pn = rep(.3, 2), a = rep(.25, 2))
  counts <- data.frame(Date = spill_on_ref$Date, S_day = 0, S_night = 0)
  wts <- data.frame(rear_type = "W", date = spill_on_ref$Date, week = 16L, weight = c(1, 3))
  a <- generate_rear_ge_draws_daynight(fit, pass_dates = spill_on_ref$Date, B = 2,
                                      daily_spill = spill_on_ref, spill_counts = counts, seed = 3)
  b <- generate_rear_ge_draws_daynight(fit, pass_dates = spill_on_ref$Date, B = 2,
                                      daily_spill = spill_on_ref, spill_counts = counts, seed = 3,
                                      daily_weights = wts)
  expect_equal(as.matrix(a[, -1]), as.matrix(b[, -1]), tolerance = 1e-10)
})
