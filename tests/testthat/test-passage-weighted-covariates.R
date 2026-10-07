passage_weight_fixture <- function() {
  data.frame(tag = c("g", "g", "s", "s", "both", "both", "transport", "load", "h"),
    rear_type = c(rep("W", 8), "H"), source = c(rep("main", 7), "transport", "main"),
    site = c("GRJ", "GRJ", "GRS", "GRS", "GRS", "GRJ", "GRJ", "GRJ", "GRS"),
    det_time = as.POSIXct(c("2025-04-14 12:00:00", "2025-04-15 12:00:00",
      "2025-04-14 12:00:00", "2025-04-15 12:00:00", "2025-04-15 01:00:00",
      "2025-04-15 02:00:00", "2025-04-14 12:00:00", "2025-04-17 12:00:00",
      "2025-04-14 12:00:00"), tz = "UTC"))
}
weekly_p_fixture <- function() data.frame(rear_type = c("W", "H"), week = 16L, p_mean = c(.25, .5))
daily_cov_fixture <- function() data.frame(Date = as.Date("2025-04-14") + 0:7,
                                           spill.per = seq(10, 80, 10), outflow = seq(100, 170, 10))

test_that("passage weights deduplicate fish, prioritize GRJ and ignore loading records", {
  w <- calc_daily_passage_weights(passage_weight_fixture(), weekly_p_fixture(), weeks = 16:17)
  expect_equal(sum(w$G), 3L)
  expect_equal(sum(w$S), 2L)
  expect_equal(w$weight[w$rear_type == "W"], c(6, 1))
  expect_equal(w$weight[w$rear_type == "H"], 2)
  expect_equal(as.character(w$date[w$rear_type == "W"]), c("2025-04-14", "2025-04-15"))
})

test_that("decay endpoints and initial windows use calendar lags and preserve flow", {
  ds <- daily_cov_fixture()[1:3, ]
  ds$spill.per <- c(10, 40, 70)
  expect_equal(apply_spill_decay(ds, decay = 0)$spill.per, ds$spill.per)
  expect_equal(apply_spill_decay(ds, decay = .5, lag_days = 2)$spill.per,
                c(10, 30, (70 + .5 * 40 + .25 * 10) / 1.75))
  expect_equal(apply_spill_decay(ds, decay = 1, lag_days = 1)$spill.per, c(10, 25, 55))
  expect_equal(apply_spill_decay(ds, decay = .8, lag_days = 0)$spill.per, ds$spill.per)
  expect_equal(apply_spill_decay(ds, decay = .5)$outflow, ds$outflow)
  reversed <- ds[3:1, ]
  expect_equal(apply_spill_decay(reversed, .5)$Date, ds$Date)
})

test_that("weekly matrices honor rear/week order and audit empty-week fallbacks", {
  w <- calc_daily_passage_weights(passage_weight_fixture(), weekly_p_fixture(), weeks = 16:17)
  ds <- daily_cov_fixture()
  out <- build_weighted_weekly_covariates(w, ds, weeks = c(17L, 16L),
    rear_levels = c("W", "H"), shared_spill = c(80, 20), shared_outflow = c(170, 110))
  expect_identical(dimnames(out$psi_spill), list(c("W", "H"), c("17", "16")))
  expect_equal(unname(out$psi_spill), rbind(c(80, 80 / 7), c(80, 10)))
  expect_equal(unname(out$psi_outflow), rbind(c(170, 710 / 7), c(170, 100)))
  expect_equal(out$audit$total_weight, c(0, 7, 0, 2))
  expect_equal(out$audit$fallback, c(TRUE, FALSE, TRUE, FALSE))
  w$weight <- 0
  zero <- build_weighted_weekly_covariates(w, ds, 16L, c("W", "H"), 30, 120)
  expect_true(all(zero$audit$fallback))
  expect_equal(unname(zero$psi_spill), matrix(30, 2, 1))
})

test_that("wrapper matches the established script calculations at configurable decay", {
  ds <- daily_cov_fixture()
  fit <- list(rear_levels = c("H", "W"), weeks = 16:17,
              lgr_spill_pct = c(40, 80), outflow = c(130, 170),
              samples = list(cbind(`p[1,1]` = c(.4, .6), `p[1,2]` = c(.5, .5),
                                   `p[2,1]` = c(.2, .3), `p[2,2]` = c(.25, .25))))
  for (decay in c(0, .5, .8)) {
    out <- prep_rear_ge_covariates(passage_weight_fixture(), ds, fit, decay = decay, lag_days = 7)
    # Independent reproduction of the old run-script lag calculation.
    expected <- vapply(seq_len(nrow(ds)), function(i) {
      k <- 0:min(7, i - 1)
      sum(decay^k * ds$spill.per[i - k]) / sum(decay^k)
    }, numeric(1))
    expect_equal(out$daily_spill$spill.per, expected)
    expect_equal(out$psi_spill["W", "16"], (6 * expected[1] + expected[2]) / 7)
    expect_equal(out$psi_outflow["W", "16"], 710 / 7)
    expect_equal(out$psi_spill["W", "17"], 80) # no PIT weight: shared fallback
    expect_equal(out$settings$decay, decay)
    expect_false(out$settings$outflow_lagged)
    expect_equal(out$p_week$p_mean, c(.5, .5, .25, .25))
  }
})

test_that("invalid inputs fail instead of silently changing the weighting", {
  ev <- passage_weight_fixture(); p <- weekly_p_fixture(); ds <- daily_cov_fixture()
  expect_error(calc_daily_passage_weights(ev, rbind(p, p)), "unique")
  p$p_mean[1] <- 0
  expect_error(calc_daily_passage_weights(ev, p), "p_mean")
  expect_error(calc_daily_passage_weights(ev, weekly_p_fixture()[1, ]), "lack p_mean")
  ev$rear_type[2] <- "H"
  expect_error(calc_daily_passage_weights(ev, weekly_p_fixture()), "more than one rear")
  expect_error(apply_spill_decay(ds, -1), "decay")
  expect_error(apply_spill_decay(ds, .5, lag_days = 1.5), "lag_days")
  expect_error(apply_spill_decay(ds[-2, ]), "consecutive")
  expect_error(apply_spill_decay(rbind(ds, ds[1, ])), "duplicate dates")
  ds$spill.per[1] <- NA_real_
  expect_error(apply_spill_decay(ds), "Daily spill")
  w <- calc_daily_passage_weights(passage_weight_fixture(), weekly_p_fixture())
  expect_error(build_weighted_weekly_covariates(w, daily_cov_fixture()[-1, ], 16:17,
                 c("H", "W"), c(40, 80), c(130, 170)), "lack daily")
  w$week[1] <- 17L
  expect_error(build_weighted_weekly_covariates(w, daily_cov_fixture(), 16:17,
                 c("H", "W"), c(40, 80), c(130, 170)), "ISO")
})

test_that("empty passage records keep finite shared covariates", {
  ev <- passage_weight_fixture()[FALSE, ]
  w <- calc_daily_passage_weights(ev, weekly_p_fixture())
  expect_equal(nrow(w), 0L)
  out <- build_weighted_weekly_covariates(w, daily_cov_fixture(), 16:17,
                                         c("H", "W"), c(40, 80), c(130, 170))
  expect_true(all(out$audit$fallback))
})
