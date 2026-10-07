# events: cleaned PIT histories, with transport dated to bypass passage.
# daily_spill: unlagged daily Date, spill.per, outflow.
# prod_fit: reference production rear fit used for weekly mean p and fallback covariates.
# Change only this value to prepare another spill-decay calibration.
SPILL_DECAY <- 0.5
covariates <- smoltEASE::prep_rear_ge_covariates(
  events, daily_spill, ge_fit = prod_fit,
  decay = SPILL_DECAY, lag_days = 7, tz = "UTC"
)

# Pass matrices to the EXISTING fitter (same ge_dn and fitting settings as before).
# Refit when changing SPILL_DECAY; use a distinct cache filename for each decay.
fit <- smoltEASE::fit_ge_rear_daynight(
  ge_dn, weeks = prod_fit$weeks, parent = prod_fit$parent,
  target_rear = "W", day_ge = "shared", p_night_offset = 0,
  psi_spill = covariates$psi_spill,
  psi_outflow = covariates$psi_outflow,
  trans_on = prod_fit$trans_on,
  delta_sd = 0.5, n_iter = 100000, n_adapt = 5000,
  n_burnin = 20000, n_chains = 4, n_thin = 10,
  seed = SEED, rhat_threshold = 1.01
)
# Retain the preparation settings alongside a cached fit.
attr(fit, "covariate_settings") <- covariates$settings

ge_W <- smoltEASE::generate_rear_ge_draws_daynight(
  fit, rear_type = "W", pass_dates = passage$SampleEndDate,
  B = 5000, daily_spill = covariates$daily_spill,
  spill_counts = spill_counts_W, shrink_k = 10, seed = SEED
)
# covariates$audit flags rear-weeks using shared references rather than PIT weights.
