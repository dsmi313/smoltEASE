% Generated from R/rear_ge_covariates.R; edit that file.
\name{prep_rear_ge_covariates}
\alias{prep_rear_ge_covariates}
\title{Prepare passage-weighted GE covariates with configurable spill decay}
\usage{
prep_rear_ge_covariates(events, daily_spill, ge_fit, decay = 0.5,
                                    lag_days = 7L, tz = "UTC")
}
\arguments{
\item{events}{Cleaned PIT passage records; see \code{\link{calc_daily_passage_weights}()}.}
\item{daily_spill}{Unlagged daily Date, spill.per, outflow.}
\item{ge_fit}{Saved reference rear fit with samples, rear_levels, weeks, lgr_spill_pct and the shared weekly outflow vector.}
\item{decay}{Spill-decay factor in [0,1]. Default 0.5.}
\item{lag_days}{Number of preceding calendar days. Default 7.}
\item{tz}{Clock-time zone for passage dates. Default "UTC".}
}
\value{List with psi_spill, psi_outflow, audit, daily_weights, p_week, daily_spill, and settings. The matrices can be passed directly to \code{\link{fit_ge_rear_daynight}()} or \code{\link{fit_ge_rear_model2}()}.}
\description{Combines \code{\link{calc_daily_passage_weights}()}, \code{\link{apply_spill_decay}()} and
\code{\link{build_weighted_weekly_covariates}()}. Extracts posterior mean weekly p
from a saved reference rear fit, normally the production fit. The reference
fit supplies rear/week order and fallback covariates, not new fitted models.
Use the returned daily_spill both for weighting and for daily GE generation.
Changing decay requires a corresponding refit; changing daily covariates
alone on an old fit does not perform a new decay calibration.}
