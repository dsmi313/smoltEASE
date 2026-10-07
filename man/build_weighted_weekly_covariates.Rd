% Generated from R/rear_ge_covariates.R; edit that file.
\name{build_weighted_weekly_covariates}
\alias{build_weighted_weekly_covariates}
\title{Build rear-weighted weekly GE covariate matrices}
\usage{
build_weighted_weekly_covariates(daily_weights, daily_spill, weeks,
                                            rear_levels, shared_spill, shared_outflow)
}
\arguments{
\item{daily_weights}{Result of \code{\link{calc_daily_passage_weights}()}, or a data frame with rear_type, date, week, weight; unique rear-date rows.}
\item{daily_spill}{Daily covariates, possibly from \code{\link{apply_spill_decay}()}.}
\item{weeks}{Unique fitted ISO weeks, in model order.}
\item{rear_levels}{Unique rear groups, in model order.}
\item{shared_spill,shared_outflow}{Finite shared weekly reference vectors, one per week. Spill must be in [0,100] and outflow nonnegative.}
}
\value{List with psi_spill, psi_outflow, and audit (rear, week, total_weight, n_days, fallback, spill.per, outflow).}
\description{Uses daily PIT passage weights to average the supplied daily spill and
outflow separately by rear and ISO week. Does not apply spill decay itself.
Rear-weeks with no positive weight keep their shared weekly reference;
the audit identifies every such fallback. Matrices are in exactly the
supplied rear_levels and weeks order, ready for psi_spill and psi_outflow.}
