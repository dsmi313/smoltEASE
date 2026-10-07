% Generated from R/rear_ge_covariates.R; edit that file.
\name{apply_spill_decay}
\alias{apply_spill_decay}
\title{Apply a configurable spill-decay window}
\usage{
apply_spill_decay(daily_spill, decay = 0.5, lag_days = 7L)
}
\arguments{
\item{daily_spill}{Data frame with Date, spill.per (percent, 0 to 100), outflow (nonnegative daily mean flow rate). One row per consecutive date.}
\item{decay}{Spill-decay factor between zero and one. Default 0.5.}
\item{lag_days}{Number of preceding calendar days; default 7.}
}
\value{Sorted daily_spill with spill.per replaced by the decayed spill. Other columns are retained. Attribute spill_decay records the settings.}
\description{For date d, computes sum(decay^k * spill[d-k]) / sum(decay^k), with
k = 0,...,lag_days. At the beginning of the supplied series only available
preceding days are used. Provide earlier dates for a full initial window.
decay = 0 gives same-day spill; decay = 1 gives an equally weighted window.
Calendar dates must be consecutive. Outflow is always left at its daily value.}
