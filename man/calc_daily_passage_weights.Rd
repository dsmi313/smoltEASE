% Generated from R/rear_ge_covariates.R; edit that file.
\name{calc_daily_passage_weights}
\alias{calc_daily_passage_weights}
\title{Daily PIT passage weights by rear group}
\usage{
calc_daily_passage_weights(events, p_week, weeks = NULL, tz = "UTC")
}
\arguments{
\item{events}{Cleaned data frame with tag, rear_type, source, site, det_time. det_time must be POSIXct; these are passage times, not loading times.}
\item{p_week}{Data frame with rear_type, week, p_mean; unique rear-week rows.}
\item{weeks}{Optional ISO weeks to keep, in model order.}
\item{tz}{Time zone used to turn passage times into calendar dates. Default "UTC" matches PTAGIS local clock strings parsed with a UTC label.}
}
\value{Data frame with rear_type, date, week, G, S, p_mean, weight.}
\description{Counts each eligible tag once, on its first main GRJ/GRS passage date.
GRJ takes route priority for a tag detected at both sites. Transported fish
must already have their bypass passage records in source = "main"; loading
records are not counted. The caller must first exclude dam-tagged fish,
adult releases, reverse movements, and unusable transport histories.
Weights are G + S / p_mean, using a fixed posterior mean weekly GRS
detection probability. Uncertainty in these plug-in weights is not propagated.}
