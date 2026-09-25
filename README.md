# smoltEASE

Smolt escapement estimation with propagated uncertainty at Lower Granite Dam.

## Overview

`smoltEASE` estimates juvenile salmonid escapement past Lower Granite Dam (LGR) and attaches defensible uncertainty to those estimates.

Juvenile abundance at LGR has traditionally been estimated with SCRAPI, which propagates sampling uncertainty through a bootstrap but treats two important quantities as if they were known exactly: guidance efficiency (GE), the proportion of fish routed to the collection system, and genetic stock composition (GSI). Both are estimated from data and both carry real uncertainty. Treating them as fixed understates the uncertainty in the final abundance estimate.

`smoltEASE` implements SCRAPI2, which models GE and GSI explicitly and propagates their uncertainty into the escapement estimate alongside the existing bootstrap. The result is an abundance estimate whose interval reflects all three sources rather than one.

This brings juvenile estimation into alignment with the EASE framework used for adult returns at the same dam, so uncertainty is handled consistently across life stages.

## What it does

- Estimates smolt escapement and stock-specific abundance via `SCRAPI2`
- Fits a Bayesian hierarchical model for guidance efficiency as a function of spill conditions
- Generates GE posterior draws and propagates them through the escapement calculation
- Incorporates GSI uncertainty via posterior draws rather than point estimates
- Remains compatible with legacy `SCOBI::SCRAPI` inputs and workflows

## Installation

```r
remotes::install_github("dsmi313/smoltEASE",upgrade = F)
```

## Usage

Worked examples are in the package vignettes.

## Methods

The estimator has three components and three corresponding sources of uncertainty.

| Component | What it contributes | How uncertainty enters |
|---|---|---|
| Trap sampling and expansion | Daily abundance from sampled counts | Bootstrap resampling |
| Guidance efficiency | Expansion from collected fish to river-wide passage | Draws from the posterior of a Bayesian hierarchical model |
| Stock composition | Allocation of total abundance to genetic stocks | Draws from the GSI posterior |

Each GE and GSI draw is pushed through the full escapement calculation, so the spread of the resulting estimates reflects Monte Carlo integration over those parameters rather than a plug-in at their point estimates.

The composed interval combines bootstrap resampling of the trap sample with draws from the GE and GSI posteriors, so it is a hybrid rather than a pure confidence or credible interval. GE draws are propagated as complete joint MCMC rows restricted to the intersection of each stratum's central 95% marginal interval. This preserves dependence among weekly GE values and the shared spill parameters while limiting the influence of extreme posterior-tail GE on the inverse expansion, but it means the propagated draws are not the full joint posterior.

Simulation testing used 20 seeds, each drawing its own hyperparameters from the real-data posterior. The multistate model recovered GE with mean absolute bias 0.017 and 0.979 coverage, with coverage of 0.968, 0.957, and 0.979 for the GRS detection and the two route-specific downstream probabilities. Against an independently fixed 300,000-smolt target, fitted-GE SCRAPI2 had a median absolute error of 3.3% and all 20 nominal 90% intervals contained the target. The known-GE arm also covered 20 of 20, so the observed over-coverage is present at the sampling-noise floor rather than introduced by estimating GE. Twenty seeds cannot distinguish correct calibration from over-coverage, and GE and abundance are estimated sequentially rather than in a single joint likelihood.

## Rear-type partial pooling

The experimental six-cell rear model tests whether known wild fish can borrow information from hatchery and unknown-rear PIT histories without forcing all rear groups to share one GE curve.

```r
rear_fit <- fit_ge_rear_model2(
  ge_data = list(H = ge_H, W = ge_W, U = ge_U),
  weeks = 13:26,
  target_rear = "W",
  delta_mode = "free",
  delta_sd = 0.5,
  rear_sd_scale = 1,
  n_iter = 50000,
  n_adapt = 5000,
  n_burnin = 20000,
  n_chains = 3,
  n_thin = 10,
  seed = 11
)

ge_draws_W <- generate_rear_ge_draws(
  rear_fit,
  rear_type = "W",
  pass_dates = passageData$SampleEndDate,
  B = 5000,
  daily_spill = lgr_daily_covariates, # Date, spill.per, outflow
  clip_to_ci = FALSE,
  seed = 12
)
```

The model estimates rear-specific weekly GE and transport probabilities while sharing the seasonal hierarchy, covariate responses, spillway detection, downstream recovery, and route offset. The GE regression includes standardized percent spill, standardized mean outflow, and their interaction. The `U` level is an observed unknown-rear classification, not a latent hatchery/wild mixture. Treat this model as exploratory until posterior predictive checks and recovery simulations support its assumptions.

## Relationship to other tools

- `SCOBI::SCRAPI`: the legacy estimator. `smoltEASE` reads the same inputs and reduces to comparable results when GE and GSI uncertainty are switched off.
- `escapeLGD`: the adult-return counterpart at the same dam. Both implement the same general estimator structure at different life stages.

## Status

Under active development. Function names and arguments may change.

## Author

David R. Smith
Anadromous Research Biologist, Nampa Research Station
Idaho Department of Fish and Game
david.smith@idfg.idaho.gov

## License

MIT. See `LICENSE`.
