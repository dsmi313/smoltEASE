# smoltEASE

**Smolt Escapement Analysis with Uncertainty at Lower Granite Dam**

---

## Overview

`smoltEASE` provides tools for estimating juvenile salmonid escapement at Lower Granite Dam (LGR) using extensions of the SCRAPI framework.

The package implements **SCRAPI2**, which improves on traditional SCRAPI by propagating additional sources of uncertainty through the estimation process:

- **Guidance Efficiency (GE) uncertainty** via Bayesian hierarchical modeling  
- **Genetic Stock Identification (GSI) uncertainty** via posterior draws  
- Full integration with bootstrap-based escapement estimation  

This approach aligns juvenile escapement estimation with the **EASE framework** used for adult returns, enabling consistent uncertainty propagation across life stages.

---

## Key Features

- Estimate smolt escapement using `SCRAPI2`
- Fit Bayesian models for GE as a function of spill conditions
- Generate posterior draws of GE for uncertainty propagation
- Incorporate GSI uncertainty into stock-specific estimates
- Maintain compatibility with legacy `SCOBI::SCRAPI` workflows

---

## Installation

Install directly from GitHub:

```r
remotes::install_github("dsmi313/smoltEASE")
