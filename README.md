# The Tail You Cannot Fit

> Extreme Value Theory and Operational Loss Capital Estimation

Author: Patrick Lefler

Published: 2026-06-17

Rendered:

## Project Introduction

> Quantifies how much a log-normal capital model understates operational loss tail risk versus a Generalized Pareto fit, using NexaCore's loss register.

## Overview

This project applies the Peaks-Over-Threshold method from Extreme Value Theory to a synthetic operational loss register for NexaCore Financial Technologies, fitting a Generalized Pareto Distribution to losses exceeding a threshold selected via mean excess and shape-parameter stability diagnostics. The resulting Value-at-Risk and Expected Shortfall estimates, at the 99th, 99.9th, and 99.95th percentiles, are compared against a log-normal benchmark fit to sub-threshold losses and against historical simulation. The intended outcome for decision-makers is a concrete, euro-denominated answer to a question most capital models never test directly: how much does the choice of loss distribution change the capital figure a CRO or capital adequacy committee relies on.

## Tech Stack

* **Language:** Polyglot (R, Python)
* **Framework:** [Quarto](https://quarto.org/)
* **Primary Libraries:** R — `evd`, `POT`, `fitdistrplus`, `tidyverse`, `ggplot2`, `patchwork`, `kableExtra`. Python — `scipy`, `numpy`, `pandas`
* **Deployment/Output:** Self-contained HTML report

## Repository Structure

```
evt-operational-loss/
├── data/                       # Synthetic loss register and all derived EVT outputs
├── scripts/
│   ├── generate_losses.py      # Synthetic loss register generation
│   └── compute_evt.R           # GPD/GEV fitting, capital comparison, diagnostic figures
├── output/                     # Diagnostic PNGs and rendered HTML
├── _brand.yml                  # Brand configuration
├── _quarto.yml                 # Project configuration
├── INSTRUCTIONS.md             # Project specification and decision log
└── index.qmd                   # Main Quarto entry point
```

## Key Findings

> GPD-derived Value-at-Risk exceeds a conventional log-normal benchmark at every confidence level tested, by 110.6% at the 99th percentile, 63.0% at the 99.9th, and 52.0% at the 99.95th. The percentage gap narrows as confidence rises, but the euro gap holds roughly steady between €117,700 and €153,500 across all three levels, meaning the undercapitalization risk does not diminish where it matters most for regulatory capital adequacy.

> The fitted shape parameter (ξ = 0.1195) places NexaCore's loss-generating process in the heavy-tailed Fréchet domain, though the standard error is wide enough that the confidence interval does not exclude zero — a direct, honest consequence of estimating tail behavior from 83 exceedances rather than thousands, and itself a finding rather than a limitation to hide.

> Historical simulation's Expected Shortfall estimate flattens at the single largest loss in the register once required confidence exceeds what the sample can support, illustrating why an empirical method without a theoretical basis for extrapolation cannot function as a genuine tail risk model.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

Patrick Lefler [https://www.linkedin.com/in/patricklefler/](https://www.linkedin.com/in/patricklefler/) | [patrick-lefler.github.io](https://patrick-lefler.github.io) | [https://substack.com/@pflefler](https://substack.com/@pflefler)
