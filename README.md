# OHARP_econ
# OHARP Economic Evaluation

A One Health modelling analysis of the economic burden of ESBL-producing Enterobacterales in Singapore over **2026–2030**.

The analysis estimates infection-related healthcare costs and productivity losses, compares hospital and community hand-hygiene scenarios, and examines uncertainty through probabilistic sensitivity analysis (PSA) and value-of-information analysis.

## Repository Contents

| File | Purpose |
|---|---|
| `total_cost.R` | Estimates baseline burden and compares hand-hygiene scenarios |
| `attributable_burden.R` | Estimates costs and deaths attributable to antimicrobial resistance |
| `voi.R` | Calculates estimation-based EVPI and EVPPI |
| `data/` | Contains model parameters, calibrated posterior samples and saved simulation draws |

Each script runs independently. No patient-level records are included.

## Getting Started

Install the required R packages:

```r
install.packages(c(
  "deSolve", "dplyr", "tidyr", "ggplot2", "scales",
  "foreach", "doSNOW", "doParallel", "purrr", "voi"
))
```

Set your working directory to the repository folder and run the analysis you need:

```r
source("total_cost.R")
source("attributable_burden.R")
source("voi.R")
```

Results are saved in the `outputs/` folder.

## Reproducing Results

By default, the total-cost and attributable-burden scripts reproduce the reported tables from the included **1,000 PSA draws**. The VOI script refits the uncertainty analysis using the saved baseline draws.

To run a new PSA, set `RUN_MODE <- "fresh"` at the top of `total_cost.R` or `attributable_burden.R`. The number of iterations and random seed can also be adjusted there.

## Interpretation

- Cost summaries are reported in **million Singapore dollars (SGD)**, with **95% uncertainty intervals**.
- Five-year totals are calculated within each PSA iteration before summarisation.
- Intervention comparisons estimate infection-related cost changes; programme implementation costs are not included.
- Resistance-attributable burden is calculated separately from intervention effects.
- VOI measures uncertainty in projected costs, not the monetary value of a policy decision.
