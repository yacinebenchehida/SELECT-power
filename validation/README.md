# Validation of rescaling

SPOT simulates long periods (more than 2,000 generations) with rescaling: each simulated generation stands for λ real ones, with selection (1+s)^λ − 1, breeding numbers N/λ, migration 1 − (1−m)^λ and sampling times t/λ (see the README, "Very long series: rescaling"). These scripts compare rescaled results with full simulations of every generation.

| File | Content |
|---|---|
| `rescale_validation.js` | Runs the same designs with full and rescaled simulations and writes `rescale_results.tsv`. |
| `rescale_limits.js` | Forces increasing scaling factors to find where rescaling breaks down; used to set the limits. |
| `make_figure.py` | Draws `rescale_validation.png` (full against rescaled power, and neutral spread of allele frequency) from the results table; the figure itself is not stored in the repository. |
| `rescale_results.tsv` | Results of the last validation run. |

Run from the repository root:

```
node validation/rescale_validation.js > validation/rescale_results.tsv   # about 15 minutes
node validation/rescale_limits.js                                         # about 2 minutes
python3 validation/make_figure.py
```

## Design of the comparison

- Four sampling points equally spaced over 5,000, 20,000 or 60,000 generations; 20, 20, 20 and 40 individuals; default sequencing settings.
- Haploid and diploid; breeding numbers 1,000, 10,000 and 100,000; s = 0.0001, 0.0005 and 0.002; starting frequency 0.01 and 0.1. Only a subset was run at 60,000 generations, to keep the full simulations affordable.
- 1,000 replicates under selection and 2,000 neutral for each estimate, with different seeds for the full and rescaled runs.
- Scaling factors chosen automatically by SPOT: λ = 3 to 30.

## Results

- **Power** (80 designs): the mean difference between rescaled and full simulations was −0.3 percentage points, and the mean absolute difference 1.2 points. 5 of the 80 differences exceeded 2 Monte Carlo standard errors (about 4 expected by chance alone) and 1 exceeded 3; the root-mean-square standardised difference was 1.13 (1 expected if the only difference is Monte Carlo noise). The largest difference (8 points) was for a haploid population of 1,000 at the edge of the limits below.
- **Neutral drift**: the spread of the neutral allele frequency at the last point agreed between full and rescaled runs, and with the theoretical expectation.
- **Speed**: rescaled runs were a median of 4 times faster in this grid, and up to 29 times faster; the gain grows with the number of generations.

## Limits

`rescale_limits.js` forced λ from 2 to 200. Results stayed within Monte Carlo error while the rescaled selection coefficient was at most about 0.025 and the rescaled breeding numbers at least about 200; beyond that (for example scaled s ≈ 0.05, or scaled N ≈ 50), power differed by 5–7 points. SPOT therefore caps λ so that scaled s ≤ 0.025 and scaled N ≥ 200. If these caps would leave more than 10,000 generations to simulate, λ is raised anyway to keep the page fast, and the page marks the result as approximate.
