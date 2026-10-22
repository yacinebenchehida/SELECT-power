# Validation

Five analyses checking that SPOT gives correct results. Each has its own page with methods, results, a figure and the commands to reproduce it; they are written so they can be adapted for the paper.

| Analysis | Question | Main result |
|---|---|---|
| [False-positive rate](false_positive_rate.md) | Without selection, how often does each test declare selection? | SPOT's drift-aware test rejects 4.5–5.4% of neutral loci at α = 0.05 (trend test; 3.7–5.2% for the first-vs-last test) across effective sizes, inheritance systems and designs; the naive test rejects up to 23%. With many candidate loci, SPOT simulates more neutral replicates, which keeps the rate at 1.0–1.2 times nominal (an earlier Gaussian fit gave up to 8 times). |
| [Agreement with theory](theory.md) | Does SPOT reproduce results known without simulation? | Power without drift matches an exact calculation (240 comparisons, within Monte Carlo error); loss of heterozygosity under drift matches theory within 0.75%; selection trajectories match the closed form (haploid) or an independent recursion (diploid, haplodiploid) to within 2 × 10⁻¹⁵. |
| [Comparison with SLiM](slim.md) | Does an independent, individual-based simulator give the same results? | Same scenarios in SLiM 5.2 (diploid, haploid, haplodiploid; 2,000 replicates each): allele-frequency distributions agree (all KS p > 0.15), power differs by at most 4 points, and SPOT's threshold gives 4.7–5.9% false positives on SLiM's neutral populations (slightly above 5% in two scenarios). |
| [Monte Carlo precision](monte_carlo_precision.md) | How much does power vary between runs? | SD of 1.4–2.3 percentage points at the default 1,000 replicates (trend test; up to 3.1 for the first-vs-last test), falling roughly as 1/√R. |
| [Rescaling](rescaling.md) | Are very long series (more than 2,000 generations) simulated accurately? | Rescaled and full simulations agree within Monte Carlo error (80 designs). |

All SPOT simulations call the same code as the page (`sim.js`). The theoretical values are computed independently in `theory_check.py`, and the independent simulations are run with SLiM (`slim_check.slim`). Figures are drawn with `svgplot.py` (Python standard library only) so they display directly on GitHub; `theory_check.py` also needs NumPy.
