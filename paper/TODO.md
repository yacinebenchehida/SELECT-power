# SPOT paper: to-do list

Working checklist for the SPOT paper. Related literature is listed in [`RELATED_WORK.md`](RELATED_WORK.md).

## 1. Software before submission

- [x] **Rescaling for very many generations** (historical DNA of bacteria and viruses): done. λ chosen automatically (about 2,000 simulated generations), capped at scaled s ≤ 0.025 and scaled N ≥ 200; validated against full simulations; summary and figure in [`validation/rescaling.md`](../validation/rescaling.md).
- [x] **Use all processor cores (Web Workers)**: done for the power curve, planner and sampling comparison; results identical to serial runs.
- [x] **Save, load and share settings; export results as CSV**: done (shareable link, settings file, CSV downloads). Useful for the paper: each figure can come with a link that reproduces it.
- [ ] Final pass on the page, README and guide after these changes; check every number quoted in the guide.
- [ ] Organism-specific advice in the README ("Choosing settings for your organism": insects, long-lived vertebrates, plants, pathogens).
- [ ] Freeze a version for the paper (currently 0.9.1; release 1.0.0 on GitHub and archive it on Zenodo to get a DOI, which the journal requires).
- [x] **Licence**: MIT (`LICENSE`), with `CITATION.cff` for GitHub's "Cite this repository" button.
- [x] **Example scenarios** in the page (bumblebees, long-lived mammal, human ancient DNA, bacterial pathogen, large diploid population, analytic baseline), version number and "How to cite", and PNG/SVG download of every chart.

## 2. Validation (main results section)

- [x] **False-positive rate**: done, see [`validation/false_positive_rate.md`](../validation/false_positive_rate.md). The drift-aware test holds 5% (trend test 4.5–5.4%) across Ne, inheritance systems and design factors; the naive test reaches up to 23%.
- [x] **Gaussian fallback for many loci fixed** (version 0.9.1): SPOT now simulates 10/α′ neutral replicates (up to 200,000), which keeps the false-positive rate at 1.0–1.2 times nominal; the earlier Gaussian fit gave up to 8 times. See the same file.
- [x] **Agreement with theory**: done, see [`validation/theory.md`](../validation/theory.md) (exact power without drift, heterozygosity loss under drift, selection trajectories).
- [x] **Monte Carlo precision**: done, see [`validation/monte_carlo_precision.md`](../validation/monte_carlo_precision.md).
- [x] **Rescaling against full simulations** across s, Ne, ploidy and number of generations: done, see [`validation/rescaling.md`](../validation/rescaling.md) (methods, results, figure and caveats, ready to adapt as a results subsection or supplement).
- [x] **Comparison with SLiM**: done, see [`validation/slim.md`](../validation/slim.md) (6 scenarios, diploid, haploid and haplodiploid, full simulation; rescaling not compared with SLiM). Small excess of false positives (5.6–5.9%) in two scenarios, cause not identified.

## 3. Worked examples (application section)

- [x] **Main example: museum bumblebee design**: done, see [`worked_examples/README.md`](../worked_examples/README.md). Detectable s by starting frequency; Ne anchored on heterozygosity loss and power against Ne; naive test false positives; where to add sampling effort, with a contrast between strong drift (Ne = 1,500) and weak drift (Ne = 100,000).
- [x] **Short contrasting example**: done (long-lived mammal, 7 generations, with and without a bottleneck, using Ne per interval).
- [x] Optional pathogen example (haploid, 17,000 generations, rescaled): done.

## 4. Introduction and statement of the gap

- [ ] Time series from museum collections, ancient DNA and monitoring are increasingly used to detect selection, and many inference methods exist (WFABC, ApproxWF, diplo-locus, He et al., and others).
- [ ] No easy way to check beforehand whether a sampling design can detect selection. Existing design work is either lab evolve-and-resequence (command-line simulators) or bespoke, scenario-specific simulations (Anchieri et al. 2026 for human aDNA; SimOutbreakSelection for epidemics). Analytical theory (Lynch & Ho 2020) has no software and ignores genotype errors.
- [ ] Drift and sampling can masquerade as selection: naive tests are anti-conservative under drift (Spitzer et al. 2020; Wiberg et al. 2017); the Black Death debate (Klunk et al. 2022; Barton et al. 2025; Vilgalys et al. 2025) is a high-profile example.
- [ ] SPOT fills the gap: instant, installation-free, drift-aware, and models the complications of historical samples (few individuals, low coverage, haploid/haplodiploid systems, changing Ne).
- [ ] Comparison table: SPOT against MimicrEE2, SimOutbreakSelection, Anchieri et al. 2026 and diplo-locus, with columns: web-based, power output, drift-aware null, genotype calling, haplodiploidy.
- [ ] Verify every reference in `RELATED_WORK.md` before citing (a few are marked "not verified").

## 5. Methods section

- [ ] Model: Wright–Fisher with selection (haploid, diploid, haplodiploid), sex-specific breeding numbers, change in population size, Ne per interval, migration.
- [ ] Sampling and sequencing: individuals, siblings, Poisson depth, genotype calling by maximum likelihood, minimum depth.
- [ ] Test and power: trend (Cochran–Armitage) statistic; threshold from neutral replicates (empirical quantile, Gaussian fit when α/m is too small); Bonferroni across candidate loci; naive test for comparison.
- [ ] Outputs: power tiles, power curve, sample-size planner.
- [ ] Assumptions and limitations (single locus, no linkage or structure, discrete generations).
- [ ] Code and data availability: GitHub URL, Zenodo DOI, licence.

## 6. Journal and submission

- [ ] Target order:
  1. *Molecular Ecology Resources* (best fit for museum, aDNA and temporal monitoring).
  2. *Methods in Ecology and Evolution* (Application article).
  3. *Ecology and Evolution* (open-access fallback).
- [ ] Check each journal's current author guidelines (article type, word limit, fees) before writing.
- [ ] Keep it short: software description, validation, worked example.
- [ ] Show the paper, README and app to colleagues and ask them to review before submission.
