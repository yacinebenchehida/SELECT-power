# SPOT paper: to-do list

Working checklist for the SPOT paper. Related literature is listed in [`RELATED_WORK.md`](RELATED_WORK.md).

## 1. Software before submission

- [ ] **Rescaling for very many generations** (historical DNA of bacteria and viruses). Simulate at most about 2,000 generations, with selection made proportionally stronger and the population proportionally smaller (Ne·s and total drift unchanged). Choose the scaling factor automatically, keep scaled s below about 0.1, say on the page when rescaling is used, and allow generation times down to about one hour (0.0001 years). Validate against full simulations (for example 20,000 generations) before release.
- [ ] **Use all processor cores (Web Workers)**. Number of threads chosen automatically; a fixed random seed per block of replicates so results do not depend on the number of threads.
- [ ] Optional: **replicate populations (1 to 100)**, evidence combined across replicates.
- [ ] Final pass on the page, README and guide after these changes; check every number quoted in the guide.
- [ ] Organism-specific advice in the README ("Choosing settings for your organism": insects, long-lived vertebrates, plants, pathogens).
- [ ] Freeze a version for the paper (GitHub release with a version number; archive it on Zenodo to get a DOI).
- [ ] Choose and add a licence before publication (for example MIT).

## 2. Validation (main results section)

- [ ] **False-positive rate is correct.** With no selection, across small and large Ne, 2 to 10 time points, haploid, diploid and haplodiploid systems, low and high depth: the drift-aware test should reject about 5% of neutral replicates, while the naive test rejects far more when drift is strong.
  - Figure: false-positive rate against Ne, one line for each test.
- [ ] **Agreement with known theory.**
  - With no drift and perfect genotypes, power should match the exact binomial (textbook) calculation.
  - The neutral variance of allele-frequency change should match p(1−p)·t/(2Ne) (plus sampling variance).
  - Haploid deterministic trajectory against p(1+s)^t / (1 − p + p(1+s)^t).
  - Table or scatter plot of SPOT against theory (points on the diagonal).
- [ ] **Monte Carlo precision.** Repeat runs with different seeds: power varies by about ±1–2 points with the default number of replicates. This justifies the defaults.
- [ ] If rescaling is added: rescaled against full simulations across s, Ne and number of generations.

## 3. Worked examples (application section)

- [ ] **Main example: museum bumblebee design.** Four decades (1930s, 1950s, 1970s, 2020s, entered as one representative year each), females only, haplodiploid.
  - Ne anchored on the heterozygosity loss observed between the first and last samples.
  - Plausible ranges of s and starting frequency.
  - Planner and power curves: for example, whether adding historical specimens helps more than enlarging the modern sample; below which Ne no sample size reaches 80% power.
  - One clear take-home message (for example: "s ≥ 0.05 on a common allele is detectable; rare alleles are not, whatever the sample size").
- [ ] **Short contrasting example**: a diploid long-lived species (long generation time) or a bottleneck between two sampling dates using Ne per interval.
- [ ] Optional: a pathogen example (haploid, many generations) if rescaling is ready.

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
