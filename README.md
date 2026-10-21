# SPOT-simulator: Selection Power Over Time

**▶ Run it in the browser: <https://yacinebenchehida.github.io/SPOT-simulator/>**

A browser-based tool for estimating the **statistical power to detect selection at a single candidate locus** from allele frequencies sampled at **two or more points in time** (for example museum specimens plus a modern sample, ancient DNA series, or any time series of population samples). Generation time is set in years, so it covers both short-lived species (several generations per year) and long-lived ones (decades per generation).

It simulates populations forward in time with **genetic drift**, **haploid, diploid or haplodiploid inheritance**, **sex-specific breeding numbers and selection**, then samples individuals and passes them through a **low-coverage sequencing and genotype-calling** step. Power is computed against a **drift-aware null**: the test threshold comes from simulated neutral replicates, not from the textbook N(0, 1) distribution.

Everything runs locally in the browser. A typical run (1,000 replicates under selection plus 2,000 neutral replicates) takes about 0.1–0.4 s.

**New to the method?** [SIMULATION_GUIDE.md](SIMULATION_GUIDE.md) explains step by step, with worked numerical examples, how drift, selection and sampling are simulated and how a power value is obtained.

## What power means, in plain words

**Power is the probability that a sampling design detects selection when selection is actually acting.** A power of 80% means that selection of the chosen strength is detected in 8 cases out of 10 and missed in 2.

How the simulator estimates it:

1. **Replay history many times.** The population's past is replayed on the computer: by default 1,000 times with the allele favoured by selection, and 2,000 times with no selection at all.
2. **Sample and sequence.** In every replay, specimens are drawn at each sampling date and "sequenced" with the same limits as real data: a limited number of specimens per date, and some genotypes read wrongly or not at all.
3. **Score each replay.** One number per replay says how clearly the allele became more common over time. A high score means a clear rise.
4. **Set the threshold with the no-selection replays.** Even without selection, allele frequencies wander up and down by chance (drift), a lot in small populations. The threshold is set so that only 5% of these chance-only replays clear it (5% is the significance level, adjustable in the settings).
5. **Count.** Power is the share of the selection replays that clear the threshold.

Why the threshold matters: a standard test that ignores drift sets the threshold too low, so it often "finds" selection where there is none. If power stays low whatever the sample size, chance drift is the limit, not the number of specimens; the sample-size planner shows this directly.

---

## Contents

- [What power means, in plain words](#what-power-means-in-plain-words)
- [Quick start](#quick-start)
- [What the simulator does, step by step](#what-the-simulator-does-step-by-step)
- [Parameters](#parameters)
  - [Inheritance and population](#1-inheritance-and-population)
  - [Selection](#2-selection)
  - [Sampling design](#3-sampling-design)
  - [Sequencing and genotype calls](#4-sequencing-and-genotype-calls)
  - [Test and simulation](#5-test-and-simulation)
  - [Presets](#presets)
- [Choosing values that are unknown](#choosing-values-that-are-unknown)
- [Outputs](#outputs)
- [How power is computed](#how-power-is-computed)
- [Model equations](#model-equations)
- [Using the simulation core from JavaScript / Node](#using-the-simulation-core-from-javascript--node)
- [Assumptions and limitations](#assumptions-and-limitations)
- [Files](#files)

---

## Quick start

1. Open the [online version](https://yacinebenchehida.github.io/SPOT-simulator/), or download `index.html` and open it in any modern browser. No installation or server is needed.
2. Set the parameters that describe **the study design**: inheritance system, years per generation, sampling years and individuals per time point (add or remove time points as needed), sequencing depth and number of candidate loci.
3. Use the **parameter power curve** to explore the parameters that are **unknown**, such as effective population size, starting frequency, selection strength, dominance and onset of selection. Viewing power over a range is more informative than guessing one value.
4. Read the first tile, **"Power, trend test across all time points, with a drift-aware threshold"**. This is the headline result.

Parameters with no reason to change can stay at their defaults. The defaults are neutral: constant population size, no migration, the same selection in both sexes, unrelated individuals, and female-only sampling.

---

## What the simulator does, step by step

For each replicate:

1. **Start** at the first sampling point with the selected allele at frequency **p₀** in both sexes.
2. **Each generation**:
   - form female zygotes from the breeding females and males (in the haploid system, each individual simply carries one allele);
   - apply **selection** in females (diploid) and in males (haploid or diploid);
   - apply **drift** by binomially sampling the next generation's breeding females and males;
   - optionally apply **population change** and **migration**.
3. **At each sampling year**, draw the chosen number of individuals:
   - diploid individuals carry 2 alleles;
   - haploid males carry 1 allele;
   - individuals can be unrelated or grouped into sibling families.
4. **Sequence** each individual: read depth ~ Poisson(mean depth), alternative-allele reads ~ Binomial(depth, expected fraction). The genotype is called by maximum likelihood, and individuals below the minimum depth are dropped.
5. **Compute** a one-sided Cochran–Armitage trend statistic **Z** from the called allele counts. There are two versions: all sampling points, and the first vs the last point only.

This is done **R times with selection** and **2R times without selection** (s = 0, everything else identical). The threshold is the (1 − α/m) quantile of the neutral Z values. **Power** is the fraction of selection replicates above that threshold.

---

## Parameters

Each entry gives the label used in the interface, the key in the JavaScript parameter object, the default value, and practical guidance.

### 1. Inheritance and population

| Interface label | Key | Default | Meaning |
|---|---|---|---|
| Haplodiploid / Diploid / Haploid | `system` | `'haplodiploid'` | Inheritance system. **Haplodiploid**: females are diploid and males are haploid, developing from unfertilised eggs, so a son inherits only his mother's allele. **Diploid**: both sexes are diploid. **Haploid**: every individual carries one copy of the locus (for example haploid organisms or a haploid life stage); there are no sexes, no dominance and no families, so the settings for males, dominance and siblings are hidden. |
| Years per generation | `genTime` | `1` | Generation time in years (between 0.05 and 100). The number of generations simulated is (last year − first year) / this value, rounded. Use 0.5 for two generations a year, 1 for one, and the average parental age for long-lived species (around 25–30 for humans). Older parameter sets that use `volt` (generations per year) are still accepted. |
| Breeding females per generation (haploid: breeding individuals) | `Nf` | `1000` | Number of females that actually reproduce each generation. This is **not** a census count of all individuals, and it is not Ne itself (see below). In social insects this is roughly the number of successful reproductive females, not the number of workers. |
| Breeding males per generation | `Nm` | `1000` | Number of males that actually reproduce each generation. Raise it relative to `Nf` to approximate multiple mating by females. Not used in the haploid system. |
| No drift (infinite population) | `noDrift` | `false` | Switches drift off. Only useful to reproduce a deterministic calculation. Leave it off for realistic power. |
| Set Ne separately for each interval | `neByBlock`, `neBlocks` | off | When ticked, an extra column "Ne to next" appears in the sampling design. Each value is the effective population size from that sampling point to the next (`neBlocks` has one value per interval, so one fewer than the number of points). Ne changes abruptly at the sampling dates, which gives simple histories such as a bottleneck between two dates followed by recovery. Breeding females and males are scaled to give exactly that Ne, keeping their ratio. This replaces the breeding numbers and the exponential change. |
| Breeding numbers at the last point relative to the first | `decline` | `1` (100%) | Exponential change in `Nf` and `Nm` between the first and last sampling year. 1 means constant size, 0.1 a ten-fold decline, 2 a doubling. Use it to ask "what if the population declined?". Otherwise keep it at 1. |
| Migrants per generation (fraction) | `mig` | `0` | Each generation, this fraction of breeding adults is replaced by migrants. **Keep it at 0** unless immigration from another population is to be modelled. |
| Allele frequency in the migrant source | `pMig` | `-1` (= p₀) | Allele frequency among migrants. `-1` means "same as p₀", a source where the allele is not under selection. Only used when `mig > 0`. |

**How breeding numbers relate to Ne.** The simulator converts `Nf` and `Nm` into an effective population size, shown under the sliders and in the Ne tile:

- haplodiploid: **Ne = 9·Nf·Nm / (4·Nm + 2·Nf)**. With equal numbers N of each sex, Ne = 1.5 N.
- diploid: **Ne = 4·Nf·Nm / (Nf + Nm)**. With equal numbers N of each sex, Ne = 2 N.

For most purposes **Ne is the only thing that matters**. If the sex-specific numbers are unknown, keep `Nf = Nm` and treat the pair as a single dial for Ne.

### 2. Selection

| Interface label | Key | Default | Meaning |
|---|---|---|---|
| Allele frequency at the first point (p₀) | `p0` | `0.01` | Frequency of the favoured allele at the first sampling year. Low values (≤ 0.01) represent a rare variant, and values around 0.1 represent common standing variation. |
| Selection coefficient s | `s` | `0.05` | Relative fitness advantage per generation. In diploids the genotype fitnesses are 1+s (AA), 1+hs (Aa) and 1 (aa). Values of 0.01–0.05 are common for natural selection. Strong selection, such as insecticide resistance, can be an order of magnitude higher. |
| Dominance in diploid individuals | `h` | `0.5` | Effect of the allele in heterozygotes: 1 = dominant, 0.5 = additive, 0 = recessive. It applies to all diploid individuals. Haploid males have no heterozygotes. If unknown, keep 0.5 and check 0 and 1 with the parameter power curve. |
| Selection in males (× s) | `sMale` | `1` | Selection strength in males as a multiple of `s`. In haplodiploids, haploid males with the allele have fitness 1 + s·sMale. **Keep it at 1** when the trait is expressed in both sexes. Set 0 if the trait is not expressed in males, or values above 1 if the effect is stronger in males. |
| Selection starts in year | `onset` | year 1 | Calendar year when selection begins. Before it, the allele only drifts. Selection acts from the generation after this year. |

### 3. Sampling design

| Interface label | Key | Default | Meaning |
|---|---|---|---|
| Time points (Year) | `years` | `[1935, 1955, 1975, 2025]` | Calendar years of the sampling points: 2 to 20 points, added or removed with **+ Add time point** and **×**. Years must increase strictly and can be negative (before year 0). The first point defines the start of the simulation, and the whole period can span at most 5,000 generations. Samples collected over a period (for example a decade) can be pooled and given one representative year, such as the median. |
| Time points (Individuals) | `n` | `[20, 20, 20, 40]` | Number of **individuals** sampled at each point (up to 1,000). The first and last points need at least 1. Intermediate points can be set to 0 to drop them; with only the first and last points sampled, the trend test equals the two-point test. If two sampling years fall in the same generation, both samples come from the same population. |
| Males in the sample | `maleFrac` | `0` | Probability that each sampled individual is male. In a haplodiploid, males are haploid and contribute **1 allele**, against 2 for a female. In a diploid, males contribute 2 alleles and differ from females only if `sMale ≠ 1`. 0 means females only. |
| Sampled individuals per family | `kSib` | `1` | 1 means every individual comes from a different family. k > 1 samples groups of k siblings sharing a mother and a father. Siblings carry correlated alleles: in haplodiploids, sisters share their father's allele exactly. This reduces the effective sample size. Raise it only if related individuals are expected in the sample. |

### 4. Sequencing and genotype calls

| Interface label | Key | Default | Meaning |
|---|---|---|---|
| Perfect genotypes | `perfect` | `false` | Ignores sequencing: true genotypes are observed directly. |
| Mean depth, earlier points | `depthEarly` | `10` | Mean read depth at the locus for individuals sampled at every point except the last. Depth per individual ~ Poisson(mean). |
| Mean depth, last point | `depthLast` | `15` | Mean read depth for individuals at the last point. Separate values let earlier samples, such as older or degraded material, have lower coverage. |
| Per-read error rate | `err` | `0.005` | Probability that a read shows the wrong allele (sequencing error plus residual damage). |
| Minimum depth to call a genotype | `minDP` | `3` | Individuals with fewer reads at the locus are treated as missing. |

Diploid genotypes are called by maximum likelihood among the three genotypes with a flat prior. At low depth, heterozygotes are often called homozygous (allelic dropout). Haploid males are called by majority of reads, with ties broken at random.

### 5. Test and simulation

| Interface label | Key | Default | Meaning |
|---|---|---|---|
| Significance level α | `alpha` | `0.05` | Per-study significance level. |
| Candidate loci tested (Bonferroni m) | `m` | `1` | Number of candidate loci to be tested. Each is judged at α/m. |
| Replicates under selection | `reps` | `1000` | Number of replicates with selection. Twice as many neutral replicates are run. More replicates give a more precise estimate, not a different answer. With 1,000 replicates, a power of 0.80 has a standard error of about ±0.013. |
| Random seed | `seed` | `1` | Makes a run exactly reproducible. |

If α/m is so small that fewer than 10 neutral replicates are expected beyond the threshold, the threshold is taken from a **Gaussian fit** to the neutral Z values instead of the empirical quantile. The status line says when this happens, and hollow markers show it in the parameter power curve.

### Presets

- **Reset to defaults** restores the neutral defaults listed above.
- **Analytic baseline** is diploid, with no drift, perfect genotypes and a dominant allele. Only binomial sampling noise remains, which reproduces a classic deterministic power calculation. It is useful as a reference point.

---

## Choosing values that are unknown

**Set these to match the study design:** inheritance system, years per generation, sampling years, n per point, sex of sampled individuals, sequencing depth, and number of candidate loci.

**Explore these with the parameter power curve:** Ne (via breeding numbers, or per interval for a known demographic history), p₀, s, h and onset year. They are rarely known, and power depends strongly on them.

**Leave these at their defaults unless there is a specific reason:** no drift (off), population change (1), Ne per interval (off), migration (0), selection in males (1), individuals per family (1), error rate, minimum depth, α, replicates and seed.

### Anchoring Ne with heterozygosity data

If the change in genome-wide heterozygosity between the first and last sampling points can be measured, it gives a plausible Ne. Neutral heterozygosity decays as

```
H_t / H_0 = (1 − 1/(2·Ne))^G          G = number of generations
Ne = 1 / ( 2 · (1 − (H_t/H_0)^(1/G)) )
```

Heterozygosity lost over **100 generations** for a range of Ne:

| Ne | 50 | 100 | 250 | 500 | 1,000 | 2,500 | 5,000 | 10,000 |
|---|---|---|---|---|---|---|---|---|
| Heterozygosity lost | 63% | 39% | 18% | 9.5% | 4.9% | 2.0% | 1.0% | 0.5% |

The tile "Neutral heterozygosity lost from the first to the last point" shows this value for the current settings, so the breeding numbers can be tuned until it matches observed data.

---

## Outputs

### Tiles

| Tile | Meaning |
|---|---|
| **Power, trend test across all time points, drift-aware threshold** | Headline result: the fraction of selection replicates whose trend Z (all points) exceeds the (1 − α/m) quantile of neutral replicates. |
| Power, two-point test, drift-aware threshold | The same, using only the first and last points. |
| False-positive rate of a naive test | The fraction of **neutral** replicates rejected by the textbook threshold z₁₋α/m, which ignores drift. Shown in red when it clearly exceeds α/m. |
| Power that naive test would report | Power using the textbook threshold. It is misleading when the false-positive rate above is inflated. |
| Expected frequency at the last point without drift | Deterministic allele frequency at the last sampling point. |
| Neutral heterozygosity lost | 1 − Π(1 − 1/(2·Ne_g)) over the generations from the first to the last point (1 − Π(1 − 1/N_g) in the haploid system). |
| Effective population size Ne | Computed from the breeding numbers (start → end if the population changes), or the lowest and highest values when Ne is set per interval. |
| Genotypes called correctly | Accuracy of heterozygote calls (or haploid calls), and the share of individuals below the minimum depth. |

### Charts

- **Allele trajectories:** 40 selection replicates, 40 neutral replicates, the deterministic expectation, and the sampling points.
- **How power is read off:** histograms of Z under neutrality and under selection, with the drift-aware threshold (solid line) and the naive threshold (dashed line). Power is the shaded area of the selection histogram beyond the solid line.
- **Parameter power curve:** shows how the overall power changes when one parameter changes. All values set by the user stay fixed except the chosen parameter, which is stepped through a range of values (for example 10, 20, 30 … 120 individuals at the last point). At each value the whole power calculation is rerun from scratch: new replays with and without selection, a new threshold, and the power of the complete design. Changing the sample size at the last point therefore gives the power of the design 20 / 20 / 20 / 10, then 20 / 20 / 20 / 20, and so on, with the earlier points kept at their current values. Blue: trend test across all points. Green: two-point test. Red: false-positive rate of the naive test. The vertical line marks the current value and the dashed line marks 80% power.
- **Sample-size planner:** a grid of power (trend test across all points, drift-aware) with rows that scale the whole sampling design (×0.5, ×1, ×2, ×4) and columns that vary one parameter (s, h, Ne or p₀). Cells reaching 80% are outlined. Each cell uses up to 500 replicates under selection and 1,000 neutral.
- **Sampling haploid males vs diploid females:** power for four sampling designs under haplodiploidy, with genotype-call accuracy shown on hover. Shown only when the haplodiploid system is selected (in a diploid, both sexes carry two alleles, so the sex of sampled individuals does not matter).

---

## How power is computed

**Test statistic.** At sampling point i (time tᵢ in generations since the first point), xᵢ copies of the selected allele are called among nᵢ alleles:

```
Z = Σᵢ xᵢ (tᵢ − t̄) / sqrt( p̄ (1 − p̄) Σᵢ nᵢ (tᵢ − t̄)² )
t̄ = Σ nᵢ tᵢ / Σ nᵢ      p̄ = Σ xᵢ / Σ nᵢ
```

This is the one-sided Cochran–Armitage trend test on allele counts. The two-point version uses the first and last points only and reduces to a one-sided two-proportion z-test.

**Why a drift-aware threshold?** With binomial sampling noise only, Z ~ N(0, 1) under the null. Drift adds variance between generations, so the real neutral distribution of Z is wider. A test using z₁₋α would then flag too many neutral loci as selected. The simulator therefore takes the threshold from neutral replicates, which experience the same drift, sampling and genotyping error:

```
critical value  c = quantile_{1 − α/m}( Z_neutral )
power           = P( Z_selection > c )
naive FPR       = P( Z_neutral > z_{1 − α/m} )
```

With real data, the equivalent is an **empirical null** from genome-wide SNPs matched on starting frequency and measured in the same samples.

**Analytic baseline, for comparison.** A deterministic calculation assumes the frequency follows a fixed path from p₀ to pₜ, with binomial sampling only:

```
Power = Σₓ Σᵧ Bin(x; n₁, p₀) · Bin(y; n₂, pₜ) · 1[ test(x, y) rejects at α ]
```

The "Analytic baseline" preset reproduces this setting within Monte Carlo error.

---

## Model equations

Per generation g, with p_f and p_m the allele frequencies in breeding females and males:

```
Female zygotes (diploid, haplodiploid):  AA = p_f·p_m   Aa = p_f(1−p_m) + p_m(1−p_f)   aa = (1−p_f)(1−p_m)
Selection in females:           weights 1+s, 1+hs, 1      →  p_f* = (AA(1+s) + ½Aa(1+hs)) / w̄
Males, haplodiploid:            allele from mother only; weights 1+s♂, 1
                                →  p_m* = p_f(1+s♂) / (p_f(1+s♂) + 1 − p_f)
Males, diploid:                 same zygote frequencies as females; weights 1+s♂, 1+hs♂, 1
Male selection:                 s♂ = s × sMale
Haploid system:                 weights 1+s, 1  →  p* = p(1+s) / (1 + p·s);   drift p' = Bin(N(g), p*) / N(g);   Ne = N
Drift:                          p_f' = Bin(2·Nf(g), p_f*) / (2·Nf(g))
                                p_m' = Bin(Nm(g), p_m*) / Nm(g)            (haplodiploid)
                                p_m' = Bin(2·Nm(g), p_m*) / (2·Nm(g))      (diploid)
Population change:              N(g) = N × decline^(g/G),  or, with per-interval Ne, N(g) scaled so that Ne(g) = Ne of the interval containing g
Migration (after drift):        p' ← (1 − mig)·p' + mig·p_source           (both sexes)
Sampling, kSib = 1:             each individual is male with probability maleFrac;
                                diploid genotypes drawn from that sex's post-selection frequencies;
                                haploid males carry allele A with probability p_m*
Sampling, kSib > 1:             family = one mother (breeding female) + one father (breeding male);
                                daughters get one random maternal and one paternal allele;
                                haploid sons get one maternal allele
Heterozygosity retained:        H(last point)/H(first point) = Π_g (1 − 1/(2·Ne(g)))
```

Selection starts in the generation after the onset year, so a sample taken at the onset year sits exactly at p₀.

---

## Using the simulation core from JavaScript / Node

`sim.js` contains the pure simulation functions, with no DOM code. They are the same functions embedded in `index.html`.

```js
const S = require('./sim.js');

const params = {
  system: 'haplodiploid',          // or 'diploid' or 'haploid'
  genTime: 1,                      // years per generation
  Nf: 1000, Nm: 1000,              // breeding females / males per generation
  noDrift: false,
  decline: 1,                      // breeding numbers at the last point relative to the first
  neByBlock: false, neBlocks: null,// optional: Ne for each interval between sampling points
  mig: 0, pMig: -1,                // migration fraction; source frequency (-1 = p0)
  p0: 0.01, s: 0.05, h: 0.5,       // starting frequency, selection, dominance
  sMale: 1,                        // selection in males, multiple of s
  onset: 1935,                     // year selection starts
  years: [1935, 1955, 1975, 2025], // sampling years (2 to 20 points)
  n: [20, 20, 20, 40],             // individuals per sampling point (same length as years)
  maleFrac: 0,                     // probability a sampled individual is male
  kSib: 1,                         // sampled individuals per family
  perfect: false,
  depthEarly: 10, depthLast: 15,   // mean depth, earlier points and last point
  err: 0.005, minDP: 3,
  alpha: 0.05, m: 1,
  reps: 1000, seed: 1
};

const r = S.runPower(params);
console.log(r.trend.powerEmp);     // power, trend test across all points, drift-aware threshold
console.log(r.two.powerEmp);       // power, first vs last point, drift-aware threshold
console.log(r.trend.fprNaive);     // naive false-positive rate
console.log(S.neFor(params));      // Ne from breeding numbers
console.log(1 - S.hetRetained(params)); // neutral heterozygosity lost
```

| Function | Returns |
|---|---|
| `runPower(params, {keep})` | `{ trend, two, alphaAdj, zNaive, sel, neu, trajS, trajN, hetAcc, hapAcc, missRate }`. `trend` and `two` each contain `crit`, `method` (`'empirical'` or `'gaussian'`), `powerEmp`, `powerNaive` and `fprNaive`. `sel` and `neu` hold the Z values. `keep` = number of trajectories to return. |
| `deterministic(params, both)` | The deterministic female frequency trajectory, or `{f, m}` for both sexes if `both` is true. |
| `neFor(params)` | Ne from `Nf`, `Nm` and `system` (`Nf` alone in the haploid system). |
| `sizeFactor(params, gens, g)` | Multiplier applied to the breeding numbers in generation `g` (exponential change or per-interval Ne). |
| `hetRetained(params)` | Expected fraction of neutral heterozygosity retained from the first to the last point. |
| `gensOf(params)` | Generation index of each sampling year. |

---

## Assumptions and limitations

- One biallelic, autosomal locus. There is no linkage, so hitchhiking and background selection are not modelled.
- Discrete, non-overlapping generations. Every generation is a breeding generation. For long-lived species with overlapping generations, the generation time and Ne absorb most of the effect. Note that an individual's alleles reflect the population at its birth, not at its sampling date, so for long-lived species the birth year (or an estimate of it) is the better date to enter.
- Selection is constant in strength from the onset year onwards and acts on viability before breeding.
- No population structure. A single panmictic population, apart from the optional migrant source.
- Population size is constant or changes exponentially. There is no random fluctuation.
- Families (`kSib > 1`) use one mother and one father, so there is no multiple paternity within a family.
- Genotype calling uses a simple per-read error model and a flat prior. Reference bias, mapping quality and damage patterns are not modelled explicitly, beyond the single error rate.
- Population size is either constant or changes exponentially between the first and last points. More complex histories, such as a bottleneck followed by recovery, need a dedicated forward simulator such as SLiM or simuPOP.
- Single-locus power only. Genome-wide scans would need a much stricter threshold. The Bonferroni m approximates testing several candidates.

---

## Files

| File | Content |
|---|---|
| `index.html` | Self-contained interactive simulator. Open it in a browser. |
| `sim.js` | Simulation core (pure functions), usable in Node or in the browser. |
| `README.md` | This documentation. |
| `SIMULATION_GUIDE.md` | Step-by-step explanation of the simulations, with worked examples. |
| `POSSIBLE_IMPROVEMENTS.md` | Ideas that are not implemented yet and may be added depending on feedback. |
| `examples/` | Short R scripts (selection only, drift only, selection + drift with the power calculation) and a Node script that reproduces the numbers in the guide. |
