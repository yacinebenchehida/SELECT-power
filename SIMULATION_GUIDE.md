# How the simulations work: a guide with worked examples

This guide explains, step by step and with real numbers, what SELECT (SELECTion power over Time) does under the hood: how allele frequencies are simulated, how specimens are sampled, how a power value is produced, and what the more technical pieces of the tool do. It complements the [README](README.md), which describes the parameters and the interface.

Every number in this guide was produced with the simulator's own code (`sim.js`) or with the short R scripts mentioned at the end. Settings are given with each example so that it can be reproduced.

---

## Contents

1. [The big picture](#1-the-big-picture)
2. [Two separate lotteries: drift and sampling](#2-two-separate-lotteries-drift-and-sampling)
3. [Building block A: sampling specimens](#3-building-block-a-sampling-specimens)
4. [Building block B: selection alone](#4-building-block-b-selection-alone)
5. [Building block C: drift alone](#5-building-block-c-drift-alone)
6. [Selection and drift together](#6-selection-and-drift-together)
7. [From replays to a power value](#7-from-replays-to-a-power-value)
8. [The "No drift" option](#8-the-no-drift-option)
9. [Why more specimens do not always help](#9-why-more-specimens-do-not-always-help)
10. [Why a drift-aware threshold is needed](#10-why-a-drift-aware-threshold-is-needed)
11. [The more technical pieces, explained intuitively](#11-the-more-technical-pieces-explained-intuitively)
12. [How each setting moves the power: reference numbers](#12-how-each-setting-moves-the-power-reference-numbers)
13. [Glossary](#13-glossary)
14. [Reproducing the examples](#14-reproducing-the-examples)

---

## 1. The big picture

A power analysis asks one question:

> **If selection of a given strength really acted on an allele, how often would a given sampling design detect it?**

The simulator answers by brute force. It replays the history of a population many times on the computer, collects virtual specimens at the sampling dates exactly as a real study would, scores each replay, and counts how often selection is detected.

```
 replay the population's history     catch specimens at the       score the       compare with replays    count
 (generation by generation)    ──►   sampling dates          ──►  observed   ──►  without selection   ──►  detections
 selection + drift                   (+ sequencing errors)        change          (the threshold)          = power
```

Two sets of replays are run, in the same world (same population size, same drift, same sampling), differing only in selection:

| Set | What it contains | Default number | Purpose |
|---|---|---|---|
| Without selection (s = 0) | drift + sampling | 2,000 | Shows how much apparent change happens **by chance**. Sets the threshold. |
| With selection (chosen s) | selection + drift + sampling | 1,000 | Shows what the data look like **when selection acts**. Power = share above the threshold. |

The rest of this guide builds this up piece by piece.

---

## 2. Two separate lotteries: drift and sampling

Two different kinds of chance are involved. They are easy to confuse, so it helps to keep them apart from the start.

| | **Drift** | **Sampling** |
|---|---|---|
| Who draws | Nature: which individuals breed | The study: which individuals end up in the tubes |
| How often | Every generation (e.g. 90 times between 1935 and 2025) | Only at the sampling dates (e.g. 4 times) |
| Size of each draw | 2N gene copies (the breeding population) | 2n gene copies (the specimens) |
| What it changes | The population's **true** allele frequency | Only the **estimate** of it |
| Removed by | An infinitely large population ("No drift") | Sequencing every individual |
| Reduced by | Larger Ne | More specimens |

**Selection is not a lottery.** It is a deterministic push: the same in every replay.

**Key consequence:** more specimens can shrink the sampling lottery, but they can never shrink the drift lottery.

---

## 3. Building block A: sampling specimens

### The idea

Suppose that in a given year the allele frequency in the population is *p*. Picking one gene copy at random from the population is like flipping a weighted coin that lands "allele" with probability *p*. Sampling *n* diploid individuals means flipping 2*n* such coins. The number of allele copies found follows a **binomial distribution**:

```
copies found ~ Binomial(2n, p)          sample frequency p̂ = copies found / 2n
```

The "hat" on p̂ marks a value **estimated from the specimens**, as opposed to the true population value *p*.

### Example 1: a rare allele, small sample

True frequency p = 0.01 (1%), 20 diploid individuals = 40 gene copies.

| Copies found | Probability | Sample frequency p̂ |
|---|---|---|
| 0 | 67% | 0% |
| 1 | 27% | 2.5% |
| 2 | 5% | 5% |
| 3 or more | 1% | ≥ 7.5% |

In two samples out of three, the allele is not seen at all, even though it is present in the population.

### Example 2: a common allele, larger sample

True frequency p = 0.16, 40 individuals = 80 gene copies. On average 12.8 copies are found, but individual samples typically give anywhere from about 7 to 19 copies (p̂ between 9% and 24%).

**Take-away:** the same population, sampled twice, gives two different p̂ values. This is sampling noise. It exists even when nothing at all is changing in the population.

---

## 4. Building block B: selection alone

### The idea

With selection and an infinitely large population (no drift), the allele frequency follows a fixed, predictable curve. Each generation:

1. Offspring genotypes form from the parents' allele frequency (Hardy–Weinberg proportions).
2. Genotypes have different fitnesses: 1+s (two copies), 1+hs (one copy), 1 (none).
3. The allele frequency among the surviving breeders is recomputed. That is the next generation's frequency.

No randomness is involved. Every replay follows exactly the same curve.

### Example 3: the selection curve

Diploid, additive (h = 0.5), s = 0.05, starting at p₀ = 0.01, 90 generations:

| Generation | 0 | 45 | 90 |
|---|---|---|---|
| True frequency | 0.0100 | 0.0298 | 0.0851 |

### Example 4: sampling from the selection curve

Now 20, 20 and 40 individuals are sampled at generations 0, 45 and 90 (40, 40 and 80 gene copies). The population curve is identical in every replay; only the specimens change. Five replays:

| Replay | Copies at gen 0 (of 40) | Copies at gen 45 (of 40) | Copies at gen 90 (of 80) |
|---|---|---|---|
| 1 | 0 | 1 | 7 |
| 2 | 1 | 0 | 10 |
| 3 | 2 | 1 | 8 |
| 4 | 0 | 0 | 4 |
| 5 | 1 | 1 | 9 |

Each sampling point is drawn **independently** from its own known frequency: Binomial(40, 0.0100), Binomial(40, 0.0298), Binomial(80, 0.0851). What was caught at one date has no influence on the next date, because catching specimens does not change the population.

Replay 4 looks much less convincing than replay 2, purely through sampling luck.

---

## 5. Building block C: drift alone

### The idea

In a finite population, the next generation's breeders carry a random draw of the current generation's gene copies. With N breeding diploid individuals (2N gene copies):

```
copies in next generation ~ Binomial(2N, p_current)          p_next = copies / 2N
```

On average p_next equals p_current, but it is never exactly equal. The gap is **drift**.

Crucially, **each generation is drawn from the previous generation**, not from the starting value. The frequency therefore performs a random walk: small random steps that accumulate.

```
p₁ = Bin(2N, p₀)/2N  →  p₂ = Bin(2N, p₁)/2N  →  p₃ = Bin(2N, p₂)/2N  →  …
```

This is why, with drift, the whole population has to be simulated generation by generation, even when only a few dates are sampled.

### How big are the steps?

One generation's step has a standard deviation of √(p(1−p) / 2N). The steps add up, so after t generations the typical distance from the start is roughly:

```
SD after t generations ≈ √( p(1−p) · t / 2Ne )
```

The spread grows with the **square root** of time: 100 generations give about 10 times the spread of a single generation, not 100 times.

### Example 5: step sizes for p = 0.30

| Ne (2Ne gene copies) | Typical step per generation | Typical distance after 100 generations |
|---|---|---|
| 10,000 (20,000) | ± 0.3 points | ± 3 points (most runs end between 0.24 and 0.36) |
| 100 (200) | ± 3 points | ± 30 points (many runs hit 0 or 1 and stop) |

### Example 6: five drift-only replays

Diploid, N = 500 breeding individuals, no selection, p₀ = 0.10, sampling 20/20/40 individuals at generations 0, 45, 90:

| Replay | True p at gen 0 | True p at gen 45 | True p at gen 90 | Copies found (40 / 40 / 80) |
|---|---|---|---|---|
| 1 | 0.100 | 0.118 | 0.183 | 3 / 5 / 14 |
| 2 | 0.100 | 0.166 | 0.189 | 2 / 3 / 16 |
| 3 | 0.100 | 0.067 | 0.048 | 6 / 3 / 3 |
| 4 | 0.100 | 0.101 | 0.122 | 3 / 5 / 12 |
| 5 | 0.100 | 0.129 | 0.099 | 3 / 5 / 7 |

Now the **true** frequencies differ between replays (drift), and the sampled copies add sampling noise on top. Replays 1 and 2 show an apparent rise from 10% to about 18–19% **without any selection**. A test that ignored drift would be fooled by them.

**Mental image:** a fan of lines, all starting at the same point, wiggling and spreading apart over time. Small Ne gives a wide, ragged fan; large Ne gives a narrow one.

---

## 6. Selection and drift together

### The idea

Each generation, two things happen in order:

```
p (this generation)
   │
   ▼  1. selection: a formula, no randomness  →  pushes p up slightly to p*
   │
   ▼  2. drift: binomial draw  →  Binomial(2N, p*)
p (next generation) = copies / 2N
```

**Selection tilts the coin; drift is the luck of flipping a limited number of coins.** Selection is not itself a binomial draw: it only sets the probability used in the drift draw.

A useful picture is a boat on a river: the current (selection) carries it steadily downstream, while the waves (drift) knock it randomly left and right.

### Example 7: three generations traced by hand

Diploid, Ne = 10,000 (20,000 gene copies), p₀ = 0.30, s = 0.05, additive:

| Generation | Start p | After selection p* (formula) | Drift draw: Binomial(20,000, p*) | Next p |
|---|---|---|---|---|
| 0 → 1 | 0.30000 | 0.30517 | 6,094 copies | 0.30470 |
| 1 → 2 | 0.30470 | 0.30992 | 6,132 copies | 0.30660 |
| 2 → 3 | 0.30660 | 0.31183 | 6,123 copies | 0.30615 |

Selection always pushes up by about +0.5 points; drift lands a little above or below (in the last row it even fell back). In this replay the population reached 0.606 at generation 50 and 0.845 at generation 100. Sampling 20, 20 and 40 individuals at generations 0, 50 and 100 then gave 8/40, 26/40 and 65/80 copies.

### Who wins: push or waves?

The push accumulates in a straight line over time; the waves partly cancel and grow only like √time. So, given enough generations, selection dominates, unless the waves are large compared with the push.

### Example 8: the same selection under two population sizes

p₀ = 0.30, s = 0.05 additive, 100 generations, 2,000 replays each (deterministic expectation: 0.83):

| Ne | Push per generation | Drift step | Where replays end |
|---|---|---|---|
| 10,000 | ≈ +0.5 points | ≈ ± 0.3 points | 95% between 0.79 and 0.87: a narrow fan hugging the selection curve |
| 100 | ≈ +0.5 points | ≈ ± 3 points | Anywhere between 0 and 1; median still 0.84, but 3% lose the allele and 10% end below 0.35, looking as if nothing happened |

**Rule of thumb:** selection dominates when Ne × s is well above 1 (here 500 versus 5).

**Rare alleles are most at risk.** The push is proportional to p(1−p), so at p = 1% it is tiny. For a while, a rare favoured allele is at the mercy of drift and can be lost before selection takes hold.

---

## 7. From replays to a power value

This is the part that turns simulations into one number.

### The score: what a real study would see

Each replay is summarised by a single score computed **only from the specimens** (the sample frequencies p̂), never from the population's true frequencies or its path. A real study cannot see the population, so the score must not see it either.

The simplest possible score is the **observed change**:

```
change = p̂(last date) − p̂(first date)
```

The tool uses a slightly better score that also uses the middle dates and the sample sizes (the trend score, section 11.1). The logic below is identical.

### The four steps

1. **Run the no-selection replays** (default 2,000) and compute each one's score.
2. **Set the threshold:** sort the scores from smallest to largest, set aside the top 5% (the 100 largest), and take the largest remaining value: the 1,900th. A score above the threshold counts as "selection detected". By construction, this rule is fooled in 5% of cases when nothing is happening: that is the significance level α = 0.05.
3. **Run the selection replays** (default 1,000) in the same world (same drift, same sampling, plus selection) and compute their scores.
4. **Count:** power = share of selection replays whose score is strictly above the threshold.

### Example 9: the whole calculation by hand (small version)

Diploid, Ne = 10,000, p₀ = 0.30, 100 generations, 20 individuals at generation 0 (40 copies) and 40 at generation 100 (80 copies). Score = observed change. To make it readable: 20 no-selection replays and 10 selection replays (s = 0.01).

**Step 1. Twenty no-selection replays, observed changes sorted:**

```
−0.100  −0.100  −0.087  −0.075  −0.075  −0.050  −0.050  −0.025  +0.012  +0.013
+0.013  +0.025  +0.037  +0.038  +0.050  +0.062  +0.062  +0.075  +0.088  +0.112
```

Nothing happened in these populations, yet the observed change ranges from −10 to +11 points: drift plus sampling luck.

**Step 2. The threshold.** 5% of 20 replays is one replay. The threshold is placed so that only one neutral replay is above it: **+0.088** (only +0.112 clears it).

**Step 3. Ten selection replays (s = 0.01):**

| Replay | True p at gen 100 | Sample p̂ (gen 0 → 50 → 100) | Observed change | Above +0.088? |
|---|---|---|---|---|
| 1 | 0.42 | 0.275 → 0.425 → 0.438 | +0.162 | yes |
| 2 | 0.44 | 0.300 → 0.250 → 0.425 | +0.125 | yes |
| 3 | 0.40 | 0.375 → 0.375 → 0.425 | +0.050 | **no (missed)** |
| 4 | 0.45 | 0.250 → 0.375 → 0.463 | +0.213 | yes |
| 5 | 0.41 | 0.300 → 0.300 → 0.537 | +0.237 | yes |
| 6 | 0.39 | 0.425 → 0.250 → 0.438 | +0.013 | **no (missed)** |
| 7 | 0.43 | 0.225 → 0.250 → 0.463 | +0.238 | yes |
| 8 | 0.42 | 0.300 → 0.275 → 0.438 | +0.138 | yes |
| 9 | 0.48 | 0.325 → 0.400 → 0.450 | +0.125 | yes |
| 10 | 0.42 | 0.300 → 0.400 → 0.525 | +0.225 | yes |

**Step 4. Count:** 8 of 10 clear the threshold, so **power = 80%**.

In replays 3 and 6 the population really rose (to 0.40 and 0.39), but the specimens caught at generation 0 happened to over-represent the allele, so the observed change looked no larger than what neutral replays produce by chance. Those are the misses.

With only 20 and 10 replays, this estimate is rough. The tool does the same with 2,000 and 1,000 replays, which pins it down much better.

### Example 10: power for several selection strengths

Same setting as Example 9 but with the full trend score and thousands of replays:

| s | True rise over 100 generations | Power |
|---|---|---|
| 0.002 | 0.30 → 0.32 | 9% |
| 0.005 | 0.30 → 0.36 | 17% |
| 0.01 | 0.30 → 0.41 | 31% |
| 0.05 | 0.30 → 0.83 | ≈ 100% |

(Example 9's 80% at s = 0.01 is higher than the 31% here only because 20 neutral replays give a very imprecise threshold. This is exactly why the tool uses thousands of replays.)

### Why 2,000 replays for the threshold but 1,000 for selection?

- The threshold sits in the **tail** of the neutral distribution, where values are sparse. With 2,000 replays, 100 of them lie above the 5% cut, enough to place it reliably. With 1,000 there would be only 50, and the threshold would wobble more from run to run.
- The selection set only estimates a **percentage**. With 1,000 replays, a power of 80% is known to about ± 1.3 points: √(0.8 × 0.2 / 1,000).
- Together, 3,000 simulated populations take about 0.1–0.4 s in the browser.

The totals can be changed with the "Replicates under selection" menu (500, 1,000, 2,000 or 5,000; the no-selection set is always twice that).

---

## 8. The "No drift" option

### What it does

"No drift" means the population is treated as so large that it never wanders by chance. The drift draw is skipped: each generation, p_next = p* exactly. It applies to **both** sets of replays at once:

| | No-selection replays | Selection replays |
|---|---|---|
| Default (drift on) | drift + sampling | selection + drift + sampling |
| "No drift" ticked | population frozen at p₀; sampling only | population follows the fixed selection curve; sampling only |

The specimens still come from the population: "infinite" only means that taking a few dozen individuals, or one generation of breeding, does not change its frequency. In practice, a population of tens of thousands behaves almost identically.

### How the threshold is set without drift

The recipe does not change; only the no-selection replays become simpler.

### Example 11: the threshold and power without drift

p₀ = 0.30; 20 individuals at the start (40 copies), 40 at the end (80 copies).

1. One no-selection replay: start Binomial(40, 0.30), say 11 copies (p̂ = 0.275); end Binomial(80, 0.30), the same 0.30 since nothing changed, say 27 copies (p̂ = 0.338). Observed change: +0.063, purely sampling luck.
2. Repeat 2,000 times, sort, set aside the top 5%: the threshold lands at about **+0.15**. Two samples this small can differ by about 15 points even when the population has not moved at all.
3. Selection replays (s = 0.01) follow the fixed curve 0.30 → 0.41: start Binomial(40, 0.30), end Binomial(80, 0.41). Their observed changes scatter around +0.11, and about **35%** clear +0.15: the power.

With drift switched back on, the no-selection population also wanders between the two dates, so the threshold rises and power drops.

### Example 12: the four combinations

The clearest way to separate the two lotteries is to switch each one off in turn. Diploid, p₀ = 0.30, s = 0.01, perfect genotypes, four sampling points:

| | Sampling: 20/20/20/40 individuals | Sampling: 500 per point (close to a census) |
|---|---|---|
| **No drift** | 34% (sampling noise only) | 100% (no noise left) |
| **Drift, Ne = 10,000** | 32% (both lotteries) | 85% (drift only) |

Reading the table:

- Top left to top right: more specimens remove the sampling noise, and power reaches 100%.
- Bottom right: even with 500 specimens per point, drift alone keeps power at 85%. **Specimens cannot remove drift.**
- Left column: at this Ne, drift adds little on top of the sampling noise of 20–40 specimens.

---

## 9. Why more specimens do not always help

### The variance split

The observed change mixes two independent sources of noise. Approximately:

```
noise² in the observed change ≈ Var_drift + Var_sampling(start) + Var_sampling(end)

Var_drift     ≈ p(1−p) · t / (2·Ne)        ← only Ne reduces it
Var_sampling  ≈ p(1−p) / (2n)              ← only more specimens reduce it
```

Power depends on how big the true change is compared with this total noise. When the drift term dominates, adding specimens shrinks a term that was already small.

### Example 13: drift versus sampling noise

p ≈ 1%, 90 generations, 40 specimens at the end (sampling noise computed at the final frequency, about 16%):

| Ne | Drift noise (SD) | Sampling noise at the end (SD) |
|---|---|---|
| 200 | ≈ 4.7 points (allele often lost) | ≈ 4 points |
| 2,000 | ≈ 1.5 points | ≈ 4 points |
| 20,000 | ≈ 0.5 points | ≈ 4 points |

### Example 14: the planner's message

| Design (individuals at the 4 points) | Default settings (haplodiploid, Ne ≈ 1,500, p₀ = 0.01, s = 0.05) | Diploid, Ne = 10,000, p₀ = 0.30, s = 0.01 |
|---|---|---|
| 10 / 10 / 10 / 20 | 62% | 20% |
| 20 / 20 / 20 / 40 | 67% | 29% |
| 40 / 40 / 40 / 80 | 70% | 45% |
| 80 / 80 / 80 / 160 | 73% | 59% |

With a rare starting allele and moderate Ne (left), quadrupling the sample barely helps: drift and the rarity of the allele are the limit. With a common allele and large Ne (right), sampling is the limit, and more specimens help a lot.

---

## 10. Why a drift-aware threshold is needed

A standard ("naive") test assumes that the only noise is sampling noise. It uses a fixed textbook threshold (1.645 for a one-sided test at 5%). When drift is present, neutral alleles wander more than that test expects, so it calls selection too often.

### Example 15: how often the naive test is fooled

Haplodiploid, p₀ = 0.30, s = 0.02, default sampling, equal numbers of breeding females and males:

| Breeding females = males | Ne | Naive test: false positives (should be 5%) | Power the naive test reports | Real (drift-aware) power |
|---|---|---|---|---|
| 100 | 150 | 28% | 63% | 24% |
| 300 | 450 | 20% | 77% | 50% |
| 1,000 | 1,500 | 12% | 86% | 74% |
| 3,000 | 4,500 | 8% | 92% | 87% |
| 30,000 | 45,000 | 6% | 95% | 95% |

In small populations, the naive test both raises false alarms and overstates power. The drift-aware threshold, taken from replays that contain drift, fixes both. The two agree only when Ne is large.

---

## 11. The more technical pieces, explained intuitively

### 11.1 The trend score (Cochran–Armitage test)

**Intuition.** Plot the sample frequency at each sampling date against time and draw the best straight line through the points. The score is the **slope of that line divided by how uncertain the slope is**. A steep, consistent rise based on many specimens gives a high score; a flat or zig-zag pattern gives a score near 0; a fall gives a negative score. Dates with more specimens weigh more.

**Formula.** At sampling point i (time tᵢ in generations), xᵢ copies of the allele are found among nᵢ gene copies:

```
Z = Σᵢ xᵢ (tᵢ − t̄) / √( p̄(1 − p̄) Σᵢ nᵢ (tᵢ − t̄)² )
t̄ = Σ nᵢ tᵢ / Σ nᵢ        p̄ = Σ xᵢ / Σ nᵢ
```

**Example 16: computing a score by hand.** Copies found: 11 of 40 at t = 0, 20 of 40 at t = 50, 50 of 80 at t = 100.

- Pooled frequency p̄ = 81 / 160 = 0.506; weighted mean time t̄ = (40·0 + 40·50 + 80·100) / 160 = 62.5.
- Numerator: 11·(0 − 62.5) + 20·(50 − 62.5) + 50·(100 − 62.5) = −687.5 − 250 + 1,875 = **937.5**.
- Denominator: √(0.506 · 0.494 · [40·62.5² + 40·12.5² + 80·37.5²]) = √(0.25 · 275,000) ≈ **262.2**.
- **Z ≈ 3.58**, far above the usual thresholds (1.6–2.0): a clear rise.

**The two-point version** uses only the first and last dates. With the same data (11/40 → 50/80) it gives Z ≈ 3.62, identical to the classic two-proportion z-test. The four-point version helps most when selection acts steadily from early on; when selection starts late, the middle dates add little (Example 22).

### 11.2 The empirical null with real data

In the simulator, the no-selection replays define the threshold. With real data, the equivalent is the set of **genome-wide SNPs**: most of them are neutral, they experienced the same drift, and they were sequenced in the same specimens with the same errors. Computing the same score for thousands of SNPs with a similar starting frequency gives a realistic distribution of chance changes. A candidate locus is then judged against the top 5% of that distribution. This handles drift, population structure and sequencing artefacts together, without needing an estimate of Ne.

### 11.3 Testing several candidate loci (Bonferroni)

Testing m loci at 5% each would produce false alarms much more often than 5% overall. The Bonferroni correction tests each locus at 5% / m instead, which raises the threshold.

**Example 17** (default settings): m = 1 → threshold 1.63, power 67%; m = 10 → threshold 2.86, power 50%; m = 100 → threshold 3.16, power 44%.

### 11.4 The bell-curve shortcut (Gaussian fallback)

With a very strict cut, too few neutral replays remain above it to place the threshold reliably (with 2,000 replays and α/m = 0.0005, only one would). The tool then fits a bell curve (normal distribution) to the 2,000 neutral scores and reads the threshold from the curve. The status line says when this happens, and the parameter power curve marks such points with hollow markers. In Example 17, m = 100 uses this shortcut.

### 11.5 Effective population size Ne

**Intuition.** Ne is the size of an idealised population that would drift as much as the real one. It is usually much smaller than the census size, because only some individuals breed, and the rarer sex limits it.

- Diploid: **Ne = 4·Nf·Nm / (Nf + Nm)**. With N females and N males, Ne = 2N.
- Haplodiploid: **Ne = 9·Nf·Nm / (4·Nm + 2·Nf)**. With N of each, Ne = 1.5N. Males carry only one copy, so they contribute fewer gene copies to the next generation and the population drifts more for the same head count.
- A very skewed breeding sex ratio reduces Ne strongly: drift is dominated by the rarer sex.

### 11.6 Heterozygosity loss: anchoring Ne with real data

Drift steadily erodes neutral genetic diversity:

```
H_t / H₀ = (1 − 1/(2Ne))^t
Ne = 1 / ( 2 · (1 − (H_t/H₀)^(1/t)) )
```

**Example 18:** heterozygosity lost over 100 generations:

| Ne | 50 | 100 | 250 | 500 | 1,000 | 2,500 | 5,000 | 10,000 |
|---|---|---|---|---|---|---|---|---|
| Heterozygosity lost | 63% | 39% | 18% | 9.5% | 4.9% | 2.0% | 1.0% | 0.5% |

Measuring genome-wide heterozygosity in the first and last samples therefore gives a data-based estimate of Ne, which can then be entered in the simulator. A loss of about 25% over 90 generations, for instance, corresponds to Ne ≈ 160.

### 11.7 Dominance and haplodiploidy

**Intuition.** In a diploid, a rare allele sits almost entirely in heterozygotes. If it is recessive (h = 0), heterozygotes look normal, so selection barely sees it. In a haplodiploid, males are haploid: every allele a male carries is expressed, so even a recessive allele is exposed to selection through males.

**Example 19:** default settings (p₀ = 0.01, s = 0.05, 90 generations):

| Dominance h | Final frequency, haplodiploid | Final frequency, diploid | Power, haplodiploid | Power, diploid |
|---|---|---|---|---|
| 0 (recessive) | 0.045 | 0.011 | 29% | 7% |
| 0.5 (additive) | 0.160 | 0.085 | 67% | 51% |
| 1 (dominant) | 0.365 | 0.333 | 89% | 95% |

(The diploid can do slightly better for dominant alleles because, with the same breeding numbers, its Ne is larger: 2N versus 1.5N.)

Setting selection in males to 0 (a trait not expressed in males) removes this advantage: power drops from 67% to 26% with default settings.

### 11.8 Genotype calling at low coverage

**Intuition.** Each specimen is "sequenced": the number of reads covering the locus is random (Poisson with the chosen mean depth), and each read shows one of the two gene copies, with a small error rate. With few reads, a heterozygote can by chance show only one allele and be called homozygous (**allelic dropout**). Specimens with fewer reads than the minimum depth are dropped.

**Example 20:** default settings, minimum depth 3, same mean depth at all points:

| Mean depth | Heterozygotes called correctly | Specimens dropped | Power |
|---|---|---|---|
| 2× | 82% | 68% | 50% |
| 3× | 85% | 42% | 59% |
| 6× | 93% | 6% | 64% |
| 10× | 97% | 0.3% | 68% |
| 20× | 99.9% | 0% | 66% |

Haploid males have no heterozygotes, so their calls are essentially always correct; but low depth still drops the same share of specimens. Above about 6–10×, more depth no longer changes power: drift and sample size are then the limits.

### 11.9 Related specimens (siblings)

Siblings share parents, so their alleles are correlated: in haplodiploids, sisters share their father's allele exactly. Ten sisters therefore carry less independent information than ten unrelated females, and the sample behaves as if it were smaller.

**Example 21** (default settings): 1 per family → 67%; 2 → 60%; 5 → 54%; 10 → 48%. The naive false-positive rate also rises (from 4.5% to about 8%).

### 11.10 Population decline, migration and timing

**Example 22** (default settings unless stated):

| Change | Power |
|---|---|
| Constant population (reference) | 67% |
| Breeding numbers fall to 30% by the last date | 58% (heterozygosity lost: 5.6% instead of 3.0%) |
| Breeding numbers fall to 10% by the last date | 44% (heterozygosity lost: 10.9%) |
| 1% migrants per generation (source at p₀) | 54% |
| 5% migrants per generation (source at p₀) | 11% |
| Selection starts in 1975 instead of 1935 | 29% (two-point test: 26%) |
| Selection starts in 2000 | 15% |

- **Decline:** drift accelerates as the population shrinks, right when the last sample is taken.
- **Migration** from a source where the allele is not favoured keeps pulling the frequency back towards p₀.
- **Late onset:** fewer generations of selection mean a smaller change; samples taken before the onset only measure the baseline.

### 11.11 Sampling haploid males (haplodiploids only)

**Example 23** (default settings, haplodiploid):

| Design | Gene copies | Power |
|---|---|---|
| Females only (20/20/20/40) | 200 | 67% |
| Males only, same number of individuals | 100 | 61% |
| Males only, same number of gene copies (40/40/40/80) | 200 | 68% |

A male gives one copy instead of two, so per specimen males carry less information. Per gene copy they are equivalent, with cleaner genotype calls. The main potential advantage of males, directly observed haplotypes without phasing, is not modelled here.

---

## 12. How each setting moves the power: reference numbers

Default settings (haplodiploid, 1,000 breeding females and 1,000 breeding males, Ne = 1,500, p₀ = 0.01, s = 0.05, h = 0.5, sampling 20/20/20/40 at 1935/1955/1975/2025, depth 10×/15×) give a four-point power of **67%** and a two-point power of **62%**. Switching only to diploid gives 51% (two-point 50%).

| Setting changed | New power | Section |
|---|---|---|
| Sample size ×0.5 / ×2 / ×4 | 62% / 70% / 73% | 9 |
| Recessive (h = 0) / dominant (h = 1) | 29% / 89% | 11.7 |
| No selection in males | 26% | 11.7 |
| Mean depth 2× / 3× | 50% / 59% | 11.8 |
| 5 siblings per family | 54% | 11.9 |
| Population falls to 10% | 44% | 11.10 |
| 5% migrants per generation | 11% | 11.10 |
| Selection starts in 1975 | 29% | 11.10 |
| 10 candidate loci (Bonferroni) | 50% | 11.3 |

---

## 13. Glossary

| Term | Plain meaning |
|---|---|
| Allele frequency (p) | Share of gene copies in the population that carry the allele. |
| Sample frequency (p̂) | Share of gene copies carrying the allele among the specimens: an estimate of p. |
| Binomial draw | Counting successes among a fixed number of independent weighted coin flips. |
| Selection coefficient (s) | Relative advantage per generation of carrying the allele (0.05 = 5% more offspring). |
| Dominance (h) | How much of the advantage a heterozygote gets: 0 = none (recessive), 0.5 = half (additive), 1 = full (dominant). |
| Drift | Random change in allele frequency caused by the finite number of breeders each generation. |
| Effective population size (Ne) | Size of an idealised population that would drift as much as the real one. |
| Replay | One simulated history of the population, from the first to the last sampling date. |
| Neutral / no-selection replay | A replay with s = 0: only chance changes the frequency. |
| Score (Z) | One number per replay summarising how clearly the sampled frequencies rise over time. |
| Threshold | The score that only 5% (α) of no-selection replays exceed. Other texts may call it the critical value or cut-off. |
| Significance level (α) | Accepted rate of false alarms: how often the threshold is cleared when nothing is happening. |
| Power | Share of selection replays whose score clears the threshold. |
| False positive | Calling selection when the allele is neutral. |
| Naive test | A test that assumes only sampling noise and ignores drift. |
| Empirical null | Distribution of scores from neutral loci (real data) or neutral replays (simulations), used to set the threshold. |
| Allelic dropout | A heterozygote called homozygous because, by chance, only one allele was read. |

---

## 14. Reproducing the examples

All scripts are in the [`examples/`](examples/) folder.

| Script | What it does | Examples |
|---|---|---|
| `selection_only_sampling.R` | Selection without drift: one deterministic trajectory, then repeated sampling of specimens from it. Plain loops, so every step is visible. | 3–4 |
| `drift_only_sampling.R` | Drift without selection: the population is re-simulated in every replay, then sampled. | 6 |
| `selection_drift_power.R` | Selection and drift together, then the full power calculation (no-selection replays → threshold → selection replays → count), with the observed change as the score. | 7–10 (same logic) |
| `guide_examples.js` | Runs the simulator's own functions (`sim.js`) for every other number in this guide. Run from the repository root with `node examples/guide_examples.js`. | 11, 12, 14, 15, 17 and 19–23, and the reference table |

Notes:

- Examples 1, 2, 5, 13, 16 and 18 are direct calculations from the formulas given next to them, and need no script.
- `guide_examples.js` uses the default parameters listed in the README, changing only the settings stated in each example, with 2,000 replicates under selection and 4,000 neutral (seed 1). Small differences (± 1–2 points) from the web page are expected, because the page uses 1,000 and 2,000 replicates by default.
- Examples 7–10 were first worked out with a minimal script of the same model; `selection_drift_power.R` implements that logic. Its numbers differ somewhat from the guide because random draws differ, and because it scores replays by the simple observed change rather than the trend score used for Example 10.
