// Compares SPOT with SLiM (validation/slim.md). Run from the repository root after validation/slim_runs.sh:
//     node validation/slim_compare.js        -> validation/slim_results.tsv (frequencies), slim_power.tsv (power) and slim_fpr.tsv (false positives)
//
// For each scenario of slim_runs.sh, SPOT simulates the same number of replicates (2,000) with the same
// parameters, using its own code (sim.js, simulateRep). Both simulators record, at generations 0, 20, 40 and 90:
//   pfs : allele frequency among breeding females after selection (haploid: all individuals)
//   x/n : allele copies among the sampled females (perfect genotypes)
// Then, separately for each simulator and with exactly the same code:
//   1. the distribution of pfs is summarised (mean, SD) and the two distributions are compared
//      (difference of means in standard errors, ratio of SDs, two-sample Kolmogorov-Smirnov test at generation 90);
//   2. the trend test Z (Cochran-Armitage across the four sampling dates, as in SPOT) is computed for every replicate;
//      the drift-aware threshold is the 95% quantile of Z in the neutral runs (as in SPOT), and power is the share of
//      selected replicates above it.
//   3. SPOT's threshold is also applied to SLiM's neutral replicates: if SPOT simulates drift correctly, this
//      false-positive rate should be close to 5%. For a precise version of this check (slim_fpr.tsv), SPOT's
//      threshold is taken from 20,000 neutral replicates and applied to all SLiM neutral replicates of the
//      scenario, including the extra runs _b, _c, _d of slim_runs.sh when present (up to 8,000).
const S = require('../sim.js');
const fs = require('fs');
const GENS = [0, 20, 40, 90], NS = [20, 20, 20, 40], R = 2000;

// same scenarios as slim_runs.sh: [selected run, matching neutral run, system, NF, NM, P0, S, H]
const SC = [
  ['diploid_N1000_s0.02', 'diploid_N1000_neutral', 'diploid', 500, 500, 0.1, 0.02, 0.5],
  ['diploid_N1000_s0.05_rec', 'diploid_N1000_neutral', 'diploid', 500, 500, 0.1, 0.05, 0],
  ['diploid_N200_s0.05', 'diploid_N200_neutral', 'diploid', 100, 100, 0.1, 0.05, 0.5],
  ['haploid_N1000_s0.02', 'haploid_N1000_neutral', 'haploid', 1000, 1000, 0.05, 0.02, 0.5],
  ['haplodiploid_N1000_s0.02', 'haplodiploid_N1000_neutral', 'haplodiploid', 500, 500, 0.1, 0.02, 0.5],
  ['haplodiploid_N2000_s0.05', 'haplodiploid_N2000_neutral', 'haplodiploid', 1000, 1000, 0.01, 0.05, 0.5],
];

// ---- reading SLiM output: columns replicate, generation, pfs, x, nall ----
function readSlim(name) {
  const reps = new Map();
  for (const line of fs.readFileSync(`validation/slim_out/${name}.tsv`, 'utf8').trim().split('\n')) {
    const [id, g, pfs, x, n] = line.split('\t').map(Number), k = GENS.indexOf(g);
    if (!reps.has(id)) reps.set(id, {pfs: [], x: [], n: []});
    const r = reps.get(id); r.pfs[k] = pfs; r.x[k] = x; r.n[k] = n;
  }
  return [...reps.values()];
}

// ---- SPOT: same scenario with its own simulation code ----
function runSpot(system, Nf, Nm, p0, s, h, seed) {
  const P = S.prep({system, neByBlock: false, neBlocks: null, genTime: 1, Nf, Nm, noDrift: false, p0, s, h, onset: 0,
    years: GENS, n: NS, perfect: true, depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3, alpha: 0.05, m: 1, sMale: 1,
    decline: 1, mig: 0, pMig: -1, kSib: 1, maleFrac: 0, reps: R, seed});
  const rng = S.makeRng(seed), out = [];
  for (let i = 0; i < R; i++) {
    const t = S.simulateRep(P, s, rng, true);
    out.push({pfs: GENS.map(g => t.traj[g]), x: t.x, n: t.n});
  }
  return out;
}

// ---- statistics (identical for both simulators) ----
function trendZ(x, n) {    // Cochran-Armitage trend statistic, as trendZ() in sim.js
  let N = 0, X = 0, T = 0;
  for (let i = 0; i < x.length; i++) { N += n[i]; X += x[i]; T += n[i] * GENS[i]; }
  const pb = X / N, tb = T / N; if (pb <= 0 || pb >= 1) return 0;
  let num = 0, den = 0;
  for (let i = 0; i < x.length; i++) { num += x[i] * (GENS[i] - tb); den += n[i] * (GENS[i] - tb) ** 2; }
  return num / Math.sqrt(pb * (1 - pb) * den);
}
const quant95 = z => { const a = [...z].sort((p, q) => p - q); return a[Math.min(a.length - 1, Math.ceil(0.95 * a.length) - 1)]; };
const mean = a => a.reduce((p, q) => p + q, 0) / a.length;
const sd = a => { const m = mean(a); return Math.sqrt(a.reduce((p, q) => p + (q - m) ** 2, 0) / (a.length - 1)); };
function ks(a, b) {        // two-sample Kolmogorov-Smirnov statistic D and asymptotic p-value
  const A = [...a].sort((p, q) => p - q), B = [...b].sort((p, q) => p - q);
  let i = 0, j = 0, D = 0;
  while (i < A.length && j < B.length) {
    const v = Math.min(A[i], B[j]);
    while (i < A.length && A[i] <= v) i++;
    while (j < B.length && B[j] <= v) j++;
    D = Math.max(D, Math.abs(i / A.length - j / B.length));
  }
  const ne = A.length * B.length / (A.length + B.length), lam = (Math.sqrt(ne) + 0.12 + 0.11 / Math.sqrt(ne)) * D;
  let p = 0; for (let k = 1; k <= 100; k++) p += 2 * (-1) ** (k - 1) * Math.exp(-2 * k * k * lam * lam);
  return [D, Math.min(1, Math.max(0, p))];
}

// ---- comparison ----
const rows = [['scenario', 'run', 'generation', 'mean_spot', 'mean_slim', 'z_mean', 'sd_spot', 'sd_slim', 'sd_ratio', 'ks_D', 'ks_p']];
const prow = [['scenario', 'system', 'NF', 'NM', 'p0', 's', 'h', 'threshold_spot', 'threshold_slim', 'power_spot', 'power_slim',
  'se_difference', 'z_power', 'fpr_spot_threshold_on_slim_neutral']];
const done = new Set(), neutralSpec = new Map();
SC.forEach(([sel, neu, system, Nf, Nm, p0, s, h], k) => {
  neutralSpec.set(neu, [system, Nf, Nm, p0, h]);
  const runs = {sel: [readSlim(sel), runSpot(system, Nf, Nm, p0, s, h, 100 + k)],
                neu: [readSlim(neu), runSpot(system, Nf, Nm, p0, 0, h, 200 + k)]};
  for (const [lab, name] of [['selected', sel], ['neutral', neu]]) {
    if (done.has(name)) continue; done.add(name);
    const [sl, sp] = runs[lab === 'selected' ? 'sel' : 'neu'];
    GENS.forEach((g, i) => {
      if (g === 0) return;
      const a = sp.map(r => r.pfs[i]), b = sl.map(r => r.pfs[i]);
      const z = (mean(a) - mean(b)) / Math.sqrt(sd(a) ** 2 / a.length + sd(b) ** 2 / b.length);
      const [D, p] = g === 90 ? ks(a, b) : ['', ''];
      rows.push([name, lab, g, mean(a).toFixed(4), mean(b).toFixed(4), z.toFixed(2), sd(a).toFixed(4), sd(b).toFixed(4),
        (sd(a) / sd(b)).toFixed(3), D === '' ? '' : D.toFixed(4), p === '' ? '' : p.toFixed(3)]);
    });
  }
  const Z = rs => rs.map(r => trendZ(r.x, r.n));
  const [zSelSlim, zSelSpot, zNeuSlim, zNeuSpot] = [Z(runs.sel[0]), Z(runs.sel[1]), Z(runs.neu[0]), Z(runs.neu[1])];
  const cSpot = quant95(zNeuSpot), cSlim = quant95(zNeuSlim);
  const pSpot = mean(zSelSpot.map(z => z > cSpot ? 1 : 0)), pSlim = mean(zSelSlim.map(z => z > cSlim ? 1 : 0));
  const se = Math.sqrt(pSpot * (1 - pSpot) / R + pSlim * (1 - pSlim) / R);
  const fpr = mean(zNeuSlim.map(z => z > cSpot ? 1 : 0));
  prow.push([sel, system, Nf, Nm, p0, s, h, cSpot.toFixed(3), cSlim.toFixed(3), pSpot.toFixed(4), pSlim.toFixed(4),
    se.toFixed(4), ((pSpot - pSlim) / se).toFixed(2), fpr.toFixed(4)]);
  console.error(sel, 'power SPOT', pSpot.toFixed(3), 'SLiM', pSlim.toFixed(3), 'FPR of SPOT threshold on SLiM neutral', fpr.toFixed(3));
});
fs.writeFileSync('validation/slim_results.tsv', rows.map(r => r.join('\t')).join('\n') + '\n');
fs.writeFileSync('validation/slim_power.tsv', prow.map(r => r.join('\t')).join('\n') + '\n');

// ---- precise check of the false-positive rate of SPOT's threshold on SLiM's neutral populations ----
const frow = [['neutral_scenario', 'system', 'NF', 'NM', 'p0', 'threshold_spot_20000', 'n_slim', 'fpr_on_slim', 'se', 'z_vs_5pc',
  'sd_Z_spot', 'sd_Z_slim', 'fpr_on_independent_spot']];
let k2 = 0;
for (const [neu, [system, Nf, Nm, p0, h]] of neutralSpec) {
  const files = [neu, neu + '_b', neu + '_c', neu + '_d'].filter(f => fs.existsSync(`validation/slim_out/${f}.tsv`));
  const zSlim = files.flatMap(f => readSlim(f).map(r => trendZ(r.x, r.n)));
  const big = []; for (let j = 0; j < 10; j++) big.push(...runSpot(system, Nf, Nm, p0, 0, h, 500 + 10 * k2 + j));
  const zSpot = big.map(r => trendZ(r.x, r.n)), c = quant95(zSpot);
  const test = []; for (let j = 0; j < 4; j++) test.push(...runSpot(system, Nf, Nm, p0, 0, h, 900 + 10 * k2 + j));
  const fprSpot = mean(test.map(r => trendZ(r.x, r.n) > c ? 1 : 0));
  const fpr = mean(zSlim.map(z => z > c ? 1 : 0)), se = Math.sqrt(0.05 * 0.95 / zSlim.length);
  frow.push([neu, system, Nf, Nm, p0, c.toFixed(3), zSlim.length, fpr.toFixed(4), se.toFixed(4), ((fpr - 0.05) / se).toFixed(2),
    sd(zSpot).toFixed(3), sd(zSlim).toFixed(3), fprSpot.toFixed(4)]);
  console.error(neu, 'SLiM neutral replicates', zSlim.length, 'FPR of SPOT threshold', fpr.toFixed(4), '(SPOT on itself', fprSpot.toFixed(4) + ')');
  k2++;
}
fs.writeFileSync('validation/slim_fpr.tsv', frow.map(r => r.join('\t')).join('\n') + '\n');
