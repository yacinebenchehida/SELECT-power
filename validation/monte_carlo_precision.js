// Monte Carlo precision: how much the power estimate varies between runs (different random seeds) for the example
// scenarios of the page, at 500, 1,000 (default), 2,000 and 5,000 replicates under selection (twice as many neutral).
// Run from the repository root: node validation/monte_carlo_precision.js > validation/mc_precision.tsv
const S = require('../sim.js');
const fs = require('fs');
const D = {system: 'haplodiploid', neByBlock: false, neBlocks: null, genTime: 1, Nf: 1000, Nm: 1000, noDrift: false,
  p0: 0.01, s: 0.05, h: 0.5, onset: 1935, years: [1935, 1955, 1975, 2025], n: [20, 20, 20, 40], perfect: false,
  depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3, alpha: 0.05, m: 1, sMale: 1, decline: 1, mig: 0, pMig: -1,
  kSib: 1, maleFrac: 0, reps: 1000, seed: 1};
// the example scenarios of the page (same values as PRESETS in parts/ui.js)
const SC = [
  ['Bumblebees', D],
  ['Long-lived mammal', Object.assign({}, D, {system: 'diploid', genTime: 20, Nf: 2500, Nm: 2500, neByBlock: true, neBlocks: [5000, 200], p0: 0.1, s: 0.2, onset: 1880, years: [1880, 1950, 2020], n: [20, 20, 40], depthEarly: 4, depthLast: 15})],
  ['Human ancient DNA', Object.assign({}, D, {system: 'diploid', genTime: 28, Nf: 5000, Nm: 5000, p0: 0.1, s: 0.005, onset: -7000, years: [-7000, -5000, -3000, -1000, 1000, 2000], n: [50, 50, 50, 50, 50, 100], depthEarly: 2, depthLast: 20, err: 0.01, minDP: 1})],
  ['Bacterial pathogen', Object.assign({}, D, {system: 'haploid', genTime: 0.01, Nf: 100000, Nm: 100000, p0: 0.05, s: 0.0002, onset: 1850, years: [1850, 1900, 1950, 2000, 2020], n: [10, 10, 15, 20, 40], depthEarly: 10, depthLast: 30})],
  ['Large diploid population', Object.assign({}, D, {system: 'diploid', genTime: 1, Nf: 1e6, Nm: 1e6, p0: 0.05, s: 0.02, onset: 1900, years: [1900, 1950, 2000, 2020], n: [30, 30, 30, 60]})],
];
const SEEDS = 40, REPS = [500, 1000, 2000, 5000];
const sd = a => { const m = a.reduce((x, y) => x + y, 0) / a.length; return Math.sqrt(a.reduce((x, y) => x + (y - m) ** 2, 0) / (a.length - 1)); };
const all = {};
console.log(['scenario', 'replicates', 'statistic', 'mean_power', 'sd_power', 'sd_binomial_only', 'sd_ratio', 'min', 'max', 'seeds'].join('\t'));
const t0 = Date.now();
for (const [lab, P] of SC) {
  all[lab] = {};
  for (const R of REPS) {
    const v = {trend: [], two: []};
    for (let i = 0; i < SEEDS; i++) { const r = S.runPower(Object.assign({}, P, {reps: R, seed: 5000 + 97 * i})); v.trend.push(r.trend.powerEmp); v.two.push(r.two.powerEmp); }
    all[lab][R] = v;
    for (const k of ['trend', 'two']) {
      const a = v[k], m = a.reduce((x, y) => x + y, 0) / a.length, b = Math.sqrt(m * (1 - m) / R);
      console.log([lab, R, k, m.toFixed(4), sd(a).toFixed(4), b.toFixed(4), (sd(a) / b).toFixed(2), Math.min(...a).toFixed(3), Math.max(...a).toFixed(3), SEEDS].join('\t'));
    }
    console.error(lab, R, ((Date.now() - t0) / 1000).toFixed(0) + 's');
  }
}
fs.writeFileSync('validation/mc_precision.json', JSON.stringify(all));
