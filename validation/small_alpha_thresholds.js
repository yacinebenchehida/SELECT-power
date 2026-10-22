// Small per-locus alpha (many candidate loci): compares ways of setting the drift-aware threshold when there are too
// few neutral replicates for an empirical quantile (alpha' x 2,000 < 10):
//   gaussian    : current SPOT rule, mean + SD x z(1 - alpha') of the 2,000 neutral Z values
//   exp_tail    : exponential tail fitted to the neutral Z values above their 95% quantile (peaks over threshold)
//   empirical_N : empirical quantile from enough neutral replicates (10 / alpha'), i.e. more simulation
// Realised false-positive rates are measured on 1,000,000 independent neutral replicates.
// Run from the repository root: node validation/small_alpha_thresholds.js > validation/small_alpha.tsv
const S = require('../sim.js');
const BASE = {system: 'haplodiploid', neByBlock: false, neBlocks: null, genTime: 1, Nf: 1000, Nm: 1000, noDrift: false,
  p0: 0.1, s: 0, h: 0.5, onset: 1935, years: [1935, 1955, 1975, 2025], n: [20, 20, 20, 40], perfect: false,
  depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3, alpha: 0.05, m: 1, sMale: 1, decline: 1, mig: 0, pMig: -1,
  kSib: 1, maleFrac: 0, reps: 1000, seed: 1};
const neutral = (P, n, seed) => { const r = S.runPower(Object.assign({}, P, {s: 0, reps: Math.ceil(n / 3), seed}));
  return {trend: Float64Array.from([...r.sel.trend, ...r.neu.trend]).sort(), two: Float64Array.from([...r.sel.two, ...r.neu.two]).sort()}; };
const empirical = (z, a) => z[Math.min(z.length - 1, Math.ceil((1 - a) * z.length) - 1)];
function gaussian(z, a) { let m = 0, v = 0; for (const x of z) m += x; m /= z.length; for (const x of z) v += (x - m) ** 2; v /= z.length - 1; return m + Math.sqrt(v) * S.qnorm(1 - a); }
function expTail(z, a) { const u = empirical(z, 0.05); let s = 0, k = 0; for (const x of z) if (x > u) { s += x - u; k++; }
  const beta = k ? s / k : 0, pu = k / z.length; return pu > 0 && beta > 0 ? u + beta * Math.log(pu / a) : u; }
const D = [['haplodiploid Ne=300', Object.assign({}, BASE, {Nf: 200, Nm: 200})], ['haplodiploid Ne=1500', BASE],
  ['haplodiploid Ne=100000', Object.assign({}, BASE, {Nf: 66667, Nm: 66667})], ['haplodiploid Ne=1500, p0=0.01', Object.assign({}, BASE, {p0: 0.01})]];
console.log(['design', 'statistic', 'alpha_per_locus', 'method', 'fpr_mean', 'fpr_sd', 'ratio_to_alpha', 'n_neutral_for_threshold'].join('\t'));
const t0 = Date.now();
for (const [lab, P] of D) {
  const test = neutral(P, 1000000, 424242);
  const small = [], big = {};
  for (let i = 0; i < 20; i++) small.push(neutral(P, 2000, 7000 + i));
  for (const a of [0.0005, 0.00005]) big[a] = [0, 1, 2, 3, 4].map(i => neutral(P, Math.ceil(10 / a), 9000 + i));
  for (const k of ['trend', 'two']) for (const a of [0.0005, 0.00005]) {
    const frac = c => { let n = 0; for (const x of test[k]) if (x > c) n++; return n / test[k].length; };
    const res = {gaussian: small.map(z => frac(gaussian(z[k], a))), exp_tail: small.map(z => frac(expTail(z[k], a))),
      empirical_N: big[a].map(z => frac(empirical(z[k], a)))};
    for (const meth in res) { const v = res[meth], m = v.reduce((x, y) => x + y, 0) / v.length, sd = Math.sqrt(v.reduce((x, y) => x + (y - m) ** 2, 0) / (v.length - 1));
      console.log([lab, k, a, meth, m.toExponential(3), sd.toExponential(2), (m / a).toFixed(2), meth === 'empirical_N' ? Math.ceil(10 / a) : 2000].join('\t')); }
  }
  console.error(lab, ((Date.now() - t0) / 1000).toFixed(0) + 's');
}
