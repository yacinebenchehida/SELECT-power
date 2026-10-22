// False-positive rate of the drift-aware test and of the naive test.
// For each design: SPOT sets its threshold from its own neutral replicates (20 independent runs, default 1,000 + 2,000
// replicates each), and the realised false-positive rate of each threshold is measured on a large, independent set of
// neutral replicates. The naive test uses the textbook threshold z(1 - alpha') on the same independent set.
// Run from the repository root:  node validation/false_positive_rate.js A > validation/fpr_A.tsv   (parts A, B, C)
// Part C can be split, e.g. `C 0,1` and `C 2,3`, to run on two processor cores.
const S = require('../sim.js');

const BASE = {system: 'haplodiploid', neByBlock: false, neBlocks: null, genTime: 1, Nf: 1000, Nm: 1000, noDrift: false,
  p0: 0.1, s: 0, h: 0.5, onset: 1935, years: [1935, 1955, 1975, 2025], n: [20, 20, 20, 40], perfect: false,
  depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3, alpha: 0.05, m: 1, sMale: 1, decline: 1, mig: 0, pMig: -1,
  kSib: 1, maleFrac: 0, reps: 1000, seed: 1};
const N_THRESHOLDS = 20;

// breeding numbers that give a chosen Ne in each system
function withNe(P, system, Ne) {
  const Q = Object.assign({}, P, {system});
  if (system === 'haploid') { Q.Nf = Ne; Q.Nm = Ne; }
  else if (system === 'diploid') { Q.Nf = Ne / 2; Q.Nm = Ne / 2; }
  else { Q.Nf = Ne / 1.5; Q.Nm = Ne / 1.5; }   // haplodiploid with Nf = Nm: Ne = 9N^2 / 6N = 1.5 N
  return Q;
}

// independent neutral Z values: every replicate of a run with s = 0 is neutral (1 + 2 per "rep")
function neutralZ(P, nWanted, seed) {
  const r = S.runPower(Object.assign({}, P, {s: 0, reps: Math.ceil(nWanted / 3), seed, alpha: 0.05, m: 1}));
  const out = {};
  for (const k of ['trend', 'two']) out[k] = Float64Array.from([...r.sel[k], ...r.neu[k]]);
  return out;
}

function evaluate(label, P, nTest) {
  const thr = [];
  for (let i = 0; i < N_THRESHOLDS; i++) thr.push(S.runPower(Object.assign({}, P, {seed: 1000 + i})));
  const test = neutralZ(P, nTest, 99991);
  const a = P.alpha / P.m, zc = S.qnorm(1 - a);
  const rows = [];
  for (const k of ['trend', 'two']) {
    const z = test[k], nT = z.length;
    const fpr = thr.map(r => { let c = 0; for (const v of z) if (v > r[k].crit) c++; return c / nT; });
    let naive = 0; for (const v of z) if (v > zc) naive++;
    const mean = fpr.reduce((x, y) => x + y, 0) / fpr.length;
    const sd = Math.sqrt(fpr.reduce((x, y) => x + (y - mean) * (y - mean), 0) / (fpr.length - 1));
    const zero = z.reduce((c, v) => c + (v === 0 ? 1 : 0), 0) / nT;
    rows.push([label, k, P.system, Math.round(S.neFor(P)), P.p0, P.years.length, S.gensOf(P).slice(-1)[0], +a.toPrecision(3),
      thr[0][k].method, mean.toExponential(4), sd.toExponential(3), Math.min(...fpr).toExponential(3), Math.max(...fpr).toExponential(3),
      (naive / nT).toExponential(4), nT, (thr.reduce((x, r) => x + r[k].crit, 0) / thr.length).toFixed(3), zc.toFixed(3), zero.toFixed(4)]);
  }
  return rows;
}

const part = process.argv[2] || 'A';
const header = ['label', 'statistic', 'system', 'Ne', 'p0', 'points', 'generations', 'alpha_per_locus', 'threshold_method',
  'fpr_drift_aware_mean', 'fpr_drift_aware_sd', 'fpr_drift_aware_min', 'fpr_drift_aware_max', 'fpr_naive', 'n_test',
  'threshold_mean', 'threshold_naive', 'share_Z_zero'];
console.log(header.join('\t'));
const t0 = Date.now();
const emit = rows => { for (const r of rows) console.log(r.join('\t')); console.error(rows[0][0], ((Date.now() - t0) / 1000).toFixed(0) + 's'); };

if (part === 'A') { // effective size x inheritance system
  for (const system of ['haploid', 'diploid', 'haplodiploid'])
    for (const Ne of [100, 300, 1000, 3000, 10000, 30000, 100000, 1000000])
      emit(evaluate(`A ${system} Ne=${Ne}`, withNe(BASE, system, Ne), 60000));
}
if (part === 'B') { // one design factor at a time, haplodiploid, Ne = 1,500 (1,000 breeding females and males)
  const ten = [1935, 1945, 1955, 1965, 1975, 1985, 1995, 2005, 2015, 2025];
  const V = [
    ['base: 4 points, 20/20/20/40, 10x/15x, p0 = 0.1', {}],
    ['2 points (1935, 2025)', {years: [1935, 2025], n: [20, 40]}],
    ['10 points, 20 each', {years: ten, n: ten.map(() => 20)}],
    ['p0 = 0.01', {p0: 0.01}],
    ['p0 = 0.5', {p0: 0.5}],
    ['depth 2x at all points', {depthEarly: 2, depthLast: 2, minDP: 1}],
    ['perfect genotypes', {perfect: true}],
    ['siblings: 4 per family', {kSib: 4}],
    ['decline to 10% by the last point', {decline: 0.1}],
    ['bottleneck: Ne 1,500 / 100 / 1,500', {neByBlock: true, neBlocks: [1500, 100, 1500]}],
    ['migration 1% per generation', {mig: 0.01}],
    ['200 individuals per point', {n: [200, 200, 200, 200]}],
  ];
  for (const [lab, o] of V) emit(evaluate('B ' + lab, Object.assign({}, BASE, o), 60000));
}
if (part === 'C') { // small per-locus alpha (Bonferroni): SPOT raises the number of neutral replicates (neutralReps in sim.js)
  const D = [['haplodiploid Ne=300', withNe(BASE, 'haplodiploid', 300)], ['haplodiploid Ne=1500', BASE],
    ['haplodiploid Ne=100000', withNe(BASE, 'haplodiploid', 100000)], ['haplodiploid Ne=1500, p0=0.01', Object.assign({}, BASE, {p0: 0.01})]];
  const pick = process.argv[3] ? process.argv[3].split(',').map(Number) : D.map((_, i) => i);   // optional subset of designs, e.g. 0,1
  for (const [lab, P] of pick.map(i => D[i])) for (const m of [1, 10, 100, 1000])
    emit(evaluate(`C ${lab} m=${m}`, Object.assign({}, P, {m}), m >= 100 ? 1000000 : 120000));
}
