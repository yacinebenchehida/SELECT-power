// Worked examples for the paper: what SPOT tells a user planning a study.
// Run from the repository root:  node worked_examples/worked_examples.js A > worked_examples/bee_s_p0.tsv   (parts A-E)
//   A: bumblebees, power against s for several starting frequencies
//   B: bumblebees, where to add samples (historical, modern, depth)
//   C: bumblebees, power against Ne for 1x, 2x and 5x the sample sizes
//   D: long-lived mammal, power with and without a bottleneck between sampling dates
//   E: bacterial pathogen (haploid, many generations, rescaled), power against s
const S = require('../sim.js');
const BEE = {system: 'haplodiploid', neByBlock: false, neBlocks: null, genTime: 1, Nf: 1000, Nm: 1000, noDrift: false,
  p0: 0.01, s: 0.05, h: 0.5, onset: 1935, years: [1935, 1955, 1975, 2025], n: [20, 20, 20, 40], perfect: false,
  depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3, alpha: 0.05, m: 1, sMale: 1, decline: 1, mig: 0, pMig: -1,
  kSib: 1, maleFrac: 0, reps: 1000, seed: 11};
// same values as the example scenarios of the page (PRESETS in parts/ui.js)
const MAMMAL = Object.assign({}, BEE, {system: 'diploid', genTime: 20, Nf: 2500, Nm: 2500, neByBlock: true, neBlocks: [5000, 200],
  p0: 0.1, s: 0.2, onset: 1880, years: [1880, 1950, 2020], n: [20, 20, 40], depthEarly: 4, depthLast: 15});
const PATHOGEN = Object.assign({}, BEE, {system: 'haploid', genTime: 0.01, Nf: 100000, Nm: 100000, p0: 0.05, s: 0.0002,
  onset: 1850, years: [1850, 1900, 1950, 2000, 2020], n: [10, 10, 15, 20, 40], depthEarly: 10, depthLast: 30});

const run = (P, extra) => { const r = S.runPower(Object.assign({}, P, extra));
  return [r.trend.powerEmp.toFixed(4), r.trend.powerNaive.toFixed(4), r.trend.fprNaive.toFixed(4), r.two.powerEmp.toFixed(4)]; };
const H = ['power_trend', 'power_naive', 'fpr_naive', 'power_first_vs_last'];
const part = process.argv[2];
const t0 = Date.now(), log = s => console.error(s, ((Date.now() - t0) / 1000).toFixed(0) + 's');

if (part === 'A') {
  console.log(['p0', 's', ...H].join('\t'));
  for (const p0 of [0.01, 0.05, 0.1, 0.3]) { for (const s of [0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2])
    console.log([p0, s, ...run(BEE, {p0, s})].join('\t')); log('p0 ' + p0); }
}
if (part === 'B') { // each option adds 60 individuals or doubles historical depth; 5,000 replicates to keep Monte Carlo error below 1 point
  // two regimes: Ne = 1,500 (drift strong over 90 generations) and Ne = 100,000 (drift weak, sampling error dominates);
  // in each, s is chosen so that the reference design has about 55-65% power
  const D = [['reference: 20/20/20/40', {}], ['+60 historical: 40/40/40/40', {n: [40, 40, 40, 40]}],
    ['+60 modern: 20/20/20/100', {n: [20, 20, 20, 100]}], ['historical depth 20x (instead of 10x)', {depthEarly: 20}],
    ['first and last only: 50/50 (same total)', {years: [1935, 2025], n: [50, 50]}]];
  const R = [[1500, 1000, [[0.01, 0.05], [0.1, 0.02], [0.3, 0.015]]], [100000, 66667, [[0.01, 0.035], [0.1, 0.015], [0.3, 0.012]]]];
  console.log(['Ne', 'design', 'p0', 's', 'individuals', ...H].join('\t'));
  for (const [Ne, N, cases] of R) for (const [p0, s] of cases) { for (const [lab, o] of D) {
    const P = Object.assign({}, BEE, {Nf: N, Nm: N, reps: 5000}, o);
    console.log([Ne, lab, p0, s, P.n.reduce((a, b) => a + b, 0), ...run(P, {p0, s})].join('\t')); } log('B ' + Ne + ' ' + p0); }
}
if (part === 'C') {
  console.log(['N_breeding_each_sex', 'Ne', 'het_retained', 'n_scale', 'n', ...H].join('\t'));
  for (const k of [1, 2, 5]) { for (const N of [100, 200, 300, 500, 1000, 2000, 5000, 20000, 100000]) {
    const n = BEE.n.map(x => x * k), P = Object.assign({}, BEE, {Nf: N, Nm: N, n, p0: 0.1, s: 0.05});
    console.log([N, Math.round(S.neFor(P)), S.hetRetained(P).toFixed(4), k, n.join('/'), ...run(P)].join('\t')); } log('C x' + k); }
}
if (part === 'D') {
  console.log(['history', 's', ...H].join('\t'));
  for (const [lab, o] of [['bottleneck: Ne 5,000 then 200', {}], ['constant Ne 5,000', {neBlocks: [5000, 5000]}], ['severe bottleneck: Ne 5,000 then 50', {neBlocks: [5000, 50]}]]) {
    for (const s of [0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5]) console.log([lab, s, ...run(MAMMAL, Object.assign({s}, o))].join('\t')); log(lab); }
}
if (part === 'E') {
  console.log(['s_per_generation', 's_per_year', ...H].join('\t'));
  for (const s of [0.00002, 0.00005, 0.0001, 0.0002, 0.0005]) { console.log([s, (s * 100).toPrecision(2), ...run(PATHOGEN, {s, reps: 500})].join('\t')); log('E ' + s); }
}
