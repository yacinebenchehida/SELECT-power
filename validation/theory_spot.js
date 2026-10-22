// SPOT side of the comparison with theory (theory_check.py computes the theoretical values independently).
// Writes validation/theory_spot.json. Run from the repository root: node validation/theory_spot.js
const S = require('../sim.js');
const fs = require('fs');
const BASE = {system: 'diploid', neByBlock: false, neBlocks: null, genTime: 1, Nf: 1000, Nm: 1000, noDrift: false,
  p0: 0.1, s: 0, h: 0.5, onset: 1935, years: [1935, 1955, 1975, 2025], n: [20, 20, 20, 40], perfect: true,
  depthEarly: 10, depthLast: 15, err: 0.005, minDP: 3, alpha: 0.05, m: 1, sMale: 1, decline: 1, mig: 0, pMig: -1,
  kSib: 1, maleFrac: 0, reps: 5000, seed: 7};
const out = {power: [], drift: [], selection: []};
const t0 = Date.now();

// (a) power with no drift and perfect genotypes: only sampling noise remains
const SYSTEMS = [['haploid', 0.5], ['diploid', 0.5], ['diploid', 0], ['diploid', 1], ['haplodiploid', 0.5]];
const DESIGNS = [['4 points', [1935, 1955, 1975, 2025], [20, 20, 20, 40]], ['2 points', [1935, 2025], [20, 40]]];
let k = 0;
for (const [system, h] of SYSTEMS) for (const p0 of [0.05, 0.2]) for (const s of [0.005, 0.01, 0.02, 0.05]) for (const [dl, years, n] of DESIGNS) {
  const P = Object.assign({}, BASE, {system, h, p0, s, years, n, noDrift: true, seed: 100 + k++});
  const r = S.runPower(P);
  out.power.push({system, h, p0, s, design: dl, years, n, reps: P.reps,
    trend: {crit: r.trend.crit, powerEmp: r.trend.powerEmp, powerNaive: r.trend.powerNaive, fprNaive: r.trend.fprNaive},
    two: {crit: r.two.crit, powerEmp: r.two.powerEmp, powerNaive: r.two.powerNaive, fprNaive: r.two.fprNaive}});
}
console.error('power', ((Date.now() - t0) / 1000).toFixed(0) + 's');

// (b) neutral drift: mean of p(1 - p) relative to p0(1 - p0), every generation, from 100,000 neutral trajectories
// overall allele frequency: haploid p; diploid (pf + pm) / 2; haplodiploid (2 pf + pm) / 3 (copies carried by each sex)
const DRIFT = [
  ['haploid, N = 100', {system: 'haploid', Nf: 100, Nm: 100}],
  ['haploid, N = 1,000', {system: 'haploid', Nf: 1000, Nm: 1000}],
  ['diploid, Nf = Nm = 50', {system: 'diploid', Nf: 50, Nm: 50}],
  ['diploid, Nf = Nm = 500', {system: 'diploid', Nf: 500, Nm: 500}],
  ['diploid, Nf = 100, Nm = 20', {system: 'diploid', Nf: 100, Nm: 20}],
  ['haplodiploid, Nf = Nm = 100', {system: 'haplodiploid', Nf: 100, Nm: 100}],
  ['haplodiploid, Nf = 300, Nm = 50', {system: 'haplodiploid', Nf: 300, Nm: 50}],
  ['diploid, Nf = Nm = 500, decline to 10%', {system: 'diploid', Nf: 500, Nm: 500, decline: 0.1}],
  ['diploid, Ne per interval 1,000 / 50 / 1,000', {system: 'diploid', Nf: 500, Nm: 500, neByBlock: true, neBlocks: [1000, 50, 1000]}],
];
for (const [lab, o] of DRIFT) {
  const P = S.prep(Object.assign({}, BASE, o, {p0: 0.3, n: [0, 0, 0, 0], perfect: true}));
  const R = S.makeRng(31), reps = 100000, G = S.gensOf(P).slice(-1)[0];
  const het = new Float64Array(G + 1), het2 = new Float64Array(G + 1);
  for (let r = 0; r < reps; r++) {
    const tr = S.simulateRep(P, 0, R, true);
    for (let g = 0; g <= G; g++) {
      const pf = tr.traj[g], pm = tr.trajM[g];
      const p = o.system === 'haploid' ? pf : o.system === 'diploid' ? (pf + pm) / 2 : (2 * pf + pm) / 3;
      het[g] += p * (1 - p); het2[g] += (p * (1 - p)) ** 2;
    }
  }
  out.drift.push({label: lab, params: o, p0: P.p0, G, reps, spotHetRetained: S.hetRetained(Object.assign({}, BASE, o)),
    het: Array.from(het, v => v / reps / (P.p0 * (1 - P.p0))),
    se: Array.from(het, (v, g) => Math.sqrt((het2[g] / reps - (v / reps) ** 2) / reps) / (P.p0 * (1 - P.p0)))});
}
console.error('drift', ((Date.now() - t0) / 1000).toFixed(0) + 's');

// (c) deterministic selection trajectories (no drift), 300 generations
for (const [system, h] of [['haploid', 0.5], ['diploid', 0.5], ['diploid', 0], ['diploid', 1], ['haplodiploid', 0.5], ['haplodiploid', 0]])
  for (const s of [0.01, 0.05, 0.2]) for (const p0 of [0.01, 0.2]) {
    const P = Object.assign({}, BASE, {system, h, s, p0, years: [0, 300], n: [0, 0], onset: 0});
    const d = S.deterministic(P, true);
    out.selection.push({system, h, s, p0, f: Array.from(d.f), m: Array.from(d.m)});
  }
// stochastic mean with a large population should follow the deterministic path (diploid, additive)
{
  const P = S.prep(Object.assign({}, BASE, {system: 'diploid', s: 0.05, p0: 0.01, Nf: 1e6, Nm: 1e6, years: [0, 300], n: [0, 0], onset: 0}));
  const R = S.makeRng(5), reps = 500, mean = new Float64Array(301);
  for (let r = 0; r < reps; r++) { const tr = S.simulateRep(P, P.s, R, true); for (let g = 0; g <= 300; g++) mean[g] += tr.traj[g] / reps; }
  out.stochasticMean = {system: 'diploid', h: 0.5, s: 0.05, p0: 0.01, N: 1e6, reps, f: Array.from(mean)};
}
fs.writeFileSync('validation/theory_spot.json', JSON.stringify(out));
console.error('done', ((Date.now() - t0) / 1000).toFixed(0) + 's');
