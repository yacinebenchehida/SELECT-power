"""Agreement of SPOT with theory. Theoretical values are computed here independently of SPOT's code, from the model
equations in the README:
  (a) power with no drift and perfect genotypes: exact enumeration of all possible allele counts at every sampling point;
  (b) neutral drift: expected heterozygosity H_t / H_0 = prod_g (1 - 1/(2 Ne_g))  (haploid: prod (1 - 1/N_g));
  (c) selection without drift: closed form for haploids, p_t = p0 (1+s)^t / (1 - p0 + p0 (1+s)^t), and the
      recursions of the README for diploids and haplodiploids.
Run from the repository root after theory_spot.js:  python3 validation/theory_check.py
Writes theory_power.tsv, theory_drift.tsv, theory_selection.tsv and theory_check.svg in validation/."""
import json, math, sys
import numpy as np
sys.path.insert(0, 'validation')
from svgplot import Fig, Ax, COL

D = json.load(open('validation/theory_spot.json'))
Z95 = 1.6448536269514722


def det(system, h, s, p0, G, smale=1.0):
    """deterministic recursion; element g = genotype/allele frequencies after selection in generation g (selection from g = 1)"""
    pf = pm = p0
    out = []
    for g in range(G + 1):
        sg = s if g > 0 else 0.0
        sm = sg * smale
        if system == 'haploid':
            ps = pf * (1 + sg) / (1 + pf * sg)
            out.append({'p': ps, 'pm': ps, 'pmf': None})
            pf = pm = ps
            continue
        AA, Aa = pf * pm, pf * (1 - pm) + pm * (1 - pf)
        aa = 1 - AA - Aa
        w = AA * (1 + sg) + Aa * (1 + h * sg) + aa
        fAA, fAa, faa = AA * (1 + sg) / w, Aa * (1 + h * sg) / w, aa / w
        pfs = fAA + fAa / 2
        if system == 'haplodiploid':
            pms = pf * (1 + sm) / (pf * (1 + sm) + 1 - pf)
        else:
            wm = AA * (1 + sm) + Aa * (1 + h * sm) + aa
            pms = (AA * (1 + sm) + Aa * (1 + h * sm) / 2) / wm
        out.append({'p': pfs, 'pm': pms, 'pmf': (faa, fAa, fAA)})
        pf, pm = pfs, pms
    return out


def count_pmf(system, state, n_ind):
    """distribution of the number of copies of the allele among n_ind sampled individuals (females if diploid)"""
    if system == 'haploid':
        p = state['p']
        k = np.arange(n_ind + 1)
        return np.array([math.comb(n_ind, int(i)) * p ** i * (1 - p) ** (n_ind - i) for i in k])
    g = np.array(state['pmf'])
    out = np.array([1.0])
    for _ in range(n_ind):
        out = np.convolve(out, g)
    return out


def z_grid(nall, t):
    """trend statistic for every combination of counts (x_1..x_K), alleles nall_i per point, times t_i"""
    K = len(nall)
    grids = np.meshgrid(*[np.arange(n + 1) for n in nall], indexing='ij')
    N = sum(nall)
    tb = sum(n * ti for n, ti in zip(nall, t)) / N
    X = sum(grids)
    num = sum(x * (ti - tb) for x, ti in zip(grids, t))
    den = sum(n * (ti - tb) ** 2 for n, ti in zip(nall, t))
    pb = X / N
    with np.errstate(divide='ignore', invalid='ignore'):
        Z = num / np.sqrt(pb * (1 - pb) * den)
    Z[(pb <= 0) | (pb >= 1)] = 0.0
    return Z


def joint(pmfs):
    P = pmfs[0]
    for q in pmfs[1:]:
        P = np.multiply.outer(P, q)
    return P


def exact_power(r):
    system, h, s, p0, years, n = r['system'], r['h'], r['s'], r['p0'], r['years'], r['n']
    t = [y - years[0] for y in years]
    G = t[-1]
    sel = det(system, h, s, p0, G)
    neu = det(system, h, 0.0, p0, G)
    per = 1 if system == 'haploid' else 2
    out = {}
    for stat, idx in (('trend', list(range(len(t)))), ('two', [0, len(t) - 1])):
        nall = [n[i] * per for i in idx]
        tt = [t[i] for i in idx]
        Z = z_grid(nall, tt)
        Ps = joint([count_pmf(system, sel[t[i]], n[i]) for i in idx])
        Pn = joint([count_pmf(system, neu[t[i]], n[i]) for i in idx])
        zs, pz = Z.ravel(), Pn.ravel()
        o = np.argsort(zs, kind='stable')
        cdf = np.cumsum(pz[o])
        crit = zs[o][np.searchsorted(cdf, 0.95 - 1e-12)]   # smallest z with P(Z <= z) >= 0.95, as SPOT's empirical quantile
        sc = r[stat]['crit']   # SPOT's own threshold, from its neutral replicates
        out[stat] = {'crit': float(crit), 'powerEmp': float(Ps[Z > crit + 1e-9].sum()), 'powerNaive': float(Ps[Z > Z95].sum()),
                     'fprNaive': float(Pn[Z > Z95].sum()), 'fprExact': float(Pn[Z > crit + 1e-9].sum()),
                     'powerAtSpot': float(Ps[Z > sc + 1e-9].sum()), 'sizeAtSpot': float(Pn[Z > sc + 1e-9].sum())}
    return out


# ---------- (a) power ----------
rows, pts = [], []
for r in D['power']:
    ex = exact_power(r)
    for stat in ('trend', 'two'):
        if r['design'] == '2 points' and stat == 'two':
            continue      # identical to the trend test with two points
        for kind, ekind in (('powerNaive', 'powerNaive'), ('powerEmp', 'powerAtSpot')):
            a, b = ex[stat][ekind], r[stat][kind]
            se = math.sqrt(max(a * (1 - a), 1e-12) / r['reps'])
            z = (b - a) / se if a * (1 - a) > 0 else 0.0
            rows.append([r['system'], r['h'], r['p0'], r['s'], r['design'], stat, kind, f'{a:.4f}', f'{b:.4f}', f'{z:.2f}',
                         f'{ex[stat]["crit"]:.3f}', f'{r[stat]["crit"]:.3f}', f'{ex[stat]["fprExact"]:.4f}', f'{ex[stat]["sizeAtSpot"]:.4f}',
                         f'{ex[stat]["powerEmp"]:.4f}'])
            pts.append((a, b, kind, stat, z, ex[stat]['crit'], r[stat]['crit'], r['system'], ex[stat]['sizeAtSpot'], ex[stat]['fprExact'], ex[stat]['powerEmp']))
with open('validation/theory_power.tsv', 'w') as f:
    f.write('\t'.join(['system', 'h', 'p0', 's', 'design', 'statistic', 'threshold', 'power_exact_at_spot_threshold', 'power_spot', 'z',
                       'crit_exact', 'crit_spot', 'size_at_exact_crit', 'size_at_spot_crit', 'power_at_exact_crit']) + '\n')
    for r in rows:
        f.write('\t'.join(map(str, r)).replace('powerNaive', 'naive z(0.95)').replace('powerEmp', 'neutral 95% quantile') + '\n')
zs = np.array([p[4] for p in pts])
d = np.array([p[1] - p[0] for p in pts])
crit_same = sum(abs(p[5] - p[6]) < 1e-6 for p in pts)
print(f'(a) power: {len(pts)} comparisons; mean SPOT - exact {d.mean() * 100:+.2f} points; max |diff| {abs(d).max() * 100:.1f} points;'
      f' |z| > 2: {(abs(zs) > 2).sum()}, > 3: {(abs(zs) > 3).sum()}; RMS z {math.sqrt((zs ** 2).mean()):.2f};'
      f' threshold identical to exact in {crit_same} of {len(pts)}')
for p in pts:
    if abs(p[4]) > 3:
        print('   large z:', p)
for kind in ('powerNaive', 'powerEmp'):
    zz = np.array([p[4] for p in pts if p[2] == kind]); dd = np.array([p[1] - p[0] for p in pts if p[2] == kind])
    print(f'    {kind}: n {len(zz)}, mean diff {dd.mean()*100:+.2f}, max |diff| {abs(dd).max()*100:.1f}, |z|>2 {(abs(zz)>2).sum()}, RMS z {math.sqrt((zz**2).mean()):.2f}')
em = [p for p in pts if p[2] == 'powerEmp']
sz = np.array([p[8] for p in em]); sx = np.array([p[9] for p in em])
print(f'    size of SPOT thresholds (exact): mean {sz.mean():.4f}, range {sz.min():.4f}-{sz.max():.4f}; size at exact 95% quantile: mean {sx.mean():.4f}, range {sx.min():.4f}-{sx.max():.4f}')
dpow = np.array([p[1] - p[10] for p in em])
print(f'    SPOT power vs exact power at the exact quantile: mean {dpow.mean()*100:+.2f}, max |diff| {abs(dpow).max()*100:.1f} points')

# ---------- (b) drift ----------
def ne(system, Nf, Nm):
    if system == 'haploid':
        return Nf
    if system == 'haplodiploid':
        return 9 * Nf * Nm / (4 * Nm + 2 * Nf)
    return 4 * Nf * Nm / (Nf + Nm)


drift_rows, drift_z = [], []
for r in D['drift']:
    o, G = r['params'], r['G']
    sysm, Nf0, Nm0 = o['system'], o['Nf'], o['Nm']
    gens = [0, 20, 40, 90]          # sampling generations of the default design (1935, 1955, 1975, 2025)
    H, theo = 1.0, [1.0]
    for g in range(G):
        if o.get('neByBlock'):
            j = 0
            while j < len(gens) - 2 and g >= gens[j + 1]:
                j += 1
            f = o['neBlocks'][j] / ne(sysm, Nf0, Nm0)
        else:
            f = o.get('decline', 1) ** (g / G)
        Nf, Nm = max(2, round(Nf0 * f)), max(2, round(Nm0 * f))
        e = ne(sysm, Nf, Nm)
        H *= (1 - 1 / e) if sysm == 'haploid' else (1 - 1 / (2 * e))
        theo.append(H)
    sim, se = r['het'], r['se']
    rel = [sim[g] / theo[g] - 1 for g in range(1, G + 1)]
    drift_rows.append((r['label'], theo, sim, r['spotHetRetained']))
    zmax = max(abs(sim[g] - theo[g]) / se[g] for g in range(1, G + 1))
    drift_z.append((r['label'], (sim[-1] - theo[-1]) / se[-1]))
    print(f'(b) {r["label"]:45s} H90/H0 sim {sim[-1]:.4f} ± {se[-1]:.4f} theory {theo[-1]:.4f} (SPOT tile {r["spotHetRetained"]:.4f});'
          f' z at generation 90 {(sim[-1] - theo[-1]) / se[-1]:+.2f}; max relative deviation {max(abs(x) for x in rel) * 100:.2f}%, max |z| over generations {zmax:.2f}')
with open('validation/theory_drift.tsv', 'w') as f:
    f.write('scenario\tgeneration\tH_ratio_simulated\tH_ratio_theory\n')
    for lab, th, sm, _ in drift_rows:
        for g in range(len(th)):
            f.write(f'{lab}\t{g}\t{sm[g]:.5f}\t{th[g]:.5f}\n')

# ---------- (c) selection ----------
maxdiff, sel_rows = 0.0, []
for r in D['selection']:
    G = len(r['f']) - 1
    if r['system'] == 'haploid':
        th = [r['p0'] * (1 + r['s']) ** g / (1 - r['p0'] + r['p0'] * (1 + r['s']) ** g) for g in range(G + 1)]
        thm = th
    else:
        dd = det(r['system'], r['h'], r['s'], r['p0'], G)
        th, thm = [x['p'] for x in dd], [x['pm'] for x in dd]
    md = max(max(abs(a - b) for a, b in zip(r['f'], th)), max(abs(a - b) for a, b in zip(r['m'], thm)))
    maxdiff = max(maxdiff, md)
    sel_rows.append((r, th, md))
print(f'(c) selection: {len(sel_rows)} trajectories (300 generations); largest |SPOT - theory| = {maxdiff:.2e}')
sm = D['stochasticMean']
dd = det('diploid', 0.5, sm['s'], sm['p0'], 300)
print(f'    stochastic mean of {sm["reps"]} replicates (N = 1e6) vs deterministic: max |diff| = '
      f'{max(abs(a - b["p"]) for a, b in zip(sm["f"], dd)):.4f}')
with open('validation/theory_selection.tsv', 'w') as f:
    f.write('system\th\ts\tp0\tmax_abs_difference_spot_vs_theory\n')
    for r, _, md in sel_rows:
        f.write(f'{r["system"]}\t{r["h"]}\t{r["s"]}\t{r["p0"]}\t{md:.2e}\n')

# ---------- figure ----------
F = Fig(1180, 440)
# (a) power scatter
A = Ax(F, 70, 50, 300, 300, (0, 1), (0, 1), title='Power without drift: SPOT vs exact', label='a')
A.frame([0, .2, .4, .6, .8, 1], [0, .2, .4, .6, .8, 1], 'exact power', 'SPOT power')
xs = [i / 100 for i in range(101)]
se = [2 * math.sqrt(x * (1 - x) / 5000) for x in xs]
A.band(xs, [x - e for x, e in zip(xs, se)], [x + e for x, e in zip(xs, se)], COL['grey'], .25)
A.line([0, 1], [0, 1], COL['muted'], 1, '4 4')
for p in pts:
    A.points([p[0]], [p[1]], COL['blue'] if p[2] == 'powerEmp' else COL['red'], 'o' if p[3] == 'trend' else 's', 3.2)
A.legend([('drift-aware threshold (SPOT\'s own)', COL['blue'], 'o'), ('naive threshold z₀.₉₅', COL['red'], 'o'),
          ('trend test, all points', COL['muted'], 'o'), ('first vs last point', COL['muted'], 's')], A.x + 8, A.y + 16)
# (b) heterozygosity
B = Ax(F, 470, 50, 300, 300, (0, 90), (0, 1), title='Neutral drift: heterozygosity retained', label='b')
B.frame([0, 20, 40, 60, 80], [0, .2, .4, .6, .8, 1], 'generation', 'H(t) / H(0)')
pal = [COL['blue'], COL['green'], COL['amber'], COL['purple'], COL['red']]
show = ['haploid, N = 100', 'diploid, Nf = Nm = 50', 'haplodiploid, Nf = Nm = 100', 'diploid, Nf = Nm = 500, decline to 10%',
        'diploid, Ne per interval 1,000 / 50 / 1,000']
leg = []
for k, lab in enumerate(show):
    th, sim = next((t, s) for l, t, s, _ in drift_rows if l == lab)
    B.line(list(range(len(th))), th, pal[k], 1.6)
    g = list(range(0, len(sim), 6))
    B.points(g, [sim[i] for i in g], pal[k], 'o', 3, True)
    leg.append((lab.replace('Ne per interval', 'Ne'), pal[k], 'line'))
B.legend(leg + [('lines: theory, circles: SPOT', COL['muted'], 'oh')], B.x + 8, B.y + B.h - 100)
# (c) haploid selection
C = Ax(F, 870, 50, 280, 300, (0, 300), (0, 1), title='Selection without drift (haploid)', label='c')
C.frame([0, 50, 100, 150, 200, 250, 300], [0, .2, .4, .6, .8, 1], 'generation', 'allele frequency')
for k, s in enumerate([0.01, 0.05, 0.2]):
    for p0 in (0.01, 0.2):
        r, th, _ = next(x for x in sel_rows if x[0]['system'] == 'haploid' and x[0]['s'] == s and x[0]['p0'] == p0)
        C.line(list(range(0, 301, 2)), th[::2], pal[k], 1.6)
        g = list(range(0, 301, 10))
        C.points(g, [r['f'][i] for i in g], pal[k], 'o', 2.6, True)
C.legend([('s = 0.01', pal[0], 'line'), ('s = 0.05', pal[1], 'line'), ('s = 0.2', pal[2], 'line')], C.X(185), C.Y(0.42))
F.text(70, 420, f'a: {len(D["power"])} designs × 2 thresholds; shaded: ±2 Monte Carlo SE (5,000 replicates). b: 100,000 neutral trajectories each, p₀ = 0.3. '
       'c: p₀ = 0.01 and 0.2; lines: closed form, circles: SPOT.', 11, color=COL['muted'])
print('figure bytes', F.save('validation/theory_check.svg'))
