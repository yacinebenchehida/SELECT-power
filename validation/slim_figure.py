"""Figure of the SLiM comparison (results from slim_compare.js).
Run from the repository root: python3 validation/slim_figure.py  ->  validation/slim_check.svg"""
import csv, math, sys
sys.path.insert(0, 'validation')
from svgplot import Fig, Ax, COL, esc

rd = lambda f: list(csv.DictReader(open('validation/' + f), delimiter='\t'))
T, P, Q = rd('slim_results.tsv'), rd('slim_power.tsv'), rd('slim_fpr.tsv')
SYS = {'diploid': (COL['blue'], 'o'), 'haploid': (COL['green'], 's'), 'haplodiploid': (COL['amber'], '^')}
sysof = lambda name: name.split('_')[0]
pct = lambda v: f'{v * 100:.0f}%'
F = Fig(1360, 450)

# (a) mean allele frequency, SPOT against SLiM (each point: one scenario and sampling generation)
a = Ax(F, 70, 50, 250, 250, (0, 0.45), (0, 0.45), title='Mean allele frequency', label='a')
a.frame([0, .1, .2, .3, .4], [0, .1, .2, .3, .4], 'SPOT', 'SLiM')
a.line([0, 0.45], [0, 0.45], COL['muted'], 1, '4 4')
for r in T:
    c, sh = SYS[sysof(r['scenario'])]
    a.points([float(r['mean_spot'])], [float(r['mean_slim'])], c, sh, 3.5, r['run'] == 'neutral')
# (b) standard deviation of allele frequency between replicates
b = Ax(F, 400, 50, 250, 250, (0, 0.3), (0, 0.3), title='Spread between replicates (SD)', label='b')
b.frame([0, .1, .2, .3], [0, .1, .2, .3], 'SPOT', 'SLiM')
b.line([0, 0.3], [0, 0.3], COL['muted'], 1, '4 4')
for r in T:
    c, sh = SYS[sysof(r['scenario'])]
    b.points([float(r['sd_spot'])], [float(r['sd_slim'])], c, sh, 3.5, r['run'] == 'neutral')
a.legend([('diploid', SYS['diploid'][0], 'o'), ('haploid', SYS['haploid'][0], 's'), ('haplodiploid', SYS['haplodiploid'][0], '^'),
          ('hollow: neutral runs', COL['muted'], 'oh')], a.x + 8, a.y + 16)
# (c) power, SPOT against SLiM, with +-2 standard errors of the difference
c_ = Ax(F, 730, 50, 250, 250, (0.1, 0.8), (0.1, 0.8), title='Power (drift-aware trend test)', label='c')
c_.frame([.2, .4, .6, .8], [.2, .4, .6, .8], 'SPOT', 'SLiM', xfmt=pct, yfmt=pct)
c_.line([0.1, 0.8], [0.1, 0.8], COL['muted'], 1, '4 4')
for r in P:
    col, sh = SYS[r['system']]
    x, y, se = float(r['power_spot']), float(r['power_slim']), float(r['se_difference'])
    c_.errbar([x], [y - 2 * se], [y + 2 * se], col)
    c_.points([x], [y], col, sh, 4)
# (d) false-positive rate of SPOT's threshold on SLiM's neutral populations, +-2 SE
Qs = Q
d = Ax(F, 1060, 50, 280, 250, (-0.5, len(Qs) - 0.5), (0.03, 0.08), title="SPOT's threshold on SLiM", label='d')
d.frame([], [.03, .04, .05, .06, .07, .08], '', 'false-positive rate (nominal 5%)', yfmt=pct)
d.hline(0.05, COL['muted'])
for i, r in enumerate(Qs):
    col, sh = SYS[r['system']]
    v, se = float(r['fpr_on_slim']), float(r['se'])
    d.errbar([i], [v - 2 * se], [v + 2 * se], col)
    d.points([i], [v], col, sh, 4)
    ne = float(r['NF']) * {'haploid': 1, 'diploid': 2, 'haplodiploid': 1.5}[r['system']]
    for j, t in enumerate([{'haplodiploid': 'haplodipl.'}.get(r['system'], r['system']), f"Ne {ne:,.0f}", f"n = {int(r['n_slim']):,}"]):
        F.text(d.X(i), d.y + d.h + 14 + 12 * j, esc(t), 10, 'middle', COL['ink'] if j < 2 else COL['muted'])
F.text(70, 395, 'Same scenarios simulated by SPOT and by SLiM 5.2 (2,000 replicate populations each; d: up to 8,000 SLiM neutral replicates, SPOT threshold from 20,000). '
       'Sampling at generations 0, 20, 40 and 90,', 11, color=COL['muted'])
F.text(70, 411, '20/20/20/40 females, perfect genotypes. a, b: allele frequency among breeding females after selection at generations 20, 40 and 90. '
       'c, d: bars show ±2 standard errors.', 11, color=COL['muted'])
print('figure bytes', F.save('validation/slim_check.svg'))
