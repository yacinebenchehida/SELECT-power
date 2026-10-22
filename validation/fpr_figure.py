"""Figure of the false-positive rate analysis (results from false_positive_rate.js).
Run from the repository root: python3 validation/fpr_figure.py  ->  validation/false_positive_rate.svg"""
import csv, sys
sys.path.insert(0, 'validation')
from svgplot import Fig, Ax, COL, esc

rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
A_, B_, C_ = rd('validation/fpr_A.tsv'), rd('validation/fpr_B.tsv'), rd('validation/fpr_C.tsv')
G_ = rd('validation/fpr_C_gaussian.tsv')   # same analysis with the earlier Gaussian-fit rule, for comparison
fl = lambda r, k: float(r[k])
F = Fig(1310, 490)
SH = {'haploid': 'o', 'diploid': 's', 'haplodiploid': '^'}

# (a) false-positive rate against Ne, trend test
a = Ax(F, 75, 50, 310, 320, (70, 1.4e6), (0, 0.2), xlog=True, title='False-positive rate against Ne', label='a')
a.frame([100, 1000, 10000, 100000, 1000000], [0, .05, .1, .15, .2], 'effective population size Ne', 'false-positive rate (α = 0.05)')
a.hline(0.05, COL['muted'])
for sysm, sh in SH.items():
    rows = [r for r in A_ if r['system'] == sysm and r['statistic'] == 'trend']
    ne = [fl(r, 'Ne') for r in rows]
    a.line(ne, [fl(r, 'fpr_naive') for r in rows], COL['red'], 1.4)
    a.points(ne, [fl(r, 'fpr_naive') for r in rows], COL['red'], sh, 4)
    m, s = [fl(r, 'fpr_drift_aware_mean') for r in rows], [fl(r, 'fpr_drift_aware_sd') for r in rows]
    a.errbar(ne, [x - y for x, y in zip(m, s)], [x + y for x, y in zip(m, s)], COL['blue'])
    a.line(ne, m, COL['blue'], 1.4)
    a.points(ne, m, COL['blue'], sh, 4)
a.legend([('naive test (ignores drift)', COL['red'], 'o-'), ('SPOT drift-aware test (±1 SD)', COL['blue'], 'o-'),
          ('haploid', COL['muted'], 'o'), ('diploid', COL['muted'], 's'), ('haplodiploid', COL['muted'], '^')], a.X(4000), a.Y(0.175))

# (b) one design factor at a time (haplodiploid, Ne = 1,500), trend test
rows = [r for r in B_ if r['statistic'] == 'trend']
lab = lambda r: r['label'][2:].replace('base: 4 points, 20/20/20/40, 10x/15x, p0 = 0.1', 'reference design').replace('x', '×').replace('p0', 'p₀')
b = Ax(F, 630, 50, 230, 320, (0, 0.25), (len(rows) - 0.5, -0.5), title='Design variations (Ne = 1,500)', label='b')
b.frame([0, .05, .1, .15, .2, .25], [], 'false-positive rate (α = 0.05)', '', grid=False)
b.line([0.05, 0.05], [len(rows) - 0.5, -0.5], COL['muted'], 1.2, '5 4')
for i, r in enumerate(rows):
    y = i
    F.o.append(f'<line x1="{b.x}" x2="{b.x + b.w}" y1="{b.Y(y):.1f}" y2="{b.Y(y):.1f}" stroke="{COL["line"]}"/>')
    F.text(b.x - 8, b.Y(y) + 4, esc(lab(r)), 11, 'end', COL['ink'])
    m, s = fl(r, 'fpr_drift_aware_mean'), fl(r, 'fpr_drift_aware_sd')
    F.o.append(f'<line x1="{b.X(m - s):.1f}" x2="{b.X(m + s):.1f}" y1="{b.Y(y):.1f}" y2="{b.Y(y):.1f}" stroke="{COL["blue"]}" stroke-width="1.4"/>')
    b.points([m], [y], COL['blue'], 'o', 4)
    b.points([fl(r, 'fpr_naive')], [y], COL['red'], 'o', 4)

# (c) small per-locus alpha: realised / nominal false-positive rate, trend test
c = Ax(F, 1000, 50, 280, 320, (3e-5, 0.08), (0.1, 1000), xlog=True, ylog=True, title='Many loci tested (Bonferroni)', label='c')
c.frame([5e-5, 5e-4, 5e-3, 5e-2], [0.1, 1, 10, 100, 1000], 'per-locus α′ = 0.05 / m', 'realised / nominal false-positive rate')
F.o.append(f'<rect x="{c.X(3e-5):.1f}" y="{c.y}" width="{c.X(0.0025) - c.X(3e-5):.1f}" height="{c.h}" fill="{COL["grey"]}" fill-opacity=".12"/>')
F.text(c.X(3.5e-5), c.Y(0.25), 'more neutral replicates', 10, color=COL['muted'])
F.text(c.X(3.5e-5), c.Y(0.25) + 12, '(20,000 and 200,000)', 10, color=COL['muted'])
c.hline(1, COL['muted'])
DS = [('haplodiploid Ne=300', 'o', 'Ne = 300'), ('haplodiploid Ne=1500', 's', 'Ne = 1,500'), ('haplodiploid Ne=100000', '^', 'Ne = 100,000'),
      ('haplodiploid Ne=1500, p0=0.01', 'd', 'Ne = 1,500, p₀ 0.01')]
for key, sh, _ in DS:
    rows = [r for r in C_ if r['label'].startswith('C ' + key + ' m=') and r['statistic'] == 'trend']
    al = [fl(r, 'alpha_per_locus') for r in rows]
    old = [r for r in G_ if r['label'].startswith('C ' + key + ' m=') and r['statistic'] == 'trend']
    ys = [max(0.1, fl(r, 'fpr_drift_aware_mean') / x) for r, x in zip(old, al)]
    c.line(al, ys, COL['grey'], 1.1, '4 3')
    c.points(al, ys, COL['grey'], sh, 3.4, True)
    for col, k in ((COL['red'], 'fpr_naive'), (COL['blue'], 'fpr_drift_aware_mean')):
        ys = [max(0.1, fl(r, k) / x) for r, x in zip(rows, al)]
        c.line(al, ys, col, 1.2)
        c.points(al, ys, col, sh, 3.8)
c.legend([(lab_, COL['muted'], sh) for _, sh, lab_ in DS] + [('earlier Gaussian fit', COL['grey'], 'dash')], c.X(0.0011), c.Y(700))
F.text(75, 440, 'Trend test across all points. Each drift-aware value is the mean realised rate of 20 thresholds, each set by SPOT from its own neutral replicates', 11, color=COL['muted'])
F.text(75, 456, '(2,000 by default; c: up to 200,000) and tested on 60,000 (c: up to 1,000,000) independent neutral replicates. Red: naive test.', 11, color=COL['muted'])
F.text(75, 472, 'Reference design: haplodiploid, 1935/1955/1975/2025, 20/20/20/40 females, 10×/15× depth, p₀ = 0.1. a: Ne varied in each inheritance system. '
       'b: one factor changed at a time.', 11, color=COL['muted'])
print('figure bytes', F.save('validation/false_positive_rate.svg'))
