# Figure: power from full simulations vs rescaled simulations (one point per design).
# Run from the repository root after rescale_validation.js: python3 validation/make_figure.py
import csv, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows = list(csv.DictReader(open('validation/rescale_results.tsv'), delimiter='\t'))
fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
cols = {'haploid': '#2B59D9', 'diploid': '#0E9F6E'}
marks = {5000: 'o', 20000: 's', 60000: '^'}
for ax, (ke, kr, title) in zip(axes, [('power_exact', 'power_rescaled', 'Power, trend test'), ('neutralSD_exact', 'neutralSD_rescaled', 'SD of neutral allele frequency at the last point')]):
    lim = 0
    for r in rows:
        x, y = float(r[ke]), float(r[kr]); lim = max(lim, x, y)
        ax.scatter(x, y, s=28, color=cols[r['system']], marker=marks[int(r['G'])], alpha=.8, edgecolor='white', linewidth=.5)
    lim = 1.0 if ke.startswith('power') else lim * 1.05
    ax.plot([0, lim], [0, lim], color='#555', lw=1, ls='--')
    if ke.startswith('power'):
        xs = [i / 100 for i in range(101)]
        se = [2 * math.sqrt(2 * p * (1 - p) / 1000) for p in xs]
        ax.fill_between(xs, [p - e for p, e in zip(xs, se)], [p + e for p, e in zip(xs, se)], color='#999', alpha=.15, lw=0)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect('equal')
    ax.set_xlabel('full simulation'); ax.set_ylabel('rescaled simulation'); ax.set_title(title, fontsize=11)
from matplotlib.lines import Line2D
h = [Line2D([], [], color=c, marker='o', ls='', label=k) for k, c in cols.items()] + [Line2D([], [], color='#444', marker=m, ls='', label=f'{g:,} generations') for g, m in marks.items()]
axes[0].legend(handles=h, fontsize=8, loc='lower right', frameon=False)
fig.suptitle('SPOT rescaling: full vs rescaled simulations (shaded: ±2 Monte Carlo SE)', fontsize=11)
fig.tight_layout(); fig.savefig('validation/rescale_validation.png', dpi=160)
d = [float(r['power_rescaled']) - float(r['power_exact']) for r in rows]
print(len(rows), 'designs; mean power difference %.4f; max |diff| %.3f' % (sum(d) / len(d), max(map(abs, d))))
