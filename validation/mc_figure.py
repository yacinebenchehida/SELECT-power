"""Figure of the Monte Carlo precision analysis (results from monte_carlo_precision.js).
Run from the repository root: python3 validation/mc_figure.py  ->  validation/monte_carlo_precision.svg"""
import csv, json, math, sys
sys.path.insert(0, 'validation')
from svgplot import Fig, Ax, COL, esc

T = [r for r in csv.DictReader(open('validation/mc_precision.tsv'), delimiter='\t') if r['statistic'] == 'trend']
J = json.load(open('validation/mc_precision.json'))
SC = list(J.keys())
pal = [COL['amber'], COL['purple'], COL['blue'], COL['red'], COL['green']]
F = Fig(1060, 440)

a = Ax(F, 75, 50, 330, 300, (400, 6500), (0.4, 5), xlog=True, ylog=True, title='Spread of power between runs', label='a')
a.frame([500, 1000, 2000, 5000], [0.5, 1, 2, 4], 'replicates under selection (2× as many neutral)', 'SD of power (percentage points)')
xs = [400, 6500]
a.line(xs, [100 * math.sqrt(0.25 / x) for x in xs], COL['grey'], 1.4, '5 4')
for k, sc in enumerate(SC):
    rows = [r for r in T if r['scenario'] == sc]
    R = [float(r['replicates']) for r in rows]
    a.line(R, [100 * float(r['sd_power']) for r in rows], pal[k], 1.6)
    a.points(R, [100 * float(r['sd_power']) for r in rows], pal[k], 'o', 4)
a.legend([(sc, pal[k], 'o-') for k, sc in enumerate(SC)] + [('sampling noise only, p = 0.5', COL['grey'], 'dash')], a.X(1350), a.Y(5) + 16)

b = Ax(F, 520, 50, 500, 300, (-0.5, len(SC) - 0.5), (0, 1), title='Power from 40 runs at the default 1,000 replicates', label='b')
b.frame([], [0, .2, .4, .6, .8, 1], '', 'power (trend test, drift-aware)', yfmt=lambda v: f'{v * 100:.0f}%')
for k, sc in enumerate(SC):
    v = J[sc]['1000']['trend']
    m = sum(v) / len(v)
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))
    for i, x in enumerate(v):
        b.points([k - 0.25 + 0.5 * ((i * 0.618) % 1)], [x], pal[k], 'o', 2.6, True)
    b.line([k - 0.32, k + 0.32], [m, m], COL['ink'], 2)
    b.errbar([k + 0.38], [m - 2 * sd], [m + 2 * sd], COL['ink'], 1.6)
    F.text(b.X(k), b.y + b.h + 18, esc({'Long-lived mammal': 'Mammal', 'Human ancient DNA': 'Human aDNA', 'Bacterial pathogen': 'Pathogen', 'Large diploid population': 'Large diploid'}.get(sc, sc)), 11, 'middle', COL['ink'])
    F.text(b.X(k), b.y + b.h + 33, f'{m * 100:.1f}% ± {sd * 100:.1f}', 11, 'middle', COL['muted'])
F.text(75, 418, 'Example scenarios of the page, 40 random seeds each. a: dashed line, binomial standard deviation √(p(1−p)/R) at p = 0.5, the largest possible. '
       'b: black bar = mean, whisker = ±2 SD; numbers = mean ± SD.', 11, color=COL['muted'])
print('figure bytes', F.save('validation/monte_carlo_precision.svg'))
