"""Figures of the worked examples (results from worked_examples.js).
Run from the repository root: python3 worked_examples/figures.py  ->  worked_examples/bumblebees.svg, contrasts.svg"""
import csv, sys
sys.path.insert(0, 'validation')
from svgplot import Fig, Ax, COL, esc

rd = lambda f: list(csv.DictReader(open('worked_examples/' + f), delimiter='\t'))
fl = float
pct = lambda v: f'{v * 100:.0f}%'
pal = [COL['blue'], COL['green'], COL['amber'], COL['purple'], COL['red']]

# ---------- Figure 1: bumblebees ----------
F = Fig(1250, 836)
A = rd('bee_s_p0.tsv')
a = Ax(F, 75, 50, 330, 290, (0.004, 0.25), (0, 1), xlog=True, title='Power against selection', label='a')
a.frame([0.005, 0.01, 0.02, 0.05, 0.1, 0.2], [0, .2, .4, .6, .8, 1], 'selection coefficient s', 'power (drift-aware trend test)', yfmt=pct)
a.hline(0.8, COL['muted'])
for k, p0 in enumerate(['0.01', '0.05', '0.1', '0.3']):
    r = [x for x in A if x['p0'] == p0]
    xs, ys = [fl(x['s']) for x in r], [fl(x['power_trend']) for x in r]
    a.line(xs, ys, pal[k]); a.points(xs, ys, pal[k])
a.legend([(f'p₀ = {p}', pal[k], 'o-') for k, p in enumerate(['0.01', '0.05', '0.1', '0.3'])], a.x + 10, a.y + 18)

C = rd('bee_ne.tsv')
b = Ax(F, 520, 50, 290, 290, (110, 2e5), (0, 1), xlog=True, title='Power against Ne (s = 0.05, p₀ = 0.1)', label='b')
b.frame([100, 1000, 10000, 100000], [0, .2, .4, .6, .8, 1], 'effective population size Ne', 'power (drift-aware trend test)', yfmt=pct)
b.hline(0.8, COL['muted'])
c = Ax(F, 920, 50, 290, 290, (110, 2e5), (0, 0.3), xlog=True, title='False positives of the naive test', label='c')
c.frame([100, 1000, 10000, 100000], [0, .05, .1, .15, .2, .25, .3], 'effective population size Ne', 'false-positive rate (nominal 5%)', yfmt=pct)
c.hline(0.05, COL['muted'])
labs = {'1': '20/20/20/40', '2': '40/40/40/80', '5': '100/100/100/200'}
for k, sc in enumerate(['1', '2', '5']):
    r = [x for x in C if x['n_scale'] == sc]
    xs = [fl(x['Ne']) for x in r]
    b.line(xs, [fl(x['power_trend']) for x in r], pal[k]); b.points(xs, [fl(x['power_trend']) for x in r], pal[k])
    c.line(xs, [fl(x['fpr_naive']) for x in r], pal[k]); c.points(xs, [fl(x['fpr_naive']) for x in r], pal[k])
b.legend([(labs[s], pal[k], 'o-') for k, s in enumerate(['1', '2', '5'])], b.X(1500), b.Y(0.42))
c.legend([(labs[s], pal[k], 'o-') for k, s in enumerate(['1', '2', '5'])], c.X(3000), c.Y(0.29))

B = rd('bee_allocation.tsv')
designs = list(dict.fromkeys(x['design'] for x in B))[1:]
d = Ax(F, 330, 450, 420, 290, (-5, 25), (len(designs) - 0.5, -0.5), title='Gain from extra effort over the reference design', label='d')
d.frame([-5, 0, 5, 10, 15, 20, 25], [], 'change in power (percentage points)', '', grid=False)
d.line([0, 0], [len(designs) - 0.5, -0.5], COL['muted'], 1.2, '5 4')
ref = {(x['Ne'], x['p0']): fl(x['power_trend']) for x in B if x['design'].startswith('reference')}
p0s = ['0.01', '0.1', '0.3']
for i, lab in enumerate(designs):
    F.o.append(f'<line x1="{d.x}" x2="{d.x + d.w}" y1="{d.Y(i):.1f}" y2="{d.Y(i):.1f}" stroke="{COL["line"]}"/>')
    F.text(d.x - 8, d.Y(i) + 4, esc(lab), 11, 'end', COL['ink'])
    for j, Ne in enumerate(['1500', '100000']):
        for k, p0 in enumerate(p0s):
            v = next(fl(x['power_trend']) for x in B if x['design'] == lab and x['Ne'] == Ne and x['p0'] == p0)
            d.points([100 * (v - ref[(Ne, p0)])], [i + (-0.17 if j == 0 else 0.17)], pal[k], 'o', 4.5, j == 1)
d.legend([(f'p₀ = {p}', pal[k], 'o') for k, p in enumerate(p0s)] + [('filled: Ne = 1,500 (strong drift)', COL['muted'], 'o'),
          ('hollow: Ne = 100,000 (weak drift)', COL['muted'], 'oh')], d.x + d.w + 20, d.y + 18)
F.text(75, 790, 'Bumblebee scenario of the page: haplodiploid, Nf = Nm = 1,000 (Ne = 1,500), samples 1935/1955/1975/2025, 20/20/20/40 females, depth 10×/15×, α = 0.05. '
       '1,000 replicates under selection per point; Monte Carlo SD about 1.5 points.', 11, color=COL['muted'])
F.text(75, 806, 'b, c: breeding numbers varied (Ne shown) with 1×, 2× and 5× the reference sample sizes. d: each option adds 60 individuals, '
       'doubles historical depth, or moves the same 100 individuals to the first and last dates;', 11, color=COL['muted'])
F.text(75, 822, 's chosen so that the reference design has 55–67% power; 5,000 replicates under selection per point (SD of a difference about 1 point).', 11, color=COL['muted'])
print('bumblebees.svg', F.save('worked_examples/bumblebees.svg'))

# ---------- Figure 2: contrasting examples ----------
F = Fig(1000, 440)
M = rd('mammal.tsv')
e = Ax(F, 75, 50, 340, 290, (0, 0.52), (0, 1), title='Long-lived mammal (7 generations)', label='a')
e.frame([0, .1, .2, .3, .4, .5], [0, .2, .4, .6, .8, 1], 'selection coefficient s (per generation)', 'power', yfmt=pct)
e.hline(0.8, COL['muted'])
hist = list(dict.fromkeys(x['history'] for x in M))
order = [hist[1], hist[0], hist[2]]
for k, h in enumerate(order):
    r = [x for x in M if x['history'] == h]
    xs = [fl(x['s']) for x in r]
    e.line(xs, [fl(x['power_trend']) for x in r], pal[k]); e.points(xs, [fl(x['power_trend']) for x in r], pal[k])
    e.line(xs, [fl(x['power_naive']) for x in r], pal[k], 1.2, '4 3')
e.legend([(h.replace('bottleneck: ', 'bottleneck, ').replace('constant', 'constant'), pal[k], 'o-') for k, h in enumerate(order)]
         + [('dashed: naive test', COL['muted'], 'dash')], e.x + 10, e.y + 18)
P = rd('pathogen.tsv')
f = Ax(F, 570, 50, 340, 290, (1.5e-5, 7e-4), (0, 1), xlog=True, title='Bacterial pathogen (17,000 generations, rescaled)', label='b')
f.frame([2e-5, 5e-5, 1e-4, 2e-4, 5e-4], [0, .2, .4, .6, .8, 1], 'selection coefficient s (per generation)', 'power', yfmt=pct)
f.hline(0.8, COL['muted'])
xs = [fl(x['s_per_generation']) for x in P]
f.line(xs, [fl(x['power_trend']) for x in P], pal[0]); f.points(xs, [fl(x['power_trend']) for x in P], pal[0])
f.line(xs, [fl(x['power_naive']) for x in P], pal[0], 1.2, '4 3')
f.legend([('drift-aware test', pal[0], 'o-'), ('naive test', pal[0], 'dash')], f.x + 10, f.y + 18)
F.text(75, 392, 'a: diploid, generation time 20 years, samples 1880/1950/2020 (20/20/40), Ne 5,000 before 1950 and as shown after; p₀ = 0.1.', 11, color=COL['muted'])
F.text(75, 408, 'b: haploid, 100 generations per year, N = 100,000, samples 1850/1900/1950/2000/2020; p₀ = 0.05.', 11, color=COL['muted'])
F.text(75, 424, 'Example scenarios of the page; 1,000 (a) and 500 (b) replicates under selection per point.', 11, color=COL['muted'])
print('contrasts.svg', F.save('worked_examples/contrasts.svg'))
