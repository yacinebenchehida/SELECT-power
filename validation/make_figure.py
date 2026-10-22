# Figure: full against rescaled simulations (one point per design), written as a small SVG.
# Run from the repository root after rescale_validation.js: python3 validation/make_figure.py
# Uses only the Python standard library.
import csv, math

rows = list(csv.DictReader(open('validation/rescale_results.tsv'), delimiter='\t'))
COL = {'haploid': '#2B59D9', 'diploid': '#0E9F6E'}
GENS = [5000, 20000, 60000]
PW, PH, TOP, LEFT, GAP = 330, 330, 70, 75, 130   # plot size and placement
W, H = LEFT + 2 * PW + GAP + 30, TOP + PH + 70
out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Helvetica,Arial,sans-serif" font-size="13">',
       f'<rect width="{W}" height="{H}" fill="#fff"/>',
       f'<text x="{W/2}" y="24" text-anchor="middle" font-size="15">SPOT rescaling: full against rescaled simulations (shaded: ±2 Monte Carlo SE)</text>']

def marker(x, y, shape, col):
    if shape == 0: return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{col}" fill-opacity=".8" stroke="#fff" stroke-width=".6"/>'
    if shape == 1: return f'<rect x="{x-4:.1f}" y="{y-4:.1f}" width="8" height="8" fill="{col}" fill-opacity=".8" stroke="#fff" stroke-width=".6"/>'
    return f'<path d="M{x:.1f} {y-5:.1f}L{x+4.5:.1f} {y+3.5:.1f}H{x-4.5:.1f}Z" fill="{col}" fill-opacity=".8" stroke="#fff" stroke-width=".6"/>'

def panel(x0, ke, kr, title, lim, step, band):
    X = lambda v: x0 + v / lim * PW
    Y = lambda v: TOP + PH - v / lim * PH
    o = [f'<text x="{x0+PW/2}" y="{TOP-10}" text-anchor="middle" font-size="14">{title}</text>']
    if band:  # ±2 SE of a difference between two estimates from 1,000 replicates each
        xs = [i / 100 for i in range(101)]
        up = [(p, min(1, p + 2 * math.sqrt(2 * p * (1 - p) / 1000))) for p in xs]
        lo = [(p, max(0, p - 2 * math.sqrt(2 * p * (1 - p) / 1000))) for p in reversed(xs)]
        o.append('<path d="M' + 'L'.join(f'{X(a):.1f} {Y(b):.1f}' for a, b in up + lo) + 'Z" fill="#999" fill-opacity=".18"/>')
    o.append(f'<line x1="{X(0)}" y1="{Y(0)}" x2="{X(lim)}" y2="{Y(lim)}" stroke="#555" stroke-dasharray="5 4"/>')
    o.append(f'<rect x="{x0}" y="{TOP}" width="{PW}" height="{PH}" fill="none" stroke="#222"/>')
    v = 0.0
    while v <= lim + 1e-9:
        lab = f'{v:.1f}' if step >= 0.1 else f'{v:.2f}'
        o.append(f'<line x1="{X(v):.1f}" y1="{Y(0)}" x2="{X(v):.1f}" y2="{Y(0)+5}" stroke="#222"/><text x="{X(v):.1f}" y="{Y(0)+20}" text-anchor="middle">{lab}</text>')
        o.append(f'<line x1="{X(0)-5}" y1="{Y(v):.1f}" x2="{X(0)}" y2="{Y(v):.1f}" stroke="#222"/><text x="{X(0)-8}" y="{Y(v)+4:.1f}" text-anchor="end">{lab}</text>')
        v += step
    o.append(f'<text x="{x0+PW/2}" y="{TOP+PH+45}" text-anchor="middle">full simulation</text>')
    o.append(f'<text transform="translate({x0-50},{TOP+PH/2}) rotate(-90)" text-anchor="middle">rescaled simulation</text>')
    for r in rows:
        o.append(marker(X(float(r[ke])), Y(float(r[kr])), GENS.index(int(r['G'])), COL[r['system']]))
    return o

lim2 = math.ceil(max(max(float(r['neutralSD_exact']), float(r['neutralSD_rescaled'])) for r in rows) * 1.05 / 0.05) * 0.05
out += panel(LEFT, 'power_exact', 'power_rescaled', 'Power, trend test', 1.0, 0.2, True)
out += panel(LEFT + PW + GAP, 'neutralSD_exact', 'neutralSD_rescaled', 'SD of neutral allele frequency at the last point', lim2, 0.05, False)
# legend inside the first panel
lx, ly = LEFT + PW - 150, TOP + PH - 112
for i, (k, c) in enumerate(COL.items()):
    out += [marker(lx, ly + 18 * i, 0, c), f'<text x="{lx+12}" y="{ly+18*i+4}">{k}</text>']
for i, g in enumerate(GENS):
    out += [marker(lx, ly + 36 + 18 * i, i, '#444'), f'<text x="{lx+12}" y="{ly+36+18*i+4}">{g:,} generations</text>']
out.append('</svg>')
open('validation/rescale_validation.svg', 'w').write('\n'.join(out) + '\n')

d = [float(r['power_rescaled']) - float(r['power_exact']) for r in rows]
print(len(rows), 'designs; mean power difference %.4f; max |diff| %.3f' % (sum(d) / len(d), max(map(abs, d))))
