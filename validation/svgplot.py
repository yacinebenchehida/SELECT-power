"""Minimal SVG plotting (Python standard library only), used for the validation figures.
Figures stay small enough to view directly on GitHub and are fully reproducible without plotting packages."""
import math

COL = {'blue': '#2B59D9', 'green': '#0E9F6E', 'red': '#D9363E', 'amber': '#C27A12', 'purple': '#8B4FD6',
       'grey': '#9AA3B8', 'ink': '#0F1730', 'muted': '#5A6480', 'line': '#E1E5EF'}


def esc(t):
    return str(t).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def fmt(v):
    if v == 0:
        return '0'
    a = abs(v)
    if a >= 10000:
        e = int(math.floor(math.log10(a) + 1e-9))
        m = v / 10 ** e
        return ('10' if abs(m - 1) < 1e-9 else f'{m:g}×10') + f'<tspan dy="-5" font-size="8">{e}</tspan>'
    if a >= 1:
        return f'{v:,.0f}' if abs(v - round(v)) < 1e-9 else f'{v:g}'
    return f'{v:g}'


class Fig:
    def __init__(self, W, H, title=None):
        self.W, self.H, self.o = W, H, [f'<rect width="{W}" height="{H}" fill="#fff"/>']
        self.defs, self.defo = {}, []
        if title:
            self.o.append(f'<text x="{W / 2}" y="22" text-anchor="middle" font-size="15" font-weight="600" fill="{COL["ink"]}">{esc(title)}</text>')

    def text(self, x, y, t, size=12, anchor='start', color=None, weight=None, rot=None):
        tr = f' transform="rotate({rot} {x:.1f} {y:.1f})"' if rot else ''
        w = f' font-weight="{weight}"' if weight else ''
        self.o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" fill="{color or COL["ink"]}"{w}{tr}>{t}</text>')

    def save(self, path):
        s = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.W}" height="{self.H}" viewBox="0 0 {self.W} {self.H}" '
             f'font-family="Helvetica,Arial,sans-serif">\n<defs>' + ''.join(self.defo) + '</defs>\n' + '\n'.join(self.o) + '\n</svg>\n')
        open(path, 'w').write(s)
        return len(s)


class Ax:
    def __init__(self, fig, x, y, w, h, xlim, ylim, xlog=False, ylog=False, title=None, label=None):
        self.f, self.x, self.y, self.w, self.h = fig, x, y, w, h
        self.xlim, self.ylim, self.xlog, self.ylog = xlim, ylim, xlog, ylog
        if title:
            fig.text(x, y - 10, esc(title), 13, weight='600')
        if label:
            fig.text(x - 52, y - 10, label, 15, weight='700')

    def _t(self, v, lim, log):
        a, b = lim
        if log:
            v, a, b = math.log10(v), math.log10(a), math.log10(b)
        return (v - a) / (b - a)

    def X(self, v):
        return self.x + self._t(v, self.xlim, self.xlog) * self.w

    def Y(self, v):
        return self.y + self.h - self._t(v, self.ylim, self.ylog) * self.h

    def frame(self, xt, yt, xlab='', ylab='', xfmt=fmt, yfmt=fmt, grid=True):
        o, c = self.f.o, COL
        for v in yt:
            yy = self.Y(v)
            if grid:
                o.append(f'<line x1="{self.x}" x2="{self.x + self.w}" y1="{yy:.1f}" y2="{yy:.1f}" stroke="{c["line"]}"/>')
            o.append(f'<text x="{self.x - 6}" y="{yy + 4:.1f}" font-size="11" text-anchor="end" fill="{c["muted"]}">{yfmt(v)}</text>')
        for v in xt:
            xx = self.X(v)
            o.append(f'<line x1="{xx:.1f}" x2="{xx:.1f}" y1="{self.y + self.h}" y2="{self.y + self.h + 4}" stroke="{c["muted"]}"/>')
            o.append(f'<text x="{xx:.1f}" y="{self.y + self.h + 17}" font-size="11" text-anchor="middle" fill="{c["muted"]}">{xfmt(v)}</text>')
        o.append(f'<rect x="{self.x}" y="{self.y}" width="{self.w}" height="{self.h}" fill="none" stroke="{c["ink"]}" stroke-width=".8"/>')
        if xlab:
            self.f.text(self.x + self.w / 2, self.y + self.h + 36, xlab, 12, 'middle', c['ink'])
        if ylab:
            self.f.text(self.x - 44, self.y + self.h / 2, ylab, 12, 'middle', c['ink'], rot=-90)

    def line(self, xs, ys, color, width=1.8, dash=None, opacity=1):
        d = 'M' + 'L'.join(f'{self.X(a):.1f} {self.Y(b):.1f}' for a, b in zip(xs, ys))
        da = f' stroke-dasharray="{dash}"' if dash else ''
        self.f.o.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}"{da} stroke-opacity="{opacity}"/>')

    def hline(self, v, color, dash='5 4', width=1.2):
        self.line([self.xlim[0], self.xlim[1]], [v, v], color, width, dash)

    def band(self, xs, lo, hi, color, opacity=.15):
        pts = [(self.X(a), self.Y(b)) for a, b in zip(xs, hi)] + [(self.X(a), self.Y(b)) for a, b in zip(reversed(xs), reversed(lo))]
        self.f.o.append('<path d="M' + 'L'.join(f'{a:.1f} {b:.1f}' for a, b in pts) + f'Z" fill="{color}" fill-opacity="{opacity}"/>')

    def marker(self, px, py, color, shape='o', size=4, hollow=False):
        """markers are defined once in <defs> and placed with <use>, which keeps files small"""
        key = (color, shape, size, hollow)
        f = self.f
        if key not in f.defs:
            mid = f'm{len(f.defs)}'
            fill, st = ('#fff', color) if hollow else (color, '#fff')
            sw = 1.5 if hollow else .6
            s = size
            if shape == 'o':
                g = f'<circle r="{s}"'
            elif shape == 's':
                g = f'<rect x="{-s}" y="{-s}" width="{2 * s}" height="{2 * s}"'
            elif shape == '^':
                g = f'<path d="M0 {-s - 1}L{s + .5} {s - .5}H{-s - .5}Z"'
            else:
                g = f'<path d="M0 {-s - 1}L{s + 1} 0L0 {s + 1}L{-s - 1} 0Z"'
            f.defs[key] = mid
            f.defo.append(f'<g id="{mid}">{g} fill="{fill}" stroke="{st}" stroke-width="{sw}"/></g>')
        return f'<use href="#{f.defs[key]}" x="{px:.1f}" y="{py:.1f}"/>'

    def points(self, xs, ys, color, shape='o', size=4, hollow=False):
        for a, b in zip(xs, ys):
            self.f.o.append(self.marker(self.X(a), self.Y(b), color, shape, size, hollow))

    def errbar(self, xs, lo, hi, color, width=1.2):
        for a, l, h in zip(xs, lo, hi):
            xx = self.X(a)
            self.f.o.append(f'<line x1="{xx:.1f}" x2="{xx:.1f}" y1="{self.Y(l):.1f}" y2="{self.Y(h):.1f}" stroke="{color}" stroke-width="{width}"/>')

    def legend(self, items, x=None, y=None, dy=17):
        """items: (label, color, kind) with kind 'line', 'dash', or a marker shape ('o', 's', '^', 'd'); 'o-' = line + marker"""
        x = self.x + 10 if x is None else x
        y = self.y + 14 if y is None else y
        for i, (lab, color, kind) in enumerate(items):
            yy = y + i * dy
            if kind in ('line', 'dash', 'o-', 's-', '^-'):
                da = ' stroke-dasharray="5 3"' if kind == 'dash' else ''
                self.f.o.append(f'<line x1="{x}" x2="{x + 20}" y1="{yy - 4}" y2="{yy - 4}" stroke="{color}" stroke-width="2"{da}/>')
                if kind.endswith('-'):
                    self.f.o.append(self.marker(x + 10, yy - 4, color, kind[0], 3.5))
            else:
                self.f.o.append(self.marker(x + 10, yy - 4, color, kind.rstrip('h'), 4, kind.endswith('h')))
            self.f.text(x + 27, yy, esc(lab), 11, color=COL['ink'])
