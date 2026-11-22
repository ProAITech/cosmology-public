"""Draft manuscript, Figure 1 (SVG, no plotting library).  Reads explore_union.json.
Left: scatter of the foreground-subtracted RM against redshift in the union of 10,449 sight lines (filled) and in the three
independent samples (open).  Right: the net slope of each independent sample, their weighted mean, and the union.
Run: .venv/bin/python -I entries/174-quaia-replication/figure.py   -> figure1.svg
"""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
r = json.load(open(os.path.join(HERE, 'explore_union.json')))
W, H = 920, 400
INK, GRID = '#222', '#bbb'
COL = dict(union='#1f4e9c', entry171='#b5651d', desi_only='#6a8f2f', quaia_new='#8a3fa0')
NAME = dict(union='all 10,449 sight lines', entry171='Milliquas (discovery), 4,867', desi_only='DESI only, 1,833', quaia_new='Quaia (new), 3,749')
mn = lambda v: ('%s' % v).replace('-', '−')
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="Helvetica, Arial, sans-serif" font-size="12">' % (W, H),
     '<rect width="%d" height="%d" fill="#fff"/>' % (W, H)]
# ---- left panel
L, T, PW, PH = 62, 28, 440, 310
xr, yr = (np.log10(1.3), np.log10(4.0)), (6.5, 13.5)
sx = lambda z: L + (np.log10(1 + z) - xr[0]) / (xr[1] - xr[0]) * PW
sy = lambda v: T + PH - (v - yr[0]) / (yr[1] - yr[0]) * PH
o.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="%s"/>' % (L, T, PW, PH, INK))
for z in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
    o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="%s"/><text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (sx(z), T + PH, sx(z), T + PH + 5, INK, sx(z), T + PH + 19, z))
for v in range(7, 14):
    o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s"/><text x="%d" y="%.1f" text-anchor="end">%d</text>' % (L - 5, sy(v), L, sy(v), INK, L - 8, sy(v) + 4, v))
    o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" stroke-width="0.4"/>' % (L, sy(v), L + PW, sy(v), GRID))
o.append('<text x="%.1f" y="%d" text-anchor="middle">redshift z (axis linear in log(1+z))</text>' % (L + PW / 2, T + PH + 38))
o.append('<text transform="translate(18,%.1f) rotate(-90)" text-anchor="middle">scatter of ΔRM, 1.4826 × MAD (rad m⁻²)</text>' % (T + PH / 2))
def series(b, col, filled, dx=0.0, shape='c'):
    for z, s, e in zip(b['zmed'], b['sigma'], b['err']):
        x, y = sx(z) + dx, sy(s)
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%s"/>' % (x, sy(s - e), x, sy(s + e), col, 1.6 if filled else 0.8))
        if shape == 'c':
            o.append('<circle cx="%.1f" cy="%.1f" r="%s" fill="%s" stroke="%s"/>' % (x, y, 5 if filled else 3.2, col if filled else '#fff', col))
        elif shape == 's':
            o.append('<rect x="%.1f" y="%.1f" width="6.4" height="6.4" fill="#fff" stroke="%s"/>' % (x - 3.2, y - 3.2, col))
        else:
            o.append('<path d="M%.1f %.1f l3.8 6.6 l-7.6 0 z" fill="#fff" stroke="%s"/>' % (x, y - 4.4, col))
series(r['parts']['desi_only']['bins'], COL['desi_only'], False, shape='t')
series(r['parts']['quaia_new']['bins'], COL['quaia_new'], False, shape='s')
series(r['parts']['entry171']['bins'], COL['entry171'], False, shape='c')
u = r['union_12_cells']
series(u['bins'], COL['union'], True)
# weighted least-squares line through the union's bins, drawn against log10(1+z_median)
xs = np.log10(1 + np.array(u['bins']['zmed'])); ys = np.array(u['bins']['sigma']); w = 1 / np.array(u['bins']['err']) ** 2
bfit = ((w.sum() * (w * xs * ys).sum() - (w * xs).sum() * (w * ys).sum()) / (w.sum() * (w * xs * xs).sum() - (w * xs).sum() ** 2))
afit = ((w * ys).sum() - bfit * (w * xs).sum()) / w.sum()
z0, z1 = u['bins']['zmed'][0], u['bins']['zmed'][-1]
o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.4" stroke-dasharray="5 3"/>' % (
    sx(z0), sy(afit + bfit * np.log10(1 + z0)), sx(z1), sy(afit + bfit * np.log10(1 + z1)), COL['union']))
ly = T + 14
for key, shape in (('union', 'C'), ('entry171', 'c'), ('quaia_new', 's'), ('desi_only', 't')):
    x = L + 14
    if shape == 'C':
        o.append('<circle cx="%d" cy="%d" r="5" fill="%s"/>' % (x, ly - 4, COL[key]))
    elif shape == 'c':
        o.append('<circle cx="%d" cy="%d" r="3.2" fill="#fff" stroke="%s"/>' % (x, ly - 4, COL[key]))
    elif shape == 's':
        o.append('<rect x="%.1f" y="%.1f" width="6.4" height="6.4" fill="#fff" stroke="%s"/>' % (x - 3.2, ly - 7.2, COL[key]))
    else:
        o.append('<path d="M%d %.1f l3.8 6.6 l-7.6 0 z" fill="#fff" stroke="%s"/>' % (x, ly - 8.4, COL[key]))
    o.append('<text x="%d" y="%d">%s</text>' % (x + 12, ly, NAME[key]))
    ly += 17
# ---- right panel: net slopes
L2, PW2 = 640, 250
xr2 = (-6.0, 10.0)
sx2 = lambda v: L2 + (v - xr2[0]) / (xr2[1] - xr2[0]) * PW2
o.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="%s"/>' % (L2, T, PW2, PH, INK))
for v in (-4, 0, 4, 8):
    o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="%s"/><text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (sx2(v), T + PH, sx2(v), T + PH + 5, INK, sx2(v), T + PH + 19, mn(v)))
o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="%s" stroke-dasharray="3 3"/>' % (sx2(0), T, sx2(0), T + PH, INK))
o.append('<text x="%.1f" y="%d" text-anchor="middle">net slope (rad m⁻² per unit log₁₀(1+z))</text>' % (L2 + PW2 / 2, T + PH + 38))
c = r['combined_independent']
rowsR = [('Milliquas (discovery)', c['net_slopes'][0], c['sds'][0], COL['entry171'], False),
         ('DESI only', c['net_slopes'][1], c['sds'][1], COL['desi_only'], False),
         ('Quaia (new), corrected', c['net_slopes'][2], c['sds'][2], COL['quaia_new'], False),
         ('two follow-ups, mean', c['replications_only']['weighted_mean'], c['replications_only']['sd'], INK, False),
         ('three samples, mean', c['weighted_mean'], c['sd'], INK, True),
         ('union, one fit', u['net'], u['null_sd'], COL['union'], True)]
for i, (name, v, e, col, filled) in enumerate(rowsR):
    y = T + 34 + i * 47
    o.append('<text x="%d" y="%.1f" text-anchor="end">%s</text>' % (L2 - 8, y + 4, name))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.6"/>' % (sx2(v - e), y, sx2(v + e), y, col))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="0.7"/>' % (sx2(v - 2 * e), y, sx2(v + 2 * e), y, col))
    o.append('<circle cx="%.1f" cy="%.1f" r="%s" fill="%s" stroke="%s"/>' % (sx2(v), y, 5 if filled else 3.6, col if filled else '#fff', col))
    o.append('<text x="%.1f" y="%.1f" font-size="11" fill="#444">%s ± %.2f</text>' % (sx2(xr2[1]) - 78, y - 8, mn('%+.2f' % v), e))
o.append('</svg>')
open(os.path.join(HERE, 'figure1.svg'), 'w', encoding='utf-8').write('\n'.join(o))
print('written figure1.svg; union line slope %.2f' % bfit)
