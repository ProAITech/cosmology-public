"""Entry 176, Figure 1 (SVG, no plotting library).  Reads results.json and explore.json.
Left: the noise-removed variance of dRM against distance from the nearest eRASS1 cluster in units of R500, with the rotated-cluster
null (mean and one standard deviation).  Right: the excess variance inside R500 for four quarters of the clusters by mass.
Run: .venv/bin/python -I entries/176-cluster-faraday-profile/figure.py   -> figure1.svg
"""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
r = json.load(open(os.path.join(HERE, 'results.json'))); e = json.load(open(os.path.join(HERE, 'explore.json')))
oh = json.load(open(os.path.join(HERE, 'other_half_profile.json')))      # Entry 179
W, H = 920, 400
INK, GRID, BLUE, RED, ORANGE = '#222', '#bbb', '#1f4e9c', '#b03a2e', '#d9822b'
mn = lambda v: ('%s' % v).replace('-', '−')
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="Helvetica, Arial, sans-serif" font-size="12">' % (W, H),
     '<rect width="%d" height="%d" fill="#fff"/>' % (W, H)]
def panel(L, T, PW, PH, xr, yr, xt, yt, xl, yl, logx=False):
    fx = (lambda v: np.log10(v)) if logx else (lambda v: v)
    sx = lambda v: L + (fx(v) - fx(xr[0])) / (fx(xr[1]) - fx(xr[0])) * PW
    sy = lambda v: T + PH - (v - yr[0]) / (yr[1] - yr[0]) * PH
    o.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="%s"/>' % (L, T, PW, PH, INK))
    for v in xt:
        o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="%s"/><text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (sx(v), T + PH, sx(v), T + PH + 5, INK, sx(v), T + PH + 19, mn(v)))
    for v in yt:
        o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s"/><text x="%d" y="%.1f" text-anchor="end">%s</text>' % (L - 5, sy(v), L, sy(v), INK, L - 8, sy(v) + 4, mn(v)))
        o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" stroke-width="0.4"/>' % (L, sy(v), L + PW, sy(v), GRID))
    o.append('<text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (L + PW / 2, T + PH + 38, xl))
    o.append('<text transform="translate(%d,%.1f) rotate(-90)" text-anchor="middle">%s</text>' % (L - 42, T + PH / 2, yl))
    return sx, sy
P = r['reported_only']['profile']
sx, sy = panel(64, 28, 430, 310, (0, 8), (0, 340), (0, 1, 2, 3, 4, 6, 8), (0, 50, 100, 150, 200, 250, 300),
               'distance from the nearest cluster, in units of R₅₀₀', 'variance of ΔRM, noise removed (rad² m⁻⁴)')
for p in P:
    o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#ddd"/>' % (sx(p['lo']), sy(p['null_mean'] + p['null_sd']), sx(p['hi']) - sx(p['lo']), sy(p['null_mean'] - p['null_sd']) - sy(p['null_mean'] + p['null_sd'])))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#777"/>' % (sx(p['lo']), sy(p['null_mean']), sx(p['hi']), sy(p['null_mean'])))
for p in P:
    xm = 0.5 * (p['lo'] + p['hi'])
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.2"/>' % (sx(p['lo']), sy(p['V']), sx(p['hi']), sy(p['V']), BLUE))
    o.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s"/>' % (sx(xm), sy(p['V']), BLUE))
    o.append('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="10.5" fill="#444">%s</text>' % (sx(xm), sy(p['V']) - 9, mn('%+.1f' % p['t'])))
# Entry 179: clusters on the other half of the sky (MCXC and PSZ2); its excess is drawn on the eRASS1 null of the same range
base = {(p['lo'], p['hi']): p['null_mean'] for p in P}
for q in oh['profile']:
    xm = 0.5 * (q['lo'] + q['hi']) + 0.12
    y = q['excess_variance'] + base.get((q['lo'], q['hi']), q['null_mean'])
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.2"/>' % (sx(xm), sy(max(y - q['null_sd'], 0)), sx(xm), sy(min(y + q['null_sd'], 340)), ORANGE))
    o.append('<rect x="%.1f" y="%.1f" width="7" height="7" fill="%s" transform="rotate(45 %.1f %.1f)"/>' % (sx(xm) - 3.5, sy(min(y, 340)) - 3.5, ORANGE, sx(xm), sy(min(y, 340))))
o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-dasharray="4 3"/>' % (sx(1), sy(0), sx(1), sy(300), RED))
o.append('<text x="%.1f" y="%.1f" fill="%s" font-size="11">R₅₀₀</text>' % (sx(1) + 5, sy(300) + 14, RED))
lx, ly = sx(3.2), sy(270)
o.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s"/><text x="%.1f" y="%.1f">9,050 eRASS1 clusters (numbers: t)</text>' % (lx, ly - 4, BLUE, lx + 10, ly))
o.append('<rect x="%.1f" y="%.1f" width="14" height="9" fill="#ddd"/><text x="%.1f" y="%.1f">the same clusters rotated in longitude, ±1σ</text>' % (lx - 7, ly + 9, lx + 10, ly + 17))
o.append('<rect x="%.1f" y="%.1f" width="7" height="7" fill="%s" transform="rotate(45 %.1f %.1f)"/><text x="%.1f" y="%.1f">1,101 MCXC + PSZ2 clusters, other half of the sky</text>' % (lx - 3.5, ly + 27.5, ORANGE, lx, ly + 31, lx + 10, ly + 35))
o.append('<text x="%.1f" y="%.1f" font-size="10.5" fill="#555">(excess ± 1σ, drawn on the eRASS1 null)</text>' % (lx + 10, ly + 49))
M = e['mass']
sx2, sy2 = panel(600, 28, 290, 310, (3, 100), (-50, 300), (5, 10, 20, 50, 100), (0, 100, 200, 300),
                 'median M₅₀₀ of the quarter (10¹³ M☉)', 'excess variance inside R₅₀₀ (rad² m⁻⁴)', logx=True)
o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (sx2(3), sy2(0), sx2(100), sy2(0), INK))
a = e['mass_scaling']['exponent']
ms = np.array([x['m500_median_1e13'] for x in M]); ex = np.array([x['ranges'][2]['excess_variance'] for x in M]); sd = np.array([x['ranges'][2]['null_sd'] for x in M])
A = (ex * ms ** a / sd ** 2).sum() / ((ms ** a) ** 2 / sd ** 2).sum()
xs = np.logspace(np.log10(4), np.log10(55), 40)
o.append('<polyline points="%s" fill="none" stroke="%s" stroke-dasharray="5 3"/>' % (' '.join('%.1f,%.1f' % (sx2(v), sy2(A * v ** a)) for v in xs), BLUE))
for x in M:
    v, s = x['ranges'][2]['excess_variance'], x['ranges'][2]['null_sd']
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.4"/>' % (sx2(x['m500_median_1e13']), sy2(v - s), sx2(x['m500_median_1e13']), sy2(v + s), BLUE))
    o.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s"/>' % (sx2(x['m500_median_1e13']), sy2(v), BLUE))
    v2, s2 = x['ranges'][3]['excess_variance'], x['ranges'][3]['null_sd']
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#888" stroke-width="1"/>' % (sx2(x['m500_median_1e13']) + 6, sy2(v2 - s2), sx2(x['m500_median_1e13']) + 6, sy2(v2 + s2)))
    o.append('<circle cx="%.1f" cy="%.1f" r="3.2" fill="#fff" stroke="#888"/>' % (sx2(x['m500_median_1e13']) + 6, sy2(v2)))
o.append('<text x="%.1f" y="%.1f" fill="%s" font-size="11">∝ M^%.2f</text>' % (sx2(30), sy2(A * 30 ** a) - 26, BLUE, a))
o.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s"/><text x="%.1f" y="%.1f">inside R₅₀₀</text>' % (sx2(3.6), sy2(285), BLUE, sx2(3.6) + 10, sy2(285) + 4))
o.append('<circle cx="%.1f" cy="%.1f" r="3.2" fill="#fff" stroke="#888"/><text x="%.1f" y="%.1f">1 to 2 R₅₀₀</text>' % (sx2(3.6), sy2(262), sx2(3.6) + 10, sy2(262) + 4))
o.append('</svg>')
open(os.path.join(HERE, 'figure1.svg'), 'w', encoding='utf-8').write('\n'.join(o))
print('written figure1.svg', [(round(x['m500_median_1e13'], 1), x['clusters'], round(x['z_median'], 3)) for x in M], 'A', round(A, 2), 'a', a)
