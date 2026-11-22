"""Entry 167 draft, Figure 1 (SVG, no plotting library).  Left: the Giant Arc's published members (Lopez et al. 2022, Fig. 5, read
in ga_members.py) in the test's tangent plane, with the 200 and 300 Mpc bands.  Right: E in the 52 null bands, the Arc's band, and
the excess the published overdensity predicts.   Run: .venv/bin/python -I entries/167-desi-arc/figure.py  -> figure1.svg
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ga_members as G
res = json.load(open(os.path.join(HERE, 'results.json')))
gm = json.load(open(os.path.join(HERE, 'ga_members.json')))
nullE = np.array([n['E'] for n in res['primary']['null']['nodes']])
E = res['primary']['arc']['E']
c13 = [k['predicted_E'] for k in gm['bands']['200']['cases'] if k['delta'] == 1.3]
call = [k['predicted_E'] for k in gm['bands']['200']['cases']]
W, H = 900, 380
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="Helvetica, Arial, sans-serif" font-size="12">' % (W, H),
     '<rect width="%d" height="%d" fill="#fff"/>' % (W, H)]
# ---- left panel: x, y in Mpc
L, T, PW, PH = 60, 30, 380, 290
xr, yr = (-650, 650), (-496, 496)
sx = lambda v: L + (v - xr[0]) / (xr[1] - xr[0]) * PW
sy = lambda v: T + PH - (v - yr[0]) / (yr[1] - yr[0]) * PH
o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#dbe6f4"/>' % (sx(-500), sy(150), sx(500) - sx(-500), sy(-150) - sy(150)))
o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#a9c4e6"/>' % (sx(-500), sy(100), sx(500) - sx(-500), sy(-100) - sy(100)))
o.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#333"/>' % (L, T, PW, PH))
for v in (-500, -250, 0, 250, 500):
    o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#333"/><text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (sx(v), T + PH, sx(v), T + PH + 5, sx(v), T + PH + 19, str(v).replace('-', '\u2212')))
for v in (-400, -200, 0, 200, 400):
    o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#333"/><text x="%d" y="%.1f" text-anchor="end">%s</text>' % (L - 5, sy(v), L, sy(v), L - 8, sy(v) + 4, str(v).replace('-', '\u2212')))
for x, y, g in zip(G.x, G.y, G.grp):
    o.append('<circle cx="%.1f" cy="%.1f" r="3.2" fill="%s" stroke="%s"/>' % (sx(x), sy(y), '#111' if g == 'B' else 'none', '#111' if g == 'B' else '#c0392b'))
o.append('<text x="%.1f" y="%d" text-anchor="middle">x, east (Mpc)</text>' % (L + PW / 2, T + PH + 38))
o.append('<text transform="translate(16,%.1f) rotate(-90)" text-anchor="middle">y, north (Mpc)</text>' % (T + PH / 2))
o.append('<text x="%d" y="%d" font-weight="bold">a</text><text x="%d" y="%d">Published members: 44 main (filled), 11 secondary (open)</text>' % (L, T - 10, L + 16, T - 10))
o.append('<text x="%.1f" y="%.1f" fill="#1f4e8c">bands 200 and 300 Mpc wide: 37 and 41 of the 44</text>' % (sx(-630), sy(-260)))
# ---- right panel: histogram of null E
L2, PW2 = 520, 350
er = (-0.4, 0.8)
ex = lambda v: L2 + (v - er[0]) / (er[1] - er[0]) * PW2
bins = np.arange(-0.4, 0.4001, 0.05)
cnt, _ = np.histogram(nullE, bins)
ymax = cnt.max() + 1
ey = lambda c: T + PH - c / ymax * PH
o.append('<rect x="%.1f" y="%d" width="%.1f" height="%d" fill="#f6d9d5"/>' % (ex(min(call)), T, ex(max(call)) - ex(min(call)), PH))
o.append('<rect x="%.1f" y="%d" width="%.1f" height="%d" fill="#e9a79f"/>' % (ex(min(c13)), T, ex(max(c13)) - ex(min(c13)), PH))
for b0, c in zip(bins[:-1], cnt):
    if c:
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#8a8f98" stroke="#fff"/>' % (ex(b0), ey(c), ex(b0 + 0.05) - ex(b0), T + PH - ey(c)))
o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#1f4e8c" stroke-width="2.5"/>' % (ex(E), T, ex(E), T + PH))
o.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#333"/>' % (L2, T, PW2, PH))
for v in (-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8):
    o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#333"/><text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (ex(v), T + PH, ex(v), T + PH + 5, ex(v), T + PH + 19, ('%+.1f' % v).replace('-', '−') if v else '0'))
for c in range(0, ymax + 1, 2):
    o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#333"/><text x="%d" y="%.1f" text-anchor="end">%d</text>' % (L2 - 5, ey(c), L2, ey(c), L2 - 8, ey(c) + 4, c))
o.append('<text x="%.1f" y="%d" text-anchor="middle">excess E in the 1,000 × 200 Mpc band</text>' % (L2 + PW2 / 2, T + PH + 38))
o.append('<text transform="translate(%d,%.1f) rotate(-90)" text-anchor="middle">null bands</text>' % (L2 - 34, T + PH / 2))
o.append('<text x="%d" y="%d" font-weight="bold">b</text><text x="%d" y="%d">52 null bands (grey), the Arc’s band (blue), prediction (red)</text>' % (L2, T - 10, L2 + 16, T - 10))
o.append('<text x="%.1f" y="%d" fill="#1f4e8c" text-anchor="end">DESI: %s</text>' % (ex(E) - 6, T + 16, ('%+.2f' % E).replace('-', '−')))
o.append('<text x="%.1f" y="%d" fill="#a93226" text-anchor="middle">published</text><text x="%.1f" y="%d" fill="#a93226" text-anchor="middle">overdensity</text>' % (ex(np.mean(call)), T + 16, ex(np.mean(call)), T + 31))
o.append('</svg>')
open(os.path.join(HERE, 'figure1.svg'), 'w').write('\n'.join(o))
print('figure1.svg', len('\n'.join(o)), 'bytes; null n', len(nullE), 'pred', min(call), max(call))
