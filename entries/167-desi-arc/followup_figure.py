"""Entry 167 follow-up, Figure 2 (SVG, no plotting library).  Left: the Giant Arc's 44 published main members in the test's tangent
plane, marked by what DESI DR1 shows on them (followup.json, part B).  Right: the band excess E with all absorbers (Entry 167) and with
strong absorbers only (part A), each against the 52 null bands (mean +- 1 and 2 sd), beside the excess the published overdensity implies
(ga_members.json).   Run: .venv/bin/python -I entries/167-desi-arc/followup_figure.py  -> figure2.svg
"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
F = json.load(open(os.path.join(HERE, 'followup.json')))
R = json.load(open(os.path.join(HERE, 'results.json')))
G = json.load(open(os.path.join(HERE, 'ga_members.json')))['bands']['200']['cases']
M = lambda v: ('%+.2f' % v).replace('-', '−') if abs(v) > 1e-9 else '0'
BLUE, GREY, RED, INK, RULE, BAND = '#2f5bb7', '#8a8a8a', '#b5462f', '#222', '#e3e3e3', '#dbe6f4'
W, H = 960, 400
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="Helvetica, Arial, sans-serif" font-size="12">' % (W, H),
     '<rect width="100%" height="100%" fill="#fff"/>']
# ---- left: members
L, T, PW, PH = 60, 40, 380, 290
xr, yr = (-650, 650), (-496, 496)
sx = lambda v: L + (v - xr[0]) / (xr[1] - xr[0]) * PW
sy = lambda v: T + PH - (v - yr[0]) / (yr[1] - yr[0]) * PH
o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>' % (sx(-500), sy(100), sx(500) - sx(-500), sy(-100) - sy(100), BAND))
o.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#333"/>' % (L, T, PW, PH))
for v in (-500, -250, 0, 250, 500):
    o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#333"/><text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (sx(v), T + PH, sx(v), T + PH + 5, sx(v), T + PH + 19, str(v).replace('-', '−')))
for v in (-400, -200, 0, 200, 400):
    o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#333"/><text x="%d" y="%.1f" text-anchor="end">%s</text>' % (L - 5, sy(v), L, sy(v), L - 8, sy(v) + 4, str(v).replace('-', '−')))
B = F['B_members']
for m in B['members']:
    x, y = sx(max(xr[0] + 5, min(xr[1] - 5, m['x_mpc']))), sy(m['y_mpc'])
    if m.get('recovered'):
        o.append('<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/>' % (x, y, BLUE))
    elif m['matched']:
        o.append('<circle cx="%.1f" cy="%.1f" r="4" fill="#fff" stroke="%s" stroke-width="1.5"/>' % (x, y, BLUE))
    else:
        o.append('<path d="M%.1f %.1fl7 7m0 -7l-7 7" stroke="%s" stroke-width="1.5"/>' % (x - 3.5, y - 3.5, GREY))
o.append('<text x="%.1f" y="%d" text-anchor="middle">x, east (Mpc)</text>' % (L + PW / 2, T + PH + 38))
o.append('<text transform="translate(16,%.1f) rotate(-90)" text-anchor="middle">y, north (Mpc)</text>' % (T + PH / 2))
o.append('<text x="%d" y="22" fill="%s" font-size="13" font-weight="bold">a  The 44 published members on DESI DR1</text>' % (L - 48, INK))
ly = T + PH - 52
o.append('<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/><text x="%.1f" y="%.1f">absorber found again, z 0.742–0.862 (%d)</text>' % (L + 14, ly, BLUE, L + 24, ly + 4, B['N_r']))
o.append('<circle cx="%.1f" cy="%.1f" r="4" fill="#fff" stroke="%s" stroke-width="1.5"/><text x="%.1f" y="%.1f">DESI quasar within 1′, no absorber listed (%d)</text>' % (L + 14, ly + 17, BLUE, L + 24, ly + 21, B['N_m'] - B['N_r']))
o.append('<path d="M%.1f %.1fl7 7m0 -7l-7 7" stroke="%s" stroke-width="1.5"/><text x="%.1f" y="%.1f">no DESI quasar within 1′ (%d)</text>' % (L + 10.5, ly + 30.5, GREY, L + 24, ly + 38, 44 - B['N_m']))
o.append('<text x="%.1f" y="%.1f" fill="#1f4e8c">test band, 1,000 × 200 Mpc</text>' % (sx(-630), sy(-130)))
# ---- right: E against the null
L2, T2, PW2 = 640, 40, 290
er = (-0.6, 0.8)
ex = lambda v: L2 + (v - er[0]) / (er[1] - er[0]) * PW2
rows = [('all absorbers (Entry 167)', R['primary']['arc']['E'], R['primary']['null']['mean'], R['primary']['null']['sd'], 'z 0.742–0.862', BLUE, False),
        ('EW ≥ 0.6 Å (decides)', F['A_ew06']['arc']['E'], F['A_ew06']['null']['mean'], F['A_ew06']['null']['sd'], 'z 0.742–0.862', BLUE, True),
        ('EW ≥ 1.0 Å (check)', F['check_ew10']['arc']['E'], F['check_ew10']['null']['mean'], F['check_ew10']['null']['sd'], 'z 0.742–0.862', BLUE, False),
        ('EW ≥ 0.6 Å, slice below', F['check_ew06_slice_below']['arc']['E'], F['check_ew06_slice_below']['null']['mean'], F['check_ew06_slice_below']['null']['sd'], 'z 0.630–0.742', GREY, False)]
rh = 62; PH2 = rh * len(rows)
pmin, pmax = min(k['predicted_E'] for k in G), max(k['predicted_E'] for k in G)
o.append('<rect x="%.1f" y="%d" width="%.1f" height="%d" fill="#f6e3de"/>' % (ex(pmin), T2, ex(pmax) - ex(pmin), rh * 3))
o.append('<text x="%.1f" y="%d" text-anchor="middle" fill="%s" font-size="11">published members imply</text>' % ((ex(pmin) + ex(pmax)) / 2, T2 - 6, RED))
for v in (-0.4, -0.2, 0, 0.2, 0.4, 0.6):
    o.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s"/>' % (ex(v), ex(v), T2, T2 + PH2, '#888' if v == 0 else RULE))
    o.append('<text x="%.1f" y="%d" text-anchor="middle" fill="#333">%s</text>' % (ex(v), T2 + PH2 + 34, M(v)))
for i, (lab, E, m, sd, sl, col, bold) in enumerate(rows):
    y = T2 + rh * (i + 0.5)
    o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="14" fill="#ececec"/>' % (ex(m - 2 * sd), y - 7, ex(m + 2 * sd) - ex(m - 2 * sd)))
    o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="14" fill="#cfcfcf"/>' % (ex(m - sd), y - 7, ex(m + sd) - ex(m - sd)))
    o.append('<circle cx="%.1f" cy="%.1f" r="%s" fill="%s"/>' % (ex(E), y, 6 if bold else 4.5, col))
    o.append('<text x="%d" y="%.1f" text-anchor="end" fill="%s"%s>%s</text>' % (L2 - 10, y, INK, ' font-weight="bold"' if bold else '', lab))
    o.append('<text x="%d" y="%.1f" text-anchor="end" fill="#666" font-size="11">%s · E %s</text>' % (L2 - 10, y + 15, sl, M(E)))
o.append('<text x="%.1f" y="%d" text-anchor="middle" fill="#333">band excess E (dark grey: null ±1 sd; light: ±2 sd)</text>' % (L2 + PW2 / 2, T2 + PH2 + 54))
o.append('<text x="%d" y="22" fill="%s" font-size="13" font-weight="bold">b  Excess in the band, against 52 bands elsewhere</text>' % (L2 - 160, INK))
o.append('</svg>')
open(os.path.join(HERE, 'figure2.svg'), 'w', encoding='utf-8').write('\n'.join(o)); print('figure2.svg written')
