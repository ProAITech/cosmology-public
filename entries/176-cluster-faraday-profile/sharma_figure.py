"""Entry 176 follow-up, Figure 2 (SVG, no plotting library).  Left: the estimator of Sharma, Ravi & Krause on this sample, xi(R) with
jackknife errors (sharma.json, frozen), with a local background (sharma_posthoc.json, P2, after the result), and the spread of 1,000
shuffles.  Right: the significance of each statistic against clusters turned in Galactic longitude, the 2 x 2 of the follow-up rule
(filled: turned over the whole sky, fixed in advance; open: turned within the eRASS1 half, after the result), and S1 itself.
Run: python3 entries/176-cluster-faraday-profile/sharma_figure.py  -> figure2.svg
"""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
F = json.load(open(os.path.join(HERE, 'sharma.json')))
P = json.load(open(os.path.join(HERE, 'sharma_posthoc.json')))
M = lambda v, f='%.0f': (f % v).replace('-', '−')
BLUE, GREY, RED, INK, RULE, ORANGE = '#2f5bb7', '#8a8a8a', '#b5462f', '#222', '#e3e3e3', '#c77c1e'
W, H = 960, 470
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="Helvetica, Arial, sans-serif" font-size="12">' % (W, H),
     '<rect width="100%" height="100%" fill="#fff"/>']
# ---- left: xi(R)
L, T, PW, PH = 78, 44, 360, 300
yr = (-600.0, 2700.0)
ex = lambda R: L + (math.log10(R) - 1.0) / 2.0 * PW
sy = lambda v: T + PH - (v - yr[0]) / (yr[1] - yr[0]) * PH
S1, P2, sh = F['S1']['bins'], P['P2']['bins'], F['S1']['shuffle']
for j, b in enumerate(S1):                                     # shuffle spread, mean +- 1 sd, per bin
    m, s = sh['xi_mean'][j], sh['xi_sd'][j]
    o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#ececec"/>' % (ex(b['lo_kpc']), sy(min(m + s, yr[1])), ex(b['hi_kpc']) - ex(b['lo_kpc']),
                                                                                    sy(max(m - s, yr[0])) - sy(min(m + s, yr[1]))))
o.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#333"/>' % (L, T, PW, PH))
for v in (-500, 0, 500, 1000, 1500, 2000, 2500):
    o.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="%s"/>' % (L, L + PW, sy(v), sy(v), '#888' if v == 0 else RULE))
    o.append('<text x="%d" y="%.1f" text-anchor="end">%s</text>' % (L - 7, sy(v) + 4, '{:,}'.format(v).replace('-', '−')))
for v, lab in ((10, '10'), (30, '30'), (100, '100'), (300, '300'), (1000, '1,000')):
    o.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="#333"/><text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (ex(v), ex(v), T + PH, T + PH + 5, ex(v), T + PH + 19, lab))
o.append('<text x="%.1f" y="%d" text-anchor="middle">projected distance from the cluster, R (kpc)</text>' % (L + PW / 2, T + PH + 38))
o.append('<text transform="translate(18,%.1f) rotate(-90)" text-anchor="middle">ξ = ⟨ΔRM²⟩ − random (rad² m⁻⁴)</text>' % (T + PH / 2))
for pts, col, dx, filled in ((S1, BLUE, 0.97, True), (P2, ORANGE, 1.03, False)):
    for b in pts:
        R = math.sqrt(b['lo_kpc'] * b['hi_kpc']) * dx
        x, y = ex(R), sy(b['xi'])
        lo, hi = max(b['xi'] - b['sigma_jk'], yr[0]), min(b['xi'] + b['sigma_jk'], yr[1])
        kept = b['kept']
        col2 = col if kept else GREY
        o.append('<line x1="%.1f" x2="%.1f" y1="%.1f" y2="%.1f" stroke="%s" stroke-width="1.3"/>' % (x, x, sy(lo), sy(hi), col2))
        if filled:
            o.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s"/>' % (x, y, col2))
        else:
            o.append('<rect x="%.1f" y="%.1f" width="8" height="8" transform="rotate(45 %.1f %.1f)" fill="#fff" stroke="%s" stroke-width="1.5"/>' % (x - 4, y - 4, x, y, col2))
ly = T + PH + 56
o.append('<circle cx="%d" cy="%d" r="4.5" fill="%s"/><text x="%d" y="%d">fixed in advance: random term from turned clusters (S/N %s)</text>' % (L + 14, ly, BLUE, L + 24, ly + 4, M(F['S1']['snr'], '%.1f')))
o.append('<rect x="%d" y="%d" width="8" height="8" transform="rotate(45 %d %d)" fill="#fff" stroke="%s" stroke-width="1.5"/><text x="%d" y="%d">after the result: local background, 2.5–10 Mpc (S/N %s)</text>' % (
    L + 10, ly + 13, L + 14, ly + 17, ORANGE, L + 24, ly + 21, M(P['P2']['snr'], '%.1f')))
o.append('<rect x="%d" y="%d" width="12" height="9" fill="#ececec"/><text x="%d" y="%d">1,000 shuffles of ΔRM, mean ± 1 sd; grey point: fewer than 30 pairs</text>' % (L + 8, ly + 30, L + 24, ly + 38))
o.append('<text x="%d" y="24" fill="%s" font-size="13" font-weight="bold">a  Their estimator on 9,050 eRASS1 clusters</text>' % (L - 60, INK))
# ---- right: significance
TB = F['two_by_two']; TP = P['two_by_two_P1']
L2, T2, PW2 = 700, 44, 230
xr = (0.0, 11.0)
sx = lambda v: L2 + (v - xr[0]) / (xr[1] - xr[0]) * PW2
rows = [('clipped scale, range A (Entry 176)', TB['rangeA_mad']['t'], TP['rangeA_mad']['t']),
        ('clipped scale, their five bins', TB['physical_mad']['snr_rot'], TP['physical_mad']['snr_rot']),
        ('unclipped mean, range A', TB['rangeA_unclipped']['t'], TP['rangeA_unclipped']['t']),
        ('unclipped mean, their five bins', TB['physical_unclipped']['snr_rot'], TP['physical_unclipped']['snr_rot'])]
rh = 38
PH2 = rh * 6
o.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#333"/>' % (L2, T2, PW2, PH2))
for v in range(0, 12, 2):
    o.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s"/>' % (sx(v), sx(v), T2, T2 + PH2, RULE if v else '#888'))
    o.append('<text x="%.1f" y="%d" text-anchor="middle">%d</text>' % (sx(v), T2 + PH2 + 18, v))
o.append('<text x="%.1f" y="%d" text-anchor="middle">significance against turned clusters (σ)</text>' % (L2 + PW2 / 2, T2 + PH2 + 38))
for i, (lab, a, b) in enumerate(rows):
    y = T2 + rh * (i + 0.5)
    o.append('<text x="%d" y="%.1f" text-anchor="end" fill="%s">%s</text>' % (L2 - 8, y + 4, INK, lab))
    o.append('<circle cx="%.1f" cy="%.1f" r="5" fill="%s"/>' % (sx(a), y, BLUE))
    o.append('<circle cx="%.1f" cy="%.1f" r="5" fill="#fff" stroke="%s" stroke-width="1.5"/>' % (sx(b), y, BLUE))
y = T2 + rh * 4.0
o.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="#999" stroke-dasharray="3 3"/>' % (L2 - 230, L2 + PW2, y, y))
y = T2 + rh * 4.5
o.append('<text x="%d" y="%.1f" text-anchor="end" fill="%s">their S/N, jackknife</text>' % (L2 - 8, y + 4, INK))
o.append('<circle cx="%.1f" cy="%.1f" r="5" fill="%s"/>' % (sx(F['S1']['snr']), y, BLUE))
o.append('<circle cx="%.1f" cy="%.1f" r="5" fill="#fff" stroke="%s" stroke-width="1.5"/>' % (sx(P['P1']['snr']), y, BLUE))
o.append('<rect x="%.1f" y="%.1f" width="8" height="8" transform="rotate(45 %.1f %.1f)" fill="#fff" stroke="%s" stroke-width="1.5"/>' % (sx(P['P2']['snr']) - 4, y - 4, sx(P['P2']['snr']), y, ORANGE))
y = T2 + rh * 5.5
o.append('<text x="%d" y="%.1f" text-anchor="end" fill="%s">Sharma et al., 7,775 ROSAT clusters</text>' % (L2 - 8, y + 4, INK))
o.append('<circle cx="%.1f" cy="%.1f" r="5" fill="%s"/>' % (sx(3.6), y, RED))
y = T2 + PH2 + 62
o.append('<circle cx="%d" cy="%.1f" r="5" fill="%s"/><text x="%d" y="%.1f">turned over the whole sky (fixed in advance)</text>' % (L2 - 120, y, BLUE, L2 - 110, y + 4))
o.append('<circle cx="%d" cy="%.1f" r="5" fill="#fff" stroke="%s" stroke-width="1.5"/><text x="%d" y="%.1f">turned within the eRASS1 half (after the result)</text>' % (L2 - 120, y + 17, BLUE, L2 - 110, y + 21))
o.append('<text x="%d" y="24" fill="%s" font-size="13" font-weight="bold">b  What the significance depends on</text>' % (L2 - 230, INK))
o.append('</svg>')
open(os.path.join(HERE, 'figure2.svg'), 'w', encoding='utf-8').write('\n'.join(o)); print('figure2.svg written')
