"""Entry 175, exploratory: the alignment and the link against the removal level f, with random-direction controls at every f.

Made after the frozen result was seen; nothing here decides anything.  Reads work/template.npz (written by qso_dust.py run), so
the dust map is not opened again.  Writes explore.json and figure1.svg.
Run from the repository root:  .venv/bin/python -I entries/175-quasar-polarization-dust/explore.py     (about 6 minutes)
Shuffles: 10,000 for each true-direction row, 1,000 for each random-direction trial.
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qso_dust as q

z = np.load(os.path.join(HERE, 'work', 'template.npz'))
ra, dec, zz, p, pa, tq, tu, tpa = z['ra'], z['dec'], z['z'], z['qp'], z['qpa'], z['q'], z['u'], z['pa']
n = len(p)
nvs = q.CFG['nvs']
prep = q.prepare(ra, dec, zz, max(nvs))
g = np.random.default_rng(q.CFG['seed'] + 123)
rows = []
for f in [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0]:
    p2, pa2 = q.correct(p, pa, tq, tu, f)
    r = q.align_test(prep, pa2 * q.D2R, nvs, 10000, q.CFG['seed'] + 1000 + int(round(100 * f)))
    c, e = q.link(pa2, tpa)
    ctl = []
    if f > 0:
        for k in range(8):
            p3, pa3 = q.correct(p, pa, tq, tu, f, phase=g.uniform(0, 2 * np.pi, n))
            ctl.append(q.align_test(prep, pa3 * q.D2R, nvs, 1000, q.CFG['seed'] + 2000 + k)['sigma'][2])
    rows.append(dict(f=f, sigma=r['sigma'], sigma_n30=r['sigma'][2], p_global=r['p_global'], S_global=r['sigma_global'], link=c, link_err=e,
                     control_n30=ctl, below_cut=int((p2 < q.CFG['pcut']).sum())))
    print('f %.2f: n=30 %.2f, global %.2f sigma (p %.4f), link %+.3f, controls %s' % (
        f, r['sigma'][2], r['sigma_global'], r['p_global'], c, [round(x, 2) for x in ctl]), flush=True)
# how the predicted foreground compares with the quasars' own polarization
ratio = np.hypot(tq, tu) / p
out = dict(rows=rows, foreground_over_quasar=dict(median=float(np.median(ratio)), above_quarter=float((ratio > 0.25).mean()),
                                                   above_half=float((ratio > 0.5).mean())))
json.dump(out, open(os.path.join(HERE, 'explore.json'), 'w'), indent=1)

# ---- figure
W, H = 900, 380
INK, GRID, BLUE, GREY, RED = '#222', '#bbb', '#1f4e9c', '#999', '#b03a2e'
mn = lambda v: ('%s' % v).replace('-', '−')
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="Helvetica, Arial, sans-serif" font-size="12">' % (W, H),
     '<rect width="%d" height="%d" fill="#fff"/>' % (W, H)]
def panel(L, T, PW, PH, xr, yr, xt, yt, xl, yl):
    sx = lambda v: L + (v - xr[0]) / (xr[1] - xr[0]) * PW
    sy = lambda v: T + PH - (v - yr[0]) / (yr[1] - yr[0]) * PH
    o.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="%s"/>' % (L, T, PW, PH, INK))
    for v in xt:
        o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="%s"/><text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (sx(v), T + PH, sx(v), T + PH + 5, INK, sx(v), T + PH + 19, mn(v)))
    for v in yt:
        o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s"/><text x="%d" y="%.1f" text-anchor="end">%s</text>' % (L - 5, sy(v), L, sy(v), INK, L - 8, sy(v) + 4, mn(v)))
        o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" stroke-width="0.4"/>' % (L, sy(v), L + PW, sy(v), GRID))
    o.append('<text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (L + PW / 2, T + PH + 38, xl))
    o.append('<text transform="translate(%d,%.1f) rotate(-90)" text-anchor="middle">%s</text>' % (L - 40, T + PH / 2, yl))
    return sx, sy
fs = [r['f'] for r in rows]
sx, sy = panel(60, 28, 380, 300, (-0.1, 3.1), (-0.5, 5.5), (0, 0.5, 1, 1.5, 2, 2.5, 3), (0, 1, 2, 3, 4, 5),
               'foreground removed, in units of Planck’s prediction (f)', 'alignment (σ)')
ctl = [r for r in rows if r['control_n30']]
top = ' '.join('%.1f,%.1f' % (sx(r['f']), sy(max(r['control_n30']))) for r in ctl)
bot = ' '.join('%.1f,%.1f' % (sx(r['f']), sy(min(r['control_n30']))) for r in reversed(ctl))
o.append('<polygon points="%s %s" fill="#ddd" stroke="none"/>' % (top, bot))
o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-dasharray="4 3"/>' % (sx(-0.1), sy(3), sx(3.1), sy(3), RED))
o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-dasharray="4 3"/>' % (sx(-0.1), sy(2), sx(3.1), sy(2), RED))
o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.8"/>' % (' '.join('%.1f,%.1f' % (sx(r['f']), sy(r['sigma_n30'])) for r in rows), BLUE))
for r in rows:
    o.append('<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/>' % (sx(r['f']), sy(r['sigma_n30']), BLUE))
o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.4" stroke-dasharray="2 3"/>' % (' '.join('%.1f,%.1f' % (sx(r['f']), sy(max(r['S_global'], -0.5))) for r in rows), INK))
o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (sx(1), sy(-0.5), sx(1), sy(5.5), GREY))
o.append('<text x="%.1f" y="%.1f" fill="#555" font-size="11">Planck’s prediction</text>' % (sx(1) + 5, sy(5.5) + 14))
o.append('<text x="%.1f" y="%.1f" fill="%s" font-size="11">3σ</text><text x="%.1f" y="%.1f" fill="%s" font-size="11">2σ</text>' % (sx(3.1) - 22, sy(3) - 4, RED, sx(3.1) - 22, sy(2) - 4, RED))
lx, ly = sx(0.08), sy(0.95)
o.append('<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/><text x="%.1f" y="%.1f">30 neighbours, true direction</text>' % (lx, ly - 4, BLUE, lx + 10, ly))
o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.4" stroke-dasharray="2 3"/><text x="%.1f" y="%.1f">global, five neighbour counts</text>' % (lx - 8, ly + 13, lx + 6, ly + 13, INK, lx + 10, ly + 17))
o.append('<rect x="%.1f" y="%.1f" width="14" height="9" fill="#ddd"/><text x="%.1f" y="%.1f">random directions (8 trials)</text>' % (lx - 8, ly + 26, lx + 10, ly + 34))
sx2, sy2 = panel(545, 28, 320, 300, (-0.1, 3.1), (-0.2, 0.3), (0, 0.5, 1, 1.5, 2, 2.5, 3), (-0.2, -0.1, 0, 0.1, 0.2, 0.3),
                 'foreground removed (f)', 'link to the foreground direction, mean cos 2Δ')
e = rows[0]['link_err']
o.append('<polygon points="%s %s" fill="#dbe6f4" stroke="none"/>' % (' '.join('%.1f,%.1f' % (sx2(r['f']), sy2(r['link'] + e)) for r in rows),
                                                                    ' '.join('%.1f,%.1f' % (sx2(r['f']), sy2(r['link'] - e)) for r in reversed(rows))))
o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (sx2(-0.1), sy2(0), sx2(3.1), sy2(0), INK))
o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.8"/>' % (' '.join('%.1f,%.1f' % (sx2(r['f']), sy2(r['link'])) for r in rows), BLUE))
for r in rows:
    o.append('<circle cx="%.1f" cy="%.1f" r="3.5" fill="%s"/>' % (sx2(r['f']), sy2(r['link']), BLUE))
o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (sx2(1), sy2(-0.2), sx2(1), sy2(0.3), GREY))
o.append('</svg>')
open(os.path.join(HERE, 'figure1.svg'), 'w', encoding='utf-8').write('\n'.join(o))
print('written explore.json, figure1.svg;', json.dumps(out['foreground_over_quasar']))
