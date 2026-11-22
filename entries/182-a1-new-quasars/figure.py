"""Entry 182: figure1.svg - corrected polarization angles in region A1: the 40 quasars of 2005 measured after A1 was defined (top) and
the 20 measured later (bottom), with A1's sector (146-46 deg) shaded and the direction of dust and stars in the patch (58 deg).
Reads results.json (this folder) and, for the 2005 quasars, Entry 175's foreground cache and the 355-quasar table:
  ../175-quasar-polarization-dust/work/template.npz, data/hutsemekers-2005/qsos.tsv (VizieR J/A+A/441/915/qsos).
Run from the repository root:  .venv/bin/python -I entries/182-a1-new-quasars/figure.py
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qso_dust as q
from astropy.coordinates import SkyCoord, FK4
import astropy.units as u

R = json.load(open(os.path.join(HERE, 'results.json')))
new = np.array([d['pa1'] for d in R['A1']['quasars']])
T = np.load(os.path.join(HERE, '..', '175-quasar-polarization-dust', 'work', 'template.npz'))
Q = q.load_quasars()
head, rows = q.read_tsv(q.QSOS)
rp = np.array([int(r[head.index('r_Pol')]) for r in rows])
ra50 = SkyCoord(ra=Q['ra'] * u.deg, dec=Q['dec'] * u.deg).transform_to(FK4(equinox='B1950')).ra.hour
m = (ra50 >= 11.25) & (ra50 < 14.5) & (Q['z'] >= 1.0) & (Q['z'] <= 2.3) & (rp >= 7) & ~np.isin(np.array(Q['name']), ['B1222+228', 'B1246-057'])
old = q.correct(Q['p'], Q['pa'], T['q'], T['u'], 1.0)[1][m]
assert m.sum() == 40
W, H = 640, 420
L, Rm, top, ph, gap = 64, 20, 30, 150, 50
pw = W - L - Rm
X = lambda a: L + pw * a / 180.0
out = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="Helvetica, Arial, sans-serif" font-size="12">' % (W, H),
       '<rect width="100%" height="100%" fill="#ffffff"/>']
for k, (ang, title) in enumerate(((old, '2005: 40 quasars measured after A1 was defined (27 in the sector, P = 0.0028)'),
                                  (new, 'Measured later: 20 quasars (9 in the sector, P = 0.57)'))):
    y0 = top + k * (ph + gap)
    h, _ = np.histogram(ang, bins=18, range=(0, 180))
    ymax = 9
    Y = lambda c: y0 + ph - ph * c / ymax
    out.append('<rect x="%.1f" y="%d" width="%.1f" height="%d" fill="#e8eefc"/>' % (X(146), y0, X(180) - X(146), ph))
    out.append('<rect x="%.1f" y="%d" width="%.1f" height="%d" fill="#e8eefc"/>' % (X(0), y0, X(46) - X(0), ph))
    for i, c in enumerate(h):
        out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#2f5bb7"/>' % (X(10 * i) + 1, Y(c), X(10) - X(0) - 2, Y(0) - Y(c)))
    u0 = len(ang) / 18.0
    out.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="#555" stroke-dasharray="4 3"/>' % (L, L + pw, Y(u0), Y(u0)))
    out.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="#b5462f" stroke-width="2"/>' % (X(58), X(58), y0, y0 + ph))
    out.append('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="#222"/>' % (L, L + pw, y0 + ph, y0 + ph))
    out.append('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="#222"/>' % (L, L, y0, y0 + ph))
    for c in range(0, ymax + 1, 3):
        out.append('<text x="%d" y="%.1f" text-anchor="end" fill="#222">%d</text>' % (L - 6, Y(c) + 4, c))
    for a in range(0, 181, 30):
        out.append('<text x="%.1f" y="%d" text-anchor="middle" fill="#222">%d°</text>' % (X(a), y0 + ph + 15, a))
    out.append('<text x="%d" y="%d" fill="#222" font-weight="bold">%s</text>' % (L, y0 - 8, title))
    out.append('<text x="18" y="%.1f" transform="rotate(-90 18 %.1f)" text-anchor="middle" fill="#222">quasars per 10°</text>' % (y0 + ph / 2, y0 + ph / 2))
out.append('<text x="%.1f" y="%d" text-anchor="middle" fill="#222">polarization angle after removing the predicted dust foreground (east of north)</text>' % (L + pw / 2, H - 8))
out.append('<text x="%.1f" y="%d" fill="#b5462f">dust and stars, 58°</text>' % (X(58) + 4, top + 14))
out.append('<text x="%.1f" y="%d" fill="#2f5bb7">A1 sector</text>' % (X(147), top + 14))
out.append('</svg>')
open(os.path.join(HERE, 'figure1.svg'), 'w').write('\n'.join(out))
print('figure1.svg written; 2005 in sector', int(((old >= 146) | (old <= 46)).sum()), 'later', int(((new >= 146) | (new <= 46)).sum()))
