"""Entry 182: region flags, Galactic latitude, the 10-arcsecond match to the 355 quasars, duplicate groups and Planck's predicted
foreground for the rows in positions.txt (read in browser memory from VizieR J/A+A/606/A101/table4 and J/A+A/620/A68/table3).
Prints the object that a1_new.js part2() takes.  Uses Entry 175's code (copied unchanged: qso_dust.py) and its data:
  data/planck-pr3-gnilc-dust/COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits (IRSA, 1,811,949,120 bytes)
  data/hutsemekers-2005/qsos.tsv (VizieR J/A+A/441/915/qsos)
Run from the repository root:  .venv/bin/python -I entries/182-a1-new-quasars/foreground.py
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qso_dust as q
from astropy.coordinates import SkyCoord, FK4
import astropy.units as u

rows = []
for line in open(os.path.join(HERE, 'positions.txt')):
    if line.startswith('#'):
        continue
    for item in line.strip().split(';'):
        if item.strip():
            i, ra, de, z = item.split(); rows.append((int(i), float(ra), float(de), float(z)))
idx = np.array([r[0] for r in rows]); ra = np.array([r[1] for r in rows]); dec = np.array([r[2] for r in rows]); z = np.array([r[3] for r in rows])
c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
ra50 = c.transform_to(FK4(equinox='B1950')).ra.hour
absb = np.abs(c.galactic.b.deg)
A1 = (ra50 >= 11.25) & (ra50 < 14.5) & (z >= 1.0) & (z <= 2.3)
A3 = (ra50 >= 21 + 20 / 60) & (ra50 < 24.0) & (z >= 0.7) & (z <= 1.5)
Q = q.load_quasars()
c355 = SkyCoord(ra=Q['ra'] * u.deg, dec=Q['dec'] * u.deg)
sep355 = c.match_to_catalog_sky(c355)[1].arcsec
match = sep355 < 10
grp = np.arange(len(rows))
for a in range(len(rows)):
    for b in range(a):
        if c[a].separation(c[b]).arcsec < 10:
            grp[a] = grp[b]; break
T, _ = q.template_at(ra, dec)
out = {int(idx[k]): [int(A1[k]), int(A3[k]), round(float(absb[k]), 2), int(match[k]), round(float(T['q'][k]), 4), round(float(T['u'][k]), 4), int(idx[grp[k]])]
       for k in range(len(rows))}
json.dump(dict(rows=len(rows), in_A1=int(A1.sum()), in_A3=int(A3.sum()), A1_highlat=int((A1 & (absb >= 30)).sum()),
               A1_highlat_new=int((A1 & (absb >= 30) & ~match).sum()), A3_highlat_new=int((A3 & (absb >= 30) & ~match).sum()),
               matched_355=int(match.sum()), duplicate_rows=int((grp != np.arange(len(rows))).sum()), fg=out),
          open(os.path.join(HERE, 'foreground.json'), 'w'))
print('rows %d; A1 %d (|b|>=30: %d, not among the 355: %d); A3 %d (new, |b|>=30: %d); matched to the 355: %d; duplicate rows %d' % (
    len(rows), A1.sum(), (A1 & (absb >= 30)).sum(), (A1 & (absb >= 30) & ~match).sum(), A3.sum(), (A3 & (absb >= 30) & ~match).sum(), match.sum(),
    (grp != np.arange(len(rows))).sum()))
print('matched to the 355:', [(int(idx[k]), round(float(sep355[k]), 1)) for k in np.nonzero(match)[0]])
print('FG=' + json.dumps({k: v for k, v in out.items() if v[0] or v[1]}, separators=(',', ':')))
