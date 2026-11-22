"""Entry 167, added after the result (not part of the frozen rule): where the Giant Arc's published members lie in the test band.

Positions were read on 21 November 2026 from the vector graphic of Figure 5 of Lopez, Clowes & Williger (2022), arXiv:2201.06875
(ar5iv asset CHMS_GA_bigger_FOV.svg), in browser memory: 44 black circles (the main group) and 11 red (the secondary group).
The axis labels are drawn as glyphs; their repetition pattern reads 200, 205, ..., 225 on x and 30, 35, 40, 45 on y, taken as
RA and Dec in degrees. With that reading the main group's RA midpoint is 214.8, beside the cluster's 215.4.
This script projects them as the frozen rule does (gnomonic about RA 215.405, Dec +38.361, times 2,837 Mpc), counts them in the
band, and turns the published overdensity into a predicted band excess.   Run: .venv/bin/python -I entries/167-desi-arc/ga_members.py
"""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = ("B 206.4 38.39;B 208.95 38.21;B 210.29 38.24;B 212.38 38.64;B 212.6 38.73;B 212.66 37.15;B 215.81 37.87;B 216.25 38.71;"
       "B 216.26 38.35;B 219.39 39;B 219.59 39.57;B 221.06 38.98;B 221.42 40.48;B 222.16 42.2;B 223.71 39.71;B 223.76 42.02;"
       "B 223.86 39.83;B 205.03 38.56;B 207.79 38.28;B 210.09 38.11;B 210.11 39.53;B 211.22 39.94;B 211.75 38.46;B 211.75 38.9;"
       "B 211.86 37.48;B 212.03 37.75;B 213.07 36.24;B 213.16 38.29;B 213.51 37.76;B 214 37.82;B 214.67 39.38;B 216.41 37.96;"
       "B 217.13 39.47;B 217.79 38.74;B 218.17 38.54;B 218.33 38.93;B 222.35 40.1;B 222.41 39.27;B 222.48 42.1;B 222.66 40.09;"
       "B 222.66 40.41;B 222.86 39.77;B 223.36 39.62;B 224.6 40.43;R 202.53 38.01;R 202.96 39.57;R 203.2 39.97;R 203.46 38.96;"
       "R 201.24 41.08;R 202.26 38.44;R 202.99 39.91;R 203.32 40.31;R 203.84 39.57;R 204.68 41.43;R 205.8 38.86")
D = np.pi / 180
rows = [r.split() for r in RAW.split(';')]
grp = np.array([r[0] for r in rows]); ra = np.array([float(r[1]) for r in rows]) * D; de = np.array([float(r[2]) for r in rows]) * D
a0, d0 = 215.405 * D, 38.361 * D
cc = np.sin(d0) * np.sin(de) + np.cos(d0) * np.cos(de) * np.cos(ra - a0)
x = np.cos(de) * np.sin(ra - a0) / cc * 2837
y = (np.cos(d0) * np.sin(de) - np.sin(d0) * np.cos(de) * np.cos(ra - a0)) / cc * 2837
res = json.load(open(os.path.join(HERE, 'results.json')))
meas = {200: (res['primary']['arc']['E'], res['primary']['null']['mean'], res['primary']['null']['sd']),
        300: (res['checks']['wide_band_300']['arc']['E'], res['checks']['wide_band_300']['null_mean'], res['checks']['wide_band_300']['null_sd'])}
out = {'main': int((grp == 'B').sum()), 'secondary': int((grp == 'R').sum()),
       'main_x_range': [float(x[grp == 'B'].min()), float(x[grp == 'B'].max())], 'main_y_range': [float(y[grp == 'B'].min()), float(y[grp == 'B'].max())], 'bands': {}}
for w in (200, 300):
    inb = (np.abs(x) <= 500) & (np.abs(y) <= w / 2)
    nB, nR = int((inb & (grp == 'B')).sum()), int((inb & (grp == 'R')).sum())
    E, mu, sd = meas[w]
    b = dict(main_in_band=nB, f=nB / 44, secondary_in_band=nR, measured_E=E, null_mean=mu, null_sd=sd, cases=[])
    for nfield in (504, 630):                              # field count: Sawala et al.; or 44 / 0.07 from the discovery's "about 7 per cent"
        nexp = nfield * 1000 * w / (1541 * 1615)
        for delta in (1.0, 1.3, 1.6):
            pred = 44 * delta / (1 + delta) * (nB / 44) / nexp
            b['cases'].append(dict(field=nfield, delta=delta, expected_in_band=round(nexp, 1), predicted_E=round(pred, 3), sigma_above=round((pred - (E - mu)) / sd, 2)))
    out['bands'][str(w)] = b
    c = [k for k in b['cases'] if k['delta'] == 1.3]
    print('band %d: %d of 44 main members (f = %.2f), %d of 11 secondary; predicted E %.2f to %.2f; %.1f to %.1f sd above the DESI value'
          % (w, nB, nB / 44, nR, min(k['predicted_E'] for k in c), max(k['predicted_E'] for k in c), min(k['sigma_above'] for k in c), max(k['sigma_above'] for k in c)))
    print('   all cases: sigma %.1f to %.1f' % (min(k['sigma_above'] for k in b['cases']), max(k['sigma_above'] for k in b['cases'])))
json.dump(out, open(os.path.join(HERE, 'ga_members.json'), 'w'), indent=1)
