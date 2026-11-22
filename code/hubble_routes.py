"""hubble_routes.py - Workbench Entry 106: the early-universe routes to the Hubble gap, scored with 2025-26 data.

Every value is a published posterior (mean +- 1 sigma, or asymmetric), with the data named. The targets are the ten one-step
direct measurements of Entry 56 (71.4 +- 1.2), SH0ES (73.17 +- 0.86) and CCHP (70.39 +- 1.94). For each route:
the tension with each target, using the error on the side facing it, and the share of the gap it closes against the
standard model fitted to the same data.
"""
import json, os, math
HERE = os.path.dirname(os.path.abspath(__file__))
TARGETS = {'direct': (71.4, 1.2), 'SH0ES': (73.17, 0.86), 'CCHP': (70.39, 1.94)}
# route, label, data, H0, err_minus, err_plus, LCDM H0 on the same data, extra numbers, source key
R = [
 ('none', 'Standard model', 'Planck + ACT DR6 + DESI DR2 + Pantheon+', 68.35, 0.27, 0.27, 68.35, 0, 'toda'),
 ('none', 'Standard model', 'Planck + ACT + SPT + DESI DR2', 68.27, 0.26, 0.26, 68.27, 0, 'lynch'),
 ('neff', 'Extra radiation (N_eff free)', 'Planck + ACT + SPT + DESI DR2 + helium + deuterium', 68.16, 0.50, 0.50, 68.27, 1, 'goldstein'),
 ('clump', 'Clumping by primordial fields', 'Planck + ACT + SPT + DESI DR2', 68.6, 0.3, 0.3, 68.27, 2, 'lynch'),
 ('clump', 'Clumping by primordial fields', 'Planck + ACT DR6 + DESI DR1', 69.1, 0.5, 0.5, 68.2, 1, 'act'),
 ('ede', 'Early dark energy', 'Planck + ACT DR6 + DESI DR2 + Pantheon+', 69.17, 0.90, 0.38, 68.35, 3, 'toda'),
 ('ede', 'Early dark energy', 'Planck + ACT DR6 + DESI DR1', 69.9, 1.5, 0.8, 68.2, 3, 'act'),
 ('ede', 'Early dark energy (profile likelihood)', 'Planck + ACT DR6 + DESI DR2 + Pantheon+', 71.0, 1.1, 1.1, 68.35, 3, 'poulin'),
 ('rec', 'Free-form recombination', 'Planck + ACT DR6 + DESI DR1', 69.6, 1.0, 1.0, 68.2, 4, 'act'),
 ('me', 'Heavier electron', 'Planck + ACT DR6 + DESI DR2 + Pantheon+', 69.41, 0.68, 0.68, 68.35, 1, 'toda'),
 ('me', 'Heavier electron', 'Planck + ACT DR6 + DESI DR1', 69.8, 1.1, 1.1, 68.2, 1, 'act'),
 ('mek', 'Heavier electron + curvature', 'Planck + ACT DR6 + DESI DR1', 71.0, 1.7, 1.7, 68.2, 2, 'act'),
 ('mek', 'Heavier electron + curvature', 'Planck + DESI DR1 + Pantheon+', 71.61, 1.00, 1.00, 68.2, 2, 'sv'),
]
out = {'targets': TARGETS, 'rows': []}
for key, lab, data, h, em, ep, h_l, k, src in R:
    row = dict(route=key, label=lab, data=data, H0=h, em=em, ep=ep, k=k, src=src)
    for t, (ht, et) in TARGETS.items():
        e = ep if h < ht else em
        row['sig_'+t] = round((ht - h)/math.hypot(e, et), 2)
    row['gap_closed_direct'] = round((h - h_l)/(71.4 - h_l), 2) if key != 'none' else 0.0
    out['rows'].append(row)
# the best route to each target, in the DESI DR2 era (posteriors only, no profile)
dr2 = [r for r in out['rows'] if 'DR2' in r['data'] and 'profile' not in r['label'] and r['route'] != 'none']
out['dr2_best'] = max(dr2, key=lambda r: r['H0'])['label']
out['dr2_range'] = [min(r['H0'] for r in dr2), max(r['H0'] for r in dr2)]
# curvature in the m_e + Omega_k fit
out['mek_act'] = dict(me=[1.022, 0.016], Ok=[-0.0031, 0.0037])
json.dump(out, open(os.path.join(HERE, 'hubble_routes.json'), 'w'), indent=1)
for r in out['rows']: print(f"{r['label']:42s} {r['data']:52s} {r['H0']:6.2f} -{r['em']} +{r['ep']}  direct {r['sig_direct']:+.2f}  SH0ES {r['sig_SH0ES']:+.2f}  CCHP {r['sig_CCHP']:+.2f}  gap {r['gap_closed_direct']}")
print(out['dr2_best'], out['dr2_range'])
