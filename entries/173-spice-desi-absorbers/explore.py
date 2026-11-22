"""Entry 173, exploratory checks made after the frozen result was seen.  Not part of the rule; nothing here decides anything.

Reads work/sample.npz (written by absorber_split.py run) and Milliquas v8, and asks which DESI sight lines carry the growth:
  (1) sight lines that Entry 171 also had (a Milliquas Type Q quasar on the same SPICE-RACS source) against DESI-only ones;
  (2) the two parts of the frozen 'good spectra' check taken one at a time (z >= 0.4; TSNR2_QSO at or above A's lower quartile);
  (3) Entry 171's own sample split by whether DESI has the quasar, and the union of the two samples.
Run from the repository root:  .venv/bin/python -I entries/173-spice-desi-absorbers/explore.py
"""
import json, os, sys
import numpy as np
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import absorber_split as a

log = lambda s: print(s, flush=True)
S = dict(np.load(os.path.join(a.WORK, 'sample.npz')))
rows, src, _ = a.load_spice()
mq = a.load_milliquas()
tree = cKDTree(a.unit(src['ra'], src['dec']))
mm = a.match(mq, rows)
# Entry 171's kept sample, rebuilt: Milliquas Type Q, |b| >= 20, >= 10 neighbours
idx = mm['src']
lat = np.abs(src['b'][idx]) >= a.CFG['bcut']
med, nn = a.foreground(src, idx[lat], a.CFG['fg_radius'], a.CFG['min_neigh'], tree)
ok = np.isfinite(med)
k = idx[lat][ok]
M = dict(z=mq['z'][mm['q'][lat][ok]], d=src['rm'][k] - med[ok], e=src['e'][k], b=src['b'][k], l=src['l'][k], src=k)
log('Entry 171 sample rebuilt: %d' % len(k))
out = {}
N = 300
inM = np.isin(S['src'], M['src'])
out['desi_also_in_entry171'] = a.analyse(a.sub(S, inM), N, 'DESI sight lines also in Entry 171', log=log)
out['desi_only'] = a.analyse(a.sub(S, ~inM), N, 'DESI-only sight lines', log=log)
cut = float(np.quantile(S['tsnr2'][S['A']], 0.25))
out['z_ge_0p4'] = a.analyse(a.sub(S, S['z'] >= 0.4), N, 'D, z >= 0.4', log=log)
out['tsnr2_ge_cut'] = a.analyse(a.sub(S, S['tsnr2'] >= cut), N, 'D, TSNR2 >= %.1f' % cut, log=log)
out['tsnr2_lt_cut'] = a.analyse(a.sub(S, S['tsnr2'] < cut), N, 'D, TSNR2 < %.1f' % cut, log=log)
inD = np.isin(M['src'], S['src'])
out['entry171_in_desi'] = a.analyse(a.sub(M, inD), N, 'Entry 171 sight lines with a DESI quasar', log=log)
out['entry171_not_in_desi'] = a.analyse(a.sub(M, ~inD), N, 'Entry 171 sight lines without a DESI quasar', log=log)
U = {key: np.r_[M[key], S[key][~inM]] for key in ('z', 'd', 'e', 'b', 'l')}
out['union'] = a.analyse(U, 1000, 'union of Entry 171 and DESI-only', log=log)
# redshift agreement on the common sight lines
zm = dict(zip(M['src'].tolist(), M['z'].tolist()))
dz = np.array([zm[s] - z for s, z in zip(S['src'][inM].tolist(), S['z'][inM].tolist())])
out['common_redshift'] = dict(n=int(inM.sum()), agree_within_0p05=float((np.abs(dz) < 0.05).mean()), median_abs=float(np.median(np.abs(dz))))
out['tsnr2_cut'] = cut
out['counts'] = dict(D=int(len(S['z'])), also_in_entry171=int(inM.sum()), desi_only=int((~inM).sum()),
                     entry171=int(len(M['z'])), entry171_in_desi=int(inD.sum()), union=int(len(U['z'])),
                     z_lt_0p4=int((S['z'] < 0.4).sum()), tsnr2_lt_cut=int((S['tsnr2'] < cut).sum()))
log(json.dumps(out['counts']))
json.dump(out, open(os.path.join(HERE, 'explore.json'), 'w'), indent=1, default=float)
