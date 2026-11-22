"""flash_compare.py - Workbench Entry 95
Large-angle E-mode polarization from the plasma code (cmbreion.py) for three histories with the same scaling as Entry 94
(first temperature peak = 5,733 uK^2): tau = 0.054 in one late step (Planck), tau = 0.09 in one late step, and tau = 0.09 made
of a late step at z = 6.95 (what the quasar and galaxy spectra allow) plus an early burst between z = 15 and 30.
Reads cmbreion_tau054.json, cmbreion_tau09.json and cmbreion_flash09.json; writes flash_compare.json."""
import json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
runs = {k: json.load(open(os.path.join(HERE, f'cmbreion_{k}.json'))) for k in ('tau054', 'tau09', 'flash09')}
ells = np.array(runs['tau054']['ells']); low = ells <= 30
out = {'ells_low': [int(l) for l in ells[low]], 'flash': runs['flash09'].get('flash'), 'z_re_flash': runs['flash09']['z_re']}
for k, r in runs.items():
    EE = np.array(r['runs']['0.965']['EE'])/16*5733.0          # Ma-Bertschinger E normalisation, then uK^2 on the first temperature peak
    e = EE[low]; L = ells[low]
    out[k] = dict(EE_low=[round(float(x), 5) for x in e], mean_2_10=round(float(e[(L >= 2) & (L <= 10)].mean()), 4),
                  mean_10_20=round(float(e[(L >= 10) & (L <= 20)].mean()), 4), mean_20_30=round(float(e[(L >= 20) & (L <= 30)].mean()), 4),
                  peak=round(float(e.max()), 4), l_peak=int(L[np.argmax(e)]))
for band in ('mean_2_10', 'mean_10_20', 'mean_20_30'):
    out['flash_over_054_' + band] = round(out['flash09'][band]/out['tau054'][band], 2)
    out['late09_over_054_' + band] = round(out['tau09'][band]/out['tau054'][band], 2)
json.dump(out, open(os.path.join(HERE, 'flash_compare.json'), 'w'), indent=1)
for k in ('tau054', 'tau09', 'flash09'): print(k, {kk: v for kk, v in out[k].items() if kk != 'EE_low'})
print({k: v for k, v in out.items() if k.startswith(('flash_over', 'late09'))})
