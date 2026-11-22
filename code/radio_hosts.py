"""Entry 88: are double radio galaxies young? Peratt's sequence against the stars of their hosts.

Peratt (1986, IEEE Trans. Plasma Sci. 14, 763) places synchrotron-emitting double radio
galaxies and quasars in the first 1e8-1e9 years of interacting galactic Birkeland currents,
with elliptical, peculiar and spiral galaxies forming over the next 1e9-5e9 years.
If so, the stars around a double radio source should be younger than about 1 Gyr, and
younger than the spirals that come later. Dn(4000) of the host measures light-weighted
stellar age without any distance law: ~1.3 for r-weighted ages of 1-3 Gyr, ~1.85 for old
ellipticals of ~10 Gyr (Kauffmann et al. 2003).
Data: FRICAT (219 FR I) and FRIICAT (122 FR II), Capetti et al. 2017, SDSS DR7 spectra, z < 0.15;
comparison: all SDSS DR18 galaxies with MPA-JHU Dn(4000) at 0.03 < z < 0.15, binned by
stellar velocity dispersion and by shape (concentration; Galaxy Zoo spiral/elliptical).
"""
import json, numpy as np, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
rng = np.random.default_rng(88)
def load(fn, cls=False):
    out = []
    for r in open(fn).read().split(';'):
        f = r.split(',')
        d = dict(name=f[0], z=float(f[1]), dn=float(f[2]), sig=float(f[3]), cr=float(f[4]), lr=float(f[5]))
        if cls: d['cls'] = f[6]
        if d['dn'] > 0: out.append(d)
    return out
fr1 = load('fricat.txt'); fr2 = load('friicat.txt', True)
def hist(fn):
    H = {}
    for r in open(fn).read().split(';'):
        sb, db, cc, n = map(int, r.split(','))
        H[(sb, db, cc)] = H.get((sb, db, cc), 0) + n
    return H
Hc = hist('sdss_hist_conc.txt'); Hz = hist('sdss_hist_zoo.txt')
DB = np.arange(16, 53)                     # Dn bins of 0.05: 0.80 ... 2.65
centres = (DB + 0.5)*0.05
def dist(H, sb, cc=None):
    return np.array([sum(H.get((sb, d, c), 0) for c in ((0, 1, 2) if cc is None else (cc,))) for d in DB], float)
def stats_from_hist(h):
    c = np.cumsum(h)/h.sum(); med = np.interp(0.5, c, (DB+1)*0.05)
    return dict(n=int(h.sum()), median=round(float(med), 3), f_old=round(float(h[centres >= 1.85].sum()/h.sum()), 3),
                f_young=round(float(h[centres < 1.40].sum()/h.sum()), 3))
def sample_stats(s):
    dn = np.array([x['dn'] for x in s])
    b = [np.median(rng.choice(dn, len(dn))) for _ in range(3000)]
    return dict(n=len(s), median=round(float(np.median(dn)), 3), median_err=round(float(np.std(b)), 3),
                f_old=round(float(np.mean(dn >= 1.85)), 3), f_young=round(float(np.mean(dn < 1.40)), 3),
                sigma_median=float(np.median([x['sig'] for x in s])), z_median=float(np.median([x['z'] for x in s])))
res = {'FRI': sample_stats(fr1), 'FRII': sample_stats(fr2),
       'FRII_LEG': sample_stats([x for x in fr2 if x.get('cls') == 'LEG']), 'FRII_HEG': sample_stats([x for x in fr2 if x.get('cls') == 'HEG'])}
# comparison at the same velocity dispersion: weight each SDSS sigma bin by the radio hosts' sigma distribution
allr = fr1 + fr2
sbins = [3, 4, 5, 6]                          # 150-200, 200-250, 250-300, 300-350 km/s (350-400 merged into 6)
def sbin(s): return min(int(s//50), 6)
w = np.array([sum(1 for x in allr if 150 <= x['sig'] < 400 and sbin(x['sig']) == b) for b in sbins], float)
res['radio_sigma_bins'] = dict(zip(['150-200', '200-250', '250-300', '300-400'], w.astype(int).tolist()))
def matched(H, cc):
    tot = np.zeros(len(DB))
    for wi, b in zip(w, sbins):
        h = dist(H, b, cc) + (dist(H, 7, cc) if b == 6 else 0)
        if h.sum() > 0: tot += wi*h/h.sum()
    return tot
comp = {'all': matched(Hc, None), 'late_conc': matched(Hc, 0), 'early_conc': matched(Hc, 2),
        'zoo_spiral': matched(Hz, 0), 'zoo_elliptical': matched(Hz, 2)}
res['sdss_matched'] = {k: stats_from_hist(v) for k, v in comp.items()}
res['sdss_counts'] = {k: int(sum(v for (sb, db, cc), v in H.items() if (k[1] is None or cc == k[1]) and 3 <= sb <= 7))
                      for k, H in (( ('c', None), Hc), (('c0', 0), Hc), (('z0', 0), Hz), (('z2', 2), Hz))} if False else None
nz = {lab: int(sum(v for (sb, db, cc), v in Hz.items() if cc == c)) for lab, c in (('zoo_spiral', 0), ('zoo_elliptical', 2))}
nz['all'] = int(sum(Hc.values())); res['sdss_numbers'] = nz
# percentile of each radio host within all SDSS galaxies of its sigma bin
pct = []
for x in allr:
    if not (150 <= x['sig'] < 400): continue
    b = sbin(x['sig']); h = dist(Hc, b) + (dist(Hc, 7) if b == 6 else 0)
    c = np.cumsum(h)/h.sum(); edges = (DB+1)*0.05
    pct.append(float(np.interp(x['dn'], np.concatenate([[DB[0]*0.05], edges]), np.concatenate([[0], c]))))
pct = np.array(pct)
res['percentile_in_sdss_same_sigma'] = dict(n=len(pct), mean=round(float(pct.mean()), 3), err=round(float(pct.std()/np.sqrt(len(pct))), 3),
                                           above_median=round(float(np.mean(pct > 0.5)), 3))
# Peratt: hosts younger than ~1 Gyr -> Dn(4000) below ~1.4 (r-weighted 1-3 Gyr sits at ~1.3)
res['peratt_hosts_young'] = dict(n_below_1p4=int(sum(1 for x in allr if x['dn'] < 1.4)), n_total=len(allr))
# radio luminosity: do the most powerful sources have younger hosts?
lr = np.array([x['lr'] for x in allr]); dn = np.array([x['dn'] for x in allr])
from scipy import stats
res['dn_vs_radio_power_spearman'] = [round(float(v), 3) for v in stats.spearmanr(lr, dn)]
# histograms for the figure, normalised
def nh(vals):
    h = np.array([np.sum((vals >= d*0.05) & (vals < (d+1)*0.05)) for d in DB], float); return (h/h.sum()).round(4).tolist()
res['fig'] = dict(bins=[round(float(d*0.05), 2) for d in DB], FRI=nh(np.array([x['dn'] for x in fr1])), FRII=nh(np.array([x['dn'] for x in fr2])),
                  spiral=(comp['zoo_spiral']/comp['zoo_spiral'].sum()).round(4).tolist(),
                  elliptical=(comp['zoo_elliptical']/comp['zoo_elliptical'].sum()).round(4).tolist(),
                  all=(comp['all']/comp['all'].sum()).round(4).tolist())
json.dump(res, open('radio_hosts.json', 'w'), indent=1)
for k, v in res.items():
    if k != 'fig': print(k, v)
