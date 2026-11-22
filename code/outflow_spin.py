"""outflow_spin.py - Workbench Entry 133
Do stars spin along the filaments they form in, as a current pinch would make them?
Carrying Peratt's geometry (discs form across the current, spins along it; Entry 73) down to stars predicts protostellar
outflows, which leave along the spin axis, running parallel to their parent filaments.

Data: projected outflow-filament angles gamma (0-90 deg) from four published regions.
  Perseus     Stephens et al. 2017, ApJ 846, 16 (MASSES/SMA + Herschel FilFinder), 57 outflows, Table 2 gamma_F
  Orion A     Feddersen et al. 2020, ApJ 896, 11 (CARMA-NRO), 45 outflows, Table 3 per lobe; lobes of one outflow averaged
  G28.37+0.07 Kong et al. 2019, ApJ 874, 104 (ALMA), 58 outflows, counts in 0-30 / 30-60 / 60-90 deg only
  Serpens Main Green et al. 2024 (JWST), class-A outflows, PA against the filament axis at 139 deg
Model (as in Stephens et al.): filament direction isotropic; outflow at 3D angle theta from it, drawn uniformly in solid angle from
  parallel 0-20 deg, perpendicular 70-90 deg, or random 0-90 deg; both projected on the sky; gamma = projected angle.
  Mixture f_par, f_perp, 1 - f_par - f_perp; flat prior on the simplex; likelihood from 10-degree bins (30-degree for G28).
"""
import json, math
import numpy as np

PER = [76, 3, 87, 86, 12, 50, 5, 4, 27, 82, 87, 50, 40, 77, 67, 89, 50, 59, 58, 28, 57, 79, 31, 43, 32, 69, 21, 35, 55, 5, 89, 2, 89, 43,
       70, 7, 57, 85, 87, 63, 59, 7, 16, 33, 50, 89, 11, 32, 61, 72, 24, 6, 38, 5, 77, 68, 32]
ORI_LOBES = {'SMZ 11': [27], 'SMZ 17': [10, 38], 'SMZ 21': [43], 'SMZ 30': [35], 'SMZ 50': [36], 'HOPS 10': [71], 'HOPS 11': [62, 86],
             'HOPS 12': [54], 'HOPS 44': [70], 'HOPS 50': [16], 'HOPS 56': [33], 'HOPS 58': [88], 'HOPS 59': [38], 'HOPS 60': [78, 83],
             'HOPS 68': [20, 24], 'HOPS 70': [74], 'HOPS 71': [59], 'HOPS 75': [17], 'HOPS 78': [59, 43], 'HOPS 81': [40, 32],
             'HOPS 84': [58, 49], 'HOPS 87': [27, 11], 'HOPS 88': [36, 35], 'HOPS 92': [43, 41], 'HOPS 96': [55, 71], 'HOPS 99': [48],
             'HOPS 157': [84], 'HOPS 158': [50, 50], 'HOPS 160': [22], 'HOPS 166': [59, 40], 'HOPS 168': [83, 69], 'HOPS 169': [49, 35],
             'HOPS 174': [85], 'HOPS 177': [66], 'HOPS 178': [76], 'HOPS 179': [60, 28], 'HOPS 181': [48, 12], 'HOPS 182': [33, 16],
             'HOPS 192': [53, 34], 'HOPS 198': [0], 'HOPS 203': [85], 'HOPS 355': [75, 84], 'HOPS 368': [61, 84], 'HOPS 370': [37, 48],
             'HOPS 383': [88]}
ORI = [float(np.mean(v)) for v in ORI_LOBES.values()]
G28_COUNTS = [7, 19, 32]                                   # 0-30, 30-60, 60-90
SER_PA = {1: 141.2, 3: 112.0, 4: 108.3, 5: 115.6, 6: 151.6, 7: 158.8, 8: 135.2, 9: 135.6, 10: 138.9, 11: 83.4, 12: 132.7, 14: 68.1,
          16: 160.8, 17: 2.7, 20: 216.8}                   # class A outflows (disk shadows left out)
SER_PA_B = {2: 129.6, 13: 123.2, 15: 156.3, 18: 76.1}
def fold(d): d = abs(d) % 180; return min(d, 180 - d)
SER = [fold(p - 139) for p in SER_PA.values()]
SER_AB = SER + [fold(p - 139) for p in SER_PA_B.values()]
assert len(PER) == 57 and len(ORI_LOBES) == 45 and sum(len(v) for v in ORI_LOBES.values()) == 67 and sum(G28_COUNTS) == 58

rng = np.random.default_rng(130)
def sample_gamma(lo, hi, n=2_000_000):
    F = rng.normal(size=(n, 3)); F /= np.linalg.norm(F, axis=1)[:, None]
    ct = rng.uniform(math.cos(math.radians(hi)), math.cos(math.radians(lo)), n); st = np.sqrt(1 - ct**2); ph = rng.uniform(0, 2*np.pi, n)
    a = np.cross(F, [0, 0, 1]); bad = np.linalg.norm(a, axis=1) < 1e-6; a[bad] = np.cross(F[bad], [0, 1, 0])
    a /= np.linalg.norm(a, axis=1)[:, None]; b = np.cross(F, a)
    O = F*ct[:, None] + a*(st*np.cos(ph))[:, None] + b*(st*np.sin(ph))[:, None]
    g = np.degrees(np.abs(np.arctan2(F[:, 1], F[:, 0]) - np.arctan2(O[:, 1], O[:, 0]))) % 180
    return np.minimum(g, 180 - g)
G = {k: sample_gamma(*r) for k, r in {'par': (0, 20), 'perp': (70, 90), 'rand': (0, 90)}.items()}
E10 = np.arange(0, 91, 10); E30 = np.array([0, 30, 60, 90])
H10 = {k: np.histogram(v, E10)[0]/len(v) for k, v in G.items()}
H30 = {k: np.histogram(v, E30)[0]/len(v) for k, v in G.items()}

def counts10(x): return np.histogram(np.clip(x, 0, 89.999), E10)[0]
DATA = {'Perseus': ('10', counts10(PER)), 'Orion A': ('10', counts10(ORI)), 'G28.37+0.07': ('30', np.array(G28_COUNTS)),
        'Serpens Main': ('10', counts10(SER))}
fp = np.linspace(0, 1, 201)
def loglike(kind, c):
    H = H10 if kind == '10' else H30
    L = np.full((len(fp), len(fp)), -np.inf)
    for i, a in enumerate(fp):
        for j, b in enumerate(fp):
            if a + b > 1 + 1e-9: continue
            p = a*H['par'] + b*H['perp'] + (1 - a - b)*H['rand']
            L[i, j] = float(np.sum(c*np.log(np.maximum(p, 1e-300))))
    return L
def summarize(L):
    P = np.exp(L - L.max()); P /= P.sum()
    mpar = P.sum(axis=1); mperp = P.sum(axis=0)
    cpar = np.cumsum(mpar); cperp = np.cumsum(mperp)
    q = lambda c, x: float(fp[np.searchsorted(c, x)])
    i, j = np.unravel_index(np.argmax(L), L.shape)
    return dict(par_med=q(cpar, 0.5), par_95=q(cpar, 0.95), P_par_gt_half=float(mpar[fp > 0.5].sum()),
                perp_med=q(cperp, 0.5), perp_lo95=q(cperp, 0.05), best=(float(fp[i]), float(fp[j])),
                dlnL_pure_par=float(L[200, 0] - L.max()), dlnL_random=float(L[0, 0] - L.max()), dlnL_pure_perp=float(L[0, 200] - L.max()))
out = {'n': {'Perseus': 57, 'Orion A': 45, 'G28.37+0.07': 58, 'Serpens Main': len(SER)}, 'regions': {}, 'hist': {}}
LL = {}
for name, (kind, c) in DATA.items():
    LL[name] = loglike(kind, c); out['regions'][name] = summarize(LL[name]); out['hist'][name] = c.tolist()
out['regions']['Serpens Main (A+B)'] = summarize(loglike('10', counts10(SER_AB)))
out['regions']['Three regions without Serpens'] = summarize(LL['Perseus'] + LL['Orion A'] + LL['G28.37+0.07'])
out['regions']['All four'] = summarize(sum(LL.values()))
out['model_hist10'] = {k: [round(float(x), 4) for x in v] for k, v in H10.items()}
out['model_hist30'] = {k: [round(float(x), 4) for x in v] for k, v in H30.items()}
out['serpens_gamma'] = [round(x, 1) for x in SER]
out['frac_gamma_lt30'] = {'Perseus': float(np.mean(np.array(PER) < 30)), 'Orion A': float(np.mean(np.array(ORI) < 30)),
                          'G28.37+0.07': 7/58, 'Serpens Main': float(np.mean(np.array(SER) < 30)),
                          'model_par': float(H30['par'][0]), 'model_rand': float(H30['rand'][0]), 'model_perp': float(H30['perp'][0])}
json.dump(out, open('outflow_spin.json', 'w'), indent=1)
for k, v in out['regions'].items(): print(f'{k:32s}', {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items()})
print(out['frac_gamma_lt30']); print(out['hist']); print(out['model_hist10'])
