"""Entry 187: is the CMB colder around late-type spirals beyond redshift 0.02?  Rule: Workbench section 187.

Per galaxy at redshift z: theta_d = 0.2 Mpc / D_A(z) (H0 70, Om 0.3; 0.25 deg at z 0.01, as Hansen et al. 2025).
  D = mean of the used pixels whose centres lie within theta_d; A = mean of the used pixels at 2-3 Mpc; C = D - A.
  A galaxy is used if at least half of the pixels of each region are inside the Planck common mask's used area.
Stack: equal-weight means of D and C over used galaxies.  Null: NROT uniform random rotations of the whole sample.
Prediction (predict.json): a top-hat cooling of depth d inside 0.2 Mpc, smoothed by the 5' beam, gives C = d * f, f = f_disc - f_ann.
Depth estimate delta = (C - null mean) / mean f; its error sd_null / mean f.
Maps: work/<map>_l<L>.i16 from prep_maps.py (micro-K x 10, -32768 masked, RING, Nside 2048).
Samples: samples/<name>.txt, one galaxy per line: Galactic unit vector x y z, redshift (decoded from the browser read-out, decode.py).
Run from the repository root:
  .venv/bin/python -I entries/187-cmb-spirals/stack.py synthetic          # injection test on random positions, no galaxy data
  .venv/bin/python -I entries/187-cmb-spirals/stack.py run <sample> <map>  # the real stack"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, json, time
import numpy as np, healpy as hp
from pathlib import Path
from multiprocessing import get_context
HERE = Path(__file__).resolve().parent; W = HERE / 'work'; NS = 2048
NROT, SEED, NPROC = 1000, 20261187, 8
PRED = json.load(open(HERE / 'predict.json'))
def DA(z, H0=70.0, Om=0.3):
    zz = np.linspace(0, z, 400); E = np.sqrt(Om * (1 + zz) ** 3 + 1 - Om)
    return 299792.458 / H0 * np.trapezoid(1 / E, zz) / (1 + z)
def fpred(z):
    i = int(np.clip(round((z - PRED['z'][0]) / (PRED['z'][1] - PRED['z'][0])), 0, len(PRED['z']) - 1))
    return PRED['top']['f_disc'][i] - PRED['top']['f_ann'][i], PRED['fade']['f_disc'][i] - PRED['fade']['f_ann'][i]
def loadmap(name):
    q = np.fromfile(W / f'{name}.i16', dtype='<i2'); used = q != -32768
    return q.astype(np.float32) / 10.0, used
G = {}
def _init(mapname, V, Z):
    G['T'], G['U'] = loadmap(mapname) if isinstance(mapname, str) else mapname
    G['V'], G['Z'] = V, Z; da = np.array([DA(z) for z in Z]); G['td'], G['t2'], G['t3'] = 0.2 / da, 2.0 / da, 3.0 / da
def stack(R=None):
    T, U, V = G['T'], G['U'], G['V']
    if R is not None: V = V @ R.T
    D, C, idx = [], [], []
    for k, v in enumerate(V):
        p = hp.query_disc(NS, v, G['td'][k]); pu = p[U[p]]
        if len(pu) < 0.5 * len(p) or len(pu) == 0: continue
        a = np.setdiff1d(hp.query_disc(NS, v, G['t3'][k]), hp.query_disc(NS, v, G['t2'][k]), assume_unique=True); au = a[U[a]]
        if len(au) < 0.5 * len(a): continue
        d = T[pu].mean(); D.append(d); C.append(d - T[au].mean()); idx.append(k)
    return np.mean(D), np.mean(C), np.array(idx)
def randrot(rng):
    u1, u2, u3 = rng.random(3); a, b = np.sqrt(1 - u1), np.sqrt(u1)
    x, y, z, w = a * np.sin(2 * np.pi * u2), a * np.cos(2 * np.pi * u2), b * np.sin(2 * np.pi * u3), b * np.cos(2 * np.pi * u3)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
def _rot(R): d, c, i = stack(R); return d, c, len(i)
def analyse(mapname, V, Z, nrot=NROT, seed=SEED):
    _init(mapname, V, Z); d0, c0, idx = stack(None)
    rng = np.random.default_rng(seed); Rs = [randrot(rng) for _ in range(nrot)]
    with get_context('fork').Pool(NPROC) as pool: out = np.array(pool.map(_rot, Rs, chunksize=4))
    f = np.array([fpred(Z[i]) for i in idx]); ftop, ffade = f[:, 0].mean(), f[:, 1].mean()
    mD, sD, mC, sC = out[:, 0].mean(), out[:, 0].std(ddof=1), out[:, 1].mean(), out[:, 1].std(ddof=1)
    eC = c0 - mC; delta, sdelta = eC / ftop, sC / ftop
    return dict(n_in=len(Z), n_used=int(len(idx)), D=d0, C=c0, null_D=[mD, sD], null_C=[mC, sC], t_D=(d0 - mD) / sD, t_C=eC / sC,
                p_C=float((out[:, 1] <= c0).mean()), p_D=float((out[:, 0] <= d0).mean()), f_top=ftop, f_fade=ffade,
                delta_top=delta, sd_delta_top=sdelta, delta_fade=eC / ffade, upper95_cooling_top=-(delta - 1.645 * sdelta),
                n_used_rot_mean=float(out[:, 2].mean()), nrot=nrot)
def decide(r, depth=-15.0):
    if r['n_used'] < 500: return 'not tested (fewer than 500 galaxies used)'
    if r['t_C'] <= -3: return 'cooling seen'
    if r['delta_top'] >= depth + 3 * r['sd_delta_top']: return 'the claimed cooling is excluded'
    return 'unresolved'
def synthetic():
    rng = np.random.default_rng(187); n = 2800
    ra = rng.uniform(120, 250, n); dec = np.degrees(np.arcsin(rng.uniform(0, np.sin(np.radians(60)), n)))
    z = (0.02 ** 3 + rng.uniform(0, 1, n) * (0.05 ** 3 - 0.02 ** 3)) ** (1 / 3)
    V = np.array(hp.Rotator(coord=['C', 'G'])(hp.ang2vec(ra, dec, lonlat=True).T)).T
    T, U = loadmap('smica_l5')
    for depth in (0.0, -15.0):
        Tm = T.copy()
        if depth:
            inj = np.zeros(hp.nside2npix(NS))
            for v, zz in zip(V, z): inj[hp.query_disc(NS, v, 0.2 / DA(zz))] += depth
            Tm = (T + hp.smoothing(inj, fwhm=np.radians(5 / 60), lmax=4000)).astype(np.float32)
        t0 = time.time(); r = analyse((Tm, U), V, z, nrot=200)
        print(f'injected {depth:+.0f} uK: used {r["n_used"]}, C {r["C"]:+.2f}, null {r["null_C"][0]:+.2f} +- {r["null_C"][1]:.2f}, t_C {r["t_C"]:+.2f}, '
              f'delta_top {r["delta_top"]:+.2f} +- {r["sd_delta_top"]:.2f}, D t {r["t_D"]:+.2f} | {decide(r)} | {time.time()-t0:.0f}s', flush=True)
if __name__ == '__main__':
    if sys.argv[1] == 'synthetic': synthetic()
    elif sys.argv[1] == 'run':
        S = np.loadtxt(HERE / 'samples' / f'{sys.argv[2]}.txt', ndmin=2); r = analyse(sys.argv[3], S[:, :3], S[:, 3])
        r['decision_minus15'] = decide(r); r['sample'], r['map'] = sys.argv[2], sys.argv[3]
        r = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else [float(x) for x in v] if isinstance(v, list) else v) for k, v in r.items()}
        print(json.dumps(r)); json.dump(r, open(HERE / f'result_{sys.argv[2]}_{sys.argv[3]}.json', 'w'), indent=1)
