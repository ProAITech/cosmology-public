"""Entry 187, before the rule: the power of a disc stack around SDSS-like spirals on the Planck SMICA map.
Uses no galaxy: random positions in an SDSS-like patch (RA 120-250, Dec 0-60), redshifts drawn as z^2 in 0.02-0.05,
and random rigid rotations of that set as the null. Statistics: D = mean T inside 0.2 Mpc (H0 = 70, flat, Om = 0.3),
C = D minus the mean T at 2-3 Mpc. SMICA I_STOKES (K_CMB) with l <= 5 removed (low-l alm from I_STOKES_INP).
Run from the repository root: .venv/bin/python -I entries/187-cmb-spirals/power.py"""
import numpy as np, healpy as hp, time, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
SM = ROOT / 'data/planck-pr3-smica/COM_CMB_IQU-smica_2048_R3.00_full.fits'
MK = ROOT / 'data/planck-pr3-masks/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits'
NS = 2048
def DA(z, H0=70.0, Om=0.3):
    zz = np.linspace(0, z, 200); E = np.sqrt(Om * (1 + zz) ** 3 + 1 - Om)
    return 299792.458 / H0 * np.trapezoid(1 / E, zz) / (1 + z)
T = hp.read_map(SM, field=0) * 1e6; Ti = hp.read_map(SM, field=1) * 1e6
m = hp.read_map(MK, field=0) > 0.5
alm = hp.map2alm(Ti, lmax=5, iter=3); low = hp.alm2map(alm, NS, lmax=5)
T = T - low
print('map ready; used fraction', m.mean(), flush=True)
rng = np.random.default_rng(187)
def sample(n):
    ra = rng.uniform(120, 250, 4 * n); s = rng.uniform(np.sin(0), np.sin(np.radians(60)), 4 * n); dec = np.degrees(np.arcsin(s))
    z = (0.02 ** 3 + rng.uniform(0, 1, n) * (0.05 ** 3 - 0.02 ** 3)) ** (1 / 3)
    return hp.ang2vec(ra[:n], dec[:n], lonlat=True), z
def stat(V, th1, th2, th3):
    D, C = [], []
    for v, a, b, c in zip(V, th1, th2, th3):
        p = hp.query_disc(NS, v, a); p = p[m[p]]
        if len(p) < 3: continue
        q = hp.query_disc(NS, v, c); q2 = hp.query_disc(NS, v, b); q = np.setdiff1d(q, q2); q = q[m[q]]
        if len(q) < 10: continue
        d = T[p].mean(); D.append(d); C.append(d - T[q].mean())
    return np.mean(D), np.mean(C), len(D)
for n in (1400, 4000):
    V, z = sample(n); da = np.array([DA(x) for x in z])
    th1, th2, th3 = 0.2 / da, 2.0 / da, 3.0 / da
    print(f'N={n}: disc radius median {np.degrees(np.median(th1))*60:.1f} arcmin', flush=True)
    t0 = time.time(); out = []
    for k in range(int(sys.argv[1]) if len(sys.argv) > 1 else 60):
        R = hp.Rotator(rot=[rng.uniform(0, 360), np.degrees(np.arcsin(rng.uniform(-1, 1))), rng.uniform(0, 360)], deg=True)
        Vr = np.array(R(V.T)).T
        out.append(stat(Vr, th1, th2, th3))
    out = np.array(out)
    print(f'  rotations {len(out)} in {time.time()-t0:.0f}s: D mean {out[:,0].mean():+.2f} sd {out[:,0].std(ddof=1):.2f} | C mean {out[:,1].mean():+.2f} sd {out[:,1].std(ddof=1):.2f} | kept {out[:,2].mean():.0f}', flush=True)
