"""Entry 187, after the frozen result (not pre-registered): the radial CMB temperature profile around the galaxies.

For each galaxy used by the frozen stack (stack.py: at least half the pixels of the 0.2 Mpc disc and of the 2-3 Mpc annulus
unmasked), the mean of the used pixels of the map in bins of projected distance 0-0.2, 0.2-0.5, 0.5-1, 1-1.5, 1.5-2 and
2-3 Mpc (H0 70, Om 0.3, as the rule).  Equal-weight mean over galaxies per bin (a galaxy enters a bin if any of its pixels
there is used).  Null: the same 1,000 rotations as stack.py (seed 20261187).  The first bin is the claim's statistic D and
first minus last is C.
Prediction: a profile of depth -15 uK (top hat inside 0.2 Mpc; or 'fade', falling linearly to 0 at 1.5 Mpc, as predict.py),
smoothed by the SMICA beam and pixel (Gaussian, sigma as predict.json) with a Hankel transform, averaged in each bin and over
the galaxies' redshifts.
Writes radial_<sample>_<map>.json, and each rotation's profile to work/radial_rot_<sample>_<map>.npy.  Run from the repository root:
  .venv/bin/python -I entries/187-cmb-spirals/radial.py P smica_l5"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, json, time
import numpy as np, healpy as hp
from scipy.special import i0e
from multiprocessing import get_context
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stack as S                                       # DA, loadmap, randrot, NS, SEED, NROT, PRED
EDGES = np.array([0.0, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0])   # Mpc
NB = len(EDGES) - 1
G = {}
def _init(mapname, V, Z):
    G['T'], G['U'] = S.loadmap(mapname)
    G['V'], G['Z'] = V, Z; G['mpc'] = np.array([1.0 / S.DA(z) for z in Z])    # radians per Mpc
def prof(R=None):
    T, U, V = G['T'], G['U'], G['V']
    if R is not None: V = V @ R.T
    acc, cnt = np.zeros(NB), np.zeros(NB); n = 0
    for k, v in enumerate(V):
        m = G['mpc'][k]
        p = hp.query_disc(S.NS, v, 3.0 * m)
        x, y, z = hp.pix2vec(S.NS, p); r = np.arccos(np.clip(x * v[0] + y * v[1] + z * v[2], -1, 1)) / m
        b = np.searchsorted(EDGES, r, side='right') - 1; b[b >= NB] = NB - 1
        u = U[p]
        nall = np.bincount(b, minlength=NB); nuse = np.bincount(b[u], minlength=NB)
        if nuse[0] < 0.5 * nall[0] or nuse[0] == 0: continue
        if nuse[-1] < 0.5 * nall[-1]: continue
        s = np.bincount(b[u], weights=T[p[u]], minlength=NB)
        ok = nuse > 0; acc[ok] += s[ok] / nuse[ok]; cnt[ok] += 1; n += 1
    return acc / np.maximum(cnt, 1), n
def predicted(Z, depth=-15.0):
    """Beam-smoothed bin means of the top-hat and fade profiles, averaged over the redshifts Z."""
    sig = S.PRED['sigma_arcmin']; out = {}
    zs, w = np.unique(np.round(Z / 0.002) * 0.002, return_counts=True)   # the predict.json grid
    for name in ('top', 'fade'):
        tot = np.zeros(NB)
        for zz, ww in zip(zs, w):
            mpc = 1 / S.DA(zz) * 180 / np.pi * 60                       # arcmin per Mpc
            r = np.linspace(0, 3.0 * mpc + 6 * sig, 1500); dr = r[1] - r[0]
            td = 0.2 * mpc
            p0 = (r < td).astype(float) if name == 'top' else np.where(r < td, 1.0, np.clip((1.5 * mpc - r) / (1.3 * mpc), 0, 1))
            # Gaussian smoothing of a radial profile: s(r) = int p(r') exp(-(r-r')^2/2s^2) i0e(r r'/s^2) r'/s^2 dr'
            RR, RP = np.meshgrid(r, r, indexing='ij')
            ker = np.exp(-(RR - RP) ** 2 / (2 * sig ** 2)) * i0e(RR * RP / sig ** 2) * RP / sig ** 2
            s = ker @ p0 * dr
            rb = r / mpc; bm = np.zeros(NB)
            for j in range(NB):
                sel = (rb >= EDGES[j]) & (rb < EDGES[j + 1]); bm[j] = np.sum(s[sel] * r[sel]) / np.sum(r[sel])
            tot += ww * bm
        out[name] = list(depth * tot / w.sum())
    return out
if __name__ == '__main__':
    sample, mapname = sys.argv[1], sys.argv[2]
    A = np.loadtxt(HERE / 'samples' / f'{sample}.txt', ndmin=2); V, Z = A[:, :3], A[:, 3]
    t0 = time.time(); _init(mapname, V, Z); obs, n = prof(None)
    rng = np.random.default_rng(S.SEED); Rs = [S.randrot(rng) for _ in range(S.NROT)]
    with get_context('fork').Pool(S.NPROC) as pool: out = pool.map(prof, Rs, chunksize=4)
    P = np.array([o[0] for o in out]); m, sd = P.mean(0), P.std(0, ddof=1)
    np.save(HERE / 'work' / f'radial_rot_{sample}_{mapname}.npy', P)      # each rotation's profile (git-ignored), for radial_fit.py
    C0 = obs[0] - obs[-1]; Cn = P[:, 0] - P[:, -1]
    r = dict(sample=sample, map=mapname, n_in=int(len(Z)), n_used=int(n), edges_mpc=list(EDGES), profile=list(obs),
             null_mean=list(m), null_sd=list(sd), excess=list(obs - m), t=list((obs - m) / sd),
             corr=np.corrcoef(P.T).round(3).tolist(), C=float(C0), null_C=[float(Cn.mean()), float(Cn.std(ddof=1))],
             predicted_minus15=predicted(Z), nrot=S.NROT, seconds=round(time.time() - t0))
    json.dump(r, open(HERE / f'radial_{sample}_{mapname}.json', 'w'), indent=1)
    print(sample, mapname, 'used', n, f'{time.time() - t0:.0f}s')
    for j in range(NB):
        print(f'  {EDGES[j]:.1f}-{EDGES[j + 1]:.1f} Mpc: {obs[j]:+6.2f}  null {m[j]:+6.2f} +- {sd[j]:5.2f}  t {(obs[j] - m[j]) / sd[j]:+5.2f}'
              f'   pred -15 top {r["predicted_minus15"]["top"][j]:+6.2f} fade {r["predicted_minus15"]["fade"][j]:+6.2f}')
    print(f'  C (bin 1 - bin 6) {C0:+.2f}, null {Cn.mean():+.2f} +- {Cn.std(ddof=1):.2f}')
