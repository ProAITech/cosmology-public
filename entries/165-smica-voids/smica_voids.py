"""Entry 165: the deeper-slab void test frozen in Entry 148, run on the Planck PR3 SMICA map.

Protocol: Entry 148 (workbench.html), fixed 11 November 2026 before any CMB map was opened.
Choices the protocol left open: README.md in this folder, fixed 13 November 2026 before any void was evaluated.
Data (URL, size, sha256 in README.md):
  data/planck-pr3-smica/COM_CMB_IQU-smica_2048_R3.00_full.fits        table COMP-MAP, columns I_STOKES, I_STOKES_INP (K_CMB)
  data/planck-pr3-masks/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits  table, column TMASK (1 = use)
Catalogue: code/voids_slab1.txt, code/voids_slab2.txt (frozen in Entry 148).
Spectrum for the simulated skies: code/cmbreion_tau054.json (the notebook's TT, as code/prereg.py), and SMICA's own.

Statistic, per void of radius theta: mean T within theta/2 minus mean T in the annulus theta..1.5 theta, used pixels only;
the stack is the equal-weight mean over the sample.

Run from the repository root:  .venv/bin/python entries/165-smica-voids/smica_voids.py [prep|null|real|all]
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, json, time
from pathlib import Path
from multiprocessing import get_context
import numpy as np
import healpy as hp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA, CODE, WORK = ROOT / 'data', ROOT / 'code', HERE / 'work'
SMICA = DATA / 'planck-pr3-smica' / 'COM_CMB_IQU-smica_2048_R3.00_full.fits'
MASK = DATA / 'planck-pr3-masks' / 'COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits'
NS, NNAT, LMAX = 256, 2048, 700
NSIM, NROT, SEED = 1000, 1000, 20261113
UNSEEN = hp.UNSEEN
NPZ = WORK / 'lowres.npz'


# ---------------------------------------------------------------- catalogue and samples
def _load(f):
    return np.array([list(map(float, r.split(','))) for r in open(f).read().strip().split(';')])


def catalogue():
    S1, S2 = _load(CODE / 'voids_slab1.txt'), _load(CODE / 'voids_slab2.txt')
    slab = np.r_[np.ones(len(S1), int), 2 * np.ones(len(S2), int)]
    l = np.r_[S1[:, 0], S2[:, 0]]; b = np.r_[S1[:, 1], S2[:, 1]]; th = np.r_[S1[:, 2], S2[:, 2]]
    flag = np.r_[-np.ones(len(S1), int), S2[:, 7].astype(int)]   # slab 2: 0 fresh, 1 overlaps a nearer void
    return dict(slab=slab, l=l, b=b, th=th, flag=flag)


def samples(c):
    s1, s2, th = c['slab'] == 1, c['slab'] == 2, c['th']
    fresh, over = s2 & (c['flag'] == 0), s2 & (c['flag'] == 1)
    q = lambda m: m & (th >= np.percentile(th[m], 75))
    return {'check_slab1_all': s1, 'check_slab1_largest_quarter': q(s1),
            'primary_fresh': fresh, 'primary_fresh_largest_quarter': q(fresh),
            'overlapping': over, 'slab2_all': s2}


# ---------------------------------------------------------------- geometry and statistic
def regions(nside, l, b, th, valid):
    """Used-pixel indices inside theta/2 and in the annulus theta..1.5 theta."""
    v0 = hp.ang2vec(l, b, lonlat=True)
    pix = hp.query_disc(nside, v0, np.radians(1.5 * th), nest=True, inclusive=False)
    v = np.array(hp.pix2vec(nside, pix, nest=True))
    ang = np.degrees(np.arccos(np.clip(v0 @ v, -1, 1)))
    ok = valid[pix]
    return pix[(ang <= th / 2) & ok], pix[(ang > th) & (ang <= 1.5 * th) & ok]


def stat(m, inner, ann):
    return m[inner].mean() - m[ann].mean() if len(inner) and len(ann) else np.nan


def md_subtract(m, valid, nvec):
    """Fit monopole and dipole on the used pixels and subtract them (same call for data and simulations)."""
    mono, dip = hp.fit_dipole(np.where(valid, m, UNSEEN), nest=True, bad=UNSEEN, gal_cut=0)
    return np.where(valid, m - (mono + dip @ nvec), UNSEEN), mono, dip


def stack(T, members):
    v = T[members]; v = v[np.isfinite(v)]
    return (float(v.mean()), int(len(v))) if len(v) else (np.nan, 0)


def cl_notebook():
    d = json.load(open(CODE / 'cmbreion_tau054.json')); run = list(d['runs'].values())[0]
    L0, TT = np.array(d['ells'], float), np.array(run['TT']) * 5733.0   # D_l in uK^2, first peak 5,733 as prereg.py
    ell = np.arange(LMAX + 1); cl = np.zeros(LMAX + 1)
    cl[2:] = np.interp(ell[2:], L0, TT) * 2 * np.pi / (ell[2:] * (ell[2:] + 1))
    return cl


# ---------------------------------------------------------------- stage 1: low-resolution map, mask, observed spectrum
def prep():
    from astropy.io import fits
    WORK.mkdir(exist_ok=True); t0 = time.time()
    with fits.open(MASK, memmap=True) as h: good = np.asarray(h[1].data.field(0), dtype=np.float32) == 1
    with fits.open(SMICA, memmap=True) as h:
        I = np.asarray(h[1].data.field('I_STOKES'), dtype=np.float64) * 1e6
        Ip = np.asarray(h[1].data.field('I_STOKES_INP'), dtype=np.float64) * 1e6
    low = hp.ud_grade(np.where(good, I, UNSEEN), NS, pess=True, order_in='NESTED', order_out='NESTED')
    valid = low > -1e20
    nvec = np.array(hp.pix2vec(NS, np.arange(hp.nside2npix(NS)), nest=True))
    low_c, mono, dip = md_subtract(low, valid, nvec)
    ipl = hp.reorder(hp.ud_grade(Ip, NS, order_in='NESTED', order_out='NESTED'), n2r=True)
    cl_obs = hp.anafast(ipl, lmax=LMAX); cl_obs[:2] = 0
    np.savez(NPZ, low=low_c, valid=valid, mono=mono, dip=dip, cl_obs=cl_obs, cl_nb=cl_notebook())
    print(f'prep: unmasked at Nside 2048 {good.mean():.4f}; used at Nside {NS} {valid.mean():.4f}; '
          f'monopole {mono:.2f} uK, dipole {np.linalg.norm(dip):.1f} uK; rms of used pixels after removal '
          f'{low_c[valid].std():.1f} uK; {time.time() - t0:.0f} s')


# ---------------------------------------------------------------- workers (spawned processes)
G = {}


def init():
    z = np.load(NPZ); c = catalogue()
    G.update(low=z['low'], valid=z['valid'], cl={'nb': z['cl_nb'], 'obs': z['cl_obs']}, c=c,
             nvec=np.array(hp.pix2vec(NS, np.arange(hp.nside2npix(NS)), nest=True)))
    G['idx'] = [regions(NS, c['l'][i], c['b'][i], c['th'][i], G['valid']) for i in range(len(c['th']))]


def sim_task(a):
    kind, i = a
    np.random.seed(SEED + (100000 if kind == 'nb' else 200000) + i)
    m = hp.reorder(hp.synfast(G['cl'][kind], NS, lmax=LMAX, pixwin=False, new=True), r2n=True)
    m, _, _ = md_subtract(m, G['valid'], G['nvec'])
    return np.array([stat(m, *ix) for ix in G['idx']]), float(m[G['valid']].std())


def rot_task(i):
    dl = np.random.default_rng(SEED + 700000 + i).uniform(0, 360); c = G['c']; T = []
    for j in range(len(c['th'])):
        T.append(stat(G['low'], *regions(NS, (c['l'][j] + dl) % 360, c['b'][j], c['th'][j], G['valid'])))
    return np.array(T), dl


def null():
    t0 = time.time(); out = {}
    with get_context('spawn').Pool(14, initializer=init) as P:
        for kind in ('nb', 'obs'):
            r = P.map(sim_task, [(kind, i) for i in range(NSIM)], chunksize=8)
            out[f'sim_{kind}'] = np.array([x[0] for x in r]); out[f'rms_{kind}'] = np.array([x[1] for x in r])
            print(f'null: {NSIM} sims ({kind}) done, {time.time() - t0:.0f} s; mean sky rms {out[f"rms_{kind}"].mean():.1f} uK')
        r = P.map(rot_task, range(NROT), chunksize=8)
        out['rot'] = np.array([x[0] for x in r]); out['rot_dl'] = np.array([x[1] for x in r])
        print(f'null: {NROT} rotations done, {time.time() - t0:.0f} s')
    np.savez(WORK / 'nulls.npz', **out)


# ---------------------------------------------------------------- stage 3: the real stack, nulls, decisions
def native_T(c, good, I, mono, dip):
    """Cross-check on the native Nside-2048 map: every unmasked pixel, same monopole and dipole."""
    res = []
    for i in range(len(c['th'])):
        l, b, th = c['l'][i], c['b'][i], c['th'][i]
        v0 = hp.ang2vec(l, b, lonlat=True)
        pix = hp.query_disc(NNAT, v0, np.radians(1.5 * th), nest=True, inclusive=False)
        v = np.array(hp.pix2vec(NNAT, pix, nest=True)); ang = np.degrees(np.arccos(np.clip(v0 @ v, -1, 1)))
        val = I[pix] - (mono + dip @ v); ok = good[pix]
        a, b2 = (ang <= th / 2) & ok, (ang > th) & (ang <= 1.5 * th) & ok
        res.append(val[a].mean() - val[b2].mean() if a.any() and b2.any() else np.nan)
    return np.array(res)


def real():
    from astropy.io import fits
    z = np.load(NPZ); nl = np.load(WORK / 'nulls.npz'); c = catalogue(); S = samples(c)
    low, valid = z['low'], z['valid']
    idx = [regions(NS, c['l'][i], c['b'][i], c['th'][i], valid) for i in range(len(c['th']))]
    T = np.array([stat(low, *ix) for ix in idx])
    n_in = np.array([len(ix[0]) for ix in idx]); n_an = np.array([len(ix[1]) for ix in idx])
    with fits.open(MASK, memmap=True) as h: good = np.asarray(h[1].data.field(0), dtype=np.float32) == 1
    with fits.open(SMICA, memmap=True) as h: I = np.asarray(h[1].data.field('I_STOKES'), dtype=np.float64) * 1e6
    Tn = native_T(c, good, I, float(z['mono']), z['dip'])
    voids = [dict(slab=int(c['slab'][i]), l=float(c['l'][i]), b=float(c['b'][i]), radius=float(c['th'][i]),
                  flag=int(c['flag'][i]), T=None if np.isnan(T[i]) else float(T[i]),
                  T_native=None if np.isnan(Tn[i]) else float(Tn[i]), n_inner=int(n_in[i]), n_annulus=int(n_an[i]))
             for i in range(len(c['th']))]
    out = dict(generated=time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime()), seed=SEED, nside=NS, lmax=LMAX,
               nsim=NSIM, nrot=NROT, sky_rms_real=float(low[valid].std()),
               sky_rms_sim_nb=float(nl['rms_nb'].mean()), sky_rms_sim_obs=float(nl['rms_obs'].mean()),
               voids=voids, samples={})
    for name, m in S.items():
        v, n = stack(T, m); vn, nn = stack(Tn, m)
        r = dict(N=int(m.sum()), N_used=n, stack=v, stack_native=vn, N_used_native=nn, dropped=int(m.sum()) - n)
        for k, key in (('sim_nb', 'sim_nb'), ('sim_obs', 'sim_obs'), ('rot', 'rot')):
            ns_ = np.array([stack(row, m)[0] for row in nl[key]]); ns_ = ns_[np.isfinite(ns_)]
            r[k] = dict(mean=float(ns_.mean()), sigma=float(ns_.std(ddof=1)), n=int(len(ns_)),
                        z=float((v - ns_.mean()) / ns_.std(ddof=1)), p_upper=float((1 + (ns_ >= v).sum()) / (1 + len(ns_))),
                        p_two=float((1 + (np.abs(ns_ - ns_.mean()) >= abs(v - ns_.mean())).sum()) / (1 + len(ns_))),
                        q16_84=[float(x) for x in np.percentile(ns_, [15.87, 84.13])])
            r[k]['null_subtracted'] = float(v - ns_.mean())
            if name in ('primary_fresh', 'check_slab1_largest_quarter', 'overlapping'):
                r[k]['values'] = [round(float(x), 2) for x in ns_]
        out['samples'][name] = r
    chk = out['samples']['check_slab1_largest_quarter']['stack']; pri = out['samples']['primary_fresh']['stack']
    out['decision'] = dict(
        step1_check_value=chk, step1=('the claim depends on the void finder (under +10 uK)' if chk < 10 else
                                      'the nearest-slab warmth is reproduced with these voids (+10 uK or more)'),
        step2_value=pri, step2=('real effect along the path (+20 uK or more)' if pri >= 20 else
                                'chance alignment, or the nearest structures only (below +15 uK)' if pri < 15 else
                                'unresolved (between +15 and +20 uK)'))
    json.dump(out, open(HERE / 'results.json', 'w'))
    print(f"sky rms (uK): real {out['sky_rms_real']:.1f}, simulated notebook {out['sky_rms_sim_nb']:.1f}, "
          f"simulated SMICA {out['sky_rms_sim_obs']:.1f}")
    print(f"{'sample':32s} {'N':>3s} {'used':>4s} {'stack':>8s} {'native':>8s} | sigma: sim_nb  sim_obs  rot | z(nb) z(rot) | p_upper(nb) p_upper(rot)")
    for k, r in out['samples'].items():
        print(f"{k:32s} {r['N']:3d} {r['N_used']:4d} {r['stack']:8.2f} {r['stack_native']:8.2f} | "
              f"{r['sim_nb']['sigma']:14.2f} {r['sim_obs']['sigma']:7.2f} {r['rot']['sigma']:5.2f} | "
              f"{r['sim_nb']['z']:5.2f} {r['rot']['z']:5.2f} | {r['sim_nb']['p_upper']:11.4f} {r['rot']['p_upper']:11.4f}")
    print('STEP 1:', out['decision']['step1']); print('STEP 2:', out['decision']['step2'])


if __name__ == '__main__':
    st = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if st in ('prep', 'all'): prep()
    if st in ('null', 'all'): null()
    if st in ('real', 'all'): real()
