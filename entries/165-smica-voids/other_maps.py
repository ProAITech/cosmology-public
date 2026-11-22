"""Entry 165, completed: the frozen deeper-slab void test repeated on Planck's NILC, SEVEM and Commander maps.

Entry 148's protocol reads "Planck PR3 SMICA temperature, common mask, monopole and dipole removed; repeated in NILC, SEVEM and
Commander".  smica_voids.py ran SMICA (13 November 2026).  This script runs the other three with the same code path: it imports
smica_voids.py and calls its functions, changing only the map file.  Choices: README.md, section "The repeat in three more maps",
fixed 21 November 2026 before any of the three maps was opened.

Data (URL, size, sha256 in README.md): data/planck-pr3-{nilc,sevem,commander}/COM_CMB_IQU-{...}_2048_R3.00_full.fits, table 1,
column I_STOKES (K_CMB); the common mask as in smica_voids.py.
Needs work/nulls.npz from `smica_voids.py null` (the simulated-sky nulls depend only on the mask and the void positions).

Run from the repository root:  .venv/bin/python entries/165-smica-voids/other_maps.py      (about 5 minutes)
"""
import sys, json, time
from pathlib import Path
import numpy as np
import healpy as hp
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import smica_voids as sv

MAPS = {m: sv.DATA / ('planck-pr3-' + m) / ('COM_CMB_IQU-%s_2048_R3.00_full.fits' % m) for m in ('nilc', 'sevem', 'commander')}


def lowres(path, good):
    from astropy.io import fits
    with fits.open(path, memmap=True) as h:
        names = [c.name for c in h[1].columns]
        I = np.asarray(h[1].data.field('I_STOKES'), dtype=np.float64) * 1e6
        order = h[1].header.get('ORDERING'); nside = h[1].header.get('NSIDE'); units = h[1].header.get('TUNIT1')
    assert order == 'NESTED' and nside == sv.NNAT, (order, nside)
    low = hp.ud_grade(np.where(good, I, sv.UNSEEN), sv.NS, pess=True, order_in='NESTED', order_out='NESTED')
    valid = low > -1e20
    nvec = np.array(hp.pix2vec(sv.NS, np.arange(hp.nside2npix(sv.NS)), nest=True))
    low_c, mono, dip = sv.md_subtract(low, valid, nvec)
    return I, low_c, valid, float(mono), dip, dict(columns=names, ordering=order, nside=nside, unit=units)


def one_map(name, path, good, c, S, nl):
    t0 = time.time()
    I, low, valid, mono, dip, meta = lowres(path, good)
    idx = [sv.regions(sv.NS, c['l'][i], c['b'][i], c['th'][i], valid) for i in range(len(c['th']))]
    T = np.array([sv.stat(low, *ix) for ix in idx])
    Tn = sv.native_T(c, good, I, mono, dip)
    rot = []
    for i in range(sv.NROT):                                   # the same shifts as smica_voids.rot_task
        dl = np.random.default_rng(sv.SEED + 700000 + i).uniform(0, 360)
        rot.append([sv.stat(low, *sv.regions(sv.NS, (c['l'][j] + dl) % 360, c['b'][j], c['th'][j], valid)) for j in range(len(c['th']))])
    rot = np.array(rot)
    out = dict(map=name, file=str(path.relative_to(sv.ROOT)), meta=meta, used_fraction_nside256=float(valid.mean()),
               monopole=mono, dipole=float(np.linalg.norm(dip)), sky_rms=float(low[valid].std()),
               T=[None if np.isnan(x) else float(x) for x in T], samples={})
    for sname, m in S.items():
        v, n = sv.stack(T, m); vn, nn = sv.stack(Tn, m)
        r = dict(N=int(m.sum()), N_used=n, stack=v, stack_native=vn)
        for key, arr in (('sim_nb', nl['sim_nb']), ('sim_obs', nl['sim_obs']), ('rot', rot)):
            ns_ = np.array([sv.stack(row, m)[0] for row in arr]); ns_ = ns_[np.isfinite(ns_)]
            r[key] = dict(mean=float(ns_.mean()), sigma=float(ns_.std(ddof=1)), z=float((v - ns_.mean()) / ns_.std(ddof=1)),
                          p_upper=float((1 + (ns_ >= v).sum()) / (1 + len(ns_))))
        out['samples'][sname] = r
    chk = out['samples']['check_slab1_largest_quarter']['stack']; pri = out['samples']['primary_fresh']['stack']
    sig = out['samples']['primary_fresh']['sim_nb']['sigma']
    out['decision'] = dict(
        step1_check_value=chk, step1=('the claim depends on the void finder (under +10 uK)' if chk < 10 else
                                      'the nearest-slab warmth is reproduced with these voids (+10 uK or more)'),
        step2_value=pri, step2=('real effect along the path (+20 uK or more)' if pri >= 20 else
                                'chance alignment, or the nearest structures only (below +15 uK)' if pri < 15 else
                                'unresolved (between +15 and +20 uK)'),
        primary_sigma_sim=sig, excludes_35uK_at_sigma=float((35.0 - pri) / sig))
    print('%-9s used %.4f, rms %.1f uK | slab-1 largest quarter %+.2f | fresh 43: %+.2f +- %.2f (native %+.2f; rot sigma %.2f) | %s | %s  [%d s]' % (
        name, out['used_fraction_nside256'], out['sky_rms'], chk, pri, sig, out['samples']['primary_fresh']['stack_native'],
        out['samples']['primary_fresh']['rot']['sigma'], out['decision']['step1'], out['decision']['step2'], time.time() - t0), flush=True)
    return out, T


def main():
    from astropy.io import fits
    with fits.open(sv.MASK, memmap=True) as h: good = np.asarray(h[1].data.field(0), dtype=np.float32) == 1
    c = sv.catalogue(); S = sv.samples(c); nl = np.load(sv.WORK / 'nulls.npz')
    sm = json.load(open(HERE / 'results.json'))
    Tsm = np.array([np.nan if v['T'] is None else v['T'] for v in sm['voids']])
    res = dict(generated=time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime()), maps={},
               smica=dict(primary_fresh=sm['samples']['primary_fresh']['stack'], check=sm['samples']['check_slab1_largest_quarter']['stack'],
                          slab2_all=sm['samples']['slab2_all']['stack'], overlapping=sm['samples']['overlapping']['stack']))
    Ts = {'smica': Tsm}
    for name, path in MAPS.items():
        res['maps'][name], Ts[name] = one_map(name, path, good, c, S, nl)
    # how closely the four methods agree, void by void, on the fresh sample
    fresh = S['primary_fresh']
    res['agreement_fresh'] = {a + '_minus_smica': dict(mean=float(np.nanmean((Ts[a] - Tsm)[fresh])), rms=float(np.nanstd((Ts[a] - Tsm)[fresh])))
                              for a in MAPS}
    allv = np.array([res['maps'][m]['samples']['primary_fresh']['stack'] for m in MAPS] + [res['smica']['primary_fresh']])
    res['primary_fresh_four_maps'] = dict(values=allv.tolist(), min=float(allv.min()), max=float(allv.max()), mean=float(allv.mean()))
    res['all_below_15uK'] = bool((allv < 15).all())
    json.dump(res, open(HERE / 'results_other_maps.json', 'w'), indent=1)
    print('fresh 43, four maps (nilc, sevem, commander, smica):', [round(x, 2) for x in allv], '| agreement', json.dumps(res['agreement_fresh']))


if __name__ == '__main__':
    main()
