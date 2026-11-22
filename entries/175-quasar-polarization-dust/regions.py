"""Entry 175, exploratory: the original alignment regions A1 and A3 of Hutsemekers et al. (2005), before and after the dust correction.

Made after the frozen result; nothing here decides anything.  Hutsemekers et al. (2005, A&A 441, 915, Sect. 3.1 and Table 1)
test two regions defined by eye in Paper I (Hutsemekers 1998) with a binomial count:
  A1: 11h15m <= RA(B1950) < 14h30m (names up to 1429), 1.0 <= z <= 2.3, angles expected in the sector [146, 46] deg (80 deg wide, through 0);
  A3: 21h20m <= RA(B1950) <= 24h00m, 0.7 <= z <= 1.5, angles expected in [84, 164] deg.
Under uniform angles a quasar falls in its sector with probability 80/180.  The 2005 paper reports, for p >= 0.6%:
A1 27/40 (new sample, P = 2.8e-3), 42/56 (full sample, 3.3e-6); A3 14/18 (new, 4.3e-3), 24/29 (full, 2.6e-5).
The 'new sample' leaves out the quasars of Paper I, from which the regions and sectors were chosen; it is the fair test.
Here: the counts with the measured angles, after removing the dust foreground Planck predicts (f = 1, as in the frozen run) and
at f* = 2.05; and, from the simulation of simulate.py (random intrinsic angles, Planck's foreground, the 0.6% cut, no cosmic
order), how many quasars dust alone puts in each sector.
Quasar table: data/hutsemekers-2005/qsos.tsv (VizieR J/A+A/441/915/qsos); r_Pol codes from J/A+A/441/915/refs: 0 = Paper I
(1998A&A...332..410H), 1-8 earlier literature, 9 = Lamy & Hutsemekers 2000, 10 = Smith et al. 2002, 11 = Sluse et al. 2005.
Foreground: work/template.npz (Planck GNILC, as in qso_dust.py).
Run from the repository root:  .venv/bin/python -I entries/175-quasar-polarization-dust/regions.py
"""
import json, os, sys
import numpy as np
from scipy.stats import binom
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qso_dust as q

REG = dict(A1=dict(ra=(11 + 15 / 60, 14 + 30 / 60), z=(1.0, 2.3), sector=(146.0, 46.0), paper=dict(new='27/40', full='42/56')),
           A3=dict(ra=(21 + 20 / 60, 24.0), z=(0.7, 1.5), sector=(84.0, 164.0), paper=dict(new='14/18', full='24/29')))
P0 = 80.0 / 180.0


def in_sector(pa, s):
    a, b = s
    return (pa >= a) | (pa <= b) if a > b else (pa >= a) & (pa <= b)


def pbin(k, n):
    return float(binom.sf(k - 1, n, P0)) if n else float('nan')


def rayleigh(pa_deg):
    """Sector-free test of the doubled angles: mean resultant length R, mean angle, and P = exp(-n R^2) (large-n form with the
    first correction term, Mardia & Jupp 2000, eq. 6.3.6)."""
    a = np.radians(2 * np.asarray(pa_deg)); n = len(a)
    C, S = np.cos(a).mean(), np.sin(a).mean(); R = np.hypot(C, S); Z = n * R * R
    P = np.exp(-Z) * (1 + (2 * Z - Z * Z) / (4 * n) - (24 * Z - 132 * Z ** 2 + 76 * Z ** 3 - 9 * Z ** 4) / (288 * n * n))
    return dict(n=n, R=float(R), mean_angle=float(np.mod(np.degrees(np.arctan2(S, C)) / 2, 180)), P=float(min(max(P, 0.0), 1.0)))


def vtest(pa_deg, mu_deg):
    """Directional test against a direction fixed in advance (the centre of the region's sector): V = sqrt(2n) mean cos 2(theta - mu),
    one-sided normal P (Mardia & Jupp 2000, sect. 6.3.3)."""
    from scipy.stats import norm
    c = np.cos(np.radians(2 * (np.asarray(pa_deg) - mu_deg))); V = np.sqrt(2 * len(c)) * c.mean()
    return dict(mean_cos=float(c.mean()), V=float(V), P=float(norm.sf(V)))


def main():
    from astropy.coordinates import SkyCoord, FK4
    import astropy.units as u
    z = np.load(os.path.join(HERE, 'work', 'template.npz'))
    ra, dec, zz, tq, tu, tpa = z['ra'], z['dec'], z['z'], z['q'], z['u'], z['pa']
    Q = q.load_quasars()
    head, rows = q.read_tsv(q.QSOS)
    rpol = np.array([int(r[head.index('r_Pol')]) for r in rows])
    names = np.array(Q['name'])
    assert np.allclose(Q['ra'], ra) and np.allclose(Q['z'], zz)
    p, pa, ep = Q['p'], Q['pa'], Q['ep']
    ra1950 = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs').transform_to(FK4(equinox='B1950')).ra.hour
    # 'new sample': quasars measured after Paper I, less the two Paper I objects re-measured (footnote 3 of the 2005 paper)
    # Reproduces all four of the paper's counts (A1 42/56 and 27/40, A3 24/29 and 14/18): new = r_Pol >= 7 (measurements published
    # after Paper I), less the two Paper I quasars re-measured (B1222+228; the footnote's 'B1246-047' is B1246-057 in the table);
    # the RA limits are read as B1950 names, so 14h29m includes 14h29m59s.
    paper1_remeasured = np.isin(names, ['B1222+228', 'B1246-057'])
    new = (rpol >= 7) & ~paper1_remeasured
    out = dict(sector_probability=P0, regions={})
    corr = {f: q.correct(p, pa, tq, tu, f) for f in (1.0, 2.05)}
    # simulation: dust alone (random intrinsic angles, Planck's foreground, the 0.6% cut)
    def draw(g, med, sig, gain):
        """As simulate.draw (copied: importing simulate.py would rerun it): one synthetic sample at the real positions."""
        n = len(ra); p_out, pa_out = np.empty(n), np.empty(n); todo = np.arange(n)
        while len(todo):
            m = len(todo)
            pi = g.lognormal(np.log(med), sig, m); th = g.uniform(0, np.pi, m)
            qo = pi * np.cos(2 * th) + gain * tq[todo] + g.normal(0, 1, m) * ep[todo]
            uo = pi * np.sin(2 * th) + gain * tu[todo] + g.normal(0, 1, m) * ep[todo]
            po = np.hypot(qo, uo); ok = po >= q.CFG['pcut']
            p_out[todo[ok]] = po[ok]; pa_out[todo[ok]] = np.mod(0.5 * np.arctan2(uo[ok], qo[ok]) / q.D2R, 180)
            todo = todo[~ok]
        return p_out, pa_out
    sims = {}
    vt_sims = {k: {'full': [], 'new': []} for k in REG}
    g2 = np.random.default_rng(175999)
    for it in range(2000):
        ps, pas = draw(g2, 0.8, 0.9, 1.0)
        for k, R in REG.items():
            m = (ra1950 >= R['ra'][0]) & (ra1950 < R['ra'][1]) & (zz >= R['z'][0]) & (zz <= R['z'][1])
            mu = np.mod((R['sector'][0] + 40.0), 180.0)
            for lab, mm in (('full', m), ('new', m & new)):
                vt_sims[k][lab].append(vtest(pas[mm], mu)['V'])
    for med in (0.5, 0.8):
        g = np.random.default_rng(175175 + int(1000 * med))
        cnt = {k: {'full': [], 'new': []} for k in REG}
        for it in range(2000):
            ps, pas = draw(g, med, 0.9, 1.0)
            for k, R in REG.items():
                m = (ra1950 >= R['ra'][0]) & (ra1950 < R['ra'][1]) & (zz >= R['z'][0]) & (zz <= R['z'][1])
                for lab, mm in (('full', m), ('new', m & new)):
                    cnt[k][lab].append(int(in_sector(pas[mm], R['sector']).sum()))
        sims[med] = cnt
    for k, R in REG.items():
        m = (ra1950 >= R['ra'][0]) & (ra1950 < R['ra'][1]) & (zz >= R['z'][0]) & (zz <= R['z'][1])
        rec = dict(definition=dict(ra_B1950_hours=R['ra'], z=R['z'], sector_deg=R['sector']), paper_2005=R['paper'])
        fg_in = in_sector(tpa[m], R['sector'])
        rec['foreground'] = dict(median_p=float(np.median(np.hypot(tq[m], tu[m]))), share_of_foreground_angles_in_sector=float(fg_in.mean()),
                                 mean_angle_deg=float(np.mod(0.5 * np.degrees(np.arctan2(np.mean(np.sin(2 * np.radians(tpa[m]))), np.mean(np.cos(2 * np.radians(tpa[m]))))), 180)))
        for lab, mm in (('full', m), ('new', m & new)):
            r = dict(n=int(mm.sum()))
            for tag, (pp, aa) in (('measured', (p, pa)), ('f1', corr[1.0]), ('fstar', corr[2.05])):
                k_in = int(in_sector(aa[mm], R['sector']).sum())
                r[tag] = dict(in_sector=k_in, P_binomial=pbin(k_in, int(mm.sum())))
                if tag != 'measured':                      # the cut applied again after the correction
                    keep = mm & (pp >= q.CFG['pcut'])
                    k2 = int(in_sector(aa[keep], R['sector']).sum())
                    r[tag + '_recut'] = dict(n=int(keep.sum()), in_sector=k2, P_binomial=pbin(k2, int(keep.sum())))
            for med, cnt in sims.items():
                v = np.array(cnt[k][lab])
                r['dust_only_parent_%.1f' % med] = dict(mean_in_sector=float(v.mean()), sd=float(v.std(ddof=1)),
                                                        share_at_or_above_measured=float((v >= r['measured']['in_sector']).mean()),
                                                        uniform_expectation=float(P0 * mm.sum()))
            r['rayleigh'] = dict(measured=rayleigh(pa[mm]), f1=rayleigh(corr[1.0][1][mm]))
            mu = np.mod((R['sector'][0] + 40.0), 180.0)
            r['vtest'] = dict(direction_deg=float(mu), measured=vtest(pa[mm], mu), f1=vtest(corr[1.0][1][mm], mu), fstar=vtest(corr[2.05][1][mm], mu))
            v = np.array(vt_sims[k][lab])
            r['vtest']['dust_only_parent_0.8'] = dict(mean_V=float(v.mean()), sd=float(v.std(ddof=1)), share_at_or_above_measured=float((v >= r['vtest']['measured']['V']).mean()))
            rec[lab] = r
        # stars in the same patch of sky (Heiles 2000): |b| >= 30 deg, beyond 300 pc, detected at 3 sigma; same RA range (B1950)
        St = q.load_stars()
        sra = SkyCoord(ra=St['ra'] * u.deg, dec=St['dec'] * u.deg, frame='icrs').transform_to(FK4(equinox='B1950')).ra.hour
        lo, hi = R['ra']
        ms = (sra >= lo) & (sra < hi) & (np.abs(St['glat']) >= 30) & (St['dist'] > 300) & (St['p'] >= 3 * St['ep'])
        dec_q = (dec[m].min(), dec[m].max())
        ms &= (St['dec'] >= dec_q[0] - 5) & (St['dec'] <= dec_q[1] + 5)
        rec['stars'] = dict(n=int(ms.sum()), share_in_sector=float(in_sector(St['pa'][ms], R['sector']).mean()) if ms.sum() else None,
                            rayleigh=rayleigh(St['pa'][ms]) if ms.sum() > 1 else None, dec_range_of_quasars=[float(dec_q[0]), float(dec_q[1])])
        out['regions'][k] = rec
        for lab in ('full', 'new'):
            r = rec[lab]
            print('%s %s: n %d | measured %d (P %.1e) | f=1 %d (P %.1e) recut %d/%d (P %.1e) | f*=2.05 %d (P %.1e) | dust only: %.1f +- %.1f (parent 0.5%%), %.1f +- %.1f (0.8%%); share >= measured %.3f, %.3f | uniform %.1f' % (
                k, lab, r['n'], r['measured']['in_sector'], r['measured']['P_binomial'], r['f1']['in_sector'], r['f1']['P_binomial'],
                r['f1_recut']['in_sector'], r['f1_recut']['n'], r['f1_recut']['P_binomial'], r['fstar']['in_sector'], r['fstar']['P_binomial'],
                r['dust_only_parent_0.5']['mean_in_sector'], r['dust_only_parent_0.5']['sd'], r['dust_only_parent_0.8']['mean_in_sector'],
                r['dust_only_parent_0.8']['sd'], r['dust_only_parent_0.5']['share_at_or_above_measured'], r['dust_only_parent_0.8']['share_at_or_above_measured'],
                r['dust_only_parent_0.5']['uniform_expectation']), flush=True)
        for lab in ('full', 'new'):
            ry = rec[lab]['rayleigh']
            print('   %s %s sector-free: measured mean %.0f deg, R %.2f, P %.1e | after f=1 mean %.0f, R %.2f, P %.1e' % (
                k, lab, ry['measured']['mean_angle'], ry['measured']['R'], ry['measured']['P'], ry['f1']['mean_angle'], ry['f1']['R'], ry['f1']['P']))
        for lab in ('full', 'new'):
            vt = rec[lab]['vtest']
            print('   %s %s V-test toward %.0f deg: measured V %.2f (P %.1e) | f=1 %.2f (P %.1e) | f* %.2f | dust only V %.2f +- %.2f, share >= measured %.4f' % (
                k, lab, vt['direction_deg'], vt['measured']['V'], vt['measured']['P'], vt['f1']['V'], vt['f1']['P'], vt['fstar']['V'],
                vt['dust_only_parent_0.8']['mean_V'], vt['dust_only_parent_0.8']['sd'], vt['dust_only_parent_0.8']['share_at_or_above_measured']))
        st = rec['stars']
        if st['n']:
            print('   stars in the %s patch: %d, %.0f%% in the sector, mean %.0f deg, R %.2f' % (k, st['n'], 100 * st['share_in_sector'], st['rayleigh']['mean_angle'], st['rayleigh']['R']))
        print('   foreground in %s: median %.2f%%, mean angle %.0f deg, %.0f%% of foreground angles in the sector' % (
            k, rec['foreground']['median_p'], rec['foreground']['mean_angle_deg'], 100 * rec['foreground']['share_of_foreground_angles_in_sector']))
    json.dump(out, open(os.path.join(HERE, 'regions.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
