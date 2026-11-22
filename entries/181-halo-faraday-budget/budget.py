"""Entry 181: how much of the Faraday growth with quasar redshift can intervening haloes make?

Entries 171-174 found the rotation-measure scatter of quasar sight lines growing with redshift. The union of 10,449 sight lines
gives a noise-removed variance slope of 63.8 +- 18.2 rad^2/m^4 per unit log10(1+z), about 19 rad^2/m^4 between z = 0.52 and 2.65.
Entries 176-180 measured what one crossing of a cluster or group does to a sight line. Inside R500 the scatter rises by about
120 rad^2/m^4, in proportion to M500^0.86. It is the same on quasars behind the cluster, and it shows nothing beyond R500. This
entry turns those measurements into a budget. It asks what variance growth the haloes crossed by a quasar's sight line should
produce, against what is seen.

No new data are read. It is a calculation from measured quantities:
  psi        the response of the statistic (1.4826 MAD)^2 - <e^2> of a large sample to a small share f of sight lines that cross
             a halo inside R500.  Measured by mixing field sources with the real sources inside R500 (Entry 180's sample of all
             sources around 10,151 clusters, and its 74 background quasars): psi = dV/df at small f.
  crossings  how often a sight line to redshift z_s passes within R500 of a halo of mass M at z < z_s - 0.05:
             (a) catalogued clusters (eRASS1, MCXC, PSZ2), counted on the real quasar positions of Entry 180;
             (b) all haloes above M_min, from the Tinker et al. (2008) mass function (M500c, with Delta_m = 500 / Omega_m(z)),
                 an Eisenstein & Hu (1998) no-wiggle power spectrum, Planck 2015 parameters (h 0.6774, Om 0.3089, Ob 0.0486,
                 sigma8 0.8159, ns 0.9667).
  weight     a crossing of mass M at redshift z counts w = (M / M_eff)^0.86 x [(1 + z_ref) / (1 + z)]^4, where M_eff and z_ref are
             those of the clusters behind which psi was measured. The (1+z)^-4 is the dilution of rotation made at z, which
             Entry 178's far half is consistent with. A variant with no dilution gives an upper bound.
Prediction: V_add(z_s) = psi x (expected sum of w over crossings), and its slope against log10(1 + z_s) between 0.52 and 2.65.

Run from the repository root:  .venv/bin/python -I entries/181-halo-faraday-budget/budget.py      (a few minutes)
"""
import json, os, sys, time
import numpy as np
from scipy.integrate import quad
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cluster_rm as c
import background as bg

H, OM, OB, S8, NS, TCMB = 0.6774, 0.3089, 0.0486, 0.8159, 0.9667, 2.7255
OL = 1 - OM
C_KMS = 299792.458
RHO_C0 = 2.775e11 * H ** 2            # Msun / Mpc^3
GROWTH = dict(slope=63.8, err=18.2, z_lo=0.52, z_hi=2.65, note='union of Entries 171-174 (Entry 177 section text): noise-removed variance per unit log10(1+z)')
EXPONENT = 0.86
E = lambda z: np.sqrt(OM * (1 + z) ** 3 + OL)


# ---------------------------------------------------------------- linear theory
def T_eh_nowiggle(k):
    """Eisenstein & Hu (1998) zero-baryon-like transfer function with baryon suppression; k in 1/Mpc."""
    th = TCMB / 2.7
    wm, fb = OM * H * H, OB / OM
    s = 44.5 * np.log(9.83 / wm) / np.sqrt(1 + 10 * (OB * H * H) ** 0.75)
    ag = 1 - 0.328 * np.log(431 * wm) * fb + 0.38 * np.log(22.3 * wm) * fb * fb
    gam = OM * H * (ag + (1 - ag) / (1 + (0.43 * k * s) ** 4))
    q = (k / H) * th * th / gam
    L0 = np.log(2 * np.e + 1.8 * q)
    C0 = 14.2 + 731.0 / (1 + 62.5 * q)
    return L0 / (L0 + C0 * q * q)


KG = np.logspace(-4, 2, 3000)
PK_SHAPE = KG ** NS * T_eh_nowiggle(KG) ** 2


def sigma_R(R, pk):
    x = KG[None, :] * np.atleast_1d(R)[:, None]
    W = 3 * (np.sin(x) - x * np.cos(x)) / x ** 3
    integrand = KG[None, :] ** 3 * pk[None, :] * W ** 2 / (2 * np.pi ** 2)
    return np.sqrt(np.trapezoid(integrand, np.log(KG), axis=1))


NORM = S8 ** 2 / sigma_R(np.array([8.0 / H]), PK_SHAPE)[0] ** 2
PK0 = NORM * PK_SHAPE


def growth(z):
    a = 1 / (1 + z)
    f = lambda aa: 1.0 / (aa * np.sqrt(OM / aa ** 3 + OL)) ** 3
    D = lambda aa: 2.5 * OM * np.sqrt(OM / aa ** 3 + OL) * quad(f, 1e-8, aa)[0]
    return D(a) / D(1.0)


TINKER = dict(logD=np.log10([200, 300, 400, 600, 800, 1200, 1600, 2400, 3200]),
              A=[0.186, 0.200, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260], a=[1.47, 1.52, 1.56, 1.61, 1.87, 2.13, 2.30, 2.53, 2.66],
              b=[2.57, 2.25, 2.05, 1.87, 1.59, 1.51, 1.46, 1.44, 1.41], c=[1.19, 1.27, 1.34, 1.45, 1.58, 1.80, 1.97, 2.24, 2.44])


def dndlnM(M, z):
    """Comoving number density per ln M (Mpc^-3) of haloes of mass M500c (Msun) at redshift z (Tinker et al. 2008)."""
    Omz = OM * (1 + z) ** 3 / E(z) ** 2
    Dm = 500.0 / Omz
    ld = np.log10(min(max(Dm, 200.0), 3200.0))
    A0, a0, b0, c0 = (np.interp(ld, TINKER['logD'], TINKER[k]) for k in ('A', 'a', 'b', 'c'))
    alpha = 10 ** (-(0.75 / np.log10(Dm / 75.0)) ** 1.2)
    zz = min(z, 2.5)
    A, a, b = A0 * (1 + zz) ** -0.14, a0 * (1 + zz) ** -0.06, b0 * (1 + zz) ** -alpha
    rho_m = OM * RHO_C0
    R = (3 * M / (4 * np.pi * rho_m)) ** (1 / 3)
    D = growth(z)
    lnM = np.log(M)
    sig = sigma_R(R, PK0) * D
    dlnsig = np.gradient(np.log(sig), lnM)
    f = A * ((sig / b) ** -a + 1) * np.exp(-c0 / sig ** 2)
    return f * rho_m / M * np.abs(dlnsig)


def r500_comoving(M, z):
    rho_c = RHO_C0 * E(z) ** 2
    return (3 * M / (4 * np.pi * 500 * rho_c)) ** (1 / 3) * (1 + z)          # Mpc, comoving


def crossing_kernel(Mmin, zgrid, z_ref, M_eff, dilution=True):
    """Per unit redshift of the halo: expected crossings dN/dz and weighted crossings dW/dz."""
    lnM = np.linspace(np.log(Mmin), np.log(3e15), 160)
    M = np.exp(lnM)
    dN, dW = [], []
    for z in zgrid:
        n = dndlnM(M, z)
        area = np.pi * r500_comoving(M, z) ** 2
        dchi = C_KMS / (100 * H) / E(z)
        w = (M / M_eff) ** EXPONENT * (((1 + z_ref) / (1 + z)) ** 4 if dilution else 1.0)
        dN.append(np.trapezoid(n * area, lnM) * dchi)
        dW.append(np.trapezoid(n * area * w, lnM) * dchi)
    return np.array(dN), np.array(dW)


def cumulative(zgrid, dX, zs, gap=0.05):
    out = []
    for z in zs:
        m = zgrid <= z - gap
        out.append(np.trapezoid(dX[m], zgrid[m]) if m.sum() > 1 else 0.0)
    return np.array(out)


# ---------------------------------------------------------------- measured inputs
def measured():
    rows, src, sinfo = c.load_spice()
    S, k = bg.sample_with_index(src)
    mq, de, qu = bg.load_milliquas(), bg.load_desi(), bg.load_quaia()
    zs_all, cat = bg.source_redshifts(rows, len(src['rm']), [mq, de, qu])
    zq = zs_all[k]
    cl, cinfo = bg.all_clusters()
    tree = cKDTree(c.gal_unit(S['l'], S['b']))
    x, which = c.scaled_distance(tree, len(S['d']), cl)
    field = ~np.isfinite(x)
    inside = x < 1
    bgq = inside & np.isfinite(zq) & (which >= 0)
    bgq[bgq] = zq[bgq] > cl['z'][which[bgq]] + 0.1
    return S, zq, cl, x, which, field, inside, bgq


def psi_of(dF, eF, dI, eI, fs=(0.005, 0.01, 0.02, 0.05), nfield=20000, ndraw=300, seed=181):
    """dV/df: the statistic of field sources with a share f replaced by sources inside R500 (drawn with replacement)."""
    g = np.random.default_rng(seed)
    V0 = []
    out = {}
    for f in fs:
        acc = []
        for r in range(ndraw):
            kf = g.integers(0, len(dF), nfield)
            nI = int(round(f * nfield))
            ki = g.integers(0, len(dI), nI)
            base = c.V_of(dF[kf], eF[kf])
            d = np.concatenate([dF[kf[nI:]], dI[ki]]); e = np.concatenate([eF[kf[nI:]], eI[ki]])
            acc.append((c.V_of(d, e) - base) / f)
        out[f] = (float(np.mean(acc)), float(np.std(acc, ddof=1) / np.sqrt(ndraw)))
    return out


def main():
    t0 = time.time()
    log = lambda s: print(s, flush=True)
    S, zq, cl, x, which, field, inside, bgq = measured()
    dF, eF = S['d'][field], S['e'][field]
    res = dict(growth=GROWTH, inputs=dict(field_sources=int(field.sum()), inside_R500=int(inside.sum()), background_quasars_inside=int(bgq.sum())))
    # 1. psi
    P_all = psi_of(dF, eF, S['d'][inside], S['e'][inside])
    P_bg = psi_of(dF, eF, S['d'][bgq], S['e'][bgq])
    res['psi'] = dict(all_sources={str(k): v for k, v in P_all.items()}, background_quasars={str(k): v for k, v in P_bg.items()})
    psi = P_all[0.01][0]; psi_bg = P_bg[0.01][0]
    log('psi (dV/df at f = 0.01): all sources inside R500 %.1f, background quasars %.1f rad2/m4' % (psi, psi_bg))
    # the clusters behind which psi was measured
    Mi, zi = cl['m500'][which[inside]] * 1e13, cl['z'][which[inside]]
    M_eff = float(np.mean(Mi ** EXPONENT) ** (1 / EXPONENT)); z_ref = float(np.median(zi))
    w_inside = (Mi / M_eff) ** EXPONENT * ((1 + z_ref) / (1 + zi)) ** 4
    res['reference'] = dict(M_eff=M_eff, z_ref=z_ref, mean_w_with_dilution=float(np.mean(w_inside)))
    log('reference: M_eff %.2e Msun, z_ref %.3f, mean weight of the measured crossings with dilution %.2f' % (M_eff, z_ref, np.mean(w_inside)))
    wnorm = float(np.mean(w_inside))
    # 2a. catalogued clusters crossed by the real quasar sight lines
    q = np.isfinite(zq)
    tree = cKDTree(c.gal_unit(S['l'][q], S['b'][q]))
    zqq = zq[q]
    Wcat = np.zeros(q.sum()); Ncat = np.zeros(q.sum())
    vcl = c.gal_unit(cl['l'], cl['b'])
    lists = tree.query_ball_point(vcl, c.chord(np.minimum(cl['th'], 60.0)))
    for j, L in enumerate(lists):
        for i in L:
            if cl['z'][j] < zqq[i] - 0.05:
                Ncat[i] += 1
                Wcat[i] += (cl['m500'][j] * 1e13 / M_eff) ** EXPONENT * ((1 + z_ref) / (1 + cl['z'][j])) ** 4
    edges = np.quantile(zqq, np.linspace(0, 1, 9))
    rows = []
    for a, b2 in zip(edges[:-1], edges[1:]):
        m = (zqq >= a) & (zqq <= b2)
        rows.append(dict(z_median=float(np.median(zqq[m])), n=int(m.sum()), crossing_share=float(np.mean(Ncat[m] > 0)),
                         V_add=float(psi * np.mean(Wcat[m]) / wnorm)))
    res['catalogued'] = rows
    zm = np.array([r['z_median'] for r in rows]); Va = np.array([r['V_add'] for r in rows])
    slope_cat = float(np.polyfit(np.log10(1 + zm), Va, 1)[0])
    res['catalogued_slope'] = slope_cat
    log('catalogued clusters on the real sight lines: crossing share %s; V_add %s; slope %.2f rad2/m4 per log10(1+z)' % (
        [round(r['crossing_share'], 4) for r in rows], [round(v, 2) for v in Va], slope_cat))
    # 2b. all haloes from the mass function
    zgrid = np.linspace(0.0, 3.0, 61)
    zs = np.array([GROWTH['z_lo'], GROWTH['z_hi']])
    dlog = np.log10(1 + GROWTH['z_hi']) - np.log10(1 + GROWTH['z_lo'])
    # validation: predicted crossing share of massive (M500 > 2e14) haloes at z < 0.3 against the catalogue on the same sight lines
    dN_v, _ = crossing_kernel(2e14, zgrid, z_ref, M_eff)
    m = zgrid <= 0.3
    pred_share = float(np.trapezoid(dN_v[m], zgrid[m]))
    mass_ok = (cl['m500'] * 1e13 > 2e14) & (cl['z'] < 0.3)
    hits = set()
    for j in np.nonzero(mass_ok)[0]:
        hits.update(lists[j])
    obs_share = len([i for i in hits if zqq[i] > 0.35]) / max(1, int((zqq > 0.35).sum()))
    res['validation'] = dict(predicted_share_M_gt_2e14_z_lt_0p3=pred_share, catalogued_share=obs_share,
                             note='sight lines to quasars at z > 0.35; the catalogues cover about half the sky, so their share is about half the prediction')
    log('validation (M500 > 2e14, z < 0.3): predicted crossing share %.4f all sky; catalogued %.4f' % (pred_share, obs_share))
    budget = []
    for Mmin in (1e14, 3e13, 1e13, 3e12, 1e12):
        for dil in (True, False):
            dN, dW = crossing_kernel(Mmin, zgrid, z_ref, M_eff, dilution=dil)
            Ncum, Wcum = cumulative(zgrid, dN, zs), cumulative(zgrid, dW, zs)
            Vadd = psi * Wcum / wnorm
            slope = float((Vadd[1] - Vadd[0]) / dlog)
            budget.append(dict(M_min=Mmin, dilution=dil, crossings=[float(v) for v in Ncum], V_add=[float(v) for v in Vadd], slope=slope,
                               share_of_growth=slope / GROWTH['slope']))
            log('M_min %.0e, dilution %s: crossings to z = 0.52 / 2.65: %.3f / %.3f; V_add %.2f / %.2f; slope %.1f = %.0f%% of the growth' % (
                Mmin, dil, Ncum[0], Ncum[1], Vadd[0], Vadd[1], slope, 100 * slope / GROWTH['slope']))
    res['mass_function'] = budget
    res['seconds'] = round(time.time() - t0)
    json.dump(res, open(os.path.join(HERE, 'results.json'), 'w'), indent=1, default=float)
    log('written results.json %d s' % res['seconds'])


if __name__ == '__main__':
    main()
