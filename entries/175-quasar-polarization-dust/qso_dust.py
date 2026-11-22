"""Entry 175: do the gigaparsec alignments of quasar optical polarizations survive Planck's dust-polarization foreground?

Entry 132 reproduced the alignments of Hutsemekers et al. in their 355 quasars, found the angles linked to the polarization of
nearby stars, and pre-registered this test: "with Planck's 353 GHz dust polarization as the foreground template, the full sample
staying above 3 sigma in 3D would reopen cosmic order".  The rule in full is in the Workbench, section 175, frozen before the
dust map was opened.

Data (raw files are git-ignored under data/):
  data/planck-pr3-gnilc-dust/COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits   Planck PR3 GNILC thermal dust at 353 GHz,
      80 arcmin resolution (Planck 2018 results XII), https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/
      component-maps/foregrounds/, 1,811,949,120 bytes, sha256 2d68599b73aa5034340c85eed89a5270b68aaa33dfd52a286b3a2c7320491f88.
      Columns Q_STOKES, U_STOKES (K_CMB); NESTED, Galactic, POLCCONV = COSMO.
  data/hutsemekers-2005/qsos.tsv   VizieR J/A+A/441/915/qsos (Hutsemekers et al. 2005), 355 quasars, from vizier.cfa.harvard.edu,
      40,046 bytes, sha256 19573d69c6e2dee160b1662d6c1d322c6601287a4eb723db9ca177e58a183692.
  data/heiles-2000/catalog.tsv     VizieR II/226/catalog (Heiles 2000), 9,286 stars, from vizier.cfa.harvard.edu, 1,981,311 bytes,
      sha256 092fb3d15ae0485de090976cf63ba467609ce7069ae123c542662a94a0c22200.  Used only to check the angle conventions.

Template.  Dust emission at 353 GHz is polarized across the field; starlight passing the same dust is polarized along it, with
P_353 / p_V = 5.4 MJy/sr (Planck intermediate results XXI).  So the optical polarization the Galaxy adds to a quasar is, in Stokes
terms in the IAU convention, (q, u) = -(Q, -U_map) x 287.45 MJy/sr per K_CMB / 5.4 MJy/sr, a fraction; the minus sign is the 90 degree
turn and -U_map turns Planck's COSMO convention into IAU.  The vector is then turned from the Galactic to the equatorial frame.

Statistic.  S_D of Jain et al., as in Entry 132 (code/qso_align.js, ported): each quasar's n nearest neighbours in comoving space,
angles parallel-transported, D_i = |mean exp 2i theta|, S_D = mean D_i; chance level from shuffling the angles over the positions;
the global p-value allows for trying n = 10, 20, 30, 50, 80.

Run from the repository root:
  .venv/bin/python -I entries/175-quasar-polarization-dust/qso_dust.py synthetic   # convention and statistic tests; reads no data
  .venv/bin/python -I entries/175-quasar-polarization-dust/qso_dust.py stars       # the convention check on Heiles's stars
  .venv/bin/python -I entries/175-quasar-polarization-dust/qso_dust.py run         # the test; writes results.json
"""
import json, os, sys, time
import numpy as np
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
DUST = os.path.join(ROOT, 'data', 'planck-pr3-gnilc-dust', 'COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits')
QSOS = os.path.join(ROOT, 'data', 'hutsemekers-2005', 'qsos.tsv')
STARS = os.path.join(ROOT, 'data', 'heiles-2000', 'catalog.tsv')
D2R = np.pi / 180
CFG = dict(nvs=[10, 20, 30, 50, 80], nshuf=10000, nshuf_control=1000, seed=20261121, disc_deg=0.5, K_to_MJysr=287.45, R_P_over_p=5.4,
           f_scan=[round(0.25 * i, 2) for i in range(0, 41)], f_allow=1.5, ncontrol=10, pcut=0.6)
MG = np.array([[-0.0548755604, -0.8734370902, -0.4838350155],
               [0.4941094279, -0.4448296300, 0.7469822445],
               [-0.8676661490, -0.1980763734, 0.4559837762]])          # J2000 equatorial -> Galactic


# ---------------------------------------------------------------- geometry
def unit(lon, lat):
    a, d = np.asarray(lon) * D2R, np.asarray(lat) * D2R
    return np.stack([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)], axis=-1)


def basis(lon, lat):
    """North and east unit vectors at (lon, lat), in that frame's Cartesian axes."""
    a, d = np.asarray(lon) * D2R, np.asarray(lat) * D2R
    N = np.stack([-np.sin(d) * np.cos(a), -np.sin(d) * np.sin(a), np.cos(d)], axis=-1)
    E = np.stack([-np.sin(a), np.cos(a), np.zeros_like(a)], axis=-1)
    return N, E


def to_galactic(ra, dec):
    v = unit(ra, dec) @ MG.T
    return np.mod(np.arctan2(v[..., 1], v[..., 0]) / D2R, 360.0), np.arcsin(np.clip(v[..., 2], -1, 1)) / D2R


def pa_gal_to_eq(pa_gal_deg, ra, dec):
    """Position angle from Galactic north through Galactic east -> from equatorial north through east, at (ra, dec)."""
    l, b = to_galactic(ra, dec)
    Ng, Eg = basis(l, b)                                   # in Galactic axes
    th = np.asarray(pa_gal_deg) * D2R
    d = (np.cos(th)[..., None] * Ng + np.sin(th)[..., None] * Eg) @ MG      # to equatorial axes (MG is orthogonal: MG.T^-1 = MG)
    Ne, Ee = basis(ra, dec)
    return np.mod(np.arctan2((d * Ee).sum(-1), (d * Ne).sum(-1)) / D2R, 180.0)


def transport_offset(ra_a, dec_a, ra_b, dec_b):
    """Turn picked up by a position angle carried from a to b along the great circle (Entry 132's offset(a, b)), radians."""
    va, vb = unit(ra_a, dec_a), unit(ra_b, dec_b)
    ax = np.cross(va, vb); sn = np.linalg.norm(ax, axis=-1)
    Na, _ = basis(ra_a, dec_a); Nb, Eb = basis(ra_b, dec_b)
    ok = sn > 1e-12
    k = np.where(ok[..., None], ax / np.where(ok, sn, 1.0)[..., None], 0.0)
    ang = np.arctan2(sn, (va * vb).sum(-1))
    c, s = np.cos(ang)[..., None], np.sin(ang)[..., None]
    t = Na * c + np.cross(k, Na) * s + k * (k * Na).sum(-1)[..., None] * (1 - c)
    return np.where(ok, np.arctan2((t * Eb).sum(-1), (t * Nb).sum(-1)), 0.0)


def comoving(z, n=400):
    z = np.asarray(z, float)
    zz = (np.arange(n) + 0.5)[None, :] * z[:, None] / n
    return 2997.92458 / 0.7 * (1 / np.sqrt(0.3 * (1 + zz) ** 3 + 0.7)).sum(1) * z / n


# ---------------------------------------------------------------- the alignment statistic
def prepare(ra, dec, z, nvmax, mode='3d'):
    u = unit(ra, dec)
    if mode == '3d':
        x = u * comoving(z)[:, None]
        dist = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=-1)
    else:
        dist = np.arccos(np.clip(u @ u.T, -1, 1))
    nb = np.argsort(dist, axis=1, kind='stable')[:, :nvmax]
    off = transport_offset(ra[nb], dec[nb], ra[:, None] * np.ones_like(nb), dec[:, None] * np.ones_like(nb))
    r = np.take_along_axis(dist, nb, axis=1)
    return dict(nb=nb, off=off, r=r)


def SD(prep, ang, nvs):
    c = np.cumsum(np.exp(2j * (ang[prep['nb']] + prep['off'])), axis=1)
    return np.array([np.abs(c[:, nv - 1]).mean() / nv for nv in nvs])


def align_test(prep, ang, nvs, nshuf, seed):
    g = np.random.default_rng(seed)
    obs = SD(prep, ang, nvs)
    sims = np.array([SD(prep, g.permutation(ang), nvs) for _ in range(nshuf)])
    mu, sd = sims.mean(0), sims.std(0, ddof=1)
    p = (1 + (sims >= obs[None, :]).sum(0)) / (nshuf + 1)
    # global p: how often a shuffle's best n is at least as extreme as the data's best (ranks within the shuffles)
    rank = (1 + (sims[None, :, :] >= sims[:, None, :]).sum(1)) / (nshuf + 1) if nshuf <= 2000 else None
    if rank is None:                                       # same thing without the N x N array
        rank = np.empty_like(sims)
        for j in range(sims.shape[1]):
            order = np.argsort(-sims[:, j], kind='stable'); rk = np.empty(nshuf); rk[order] = np.arange(1, nshuf + 1)
            rank[:, j] = (1 + rk) / (nshuf + 1)            # 1 + number of shuffles at least as large (itself counted)
    pglob = float((1 + (rank.min(1) <= p.min()).sum()) / (nshuf + 1))
    return dict(nvs=list(nvs), S_D=obs.tolist(), sigma=((obs - mu) / sd).tolist(), p=p.tolist(), p_global=pglob,
                sigma_global=float(norm.isf(pglob)), nshuf=nshuf, median_radius_Mpc=np.median(prep['r'][:, [nv - 1 for nv in nvs]], axis=0).tolist())


# ---------------------------------------------------------------- the foreground template
def optical_from_map(Qmap, Umap):
    """Planck Q, U (K_CMB, COSMO) -> optical q, u in percent, IAU convention, Galactic frame."""
    k = -CFG['K_to_MJysr'] / CFG['R_P_over_p'] * 100.0
    return k * Qmap, k * (-Umap)


def template_at(ra, dec, maps=None):
    import healpy as hp
    if maps is None:
        from astropy.io import fits
        with fits.open(DUST, memmap=True) as h:
            Q = np.asarray(h[1].data.field('Q_STOKES'), dtype=np.float64); U = np.asarray(h[1].data.field('U_STOKES'), dtype=np.float64)
            bad = float(h[1].header.get('BAD_DATA', -1.6375e30))
        maps = (Q, U, bad)
    Q, U, bad = maps
    nside = hp.npix2nside(len(Q))
    l, b = to_galactic(ra, dec)
    qg, ug, npix = np.empty(len(l)), np.empty(len(l)), np.empty(len(l), int)
    for i in range(len(l)):
        pix = hp.query_disc(nside, hp.ang2vec(l[i], b[i], lonlat=True), CFG['disc_deg'] * D2R, nest=True)
        ok = (Q[pix] > -1e29) & (U[pix] > -1e29)
        qo, uo = optical_from_map(Q[pix][ok].mean(), U[pix][ok].mean())
        qg[i], ug[i], npix[i] = qo, uo, int(ok.sum())
    p = np.hypot(qg, ug)
    pa_g = np.mod(0.5 * np.arctan2(ug, qg) / D2R, 180.0)
    pa_e = pa_gal_to_eq(pa_g, ra, dec)
    return dict(p=p, pa_gal=pa_g, pa=pa_e, q=p * np.cos(2 * pa_e * D2R), u=p * np.sin(2 * pa_e * D2R), npix=npix, l=l, b=b), maps


def correct(p, pa_deg, fq, fu, f, phase=None):
    """Subtract f times the foreground (q, u); phase turns the foreground vector in Stokes space (random-direction control)."""
    if phase is not None:
        fq, fu = fq * np.cos(phase) - fu * np.sin(phase), fq * np.sin(phase) + fu * np.cos(phase)
    q = p * np.cos(2 * pa_deg * D2R) - f * fq; u = p * np.sin(2 * pa_deg * D2R) - f * fu
    return np.hypot(q, u), np.mod(0.5 * np.arctan2(u, q) / D2R, 180.0)


def link(pa_deg, pa_fg_deg):
    c = np.cos(2 * (pa_deg - pa_fg_deg) * D2R)
    return float(c.mean()), float(1 / np.sqrt(2 * len(c)))


def f_star(p, pa, T):
    """Smallest f on the grid's first sign change of the link between corrected angles and the foreground direction."""
    cs = [link(correct(p, pa, T['q'], T['u'], f)[1], T['pa'])[0] for f in CFG['f_scan']]
    for i in range(1, len(cs)):
        if cs[i - 1] > 0 >= cs[i]:
            f0, f1 = CFG['f_scan'][i - 1], CFG['f_scan'][i]
            return float(f0 + (f1 - f0) * cs[i - 1] / (cs[i - 1] - cs[i])), cs
    return None, cs


def decide(S1, fs, Sfs):
    """Section 175, item 6."""
    if S1 < 2:
        return 'Galactic: the foreground Planck predicts removes the alignment (S(1) < 2 sigma)'
    if S1 < 3:
        return 'unresolved: with the foreground Planck predicts removed the alignment stands between 2 and 3 sigma'
    if fs is None or Sfs >= 3:
        return 'the alignment survives the Planck foreground at 3 sigma or more, at f = 1 and at the level that ends the link: cosmic order is reopened'
    if fs <= CFG['f_allow']:
        return 'Galactic within the template\'s uncertainty: the alignment survives f = 1 but not the level that ends the link, f* <= 1.5'
    return 'unresolved: the alignment survives the foreground Planck predicts and is ended only by a larger one (f* > 1.5)'


# ---------------------------------------------------------------- tables
def sex(s, hours=False):
    s = s.strip(); sg = -1.0 if s.startswith('-') else 1.0
    a, b, c = [float(x) for x in s.lstrip('+-').split()]
    return sg * (a + b / 60 + c / 3600) * (15.0 if hours else 1.0)


def read_tsv(path):
    rows = [l.rstrip('\n').split('\t') for l in open(path, encoding='utf-8', errors='replace') if l.strip() and not l.startswith('#')]
    head = [h.strip() for h in rows[0]]
    return head, [r for r in rows[3:] if len(r) >= len(head) - 1]


def load_quasars():
    head, rows = read_tsv(QSOS)
    ix = {k: head.index(k) for k in ('Name', 'z', 'Pol', 'e_Pol', 'PPA', 'e_PPA', 'RAJ2000', 'DEJ2000')}
    f = lambda r, k: float(r[ix[k]])
    return dict(name=[r[ix['Name']].strip() for r in rows], z=np.array([f(r, 'z') for r in rows]), p=np.array([f(r, 'Pol') for r in rows]),
                ep=np.array([f(r, 'e_Pol') for r in rows]), pa=np.array([f(r, 'PPA') for r in rows]),
                ra=np.array([sex(r[ix['RAJ2000']], True) for r in rows]), dec=np.array([sex(r[ix['DEJ2000']]) for r in rows]))


def load_stars():
    head, rows = read_tsv(STARS)
    ix = {k: head.index(k) for k in ('RAJ2000', 'DEJ2000', 'Pol', 'e_Pol', 'PA', 'PAgal', 'GLAT', 'Dist')}
    out = {k: [] for k in ('ra', 'dec', 'p', 'ep', 'pa', 'pagal', 'glat', 'dist')}
    for r in rows:
        try:
            v = dict(ra=sex(r[ix['RAJ2000']], True), dec=sex(r[ix['DEJ2000']]), p=float(r[ix['Pol']]), ep=float(r[ix['e_Pol']]),
                     pa=float(r[ix['PA']]), pagal=float(r[ix['PAgal']]), glat=float(r[ix['GLAT']]), dist=float(r[ix['Dist']]))
        except (ValueError, IndexError):
            continue
        for k in out:
            out[k].append(v[k])
    return {k: np.array(v) for k, v in out.items()}


# ---------------------------------------------------------------- stages
def synthetic():
    import healpy as hp
    out = {}
    g = np.random.default_rng(5)
    # 1. frame rotation: Galactic north at the equatorial pole region and a round trip through vectors
    ra, dec = g.uniform(0, 360, 2000), np.arcsin(g.uniform(-1, 1, 2000)) / D2R
    l, b = to_galactic(ra, dec)
    # the direction to the north Galactic pole has Galactic PA 0; its equatorial PA from spherical trigonometry
    a0, d0 = 192.85948, 27.12825
    pa_ngp = np.mod(np.arctan2(np.sin((a0 - ra) * D2R), np.cos(dec * D2R) * np.tan(d0 * D2R) - np.sin(dec * D2R) * np.cos((a0 - ra) * D2R)) / D2R, 180)
    d = pa_gal_to_eq(np.zeros(2000), ra, dec) - pa_ngp
    out['frame_rotation_max_error_deg'] = float(np.abs((d + 90) % 180 - 90).max())
    # 2. map convention: a field along PA_B (IAU, Galactic) makes emission polarized at PA_B + 90; the code must return PA_B
    errs = []
    for paB in (0.0, 45.0, 90.0, 130.0):
        em = (paB + 90.0) * D2R
        Qi, Ui = np.cos(2 * em), np.sin(2 * em)            # IAU emission Stokes
        qo, uo = optical_from_map(Qi, -Ui)                 # the map stores COSMO: U_map = -U_IAU
        back = np.mod(0.5 * np.arctan2(uo, qo) / D2R, 180)
        errs.append(float(abs((back - paB + 90) % 180 - 90)))
    out['map_convention_max_error_deg'] = max(errs)
    out['unit_check_percent_per_10uK'] = float(np.hypot(*optical_from_map(10e-6, 0.0)))
    # 3. transport: along a meridian the offset is zero; along the equator it is zero; across the pole it is not
    out['transport_meridian'] = float(transport_offset(np.array([30.0]), np.array([10.0]), np.array([30.0]), np.array([50.0]))[0])
    out['transport_equator'] = float(transport_offset(np.array([30.0]), np.array([0.0]), np.array([80.0]), np.array([0.0]))[0])
    out['transport_same_dec_60'] = float(transport_offset(np.array([0.0]), np.array([60.0]), np.array([40.0]), np.array([60.0]))[0] / D2R)
    # 4. statistic: random angles give no signal; a foreground added to weakly polarized quasars gives one; removing it at f = 1 ends it
    n = 355
    ra, dec = g.uniform(0, 360, n), np.arcsin(g.uniform(0.5, 1, n)) / D2R * g.choice([-1, 1], n)
    z = g.uniform(0.3, 2.5, n)
    prep = prepare(ra, dec, z, 80)
    p_int, pa_int = g.uniform(0.3, 2.0, n), g.uniform(0, 180, n)
    r0 = align_test(prep, pa_int * D2R, CFG['nvs'], 500, 1)
    fg_pa = pa_gal_to_eq(np.full(n, 90.0), ra, dec)        # a field along Galactic longitude everywhere
    fq, fu = 0.5 * np.cos(2 * fg_pa * D2R), 0.5 * np.sin(2 * fg_pa * D2R)
    p_obs, pa_obs = correct(p_int, pa_int, fq, fu, -1.0)   # add the foreground
    r1 = align_test(prep, pa_obs * D2R, CFG['nvs'], 500, 2)
    T = dict(q=fq, u=fu, pa=fg_pa)
    fs, cs = f_star(p_obs, pa_obs, T)
    p2, pa2 = correct(p_obs, pa_obs, fq, fu, 1.0)
    r2 = align_test(prep, pa2 * D2R, CFG['nvs'], 500, 3)
    out['statistic'] = dict(random=[round(x, 2) for x in r0['sigma']], random_p_global=r0['p_global'],
                            with_foreground=[round(x, 2) for x in r1['sigma']], with_foreground_p_global=r1['p_global'],
                            link_before=round(cs[0], 3), f_star=fs, removed_at_f1=[round(x, 2) for x in r2['sigma']],
                            removed_p_global=r2['p_global'], angles_restored_max_error_deg=float(np.abs((pa2 - pa_int + 90) % 180 - 90).max()))
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(HERE, 'synthetic.json'), 'w'), indent=1)


def stars():
    """Convention check on real data that are not the quasars: Planck's predicted optical angle against Heiles's stars."""
    S = load_stars()
    out = {}
    d = pa_gal_to_eq(S['pagal'], S['ra'], S['dec']) - S['pa']
    d = (d + 90) % 180 - 90
    out['catalogue_galactic_to_equatorial'] = dict(n=int(len(d)), median_abs_deg=float(np.median(np.abs(d))), within_2deg=float((np.abs(d) < 2).mean()))
    sel = (np.abs(S['glat']) >= 30) & (S['dist'] >= 300) & (S['p'] >= 3 * S['ep']) & (S['ep'] > 0)
    T, maps = template_at(S['ra'][sel], S['dec'][sel])
    c, e = link(S['pa'][sel], T['pa'])
    out['high_latitude'] = dict(n=int(sel.sum()), mean_cos2delta=c, error=e, median_p_star=float(np.median(S['p'][sel])),
                                median_p_predicted=float(np.median(T['p'])), median_ratio=float(np.median(S['p'][sel] / T['p'])))
    sel2 = (np.abs(S['glat']) < 30) & (S['dist'] >= 1000) & (S['p'] >= 5 * S['ep']) & (S['ep'] > 0)
    T2, _ = template_at(S['ra'][sel2], S['dec'][sel2], maps)
    c2, e2 = link(S['pa'][sel2], T2['pa'])
    out['low_latitude_distant'] = dict(n=int(sel2.sum()), mean_cos2delta=c2, error=e2, median_ratio=float(np.median(S['p'][sel2] / T2['p'])))
    out['passes'] = bool(c >= 0.2 and out['catalogue_galactic_to_equatorial']['within_2deg'] > 0.9)
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(HERE, 'stars_check.json'), 'w'), indent=1)
    return out


def run():
    t0 = time.time()
    chk = json.load(open(os.path.join(HERE, 'stars_check.json')))
    if not chk['passes']:
        raise SystemExit('the convention check on stars has not passed; fix and list the fix before any quasar is corrected')
    Q = load_quasars()
    n = len(Q['z'])
    nvs = CFG['nvs']
    prep = prepare(Q['ra'], Q['dec'], Q['z'], max(nvs))
    T, maps = template_at(Q['ra'], Q['dec'])
    res = dict(generated=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), config=CFG, n=n, stars_check=chk,
               template=dict(p_median=float(np.median(T['p'])), p_quartiles=[float(x) for x in np.percentile(T['p'], [25, 75])],
                             p_max=float(T['p'].max()), quasar_p_median=float(np.median(Q['p'])), min_pixels=int(T['npix'].min())))
    log = lambda s: print(s, flush=True)
    log('template (percent): median %.3f, quartiles %.3f-%.3f, max %.3f; quasars median %.2f' % (
        res['template']['p_median'], *res['template']['p_quartiles'], res['template']['p_max'], res['template']['quasar_p_median']))
    fs, cs = f_star(Q['p'], Q['pa'], T)
    res['link'] = dict(f=CFG['f_scan'], mean_cos2delta=cs, error=float(1 / np.sqrt(2 * n)), f_star=fs)
    log('link at f = 0: %+.3f +- %.3f; f* = %s' % (cs[0], res['link']['error'], fs))
    levels = [('f0', 0.0), ('f1', 1.0)] + ([('fstar', fs)] if fs is not None else [])
    res['alignment'] = {}
    for name, f in levels:
        p2, pa2 = correct(Q['p'], Q['pa'], T['q'], T['u'], f)
        r = align_test(prep, pa2 * D2R, nvs, CFG['nshuf'], CFG['seed'] + int(round(100 * f)))
        r['f'] = f; r['below_cut_after'] = int((p2 < CFG['pcut']).sum())
        res['alignment'][name] = r
        log('%s (f = %.2f): sigma %s, p_global %.5f (%.2f sigma)' % (name, f, [round(x, 2) for x in r['sigma']], r['p_global'], r['sigma_global']))
    S1 = res['alignment']['f1']['sigma_global']
    Sfs = res['alignment']['fstar']['sigma_global'] if fs is not None else None
    res['decision'] = dict(S_0=res['alignment']['f0']['sigma_global'], S_1=S1, f_star=fs, S_fstar=Sfs, outcome=decide(S1, fs, Sfs))
    log('decision %s' % json.dumps(res['decision']))
    rep = {}
    g = np.random.default_rng(CFG['seed'] + 9)
    for name, f in levels[1:]:
        trials = []
        for k in range(CFG['ncontrol']):
            p2, pa2 = correct(Q['p'], Q['pa'], T['q'], T['u'], f, phase=g.uniform(0, 2 * np.pi, n))
            r = align_test(prep, pa2 * D2R, nvs, CFG['nshuf_control'], CFG['seed'] + 500 + k)
            trials.append(dict(sigma=r['sigma'], sigma_global=r['sigma_global']))
        rep['random_direction_' + name] = dict(f=f, sigma_n30=[round(t['sigma'][2], 2) for t in trials],
                                               sigma_global=[round(t['sigma_global'], 2) for t in trials])
        log('control %s: n = 30 sigma %s' % (name, rep['random_direction_' + name]['sigma_n30']))
    p2, pa2 = correct(Q['p'], Q['pa'], T['q'], T['u'], 1.0)
    keep = p2 >= CFG['pcut']
    prep_k = prepare(Q['ra'][keep], Q['dec'][keep], Q['z'][keep], max(nvs))
    rep['f1_cut_reapplied'] = dict(n=int(keep.sum()), **align_test(prep_k, pa2[keep] * D2R, nvs, CFG['nshuf_control'] * 3, CFG['seed'] + 71))
    prep2 = prepare(Q['ra'], Q['dec'], Q['z'], max(nvs), mode='2d')
    for name, f in levels:
        p2, pa2 = correct(Q['p'], Q['pa'], T['q'], T['u'], f)
        rep['sky_only_' + name] = align_test(prep2, pa2 * D2R, nvs, CFG['nshuf_control'] * 3, CFG['seed'] + 80 + int(round(100 * f)))
    rep['link_by_polarization'] = []
    for lo, hi in ((0.6, 1.0), (1.0, 1.6), (1.6, 3.0), (3.0, 99.0)):
        m = (Q['p'] >= lo) & (Q['p'] < hi)
        c, e = link(Q['pa'][m], T['pa'][m])
        rep['link_by_polarization'].append(dict(p_lo=lo, p_hi=hi, n=int(m.sum()), mean_cos2delta=c, error=e))
    res['reported_only'] = rep
    res['seconds'] = round(time.time() - t0)
    json.dump(res, open(os.path.join(HERE, 'results.json'), 'w'), indent=1, default=float)
    os.makedirs(os.path.join(HERE, 'work'), exist_ok=True)
    np.savez(os.path.join(HERE, 'work', 'template.npz'), **{k: v for k, v in T.items()}, qp=Q['p'], qpa=Q['pa'], ra=Q['ra'], dec=Q['dec'], z=Q['z'])
    log('written results.json %d s' % res['seconds'])


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    {'synthetic': synthetic, 'stars': stars, 'run': run}.get(cmd, lambda: print(__doc__))()
