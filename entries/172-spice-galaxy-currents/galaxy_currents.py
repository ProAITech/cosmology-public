"""Entry 172: a current through nearby galaxies along their spin?  Entry 61's galaxy-by-galaxy fit, repeated on SPICE-RACS DR2.

Entry 61 (code/possum_through.js, results code/possum_results.json) fitted, around each galaxy with a known spin, raw rotation
measures between 30 and 200 kpc with a quadratic foreground plus the two-lobe pattern of a current along the spin axis through a
halo of gas n_e = N0 (r/kpc)^-1.5 (N0 = 3.68e-2 cm^-3, the "standard halo"); errors from the same fit at blank places at random
angles.  Seven southern galaxies in POSSUM's first fields gave +0.22 +- 0.13 x 1e18 A, all positive (a lean of about 2 sigma).
This script repeats that fit on SPICE-RACS DR2, which covers the whole sky south of Dec +49 deg, for every galaxy of known spin
within reach.  The rule is frozen in the Workbench (section 172) before any rotation measure near a galaxy is fitted.

Data:
  data/spice-racs-dr2/spice-racs.dr2.fits.gz  SPICE-RACS DR2 (Thomson et al. 2026, arXiv:2605.16917), CSIRO DAP csiro:64891,
      4,968,301,262 bytes, sha256 a1a2f3c1d12a8a8b3aa8526d74b8e0115054f610e2e4186108f0f068005d0053.  Entry 166's quality cuts
      (snr_flag, leakage_flag, channel_flag, stokesI_fit_flag all False; finite rm, rm_err > 0); rows merged per source_id by
      inverse-variance weighted mean at the first row's position.  extract() streams the table once (about 1 minute).
  galaxies.txt  LVHIS rotation fits (Koribalski et al. 2018, VizieR J/MNRAS/478/1611 tables 9 and 2) with UNGC distances
      (Karachentsev et al. 2013, VizieR J/AJ/145/101), and PHANGS kinematic fits within 12 Mpc (Lang et al. 2020,
      VizieR J/ApJ/897/122 table1); read from VizieR in browser memory and typed in, as in Entry 61.

Run from the repository root:
  .venv/bin/python -I entries/172-spice-galaxy-currents/galaxy_currents.py synthetic   # injection tests on a synthetic sky
  .venv/bin/python -I entries/172-spice-galaxy-currents/galaxy_currents.py extract
  .venv/bin/python -I entries/172-spice-galaxy-currents/galaxy_currents.py run         # writes results.json
"""
import gzip, json, os, sys, time
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SPICE = os.path.join(ROOT, 'data', 'spice-racs-dr2', 'spice-racs.dr2.fits.gz')
WORK = os.path.join(HERE, 'work')
DEG = np.pi / 180
KRM = 0.812 * 6.4815e4
N0, BCORE = 3.68e-2, 10.0
O = dict(bin=30.0, bout=200.0, order=2, clip=3.0, scale=100.0)
CUTS = dict(imin=40.0, bmin=20.0, Dmin=1.5, Dmax=80.0, nmin=30, qmin=0.10)
NULL = dict(frac=0.5, nang=6, max_places=400, bmin=20.0, seed=1234)
STACK = dict(ndraw=200000, seed=1357)
LEAN = 0.20            # Entries 60-61, north and south together: +0.20 +- 0.11 (1e18 A)
ENTRY61 = ['NGC300', 'ESO154-G023', 'NGC1311', 'IC1959', 'IC1954', 'NGC685', 'NGC1672']

MG = np.array([[-0.0548755604, -0.8734370902, -0.4838350155],
               [0.4941094279, -0.4448296300, 0.7469822445],
               [-0.8676661490, -0.1980763734, 0.4559837762]])


def unit(ra, dec):
    a, d = np.asarray(ra, float) * DEG, np.asarray(dec, float) * DEG
    return np.stack([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)], axis=-1)


def gal_b(ra, dec):
    v = unit(ra, dec) @ MG.T
    return np.arcsin(np.clip(v[..., 2], -1, 1)) / DEG


def chord(deg):
    return 2 * np.sin(0.5 * deg * DEG)


# ---------------------------------------------------------------- the template (Entry 57's, as coded in Entry 61)
_TABS = {}


def ftab(i_deg):
    key = int(max(1, min(90, round(i_deg))))
    if key in _TABS:
        return _TABS[key]
    si, ci = np.sin(key * DEG), np.cos(key * DEG)
    tab = np.empty(180)
    for a in range(180):
        th = (a + 0.5) * 0.5 * DEG
        ct, st = np.cos(th), np.sin(th)
        n = 40000 if a >= 170 else 4000
        h = np.pi / n
        ph = -np.pi / 2 + (np.arange(n) + 0.5) * h
        cp, sp = np.cos(ph), np.sin(ph)
        q = st * ci * cp - sp * si
        tab[a] = np.sum(cp ** 1.5 * ct * si / (ct * ct * cp * cp + q * q)) * h
    _TABS[key] = tab
    return tab


def fval(tab, x, y):
    th = np.arctan2(np.abs(y), np.abs(x)) / DEG
    return np.where(x >= 0, 1.0, -1.0) * tab[np.minimum(179, np.floor(th * 2).astype(int))]


def t_halo(tab, x, y):
    b = np.maximum(BCORE, np.hypot(x, y))
    return KRM * N0 * b ** -1.5 * fval(tab, x, y)


def t_step(x):
    return 1.65 * np.where(x >= 0, 1.0, -1.0)


def rot(PA, e, n):
    s, c = np.sin(PA * DEG), np.cos(PA * DEG)
    return e * s + n * c, e * c - n * s


# ---------------------------------------------------------------- data
def read_header(f):
    cards = []
    while True:
        blk = f.read(2880)
        if len(blk) < 2880:
            raise SystemExit('short read in header')
        for i in range(0, 2880, 80):
            c = blk[i:i + 80].decode('ascii', errors='replace')
            cards.append(c)
            if c.startswith('END     '):
                return cards


def card(cards, key):
    for c in cards:
        if c[:8].strip() == key:
            return c[10:].split('/')[0].strip().strip("'").strip()
    return None


def table_dtype(h):
    nf = int(card(h, 'TFIELDS'))
    fields, off = {}, 0
    for i in range(1, nf + 1):
        n, t = card(h, 'TTYPE%d' % i), card(h, 'TFORM%d' % i)
        if t.endswith('A'):
            code, size = 'S%d' % int(t[:-1] or 1), int(t[:-1] or 1)
        else:
            code, size = {'D': ('>f8', 8), 'E': ('>f4', 4), 'K': ('>i8', 8), 'J': ('>i4', 4), 'L': ('S1', 1)}[t]
        fields[n] = (code, off)
        off += size
    rowbytes = int(card(h, 'NAXIS1'))
    assert off == rowbytes, (off, rowbytes)
    return np.dtype({'names': list(fields), 'formats': [v[0] for v in fields.values()], 'offsets': [v[1] for v in fields.values()],
                     'itemsize': rowbytes}), rowbytes


def extract():
    t0 = time.time()
    keep = {k: [] for k in ('ra', 'dec', 'rm', 'rm_err', 'source_id')}
    with gzip.open(SPICE, 'rb') as f:
        read_header(f)
        h = read_header(f)
        dt, rowbytes = table_dtype(h)
        nrows, done = int(card(h, 'NAXIS2')), 0
        while done < nrows:
            n = min(100_000, nrows - done)
            buf = f.read(n * rowbytes)
            if len(buf) != n * rowbytes:
                raise SystemExit('short read at row %d' % done)
            a = np.frombuffer(buf, dtype=dt)
            ok = np.ones(n, bool)
            for k in ('snr_flag', 'leakage_flag', 'channel_flag', 'stokesI_fit_flag'):
                ok &= a[k] == b'F'
            rm, e = a['rm'].astype('<f8'), a['rm_err'].astype('<f8')
            ok &= np.isfinite(rm) & np.isfinite(e) & (e > 0)
            for k in ('ra', 'dec'):
                keep[k].append(a[k][ok].astype('<f8'))
            keep['rm'].append(rm[ok]); keep['rm_err'].append(e[ok]); keep['source_id'].append(a['source_id'][ok].copy())
            done += n
    os.makedirs(WORK, exist_ok=True)
    np.savez(os.path.join(WORK, 'spice_rm.npz'), **{k: np.concatenate(v) for k, v in keep.items()})
    print('extract done, %d s' % (time.time() - t0))


def load_sources():
    z = np.load(os.path.join(WORK, 'spice_rm.npz'))
    sid, rm, e = z['source_id'], z['rm'], z['rm_err']
    u, inv = np.unique(sid, return_inverse=True)
    w = 1 / e ** 2
    sw = np.bincount(inv, weights=w)
    first = np.zeros(len(u), dtype=int); first[inv[::-1]] = np.arange(len(sid))[::-1]
    return dict(ra=z['ra'][first], dec=z['dec'][first], rm=np.bincount(inv, weights=w * rm) / sw, e=1 / np.sqrt(sw))


def load_galaxies():
    G, seen = [], set()
    for line in open(os.path.join(HERE, 'galaxies.txt')):
        if line.startswith('#') or not line.strip():
            continue
        p = line.split()
        if p[0] == 'PH':
            name, ra, de, D, i, PA = p[1], *map(float, p[2:7])
            G.append(dict(name=name, ra=ra, de=de, D=D, i=i, PA=PA, src='PHANGS', first=True))
        else:
            name, first = p[0], p[1] == '1'
            ra, de, Dl, Du = float(p[2]), float(p[3]), float(p[4]), p[5]
            G.append(dict(name=name, ra=ra, de=de, D=float(Du) if Du != 'NA' else Dl, i=float(p[7]), PA=float(p[10]),
                          src='LVHIS', first=first, flag_i=p[8], flag_PA=p[9]))
    lv = {g['name'] for g in G if g['src'] == 'LVHIS'}
    alias = {'NGC0628': 'NGC628'}
    for g in G:                               # a galaxy in both lists keeps its LVHIS fit
        g['dup'] = g['src'] == 'PHANGS' and alias.get(g['name'], g['name']) in lv
        g['b'] = float(gal_b(g['ra'], g['de']))
    return G


# ---------------------------------------------------------------- the fit (Entry 61's PT.fit)
class Sky:
    def __init__(self, src):
        self.src = src
        self.v = unit(src['ra'], src['dec'])
        self.tree = cKDTree(self.v)

    def around(self, ra, de, D, rmax_kpc):
        """Sources within rmax of (ra, de) at distance D (Mpc): gnomonic offsets in kpc east and north, and RM."""
        rd = rmax_kpc / (D * 1000.0) / DEG
        idx = np.array(self.tree.query_ball_point(unit(ra, de), chord(rd)), dtype=int)
        if len(idx) == 0:
            return np.zeros(0), np.zeros(0), np.zeros(0)
        a0, d0 = ra * DEG, de * DEG
        a, d = self.src['ra'][idx] * DEG, self.src['dec'][idx] * DEG
        cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(a - a0)
        xi = np.cos(d) * np.sin(a - a0) / cosc
        eta = (np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(a - a0)) / cosc
        sc = D * 1000.0
        return xi * sc, eta * sc, self.src['rm'][idx]


def design(e, n, i, PA, o, tmpl='halo'):
    u, v = e / o['scale'], n / o['scale']
    cols = [u ** a * v ** c for a in range(o['order'] + 1) for c in range(o['order'] + 1 - a)]
    x, y = rot(PA, e, n)
    t = t_step(x) if tmpl == 'step' else t_halo(ftab(i), x, y)
    cols.append(-t)
    return np.stack(cols, axis=1), x, y, t


def fit(e, n, rm, i, PA, o=None, tmpl='halo', inj=0.0):
    o = dict(O, **(o or {}))
    b = np.hypot(e, n)
    m = (b >= o['bin']) & (b < o['bout'])
    e, n, r = e[m], n[m], rm[m].copy()
    nall = int(m.sum())
    A, x, y, t = design(e, n, i, PA, o, tmpl)
    if inj:
        r = r - inj * t
    k = A.shape[1]
    if nall < k + 12:
        return dict(n=0, nAll=nall, A=np.nan, err=np.nan)
    med = np.sort(r)[nall >> 1]
    mad = 1.4826 * np.sort(np.abs(r - med))[nall >> 1]
    keep = np.abs(r - med) < 6 * mad
    for _ in range(8):
        Ak, rk = A[keep], r[keep]
        ATA = Ak.T @ Ak
        sol = np.linalg.solve(ATA, Ak.T @ rk)
        cov = np.linalg.inv(ATA)
        res = r - A @ sol
        rms = np.sqrt(np.sum(res[keep] ** 2) / (keep.sum() - k))
        k2 = keep & (np.abs(res) < o['clip'] * rms)
        if k2.sum() == keep.sum():
            break
        keep = k2
    return dict(n=int(keep.sum()), nAll=nall, A=float(sol[-1]), err=float(rms * np.sqrt(cov[-1, -1])), rms=float(rms))


def cover(e, n, PA, o=None):
    o = dict(O, **(o or {}))
    b = np.hypot(e, n)
    m = (b >= o['bin']) & (b < o['bout'])
    x, y = rot(PA, e[m], n[m])
    q = np.bincount((x > 0).astype(int) + 2 * (y > 0).astype(int), minlength=4)
    return int(m.sum()), float(q.min() / max(1, m.sum()))


def places(sky, g, nref, G_all, o=None, tmpl='halo', rng=None):
    """Blank places: a grid over the catalogue's sky with spacing step x the ring's outer radius, |b| >= 20 deg, away from every
    sample galaxy by 1.5 ring radii, holding at least half the galaxy's sight lines, every sky quadrant >= 10% of them;
    at most 400 places (a fixed-seed random subset), 6 random angles each."""
    o = dict(O, **(o or {}))
    rd = o['bout'] / (g['D'] * 1000.0) / DEG
    step = 0.45 if rd > 3 else (0.7 if rd > 1.5 else 1.2)
    rng = rng or np.random.default_rng(NULL['seed'])
    cand = []
    de = -90.0 + 0.5 * step * rd
    while de <= 49.0:
        c = max(np.cos(de * DEG), 1e-3)
        ras = np.arange(0.0, 360.0, step * rd / c)
        cand += [(ra, de) for ra in ras]
        de += step * rd
    cand = np.array(cand)
    bb = gal_b(cand[:, 0], cand[:, 1])
    ok = np.abs(bb) >= NULL['bmin']
    vc = unit(cand[:, 0], cand[:, 1])
    for h in G_all:
        rh = max(rd, o['bout'] / (h['D'] * 1000.0) / DEG)
        ok &= (vc @ unit(h['ra'], h['de'])) < np.cos(1.5 * rh * DEG)
    cand = cand[ok]
    order = rng.permutation(len(cand))
    out, used = [], 0
    for j in order:
        ra, de = cand[j]
        e, n, rm = sky.around(ra, de, g['D'], o['bout'] * 1.02)
        b = np.hypot(e, n)
        m = (b >= o['bin']) & (b < o['bout'])
        nin = int(m.sum())
        if nin < NULL['frac'] * nref:
            continue
        qd = np.bincount((e[m] > 0).astype(int) + 2 * (n[m] > 0).astype(int), minlength=4)
        if qd.min() < 0.1 * nin:
            continue
        for _ in range(NULL['nang']):
            f = fit(e, n, rm, g['i'], rng.uniform(0, 360), o, tmpl)
            if np.isfinite(f['A']):
                out.append(f['A'])
        used += 1
        if used >= NULL['max_places']:
            break
    return np.array(out), used


def rob(a):
    m = np.median(a)
    return float(1.4826 * np.median(np.abs(a - m)))


def stack_null(R, seed=None, ndraw=None):
    seed = STACK['seed'] if seed is None else seed
    ndraw = ndraw or STACK['ndraw']
    w = np.array([1 / r['rob'] ** 2 for r in R])
    W = w.sum()
    obs = float(np.sum(w * np.array([r['A'] for r in R])) / W)
    g = np.random.default_rng(seed)
    tot, npos = np.zeros(ndraw), np.zeros(ndraw, int)
    for wi, r in zip(w, R):
        d = r['nulls'][g.integers(0, len(r['nulls']), ndraw)]
        tot += wi * d
        npos += d > 0
    s = np.sort(tot / W)
    return dict(n=len(R), A=obs, sd=float(s.std()), mean_null=float(s.mean()),
                p_upper=float((s >= obs).mean()), p_lower=float((s <= obs).mean()),
                lo95=float(obs - s[int(0.975 * ndraw)]), hi95=float(obs - s[int(0.025 * ndraw)]),
                upper95_one_sided=float(obs - s[int(0.05 * ndraw)]), p_all_positive=float((npos == len(R)).mean()),
                names=[r['name'] for r in R])


def analyse_galaxy(sky, g, G_sample, o=None, tmpl='halo', with_null=True, log=print):
    o = dict(O, **(o or {}))
    e, n, rm = sky.around(g['ra'], g['de'], g['D'], o['bout'] * 1.02)
    nin, qmin = cover(e, n, g['PA'], o)
    out = dict(name=g['name'], src=g['src'], D=g['D'], i=g['i'], PA=g['PA'], b=round(g['b'], 1), n=nin, qmin=round(qmin, 3))
    if nin < CUTS['nmin']:
        out['status'] = 'fewer than %d sight lines' % CUTS['nmin']
        return out
    f = fit(e, n, rm, g['i'], g['PA'], o, tmpl)
    out.update(A=f['A'], formal=f['err'], n_fit=f['n'])
    if qmin < CUTS['qmin']:
        out['status'] = 'poorly covered'
    else:
        out['status'] = 'fitted'
    if with_null and np.isfinite(f['A']):
        t0 = time.time()
        nulls, used = places(sky, g, nin, G_sample, o, tmpl, np.random.default_rng(NULL['seed']))
        out.update(places=used, n_null=int(len(nulls)), sd=float(nulls.std()) if len(nulls) else np.nan,
                   rob=rob(nulls) if len(nulls) > 10 else np.nan,
                   pct=float((nulls < f['A']).mean()) if len(nulls) else np.nan)
        out['nulls'] = nulls
        log('%-12s D %5.2f i %4.1f b %+5.1f n %4d q %.2f A %+.3f formal %.3f rob %.3f places %d (%.0f s)' %
            (g['name'], g['D'], g['i'], g['b'], nin, qmin, f['A'], f['err'], out['rob'], used, time.time() - t0))
    return out


def passes(g):
    return (g['i'] >= CUTS['imin'] and abs(g['b']) >= CUTS['bmin'] and CUTS['Dmin'] <= g['D'] <= CUTS['Dmax'])


# ---------------------------------------------------------------- synthetic test (no data read)
def synthetic():
    g0 = np.random.default_rng(7)
    nsrc = 150000
    u = g0.uniform(np.sin(-90 * DEG), np.sin(49 * DEG), nsrc)
    dec, ra = np.arcsin(u) / DEG, g0.uniform(0, 360, nsrc)
    bb = gal_b(ra, dec)
    l = ra  # a smooth pattern is enough here
    rm = 20 * np.sin(l * DEG * 2 + 0.3) * np.cos(dec * DEG) + 10 / (np.abs(np.sin(bb * DEG)) + 0.3) + g0.normal(0, 8, nsrc)
    src = dict(ra=ra, dec=dec, rm=rm, e=np.full(nsrc, 2.0))
    G = [dict(name='S1', ra=13.7, de=-37.7, D=2.15, i=40, PA=290, src='syn', b=float(gal_b(13.7, -37.7))),
         dict(name='S2', ra=44.2, de=-54.6, D=5.0, i=75, PA=218, src='syn', b=float(gal_b(44.2, -54.6)))]
    out = {}
    for inj in (0.0, 0.5):
        rm2 = rm.copy()
        if inj:
            for g in G:          # add a current along the spin of inj x 1e18 A around each synthetic galaxy
                sky0 = Sky(src)
                e, n, _ = sky0.around(g['ra'], g['de'], g['D'], 1000.0)
                rd = 1000.0 / (g['D'] * 1000) / DEG
                idx = np.array(sky0.tree.query_ball_point(unit(g['ra'], g['de']), chord(rd)), dtype=int)
                x, y = rot(g['PA'], e, n)
                rm2[idx] += -inj * t_halo(ftab(g['i']), x, y)
        sky = Sky(dict(src, rm=rm2))
        for g in G:
            r = analyse_galaxy(sky, g, G)
            out['%s_inj%.1f' % (g['name'], inj)] = dict(A=round(r['A'], 3), rob=round(r['rob'], 3), n=r['n'], places=r['places'])
    print(json.dumps(out, indent=1))


def run(only=None):
    t0 = time.time()
    src = load_sources()
    sky = Sky(src)
    G = load_galaxies()
    use = [g for g in G if g['first'] and not g['dup']]
    sample = [g for g in use if passes(g)]
    res = dict(generated=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), fit=O, cuts=CUTS, null=NULL, stack=STACK,
               sources=int(len(src['rm'])), galaxies_listed=len(use), galaxies_passing_cuts=len(sample))
    print('sources %d; galaxies %d, passing the cuts %d' % (len(src['rm']), len(use), len(sample)), flush=True)
    R = []
    for g in use:
        r = analyse_galaxy(sky, g, use, log=lambda s: print(s, flush=True))
        r['passes_cuts'] = passes(g)
        R.append(r)
    good = [r for r in R if r['status'] == 'fitted' and r['passes_cuts'] and np.isfinite(r.get('rob', np.nan))]
    res['primary'] = stack_null(good) if good else None
    new = [r for r in good if r['name'] not in ENTRY61]
    res['without_entry61'] = stack_null(new) if new else None
    allfit = [r for r in R if r['status'] == 'fitted' and np.isfinite(r.get('rob', np.nan))]
    res['all_well_covered'] = stack_null(allfit) if allfit else None
    if len(good) < 5:
        res['decision'] = 'not tested (fewer than 5 galaxies pass every cut)'
    elif res['primary']:
        P = res['primary']
        if P['p_upper'] <= 0.00135 and res['without_entry61'] and res['without_entry61']['p_upper'] <= 0.05:
            out = 'current along the spin (stack p <= 0.00135, and p <= 0.05 without Entry 61 galaxies)'
        elif P['p_lower'] <= 0.00135:
            out = 'current against the spin (stack p <= 0.00135, lower tail)'
        elif P['upper95_one_sided'] < LEAN:
            out = 'the lean of Entries 60-61 is not confirmed: one-sided 95%% upper bound %.3f < +0.20' % P['upper95_one_sided']
        else:
            out = 'unresolved'
        res['decision'] = out
        print('decision:', out, json.dumps({k: P[k] for k in ('n', 'A', 'sd', 'p_upper', 'upper95_one_sided')}), flush=True)
    # reported only: other settings, as Entry 61 reported them, weighted with the primary blank-sky weights (no new nulls)
    byname = {g['name']: g for g in use}
    second = {g['name']: g for g in G if not g['first']}
    other = {}
    if good:
        w = np.array([1 / r['rob'] ** 2 for r in good])
        for label, o2, alt in (('order1', {'order': 1}, None), ('order3', {'order': 3}, None), ('from50', {'bin': 50.0}, None),
                               ('to300', {'bout': 300.0}, None), ('outer_100_200', {'bin': 100.0}, None),
                               ('turned90', None, 'turn'), ('second_fits', None, 'second'), ('step_template', None, 'step')):
            vals = []
            for r in good:
                g = dict(byname[r['name']])
                tm = 'halo'
                if alt == 'turn':
                    g['PA'] = (g['PA'] + 90) % 360
                elif alt == 'second' and r['name'] in second:
                    g = dict(second[r['name']], b=g['b'])
                elif alt == 'step':
                    tm = 'step'
                oo = dict(O, **(o2 or {}))
                e, n, rm = sky.around(g['ra'], g['de'], g['D'], oo['bout'] * 1.02)
                vals.append(fit(e, n, rm, g['i'], g['PA'], oo, tm)['A'])
            vals = np.array(vals)
            ok = np.isfinite(vals)
            other[label] = dict(A=float(np.sum(w[ok] * vals[ok]) / w[ok].sum()), n=int(ok.sum()),
                                n_positive=int((vals[ok] > 0).sum()))
    res['other_settings_weighted_as_primary'] = other
    res['galaxies'] = [{k: v for k, v in r.items() if k != 'nulls'} for r in R]
    for r in res['galaxies']:
        for k, v in list(r.items()):
            if isinstance(v, float):
                r[k] = round(v, 4)
    res['seconds'] = round(time.time() - t0)
    json.dump(res, open(os.path.join(HERE, 'results.json'), 'w'), indent=1, default=float)
    print('written results.json, %d s' % res['seconds'], flush=True)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    {'extract': extract, 'synthetic': synthetic, 'run': run}.get(cmd, lambda: print(__doc__))()
