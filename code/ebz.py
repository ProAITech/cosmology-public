"""Entry 154: does a star's metal content change its brightness at fixed mass?
Eker et al. (2018) eclipsing binaries (masses, radii, Teff) with literature or survey [Fe/H]:
A = APOGEE DR17, P = PASTEL (median of entries), L = LAMOST DR7 LRS, G = Geneva-Copenhagen (Casagrande 2011); mean of available.
Prediction: Tout et al. (1996) ZAMS L, R, Teff as functions of mass and Z (Z = 0.02 x 10^[Fe/H]).
Fusion: residuals from the solar-Z ZAMS follow the Z-shift (slope 1). External supply: no metal dependence (slope 0)."""
import json, sys, numpy as np
exec(open(sys.argv[2] if len(sys.argv) > 2 else 'tout.py').read().split("if __name__")[0])
rows = [r.split(',') for r in open(sys.argv[1] if len(sys.argv) > 1 else 'eb_feh.txt').read().strip().split(';')]
stars = []
for r in rows:
    name = r[0]; fe = float(r[9]); src = r[10]
    for k in (0, 1):
        M, Rr, T, eT = (float(x) for x in r[1+4*k:5+4*k])
        stars.append(dict(sys=name, M=M, R=Rr, T=T, eT=eT, fe=fe, src=src, Tp=float(r[3])))
for s in stars:
    Z = 0.02*10**s['fe']; M = s['M']
    lL = 2*np.log10(s['R'])+4*np.log10(s['T']/5772)
    s['rL'] = lL-np.log10(L(M)); s['rR'] = np.log10(s['R']/R(M)); s['rT'] = np.log10(s['T']/Teff(M))
    s['pL'] = np.log10(L(M, Z)/L(M)); s['pR'] = np.log10(R(M, Z)/R(M)); s['pT'] = np.log10(Teff(M, Z)/Teff(M))
    s['e'] = s['rR']-s['pR']
def sel(s, Tmax=7500, Mlo=0.5, Mhi=2.5, emax=0.3):
    return s['T'] < Tmax and s['Tp'] < Tmax and Mlo < s['M'] < Mhi and s['e'] < emax and -0.8 < s['fe'] <= FEMAX
FEMAX = 0.176   # Tout et al. fits hold for Z = 0.0001-0.03, i.e. [Fe/H] <= +0.176 with Z_sun = 0.02; beyond it the polynomial turns back
def fit(S, y, p):
    X = np.array([[s[p], s['e'], 1.0] for s in S]); Y = np.array([s[y] for s in S])
    c, *_ = np.linalg.lstsq(X, Y, rcond=None); res = Y-X@c
    return c, res
def boot(S, y, p, n=20000, seed=154):
    names = sorted(set(s['sys'] for s in S)); by = {nm: [s for s in S if s['sys'] == nm] for nm in names}
    rng = np.random.default_rng(seed); A = []
    for _ in range(n):
        pick = rng.choice(len(names), len(names)); SS = [s for i in pick for s in by[names[i]]]
        try: A.append(fit(SS, y, p)[0][0])
        except Exception: pass
    A = np.array(A); return float(np.std(A)), float(np.percentile(A, 2.5)), float(np.percentile(A, 97.5)), float((A <= 0).mean())
out = {}
S = [s for s in stars if sel(s)]
out['n_stars'] = len(S); out['n_sys'] = len(set(s['sys'] for s in S))
out['fe_range'] = [min(s['fe'] for s in S), max(s['fe'] for s in S)]; out['fe_sd'] = float(np.std([s['fe'] for s in S]))
out['M_range'] = [min(s['M'] for s in S), max(s['M'] for s in S)]
for y, p in (('rT', 'pT'), ('rL', 'pL')):
    c, res = fit(S, y, p); e, lo, hi, pneg = boot(S, y, p)
    out[y] = dict(a=round(float(c[0]), 3), b_evol=round(float(c[1]), 3), c=round(float(c[2]), 4), ea=round(e, 3), ci95=[round(lo, 3), round(hi, 3)],
                  frac_boot_le0=pneg, rms=round(float(np.sqrt(np.mean(res**2))), 4),
                  pred_range=[round(min(s[p] for s in S), 4), round(max(s[p] for s in S), 4)])
# robustness: variations of the cuts
var = {}
for lab, kw in (('Tmax7000', dict(Tmax=7000)), ('Tmax8000', dict(Tmax=8000)), ('e<0.2', dict(emax=0.2)), ('e<0.4', dict(emax=0.4)), ('M0.7-1.8', dict(Mlo=0.7, Mhi=1.8))):
    SS = [s for s in stars if sel(s, **kw)]
    c, _ = fit(SS, 'rT', 'pT'); e, lo, hi, pneg = boot(SS, 'rT', 'pT', n=4000)
    var[lab] = dict(n=len(SS), a=round(float(c[0]), 2), ea=round(e, 2))
# leave one source out
for src in 'GLPA':
    SS = [s for s in S if src not in s['src'] or len(s['src']) > 1]
    if len(set(s['sys'] for s in SS)) < 10: continue
    c, _ = fit(SS, 'rT', 'pT'); e, *_ = boot(SS, 'rT', 'pT', n=4000)
    var['without_'+src] = dict(n=len(SS), a=round(float(c[0]), 2), ea=round(e, 2))
out['variants'] = var
# metal-rich vs metal-poor halves, Teff residual after removing the fitted evolution term only
c, _ = fit(S, 'rT', 'pT'); bev = c[1]
adj = np.array([s['rT']-bev*s['e'] for s in S]); fe = np.array([s['fe'] for s in S]); pT = np.array([s['pT'] for s in S])
lo = fe < np.median(fe)
out['halves'] = dict(fe_split=float(np.median(fe)), poor_mean_fe=round(float(fe[lo].mean()), 3), rich_mean_fe=round(float(fe[~lo].mean()), 3),
                     dT_obs=round(float(adj[lo].mean()-adj[~lo].mean()), 4), dT_pred=round(float(pT[lo].mean()-pT[~lo].mean()), 4))
print(json.dumps(out, indent=1))
fig = dict(pts=[dict(sys=s['sys'], M=s['M'], fe=s['fe'], src=s['src'], pT=round(float(s['pT']), 4), adjT=round(float(s['rT']-bev*s['e']-c[2]), 4), eT=round(float(s['eT']/s['T']/np.log(10)), 4)) for s in S])
json.dump(dict(summary=out, fig=fig), open('ebz.json', 'w'), indent=1)
# binned by [Fe/H], system-averaged, for the figure
sysd = {}
for s in S: sysd.setdefault(s['sys'], []).append(s)
SA = [dict(sys=k, fe=v[0]['fe'], adj=float(np.mean([x['rT']-bev*x['e']-c[2] for x in v])), pT=float(np.mean([x['pT'] for x in v])), M=float(np.mean([x['M'] for x in v]))) for k, v in sysd.items()]
bins = [(-0.6, -0.15), (-0.15, -0.05), (-0.05, 0.05), (0.05, 0.3)]
B = []
for lo_, hi_ in bins:
    v = [x for x in SA if lo_ <= x['fe'] < hi_]
    if len(v) < 2: continue
    a_ = np.array([x['adj'] for x in v])
    B.append(dict(lo=lo_, hi=hi_, n=len(v), fe=round(float(np.mean([x['fe'] for x in v])), 3), adj=round(float(a_.mean()), 4), e=round(float(a_.std(ddof=1)/np.sqrt(len(v))), 4), pT=round(float(np.mean([x['pT'] for x in v])), 4)))
Mmed = float(np.median([s['M'] for s in S]))
curves = {str(m): [[round(f, 2), round(float(np.log10(Teff(m, 0.02*10**f)/Teff(m))), 4)] for f in np.arange(-2.2, 0.51, 0.05)] for m in (0.8, 1.0, 1.3)}
gc = {f'{m}_{fe}': dict(T_solar=int(Teff(m)), T=int(Teff(m, 0.02*10**fe)), L_ratio=round(float(L(m, 0.02*10**fe)/L(m)), 2)) for m in (0.8, 1.0) for fe in (-1.0, -1.5, -2.0)}
print(json.dumps(dict(bins=B, Mmed=Mmed, gc=gc), indent=1))
d = json.load(open('ebz.json')); d['bins'] = B; d['Mmed'] = Mmed; d['curves'] = curves; d['gc'] = gc
d['gspphot_vs_lit'] = dict(n=35, r=0.29, offset=-0.29, scatter=0.25)
d['sys_pts'] = [dict(sys=x['sys'], fe=x['fe'], M=round(x['M'], 3), pT=round(x['pT'], 4), adj=round(x['adj'], 4)) for x in SA]
json.dump(d, open('ebz.json', 'w'), indent=1)
