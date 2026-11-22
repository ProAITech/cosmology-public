"""Entry 158: metals at fixed mass, with masses from stellar oscillations instead of orbits.
seis.csv: per star 'numax*10 (muHz), dnu*100 (muHz), Teff (K), [Fe/H]*100, source' with L = Lund et al. 2017 LEGACY (spectroscopic [Fe/H]),
S = Serenelli et al. 2017 APOKASC dwarfs and subgiants (APOGEE [M/H]). Six rows of the Serenelli tables carry identical placeholder values
(numax 690.176, dnu 39.276, Teff 6188, [M/H] 0.15) and are dropped.
Masses and radii from the scaling relations (numax_sun = 3090, dnu_sun = 135.1 muHz, Teff_sun = 5772 K); no stellar-evolution models.
Test as in Entry 154: r_T = log T/T0(M) = a p_T + b e + c, where p_T is the zero-age shift metals predict (Tout et al. 1996) and
e = log R/R0(M) - p_R the expansion with age. Fusion: a = 1. Outside supply: a = 0."""
import json, sys, numpy as np
exec(open(sys.argv[2] if len(sys.argv) > 2 else 'tout.py').read().split("if __name__")[0])
r = [x.split(',') for x in open(sys.argv[1] if len(sys.argv) > 1 else 'seis.csv').read().strip().split(';')]
nm = np.array([int(x[0]) for x in r])/10; dn = np.array([int(x[1]) for x in r])/100; T = np.array([int(x[2]) for x in r]).astype(float)
fe = np.array([int(x[3]) for x in r])/100; src = np.array([x[4] for x in r])
ph = (np.abs(nm-690.2) < 0.1) & (np.abs(dn-39.28) < 0.01)
out = {'n_rows': len(r), 'placeholders_dropped': int(ph.sum())}
M = (nm/3090)**3*(dn/135.1)**-4*(T/5772)**1.5
Rr = (nm/3090)*(dn/135.1)**-2*(T/5772)**0.5
Z = 0.02*10**fe
rT = np.log10(T/Teff(M)); rR = np.log10(Rr/R(M)); pT = np.log10(Teff(M, Z)/Teff(M)); pR = np.log10(R(M, Z)/R(M)); e = rR-pR
rL = np.log10(Rr**2*(T/5772)**4/L(M)); pL = np.log10(L(M, Z)/L(M))
def sel(Tmax=7000, Mlo=0.7, Mhi=1.8, emax=0.3, femax=1.0, which=None):
    s = ~ph & (T < Tmax) & (M > Mlo) & (M < Mhi) & (e < emax) & (fe >= -femax) & (fe <= 0.176)   # Tout fits valid to Z = 0.03
    if which: s &= (src == which)
    return s
def fit(s, y=rT, p=pT):
    X = np.vstack([p[s], e[s], np.ones(s.sum())]).T; c, *_ = np.linalg.lstsq(X, y[s], rcond=None); return c, y[s]-X@c
def boot(s, y=rT, p=pT, n=10000, seed=158):
    idx = np.where(s)[0]; rng = np.random.default_rng(seed); A = []
    for _ in range(n):
        b = rng.choice(idx, len(idx)); m = np.zeros(len(r), bool); X = np.vstack([p[b], e[b], np.ones(len(b))]).T
        c, *_ = np.linalg.lstsq(X, y[b], rcond=None); A.append(c[0])
    A = np.array(A); return float(A.std()), float(np.percentile(A, 2.5)), float(np.percentile(A, 97.5)), float((A <= 0).mean())
s = sel()
out['n'] = int(s.sum()); out['n_L'] = int((s & (src == 'L')).sum()); out['n_S'] = int((s & (src == 'S')).sum())
out['fe_range'] = [float(fe[s].min()), float(fe[s].max())]; out['fe_5_95'] = [float(np.percentile(fe[s], 5)), float(np.percentile(fe[s], 95))]
out['M_median'] = round(float(np.median(M[s])), 2); out['M_16_84'] = [round(float(np.percentile(M[s], 16)), 2), round(float(np.percentile(M[s], 84)), 2)]
for lab, y, p in (('T', rT, pT), ('L', rL, pL)):
    c, res = fit(s, y, p); sd, lo, hi, fz = boot(s, y, p)
    out['resp_'+lab] = dict(a=round(float(c[0]), 3), ea=round(sd, 3), ci95=[round(lo, 2), round(hi, 2)], b_evol=round(float(c[1]), 3), c=round(float(c[2]), 4), frac_le0=fz, z_from0=round(float(c[0])/sd, 1), z_from1=round((float(c[0])-1)/sd, 1), rms=round(float(res.std()), 4))
var = {}
for lab, kw in (('emax0.2', dict(emax=0.2)), ('emax0.4', dict(emax=0.4)), ('Tmax6500', dict(Tmax=6500)), ('M0.8-1.5', dict(Mlo=0.8, Mhi=1.5)), ('fe>-0.6', dict(femax=0.6)), ('LEGACY only', dict(which='L')), ('APOKASC only', dict(which='S'))):
    ss = sel(**kw); c, _ = fit(ss); sd, *_ = boot(ss, n=3000)
    var[lab] = dict(n=int(ss.sum()), a=round(float(c[0]), 2), ea=round(sd, 2))
out['variants'] = var
c, _ = fit(s); bev = c[1]; adj = rT-bev*e-c[2]
q = np.percentile(fe[s], [0, 25, 50, 75, 100]); bins = []
for i in range(4):
    m = s & (fe >= q[i]) & ((fe < q[i+1]) if i < 3 else (fe <= q[i+1]))
    bins.append(dict(lo=round(float(q[i]), 2), hi=round(float(q[i+1]), 2), n=int(m.sum()), fe=round(float(fe[m].mean()), 3), adj=round(float(adj[m].mean()), 4), e=round(float(adj[m].std(ddof=1)/np.sqrt(m.sum())), 4), pT=round(float(pT[m].mean()), 4)))
out['quartiles'] = bins
lo = s & (fe < np.median(fe[s])); hi = s & (fe >= np.median(fe[s]))
out['halves'] = dict(poor_fe=round(float(fe[lo].mean()), 3), rich_fe=round(float(fe[hi].mean()), 3), dT_obs=round(float(adj[lo].mean()-adj[hi].mean()), 4), dT_pred=round(float(pT[lo].mean()-pT[hi].mean()), 4),
                     err=round(float(np.hypot(adj[lo].std(ddof=1)/np.sqrt(lo.sum()), adj[hi].std(ddof=1)/np.sqrt(hi.sum()))), 4))
# a mass bias with metallicity that would fake the whole response: d log M / d[Fe/H]
dlT_dlM = float(np.median((np.log10(Teff(M[s]*1.01))-np.log10(Teff(M[s])))/np.log10(1.01)))
slope_pT_fe = float(np.polyfit(fe[s], pT[s], 1)[0])
out['mass_bias_to_fake'] = dict(dlogT0_dlogM=round(dlT_dlM, 2), dpT_dfe=round(slope_pT_fe, 3), dlogM_dfe_needed=round(-slope_pT_fe*out['resp_T']['a']/dlT_dlM, 3))
print(json.dumps(out, indent=1))
pts = [[round(float(fe[i]), 2), round(float(pT[i]), 4), round(float(adj[i]), 4), round(float(M[i]), 2), str(src[i])] for i in np.where(s)[0]]
json.dump(dict(summary=out, pts=pts), open('seisz.json', 'w'))
