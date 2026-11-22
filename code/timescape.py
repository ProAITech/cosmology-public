"""Workbench Entry 128: dark energy as a clock effect? The timescape cosmology against DESI's ruler and the supernovae.

Timescape (Wiltshire 2007-2009) keeps expansion but drops dark energy: in a universe of empty voids and dense walls, clocks in walls
(where we and all galaxies sit) run slower than the volume average, and the apparent acceleration is that difference. Its tracker solution
has one shape parameter, the present void fraction f_v0 (Wiltshire 2009, eqs. 38-41, B2, B8):
   f_v(t) = 3 f_v0 Hb0 t / [3 f_v0 Hb0 t + (1 - f_v0)(2 + f_v0)]
   1 + z  = (2 + f_v) f_v^(1/3) / [3 f_v0^(1/3) Hb0 t]
   H(t)   = (4 f_v^2 + f_v + 4) / (6 t)               dressed Hubble parameter
   d_A    = t^(2/3) Int_t^t0 2 dt' / [(2 + f_v) t'^(2/3)],   D = (1 + z) d_A,   d_L = (1 + z) D
Seifert et al. (2025) report very strong Bayesian evidence for timescape over LCDM in Pantheon+ supernovae. This script tests the same
model's shape against DESI DR2's 13 ruler numbers (normalisation c/(Hb0 r_d) free, so the sound horizon and H0 do not enter), against
the radial-to-transverse ratio D_M/D_H at each redshift (free of any normalisation), and against Pantheon+ with diagonal errors.
Caveat (Wiltshire 2009): DESI converts redshifts to distances with a fiducial FLRW model; for timescape that conversion is approximate."""
import os, sys, json, math
import numpy as np
from scipy import integrate, optimize
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = HERE if os.path.exists(os.path.join(HERE, 'cosmo.py')) else os.path.abspath(os.path.join(HERE, '..', 'artifact-files', '4379c2ac-1195-4a41-978e-eb44b945a265', 'code'))
sys.path.insert(0, CODE)
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()): import cosmo as c
ZV, YB, CBI, CB, ZB = c.ZV, c.YB, c.CBI, c.CB, c.ZB
def ts_funcs(fv0):
    """returns functions of z: D(z) and D_H(z) in units of c/Hb0, and the dressed H0/Hb0"""
    t0 = (2 + fv0)/3; cc = (1 - fv0)*(2 + fv0)
    t = np.geomspace(1e-5*t0, t0, 40000)
    fv = 3*fv0*t/(3*fv0*t + cc)
    z = (2 + fv)*fv**(1/3)/(3*fv0**(1/3)*t) - 1
    integ = 2/((2 + fv)*t**(2/3))
    I = integrate.cumulative_trapezoid(integ, t, initial=0); I = I[-1] - I     # int_t^t0
    dA = t**(2/3)*I; D = (1 + z)*dA
    H = (4*fv**2 + fv + 4)/(6*t)
    o = np.argsort(z); zs, Ds, Hs = z[o], D[o], H[o]
    H0d = (4*fv0**2 + fv0 + 4)/(2*(2 + fv0))
    return (lambda zz: np.interp(zz, zs, Ds)), (lambda zz: 1/np.interp(zz, zs, Hs)), H0d
def lcdm_funcs(om):
    zg = np.linspace(0, 3, 6001); E = np.sqrt(om*(1+zg)**3 + 1 - om)
    Dc = integrate.cumulative_trapezoid(1/E, zg, initial=0)
    return (lambda zz: np.interp(zz, zg, Dc)), (lambda zz: 1/np.sqrt(om*(1+zz)**3 + 1 - om)), 1.0
def bao_vec(Df, DHf):
    dv = (ZV*Df(ZV)**2*DHf(ZV))**(1/3)
    return np.r_[dv, Df(ZB), DHf(ZB)]
def chi_bao(Df, DHf):
    v = bao_vec(Df, DHf); A = (v @ CBI @ YB)/(v @ CBI @ v); r = YB - A*v
    return float(r @ CBI @ r), A
# radial-to-transverse ratio at each of the six two-dimensional bins
iM = np.arange(1, 7); iH = np.arange(7, 13)
FAP = YB[iM]/YB[iH]
eF = np.array([FAP[k]*math.sqrt(CB[iM[k], iM[k]]/YB[iM[k]]**2 + CB[iH[k], iH[k]]/YB[iH[k]]**2 - 2*CB[iM[k], iH[k]]/(YB[iM[k]]*YB[iH[k]])) for k in range(6)])
def chi_fap(Df, DHf):
    m = Df(ZB)/DHf(ZB); return float(np.sum(((FAP - m)/eF)**2)), m
# supernovae: Pantheon+ rows, repeats averaged, diagonal errors
d = np.loadtxt(os.path.join(CODE, 'pantheon_plus_rows.txt')); zr, mr, er = d.T
zu = np.unique(zr); Z, M, E = [], [], []
for z0 in zu:
    q = zr == z0; w = 1/er[q]**2; Z.append(z0); M.append(np.sum(w*mr[q])/w.sum()); E.append(er[q].min())
Z, M, E = map(np.array, (Z, M, E)); W = 1/E**2
def chi_sn(Df, sel=None):
    s = np.ones(len(Z), bool) if sel is None else sel
    r = M[s] - 5*np.log10((1 + Z[s])*Df(Z[s])); a = np.sum(W[s]*r)/W[s].sum(); return float(np.sum(W[s]*(r - a)**2))
out = {}
def scan(maker, grid, f):
    vals = [f(*maker(p)[:2]) for p in grid]; vals = [v[0] if isinstance(v, tuple) else v for v in vals]
    i = int(np.argmin(vals)); p0 = grid[i]
    r = optimize.minimize_scalar(lambda p: (lambda v: v[0] if isinstance(v, tuple) else v)(f(*maker(p)[:2])), bounds=(grid[max(i-2, 0)], grid[min(i+2, len(grid)-1)]), method='bounded')
    c0 = r.fun
    # 1-sigma range from the grid
    ok = [g for g, v in zip(grid, vals) if v <= c0 + 1]
    return float(r.x), float(c0), (min(ok), max(ok)) if ok else (None, None)
FV = np.linspace(0.30, 0.95, 131); OM = np.linspace(0.15, 0.50, 141)
for lab, f in (('BAO (13 numbers, scale free)', lambda D, H: chi_bao(D, H)), ('D_M/D_H ratios (6)', lambda D, H: chi_fap(D, H)),
               ('Pantheon+ (diagonal)', lambda D, H: chi_sn(D)), ('Pantheon+ beyond z = 0.075', lambda D, H: chi_sn(D, Z > 0.075))):
    ts = scan(ts_funcs, FV, f); lc = scan(lcdm_funcs, OM, f)
    out[lab] = dict(ts=ts, lcdm=lc)
    print(f"{lab:32s} timescape: f_v0 = {ts[0]:.3f} ({ts[2][0]:.3f}-{ts[2][1]:.3f}), chi2 {ts[1]:7.2f} | LCDM: Om = {lc[0]:.3f} ({lc[2][0]:.3f}-{lc[2][1]:.3f}), chi2 {lc[1]:7.2f} | LCDM minus timescape {lc[1]-ts[1]:+.2f}")
# joint BAO + SN with one shape parameter each
def joint(maker, p): D, H, _ = maker(p); return chi_bao(D, H)[0] + chi_sn(D)
tsj = optimize.minimize_scalar(lambda p: joint(ts_funcs, p), bounds=(0.3, 0.95), method='bounded')
lcj = optimize.minimize_scalar(lambda p: joint(lcdm_funcs, p), bounds=(0.15, 0.5), method='bounded')
out['joint'] = dict(ts=(float(tsj.x), float(tsj.fun)), lcdm=(float(lcj.x), float(lcj.fun)))
print(f"joint BAO + supernovae: timescape f_v0 = {tsj.x:.3f}, chi2 {tsj.fun:.2f}; LCDM Om = {lcj.x:.3f}, chi2 {lcj.fun:.2f}; LCDM minus timescape {lcj.fun - tsj.fun:+.2f}")
for p, lab in ((tsj.x, 'ts'),):
    D, H, h0d = ts_funcs(p); cb, A = chi_bao(D, H); print(f"   at the joint f_v0: BAO chi2 {cb:.2f}, SN chi2 {chi_sn(D):.2f}; dressed/bare H0 = {h0d:.3f}")
    out['joint_parts_ts'] = (cb, chi_sn(D))
D, H, _ = lcdm_funcs(lcj.x); out['joint_parts_lcdm'] = (chi_bao(D, H)[0], chi_sn(D)); print(f"   LCDM parts: BAO {out['joint_parts_lcdm'][0]:.2f}, SN {out['joint_parts_lcdm'][1]:.2f}")
# the ratio curves for the figure, and the published best fit f_v0 = 0.762 (Wiltshire 2009) and Seifert et al.
zz = np.linspace(0.05, 2.6, 120)
curves = {}
for lab, (D, H, _) in (('timescape BAO fit', ts_funcs(out['BAO (13 numbers, scale free)']['ts'][0])), ('timescape SN fit', ts_funcs(out['Pantheon+ (diagonal)']['ts'][0])),
                       ('LCDM BAO fit', lcdm_funcs(out['BAO (13 numbers, scale free)']['lcdm'][0]))):
    curves[lab] = (D(zz)/H(zz)/zz).round(5).tolist()          # D_M/(z D_H): flattens the trend
out['curves'] = dict(z=zz.round(3).tolist(), **curves)
out['fap'] = dict(z=ZB.tolist(), y=(FAP/ZB).tolist(), e=(eF/ZB).tolist())
# dressed vs bare for the BAO-preferred f_v0
for fv in (out['BAO (13 numbers, scale free)']['ts'][0], out['Pantheon+ (diagonal)']['ts'][0], 0.762):
    D, H, h0d = ts_funcs(fv); print(f"   f_v0 = {fv:.3f}: dressed H0 / bare = {h0d:.3f}; dressed matter density {0.5*(1-fv)*(2+fv):.3f} (bare, tracker)")
# DESI DR2 Lyman-alpha Alcock-Paczynski (2026): D_M/r_d = 39.32 +- 0.33, D_H/r_d = 8.600 +- 0.066 at z = 2.33; AP ratio to about 1%
FLy = 39.32/8.600; eLy = 0.01*FLy
out['lya'] = (FLy, eLy)
zf = np.linspace(0.3, 2.5, 111); Dl, Hl, _ = lcdm_funcs(out['BAO (13 numbers, scale free)']['lcdm'][0]); base = Dl(zf)/Hl(zf)
rel = {}
for lab, fv in (('ts_bao', out['BAO (13 numbers, scale free)']['ts'][0]), ('ts_sn', out['Pantheon+ (diagonal)']['ts'][0]), ('ts_762', 0.762)):
    D, H, _ = ts_funcs(fv); rel[lab] = (D(zf)/H(zf)/base).round(5).tolist()
    m = D(2.33)/H(2.33); out[f'lya_{lab}'] = (float(m), float((FLy - m)/eLy))
    print(f"   Lya AP at z = 2.33: timescape f_v0 = {fv:.3f} gives {m:.3f}, data {FLy:.3f} ± {eLy:.3f}: {(FLy-m)/eLy:+.1f} sigma")
bz = Dl(ZB)/Hl(ZB); bl = Dl(2.33)/Hl(2.33)
print(f"   Lya AP: LCDM (BAO fit) gives {bl:.3f}: {(FLy-bl)/eLy:+.1f} sigma")
out['rel'] = dict(z=zf.round(3).tolist(), **rel, data=[[float(z0), float(FAP[k]/bz[k]), float(eF[k]/bz[k])] for k, z0 in enumerate(ZB)], lya=[2.33, FLy/bl, eLy/bl])
from scipy.stats import chi2 as C2
for lab in ('BAO (13 numbers, scale free)', 'D_M/D_H ratios (6)'):
    dof = 11 if lab.startswith('BAO') else 5
    out[lab]['p_ts'] = float(C2.sf(out[lab]['ts'][1], dof)); out[lab]['p_lcdm'] = float(C2.sf(out[lab]['lcdm'][1], dof))
    print(f"   {lab}: p(timescape) = {out[lab]['p_ts']:.1e}, p(LCDM) = {out[lab]['p_lcdm']:.2f}")
json.dump(out, open(os.path.join(HERE, 'timescape.json'), 'w'), indent=1)
