"""Entry 150: our own gravitational-redshift stack from public SDSS groups (Tempel et al. 2014, VizieR J/A+A/566/A1),
and the predictions it should be compared with. The stacking was run in the browser (grz_stack.js); its numbers are in grz_tempel_raw.json.
Predictions: Entry 122's NFW model (isotropic orbits, galaxies trace mass), here at the sample's median NFW mass and redshift,
bins 0-1.1, 1.1-2.1, 2.1-4.5, 4.5-6 Mpc weighted 1:1:3:3. Kinematic terms: transverse Doppler + light cone +2.5 sigma_los^2/c,
brightest-galaxy terms +2.2 km/s scaled as M^(2/3), surface-brightness selection -K sigma_los^2/c with K = 10 (flux-limited)
or 0 (satellites well above the flux limit). Plasma-only clusters: gravitational term x 0.157."""
import json, math
import numpy as np
G = 4.30091e-9; C = 299792.458; h = 0.677; Om = 0.31; FB = 0.157
R = json.load(open('grz_tempel_raw.json'))
z = R['zmed']; E2 = Om*(1+z)**3 + 1 - Om; rho_c = 2.775e11*h*h*E2
EDG = [(0, 1.1), (1.1, 2.1), (2.1, 4.5), (4.5, 6.0)]; WTS = np.array([1, 1, 3, 3.])
def nfw(M200, c=5.0, rmax=8.0):
    R200 = (3*M200/(4*math.pi*200*rho_c))**(1/3); rs = R200/c; m = math.log(1+c) - c/(1+c); A = G*M200/m
    r = np.geomspace(1e-3, rmax, 4000); x = r/rs; rho = 1/(x*(1+x)**2)
    Mr = A/G*(np.log(1+x) - x/(1+x)); Phi = -A*np.log(1+x)/r; Phi0 = -A/rs
    integrand = rho*G*Mr/r**2; cum = np.concatenate([[0], np.cumsum(0.5*(integrand[1:]+integrand[:-1])*np.diff(r))])
    return dict(R200=R200, r=r, rho=rho, Phi=Phi, Phi0=Phi0, sig2=(cum[-1]-cum)/rho)
def project(H, Rp):
    zz = np.linspace(0, math.sqrt(max(6.0**2-Rp**2, 1e-6)), 1500); rr = np.sqrt(Rp**2+zz**2)
    rho = np.interp(rr, H['r'], H['rho']); Phi = np.interp(rr, H['r'], H['Phi']); s2 = np.interp(rr, H['r'], H['sig2']); W = np.trapezoid(rho, zz)
    return np.trapezoid(rho*Phi, zz)/W, np.trapezoid(rho*s2, zz)/W
def binavg(H, a, b):
    # average over projected radius weighted by projected number in the annulus
    Rs = np.linspace(max(a, 0.05), b, 40); vals = []; wts = []
    for Rp in Rs:
        pb, s2 = project(H, Rp); zz = np.linspace(0, math.sqrt(max(6.0**2-Rp**2, 1e-6)), 400); rr = np.sqrt(Rp**2+zz**2)
        Sig = np.trapezoid(np.interp(rr, H['r'], H['rho']), zz); vals.append((-(pb-H['Phi0'])/C, s2)); wts.append(Sig*Rp)
    v = np.array(vals); w = np.array(wts); return float(np.sum(w*v[:, 0])/w.sum()), float(np.sum(w*v[:, 1])/w.sum())
pred = {}
for lab, M in (('median', R['Mmed_1e12']*1e12), ('upper quartile', R['Mq75']*1e12), ('2e14', 2e14)):
    H = nfw(M); gb = []; sb = []
    for a, b in EDG: g_, s2 = binavg(H, a, b); gb.append(g_); sb.append(s2)
    gb = np.array(gb); sb = np.array(sb); bcg = 2.2*(M/2e14)**(2/3)
    out = {'M': M, 'R200': H['R200'], 'gr_bins': gb.round(2).tolist(), 'sig_los_bins': np.sqrt(sb).round(0).tolist()}
    for K in (0.0, 10.0):
        kin = (2.5-K)*sb/C + bcg
        for nm, f in (('dark matter', 1.0), ('plasma only', FB)):
            tot = f*gb + kin; out[f'{nm}|K{K:g}'] = dict(bins=tot.round(2).tolist(), integ=float(np.sum(WTS*tot)/WTS.sum()))
    pred[lab] = out
    print(lab, f"M={M:.2e} R200={H['R200']:.2f}", 'gr', out['gr_bins'], 'sig', out['sig_los_bins'])
    for k in ('dark matter|K0', 'plasma only|K0', 'dark matter|K10', 'plasma only|K10'): print('   ', k, out[k]['bins'], round(out[k]['integ'], 1))
# comparison with the measurements
S = R['summ']
print('measured integrated (1:1:3:3):', {m: (S[m]['I'], S[m]['eI']) for m in S})
cmp = {}
for m, Kuse in (('0', 10.0), ('0.5', 10.0), ('1', 0.0), ('1.5', 0.0)):
    I, e = S[m]['I'], S[m]['eI']; P = pred['median']
    d = P[f'dark matter|K{Kuse:g}']['integ']; p = P[f'plasma only|K{Kuse:g}']['integ']
    cmp[m] = dict(I=I, e=e, K=Kuse, dm=d, pl=p, z_dm=(I-d)/e, z_pl=(I-p)/e)
    print(f'margin {m}: {I:+.1f} ± {e:.1f}; K={Kuse:g}: DM {d:+.1f} ({(I-d)/e:+.1f}σ), plasma {p:+.1f} ({(I-p)/e:+.1f}σ)')
json.dump(dict(raw=R, pred=pred, cmp=cmp), open('grz_tempel.json', 'w'), indent=1)
