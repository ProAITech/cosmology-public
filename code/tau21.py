"""Entry 152: the 21-cm optical depth of the radio CMB (Holder et al. 2026; Fronenberg 2026) against the notebook's
primordial-field heating (Entry 138) and the radio background (Entry 75).
tau_HI = 0.00847 C (1+z)^1.5 x_HI / T_s (Fronenberg eq. 1, C ~ 1 for Planck), T_b = tau/(1+z) (T_s - T_R) (eq. 3).
Entry 138's depth grid gives the trough depth at z = 17 for each field strength with full Lyman-alpha coupling (T_s = T_k)
and no other heating; we invert it to T_s and then to tau. Array sensitivities from Fronenberg's Table 1:
Delta alpha = (1 - e^-tau)/sigma_null for each scenario and model."""
import json, math, numpy as np
P = json.load(open('../pmf21/pmf_21cm.json'))
z = 17.0; T0 = 2.7255; TR0 = T0*(1+z)
std = P['standard']; Ts_std = std['Tk_17.0']; dT_std = std['dTb17']
F = dT_std/(1 - TR0/Ts_std)                       # calibrates Entry 138's 21-cm prefactor at z = 17
tau = lambda Ts, zz=z: 0.00847*(1+zz)**1.5/Ts
G = P['depth_grid']; B = np.array(G['B'])
out = dict(z=z, Ts_std=Ts_std, tau_std=tau(Ts_std), F=F, curves={})
for v in ('S', 'A', 'C'):
    dT = np.array(G['cd_'+v]); Ts = TR0/(1 - dT/F)
    out['curves'][v] = dict(B=B.tolist(), Ts=Ts.round(3).tolist(), tau=[round(tau(t), 5) for t in Ts])
# Fronenberg Table 1 (sigma_null), models EDGES tau 0.197 and fiducial tau 0.04 (z = 17.2 and 12.25)
T1 = {'HERA': (0.51, 0.14, 0.25), 'Threshold': (1.21, 0.34, 0.59), 'Dedicated': (1.92, 0.54, 0.97), 'Robust': (4.23, 1.19, 1.77)}
sig = {k: dict(edges=(1-math.exp(-0.197))/v[0], fid17=(1-math.exp(-0.04))/v[1]) for k, v in T1.items()}
out['sigma_alpha'] = sig
# bounds on B if the measured tau equals the no-heating value (adiabatic gas, full coupling), at 2 sigma
def Bbound(sig_tau, v='S', tau_true=None):
    t0 = out['tau_std'] if tau_true is None else tau_true
    tmin = t0 - 2*sig_tau
    if tmin <= 0: return None
    c = out['curves'][v]; tt = np.array(c['tau'])
    # tau decreases with B; find B where tau = tmin
    if tt.min() > tmin: return None
    i = np.where(tt < tmin)[0][0]
    return float(np.exp(np.interp(tmin, [tt[i], tt[i-1]], [np.log(B[i]), np.log(B[i-1])])))
cases = {}
for nm, s in (('Robust', sig['Robust']['fid17']), ('Robust x10 (assumed)', sig['Robust']['fid17']/10), ('Robust x30 (assumed)', sig['Robust']['fid17']/30)):
    cases[nm] = dict(sig_tau=s, B_S=Bbound(s, 'S'), B_A=Bbound(s, 'A'), B_C=Bbound(s, 'C'))
out['bounds_noheat'] = cases
# with cosmic-dawn heating as in Fronenberg's fiducial model (tau 0.04 at z = 17.2): the field can at most supply that heat
Ts_fid = tau(1.0, 17.2)/0.04
out['Ts_fid'] = Ts_fid
c = out['curves']['S']; out['B_for_Tfid'] = {v: float(np.exp(np.interp(Ts_fid, out['curves'][v]['Ts'], np.log(B)))) for v in ('S', 'A', 'C')}
# radio background an EDGES-depth trough (-500 mK at z = 17.2) needs, for a given measured tau
zE = 17.2; TCMB = T0*(1+zE)
def TR_needed(t, Tb=-0.5):
    Ts = tau(1.0, zE)/t; return Ts - Tb*(1+zE)/t
grid = np.geomspace(0.02, 0.25, 40)
out['radio'] = dict(tau=grid.round(4).tolist(), TR_over_TCMB=[round(TR_needed(t)/TCMB, 3) for t in grid])
marks = {}
for lab, t in (('EDGES cold gas', 0.197), ('no heating', tau(Ts_std, zE)), ('fiducial heated', 0.04)):
    TR = TR_needed(t); exc_today = (TR - TCMB)/(1+zE)
    marks[lab] = dict(tau=t, Ts=tau(1.0, zE)/t, TR=TR, ratio=TR/TCMB, excess_today_K=exc_today)
# today's extragalactic radio background at 78 MHz from Entry 75: 30.4 K at 310 MHz, index -2.66
Tbg78 = 30.4*(310/78.0)**2.66
for m in marks.values(): m['frac_of_bg78'] = m['excess_today_K']/Tbg78
out['marks'] = marks; out['Tbg78'] = Tbg78
print(json.dumps({k: v for k, v in out.items() if k not in ('curves', 'radio')}, indent=1))
json.dump(out, open('tau21.json', 'w'))
