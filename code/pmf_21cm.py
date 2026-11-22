"""pmf_21cm.py - Workbench Entry 138
Does a primordial magnetic field strong enough to matter leave the 21-cm sky dark?
A field left from the early universe heats the neutral gas after recombination in two ways (Sethi & Subramanian 2005;
Chluba et al. 2015, MNRAS 451, 2244, whose equations are used here):
  ambipolar diffusion  G_am = (1 - x)/(gamma x rho_b^2) * rho_B^2 * l_d^-2 * f_L(n_B + 3)
                       gamma = <sigma v>_HH+ / (2 m_H), <sigma v> = 6.49e-10 (T/K)^0.375 cm^3/s
                       f_L(y) = 0.8313 (1 - 1.020e-2 y) y^1.105 ; l_d = a / k_D ; k_D = 286.91 (nG/B) Mpc^-1
  decaying turbulence  G_tu = 1.5 m H rho_B / u ,  u = L + 1.5 ln[(1+z_i)/(1+z)],  m = 2(n_B+3)/(n_B+5),
                       L = ln(1 + 14.8 (B/nG)^-1 (k_D/Mpc^-1)^-1) = 0.050, z_i = 1088. Chluba et al. print the ratio as t_i/t_d;
                       in Sethi & Subramanian's decay law, E_B = E_B0 / [1 + ln(t/t_i)/ln(1 + t_d/t_i)]^m, it is t_d/t_i = k_J/k_D,
                       the Alfven crossing time of the damping scale over the age (about 1/17 for 1 nG at z_i, checked by hand).
                       Either way the number in the logarithm is 0.0516.
  (G_tu written on the current field energy; with turbulence alone it is exactly the published rate.)
Unlike the CMB papers, which ignore it, the field's own energy is drained by the heat it gives up:
  d(rho_B a^4)/dt = -(G_am + G_tu) a^4.
B0 is the comoving rms field at z_i before any of this heating.
Variants: S (main) = the field's spectrum in 64 shells, k^(n_B+3) exp(-k^2/k_D^2) per ln k, from 3.3 k_D down to 1.5e-8 k_D.
            Ambipolar heating = pre * rho_B,total * sum_k rho_k k^2, calibrated to the published rate at the start, and each shell
            gives up heat in proportion to rho_k k^2: small scales drain first, large scales survive.
          A = one shell at the initial damping scale: the whole field drains at the small-scale rate (most heat, earliest).
          B = one shell whose damping scale follows the decaying field (k_D = 286.91/B_now).
          C = S with the ambipolar rate divided by (8 pi)^2, the normalization Chluba et al. had to apply to reproduce Kunze & Komatsu (2014).
Hydrogen recombination: Peebles three-level atom with the RECFAST fudge 1.14 (Seager et al. 1999), alpha_B on the gas temperature,
photo-ionization on the CMB temperature; collisional ionization, excitation and ionization cooling from Cen (1992).
Gas temperature: dT/dt = -2HT + Compton + heating/(1.5 k n_tot)  (Chluba et al. 2015, eq. for T_e).
21 cm: dT_b = 27 x_HI (1 - T_g/T_S) sqrt(0.15/(Om h^2) (1+z)/10) (Ob h^2/0.023) mK (Furlanetto, Oh & Briggs 2006).
  Cosmic dawn (z <= 25): Lyman-alpha coupling taken as complete, T_S = T_k: the deepest trough the gas temperature allows.
  Dark ages: collisional coupling only, kappa_HH = 3.1e-11 T^0.357 exp(-32/T), kappa_eH from Kuhlen, Madau & Montgomery (2006).
No X-ray, Lyman-alpha or structure heating: the field is the only heat source, so each result is the field's effect alone.
"""
import json, math, os, time
import numpy as np
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))

# --- constants (cgs) ---
kB = 1.380649e-16; mH = 1.6735e-24; me = 9.1093837e-28; c = 2.99792458e10; hP = 6.62607015e-27
sT = 6.6524587e-25; aR = 7.5657e-15; Mpc = 3.08568e24; eV = 1.602176634e-12
# --- cosmology (Planck 2018) ---
h = 0.674; Om = 0.315; obh2 = 0.02237; Y = 0.245; T0 = 2.7255
H0 = 100*h*1e5/Mpc; Or = 4.18e-5/h**2; OL = 1 - Om - Or
rho_b0 = obh2*1.87847e-29; nH0 = (1 - Y)*rho_b0/mH; fHe = Y/(3.9715*(1 - Y))
def Hz(z): return H0*math.sqrt(Om*(1 + z)**3 + Or*(1 + z)**4 + OL)
def Tg(z): return T0*(1 + z)
# --- hydrogen atom ---
E_ion = 13.605693*eV; E2 = E_ion/4; E21 = E_ion - E2; L2s = 8.2245809; lam_a = 121.5670e-7; Ffudge = 1.14
def alphaB(T):
    t = T/1e4; return Ffudge*1e-13*4.309*t**-0.6166/(1 + 0.6703*t**0.5300)
def betaB(Tr):
    return alphaB(Tr)*(2*math.pi*me*kB*Tr/hP**2)**1.5*math.exp(-E2/(kB*Tr))
def saha_x(z):
    Tr = Tg(z); n = nH0*(1 + z)**3
    S = (2*math.pi*me*kB*Tr/hP**2)**1.5*math.exp(-E_ion/(kB*Tr))/n
    return (-S + math.sqrt(S*S + 4*S))/2

zi, dzi = 1088.0, 0.05*1088.0
def fL(y): return 0.8313*(1 - 1.020e-2*y)*y**1.105

LNK = np.arange(1.2, -18.0, -0.3)                                   # ln(k/k_D) bins for the spectral variant
def rhs_factory(B0, nB, variant):
    Bc0 = B0*1e-9; Ec0 = Bc0**2/(8*math.pi)                        # comoving field energy, erg/cm^3
    m = 2*(nB + 3)/(nB + 5); L = math.log(1 + 14.8/286.91); fl = fL(nB + 3)
    kD0 = 286.91/B0 if B0 > 0 else 1.0                               # comoving, 1/Mpc
    spectral = variant in ('S', 'C')
    amp = 1/(8*math.pi)**2 if variant == 'C' else 1.0
    if spectral:
        q = np.exp(LNK); w = q**(nB + 3)*np.exp(-q*q); w /= w.sum()
        e0 = Ec0*w; k2 = (q*kD0)**2                                    # comoving k^2, 1/Mpc^2
        calib = fl*kD0**2/float((w*k2).sum())                          # equals the published rate at the start
    else:
        e0 = np.array([Ec0])
    def rhs(z, y):
        x, T = min(max(y[0], 1e-12), 1.0), max(y[1], 1e-3); e = np.maximum(y[2:], 0.0)
        n = nH0*(1 + z)**3; H = Hz(z); Tr = Tg(z)
        # recombination (Peebles) + collisional ionization
        K = lam_a**3/(8*math.pi*H); b = betaB(Tr)
        C = (1 + K*L2s*n*(1 - x))/(1 + K*(L2s + b)*n*(1 - x))
        dxdt = -C*(alphaB(T)*n*x*x - b*(1 - x)*math.exp(-E21/(kB*Tr)))
        ccol = 5.85e-11*math.sqrt(T)/(1 + math.sqrt(T/1e5))*math.exp(-157809.1/T)
        dxdt += ccol*n*x*(1 - x)
        # field heating
        Et = float(e.sum()); Gam = Gtu = 0.0; de = np.zeros_like(e)
        if B0 > 0 and Et > 0:
            S = 0.5*(1 - math.tanh((z - zi)/dzi)); rb = rho_b0*(1 + z)**3
            gam = 6.49e-10*T**0.375/(2*mH); pre = amp*(1 - x)/(gam*x*rb*rb)*S
            u = L + 1.5*math.log(max((1 + zi)/(1 + z), 1.0)); ftu = 1.5*m*H/u*S       # fractional turbulent decay rate
            if spectral:
                ek2 = e*k2; sk = float(ek2.sum())
                Gam = pre*calib*Et*sk*(1 + z)**10/Mpc**2
                de = -(Gam*ek2/sk + ftu*e*(1 + z)**4)/(1 + z)**4 if sk > 0 else -ftu*e
            else:
                rB = Et*(1 + z)**4
                kD = 286.91/(math.sqrt(8*math.pi*Et)*1e9) if variant == 'B' else kD0
                Gam = pre*rB*rB*(kD*(1 + z)/Mpc)**2*fl
                de = np.array([-(Gam + ftu*rB)/(1 + z)**4])
            Gtu = ftu*Et*(1 + z)**4
        # temperature
        ntot = n*(1 + fHe + x)
        comp = 8*sT*aR*Tr**4*x/(3*me*c*(1 + fHe + x))*(Tr - T)
        cool = (7.5e-19*math.exp(-118348/T) + 1.27e-21*math.sqrt(T)*math.exp(-157809.1/T))/(1 + math.sqrt(T/1e5))*(n*x)*(n*(1 - x))
        dTdt = -2*H*T + comp + (Gam + Gtu - cool)/(1.5*kB*ntot)
        f = -1/(H*(1 + z))                                            # dt/dz
        return np.concatenate([[dxdt*f, dTdt*f], de*f])
    return rhs, Ec0, e0

ZGRID = np.concatenate([np.arange(1600, 200, -10.0), np.arange(200, 30, -1.0), np.arange(30, 7.9, -0.1)])
def run(B0, nB=-2.9, variant='S'):
    rhs, Ec0, e0 = rhs_factory(B0, nB, variant)
    z0 = 1600.0
    atol = np.concatenate([[1e-12, 1e-6], np.full(len(e0), max(Ec0, 1e-60)*1e-10)])
    sol = solve_ivp(rhs, (z0, 7.9), np.concatenate([[saha_x(z0), Tg(z0)], e0]), method='LSODA', t_eval=ZGRID, rtol=1e-7,
                    atol=atol, max_step=5.0)
    assert sol.success, sol.message
    E = sol.y[2:].sum(axis=0)/Ec0 if Ec0 > 0 else np.ones_like(sol.t)
    return dict(z=sol.t, x=sol.y[0], T=sol.y[1], E=E)

# --- 21 cm ---
A10 = 2.85e-15; Tstar = 0.0682
def kHH(T): return 3.1e-11*T**0.357*math.exp(-32/T)
def keH(T):
    lt = math.log10(min(T, 1e4)); return 10**(-9.607 + 0.5*lt*math.exp(-lt**4.5/1800))
pref = lambda z: 27*math.sqrt(0.15/(Om*h*h)*(1 + z)/10)*(obh2/0.023)
def dTb(z, x, T, coupled):
    Tr = Tg(z)
    if coupled: Ts = T
    else:
        n = nH0*(1 + z)**3
        xc = Tstar/(A10*Tr)*(n*(1 - x)*kHH(T) + n*x*keH(T))
        Ts = (1 + xc)/(1/Tr + xc/T)
    return pref(z)*(1 - x)*(1 - Tr/Ts)
def at(r, z, k): return float(np.interp(-z, -r['z'], r[k]))
def summary(r):
    out = {}
    for zz in (17.0, 20.0, 15.0, 10.4, 7.9):
        out[f'Tk_{zz}'] = at(r, zz, 'T'); out[f'xe_{zz}'] = at(r, zz, 'x')
    out['dTb17'] = dTb(17, at(r, 17, 'x'), at(r, 17, 'T'), True)
    out['dTb20'] = dTb(20, at(r, 20, 'x'), at(r, 20, 'T'), True)
    da = [(zz, dTb(zz, at(r, zz, 'x'), at(r, zz, 'T'), False)) for zz in np.arange(30, 200, 1.0)]
    zmin, vmin = min(da, key=lambda p: p[1])
    out['dark_min_z'] = float(zmin); out['dark_min_mK'] = float(vmin)
    out['dTb85'] = dTb(85, at(r, 85, 'x'), at(r, 85, 'T'), False)
    out['field_left_17'] = at(r, 17, 'E'); out['field_left_85'] = at(r, 85, 'E')
    return out

def solve_B(target, nB, variant, key, lo=1e-4, hi=3.0):
    """Field at which `key` (a depth in mK, negative) has risen to `target` mK."""
    flo = summary(run(lo, nB, variant))[key] - target; fhi = summary(run(hi, nB, variant))[key] - target
    if flo*fhi > 0: return None
    for _ in range(26):
        mid = math.sqrt(lo*hi); fm = summary(run(mid, nB, variant))[key] - target
        if fm*flo > 0: lo, flo = mid, fm
        else: hi, fhi = mid, fm
        if hi/lo < 1.01: break
    return math.sqrt(lo*hi)

if __name__ == '__main__':
    t0 = time.time(); R = {}
    base = run(0.0); S0 = summary(base)
    R['standard'] = S0
    R['check'] = dict(xe_1089=at(base, 1089, 'x'), xe_200=at(base, 200, 'x'), xe_50=at(base, 50, 'x'), xe_20=at(base, 20, 'x'),
                      Tk_200=at(base, 200, 'T'), Tk_100=at(base, 100, 'T'), Tk_50=at(base, 50, 'T'), Tk_17=at(base, 17, 'T'), Tk_10p4=at(base, 10.4, 'T'), Tg_17=Tg(17))
    print('standard', json.dumps(R['check'])); print('standard 21cm', S0['dTb17'], S0['dark_min_z'], S0['dark_min_mK'], S0['dTb85'])
    Bs = [0.0045, 0.01, 0.03, 0.045, 0.1, 0.27, 0.69]
    R['table'] = {}
    for nB, vs in ((-2.9, ('S', 'A', 'B', 'C')), (-2.95, ('S',)), (-2.0, ('S',)), (2.0, ('S',))):
        for v in vs:
            for B0 in Bs:
                s = summary(run(B0, nB, v)); R['table'][f'{nB}|{v}|{B0}'] = s
                print(f'nB={nB:+.1f} {v} B0={B0:<6} Tk17={s["Tk_17.0"]:9.2f} dTb17={s["dTb17"]:8.1f} dark={s["dark_min_mK"]:7.2f}@{s["dark_min_z"]:.0f} '
                      f'dTb85={s["dTb85"]:7.2f} left17={s["field_left_17"]:.3f} left85={s["field_left_85"]:.3f} xe17={s["xe_17.0"]:.2e} Tk10.4={s["Tk_10.4"]:.2f}')
    # thresholds: the field that halves the standard cosmic-dawn depth, brings it to -100 mK, erases it; halves the dark-ages trough
    R['thresholds'] = {}
    for nB, vs in ((-2.9, ('S', 'A', 'B', 'C')), (-2.95, ('S',)), (-2.0, ('S',)), (2.0, ('S',))):
        for v in vs:
            th = dict(half17=solve_B(S0['dTb17']/2, nB, v, 'dTb17'), zero17=solve_B(0.0, nB, v, 'dTb17'),
                      m100_17=solve_B(-100.0, nB, v, 'dTb17'), half_dark85=solve_B(S0['dTb85']/2, nB, v, 'dTb85'))
            R['thresholds'][f'{nB}|{v}'] = th; print('thr', nB, v, th)
    # the one-field corner of Entries 14 and 115: 45 pG on 1 Mpc, nearly scale-invariant. rms B0 = B_1Mpc (k_D/1 Mpc^-1)^((n+3)/2), k_D = 286.91/B0
    R['lofar'] = {}
    for nB in (-2.9, -2.95):
        B0 = 0.045
        for _ in range(60): B0 = 0.045*(286.91/B0)**((nB + 3)/2)
        for v in ('S', 'A', 'C'):
            s = summary(run(B0, nB, v)); R['lofar'][f'{nB}|{v}'] = dict(B0=B0, dTb17=s['dTb17'], dTb85=s['dTb85'], Tk17=s['Tk_17.0'], Tk10p4=s['Tk_10.4'])
    print('lofar', json.dumps(R['lofar']))
    # Entry 105's preferred clumping field, 5-10 pG (one-number fit 4.2 pG), nearly scale-invariant
    R['h0fields'] = {str(b): summary(run(b, -2.9, 'S'))['dTb17'] for b in (0.0042, 0.005, 0.01)}
    print('h0fields', R['h0fields'])
    # the reionization-limit fields of Entry 131 (rms): n_B = -2 at 0.27 nG, n_B = +2 at 0.18 nG
    R['reionlim'] = {}
    for nB, b0 in ((-2.0, 0.27), (2.0, 0.18)):
        s = summary(run(b0, nB, 'S')); R['reionlim'][f'{nB}|{b0}'] = dict(dTb17=s['dTb17'], dTb85=s['dTb85'])
    # where the field's energy goes: share left at z = 1000, 150 (gas leaves the CMB's grip), 17
    R['energy_path'] = {}
    for nB in (-2.9, -2.0, 2.0):
        r = run(0.068, nB, 'S'); R['energy_path'][str(nB)] = {str(zz): at(r, zz, 'E') for zz in (1000, 600, 300, 150, 85, 17)}
    print('reionlim', R['reionlim']); print('energy', json.dumps(R['energy_path']))
    # HERA (95%): kinetic temperature > 3.2 K and spin temperature > 4.7 K at z = 10.4; the field that alone lifts the gas above each
    R['hera'] = {}
    for Tlim in (3.2, 4.7):
        for v in ('S', 'A', 'C'):
            lo, hi = 1e-4, 3.0
            for _ in range(40):
                mid = math.sqrt(lo*hi)
                if at(run(mid, -2.9, v), 10.4, 'T') > Tlim: hi = mid
                else: lo = mid
            R['hera'][f'{Tlim}|{v}'] = math.sqrt(lo*hi)
    print('HERA field', R['hera'])
    # curves for the figures (n_B = -2.9, main variant)
    R['curves'] = {}
    for B0 in (0.0, 0.01, 0.03, 0.1, 0.3):
        r = run(B0); m = r['z'] <= 300
        R['curves'][str(B0)] = dict(z=[round(float(a), 2) for a in r['z'][m][::3]], T=[round(float(a), 4) for a in r['T'][m][::3]])
    grid = [float(v) for v in np.round(np.logspace(-3, math.log10(0.7), 25), 5)]
    R['depth_grid'] = dict(B=grid)
    for lab, nB, v, key in (('cd_S', -2.9, 'S', 'dTb17'), ('cd_A', -2.9, 'A', 'dTb17'), ('cd_C', -2.9, 'C', 'dTb17'),
                            ('cd_S_m2', -2.0, 'S', 'dTb17'), ('cd_S_p2', 2.0, 'S', 'dTb17'),
                            ('da_S', -2.9, 'S', 'dTb85'), ('da_A', -2.9, 'A', 'dTb85'), ('da_C', -2.9, 'C', 'dTb85')):
        R['depth_grid'][lab] = [round(summary(run(b, nB, v))[key], 3) for b in grid]
    R['runtime_s'] = round(time.time() - t0, 1)
    json.dump(R, open(os.path.join(HERE, 'pmf_21cm.json'), 'w'), indent=1)
    print('done', R['runtime_s'])
