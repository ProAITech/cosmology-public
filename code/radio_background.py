"""Workbench Entry 75: Peratt's synchrotron microwave background, and whether current-carrying cosmic filaments could supply the
measured excess radio synchrotron background (RSB).  cgs units."""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
k = 1.380649e-16; c = 2.99792458e10; me_c2 = 8.1871e-7; h = 6.62607e-27; sigT = 6.6524e-25; eV = 1.602177e-12
trap = getattr(np, 'trapezoid', None) or np.trapz
out = {}
# 1. self-absorbed synchrotron cannot be a 2.7 K blackbody
Tb_min = me_c2 / (3 * k)
print(f"1. Self-absorbed synchrotron: brightness temperature ~ gamma m_e c^2 / 3k >= {Tb_min:.1e} K for any relativistic electron (gamma >= 1)")
print("   and its optically thick spectrum rises as nu^2.5, against the nu^2 of the measured 2.725 K Rayleigh-Jeans tail.")
out["Tb_min"] = Tb_min
# 2. energy density of the radio synchrotron background, T = 30.4 K (nu / 310 MHz)^-2.66 (Singal et al. 2023 report)
nu = np.logspace(np.log10(22e6), np.log10(10e9), 4000)
T = 30.4 * (nu / 310e6) ** -2.66
u_rsb = trap(8 * math.pi * k * T * nu**2 / c**3, nu)
u_cmb = 4.17e-13
print(f"\n2. RSB energy density 22 MHz-10 GHz: {u_rsb:.1e} erg/cm3; three quarters of it is the excess over known sources ({0.75*u_rsb:.1e})")
print(f"   CMB: {u_cmb:.2e} erg/cm3; the RSB is {u_rsb/u_cmb:.0e} of it")
out["u_rsb"] = u_rsb
# 3. inverse Compton companion of synchrotron in weak fields: P_IC / P_sync = U_CMB / U_B
u_gamma = {"Fermi >100 MeV (EGB 1.1e-5 ph/cm2/s/sr x ~0.3 GeV)": 4 * math.pi / c * 1.1e-5 * 0.3e9 * eV,
           "COMPTEL 1-30 MeV (~1e-2 MeV/cm2/s/sr per ln E)": 4 * math.pi / c * 1e-2 * 1e6 * eV * math.log(30)}
for kk, v in u_gamma.items(): print(f"   observed {kk}: {v:.1e} erg/cm3")
u_obs = sum(u_gamma.values())
print("\n3. Electrons that make synchrotron in a field B also scatter the CMB: inverse-Compton power / synchrotron power = U_CMB / U_B")
rows = []
for label, B in (("our filament-current limit, 3e17 A at 1 Mpc", 2e-7 * 3e17 / 3.0857e22 * 1e4), ("LOFAR filament field, 15 nG", 15e-9), ("LOFAR 2022 estimate, 32 nG", 32e-9),
                 ("radio bridge A399-A401, 0.3 uG", 0.3e-6), ("galaxy disc, 5 uG", 5e-6)):
    UB = B * B / (8 * math.pi); ratio = u_cmb / UB
    gam = math.sqrt(1.4e9 / (4.2e6 * B / 1e-6 * 1e-6 * 1e6)) if False else math.sqrt(1.4e9 / (4.2 * B / 1e-6))
    E_ic = 4 / 3 * gam**2 * 6.3e-4      # eV, mean CMB photon energy 6.3e-4 eV
    u_ic = 0.75 * u_rsb * ratio
    rows.append(dict(label=label, B=B, ratio=ratio, gamma=gam, E_ic_eV=E_ic, u_ic=u_ic, over_obs=u_ic / u_obs))
    print(f"   {label:44s} B = {B*1e9:8.2f} nG: IC/sync = {ratio:.1e}; 1.4 GHz electrons gamma {gam:.1e}, IC photons ~{E_ic/1e6:.2g} MeV;"
          f" IC energy {u_ic:.1e} erg/cm3 = {u_ic/u_obs:.0e} x the observed MeV-GeV background")
out["rows"] = rows; out["u_obs"] = u_obs
Bmin = math.sqrt(8 * math.pi * u_cmb * 0.75 * u_rsb / u_obs)
print(f"   to keep the IC below the whole observed MeV-GeV background the emitting field must exceed {Bmin*1e6:.2f} uG")
out["Bmin"] = Bmin
json.dump(out, open(os.path.join(HERE, "radio_background.json"), "w"), indent=1)
