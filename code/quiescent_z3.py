"""quiescent_z3.py - Workbench Entry 63
Entry 62 audited the line-width test against the acceleration relation of nearby galaxies and said where the audit
should bite hardest: quiescent galaxies beyond redshift 3, which are denser than those at redshift 1 to 2, so that a
static distance law cannot move them into the regime where unseen mass is normal. This script runs that case.

Sample: quiescent galaxies with a stellar velocity dispersion and a size at z = 2.8-4.7.
  Forrest et al. 2022 (ApJ 938, 109, table 1): eight of their own, four from Esdaile et al. 2021, one from Saracco et al. 2020
  (a fourteenth, Tanaka et al. 2019, has only an upper limit on its size and is left out), and GS-9209 from Carnall et al. 2023.
  Dispersions are corrected to the half-light radius; radii are at rest-frame 5000 A (GS-9209: 1.2 um); H0 = 70, Omega_m = 0.3.
Dynamical mass: as published. Forrest et al. (2022) use 5 sigma^2 r / G; a plain recomputation from the listed dispersion and radius
agrees within 0.1 dex except for COS-DR3-111740 (+0.33) and 3D-EGS-31322 (-0.39), whose listed values are kept.
Nearby anchor, as in Entries 45 and 62: stars are 0.59 of the dynamical mass in SDSS quiescent galaxies above 1e11 Msun
(van de Sande et al. 2013), at a radius of 5.4 kpc (Shen et al. 2003). That anchor used a virial factor that depends on the
light profile (5.2 to 6.6 in their distant sample); a flat 5 makes the distant dynamical masses, and so the static excess, a little smaller.
"""
import json, math
import numpy as np
from scipy import integrate
C, H0, OM = 299792.458, 70.0, 0.3
CH = C / H0
zg = np.linspace(0, 12, 12001); _dc = integrate.cumulative_trapezoid(1 / np.sqrt(OM * (1 + zg) ** 3 + 1 - OM), zg, initial=0)
DM = lambda z: CH * float(np.interp(z, zg, _dc))
LAWS = {
 "standard":                      (lambda z: DM(z) / (1 + z),        lambda z: DM(z) * (1 + z)),
 "tired light":                   (lambda z: CH * math.log(1 + z),   lambda z: CH * math.log(1 + z) * math.sqrt(1 + z)),
 "tired light, stretched time":   (lambda z: CH * math.log(1 + z),   lambda z: CH * math.log(1 + z) * (1 + z)),
 "Lerner's straight line":        (lambda z: CH * z,                 lambda z: CH * z * math.sqrt(1 + z)),
}
G = 4.30091e-3; ACC = 1e6 / 3.0857e16; A0 = 1.20e-10
nu = lambda gb: 1 / (1 - math.exp(-math.sqrt(gb / A0)))
# name, z, log M*, log M_dyn (published), sigma_e [km/s], error (mean of the two sides), r_e [kpc], error, source
GAL = [("COS-DR1-99209", 2.9834, 11.22, 11.31, 401, 74, 1.08, 0.24, "F22"), ("COS-DR3-84674", 3.0094, 11.25, 11.33, 442, 137, 0.95, 0.23, "F22"), ("COS-DR3-111740", 2.7988, 10.98, 11.02, 467, 117, 0.89, 0.33, "F22"),
       ("COS-DR3-201999", 3.1313, 11.40, 11.28, 271, 57, 2.26, 0.31, "F22"), ("COS-DR3-202019", 3.1326, 11.67, 12.00, 345, 102, 7.54, 1.16, "F22"), ("XMM-VID1-2075", 3.4520, 11.52, 11.49, 379, 69, 1.85, 0.16, "F22"),
       ("XMM-VID3-1120", 3.4919, 11.47, 11.54, 419, 111, 1.71, 0.19, "F22"), ("XMM-VID3-2457", 3.4892, 11.26, 11.49, 396, 86, 1.71, 0.22, "F22"),
       ("ZF-COS-20115", 3.715, 11.06, 10.86, 283, 52, 0.66, 0.08, "E21"), ("3D-EGS-40032", 3.219, 11.31, 11.41, 275, 56, 2.40, 0.19, "E21"), ("3D-EGS-18996", 3.239, 10.99, 10.56, 196, 48, 0.63, 0.05, "E21"), ("3D-EGS-31322", 3.434, 10.99, 10.85, 201, 119, 0.61, 0.05, "E21"),
       ("C1-23152", 3.352, 11.30, 11.34, 409, 60, 1.0, 0.1, "S20"), ("GS-9209", 4.658, 10.58, 10.30, 247, 16, 0.215, 0.020, "C23")]
shen_early = lambda m: 2.88e-6 * m ** 0.56
m_loc = 10 ** 11.2; r_loc = shen_early(m_loc); z_loc = 0.06
def anchor(law):
    da, dl = LAWS[law]; fA, fL2 = da(z_loc) / LAWS["standard"][0](z_loc), (dl(z_loc) / LAWS["standard"][1](z_loc)) ** 2
    gb = G * m_loc * fL2 / 2 / (r_loc * fA * 1e3) ** 2 * ACC; ratio = (1 / 0.59) * fA / fL2
    return gb, nu(gb), ratio, ratio / nu(gb)
rng = np.random.default_rng(63); out = {"n": len(GAL), "laws": {}, "galaxies": []}
print("%d galaxies, redshift %.2f to %.2f, median %.2f; median dispersion %.0f km/s, radius %.2f kpc, log stellar mass %.2f" % (len(GAL), min(g[1] for g in GAL), max(g[1] for g in GAL), np.median([g[1] for g in GAL]), np.median([g[4] for g in GAL]), np.median([g[6] for g in GAL]), np.median([g[2] for g in GAL])))
print("law                           radius   g_bar/a0   dyn/stars    over nearby, same mass      relation allows    left over                  galaxies above nearby   stars > dynamical")
for law in LAWS:
    da, dl = LAWS[law]; gbl, nul, ratl, call = anchor(law); a = []; b = []; rr = []; gg = []; ds = []; al = []; rows = []
    for nm, z, lms, lmd, sig, esig, re_, ere, src in GAL:
        fA, fL2 = da(z) / LAWS["standard"][0](z), (dl(z) / LAWS["standard"][1](z)) ** 2
        r = re_ * fA; ms = 10 ** lms * fL2; md = 10 ** lmd * fA; gb = G * ms / 2 / (r * 1e3) ** 2 * ACC
        a.append(math.log10(md / ms / ratl)); b.append(math.log10(md / ms / nu(gb) / call)); rr.append(r); gg.append(gb / A0); ds.append(md / ms); al.append(math.log10(nu(gb) / nul))
        rows.append([nm, z, round(math.log10(md / ms / nu(gb) / call), 3), round(math.sqrt((2 * esig / sig) ** 2 + (ere / re_) ** 2) / math.log(10), 3)])
    a, b = np.array(a), np.array(b); boot = [np.median(rng.choice(b, len(b))) for _ in range(4000)]
    q = dict(r=float(np.median(rr)), gbar_over_a0=float(np.median(gg)), dyn_over_star=float(np.median(ds)), same_mass=[float(a.mean()), float(a.std(ddof=1) / len(a) ** 0.5)], allowed=float(np.mean(al)),
             same_acc=[float(b.mean()), float(b.std(ddof=1) / len(b) ** 0.5)], same_acc_median=[float(np.median(b)), float(np.std(boot))], n_above=int((b > 0).sum()), n_star_gt_dyn=int(sum(d < 1 for d in ds)), rows=rows)
    out["laws"][law] = q
    print("%-28s  %5.1f    %6.1f     %5.2f       x%5.2f (log %+.2f +- %.2f)    x%.2f              x%5.2f (log %+.2f +- %.2f)    %2d of %d               %d" %
          (law, q["r"], q["gbar_over_a0"], q["dyn_over_star"], 10 ** a.mean(), a.mean(), q["same_mass"][1], 10 ** q["allowed"], 10 ** b.mean(), b.mean(), q["same_acc"][1], q["n_above"], len(GAL), q["n_star_gt_dyn"]))
# without the two with the largest radii or errors, and by source
for lab, keep in (("without COS-DR3-202019 (7.5 kpc)", lambda g: g[0] != "COS-DR3-202019"), ("Forrest et al. only", lambda g: g[8] == "F22"), ("others only", lambda g: g[8] != "F22"), ("dispersion known to 25% or better", lambda g: g[5] / g[4] <= 0.25), ("plain 5 sigma^2 r / G within 0.1 dex of the listed mass", lambda g: abs(math.log10(5 * g[4] ** 2 * g[6] * 1e3 / G) - g[3]) <= 0.1)):
    line = []
    for law in LAWS:
        v = [r_[2] for r_, g in zip(out["laws"][law]["rows"], GAL) if keep(g)]; line.append("%+.2f" % np.mean(v))
    print("  %-36s (n=%2d)  left over, log: %s" % (lab, sum(keep(g) for g in GAL), "  ".join(line))); out.setdefault("checks", {})[lab] = line
# a profile-dependent virial factor of 6 in place of 5 raises every dynamical mass by 0.08 dex
print("  a virial factor of 6 in place of 5 adds %.2f dex to every row" % math.log10(6 / 5))
# growth of the static factor with redshift: size factor over luminosity factor
print("\nsize over luminosity factor, static over standard:"); out["F"] = {}
for law in list(LAWS)[1:]:
    da, dl = LAWS[law]; v = [(z, (da(z) / LAWS["standard"][0](z)) / (dl(z) / LAWS["standard"][1](z)) ** 2) for z in (1.0, 1.3, 2.0, 3.3, 4.7)]
    out["F"][law] = v; print("  %-28s " % law + "  ".join("z=%.1f: x%.1f" % t for t in v))
json.dump(out, open("quiescent_z3.json", "w"), indent=1, default=float)
