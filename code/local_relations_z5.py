"""local_relations_z5.py - Workbench Entry 62
In a universe that is static and does not evolve, galaxies at any redshift are the same population as nearby ones,
so they must lie on the relations that nearby galaxies obey. Entry 59 tested static distance laws with line widths
and concluded that stars and gas fall "short of the mass the line widths require". That assumed the baryons should
supply the mass. Nearby galaxies do not work that way at low surface density: they follow one relation between the
acceleration the baryons give and the acceleration observed (McGaugh, Lelli & Schombert 2016). This script re-reads
two samples under four distance laws and asks, relation by relation, which law keeps them on the local relations.

 Sample A  ALPINE (Bethermin et al. 2020, VizieR J/A+A/643/A2/catalog): 45 galaxies at z = 4.4-5.9 with a [C II]
           line width, a [C II] size (Fujimoto et al. 2020) and a [C II] luminosity. Gas mass = 30 x L[CII]
           (Zanella et al. 2018), so the gas is seen in a line and not inferred from the star formation rate.
           Catalogue values use H0 = 70, Omega_m = 0.3. Mergers (class 2) are left out of every test.
 Sample B  the seventeen JWST galaxies of Entry 59 (z = 4.0-7.6), gas inferred from the Kennicutt law.

 Tests     1. radial-acceleration relation at the half-light radius
           2. size against nearby late-type galaxies of the same stellar mass (Shen et al. 2003)
           3. star formation law: surface density of star formation against that of gas (Kennicutt 1998) - sample A only
           4. star formation rate over stellar mass against the nearby star-forming sequence (Salim et al. 2007) -
              this one does not depend on the distance law at all
"""
import json, math
import numpy as np
from scipy import integrate

C, H0, OM = 299792.458, 70.0, 0.3
CH = C / H0
zg = np.linspace(0, 12, 12001); _dc = integrate.cumulative_trapezoid(1 / np.sqrt(OM * (1 + zg) ** 3 + 1 - OM), zg, initial=0)
DM = lambda z: CH * float(np.interp(z, zg, _dc))
LAWS = {   # angular-size distance and luminosity distance, Mpc
 "standard":                      (lambda z: DM(z) / (1 + z),        lambda z: DM(z) * (1 + z)),
 "tired light":                   (lambda z: CH * math.log(1 + z),   lambda z: CH * math.log(1 + z) * math.sqrt(1 + z)),
 "tired light, stretched time":   (lambda z: CH * math.log(1 + z),   lambda z: CH * math.log(1 + z) * (1 + z)),
 "Lerner's straight line":        (lambda z: CH * z,                 lambda z: CH * z * math.sqrt(1 + z)),
}
G = 4.30091e-3                      # pc (km/s)^2 / Msun
ACC = 1e6 / 3.0857e16               # (km/s)^2 / pc  ->  m / s^2
GDAG = 1.20e-10                     # m / s^2
ARC = math.pi / 180 / 3600
shen = lambda m: 0.1 * m ** 0.14 * (1 + m / 3.98e10) ** 0.25           # kpc, late types
rar = lambda gb: gb / (1 - np.exp(-np.sqrt(gb / GDAG)))                 # observed acceleration for a baryonic one
ks = lambda sg: 2.5e-4 * sg ** 1.4                                      # Msun/yr/kpc^2 for Msun/pc^2
ms_local = lambda m: 10 ** (-0.35 * (math.log10(m) - 10) - 9.83)        # per yr

# ---- Sample A: name, z, FWHM [km/s], err, [C II] half-light radius [kpc, standard], err, size flag (0 reliable),
#      L[CII] [1e8 Lsun, standard], log M*, log SFR, class (1 rotator, 2 merger, 3 dispersion-dominated)
ALP = [["CG32",4.4105,279,13,1.94,0.3,0,8.19,9.75,2.36,2],["DC308643",4.5253,406,40,1.25,0.55,0,5.69,9.77,1.83,2],["DC351640",5.7058,132,13,2.05,0.95,0,2.35,9.67,1.2,3],["DC372292",5.1364,289,14,0.74,0.46,0,3.93,9.44,1.45,2],["DC396844",4.5424,287,13,2.56,0.33,0,11.53,9.86,1.74,1],["DC403030",4.5594,168,14,3.23,0.95,0,4.65,9.68,1.56,2],["DC416105",5.6309,202,14,1.02,0.86,1,1.34,9.45,1.08,1],["DC417567",5.67,310,16,2.07,0.58,0,3.13,9.81,1.87,2],["DC422677",4.4381,233,14,1.1,0.5,1,4.22,9.85,1.9,2],["DC432340",4.4045,151,13,0.72,6.26,1,3.57,10.1,1.79,2],["DC434239",4.4883,497,13,0.62,24.89,1,14.07,10.32,2,2],["DC454608",4.5834,232,15,2.95,0.72,0,6.44,9.79,1.67,2],["DC488399",5.6704,303,13,1.32,0.16,0,10.79,10.2,1.67,3],["DC493583",4.5134,198,13,1.89,0.51,0,4.32,9.61,1.36,2],["DC494057",5.5448,217,13,2.48,0.25,0,7.2,10.15,1.62,1],["DC494763",5.2337,253,14,1.12,0.36,0,4.87,9.53,1.3,3],
 ["DC519281",5.5759,282,14,0.36,1.28,1,5.36,9.87,1.4,2],["DC536534",5.6886,621,21,4.16,1.04,1,7.86,10.36,1.71,2],["DC539609",5.1818,287,14,1.65,0.43,0,4.91,9.38,1.48,1],["DC552206",5.5016,366,13,3.41,0.35,1,15.25,10.58,2.24,2],["DC627939",4.5341,252,13,2.12,0.31,0,7.24,9.98,1.71,2],["DC630594",4.4403,260,13,1.7,0.34,0,6.23,9.77,1.5,3],["DC683613",5.542,216,13,1.82,0.33,0,7.97,10.17,1.77,3],["DC709575",4.4121,254,16,1.34,1.62,1,2.99,9.68,1.42,1],["DC733857",4.5445,226,13,1.8,0.58,1,4.82,9.54,1.59,3],["DC773957",5.6773,344,13,2.68,0.54,0,4.18,10,1.43,2],["DC818760",4.5613,276,13,2.59,0.16,1,43,10.63,2.68,2],["DC834764",4.5058,254,17,3.17,0.89,0,5.63,9.97,1.36,3],["DC842313",4.5537,250,15,0.79,8.82,1,4.54,10.84,2.28,2],["DC848185",5.2931,275,13,3.47,0.25,1,16.02,10.37,2.46,3],["DC873321",5.1542,201,13,0.67,17.97,1,9.51,9.97,1.96,2],["DC873756",4.5457,526,13,2.36,0.11,0,36.22,10.25,0.73,2],
 ["DC880016",4.5415,274,14,2.41,0.5,0,5.53,9.75,1.5,3],["DC881725",4.5777,198,13,2.26,0.33,0,6.87,9.96,1.94,1],["VC5100537582",4.5501,206,15,1.73,0.53,0,4.43,9.76,1.17,3],["VC5100541407",4.563,177,13,4.85,0.59,1,11.8,10.12,1.54,2],["VC5100559223",4.5627,143,13,2.86,0.77,0,4.39,9.73,1.52,3],["VC5100822662",4.5205,208,13,2.59,0.37,0,7.86,10.17,1.82,2],["VC5100969402",4.5785,291,14,1.62,0.33,0,5.23,10,1.29,3],["VC5100994794",4.5802,230,13,1.86,0.32,0,5.57,9.73,1.45,3],["VC5101218326",4.5739,229,13,2.37,0.15,0,18.25,11.01,1.06,3],["VC510596653",4.5681,62,13,1.79,0.76,0,1.73,None,None,3],["VC510786441",4.4635,224,13,2.94,0.52,0,6.6,9.99,1.6,2],["VC5110377875",4.5505,234,13,2.63,0.23,0,17.14,10.17,1.99,1],["VC5180966608",4.5296,243,13,5.1,0.42,1,13.72,10.82,2.15,2]]
# whole catalogue, as computed in the browser run: 117 galaxies with a stellar mass
ALP_ALL = dict(n=117, log_ssfr=[-8.56, -8.27, -7.96], log_mstar=[9.16, 9.68, 10.17], z=[4.42, 4.58, 5.72], n_cii=74, log_tdep=[8.29, 8.63, 8.95], fgas=[0.48, 0.70, 0.81], classes={1: 9, 2: 31, 3: 15, 4: 8, 5: 11})
# ---- Sample B: name, z, r_e [arcsec], sigma [km/s], log M_dyn (published, standard), log M*, SFR  (Entry 59)
JW = [("JADES-NS-00016745",5.566,0.32,55.0,10.23,8.80,12.6),("JADES-NS-00019606",5.890,0.117,39.1,9.17,7.54,1.46),("JADES-NS-00022251",5.799,0.122,39.0,9.23,8.21,6.2),("JADES-NS-00047100",7.432,0.130,71.0,9.85,8.53,5.1),
 ("JADES-NS-10016374",5.504,0.084,62.0,9.45,7.86,2.19),("JADES-NS-20086025",7.263,0.36,25.0,10.31,8.85,24.0),("1871-10",7.5990,0.14,65.2,9.54,7.96,6.0),("1871-11",7.6100,0.17,60.9,9.77,9.07,9.7),("1871-12",7.5020,0.07,67.3,9.41,8.31,19.9),
 ("1871-29",5.0191,0.13,72.0,9.96,9.24,10.7),("1871-40",7.0948,0.14,60.6,9.73,8.47,9.2),("1871-70",7.0302,0.12,41.3,9.42,8.85,1.7),("1871-105",6.5669,0.07,49.9,9.34,7.81,2.8),("1871-451",4.4615,0.08,40.7,9.37,7.99,1.7),
 ("1871-545",4.0416,0.12,96.4,10.11,9.41,66.7),("1871-912",6.7161,0.18,44.4,9.68,9.69,27.1),("1871-1032",7.0889,0.11,37.3,9.41,9.07,6.3)]

rng = np.random.default_rng(62)
def med(v):
    v = np.array([x for x in v if x is not None and np.isfinite(x)], float)
    b = np.array([np.median(rng.choice(v, len(v))) for _ in range(4000)])
    return dict(median=float(np.median(v)), err=float(b.std()), p16=float(np.percentile(v, 16)), p84=float(np.percentile(v, 84)), n=int(len(v)))

def alpine(rows, law, vfac=2.5, alpha=30.0):
    da, dl = LAWS[law]; out = []
    for nm, z, fw, efw, r0, er0, fl, l8, lm, lsfr, cl in rows:
        if lm is None: continue
        fA, fL2 = da(z) / LAWS["standard"][0](z), (dl(z) / LAWS["standard"][1](z)) ** 2
        r = r0 * fA; ms = 10 ** lm * fL2; mg = alpha * l8 * 1e8 * fL2; sfr = 10 ** lsfr * fL2; sig = fw / 2.3548
        gbar = G * (ms + mg) / 2 / (r * 1e3) ** 2 * ACC; gobs = vfac * sig ** 2 / (r * 1e3) * ACC
        sg = mg / 2 / (math.pi * (r * 1e3) ** 2); ssf = sfr / 2 / (math.pi * r ** 2)
        out.append(dict(name=nm, z=z, r=r, mstar=ms, mgas=mg, sfr=sfr, sigma=sig, gbar=gbar, gobs=gobs, d_rar=math.log10(gobs / rar(gbar)), size_ratio=r / shen(ms), sig_gas=sg, sig_sfr=ssf,
                        d_ks=math.log10(ssf / ks(sg)), bar_over_dyn=(ms + mg) / (5 * sig ** 2 * r * 1e3 / G), d_btf=math.log10((ms + mg) / (47 * (vfac * sig ** 2) ** 2)), tdep=mg / sfr, ssfr=sfr / ms, d_ms=math.log10(sfr / ms / ms_local(ms))))
    return out

def jwst(law, gas=True):
    da, dl = LAWS[law]; out = []
    for nm, z, re_as, sig, lmd, lms, sfr0 in JW:
        fA, fL2 = da(z) / LAWS["standard"][0](z), (dl(z) / LAWS["standard"][1](z)) ** 2
        r = re_as * ARC * da(z) * 1e3; ms = 10 ** lms * fL2; md = 10 ** lmd * fA; sfr = sfr0 * fL2
        mg = (sfr / (2 * math.pi * r ** 2) / 2.5e-4) ** (1 / 1.4) * 2 * math.pi * (r * 1e3) ** 2 if gas else 0.0
        gbar = G * (ms + mg) / 2 / (r * 1e3) ** 2 * ACC; gobs = G * md / 2 / (r * 1e3) ** 2 * ACC
        out.append(dict(name=nm, z=z, r=r, mstar=ms, mgas=mg, gbar=gbar, gobs=gobs, d_rar=math.log10(gobs / rar(gbar)), size_ratio=r / shen(ms), bar_over_dyn=(ms + mg) / md, d_ms=math.log10(sfr / ms / ms_local(ms))))
    return out

res = {"laws": list(LAWS), "alpine": {}, "jwst": {}, "alpine_all": ALP_ALL}
sel = {"non-mergers, reliable size": [r for r in ALP if r[10] != 2 and r[6] == 0], "non-mergers, any size": [r for r in ALP if r[10] != 2], "all 45": ALP}
print("Sample A (ALPINE): %d with a size; %d not mergers; %d of those with a reliable size" % (len(ALP), len(sel["non-mergers, any size"]), len(sel["non-mergers, reliable size"])))
for sname, rows in sel.items():
    res["alpine"][sname] = {}
    print("\n  %s" % sname)
    print("  law                           radius kpc   log M*   log gas   RAR offset dex    size / local     star-formation law dex   baryons/dynamical   surface density of gas")
    for law in LAWS:
        R = alpine(rows, law)
        s = {k: med([r[k] for r in R]) for k in ("r", "d_rar", "size_ratio", "d_ks", "bar_over_dyn", "sig_gas", "sig_sfr", "gbar", "gobs", "d_btf", "tdep", "d_ms")}
        s["log_mstar"] = med([math.log10(r["mstar"]) for r in R]); s["log_mgas"] = med([math.log10(r["mgas"]) for r in R]); s["log_size_ratio"] = med([math.log10(r["size_ratio"]) for r in R])
        res["alpine"][sname][law] = s
        print("  %-28s  %6.1f      %5.2f    %5.2f     %+5.2f +- %.2f    x%5.1f           %+5.2f +- %.2f            %5.3f             %7.1f Msun/pc2" %
              (law, s["r"]["median"], s["log_mstar"]["median"], s["log_mgas"]["median"], s["d_rar"]["median"], s["d_rar"]["err"], s["size_ratio"]["median"], s["d_ks"]["median"], s["d_ks"]["err"], s["bar_over_dyn"]["median"], s["sig_gas"]["median"]))
main = "non-mergers, reliable size"
# points for the figure: (log g_bar, log g_obs) per galaxy per law
res["plane_alpine"] = {law: [[round(math.log10(r["gbar"]), 3), round(math.log10(r["gobs"]), 3)] for r in alpine(sel[main], law)] for law in LAWS}
res["ks_alpine"] = {law: [[round(math.log10(r["sig_gas"]), 3), round(math.log10(r["sig_sfr"]), 3)] for r in alpine(sel[main], law)] for law in LAWS}
# sensitivity of the acceleration test to how the line width is turned into a circular speed, and to the gas conversion
print("\n  Sensitivity, %s, RAR offset [dex] by law:" % main)
res["sensitivity"] = {}
for lab, kw in (("V^2 = 2.5 sigma^2, alpha = 30 (adopted)", {}), ("V^2 = 1.85 sigma^2 (disc seen at 60 deg)", dict(vfac=1.85)), ("V^2 = 3 sigma^2", dict(vfac=3.0)), ("alpha = 10", dict(alpha=10.0)), ("alpha = 100", dict(alpha=100.0))):
    v = [float(np.median([r["d_rar"] for r in alpine(sel[main], law, **kw)])) for law in LAWS]; res["sensitivity"][lab] = v
    print("    %-42s %s" % (lab, "  ".join("%+5.2f" % x for x in v)))
t = med([r["tdep"] for r in alpine(sel[main], "standard")]); res["tdep_main_myr"] = t["median"] / 1e6
print("\n  Gas over star formation rate (no distance enters): median %.0f Myr in this subsample; %.0f Myr for all 74 detections" % (t["median"] / 1e6, 10 ** ALP_ALL["log_tdep"][1] / 1e6))
print("  The Kennicutt law gives that time at a gas surface density of %.0f Msun/pc2; at 3 Msun/pc2 it gives %.1f Gyr" % (((t["median"] / 1e6) / 4000) ** (-1 / 0.4), 3 / ks(3) / 1e3))

print("\nSample B (JWST, Entry 59), gas from the Kennicutt law and, in brackets, stars only")
for law in LAWS:
    R, R0 = jwst(law), jwst(law, gas=False)
    s = {k: med([r[k] for r in R]) for k in ("r", "d_rar", "size_ratio", "bar_over_dyn", "gbar", "gobs", "d_ms")}; s0 = med([r["d_rar"] for r in R0]); s["d_rar_stars_only"] = s0
    s["log_size_ratio"] = med([math.log10(r["size_ratio"]) for r in R]); res["jwst"][law] = s
    print("  %-28s radius %5.2f kpc   RAR offset %+5.2f +- %.2f  (%+5.2f)   size / local x%5.1f   baryons/dynamical %.3f" % (law, s["r"]["median"], s["d_rar"]["median"], s["d_rar"]["err"], s0["median"], s["size_ratio"]["median"], s["bar_over_dyn"]["median"]))
res["plane_jwst"] = {law: [[round(math.log10(r["gbar"]), 3), round(math.log10(r["gobs"]), 3)] for r in jwst(law)] for law in LAWS}
# ---- Sample C: the 71 quiescent galaxies at z = 0.8-2.2 of Entry 45 (van de Sande et al. 2013, table 4) ---------------
# Entry 45 compared dynamical over stellar mass with nearby galaxies of the same mass. A static law makes these galaxies
# larger, so the fair nearby comparison is at the same acceleration. Nearby anchor: stars/dynamical = 0.59 for SDSS quiescent
# galaxies above 1e11 Msun (same paper), which have half-light radii of about 5.4 kpc (Shen et al. 2003, early types, at 10^11.2).
import os
HERE = os.path.dirname(os.path.abspath(__file__))
V = [l.split() for l in open(os.path.join(HERE, "vds2013_table4.txt")) if not l.startswith("#")]
V = [(int(r[0]), r[1]) + tuple(float(x) for x in r[2:]) for r in V]
V = [r for r in V if not (r[8] - r[7] > 0.6) and not (r[0] == 7 and r[1] == "19627")]          # the two rows Entry 45 dropped
shen_early = lambda m: 2.88e-6 * m ** 0.56
m_loc = 10 ** 11.2; r_loc = shen_early(m_loc); z_loc = 0.06
def anchor(law):                                 # the nearby galaxies are re-read with the law too (a 4 to 14% change at z = 0.06)
    da, dl = LAWS[law]; fA, fL2 = da(z_loc) / LAWS["standard"][0](z_loc), (dl(z_loc) / LAWS["standard"][1](z_loc)) ** 2
    gb = G * m_loc * fL2 / 2 / (r_loc * fA * 1e3) ** 2 * ACC; nu = float(rar(gb) / gb); ratio = (1 / 0.59) * fA / fL2
    return gb, nu, ratio, ratio / nu
gb_loc, nu_loc, ratio_loc, cal = anchor("standard")
print("\nSample C (Entry 45): %d quiescent galaxies. Nearby anchor: radius %.1f kpc, baryonic acceleration %.1f x a0, relation predicts dynamical/stellar %.2f, measured %.2f: calibration x%.2f" % (len(V), r_loc, gb_loc / GDAG, nu_loc, ratio_loc, cal))
res["quiescent"] = {"n": len(V), "anchor": dict(r=r_loc, gbar_over_a0=gb_loc / GDAG, nu=nu_loc, cal=cal)}
print("  law                           radius kpc  g_bar/a0   dynamical/stellar over nearby (Entry 45)   relation allows   left over")
for law in LAWS:
    da, dl = LAWS[law]; a = []; b = []; rr = []; gg = []; nn = []; gbl, nul, ratl, call = anchor(law)
    for ref, idn, z, re_, n_, sig, kn, lmd, lms in V:
        fA, fL2 = da(z) / LAWS["standard"][0](z), (dl(z) / LAWS["standard"][1](z)) ** 2
        r = re_ * fA; ms = 10 ** lms * fL2; md = 10 ** lmd * fA; gb = G * ms / 2 / (r * 1e3) ** 2 * ACC; nu = float(rar(gb) / gb)
        a.append(math.log10(md / ms / ratl)); b.append(math.log10(md / ms / nu / call)); rr.append(r); gg.append(gb / GDAG); nn.append(math.log10(nu / nul))
    a, b, nn = np.array(a), np.array(b), np.array(nn)
    q = dict(r=float(np.median(rr)), gbar_over_a0=float(np.median(gg)), same_mass=[float(a.mean()), float(a.std(ddof=1) / len(a) ** 0.5)], allowed=float(nn.mean()), same_acc=[float(b.mean()), float(b.std(ddof=1) / len(b) ** 0.5)])
    res["quiescent"][law] = q
    print("  %-28s  %5.1f      %5.1f      x%.2f (log %+.2f +- %.2f)                  x%.2f             x%.2f (log %+.2f +- %.2f)" % (law, q["r"], q["gbar_over_a0"], 10 ** a.mean(), a.mean(), q["same_mass"][1], 10 ** nn.mean(), 10 ** b.mean(), b.mean(), q["same_acc"][1]))
# the distance-free test
m = 10 ** ALP_ALL["log_mstar"][1]; loc = ms_local(m)
print("\nStar formation rate over stellar mass, all 117 ALPINE galaxies: median %.1f per Gyr; nearby star-forming sequence at that mass %.2f per Gyr: x%.0f" % (10 ** ALP_ALL["log_ssfr"][1] * 1e9, loc * 1e9, 10 ** ALP_ALL["log_ssfr"][1] / loc))
res["ssfr"] = dict(alpine_per_gyr=10 ** ALP_ALL["log_ssfr"][1] * 1e9, local_per_gyr=loc * 1e9, ratio=10 ** ALP_ALL["log_ssfr"][1] / loc)
for law in ("tired light",):
    fL2 = (LAWS[law][1](4.58) / LAWS["standard"][1](4.58)) ** 2; print("   with %s masses (x%.2f) the nearby sequence at that mass is %.2f per Gyr: x%.0f" % (law, fL2, ms_local(m * fL2) * 1e9, 10 ** ALP_ALL["log_ssfr"][1] / ms_local(m * fL2)))
    res["ssfr"]["ratio_tired_light"] = 10 ** ALP_ALL["log_ssfr"][1] / ms_local(m * fL2)
res["rar_curve"] = [[x, round(float(np.log10(rar(10 ** x))), 3)] for x in np.arange(-13.0, -8.49, 0.1)]
res["factors_z4p6"] = {law: [LAWS[law][0](4.58) / LAWS["standard"][0](4.58), (LAWS[law][1](4.58) / LAWS["standard"][1](4.58)) ** 2] for law in LAWS}
json.dump(res, open("local_relations_z5.json", "w"), indent=1, default=float)
