"""age_bias_referee.py - Workbench Entry 58
A 2026 paper (Nguyen, arXiv:2608.00806) reports that supernovae corrected for a progenitor-age bias (Son et al. 2025)
prefer the logarithmic law d_L = (c/H0)(1+z) ln(1+z) over flat LCDM. This script

 1. repeats the supernova comparison on Pantheon+ (1,590 light curves, diagonal errors) and on the DES five-year
    sample with the 2025 calibration (1,820 supernovae, binned in the browser run of des_dovekie_fits.js);
 2. shows that the correction is close to the difference between the two laws;
 3. tests the logarithmic law on the clustering ruler (DESI second release), which no supernova correction touches;
 4. asks what the correction does if expansion is left free: ruler + CMB summary + corrected supernovae, and the
    Hubble constant that follows;
 5. turns the question round: how much age correction do ruler and CMB allow the supernovae?
"""
import sys, os, json, math, time
HERE = os.path.dirname(os.path.abspath(__file__))
PP = os.path.join(HERE, "..", "pplus") if os.path.isdir(os.path.join(HERE, "..", "pplus")) else HERE
sys.path.insert(0, PP)
import numpy as np
from scipy import optimize, integrate
import cosmo as c
from cosmo import Model
from dark_energy_fit import sig

T0 = time.time()
out = {}
abc = lambda z: 0.183 * (1 - np.exp(-2.2 * np.asarray(z, float)))      # Nguyen's fit to Son et al.'s correction curve

# ---------------- supernova tables --------------------------------------------------------------------------
p = os.path.join(PP, "pantheon_plus_rows.txt")
if not os.path.exists(p): p = os.path.join(HERE, "pantheon_plus_rows.txt")
zr, mr, er = np.loadtxt(p).T
zu = np.unique(zr)
ZP = zu; MP = np.array([np.sum(mr[zr == a] / er[zr == a] ** 2) / np.sum(1 / er[zr == a] ** 2) for a in zu]); EP = np.array([er[zr == a].min() for a in zu])
# DES five-year, Dovekie calibration: redshift, distance modulus, error, number of supernovae; first 4 rows are the
# 197 nearby supernovae from other telescopes, the other 26 the 1,623 from the DES camera (des_dovekie_fits.js, DV.binned)
DES = np.array([[0.02762,35.4271,0.0198,40],[0.03435,35.9457,0.0141,65],[0.04758,36.6566,0.0142,51],[0.0679,37.4414,0.0152,41],
 [0.10001,38.2576,0.0433,10],[0.14763,39.227,0.0326,17],[0.18602,39.7396,0.0282,21],[0.21747,40.1285,0.0216,39],[0.2431,40.4151,0.0203,42],
 [0.27399,40.7218,0.0192,51],[0.30436,40.9774,0.0179,69],[0.33491,41.2251,0.0146,108],[0.36295,41.4047,0.0163,85],[0.39513,41.6204,0.0218,72],
 [0.42716,41.8384,0.0183,81],[0.45486,42.0227,0.0167,93],[0.48539,42.1702,0.0204,73],[0.51568,42.3075,0.0173,88],[0.54647,42.461,0.0164,101],
 [0.57419,42.6115,0.0177,96],[0.60531,42.7518,0.0181,97],[0.64253,42.92,0.0148,128],[0.67897,43.0467,0.0157,95],[0.72136,43.2125,0.0244,50],
 [0.75993,43.3637,0.024,54],[0.80241,43.4811,0.0225,54],[0.85153,43.6092,0.0364,38],[0.9079,43.8216,0.0413,26],[0.97735,43.9898,0.0386,21],
 [1.07141,44.2148,0.0426,14]])
ZD, MD, ED = DES[:, 0], DES[:, 1], DES[:, 2]
SAMPLES = {"Pantheon+": (ZP, MP, EP), "DES": (ZD, MD, ED), "DES, its own camera only": (ZD[4:], MD[4:], ED[4:])}
print("Pantheon+: %d supernovae after averaging repeats; DES: %d bins holding %d supernovae" % (len(ZP), len(ZD), DES[:, 3].sum()))

zg = np.linspace(0, 2.6, 5201)
def dc_lcdm(om, z):
    return np.interp(z, zg, integrate.cumulative_trapezoid(1 / np.sqrt(om * (1 + zg) ** 3 + 1 - om), zg, initial=0))
def chi_shape(z, m, e, dc, A=0.0):
    """chi2 of mu = 5 log10((1+z) dc) + constant, after taking A x the age correction off the magnitudes"""
    w = 1 / e ** 2; d = m - A * abc(z) - 5 * np.log10((1 + z) * dc); a = np.sum(w * d) / w.sum(); return float(np.sum(w * (d - a) ** 2))
def fit_om(z, m, e, A=0.0):
    f = lambda om: chi_shape(z, m, e, dc_lcdm(om, z), A); r = optimize.minimize_scalar(f, bounds=(0.05, 0.95), method="bounded", options={"xatol": 1e-6})
    lo = optimize.brentq(lambda o: f(o) - r.fun - 1, 0.03, r.x); hi = optimize.brentq(lambda o: f(o) - r.fun - 1, r.x, 0.97); return float(r.x), float((hi - lo) / 2), float(r.fun)

# ---------------- 1. supernovae alone ---------------------------------------------------------------------------
print("\n1. Supernovae alone: flat LCDM against the logarithmic law (chi2 of the law minus chi2 of LCDM)")
out["sn"] = {}
for nm, (z, m, e) in SAMPLES.items():
    row = {}
    for lab, A in (("before", 0.0), ("after", 1.0)):
        om, eom, c0 = fit_om(z, m, e, A); cl = chi_shape(z, m, e, np.log(1 + z), A); cz = chi_shape(z, m, e, z / np.sqrt(1 + z), A)
        row[lab] = dict(om=om, om_err=eom, chi_lcdm=c0, d_log=cl - c0, d_lerner=cz - c0)
        print("   %-26s %-6s the correction: LCDM matter share %.3f +- %.3f (chi2 %.1f); logarithmic law %+.1f; Lerner's straight line %+.1f" % (nm, lab, om, eom, c0, cl - c0, cz - c0))
    # the amplitude the logarithmic law itself asks for
    w = 1 / e ** 2; y = m - 5 * np.log10((1 + z) * np.log(1 + z)); x = abc(z); S1, Sx, Sy, Sxx, Sxy = w.sum(), (w * x).sum(), (w * y).sum(), (w * x * x).sum(), (w * x * y).sum(); den = S1 * Sxx - Sx ** 2
    row["A_log"] = [float((S1 * Sxy - Sx * Sy) / den), float(math.sqrt(S1 / den))]
    print("      for the logarithmic law to fit, the correction must be %.2f +- %.2f of the published one" % tuple(row["A_log"]))
    out["sn"][nm] = row
print("   (unbinned DES in the browser run: logarithmic law +208.6 before, -12.1 after; matter share 0.349 and 0.512.")
print("    The paper, with full error matrices: -11.9 for Pantheon+, -2.0 for DES after; +197.5 for the two together before.)")
j = sum(out["sn"][k]["after"]["d_log"] for k in ("Pantheon+", "DES")); jb = sum(out["sn"][k]["before"]["d_log"] for k in ("Pantheon+", "DES"))
print("   two samples added: %+.1f before, %+.1f after" % (jb, j)); out["sn_joint"] = [jb, j]

# ---------------- 2. the correction against the difference between the laws --------------------------------------
print("\n2. The correction, and the standard curve minus the logarithmic law (both set to zero at z = 0.03)")
zz = np.array([0.03, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0])
rows = []
for om in (0.30, 0.315, 0.35):
    d = 5 * np.log10(dc_lcdm(om, zz) / np.log(1 + zz)); d -= d[0]; rows.append(d)
a = abc(zz) - abc(0.03)
for i, z in enumerate(zz):
    print("   z %.2f: correction %.3f | LCDM minus logarithmic law: matter 0.30 %.3f, 0.315 %.3f, 0.35 %.3f" % (z, a[i], rows[0][i], rows[1][i], rows[2][i]))
out["degeneracy"] = dict(z=zz.tolist(), abc=a.tolist(), d030=rows[0].tolist(), d0315=rows[1].tolist(), d035=rows[2].tolist())
zf = np.linspace(0.03, 1.1, 200)
for om in (0.30, 0.315, 0.35):
    d = 5 * np.log10(dc_lcdm(om, zf) / np.log(1 + zf)); d -= d[0]; r = d - (abc(zf) - abc(0.03)); r -= r.mean()
    print("   matter share %.3f: the two differ by at most %.3f mag between z = 0.03 and 1.1, once levelled" % (om, np.abs(r).max()))
out["degeneracy_max"] = float(np.abs(r).max())
def gap(om, z): return float(5 * np.log10(dc_lcdm(om, np.array([z]))[0] / np.log(1 + z)) - 5 * np.log10(dc_lcdm(om, np.array([0.03]))[0] / np.log(1.03)))
curve = [[round(float(z), 3), round(float(abc(z) - abc(0.03)), 4), round(gap(0.30, z), 4), round(gap(0.35, z), 4)] for z in np.r_[0.03, np.arange(0.05, 2.01, 0.05)]]
out["curve"] = curve

# ---------------- 3. the ruler ------------------------------------------------------------------------------------
print("\n3. The clustering ruler (DESI second release, 13 numbers)")
ZB, ZV, YB, CBI, CB = c.ZB, c.ZV, c.YB, c.CBI, c.CB
def ruler_chi(dm, dh, K):
    t = K * np.r_[(ZV * dm(ZV) ** 2 * dh(ZV)) ** (1 / 3), dm(ZB), dh(ZB)]; r = YB - t; return float(r @ CBI @ r)
def ruler_fit(dm, dh):
    r = optimize.minimize_scalar(lambda K: ruler_chi(dm, dh, K), bounds=(1, 100), method="bounded", options={"xatol": 1e-7}); return float(r.fun), float(r.x)
clog, Klog = ruler_fit(lambda z: np.log(1 + z), lambda z: 1 / (1 + z))
def lcdm_ruler(om):
    return ruler_fit(lambda z: dc_lcdm(om, np.atleast_1d(z)) if np.ndim(z) else dc_lcdm(om, np.array([z]))[0], lambda z: 1 / np.sqrt(om * (1 + z) ** 3 + 1 - om))
r = optimize.minimize_scalar(lambda om: lcdm_ruler(om)[0], bounds=(0.1, 0.6), method="bounded"); omb, cb = float(r.x), float(r.fun)
lo = optimize.brentq(lambda o: lcdm_ruler(o)[0] - cb - 1, 0.1, omb); hi = optimize.brentq(lambda o: lcdm_ruler(o)[0] - cb - 1, omb, 0.6)
print("   logarithmic law, one free scale: chi2 %.1f for 13 numbers;  flat LCDM, matter share %.3f +- %.3f: chi2 %.1f" % (clog, omb, (hi - lo) / 2, cb))
out["ruler"] = dict(chi_log=clog, chi_lcdm=cb, om=omb, om_err=(hi - lo) / 2)
# the shape ratio alone: distance across over distance along, which needs no ruler length at all
ap = []
for i, z in enumerate(ZB):
    dm, dh = YB[1 + i], YB[7 + i]; sm, sh, cv = math.sqrt(CB[1 + i, 1 + i]), math.sqrt(CB[7 + i, 7 + i]), CB[1 + i, 7 + i]
    F = dm / dh; eF = F * math.sqrt((sm / dm) ** 2 + (sh / dh) ** 2 - 2 * cv / (dm * dh)); Flog = (1 + z) * math.log(1 + z)
    Fl = dc_lcdm(omb, np.array([z]))[0] * math.sqrt(omb * (1 + z) ** 3 + 1 - omb); ap.append([float(z), F, eF, Flog, Fl, (F - Flog) / eF, (F - Fl) / eF])
    print("   z %.3f: across/along %.3f +- %.3f | logarithmic law %.3f (%+.1f sigma) | LCDM %.3f (%+.1f sigma)" % (z, F, eF, Flog, (F - Flog) / eF, Fl, (F - Fl) / eF))
c_ap_log = sum(a[5] ** 2 for a in ap); c_ap_l = sum(a[6] ** 2 for a in ap)
print("   shape ratio alone, six redshifts, nothing fitted: logarithmic law chi2 %.1f; LCDM at the ruler's matter share %.1f" % (c_ap_log, c_ap_l))
out["ap"] = dict(rows=ap, chi_log=c_ap_log, chi_lcdm=c_ap_l)
# one tracer left out at a time
drop = {}
for i, nm in enumerate(c.BAO_NAMES):
    k = c.bao_mask(i); Ci = np.linalg.inv(CB[np.ix_(k, k)])
    def f(K, k=k, Ci=Ci):
        t = K * np.r_[(ZV * np.log(1 + ZV) ** 2 / (1 + ZV)) ** (1 / 3), np.log(1 + ZB), 1 / (1 + ZB)]; r = (YB - t)[k]; return float(r @ Ci @ r)
    drop[nm] = float(optimize.minimize_scalar(f, bounds=(1, 100), method="bounded").fun)
print("   logarithmic law with one tracer left out: " + "; ".join("%s %.0f" % (a, b) for a, b in drop.items()))
out["ruler_drop"] = drop
# the same for flat LCDM without the highest-redshift tracer (the Lyman-alpha forest at z = 2.33)
k = c.bao_mask(6); Ci = np.linalg.inv(CB[np.ix_(k, k)])
def lcdm_nolya(om):
    dmv = dc_lcdm(om, np.r_[ZV, ZB]); E = lambda z: np.sqrt(om * (1 + z) ** 3 + 1 - om); shape = np.r_[(ZV * dmv[0] ** 2 / E(ZV)) ** (1 / 3), dmv[1:], 1 / E(ZB)]
    return float(optimize.minimize_scalar(lambda K: float(((YB - K * shape)[k]) @ Ci @ ((YB - K * shape)[k])), bounds=(1, 100), method="bounded").fun)
rr = optimize.minimize_scalar(lcdm_nolya, bounds=(0.1, 0.6), method="bounded")
print("   without the z = 2.33 tracer (11 numbers): logarithmic law %.1f; flat LCDM %.1f (matter share %.3f)" % (drop["Lya 2.33"], rr.fun, rr.x))
out["ruler_nolya"] = dict(log=drop["Lya 2.33"], lcdm=float(rr.fun), om=float(rr.x))
# what the law predicts at z = 2.33 once its scale is set by the six tracers below it
Kl = optimize.minimize_scalar(lambda K: float(((YB - K * np.r_[(ZV * np.log(1 + ZV) ** 2 / (1 + ZV)) ** (1 / 3), np.log(1 + ZB), 1 / (1 + ZB)])[k]) @ Ci @ ((YB - K * np.r_[(ZV * np.log(1 + ZV) ** 2 / (1 + ZV)) ** (1 / 3), np.log(1 + ZB), 1 / (1 + ZB)])[k])), bounds=(1, 100), method="bounded").x
print("   scale from the six lower tracers: at z = 2.33 the law predicts across %.2f (measured %.2f +- %.2f), along %.3f (measured %.3f +- %.3f)"
      % (Kl * math.log(3.33), YB[6], math.sqrt(CB[6, 6]), Kl / 3.33, YB[12], math.sqrt(CB[12, 12])))
pred = Kl * np.r_[(ZV * np.log(1 + ZV) ** 2 / (1 + ZV)) ** (1 / 3), np.log(1 + ZB), 1 / (1 + ZB)]; sd = np.sqrt(np.diag(CB))
out["ruler_ratio"] = dict(z=[float(ZV)] + ZB.tolist() + ZB.tolist(), kind=["both"] + ["across"] * 6 + ["along"] * 6, ratio=(YB / pred).tolist(), err=(sd / pred).tolist())
print("   measured over the law's prediction: " + ", ".join("%.3f+-%.3f" % (a, b) for a, b in zip(YB / pred, sd / pred)))
# eBOSS (du Mas des Bourboux et al. 2020), an earlier survey at the same redshift: along 8.99 +- 0.19, across 37.5 +- 1.1
print("   earlier SDSS measurement at z = 2.33: along 8.99 +- 0.19 is %.1f sigma from the law's %.3f; across 37.5 +- 1.1 is %.1f sigma" % ((Kl / 3.33 - 8.99) / 0.19, Kl / 3.33, (Kl * math.log(3.33) - 37.5) / 1.1))
out["eboss"] = dict(dh=8.99, dh_err=0.19, sigma=float((Kl / 3.33 - 8.99) / 0.19))
out["lya_pred"] = dict(dm=float(Kl * math.log(3.33)), dm_obs=float(YB[6]), dm_err=float(math.sqrt(CB[6, 6])), dh=float(Kl / 3.33), dh_obs=float(YB[12]), dh_err=float(math.sqrt(CB[12, 12])))
# ruler and supernovae together under the logarithmic law, after the correction
for nm in ("Pantheon+", "DES"):
    z, m, e = SAMPLES[nm]; tot_log = clog + chi_shape(z, m, e, np.log(1 + z), 1.0)
    f = lambda om: lcdm_ruler(om)[0] + chi_shape(z, m, e, dc_lcdm(om, z), 1.0); rr = optimize.minimize_scalar(f, bounds=(0.1, 0.7), method="bounded")
    print("   ruler + corrected %s: logarithmic law %.1f; flat LCDM %.1f (matter share %.3f): law minus LCDM %+.1f" % (nm, tot_log, rr.fun, rr.x, tot_log - rr.fun))
    out["ruler_" + nm] = dict(log=tot_log, lcdm=float(rr.fun), om=float(rr.x))

# ---------------- 4. expansion left free: ruler + CMB summary + supernovae ---------------------------------------
print("\n4. Ruler + CMB summary + supernovae, with and without the correction")
def sn_chi(z, m, e, A):
    w = 1 / e ** 2; mm = m - A * abc(z)
    def f(mod):
        r = mm - 5 * np.log10((1 + z) * mod.DM(z)); a = np.sum(w * r) / w.sum(); return float(np.sum(w * (r - a) ** 2))
    return f
def unpack(x, kind):
    ob, obc, h = x[0] * 0.02223, x[1] * 0.142, x[2] * 0.68; w0, wa = (-1.0, 0.0) if kind == "lcdm" else (x[3], x[4]); return ob, obc, h, w0, wa
def total(x, kind, parts, hfix=None):
    ob, obc, h, w0, wa = unpack(x, kind)
    if hfix is not None: h = hfix
    if not (0.015 < ob < 0.03 and 0.08 < obc < 0.25 and 0.5 < h < 0.9 and -3 < w0 < 1 and -4 < wa < 2 and w0 + wa < 0): return 1e9
    mod = Model(ob, obc, h, w0, wa)
    if mod.ode <= 0: return 1e9
    return sum(pf(mod) for pf in parts)
def best(kind, parts, starts=None, hfix=None):
    seeds = starts or ([[1, 1, 1]] if kind == "lcdm" else [[1, 1, 1, -1.0, 0.0], [1, 1.0, 0.97, -0.7, -1.0], [1, 1.0, 0.94, -0.4, -1.8], [1, 1.0, 0.92, -0.1, -2.6]])
    b = None
    for s in seeds:
        r = optimize.minimize(total, s, args=(kind, parts, hfix), method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-5, maxiter=5000, maxfev=5000))
        r = optimize.minimize(total, r.x, args=(kind, parts, hfix), method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-6, maxiter=5000, maxfev=5000))
        if b is None or r.fun < b.fun: b = r
    return b
def h_profile(kind, parts, b):
    """68% range of the Hubble constant by stepping it away from the best fit"""
    h0 = unpack(b.x, kind)[2]; res = {}
    for sgn in (-1, 1):
        x = b.x.copy(); h = h0; prev = (h0, 0.0)
        for k in range(40):
            h = h0 + sgn * 0.004 * (k + 1); r = optimize.minimize(total, x, args=(kind, parts, h), method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-5, maxiter=3000, maxfev=3000)); x = r.x; d = r.fun - b.fun
            if d >= 1: res[sgn] = prev[0] + (h - prev[0]) * (1 - prev[1]) / (d - prev[1]); break
            prev = (h, d)
        else: res[sgn] = h
    return 100 * h0, 100 * (res[1] - h0), 100 * (h0 - res[-1])
bao = lambda mod: c.chi_bao(mod); cmb = lambda mod: c.chi_cmb(mod)
def q0(P):
    ob, obc, h, w0, wa = P; ode = (h * h - obc - c.OG_H2 - c.OMNU0) / h ** 2; return 0.5 * (1 + 3 * w0 * ode)
out["free"] = {}
bL0 = best("lcdm", [bao, cmb]); bW0 = best("w0wa", [bao, cmb]); P = unpack(bW0.x, "w0wa")
print("   ruler + CMB alone: constant dark energy chi2 %.2f, H0 %.2f; changing chi2 %.2f (gain %.1f = %.1f sigma), w0 %+.2f, wa %+.2f, H0 %.1f, q0 %+.2f"
      % (bL0.fun, 100 * unpack(bL0.x, "lcdm")[2], bW0.fun, bL0.fun - bW0.fun, sig(bL0.fun - bW0.fun), P[3], P[4], 100 * P[2], q0(P)))
H, up, dn = h_profile("w0wa", [bao, cmb], bW0); print("      Hubble constant under changing dark energy, ruler + CMB alone: %.1f +%.1f -%.1f" % (H, up, dn))
out["free"]["none_H"] = [H, up, dn]
out["free"]["none"] = dict(chi_l=float(bL0.fun), chi_w=float(bW0.fun), gain=float(bL0.fun - bW0.fun), w0=P[3], wa=P[4], H0=100 * P[2], q0=q0(P), H0_l=100 * unpack(bL0.x, "lcdm")[2])
for nm in ("Pantheon+", "DES"):
    z, m, e = SAMPLES[nm]
    for lab, A in (("before", 0.0), ("after", 1.0)):
        parts = [bao, cmb, sn_chi(z, m, e, A)]; bL = best("lcdm", parts); bW = best("w0wa", parts); P = unpack(bW.x, "w0wa"); PL = unpack(bL.x, "lcdm")
        H, up, dn = h_profile("w0wa", parts, bW); HL, upL, dnL = h_profile("lcdm", parts, bL); g = bL.fun - bW.fun
        print("   %-10s %-6s: constant chi2 %.1f, H0 %.2f +%.2f -%.2f | changing chi2 %.1f, gain %.1f = %.1f sigma, w0 %+.3f, wa %+.3f, H0 %.2f +%.2f -%.2f, q0 %+.2f"
              % (nm, lab, bL.fun, HL, upL, dnL, bW.fun, g, sig(g), P[3], P[4], H, up, dn, q0(P)))
        out["free"][nm + " " + lab] = dict(chi_l=float(bL.fun), chi_w=float(bW.fun), gain=float(g), sigma=sig(g), w0=P[3], wa=P[4], H0=H, H0_up=up, H0_dn=dn, q0=q0(P), H0_l=HL, H0_l_up=upL, H0_l_dn=dnL,
                                           om_l=float((PL[1] + c.OMNU0) / PL[2] ** 2), sn_l=float(parts[2](Model(*PL))), sn_w=float(parts[2](Model(*P))))
    sys.stdout.flush()
# the referee of Entry 56
H10, E10 = 71.40, 1.20
print("   against ten one-step measurements of the Hubble constant, %.1f +- %.1f:" % (H10, E10))
for k, v in out["free"].items():
    if k in ("none", "none_H"): continue
    eh = (v["H0_up"] + v["H0_dn"]) / 2; s = (H10 - v["H0"]) / math.hypot(E10, eh); v["referee"] = s
    print("      %-18s changing dark energy gives %.1f +- %.1f: %.1f sigma below" % (k, v["H0"], eh, s))

# ---------------- 5. how much correction do ruler and CMB allow? ------------------------------------------------
print("\n5. The amplitude of the correction, fitted with everything else (1 = as published, 0 = none)")
out["A"] = {}
def A_profile(kind, z, m, e, with_cmb=True):
    grid = np.arange(-1.0, 2.01, 0.25); vals = []; x = None
    for A in grid:
        parts = [bao] + ([cmb] if with_cmb else []) + [sn_chi(z, m, e, A)]
        if not with_cmb: parts.append(lambda mod: ((mod.ob - 0.02223) / 0.0005) ** 2)      # ruler alone cannot fix the baryons
        b = best(kind, parts, starts=None if x is None else [x]) if x is not None else best(kind, parts); x = b.x; vals.append(b.fun)
    vals = np.array(vals); pc = np.polyfit(grid, vals, 2); A0 = -pc[1] / (2 * pc[0]); eA = 1 / math.sqrt(pc[0]); return float(A0), float(eA), vals.tolist()
for nm in ("Pantheon+", "DES"):
    z, m, e = SAMPLES[nm]
    for kind, lab, wc in (("lcdm", "constant dark energy, ruler + supernovae", False), ("lcdm", "constant dark energy, ruler + CMB + supernovae", True), ("w0wa", "changing dark energy, ruler + CMB + supernovae", True)):
        A0, eA, vals = A_profile(kind, z, m, e, wc); print("   %-10s %-48s A = %+.2f +- %.2f" % (nm, lab, A0, eA)); out["A"][nm + " | " + lab] = [A0, eA]; sys.stdout.flush()

# the same with the ten one-step Hubble constants added as a measurement
href = lambda mod: ((100 * mod.h - H10) / E10) ** 2
def A_profile_h(z, m, e):
    grid = np.arange(-1.5, 1.51, 0.25); vals = []; x = None
    for A in grid:
        parts = [bao, cmb, sn_chi(z, m, e, A), href]; b = best("w0wa", parts, starts=None if x is None else [x, [1, 1.0, 0.97, -0.7, -1.0], [1, 1.0, 1.0, -1.0, 0.0]]); x = b.x; vals.append(b.fun)
    vals = np.array(vals); pc = np.polyfit(grid, vals, 2); return float(-pc[1] / (2 * pc[0])), float(1 / math.sqrt(pc[0])), float(np.interp(1.0, grid, vals) - vals.min()), float(np.interp(0.0, grid, vals) - vals.min())
for nm in ("Pantheon+", "DES"):
    z, m, e = SAMPLES[nm]; A0, eA, d1, d0 = A_profile_h(z, m, e)
    print("   %-10s changing dark energy, ruler + CMB + supernovae + the ten Hubble constants: A = %+.2f +- %.2f (chi2 rises by %.1f at A = 1, %.1f at A = 0)" % (nm, A0, eA, d1, d0)); out["A"][nm + " | with H0"] = [A0, eA, d1, d0]

out["seconds"] = time.time() - T0
json.dump(out, open(os.path.join(HERE, "age_bias_referee.json"), "w"), indent=1)
print("\ndone in %.0f s" % (time.time() - T0))
