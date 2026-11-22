"""ccc_tl_test.py - Workbench Entry 65
The CCC+TL model (Gupta 2023, MNRAS 524, 3385): part of the redshift comes from expansion with "covarying coupling
constants", part from tired light:  1 + z = (1 + z_x)(1 + z_t).  It has two fitted numbers, both from supernovae:
H_c = 59.51 km/s/Mpc and alpha / H_c = -0.7997.  Everything else follows:
   f(t) = exp[alpha (t - t0)]                      the factor by which the constants differed
   (H + alpha)^2 = (H_c + alpha)^2 (1+z_x)^3 / f   flat, matter only
   which integrates to  (1+z_x)^(-3/2) = C y^3 + K y,  y = f^(-1/2),  K = 1.5 (H_c + alpha)/alpha,  C = 1 - K
   H_t = (H_c + alpha)/2 x (3 + alpha/H_c) = 13.11  tired-light rate, the paper's choice; H0 = H_c + H_t = 72.62
   distance d_p = c Int dz_x / H = (c / H_t) ln(1 + z_t)      which fixes the split at each z
   d_L = d_p (1+z_x) sqrt(1+z_t),   d_A = d_p f / (1+z_x)
With no free number left, the model makes fixed predictions for tests that do not use supernova brightness:
  1. the shape test: across over along for the clustering ruler, D_M / D_H = d_p (1+z) [H(z_x)/(1+z_x) + H_t] / c   (f cancels)
  2. the ruler itself, with one free length
  3. supernova clocks: only the expansion part stretches time
  4. line widths of quiescent galaxies (Entries 45, 62, 63), read at face value with today's constants
"""
import json, math, os
import numpy as np
from scipy import integrate, optimize
HERE = os.path.dirname(os.path.abspath(__file__))
C = 299792.458
HC, AR = 59.51, -0.7997
AL = AR * HC; K = 1.5 * (1 + AR) / AR; CC = 1 - K; HT = (HC + AL) / 2 * (3 + AR); H0 = HC + HT
def yof(zx):                                   # y = f^(-1/2) from the cubic
    t = (1 + zx) ** -1.5; ymin = math.sqrt(-K / CC)
    return optimize.brentq(lambda y: CC * y ** 3 + K * y - t, ymin, 1.0 + 1e-12)
Hx = lambda zx: (HC + AL) * (1 + zx) ** 1.5 * yof(zx) - AL
_zx = np.concatenate([np.linspace(0, 1, 2001), np.linspace(1, 30, 5801)[1:]]); _dp = integrate.cumulative_trapezoid([C / Hx(z) for z in _zx], _zx, initial=0)
dp_x = lambda zx: float(np.interp(zx, _zx, _dp))
def split(z):                                  # expansion part z_x for an observed z
    g = lambda zx: math.log(1 + zx) + HT * dp_x(zx) / C - math.log(1 + z)
    return optimize.brentq(g, 0, z) if z > 0 else 0.0
def model(z):
    zx = split(z); zt = (1 + z) / (1 + zx) - 1; d = dp_x(zx); y = yof(zx); f = y ** -2
    return dict(z=z, zx=zx, zt=zt, dp=d, f=f, dL=d * (1 + zx) * math.sqrt(1 + zt), dA=d * f / (1 + zx), Heff=(1 + z) * (Hx(zx) / (1 + zx) + HT))
# standard model
def lcdm(om, h0=70.0):
    zg = np.linspace(0, 30, 30001); dc = integrate.cumulative_trapezoid(1 / np.sqrt(om * (1 + zg) ** 3 + 1 - om), zg, initial=0) * C / h0
    DM = lambda z: float(np.interp(z, zg, dc)); return DM, (lambda z: C / h0 / math.sqrt(om * (1 + z) ** 3 + 1 - om))
out = {"params": dict(Hc=HC, alpha_over_Hc=AR, alpha=AL, Ht=HT, H0=H0, K=K, C=CC)}

# ---- checks against the paper -------------------------------------------------------------------------------------------
age = -2 / AL * (-math.log(math.sqrt(-K / CC))) * 977.8
print("H_t = %.2f, H0 = %.2f km/s/Mpc; age %.2f Gyr (paper: 13.11, 72.62, 26.70); f at the beginning %.2f" % (HT, H0, age, CC / -K))
DMs, DHs = lcdm(0.3508, 72.84)
print("size against the paper's standard model (matter 0.3508):", "  ".join("z=%g: x%.1f" % (z, model(z)["dA"] / (DMs(z) / (1 + z))) for z in (1, 2, 5, 10, 15)), " (paper: 5.6, 9.5, 12.8 at 5, 10, 15)")
print("tired-light over expansion part:", "  ".join("z=%g: %.0f%%" % (z, 100 * model(z)["zt"] / model(z)["zx"]) for z in (0.01, 1, 10, 20)), " (paper: 22%, -, 20%, 18%)")
out["checks"] = dict(age=age, size_ratio={str(z): model(z)["dA"] / (DMs(z) / (1 + z)) for z in (1, 2, 5, 10, 15)}, zt_over_zx={str(z): model(z)["zt"] / model(z)["zx"] for z in (0.01, 1, 10, 20)})

# ---- 0. supernovae: does this implementation fit Pantheon+ as the paper says? -------------------------------------------
P = np.loadtxt(os.path.join(HERE, "pantheon_plus_rows.txt")); zs, ms, es = P[:, 0], P[:, 1], P[:, 2]
def sn_chi(mu):
    w = 1 / es ** 2; off = np.sum(w * (ms - mu)) / np.sum(w); return float(np.sum(w * (ms - mu - off) ** 2))
mu_c = np.array([5 * math.log10(model(z)["dL"]) for z in zs])
def mu_l(om): DM, _ = lcdm(om); return np.array([5 * math.log10(DM(z) * (1 + z)) for z in zs])
r = optimize.minimize_scalar(lambda om: sn_chi(mu_l(om)), bounds=(0.1, 0.6), method="bounded")
print("\n0. Pantheon+ (%d supernovae, diagonal errors, zero point free): CCC+TL chi2 %.1f; standard %.1f at matter share %.3f" % (len(zs), sn_chi(mu_c), r.fun, r.x))
out["sn"] = dict(n=int(len(zs)), ccc=sn_chi(mu_c), lcdm=float(r.fun), om=float(r.x))

# ---- 1 and 2. the clustering ruler (DESI second release) ---------------------------------------------------------------
DR2 = [('V',0.295,7.942,0.075),('MH',0.510,13.588,0.167,21.863,0.425,-0.459),('MH',0.706,17.351,0.177,19.455,0.330,-0.404),('MH',0.934,21.576,0.152,17.641,0.193,-0.416),
       ('MH',1.321,27.601,0.318,14.176,0.221,-0.434),('MH',1.484,30.512,0.760,12.817,0.516,-0.500),('MH',2.330,38.988,0.531,8.632,0.101,-0.431)]
ccc = ((lambda z: model(z)["dp"]), (lambda z: C / model(z)["Heff"]))
def chi(data, mod, Ks):
    dm, dh = mod; c2 = 0.0
    for d in data:
        if d[0] == 'V': c2 += ((d[2] - Ks * (d[1] * dm(d[1]) ** 2 * dh(d[1])) ** (1 / 3)) / d[3]) ** 2
        else:
            a = (d[2] - Ks * dm(d[1])) / d[3]; b = (d[4] - Ks * dh(d[1])) / d[5]; rr = d[6]; c2 += (a * a + b * b - 2 * rr * a * b) / (1 - rr * rr)
    return c2
def best(data, mod):
    r = optimize.minimize_scalar(lambda k: chi(data, mod, k), bounds=(1e-4, 1.0), method="bounded", options={"xatol": 1e-10}); return float(r.fun), float(r.x)
c_ccc, k_ccc = best(DR2, ccc)
rl = optimize.minimize_scalar(lambda om: best(DR2, lcdm(om))[0], bounds=(0.1, 0.7), method="bounded"); c_l, k_l = best(DR2, lcdm(rl.x))
print("\n2. Ruler with one free length, 13 DESI numbers: CCC+TL chi2 %.1f (ruler %.1f Mpc); standard %.1f (matter share %.3f)" % (c_ccc, 1 / k_ccc, c_l, rl.x))
print("1. Shape only, across over along (no ruler length, no constants):")
rows = []; c_shape = 0.0; c_shape_l = 0.0; DMl, DHl = lcdm(rl.x)
for d in DR2:
    if d[0] != 'MH': continue
    z = d[1]; F = d[2] / d[4]; eF = F * math.sqrt((d[3] / d[2]) ** 2 + (d[5] / d[4]) ** 2 - 2 * d[6] * d[3] * d[5] / (d[2] * d[4]))
    m = model(z); pc = m["dp"] * m["Heff"] / C; pl = DMl(z) / DHl(z); c_shape += ((F - pc) / eF) ** 2; c_shape_l += ((F - pl) / eF) ** 2
    rows.append(dict(z=z, obs=F, err=eF, ccc=pc, lcdm=pl, pull_ccc=(F - pc) / eF, pull_lcdm=(F - pl) / eF, zx=m["zx"], zt=m["zt"]))
    print("   z = %.3f: measured %.3f +- %.3f; CCC+TL %.3f (%+.1f sigma); standard %.3f (%+.1f sigma)" % (z, F, eF, pc, (F - pc) / eF, pl, (F - pl) / eF))
print("   six ratios: CCC+TL chi2 %.1f with nothing free; standard %.1f with the matter share fitted to the ruler" % (c_shape, c_shape_l))
# the 2026 full-shape forest point (Entry 58): along 8.600 +- 0.066, across 39.32 +- 0.33, correlation taken as -0.43
F26 = 39.32 / 8.600; e26 = F26 * math.sqrt((0.33 / 39.32) ** 2 + (0.066 / 8.600) ** 2 + 2 * 0.43 * 0.33 * 0.066 / (39.32 * 8.600)); m = model(2.33)
print("   DESI 2026 forest, z = 2.33: %.3f +- %.3f; CCC+TL %.3f (%+.1f sigma); standard %.3f (%+.1f sigma)" % (F26, e26, m["dp"] * m["Heff"] / C, (F26 - m["dp"] * m["Heff"] / C) / e26, DMl(2.33) / DHl(2.33), (F26 - DMl(2.33) / DHl(2.33)) / e26))
out["ruler"] = dict(ccc=c_ccc, ruler_mpc=1 / k_ccc, lcdm=c_l, om=float(rl.x)); out["shape"] = dict(rows=rows, ccc=c_shape, lcdm=c_shape_l, forest2026=dict(obs=F26, err=e26, ccc=m["dp"] * m["Heff"] / C, lcdm=DMl(2.33) / DHl(2.33)))
# what share of tired light would the shape allow if H_t were free (alpha/H_c fixed)? -> see Entry 15 for the general case

# ---- 3. clocks -----------------------------------------------------------------------------------------------------------
print("\n3. Clocks: stretch exponent b in (1+z)^b if only the expansion part stretches time")
bs = {z: math.log(1 + model(z)["zx"]) / math.log(1 + z) for z in (0.1, 0.3, 0.5, 0.8, 1.0, 1.2)}
print("   " + "  ".join("z=%.1f: %.3f" % kv for kv in bs.items()) + "   measured (DES, 1,504 supernovae): 1.003 +- 0.005 +- 0.010")
bm = float(np.mean([bs[z] for z in (0.3, 0.5, 0.8, 1.0)])); print("   mean over z = 0.3-1.0: %.3f, which is %.0f sigma from 1.003" % (bm, (1.003 - bm) / math.hypot(0.005, 0.010)))
out["clocks"] = dict(b=bs, mean=bm, sigma=(1.003 - bm) / math.hypot(0.005, 0.010))

# ---- 4. line widths of quiescent galaxies, at face value -----------------------------------------------------------------
G = 4.30091e-3; ACC = 1e6 / 3.0857e16; A0 = 1.20e-10
nu = lambda gb: 1 / (1 - math.exp(-math.sqrt(gb / A0)))
DM70, _ = lcdm(0.3, 70.0); dA_s = lambda z: DM70(z) / (1 + z); dL_s = lambda z: DM70(z) * (1 + z)
def factors(z): m = model(z); return m["dA"] / dA_s(z) * (H0 / 70.0), (m["dL"] / dL_s(z) * (H0 / 70.0)) ** 2, m["f"]   # distances put on H0 = 70 like the catalogues
m_loc = 10 ** 11.2; r_loc = 2.88e-6 * m_loc ** 0.56
def anchor(ccctl):
    fA, fL2, f = factors(0.06) if ccctl else (1, 1, 1); gb = G * m_loc * fL2 / 2 / (r_loc * fA * 1e3) ** 2 * ACC; ratio = (1 / 0.59) * fA / fL2; return nu(gb), ratio / nu(gb), ratio
V = [l.split() for l in open(os.path.join(HERE, "vds2013_table4.txt")) if not l.startswith("#")]
V = [(float(r[2]), float(r[3]), float(r[7]), float(r[8])) for r in V if not (float(r[8]) - float(r[7]) > 0.6) and not (r[0] == "7" and r[1] == "19627")]   # z, r_e, log Mdyn, log M*
Z3 = [(2.9834,1.08,11.31,11.22),(3.0094,0.95,11.33,11.25),(2.7988,0.89,11.02,10.98),(3.1313,2.26,11.28,11.40),(3.1326,7.54,12.00,11.67),(3.4520,1.85,11.49,11.52),(3.4919,1.71,11.54,11.47),(3.4892,1.71,11.49,11.26),
      (3.715,0.66,10.86,11.06),(3.219,2.40,11.41,11.31),(3.239,0.63,10.56,10.99),(3.434,0.61,10.85,10.99),(3.352,1.0,11.34,11.30),(4.658,0.215,10.30,10.58)]
print("\n4. Quiescent galaxies, dynamical over stellar mass against nearby galaxies at the same acceleration (as Entries 62 and 63)")
out["widths"] = {}
for lab, S in (("71 at z = 0.8-2.2", V), ("14 at z = 2.8-4.7", Z3)):
    for ccctl in (False, True):
        nul, cal, ratl = anchor(ccctl); b = []; ff = []; FF = []
        for z, re_, lmd, lms in S:
            fA, fL2, f = factors(z) if ccctl else (1, 1, 1); r = re_ * fA; ms = 10 ** lms * fL2; md = 10 ** lmd * fA; gb = G * ms / 2 / (r * 1e3) ** 2 * ACC
            b.append(math.log10(md / ms / nu(gb) / cal)); ff.append(f); FF.append(fA / fL2)
        b = np.array(b); print("   %-18s %-9s left over x%.2f (log %+.2f +- %.2f); size-over-luminosity factor x%.1f; the model's f there %.2f" % (lab, "CCC+TL" if ccctl else "standard", 10 ** b.mean(), b.mean(), b.std(ddof=1) / len(b) ** 0.5, np.median(FF), np.median(ff)))
        out["widths"].setdefault(lab, {})["ccctl" if ccctl else "standard"] = dict(mean=float(b.mean()), err=float(b.std(ddof=1) / len(b) ** 0.5), F=float(np.median(FF)), f=float(np.median(ff)))
out["table"] = [dict(z=z, **{k: model(z)[k] for k in ("zx", "zt", "f", "dp", "dL", "dA", "Heff")}) for z in (0.1, 0.5, 1, 1.3, 2.33, 3.3, 5, 10)]
print("\n   z      z_x     z_t     f      d_A [Mpc]  d_A / standard(0.3, H0=70)")
for t in out["table"]: print("   %-5g  %.3f   %.3f   %.2f   %8.0f   x%.2f" % (t["z"], t["zx"], t["zt"], t["f"], t["dA"], t["dA"] / dA_s(t["z"])))
json.dump(out, open(os.path.join(HERE, "ccc_tl_test.json"), "w"), indent=1, default=float)

# ---- 5. is it the tired-light part or the expansion law? A scan of the family: alpha/H_c and H_t/H_c both free ----------
def family(ar, ht):                             # distances in units of c/H_c; returns functions of observed z
    K_ = 1.5 * (1 + ar) / ar if ar != 0 else None
    def y_(zx):
        if ar == 0: return 1.0
        t = (1 + zx) ** -1.5; C_ = 1 - K_; return optimize.brentq(lambda y: C_ * y ** 3 + K_ * y - t, math.sqrt(-K_ / C_), 1.0 + 1e-12)
    hx = lambda zx: (1 + ar) * (1 + zx) ** 1.5 * y_(zx) - ar
    zz = np.concatenate([np.linspace(0, 1, 801), np.linspace(1, 6, 1001)[1:]]); dd = integrate.cumulative_trapezoid([1 / hx(z) for z in zz], zz, initial=0)
    dpx = lambda zx: float(np.interp(zx, zz, dd))
    def at(z):
        zx = optimize.brentq(lambda q: math.log(1 + q) + ht * dpx(q) - math.log(1 + z), 0, z) if ht > 0 else z
        d = dpx(zx); return zx, d, (1 + z) * (hx(zx) / (1 + zx) + ht), d * (1 + zx) * math.sqrt((1 + z) / (1 + zx))
    return at
zb = np.array([0.012, 0.02, 0.03, 0.045, 0.07, 0.1, 0.15, 0.2, 0.26, 0.33, 0.42, 0.52, 0.65, 0.8, 1.0, 1.3, 1.7, 2.3]); sn_b = []
for a, b in zip(np.r_[0.01, (zb[:-1] + zb[1:]) / 2], np.r_[(zb[:-1] + zb[1:]) / 2, 2.4]):
    q = (zs >= a) & (zs < b)
    if q.sum(): w = 1 / P[q, 2] ** 2; sn_b.append((float(np.sum(w * zs[q]) / w.sum()), float(np.sum(w * P[q, 1]) / w.sum()), float(w.sum() ** -0.5)))
sn_b = np.array(sn_b)
def scan_point(ar, ht):
    at = family(ar, ht); c_ap = 0.0
    for d in DR2:
        if d[0] != 'MH': continue
        z = d[1]; F = d[2] / d[4]; eF = F * math.sqrt((d[3] / d[2]) ** 2 + (d[5] / d[4]) ** 2 - 2 * d[6] * d[3] * d[5] / (d[2] * d[4])); zx, dd_, he, dl = at(z); c_ap += ((F - dd_ * he) / eF) ** 2
    mu = np.array([5 * math.log10(at(z)[3]) for z in sn_b[:, 0]]); w = 1 / sn_b[:, 2] ** 2; off = np.sum(w * (sn_b[:, 1] - mu)) / w.sum()
    return c_ap, float(np.sum(w * (sn_b[:, 1] - mu - off) ** 2))
print("\n5. The family with both numbers free (shape chi2 for 6 ratios; supernova chi2 for %d bins):" % len(sn_b))
mu_s = np.array([5 * math.log10(lcdm(0.349)[0](z) * (1 + z)) for z in sn_b[:, 0]]); w = 1 / sn_b[:, 2] ** 2; off = np.sum(w * (sn_b[:, 1] - mu_s)) / w.sum(); sn_std = float(np.sum(w * (sn_b[:, 1] - mu_s - off) ** 2))
grid = []; 
for ar in (-0.95, -0.9, -0.85, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1):
    for ht in (0.0, 0.05, 0.1, 0.15, 0.22, 0.3, 0.4):
        c_ap, c_sn = scan_point(ar, ht); grid.append((ar, ht, c_ap, c_sn))
grid = np.array(grid); ib = grid[:, 2].argmin(); js = grid[:, 3].argmin(); tot = grid[:, 2] + grid[:, 3]; kt = tot.argmin()
print("   binned supernovae, standard model: %.1f" % sn_std)
print("   published point (alpha/H_c -0.80, H_t/H_c 0.22): shape %.1f, supernovae %.1f" % scan_point(AR, HT / HC))
print("   best shape:        alpha/H_c %.2f, H_t/H_c %.2f: shape %.1f, supernovae %.1f" % tuple(grid[ib]))
print("   best supernovae:   alpha/H_c %.2f, H_t/H_c %.2f: shape %.1f, supernovae %.1f" % tuple(grid[js]))
print("   best of the sum:   alpha/H_c %.2f, H_t/H_c %.2f: shape %.1f, supernovae %.1f" % tuple(grid[kt]))
print("   no tired light at all, best alpha: ", min(((c[2], c[3], c[0]) for c in grid if c[1] == 0.0), key=lambda t: t[0] + t[1]))
out["family"] = dict(sn_bins=int(len(sn_b)), sn_std=sn_std, published=scan_point(AR, HT / HC), best_shape=grid[ib].tolist(), best_sn=grid[js].tolist(), best_sum=grid[kt].tolist(), grid=grid.tolist())
json.dump(out, open(os.path.join(HERE, "ccc_tl_test.json"), "w"), indent=1, default=float)
fine = [(ar, ht) + scan_point(ar, ht) for ar in np.arange(-0.60, -0.349, 0.025) for ht in (0.0, 0.02, 0.05)]
bf = min(fine, key=lambda t: t[2] + t[3]); b0 = min((t for t in fine if t[1] == 0.05), key=lambda t: t[2] + t[3])
print("   finer grid: best sum at alpha/H_c %.3f, H_t/H_c %.2f: shape %.1f, supernovae %.1f, sum %.1f (standard: %.1f + %.1f = %.1f)" % (bf + (bf[2] + bf[3], c_shape_l, sn_std, c_shape_l + sn_std)))
print("   with H_t/H_c = 0.05 (a 5%% tired-light share): best sum %.1f at alpha/H_c %.3f" % (b0[2] + b0[3], b0[0]))
out["family"]["fine_best"] = list(bf); out["family"]["fine_ht005"] = list(b0); out["family"]["std_sum"] = c_shape_l + sn_std
json.dump(out, open(os.path.join(HERE, "ccc_tl_test.json"), "w"), indent=1, default=float)
