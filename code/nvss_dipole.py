# Entry 163: the NVSS radio count dipole measured here. Fits ran in the browser (nvss_dipole.js); this records them and
# computes the kinematic expectation, ratios and excess significance with errors scaled for clustering (overdispersion 1.297).
import json, math
beta = 369.82 / 299792.458; alpha = 0.75; F = 1.297; sF = math.sqrt(F)
# flux series, |b| > 15, Dec > -30: S, N, amp%, eamp%, RA, Dec, angle to CMB, CMB-direction component %, its error %, count slope x
flux = [l.split() for l in """15 274216 1.86 0.325 132 -9.3 35.7 1.511 0.292 0.877
20 213666 2.056 0.378 135.1 -18.1 33.9 1.707 0.33 0.928
25 174264 2.172 0.426 138.9 -23 32.2 1.838 0.366 0.965
30 146741 2.224 0.459 138.6 -21.3 31.7 1.891 0.398 0.996
40 110483 2.717 0.616 132.2 -43.3 48 1.819 0.459 1.057
60 71936 3.196 0.821 158.4 -59.1 52.7 1.936 0.568 1.15""".splitlines()]
# mask variants at S >= 20 mJy: Dec min, |b| cut, N, amp, eamp, RA, Dec, angle, CMB component, error
masks = [l.split() for l in """-40 10 259405 1.88 0.334 140.1 -19.8 29.9 1.631 0.306
-40 15 233875 2.331 0.361 138.6 -23.1 32.5 1.967 0.316
-40 20 208108 2.173 0.381 144.8 -22.6 27.3 1.932 0.331
-30 10 236884 1.594 0.346 136.4 -13.4 31.7 1.356 0.32
-30 15 213666 2.056 0.378 135.1 -18.1 33.9 1.707 0.33
-30 20 190589 1.907 0.394 144.4 -18 25.5 1.721 0.345
-20 10 210420 1.694 0.365 137.5 -8.9 30.2 1.463 0.343
-20 15 189491 2.147 0.394 134.8 -10.7 33 1.801 0.355
-20 20 168816 1.969 0.403 142.3 -3.9 25.7 1.775 0.371
-10 10 180550 1.651 0.454 131.1 -21.2 38.3 1.295 0.381
-10 15 161898 2.16 0.48 129.1 -19.6 39.7 1.661 0.396
-10 20 143687 1.777 0.472 140.3 -12.3 27.8 1.573 0.415""".splitlines()]
quad = {"-40": [2.472, 0.406, 138.4, -25.8, 33.9], "-30": [2.126, 0.425, 134.1, -14.7, 34.1], "-10": [3.2, 1.039, 118.3, 21.9, 56.5]}
dens = [[-30, 9.6], [-20, 9.47], [-10, 9.48], [0, 9.43], [10, 9.56], [20, 9.53], [30, 9.5], [40, 9.51], [50, 9.55], [60, 9.43], [70, 9.38], [80, 9.47]]
rows = []
for r in flux:
    S, N = float(r[0]), int(r[1]); amp, eamp, ra, de, ang, pc, epc, x = map(float, r[2:])
    kin = (2 + x * (1 + alpha)) * beta * 100
    rows.append({"S": S, "N": N, "amp": amp, "eamp": round(eamp * sF, 3), "ra": ra, "dec": de, "ang": ang, "pc": pc, "epc": round(epc * sF, 3), "x": x,
                 "kin": round(kin, 3), "ratio": round(pc / kin, 2), "eratio": round(epc * sF / kin, 2), "excess_sig": round((pc - kin) / (epc * sF), 2), "added": round(pc - kin, 2)})
base = rows[1]
var = [{"dmin": float(m[0]), "bcut": float(m[1]), "N": int(m[2]), "amp": float(m[3]), "ra": float(m[5]), "dec": float(m[6]), "ang": float(m[7]), "pc": float(m[8])} for m in masks]
out = {"beta": beta, "alpha": alpha, "overdispersion": F, "flux": rows, "masks": var, "quad": quad, "density_by_dec": dens, "baseline": base,
       "mask_range": {"amp": [min(v["amp"] for v in var), max(v["amp"] for v in var)], "pc": [min(v["pc"] for v in var), max(v["pc"] for v in var)], "ang": [min(v["ang"] for v in var), max(v["ang"] for v in var)]}}
json.dump(out, open("nvss_dipole.json", "w"), indent=1)
print(json.dumps(base), json.dumps(out["mask_range"]))
for r in rows: print(r["S"], r["kin"], r["ratio"], r["eratio"], r["excess_sig"], r["added"])
