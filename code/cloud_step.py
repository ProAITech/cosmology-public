# Entry 161: the rotation-measure step a current along a molecular cloud would leave, from the NVSS grid
# (Taylor, Stil & Sunstrum 2009, VizieR J/ApJ/702/1230, 37,543 sources). Fits ran in the browser (cloud_step.js);
# this script records the results, computes significances and converts steps into currents.
import json, math
prof = {"ori": "-8 -17.9 6.2 24;-4.5 -28.3 6.2 22;-2.25 -29.1 9.1 10;-0.75 -5.1 12.8 5;0.75 15.6 14.4 4;2.25 10.3 7.6 15;4.5 25.2 7.6 15;8 30.7 6.2 22",
        "per": "-8 -37.4 8 25;-4.5 -34.1 11.9 11;-2.25 -42.2 10 16;-0.75 -4.4 11.8 11;0.75 37 12.1 11;2.25 23.7 13.4 9;4.5 42.3 8.7 21;8 26.5 7.3 30"}
nullz = {"ori": "-4.3 1.9 7.1 0.1 3.7 0.3 -0.6 0.7 -1.2 -0.4 -3.4 1.0 -2.3 -7.4 -0.1 -0.9 -1.7 1.3 1.5 -0.6 -1.2 -2.7 -2.1 -0.5 0.9 -1.6 0.7 0.8 0.2 1.8 3.6 -4.2 0.3 -1.3 7.0 -0.6 -0.1 1.9 -2.0 1.0 -1.4 0.7 1.4 -3.8 -0.4 -1.4 -2.9 -1.0 -0.5 -0.7 -0.4",
         "per": "-3.8 -10.0 -4.7 4.3 -0.6 -0.3 5.6 4.5 1.8 -1.6 1.1 2.0 -0.1 1.7 3.6 -1.0 0.8 -1.1 0.0 0.3 -2.4 -1.8 4.0 -3.0 -1.2 4.7 1.3 -3.4 -2.3 5.6 0.1 -4.1 -1.9 3.1 0.3 5.2 0.2 -2.1 1.1 1.3 0.4 3.2 -6.0 -1.3 -2.3 -1.4 -4.0 -3.4 6.4 0.7 -3.5 -0.5 3.3 5.7 -0.8 -3.1 -1.8 3.3 1.9 -2.4 1.4 1.5 -0.6 -1.1 0.6 0.0"}
fits = {  # weighted fits: RM = c0 + c1 d + c2 s + A_in sgn(d)[|d|<1.5] + A_out sgn(d)[|d|>=1.5], |d| <= 10 deg
 "ori": {"name": "Orion A", "spine_lb": [[208.6, -19.1], [214.5, -19.6]], "dist_pc": 420, "N": 117, "sig_int": 26.9, "Ain": [9.7, 9.7], "Aout": [23.8, 6.6], "near": [21.9, 6.4], "far": [30.8, 14.9]},
 "per": {"name": "Perseus", "spine_lb": [[158.35, -20.5], [160.57, -16.8]], "dist_pc": 300, "N": 134, "sig_int": 38.1, "Ain": [20.2, 8.6], "Aout": [34.7, 9.5], "near": [28.1, 7.8], "far": [34.0, 20.3]}}
K = 760.5  # rad/m^2 per side for I = 4.6e15 A and n_e = 1 cm^-3: 0.812 * mu0 I / 2 in uG pc
out = {"fits": fits, "profiles": {}, "null": {}, "current": {}}
zs_all = {}
for k in ("ori", "per"):
    out["profiles"][k] = [[float(x) for x in r.split()] for r in prof[k].split(";")]
    z = [float(x) for x in nullz[k].split()]; zs_all[k] = z
    zr = fits[k]["Aout"][0] / fits[k]["Aout"][1]; rms = math.sqrt(sum(v * v for v in z) / len(z))
    p = sum(abs(v) >= abs(zr) for v in z) / len(z)
    eA = fits[k]["Aout"][1] * rms  # step error calibrated on the sky
    I = {ne: fits[k]["Aout"][0] / (K * ne) * 4.6e15 for ne in (0.01, 0.03, 0.05)}
    I95 = (fits[k]["Aout"][0] + 1.96 * eA) / (K * 0.03) * 4.6e15
    out["null"][k] = {"n": len(z), "z_real": round(zr, 2), "rms_null_z": round(rms, 2), "p": round(p, 3), "A_err_sky": round(eA, 1), "z_sky": round(fits[k]["Aout"][0] / eA, 2)}
    out["current"][k] = {"I_ne0.01": I[0.01], "I_ne0.03": I[0.03], "I_ne0.05": I[0.05], "I95_ne0.03": I95,
                         "pc_per_deg": round(fits[k]["dist_pc"] * math.pi / 180, 1)}
zr2 = sum((fits[k]["Aout"][0] / fits[k]["Aout"][1]) ** 2 for k in ("ori", "per"))
cnt = sum(a * a + b * b >= zr2 for a in zs_all["ori"] for b in zs_all["per"]); tot = len(zs_all["ori"]) * len(zs_all["per"])
pf = [out["null"][k]["p"] for k in ("ori", "per")]; X = -2 * sum(math.log(v) for v in pf); fisher = math.exp(-X / 2) * (1 + X / 2)
out["combined"] = {"p_sumz2": round(cnt / tot, 3), "p_fisher": round(fisher, 3)}
json.dump(out, open("cloud_step.json", "w"), indent=1)
print(json.dumps({"null": out["null"], "combined": out["combined"], "current": {k: {a: (f"{b:.2e}" if isinstance(b, float) and b > 1e10 else b) for a, b in v.items()} for k, v in out["current"].items()}}))
