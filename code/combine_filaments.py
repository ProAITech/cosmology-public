"""Workbench Entry 68: the node-oriented filament current from two independent rotation-measure grids, in units of 1e17 A
(n_e = 1e-5 cm^-3 within 1 Mpc of the spine).  LoTSS numbers from Entries 13 and 15; SPICE-RACS numbers from filament_spice.js
(results transcribed below, as run on 11 November 2026)."""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
# LoTSS (Entry 13): B = 0.044 +- 0.062 rad/m2; 1e17 A shifts B by 2.5 sigma when injected after the published foreground subtraction.
resp_L = 2.5 * 0.062
L_full = (0.044 / resp_L, 0.062 / resp_L)                  # no loss of signal in the published foreground
L_half = (0.044 / (resp_L / 2), 0.062 / (resp_L / 2))      # half the signal lost (Entry 15's allowance)
# SPICE-RACS DR1: B with the foreground rebuilt from the same sources, response measured through that rebuild
spice = dict(Rh=0.7, Rfg=1.0, sources=3915, withResid=3685, crossings=1768, crossSources=743, filaments=406, medZ=0.0335,
             B=-0.2101, B_sd_flips=0.187, B_sd_ctrl=0.163, respB=0.1445,
             pairs15=dict(pairs=2572, P=-1.025, P_sd_flips=0.5015, respP=0.1743), pairs10=dict(pairs=1404, P=-0.357, P_sd_flips=0.289, P_sd_ctrl=0.408, respP=0.172),
             Q=dict(Q=-7.555, mean0=0.111, sd=5.08, p=0.962, lim=9.0),
             grid=[[0.7,0.7,1.0,1768,406,734,-0.21,0.193,0.145,2.67,-0.0754,0.299,0.173,3.46],[0.7,1.5,0.7,1768,406,2572,-0.254,0.243,0.132,3.67,-1.03,0.5,0.174,5.73],
                   [0.7,1.5,1.5,1768,406,2572,-0.0467,0.216,0.153,2.82,-1.03,0.5,0.174,5.73],[1.0,1.5,1.0,2524,465,3963,-0.166,0.161,0.143,2.26,-0.656,0.418,0.178,4.7],[0.5,1.5,1.0,1252,348,1555,-0.19,0.224,0.146,3.07,-1.41,0.645,0.171,7.53]],
             controls=[[-50,1269,962,0.171,0.997],[-38,1520,1118,0.219,-0.0141],[-26,1761,1317,0.11,0.349],[-14,1930,1583,0.158,0.0128],[14,1518,1040,0.171,0.526],[26,1314,896,-0.0911,0.13],[38,783,501,0.0981,0.475],[50,382,192,0.535,-0.457]],
             thirds=dict(z_edges=[0.0237,0.0535],near=[589,-0.283,0.234],mid=[589,0.115,0.375],far=[590,-0.457,0.358]))
sd_S = max(spice['B_sd_flips'], spice['B_sd_ctrl'])
S = (spice['B'] / spice['respB'], sd_S / spice['respB'])
P = (spice['pairs10']['P'] / spice['pairs10']['respP'], spice['pairs10']['P_sd_ctrl'] / spice['pairs10']['respP'])
def comb(*xs):
    w = [1 / s ** 2 for _, s in xs]; m = sum(wi * x for wi, (x, _) in zip(w, xs)) / sum(w); return m, 1 / math.sqrt(sum(w))
rows = {"LoTSS, Entry 13 (no signal lost)": L_full, "LoTSS, half the signal lost": L_half, "SPICE-RACS DR1, single sources": S, "SPICE-RACS DR1, pairs within 1°": P,
        "both grids (LoTSS half lost)": comb(L_half, S), "both grids (LoTSS none lost)": comb(L_full, S)}
print("Current toward the denser end, 1e17 A (n_e = 1e-5 cm^-3); 95% limit = max(estimate, 0) + 2 sd")
out = {}
for k, (m, s) in rows.items():
    lim = max(m, 0) + 2 * s; print(f"  {k:38s} {m:+5.2f} +- {s:4.2f}   limit {lim:4.1f}   (at n_e = 3e-6: {lim/0.3:4.1f})"); out[k] = dict(est=m, sd=s, lim=lim, lim_ne3e6=lim / 0.3)
print("  a field of 11-15 nG at 1 Mpc, if wrapped by a current: I = 2 pi r B / mu0 = %.1f to %.1f (1e17 A), whatever the gas density" % (2*math.pi*3.0857e22*11e-13/(4e-7*math.pi)/1e17, 2*math.pi*3.0857e22*15e-13/(4e-7*math.pi)/1e17))
json.dump(dict(spice=spice, rows=out), open(os.path.join(HERE, "filament_spice.json"), "w"), indent=1)
