"""Workbench Entry 67: CCC+TL's nucleosynthesis claim, and its second branch on the shape test.
Gupta & Samaras (2026, arXiv:2604.03381) slow the expansion at nucleosynthesis by f_max ~ 3 and divide every reaction rate and the
neutron decay rate by the same factor, with the temperature-redshift relation unchanged.  Rescaling every rate and the expansion by one
factor is a change of time unit, so the yields cannot change; what is left to test is the baryon density, and the assumption that the
weak rates rescale.  Uses bbn_net.py (Entry 53) and ccc_tl_test.py (Entry 65) from this folder."""
import os, sys, io, json, math, contextlib, importlib.util
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bbn_net as B

OBS = dict(Yp=(0.2449, 0.0040), D=(2.527e-5, 0.030e-5))        # Aver et al. 2021; Cooke, Pettini & Steidel 2018
ETA_PLANCK = 6.10                                               # eta x 1e10 from the CMB baryon density
NUC = ['benp', 'ddn', 'ddp', 'dpg', 'hag', 'hdp', 'hent', 'lipa', 'npdg', 'tag', 'tdn']
out = {}
def show(name, r, note=""):
    zY = (r['Yp'] - OBS['Yp'][0]) / OBS['Yp'][1]
    zD = (r['D'] - OBS['D'][0]) / math.hypot(OBS['D'][1], 0.03 * r['D'])     # 3% allowance for the reaction rates
    print(f"  {name:58s} Yp {r['Yp']:.4f} ({zY:+5.1f} sd)   D/H {r['D']*1e5:5.2f}e-5 ({zD:+6.1f} sd)   7Li/H {r['Li7']*1e10:5.2f}e-10  {note}")
    out[name] = dict(Yp=r['Yp'], D=r['D'], Li7=r['Li7'], zY=zY, zD=zD)

print("Light elements under CCC+TL's assumptions (network of Entry 53; observed Yp 0.2449 +- 0.0040, D/H 2.527 +- 0.030e-5)")
print(" Every rate and the expansion slowed by the same factor is a change of time unit: the standard yields, by construction.")
show("1. standard, eta10 = 6.10 (also CCC+TL with every rate rescaled)", B.run(eta10=ETA_PLANCK))
eta_low = 0.56 * ETA_PLANCK
show(f"2. CCC+TL supernova branch: baryons 56% of standard, eta10 = {eta_low:.2f}", B.run(eta10=eta_low))
for fm in (2.0, 3.0):
    # In units where the expansion is standard, a neutron lifetime that is NOT slowed with the rest is a lifetime shorter by f_max;
    # in the network the same normalisation sets the n <-> p rates, so this is 'weak physics keeps its laboratory rates'.
    show(f"3. f_max = {fm:.0f}, weak rates not rescaled (lifetime 880 s in local time)", B.run(eta10=ETA_PLANCK, tau=B.TAUN / fm))
    show(f"4. f_max = {fm:.0f}, neither weak nor nuclear rates rescaled", B.run(eta10=ETA_PLANCK, tau=B.TAUN / fm, scale={k: fm for k in NUC}))
# which baryon density would the supernova branch need for deuterium? bisect eta10
lo, hi = 3.0, 9.0
for _ in range(22):
    m = 0.5 * (lo + hi); d = B.run(eta10=m)['D']
    lo, hi = (m, hi) if d > OBS['D'][0] else (lo, m)
print(f"  deuterium fixes eta10 = {0.5*(lo+hi):.2f} under the rescaling ansatz, {0.5*(lo+hi)/ETA_PLANCK:.2f} of the CMB value")
out['eta_from_D'] = 0.5 * (lo + hi)

# ---- the second branch on the shape test (Entry 65's code) ----
spec = importlib.util.spec_from_file_location("cct", os.path.join(HERE, "ccc_tl_test.py")); cct = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(cct)
print("\nThe clustering ruler's shape (6 DESI DR2 ratios) and binned supernovae (18 bins), as in Entry 65:")
for name, ar in (("supernova branch, alpha/H_c = -0.80", -0.7997), ("recombination branch, alpha/H_c = -0.75", -0.75)):
    ht = (1 + ar) * (3 + ar) / 2                         # H_t / H_c, as the model fixes it
    c_ap, c_sn = cct.scan_point(ar, ht)
    print(f"  {name}: H_t/H_c = {ht:.3f}; shape chi2 {c_ap:.1f}, supernova chi2 {c_sn:.1f}   (standard model: 6.0 and 13.0)")
    out[name] = dict(ar=ar, ht=ht, shape=c_ap, sn=c_sn)
json.dump(out, open(os.path.join(HERE, "ccctl_bbn.json"), "w"), indent=1)
