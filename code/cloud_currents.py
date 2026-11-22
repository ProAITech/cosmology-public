"""Workbench Entry 73: if the line-of-sight field reversals across Orion A, California and Perseus were the toroidal field of an
axial current, how large would the current be, is it physically allowed, and what measurement would tell a current from a bent
ambient field?"""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
mu0 = 4e-7 * math.pi; pc = 3.0857e16; e = 1.602e-19
def current(B_uG, r_pc): return 2 * math.pi * r_pc * pc * B_uG * 1e-10 / mu0
out = {"cases": []}
print("Axial current that would make a toroidal field of strength B at distance r from the spine (I = 2 pi r B / mu0)")
for name, B, r in (("Perseus, |B_los| ~ 90-136 uG (Tahani+2022)", 110, 2.0), ("Perseus, same field at 5 pc", 110, 5.0),
                   ("Orion A, OMC-B region -129 uG (Tahani+2018)", 129, 1.0), ("Orion A, 100 uG at 3 pc", 100, 3.0)):
    I = current(B, r); print(f"   {name:46s} r = {r:.0f} pc: I = {I:.1e} A"); out["cases"].append(dict(name=name, B=B, r=r, I=I))
print("   for scale: the Milky Way halo circuit 9e16 A (Entry 18); a radio jet ~3e18 A")
# can weakly ionised cloud gas carry it?
B, r = 110e-10, 2.0 * pc
j = B / (mu0 * r)
for nH2, xe in ((1e2, 1e-6), (1e3, 1e-7), (1e4, 3e-8)):
    ne = 2 * nH2 * xe * 1e6
    print(f"   current density {j:.1e} A/m2; n_H2 = {nH2:.0e} cm-3, x_e = {xe:.0e}: electron drift {1e3*j/(ne*e):.1f} mm/s")
out["j"] = j
# what distinguishes: profile of B_los with distance x from the spine
print("\nWhat tells them apart: |B_los| against distance from the spine (toroid: falls as 1/x outside the current and is antisymmetric about the spine;")
print("bent ambient field: confined to the bend near the cloud edge, one side tied to the shell that bent it, returning to the ambient sign far out)")
# sightlines: Taylor et al. 2009 NVSS RMs ~1 per deg2; SPICE-RACS DR2 ~7 per deg2 (goodRM 6.7)
area = 8 * 6      # deg2: Orion A and its surroundings out to ~3 widths each side (~8 deg by 6 deg)
for name, dens in (("NVSS (Taylor+2009)", 1.0), ("SPICE-RACS DR2", 6.7)):
    print(f"   {name:20s} {dens:.1f} per deg2 -> about {area*dens:.0f} sightlines over {area} deg2 around Orion A")
out["sightlines"] = {"NVSS": area * 1.0, "SPICE_DR2": area * 6.7}
json.dump(out, open(os.path.join(HERE, "cloud_currents.json"), "w"), indent=1)
