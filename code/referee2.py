"""Workbench Entry 56: the Hubble constant as referee, updated.
Entry 37 compared what each reading of the DESI-CMB mismatch predicts for the Hubble constant with four direct measurements.
Since then: the lens time-delay result has a third less error (TDCOSMO 2026), a community distance-ladder value has been published
(73.50 +- 0.81), and a compilation lists every measurement that uses neither the ladder nor the sound horizon.
Here: the same comparison, with the measurements that need neither a ladder nor the sound horizon taken one per method."""
import os, json
import numpy as np
from scipy.stats import chi2
HERE=os.path.dirname(os.path.abspath(__file__)); out={}
read=[["constant dark energy, standard sound horizon",68.17,0.28],["constant dark energy, sound horizon 0.6% shorter (Entry 36)",69.28,0.53],
      ["changing dark energy, ruler + CMB",63.6,2.1],["changing dark energy, + Pantheon+",67.51,0.59],["changing dark energy, + Union3",65.91,0.84],["changing dark energy, + DES 5-year",66.74,0.56]]
ladder=["distance ladder, community value",73.50,0.81,0.81]
# one step from physics or geometry to a distance: no ladder, no sound horizon.  name, value, +err, -err, precise?
one=[["lens time delays, 11 lenses (TDCOSMO 2026)",71.1,2.4,2.3,1],["megamasers",73.5,3.0,2.9,1],["type II supernovae, expanding photosphere",74.9,2.7,2.7,1],
     ["supernova Refsdal, time delay behind a cluster",66.6,4.1,3.3,0],["kilonova of 2017, expanding photosphere",67.0,3.6,3.6,0],["the same event, gravitational waves and afterglow",75.46,5.34,5.39,0],
     ["cosmic chronometers",66.7,5.3,5.3,0],["supernova H0pe, time delay",71.8,9.2,8.1,0],["supernova Encore, time delay",66.9,11.2,8.1,0],["dark sirens",79.0,9.0,9.0,0]]
def mean(rows):
    v=np.array([r[1] for r in rows]); s=np.array([(r[2]+r[3])/2 for r in rows]); w=1/s**2; m=(w*v).sum()/w.sum(); return float(m),float(1/np.sqrt(w.sum())),float((w*(v-m)**2).sum())
def z(row,h,e):
    s=row[3] if h<row[1] else row[2]; return (row[1]-h)/np.hypot(s,e)
print("one-step measurements (no ladder, no sound horizon):")
for r in one: print(f"   {r[0]:50s} {r[1]:6.2f} +{r[2]:.1f} -{r[3]:.1f}")
m3=mean([r for r in one if r[4]]); m7=mean([r for r in one if not r[4]]); m10=mean(one)
print(f"   the three with errors of 3 or less: {m3[0]:.2f} +- {m3[1]:.2f} (chi2 {m3[2]:.1f} for 2)")
print(f"   the other seven: {m7[0]:.2f} +- {m7[1]:.2f} (chi2 {m7[2]:.1f} for 6); the two groups differ by {(m3[0]-m7[0])/np.hypot(m3[1],m7[1]):.1f} sigma")
print(f"   all ten: {m10[0]:.2f} +- {m10[1]:.2f} (chi2 {m10[2]:.1f} for 9, chance {chi2.sf(m10[2],9):.2f})")
print(f"   Entry 37 had lenses at 71.6 (+3.9, -3.3) and masers at 73.9 +- 3.0: the lens error has fallen by {100*(1-2.35/3.6):.0f}%")
out['one']=one; out['means']=dict(three=m3,seven=m7,ten=m10)
print("\neach reading against the measurements (sigma below; chi2 against the ten):")
print(f"   {'reading':58s} predicts        below all ten   below the three   below the ladder   chi2 (10)")
out['read']=[]
for nm,h,e in read+[["the distance ladder itself",73.50,0.81]]:
    a=(m10[0]-h)/np.hypot(m10[1],e); b=(m3[0]-h)/np.hypot(m3[1],e); c=(73.50-h)/np.hypot(0.81,e) if nm[:3]!="the" else 0.0; c2=sum(z(r,h,e)**2 for r in one)
    print(f"   {nm:58s} {h:5.2f} +- {e:.2f}   {a:+5.1f}          {b:+5.1f}            {c:+5.1f}            {c2:5.1f}")
    out['read'].append([nm,h,e,float(a),float(b),float(c),float(c2)])
# within the changing-dark-energy model the lens distances themselves are refitted
print(f"\ninside the changing-dark-energy model: lenses give 70.6 +- 2.2 (with the ruler and DES supernovae); the model's own CMB-anchored value is 66.74 +- 0.56: apart by {(70.6-66.74)/np.hypot(2.2,0.56):.1f} sigma")
m10b=mean([[one[0][0],70.6,2.2,2.2,1]]+one[1:]); print(f"   all ten with the lenses at 70.6: {m10b[0]:.2f} +- {m10b[1]:.2f}; above 66.74 by {(m10b[0]-66.74)/np.hypot(m10b[1],0.56):.1f} sigma")
out['w0wa']=dict(lens=70.6,sep=float((70.6-66.74)/np.hypot(2.2,0.56)),ten=m10b[0],ten_err=m10b[1],ten_sep=float((m10b[0]-66.74)/np.hypot(m10b[1],0.56)))
# the sound horizon each value implies, through the ruler: h x r_d = 101.5 Mpc (Entry 14)
print("\nsound horizon implied through the clustering ruler (h r_d = 101.5 Mpc); standard value 147.05 Mpc"); out['rd']=[]
print("   standard physics: 147.05 Mpc.  Ruler and CMB ask for -0.57 +- 0.27% (Entry 36).")
for nm,h,e in (("early-plasma routes, at their cap (Entry 14)",70.0,0.0),("all ten one-step measurements",m10[0],m10[1]),("the three precise ones",m3[0],m3[1]),("lenses alone",71.1,2.35),("the distance ladder",73.50,0.81)):
    rd=10150/h; er=rd*e/h; print(f"   {nm:46s} H0 {h:5.2f}: {rd:6.1f} +- {er:3.1f} Mpc, {100*(rd/147.05-1):+.1f} +- {100*er/147.05:.1f}%"); out['rd'].append([nm,float(h),float(e),float(rd),float(er),float(100*(rd/147.05-1)),float(100*er/147.05)])
print(f"   published by the lens team from its own combination: 140.7 +- 4.3 Mpc, {100*(140.7/147.05-1):+.1f} +- {100*4.3/147.05:.1f}%")
# a wider net: methods that lean on a model of something else (the light between galaxies, the gas between galaxies, a cluster's mass)
wide=one+[["gamma-ray absorption by background light",61.9,2.9,2.4,0],["fast radio bursts",68.81,5.64,5.64,0],["quasars lensed by clusters, analysis A",74.1,8.0,8.0,0],["quasars lensed by clusters, analysis B",59.1,3.6,3.5,0]]
mw=mean(wide); print(f"\nwith four model-dependent results added (gamma-ray absorption 61.9, fast radio bursts 68.8, cluster-lensed quasars 74.1 and 59.1): {mw[0]:.2f} +- {mw[1]:.2f} (chi2 {mw[2]:.1f} for 13, chance {chi2.sf(mw[2],13):.3f})")
print("   published compilation of 58 ladder-free, horizon-free results: 69.37 +- 0.34; 69.52 +- 0.42 allowing for shared inputs; its three groups 67.61 +- 0.96, 71.03 +- 0.69, 69.07 +- 0.44")
for nm,h,e in read: print(f"   {nm:58s}: {(mw[0]-h)/np.hypot(mw[1],e):+.1f} sigma below the wider mean; {(69.52-h)/np.hypot(0.42,e):+.1f} below the compilation")
out['wide']=dict(mean=mw,rows=[[nm,float((mw[0]-h)/np.hypot(mw[1],e)),float((69.52-h)/np.hypot(0.42,e))] for nm,h,e in read])
json.dump(out,open(os.path.join(HERE,'referee2.json'),'w'),indent=1)
