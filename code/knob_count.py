"""Workbench Entry 48: counting the numbers. Published and our own fits that ease the mismatch between DESI's ruler and the CMB,
either by changing recombination (a shorter sound horizon) or by changing dark energy.
For each: the gain in chi-squared over the standard model, the number of extra numbers fitted, the chance of that gain if the extra
numbers did nothing, and the Hubble constant.  Then a first-order look at what a direct Hubble constant adds to each."""
import json, os
import numpy as np
from scipy import stats
HERE=os.path.dirname(os.path.abspath(__file__))
sig=lambda c,k: float(stats.norm.isf(stats.chi2.sf(c,k)/2))
# side, label, data, extra numbers, gain, H0, error, H0 of the standard model in the same data, its error
F=[("rec","one factor on the sound horizon (Entry 36)","DESI 2 + CMB summary",1,4.58,69.28,0.53,68.17,0.28),
   ("rec","clumping strength (Entry 36)","DESI 2 + CMB summary",1,4.59,69.2,0.53,68.17,0.28),
   ("rec","field strength, simulated recombination","Planck + DESI 1",1,4.7,68.52,0.58,67.88,0.37),
   ("rec","field strength, simulated recombination","Planck + ACT + DESI 1",1,5.3,None,None,None,None),
   ("rec","field strength, simulated recombination","Planck + SPT + DESI 1",1,5.2,None,None,None,None),
   ("rec","four-number recombination","Planck + DESI 2",4,9.0,69.29,0.63,68.19,0.28),
   ("rec","four-number recombination","Planck + DESI 2 + Pantheon+",4,6.6,69.10,0.62,68.11,0.28),
   ("rec","four-number recombination","Planck + ACT + DESI 2 + Pantheon+",4,7.41,69.13,0.44,68.35,0.26),
   ("rec","seven-number recombination","three CMB experiments + DESI 2",7,14.2,69.48,0.65,67.19,0.38),
   ("de","changing dark energy (Entry 36)","DESI 2 + CMB summary",2,8.1,63.6,2.1,68.17,0.28),
   ("de","changing dark energy","Planck + DESI 2",2,10.8,63.9,1.9,68.19,0.28),
   ("de","changing dark energy","Planck + DESI 2 + Pantheon+",2,8.5,67.51,0.59,68.11,0.28),
   ("de","changing dark energy","Planck + ACT + DESI 2 + Pantheon+",2,7.96,67.64,0.60,68.35,0.26),
   ("de","changing dark energy","three CMB experiments + DESI 2",2,13.5,None,None,None,None)]
DIRECT={"three direct methods, no Cepheids":(71.45,1.48),"all four direct methods":(72.52,0.85)}
out=[]
print(f"{'fit':44s} {'data':36s} k  gain  chance   sigma  gain-2k   H0")
for side,lab,data,k,g,h0,e,hs,es in F:
    r=dict(side=side,label=lab,data=data,k=k,gain=g,p=float(stats.chi2.sf(g,k)),sigma=sig(g,k),aic=g-2*k,H0=h0,eH0=e,direct={})
    line=f"{lab:44s} {data:36s} {k}  {g:5.2f} {r['p']:.4f}  {r['sigma']:.2f}   {g-2*k:+5.1f}   "+(f"{h0} +- {e}" if h0 else "-")
    if h0:
        for nm,(hd,ed) in DIRECT.items():
            cm=(hd-h0)**2/(ed**2+e**2); cs=(hd-hs)**2/(ed**2+es**2); tot=g+cs-cm
            r["direct"][nm]=dict(extra=float(cs-cm),total=float(tot),sigma=sig(max(tot,0.01),k)); line+=f" | {nm}: {cs-cm:+.1f} -> {tot:.1f} ({sig(max(tot,0.01),k):.1f} sigma)"
    print(line); out.append(r)
print("\npublished fits with the Cepheid ladder included as a prior (their own numbers): four-number recombination gains 24.4 (%.1f sigma for 4), changing dark energy 14 (%.1f sigma for 2); field strength gains 15.25 (%.1f sigma for 1)"%(sig(24.4,4),sig(14,2),sig(15.25,1)))
print("what extra freedom buys in recombination: from one number (4.6-5.3) to four (6.6-9.0) to seven (14.2): %.1f to %.1f per added number"%((6.6-5.3)/3,(14.2-4.6)/6))
print("thresholds of gain for 2 and 3 sigma, by number of extra numbers:",{k:(round(float(stats.chi2.isf(2*stats.norm.sf(2),k)),1),round(float(stats.chi2.isf(2*stats.norm.sf(3),k)),1)) for k in range(1,8)})
json.dump(dict(fits=out,thr={k:[float(stats.chi2.isf(2*stats.norm.sf(s),k)) for s in (1,2,3)] for k in range(1,8)}),open(os.path.join(HERE,"knob_count.json"),"w"),indent=1)
