"""Workbench Entry 53, robustness of the tail scan: flat S-factors instead of the adopted shapes, and the lithium cut that deuterium allows."""
import os, json, sys
import numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import bbn_net as bn, tail_scan as ts
HERE=os.path.dirname(os.path.abspath(__file__))
S=np.array(ts.out['scan']); base=ts.out['base']; q=S[:,0]; d=S[:,2]/base['D']; li=S[:,4]/base['Li7']
print("\nfrom the scan:")
for target in (1.05,1.10):
    qq=np.interp(target,d,q); print(f"  deuterium up by {100*(target-1):.0f}%: q = {qq:.4f}, lithium x{np.interp(qq,q,li):.2f}")
qq=np.interp(-1/3,-li,q); print(f"  lithium down to one third: q = {qq:.4f}, deuterium x{np.interp(qq,q,d):.2f}")
i=list(q).index(1.01); print(f"  slope near the thermal case: d ln D / d ln Li = {np.log(d[i])/np.log(li[i]):.2f}")
saved={k:v for k,v in ts.PAIR.items()}
for k in ts.PAIR: ts.PAIR[k]=(saved[k][0],saved[k][1],lambda E: 1.0+0*E)
r=bn.run(tail=ts.ratio_table(1.0755)[0]); print(f"flat S-factors, q = 1.0755: D x{r['D']/base['D']:.3f}  3He x{r['He3']/base['He3']:.3f}  7Li x{r['Li7']/base['Li7']:.3f}   (adopted shapes: x1.439, x1.355, x0.296)")
json.dump(dict(flat=[float(r['D']/base['D']),float(r['He3']/base['He3']),float(r['Li7']/base['Li7'])]),open(os.path.join(HERE,'tail_robust.json'),'w'))
