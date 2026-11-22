"""Workbench Entry 52, check: what happens to a random chain of about 44 when the linking length is changed?
The discovery paper reports that at 85 Mpc the Giant Arc is 'still mostly detected' in two parts, that shorter lengths give no structures,
and that 105 Mpc joins almost everything.  Here: random slices whose largest group at 95 Mpc has 40 to 48 members, relinked."""
import os, json
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
HERE=os.path.dirname(os.path.abspath(__file__)); rng=np.random.default_rng(524); BOX=np.array([1541.0,1615.0,338.0]); NP=504
def groups(pts,link):
    n=len(pts); pr=cKDTree(pts).query_pairs(link,output_type="ndarray")
    if len(pr)==0: return np.arange(n),np.ones(n,int)
    g=coo_matrix((np.ones(len(pr)),(pr[:,0],pr[:,1])),shape=(n,n)); k,lab=connected_components(g,directed=False); return lab,np.bincount(lab)
rows=[]; tried=0
while len(rows)<1500:
    tried+=1; pts=rng.random((NP,3))*BOX; lab,cnt=groups(pts,95.0)
    if not 40<=cnt.max()<=48: continue
    mem=np.where(lab==np.argmax(cnt))[0]; n=len(mem); r=[n]
    for link in (85.0,75.0,65.0):
        l2,c2=groups(pts[mem],link); s=np.sort(c2)[::-1]; r+=[s[0]/n,(s[0]+(s[1] if len(s)>1 else 0))/n,s[0]]
    l3,c3=groups(pts,105.0); r+=[c3[l3[mem[0]]],c3.max()]; l4,c4=groups(pts,85.0); r+=[c4.max()]
    rows.append(r)
R=np.array(rows,float); out={}
print(f"{len(R)} random slices whose largest group at 95 Mpc has 40 to 48 members ({100*len(R)/tried:.0f}% of slices)")
for i,link in enumerate((85,75,65)):
    a=R[:,1+3*i]; b=R[:,2+3*i]; c=R[:,3+3*i]
    print(f"  relinked at {link} Mpc: largest piece keeps a median {100*np.median(a):.0f}% of the members ({np.median(c):.0f} points; 16-84% {np.percentile(c,16):.0f}-{np.percentile(c,84):.0f}); the two largest pieces together {100*np.median(b):.0f}%; two pieces hold more than half in {100*np.mean(b>0.5):.0f}% of cases")
    out[f"l{link}"]=dict(top1=float(np.median(a)),top2=float(np.median(b)),n1=float(np.median(c)),n1_16=float(np.percentile(c,16)),n1_84=float(np.percentile(c,84)),half=float(np.mean(b>0.5)))
print(f"  at 105 Mpc the group containing it has a median {np.median(R[:,10]):.0f} members ({100*np.median(R[:,10])/NP:.0f}% of the slice); the largest group in the slice {np.median(R[:,11]):.0f}")
print(f"  at 85 Mpc the largest group anywhere in the slice: median {np.median(R[:,12]):.0f}")
out["l105"]=dict(n=float(np.median(R[:,10])),nmax=float(np.median(R[:,11])))
json.dump(out,open(os.path.join(HERE,"breakup.json"),"w"),indent=1)
