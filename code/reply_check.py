"""Workbench Entry 52, check: the reply's 'nowhere to be seen', on random points.
Lopez & Clowes (2025) searched the whole 2800 Mpc cube of the simulation, thinned to the density of the Giant Arc field, with a linking
length of 65 Mpc (the mean nearest-neighbour distance of random points at that density, 0.554 n^(-1/3)), and kept groups of 10 or more.
They found no gigaparsec structures: the largest group in any of 100 samples had 19 members and few passed 370 Mpc.
The Giant Arc itself was found with a linking length of 95 Mpc.  Here: random points at the same density, at both lengths."""
import os, json
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist
HERE=os.path.dirname(os.path.abspath(__file__)); rng=np.random.default_rng(523); nbar=504/(1541.0*1615.0*338.0); out={}
print(f"density {nbar*1e9:.0f} per Gpc3; mean nearest-neighbour distance of random points {0.554*nbar**(-1/3):.0f} Mpc")
def groups(pts,link):
    n=len(pts); pr=cKDTree(pts).query_pairs(link,output_type="ndarray"); g=coo_matrix((np.ones(len(pr)),(pr[:,0],pr[:,1])),shape=(n,n)); k,lab=connected_components(g,directed=False); return lab,np.bincount(lab)
# (a) the whole cube, as in the reply
side=2800.0; npts=int(round(nbar*side**3)); print(f"\ncube of {side:.0f} Mpc: {npts} random points, 100 samples")
for link in (65.0,95.0):
    big=[]; n370=[]; n1000=[]; longest=[]
    for it in range(100 if link==65.0 else 20):
        pts=rng.random((npts,3))*side; lab,cnt=groups(pts,link); big.append(cnt.max()); c370=0; c1000=0; lmax=0
        for c in np.where(cnt>=10)[0]:
            p=pts[lab==c]
            if len(p)>3000: L=np.ptp(p,axis=0).max()          # one web: the span along an axis is enough
            else: L=pdist(p).max()
            c370+=L>370; c1000+=L>1000; lmax=max(lmax,L)
        n370.append(c370); n1000.append(c1000); longest.append(lmax)
    print(f"  linking {link:.0f} Mpc (n l^3 = {nbar*link**3:.2f}): largest group median {np.median(big):.0f} members, largest in any sample {np.max(big)}; groups of 10 or more longer than 370 Mpc: {np.mean(n370):.1f} per sample; longer than 1000 Mpc: {np.mean(n1000):.2f} per sample; longest {np.max(longest):.0f} Mpc")
    out[f"cube{int(link)}"]=dict(med=float(np.median(big)),max=int(np.max(big)),n370=float(np.mean(n370)),n1000=float(np.mean(n1000)),longest=float(np.max(longest)))
print("  the reply, for the simulation at 65 Mpc: largest group 19 members; only a few longer than 370 Mpc; none near 1000 Mpc")
# (b) the Arc's own slice at 65 Mpc
BOX=np.array([1541.0,1615.0,338.0]); big=[]
for it in range(5000):
    lab,cnt=groups(rng.random((504,3))*BOX,65.0); big.append(cnt.max())
big=np.array(big); print(f"\nthe Arc's slice at 65 Mpc: largest group median {np.median(big):.0f}, 99th percentile {np.percentile(big,99):.0f}, largest in 5000 slices {big.max()}; 44 or more in {100*np.mean(big>=44):.2f}%")
out["slice65"]=dict(med=float(np.median(big)),p99=float(np.percentile(big,99)),max=int(big.max()))
json.dump(out,open(os.path.join(HERE,"reply_check.json"),"w"),indent=1)
