"""Workbench Entry 52, check: can plain random points reproduce the simulation paper's own numbers?
Sawala et al. (2025) cut slabs of 2800 x 2800 x 338 Mpc holding 1588 points, link at 95 Mpc, take the extent of a group as its largest
pairwise distance, keep groups with axis ratio b/a < 0.5, drop groups within one linking length of a face, and report that about 19% of
samples hold a group longer than 770 Mpc and about 11% one longer than 841 Mpc.  They state that random points behave like the simulation.
Here: the same cuts on random points, with the two readings of 'a face' (the four sides of the slab, or all six faces)."""
import os, json
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist
HERE=os.path.dirname(os.path.abspath(__file__))
BOX=np.array([2800.0,2800.0,338.0]); NP=1588; LINK=95.0; rng=np.random.default_rng(522); NR=3000
keys=("all","elong","elong_sides","elong_allfaces","sides"); best={k:np.zeros(NR) for k in keys}
for it in range(NR):
    pts=rng.random((NP,3))*BOX; pr=cKDTree(pts).query_pairs(LINK,output_type="ndarray")
    g=coo_matrix((np.ones(len(pr)),(pr[:,0],pr[:,1])),shape=(NP,NP)); k,lab=connected_components(g,directed=False); cnt=np.bincount(lab)
    for c in np.where(cnt>=8)[0]:                        # a chain needs at least 8 links of 95 Mpc to pass 700 Mpc
        p=pts[lab==c]; L=pdist(p).max()
        if L<700: continue
        ev=np.sort(np.linalg.eigvalsh(np.cov(p.T)))[::-1]; ba=np.sqrt(ev[1]/ev[0])
        side_ok=(p[:,:2].min()>LINK) and (p[:,:2].max()<BOX[0]-LINK); all_ok=side_ok and (p[:,2].min()>LINK) and (p[:,2].max()<BOX[2]-LINK)
        for key,ok in (("all",True),("elong",ba<0.5),("elong_sides",ba<0.5 and side_ok),("elong_allfaces",ba<0.5 and all_ok),("sides",side_ok)):
            if ok: best[key][it]=max(best[key][it],L)
out={}
print(f"{NR} random slabs of 1588 points in 2800 x 2800 x 338 Mpc, linked at 95 Mpc; share with a group longer than 770 / 841 Mpc (largest pairwise distance)")
lab={"all":"every group","elong":"axis ratio b/a < 0.5","elong_sides":"b/a < 0.5 and clear of the four sides","elong_allfaces":"b/a < 0.5 and clear of all six faces","sides":"clear of the four sides, any shape"}
for k in keys:
    a=np.mean(best[k]>=770); b=np.mean(best[k]>=841); out[k]=[float(a),float(b)]; print(f"  {lab[k]:42s}: {100*a:5.1f}%   {100*b:5.1f}%")
print("  simulation paper, for its mock samples    :  about 19%   about 11%")
json.dump(out,open(os.path.join(HERE,"slab_check.json"),"w"),indent=1)
