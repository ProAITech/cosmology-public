"""Workbench Entry 52: how often do random points make a 'Giant Arc'?
The Giant Arc (Lopez et al. 2022) is a friends-of-friends group of 44 Mg II absorbers, 770 Mpc long, found with a linking length of 95 Mpc
among 504 absorbers in a slice 1541 x 1615 x 338 Mpc (comoving; numbers as given by Sawala et al. 2025).
Two papers disagree on whether the standard model makes such things, and both report that the simulation they use looks like random points
at this density.  So the question can be asked of random points directly: scatter 504 points in that slice, link them the same way,
and see how large the largest group is."""
import os, json
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
HERE=os.path.dirname(os.path.abspath(__file__))
BOX=np.array([1541.0,1615.0,338.0]); NP=504; rng=np.random.default_rng(52)
def largest(pts,link):
    pr=cKDTree(pts).query_pairs(link,output_type="ndarray")
    if len(pr)==0: return 1,0.0,0.0,2,0.0
    n=len(pts); g=coo_matrix((np.ones(len(pr)),(pr[:,0],pr[:,1])),shape=(n,n)); k,lab=connected_components(g,directed=False)
    cnt=np.bincount(lab); best=(0,0,0)
    # membership of the largest group, and the longest group with at least 20 members
    out_n=cnt.max(); ext_of_big=0.0; long_ext=0.0; long_n=0
    for c in np.where(cnt>=20)[0]:
        p=pts[lab==c]; d=np.sqrt(((p[:,None,:2]-p[None,:,:2])**2).sum(-1)).max()            # extent projected on the sky plane
        if cnt[c]==out_n: ext_of_big=max(ext_of_big,d)
        if d>long_ext: long_ext,long_n=d,cnt[c]
    return int(out_n),float(ext_of_big),float(long_ext),int(long_n),0.0
def run(link,nreal,box=BOX,npts=NP):
    R=np.array([largest(rng.random((npts,3))*box,link)[:4] for _ in range(nreal)]); return R
out={"box":BOX.tolist(),"n":NP}
dens=NP/np.prod(BOX); print(f"density {dens*1e9:.0f} per Gpc3; mean spacing {dens**(-1/3):.0f} Mpc; linking length 95 Mpc is {95*dens**(1/3):.2f} of it; n l^3 = {dens*95**3:.3f} (random points link up across all space near 0.65)")
NR=20000; R=run(95.0,NR); nmax,ext_big,long_ext,long_n=R.T
print(f"\n{NR} random slices, linking length 95 Mpc:")
print(f"  largest group: median {np.median(nmax):.0f} members, 16-84% {np.percentile(nmax,16):.0f}-{np.percentile(nmax,84):.0f}, 95th percentile {np.percentile(nmax,95):.0f}, 99th {np.percentile(nmax,99):.0f}")
p44=np.mean(nmax>=44); p770=np.mean(long_ext>=770); pboth=np.mean((nmax>=44)&(ext_big>=770)); p841=np.mean(long_ext>=841)
print(f"  a group of 44 or more: {100*p44:.1f}% of slices;  a group of 20 or more longer than 770 Mpc: {100*p770:.1f}%; longer than 841: {100*p841:.1f}%; 44 or more AND longer than 770: {100*pboth:.2f}%")
print(f"  longest group of 20 or more: median {np.median(long_ext[long_ext>0]):.0f} Mpc; slices with any such group {100*np.mean(long_ext>0):.0f}%")
out["l95"]=dict(nreal=NR,median=float(np.median(nmax)),p16=float(np.percentile(nmax,16)),p84=float(np.percentile(nmax,84)),p95=float(np.percentile(nmax,95)),p99=float(np.percentile(nmax,99)),p44=float(p44),p770=float(p770),p841=float(p841),pboth=float(pboth))
out["hist"]=[int(x) for x in np.bincount(nmax.astype(int),minlength=140)[:140]]
print("\nsensitivity to the linking length (5,000 slices each):"); out["link"]=[]
for link in (80,85,90,95,100,105,110):
    r=run(float(link),5000); a=np.mean(r[:,0]>=44); b=np.mean((r[:,0]>=44)&(r[:,1]>=770)); print(f"  {link:3d} Mpc: n l^3 = {dens*link**3:.2f}; median largest group {np.median(r[:,0]):3.0f}; 44 or more in {100*a:5.1f}% ; 44 or more and longer than 770 Mpc in {100*b:5.2f}%")
    out["link"].append([link,float(dens*link**3),float(np.median(r[:,0])),float(a),float(b)])
# the larger slabs of the simulation paper: 1588 points in 2800 x 2800 x 338
r=run(95.0,3000,box=np.array([2800.0,2800.0,338.0]),npts=1588)
print(f"\nslabs as in the simulation paper (1588 points in 2800 x 2800 x 338 Mpc), random points: a group of 20 or more longer than 770 Mpc in {100*np.mean(r[:,2]>=770):.0f}% (the paper: about 19%); longer than 841 in {100*np.mean(r[:,2]>=841):.0f}% (paper: about 11%); median longest {np.median(r[:,2][r[:,2]>0]):.0f} Mpc; 44 or more members in {100*np.mean(r[:,0]>=44):.0f}%")
out["bigslab"]=dict(p770=float(np.mean(r[:,2]>=770)),p841=float(np.mean(r[:,2]>=841)),p44=float(np.mean(r[:,0]>=44)))
# look-elsewhere: how many such slices does the absorber catalogue hold?
tot=63876; print(f"\nthe catalogue holds {tot} absorbers, {tot/NP:.0f} times the slice's 504")
for k in (10,30,100):
    print(f"  if {k} slices of this kind were open to inspection, the chance of at least one with 44 or more and longer than 770 Mpc: {100*(1-(1-pboth)**k):.0f}%; with 44 or more: {100*(1-(1-p44)**k):.0f}%")
out["look"]=[[k,float(1-(1-pboth)**k),float(1-(1-p44)**k)] for k in (10,30,100)]
json.dump(out,open(os.path.join(HERE,"random_arcs.json"),"w"),indent=1)
