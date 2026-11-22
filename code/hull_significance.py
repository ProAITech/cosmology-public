"""Workbench Entry 52, second part: the published 'convex hull' significance, applied to random points.
The method is that of Clowes et al. (2012), used for the Giant Arc by Lopez et al. (2022):
  - pad each member of a group to a sphere of radius half the mean link of the group's minimal spanning tree;
  - take the volume of the convex hull of those spheres (CHMS volume);
  - compare it with the CHMS volumes of sets of N points scattered in a cube at the mean density.
    Significance = (mean random volume - group volume) / scatter of random volumes; overdensity = mean random volume / group volume - 1.
The random sets in that comparison are simply scattered.  The group was found by linking close points.  Here the group is ALSO taken
from random points, found the same way, so that like is compared with like."""
import os, json
import numpy as np
from scipy.spatial import cKDTree, ConvexHull
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform
HERE=os.path.dirname(os.path.abspath(__file__))
BOX=np.array([1541.0,1615.0,338.0]); NP=504; LINK=95.0; rng=np.random.default_rng(521); nbar=NP/np.prod(BOX)
def chms(p):
    """Mean MST link and exact volume of the convex hull of equal spheres (Steiner: V + S r + M r^2 + 4 pi r^3 / 3)."""
    d=squareform(pdist(p)); e=minimum_spanning_tree(d).data; r=e.mean()/2
    h=ConvexHull(p); nrm=h.equations[:,:3]; M=0.0
    for k in range(3):
        j=h.neighbors[:,k]; oth=[m for m in range(3) if m!=k]; a=p[h.simplices[:,oth[0]]]; b=p[h.simplices[:,oth[1]]]
        M+=(np.linalg.norm(a-b,axis=1)*np.arccos(np.clip((nrm*nrm[j]).sum(1),-1,1))).sum()
    M/=4.0                                                   # each edge met twice; mean-curvature integral is half of sum(L*theta)
    return float(e.mean()), float(h.volume+h.area*r+M*r*r+4*np.pi*r**3/3)
# check of the Steiner sum on a unit cube padded by 0.5: 1 + 6*0.5 + 3*pi*0.25 + 4*pi*0.125/3
cube=np.array([[x,y,z] for x in (0,1) for y in (0,1) for z in (0,1)],float); h=ConvexHull(cube); nrm=h.equations[:,:3]; M=0
for k in range(3):
    j=h.neighbors[:,k]; oth=[m for m in range(3) if m!=k]; M+=(np.linalg.norm(cube[h.simplices[:,oth[0]]]-cube[h.simplices[:,oth[1]]],axis=1)*np.arccos(np.clip((nrm*nrm[j]).sum(1),-1,1))).sum()
assert abs(M/4-3*np.pi)<1e-9
# reference: CHMS volumes of N scattered points in a cube at the mean density
NG=np.array([5,6,7,8,10,12,15,20,25,30,36,44,52,60,72,85,100,120,150,190]); NSET=1000; ref=[]
for n in NG:
    side=(n/nbar)**(1/3); v=np.array([chms(rng.random((n,3))*side)[1] for _ in range(NSET)]); ref.append([n,v.mean()/(n/nbar),v.std()/v.mean()])
ref=np.array(ref); mu_of=lambda n: np.interp(np.log(n),np.log(ref[:,0]),ref[:,1])*n/nbar; sd_of=lambda n: np.interp(np.log(n),np.log(ref[:,0]),ref[:,2])*mu_of(n)
print("scattered sets of N points in a cube at the mean density: mean CHMS volume / cube volume, and relative scatter")
for n,a,b in ref: print(f"  N={n:4.0f}: {a:.3f}  {b:.3f}")
i44=list(NG).index(44); GA_N=44; GA_D=0.9
ga_sig_corr=(1-1/(1+GA_D))/ref[i44,2]                         # overdensity defined against the mean random volume
ga_sig_unc=(ref[i44,1]-1/(1+GA_D))/(ref[i44,1]*ref[i44,2])    # overdensity defined against the mean density
print(f"\nGiant Arc as quoted (44 members, overdensity 0.9): {ga_sig_corr:.1f} sigma by this statistic if the overdensity is against the mean random volume, {ga_sig_unc:.1f} if against the mean density (published: 4.5)")
# the same statistic for groups found by linking in random slices
NR=4000; big=[]; allg=[]; anyhit=[]
for it in range(NR):
    pts=rng.random((NP,3))*BOX; pr=cKDTree(pts).query_pairs(LINK,output_type="ndarray")
    g=coo_matrix((np.ones(len(pr)),(pr[:,0],pr[:,1])),shape=(NP,NP)); k,lab=connected_components(g,directed=False); cnt=np.bincount(lab); cmax=np.argmax(cnt); best=-9
    for c in np.where(cnt>=10)[0]:
        p=pts[lab==c]; n=len(p); edge,V=chms(p); s=(mu_of(n)-V)/sd_of(n); dl=mu_of(n)/V-1; L3=pdist(p).max(); L2=pdist(p[:,:2]).max(); best=max(best,s)
        row=[n,edge,V,dl,s,L3,L2,n/(nbar*V)-1]; allg.append(row)
        if c==cmax: big.append(row)
    anyhit.append(best)
B=np.array(big); A=np.array(allg); anyhit=np.array(anyhit); N,edge,V,dl,S,L3,L2,dunc=B.T
print(f"\nlargest linked group in each of {NR} random slices:")
print(f"  members median {np.median(N):.0f}; mean link {np.median(edge):.0f} Mpc; CHMS volume median {np.median(V)/1e9:.3f} Gpc3")
print(f"  overdensity: median {np.median(dl):.2f}, 16-84% {np.percentile(dl,16):.2f}-{np.percentile(dl,84):.2f}; the Arc's 0.9 is at the {100*np.mean(dl<GA_D):.0f}th percentile")
print(f"  significance: median {np.median(S):.1f} sigma, 16-84% {np.percentile(S,16):.1f}-{np.percentile(S,84):.1f}; at least 4.5 in {100*np.mean(S>=4.5):.0f}% of slices; at least 3 in {100*np.mean(S>=3):.0f}%")
q=N>=44; print(f"  slices whose largest group has 44 or more ({100*q.mean():.1f}%): overdensity median {np.median(dl[q]):.2f}; significance median {np.median(S[q]):.1f}; at least 4.5 sigma in {100*np.mean(S[q]>=4.5):.0f}% of them")
print(f"  44 or more AND at least 4.5 sigma: {100*np.mean(q&(S>=4.5)):.1f}% of slices; AND projected length over 770 Mpc as well: {100*np.mean(q&(S>=4.5)&(L2>=770)):.2f}%")
print(f"  any group of 10 or more in the slice reaching 4.5 sigma: {100*np.mean(anyhit>=4.5):.0f}% of slices")
lg=A[A[:,5]>700]; print(f"\nall groups longer than 700 Mpc ({len(lg)} in {NR} slices): overdensity median {np.median(lg[:,3]):.2f}; the Arc's 0.9 is at the {100*np.mean(lg[:,3]<GA_D):.0f}th percentile (simulation paper: 89th for random points); significance median {np.median(lg[:,4]):.1f}, at least 4.5 in {100*np.mean(lg[:,4]>=4.5):.0f}%")
m44=A[(A[:,0]>=40)&(A[:,0]<=48)]; print(f"groups of 40 to 48 members ({len(m44)}): overdensity median {np.median(m44[:,3]):.2f}, significance median {np.median(m44[:,4]):.1f}, at least 4.5 in {100*np.mean(m44[:,4]>=4.5):.0f}%; the Arc's 0.9 at the {100*np.mean(m44[:,3]<GA_D):.0f}th percentile")
from scipy.stats import norm
print(f"\nas a one-sided chance per slice: 44 or more and 4.5 sigma, {100*np.mean(q&(S>=4.5)):.2f}% = {norm.isf(np.mean(q&(S>=4.5))):.1f} sigma; with length over 770 Mpc too, {100*np.mean(q&(S>=4.5)&(L2>=770)):.2f}% ({int(np.sum(q&(S>=4.5)&(L2>=770)))} slices) = {norm.isf(np.mean(q&(S>=4.5)&(L2>=770))):.1f} sigma")
print(f"largest group at least 4.5 sigma: {100*np.mean(S>=4.5):.1f}%; outside -1..9: {int(np.sum((S<-1)|(S>=9)))}; overdensity against the mean density, long groups: Arc's 0.9 at the {100*np.mean(lg[:,7]<GA_D):.0f}th percentile; largest groups: median {np.median(dunc):.2f}")
for kk in (10,30,127):
    print(f"  {kk} slices open to inspection: at least one with all three in {100*(1-(1-np.mean(q&(S>=4.5)&(L2>=770)))**kk):.0f}%; with 44 or more and 4.5 sigma in {100*(1-(1-np.mean(q&(S>=4.5)))**kk):.0f}%")
out=dict(ref=ref.tolist(),ga_sig_corr=float(ga_sig_corr),ga_sig_unc=float(ga_sig_unc),med_delta=float(np.median(dl)),pct_delta=float(np.mean(dl<GA_D)),med_S=float(np.median(S)),p16=float(np.percentile(S,16)),p84=float(np.percentile(S,84)),
         frac45=float(np.mean(S>=4.5)),frac44=float(q.mean()),cond_med=float(np.median(S[q])),cond45=float(np.mean(S[q]>=4.5)),both=float(np.mean(q&(S>=4.5))),all3=float(np.mean(q&(S>=4.5)&(L2>=770))),anyhit=float(np.mean(anyhit>=4.5)),
         long_pct=float(np.mean(lg[:,3]<GA_D)),long_med=float(np.median(lg[:,3])),m44_med_S=float(np.median(m44[:,4])),m44_45=float(np.mean(m44[:,4]>=4.5)),m44_pct=float(np.mean(m44[:,3]<GA_D)),
         hist=np.histogram(S,bins=np.arange(-1,9.01,0.5))[0].tolist(),hist_lo=-1,hist_step=0.5,nsig=NR,outside=int(np.sum((S<-1)|(S>=9))),
         frac45_exact=float(np.mean(S>=4.5)),frac3=float(np.mean(S>=3)),long_pct_unc=float(np.mean(lg[:,7]<GA_D)),p_all3_n=int(np.sum(q&(S>=4.5)&(L2>=770))))
json.dump(out,open(os.path.join(HERE,"hull_significance.json"),"w"),indent=1)
