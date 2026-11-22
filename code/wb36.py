"""Workbench Entry 114: why 36 wide binaries with 3D velocities can give a 60% gravity boost or none.

Two 2026 analyses of the same 36 nearby wide binaries (radial velocities to < 100 m/s, hidden companions screened) disagree:
Chae et al. find G_eff/G = gamma = 1.60 (+0.17, -0.14); Saad & Ting, with the semi-major axis free, find 1.12 (0.90-1.38), and recover 1.56 when
the semi-major axis is replaced by a geometric de-projection of the projected separation.
This script asks how much a choice of that kind can move gamma when gravity is exactly Newtonian.

Method. A library of Keplerian orbits with a = GM = 1: eccentricity from f(e) ~ e^alpha, uniform mean anomaly, isotropic orientation.
For each: projected separation x = s/a and 3D relative speed u = v/sqrt(GM/a). The observable for a pair is v~ = v / sqrt(GM/s) = u sqrt(x).
If semi-major axes follow dN/da ~ a^-k, the distribution of v~ at fixed s is the library's, weighted by x^(k-1), and does not depend on s
(scale-free). So each pair's likelihood is f(v~/sqrt(gamma))/sqrt(gamma), with f fixed by (k, alpha). Mock samples of 36 Newtonian pairs
(gamma = 1) are drawn with a true (k, alpha), velocity errors of 7% and mass errors of 5%, and fitted with other (k, alpha), and with a
geometric short-cut: each pair placed at r = kappa s on a circular orbit, so gamma_i = v~^2 kappa."""
import os, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); rng=np.random.default_rng(114)
NL=1_500_000
def library(alpha):
    e=rng.random(NL)**(1/(alpha+1)); M=2*np.pi*rng.random(NL); E=M.copy()
    for _ in range(60): E=E-(E-e*np.sin(E)-M)/(1-e*np.cos(E))
    xo=np.cos(E)-e; yo=np.sqrt(1-e*e)*np.sin(E); r=1-e*np.cos(E)
    cz=2*rng.random(NL)-1; ph=2*np.pi*rng.random(NL); sz=np.sqrt(1-cz*cz)
    los=np.stack([sz*np.cos(ph),sz*np.sin(ph),cz])        # line of sight in the orbital frame (orbit in the x-y plane)
    dz=xo*los[0]+yo*los[1]; s=np.sqrt(np.maximum(r*r-dz*dz,1e-12))
    u=np.sqrt(2/r-1)
    return dict(x=s,u=u,r=r)
LIB={a:library(a) for a in (0.0,1.0,1.3,2.0)}
GRID=np.linspace(0,3.2,321)
def f_of(alpha,k):
    L=LIB[alpha]; vt=L['u']*np.sqrt(L['x']); w=L['x']**(k-1)
    h,_=np.histogram(vt,bins=GRID,weights=w,density=True); return h
def draw(alpha,k,n):
    L=LIB[alpha]; w=L['x']**(k-1); w=w/w.sum(); i=rng.choice(NL,size=n,p=w)
    return L['u'][i]*np.sqrt(L['x'][i]), L['r'][i]/L['x'][i]
mid=0.5*(GRID[1:]+GRID[:-1]); dg=GRID[1]-GRID[0]
def fit(vt,h,sig=0.07):
    # likelihood with Gaussian velocity error sig (fractional) and 5% mass error (2.5% on v~), on a gamma grid
    gs=np.linspace(0.5,3.0,251); s=np.hypot(sig,0.025)
    ll=[]
    for g in gs:
        # model density of observed v~: f(t/sqrt g)/sqrt g convolved with N(t, s t)
        t=mid*np.sqrt(g); dens=h/np.sqrt(g)
        P=np.array([np.sum(dens*dg*np.sqrt(g)*np.exp(-0.5*((v-t)/(s*t+1e-9))**2)/(np.sqrt(2*np.pi)*(s*t+1e-9))) for v in vt])
        ll.append(np.sum(np.log(np.maximum(P,1e-300))))
    ll=np.array(ll); i=int(np.argmax(ll)); return gs[i], gs, ll
out={"note":"gamma = G_eff/G fitted to mock Newtonian samples of 36 pairs"}
# median r/s for the geometric short-cut
for a in (1.0,1.3):
    L=LIB[a]; rs=L['r']/L['x']; out[f"r_over_s_alpha{a}"]=dict(median=float(np.median(rs)),mean_isotropic=4/np.pi)
TRUE=(1.3,1.6)                                                # superthermal wide-binary eccentricities; dN/da ~ a^-1.6
NM=300; out["cases"]={}
cases=[("matched: alpha 1.3, k 1.6",1.3,1.6),("eccentricity thermal (1.0)",1.0,1.6),("eccentricity uniform (0)",0.0,1.6),("eccentricity 2.0",2.0,1.6),
       ("separations log-flat (k 1.0)",1.3,1.0),("separations k 0.5",1.3,0.5),("separations k 2.2",1.3,2.2)]
H={(a,k):f_of(a,k) for _,a,k in cases}
mocks=[draw(*TRUE,36) for _ in range(NM)]
for lab,a,k in cases:
    gh=[]
    for vt,_ in mocks:
        noisy=vt*(1+0.07*rng.standard_normal(36))*(1+0.025*rng.standard_normal(36))
        gh.append(fit(noisy,H[(a,k)])[0])
    gh=np.array(gh); out["cases"][lab]=dict(mean=float(gh.mean()),median=float(np.median(gh)),sd=float(gh.std()),p16=float(np.percentile(gh,16)),p84=float(np.percentile(gh,84)))
    print(f"   {lab:34s} gamma = {np.median(gh):.3f}  (16-84%: {np.percentile(gh,16):.2f}-{np.percentile(gh,84):.2f})")
# semi-major axis tied to the projected separation (a = kappa s), eccentricity, phase and orientation still free:
# the likelihood is the library's v~ distribution among orbits with s/a = 1/kappa
for kap,lab in ((1.0,"a = s"),(4/np.pi,"a = (4/pi) s"),(np.sqrt(3/2),"a = sqrt(3/2) s")):
    L=LIB[1.3]; sel=np.abs(L['x']-1/kap)<0.015; vt=L['u'][sel]*np.sqrt(L['x'][sel]); hfix,_=np.histogram(vt,bins=GRID,density=True)
    gh=[]
    for v36,_ in mocks:
        noisy=v36*(1+0.07*rng.standard_normal(36))*(1+0.025*rng.standard_normal(36)); gh.append(fit(noisy,hfix)[0])
    gh=np.array(gh); out["cases"]["tied "+lab]=dict(median=float(np.median(gh)),p16=float(np.percentile(gh,16)),p84=float(np.percentile(gh,84)),n_lib=int(sel.sum()))
    print(f"   semi-major axis tied, {lab:16s} gamma = {np.median(gh):.3f}  (16-84%: {np.percentile(gh,16):.2f}-{np.percentile(gh,84):.2f})")
# how often is the true separation within 25% of the projected one (isotropic), the margin that decides a 1.12 'escaper'
for a in (1.3,):
    L=LIB[a]; w=L['x']**(TRUE[1]-1); rs=L['r']/L['x']
    out["p_r_within_25pct"]=float(np.sum(w*(rs<1.25))/np.sum(w)); print(f"   P(r/s < 1.25) with the true priors: {out['p_r_within_25pct']:.2f}")
# geometric short-cut: r = kappa s, circular: v^2 = G M / r  ->  gamma_i = v~^2 * kappa; combine as the median and as the mean
for kap,lab in ((1.0,"r = s"),(4/np.pi,"r = (4/pi) s, isotropic mean"),(np.sqrt(3/2),"r = sqrt(3/2) s, rms")):
    gm=[];gme=[]
    for vt,_ in mocks:
        noisy=vt*(1+0.07*rng.standard_normal(36)); gi=noisy**2*kap; gm.append(np.median(gi)); gme.append(np.mean(gi))
    out["cases"]["shortcut "+lab]=dict(median=float(np.median(gm)),mean=float(np.median(gme)),sd=float(np.std(gm)))
    print(f"   shortcut, {lab:28s} gamma (median of pairs) = {np.median(gm):.3f} ± {np.std(gm):.3f};  mean of pairs {np.median(gme):.3f}")
# 'Newtonian-unbound' pairs: how many of 36 exceed v_esc at an estimated r, with 7% velocity error, when r is taken as kappa*s or the true r
for kap,lab in ((None,"true r"),(1.0,"r = s"),(4/np.pi,"r = (4/pi) s")):
    cnt=[]
    for vt,rs in mocks:
        noisy=vt*(1+0.07*rng.standard_normal(36)); rr=rs if kap is None else kap   # r/s
        cnt.append(np.sum(noisy>np.sqrt(2/rr)))
    out[f"unbound_{lab}"]=dict(mean=float(np.mean(cnt)),p_ge4=float(np.mean(np.array(cnt)>=4)))
    print(f"   apparent escapers among 36 with {lab:14s}: mean {np.mean(cnt):.2f}; P(>= 4) = {np.mean(np.array(cnt)>=4):.3f}")
# hidden companions: a fraction fc of pairs carries an unseen star that adds 0-40% to the mass (not in the mass estimate) and moves the
# measured star by a kick of size kv (in units of sqrt(GM/s)), random direction. Fitted with the matched priors.
def contaminate(v36,rs,fc,kv):
    n=len(v36); c=rng.random(n)<fc; vt=v36.copy()
    boost=np.sqrt(1+0.4*rng.random(n)); vt[c]=vt[c]*boost[c]           # true orbit speed with the extra mass
    if kv>0:
        d=rng.standard_normal((n,3)); d/=np.linalg.norm(d,axis=1)[:,None]; ang=rng.standard_normal((n,3)); ang/=np.linalg.norm(ang,axis=1)[:,None]
        vec=vt[:,None]*ang; vec[c]+=kv*np.exp(0.5*rng.standard_normal((c.sum(),1)))*d[c]; vt=np.linalg.norm(vec,axis=1)
    return vt
out["companions"]={}; hm=H[(1.3,1.6)]
for fc in (0.1,0.2,0.3):
    for kv in (0.0,0.3,0.6):
        gh=[];esc=[]
        for v36,rs in mocks[:150]:
            vt=contaminate(v36,rs,fc,kv); noisy=vt*(1+0.07*rng.standard_normal(36))*(1+0.025*rng.standard_normal(36))
            gh.append(fit(noisy,hm)[0]); esc.append(np.sum(noisy>np.sqrt(2/rs)))
        gh=np.array(gh); esc=np.array(esc); key=f"fc{fc}_kv{kv}"
        out["companions"][key]=dict(fc=fc,kv=kv,gamma=float(np.median(gh)),p16=float(np.percentile(gh,16)),p84=float(np.percentile(gh,84)),esc=float(esc.mean()),p_esc4=float(np.mean(esc>=4)))
        print(f"   hidden companions {fc:.0%}, kick {kv:.1f} v_c: gamma = {np.median(gh):.2f} ({np.percentile(gh,16):.2f}-{np.percentile(gh,84):.2f}); escapers {esc.mean():.1f}, P(>=4) {np.mean(esc>=4):.2f}")
json.dump(out,open(os.path.join(HERE,"wb36.json"),"w"),indent=1)
