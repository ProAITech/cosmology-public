"""Workbench Entry 29: how many deuterons does a cosmic ray make when it stops in dense hydrogen and helium?
The helium-from-stars model (Lerner's GOLE paper, section 4.3) says nearly 0.9 deuterons per cosmic-ray proton, by this route:
p + p -> p + n + pi+ makes a neutron, and in dense gas the neutron slows down and is captured by a proton, n + p -> d + gamma.
The paper gives no calculation. This is one: a Monte Carlo of the cascade a cosmic ray starts in a thick target.

Followed: protons, neutrons, deuterons, and incoming alpha particles. Processes:
  protons    slow by ionisation (Bethe); p-p elastic; p-p inelastic above 280 MeV (pion production; a neutron comes out in about 4 of 5, or a deuteron
             directly, p p -> d pi+); p + 4He break-up (frees neutrons and deuterons).
  neutrons   n-p elastic (the slowing-down); n-p inelastic; n + 4He break-up; once slow, capture on hydrogen (a deuteron) or on 3He (lost).
  deuterons  slow by ionisation; break-up on protons and helium, which frees the neutron again.
Cross sections (mb; kinetic energy T in MeV):
  p-p inelastic: Kafexhiu et al. (2014) formula. p p -> d pi+: table through the measured peak, 3.15 mb near 580 MeV.
  n-p, p-p elastic and n-p inelastic: tables of standard free-nucleon values (see below). 4He + p channels: Glagolev et al. (bubble chamber, 620 MeV per nucleon):
  3He p n 23.0, 3H p p 17.1, d d p 3.2, d p p n 22.0, p p p n n 11.5, pion channels 12.7; total inelastic 89.5.
Everything not measured directly (how energy is shared, low-energy helium channels, deuteron break-up) is varied to show how much it matters."""
import numpy as np, json, sys
me=0.510999; mp=938.272; mn=939.565; md=1875.613; mpi=139.57; K=0.307075; NA=6.02214e23
rng=np.random.default_rng(11)
T_el=np.array([1,2,5,10,20,30,50,100,150,200,300,400,600,1000,2000,5000,1e5])
s_np_el=np.array([4250,2900,1630,940,480,300,165,73,50,42,35,32,30,24,17,12,7.5])
s_pp_el=np.array([4250,2900,1630,380,150,95,55,33,26,23.5,23,24,24.5,24,20,13,7.5]); s_pp_el[:3]=[0,0,0]
T_d=np.array([279.7,288,320,350,400,450,500,550,580,620,650,700,750,800,900,1000,1200,1500,2000,5000,1e5]); s_dpi=np.array([0,0.02,0.15,0.4,0.9,1.7,2.5,3.0,3.15,3.0,2.7,2.0,1.5,1.1,0.6,0.35,0.17,0.08,0.04,0.01,0.0])
T_he=np.array([23,30,45,60,100,200,400,1e5]); s_he=np.array([0,40,105,110,110,100,90,90])
e_bu=np.array([3.4,5,8,12,20,30,50,100,200,500,1000,1e5]); s_bu=np.array([0,60,140,170,160,135,110,82,66,57,55,55])
def s_inel_pp(T):
    T=np.asarray(T,float); x=np.maximum(T/279.7,1.0); L=np.log(x); return np.where(T>279.7,(30.7-0.96*L+0.18*L*L)*(1-x**-1.9)**3,0.0)
r_np=lambda T: 0.45+0.5*np.clip((T-600)/1400,0,1)
lin=lambda x,xs,ys: np.interp(np.log(x),np.log(xs),ys)
class Medium:
    def __init__(s,Y=0.25,plasma=1.0,f_inel=1.0,f_he=1.0,f_bu=1.0,pn_frac=0.82,rec_mean=80.0,lead=(0.5,0.95)):
        s.Y=Y; X=1-Y; s.nH=X*NA/1.008; s.nHe=Y*NA/4.0026; s.pl=plasma; s.fi=f_inel; s.fh=f_he; s.fb=f_bu; s.pn=pn_frac; s.rec=rec_mean; s.lead=lead
        s.comp=[(X,1/1.008,19.2e-6),(Y,2/4.0026,41.8e-6)]; s.Tcut=23.0 if Y>0 else 279.7
        s.grid=np.geomspace(5,2e5,900)
        s.Sg=s.S(s.grid); s.mu_p=s.mu_proton(s.grid); s.Lp=np.concatenate([[0],np.cumsum(0.5*(s.mu_p[1:]/s.Sg[1:]+s.mu_p[:-1]/s.Sg[:-1])*np.diff(s.grid))])
        eg=s.grid; s.mu_d=(s.nH+4*s.nHe)*lin(np.maximum(eg/2,3.4),e_bu,s_bu)*1e-27*s.fb*(eg/2>3.4); Sd=s.S(eg/2)     # deuteron of kinetic energy eg: same speed as a proton of eg/2
        s.Ld=np.concatenate([[0],np.cumsum(0.5*(s.mu_d[1:]/Sd[1:]+s.mu_d[:-1]/Sd[:-1])*np.diff(eg))])
        mu_a=(s.nH*lin(np.maximum(eg,23),T_he,s_he)*s.fh+s.nHe*270.0*(eg>20))*1e-27; s.La=np.concatenate([[0],np.cumsum(0.5*(mu_a[1:]/s.Sg[1:]+mu_a[:-1]/s.Sg[:-1])*np.diff(eg))])   # alpha, per nucleon energy
    def S(s,T):
        g=1+np.asarray(T,float)/mp; b2=1-1/g**2; return s.pl*sum(w*K*za/b2*(np.log(2*me*b2*g*g/I)-b2) for w,za,I in s.comp if w>0)
    def parts(s,T,neutron):
        el=(lin(T,T_el,s_np_el) if neutron else lin(T,T_el,s_pp_el))*(T>(1.0 if neutron else 25.0))
        inel=s_inel_pp(T)*s.fi*(r_np(T) if neutron else 1.0); he=lin(max(T,23.0),T_he,s_he)*s.fh if (T>23 and s.Y>0) else 0.0
        return s.nH*el*1e-27, s.nH*float(inel)*1e-27, s.nHe*he*1e-27
    def mu_proton(s,T):
        el=lin(T,T_el,s_pp_el)*(T>25.0); he=lin(np.maximum(T,23.0),T_he,s_he)*s.fh*(T>23)
        return (s.nH*(el+s_inel_pp(T)*s.fi)+s.nHe*he)*1e-27
def d_pi_energy(T):                                   # lab kinetic energy of the deuteron in p p -> d pi+, angular distribution 1/3 + cos^2
    s=4*mp*mp+2*mp*T; rs=np.sqrt(s); Ed=(s+md*md-mpi*mpi)/(2*rs); pd=np.sqrt(max(Ed*Ed-md*md,0)); g=(T+2*mp)/rs; b=np.sqrt(1-1/g**2)
    while True:
        c=rng.uniform(-1,1)
        if rng.uniform(0,4/3)<1/3+c*c: break
    return max(g*(Ed+b*pd*c)-md,1.0)
HE=[("A",23.0),("B",17.1),("C",3.2),("D",22.0),("E",11.5),("F1",2.0),("F2",3.0),("F3",0.65),("F4",1.1),("F5",0.65),("F6",1.3),("F7",4.0)]
def he_channel(T):
    if T<150:
        pick=35*np.exp(-(T-35)/40); tot=float(lin(max(T,23.0),T_he,s_he))
        if rng.random()<min(pick/max(tot,1e-9),1.0): return "P"
        return rng.choice(["A","B","D","E"],p=[0.45,0.33,0.17,0.05])
    names=[n for n,_ in HE]; w=np.array([v for _,v in HE]); 
    if T<350: w=w*np.array([1,1,1,1,1,0,0,0,0,0,0,0])
    return rng.choice(names,p=w/w.sum())
def cascade(M,kind,T0):
    """returns (slow neutrons, deuterons left intact, 3He made, 3H made) for one primary"""
    nth=0; dsv=0; he3=0; h3=0; stack=[(kind,T0)]
    def push(k,T):
        nonlocal nth
        if k=="n":
            if T<=M.Tcut: nth+=1
            else: stack.append((k,T))
        elif k=="p":
            if T>M.Tcut: stack.append((k,T))
    while stack:
        k,T=stack.pop()
        if k=="d":
            L0=np.interp(T,M.grid,M.Ld); xi=rng.exponential()
            if xi>=L0: dsv+=1; continue
            T1=np.interp(L0-xi,M.Ld,M.grid); push("n",T1/2); push("p",T1/2); continue
        if k=="a":                                  # alpha particle, T = kinetic energy per nucleon
            L0=np.interp(T,M.grid,M.La); xi=rng.exponential()
            if xi>=L0: continue
            T1=np.interp(L0-xi,M.La,M.grid); ch=he_channel(T1); P3=1-np.exp(-np.interp(T1,M.grid,M.La)*65/90)   # chance that a mass-3 fragment breaks up too
            def frag3(nn,npr):
                nonlocal he3,h3
                if rng.random()<P3*(0.75 if npr==2 else 3.0 if P3*3<1 else 1/P3):
                    for _ in range(nn): push("n",T1)
                    for _ in range(npr): push("p",T1)
                elif npr==2: he3+=1
                else: h3+=1
            if ch=="P": stack.append(("d",2*T1)); frag3(1,2)
            elif ch=="A": push("n",T1); frag3(1,2)
            elif ch=="B": push("p",T1); frag3(2,1)
            elif ch=="C": stack.append(("d",2*T1)); stack.append(("d",2*T1))
            elif ch=="D": stack.append(("d",2*T1)); push("n",T1); push("p",T1)
            elif ch=="E": push("n",T1); push("n",T1); push("p",T1); push("p",T1)
            elif ch in ("F1","F7"): stack.append(("a",max(T1-60,1)))
            elif ch=="F2": push("n",T1); frag3(2,1)
            elif ch=="F3": stack.append(("d",2*T1)); stack.append(("d",2*T1))
            elif ch=="F4": stack.append(("d",2*T1)); push("n",T1); push("n",T1)
            elif ch=="F5": push("n",T1); push("n",T1); push("n",T1); push("p",T1)
            elif ch=="F6": stack.append(("d",2*T1)); frag3(1,2)
            continue
        neutron=(k=="n")
        if neutron: T1=T
        else:
            L0=np.interp(T,M.grid,M.Lp); xi=rng.exponential()
            if xi>=L0: continue
            T1=np.interp(L0-xi,M.Lp,M.grid)
            if T1<=M.Tcut: continue
        a,b,c=M.parts(T1,neutron); tot=a+b+c
        if tot<=0:
            if neutron: nth+=1
            continue
        u=rng.random()*tot
        if u<a:                                      # elastic on a proton
            Tr=rng.uniform(0,T1) if T1<2*M.rec else min(rng.exponential(M.rec),T1/2)
            if neutron and T1<400 and rng.random()<0.5: push("n",Tr); push("p",T1-Tr)     # backward (charge-exchange) scattering: the neutron comes out slow
            else: push(k,T1-Tr); push("p",Tr)
        elif u<a+b:                                  # pion production
            sd=float(lin(max(T1,279.7),T_d,s_dpi))*(0.5 if neutron else 1.0); si=float(s_inel_pp(T1))*M.fi*(r_np(T1) if neutron else 1.0)
            if rng.random()<sd/max(si,1e-9): stack.append(("d",d_pi_energy(T1))); continue
            fpi=rng.uniform(0.2,0.4); Kn=max(T1*(1-fpi),1.0); x=rng.uniform(*M.lead); Ta,Tb=Kn*x,Kn*(1-x)
            if neutron:
                v=rng.random(); kinds=("n","n") if v<0.2 else ("p","p") if v<0.4 else (("n","p") if rng.random()<0.5 else ("p","n"))
            else:
                if T1<1500: v=rng.random(); kinds=(("n","p") if rng.random()<0.5 else ("p","n")) if v<M.pn else ("p","p")
                else: v=rng.random(); kinds=("p","p") if v<0.30 else ("n","n") if v<0.40 else (("n","p") if rng.random()<0.5 else ("p","n"))
            push(kinds[0],Ta); push(kinds[1],Tb)
        else:                                        # break-up of a helium nucleus
            ch=he_channel(T1); Tr=min(rng.exponential(60.0),0.6*T1); Tp=max(T1-28-Tr,1.0); o="p" if neutron else "n"   # o: the other kind of nucleon
            if ch=="P": dsv+=1; (he3,h3)=(he3,h3+1) if neutron else (he3+1,h3)
            elif ch=="A": push(k,Tp); push(o,Tr); (he3,h3)=(he3,h3+1) if neutron else (he3+1,h3)
            elif ch=="B": push(k,Tp); push(k,Tr); (he3,h3)=(he3+1,h3) if neutron else (he3,h3+1)
            elif ch=="C": push(k,Tp); dsv+=2
            elif ch=="D":
                push(k,Tp); dsv+=1
                if rng.random()<0.5: push(k,Tr); push(o,8.0)
                else: push(o,Tr); push(k,8.0)
            elif ch=="E": push(k,Tp); push(k,Tr); push(k,8.0); push(o,8.0); push(o,8.0)
            elif ch=="F1": push(o,max(Tp-200,1.0))
            elif ch=="F2": push(k,max(Tp-200,1.0)); push(o,Tr); (he3,h3)=(he3+1,h3) if neutron else (he3,h3+1)
            elif ch=="F3": push(o,max(Tp-200,1.0)); dsv+=2
            elif ch=="F4": push(o,max(Tp-200,1.0)); push(o,Tr); dsv+=1
            elif ch=="F5": push(o,max(Tp-200,1.0)); push(o,Tr); push(o,8.0); push(k,8.0); push(k,8.0)
            elif ch=="F6": dsv+=1
            else: push(k,max(Tp-200,1.0))
    return nth,dsv,he3,h3
def yield_curve(M,kind,energies,n=1500):
    out=[]
    for T in energies:
        r=np.array([cascade(M,kind,T) for _ in range(n)]); out.append(r.mean(0))
    return np.array(out)
if __name__=="__main__":
    quick=len(sys.argv)>1 and sys.argv[1]=="quick"
    E=np.array([50,100,200,300,400,500,600,800,1000,1500,2000,3000,5000,10000,20000.]); nper=400 if quick else (int(sys.argv[1]) if len(sys.argv)>1 else 2500)
    spec=lambda T: (T+980.0)**-2.7                         # mean 1.4 GeV per nucleon, as in Entry 25
    Tf=np.geomspace(10,1e5,400); w=spec(Tf)*np.gradient(Tf); w/=w.sum()
    def averaged(M,kind="p"):
        yc=yield_curve(M,kind,E,nper); avg=[float((w*np.interp(np.log(Tf),np.log(E),yc[:,j],left=0)).sum()) for j in range(4)]; return yc,avg
    out={"E":E.tolist(),"cases":{}}
    cases=[("hydrogen only",dict(Y=0.0)),("hydrogen and helium (Y = 0.25)",dict(Y=0.25)),("ionised gas, stops twice as fast",dict(Y=0.25,plasma=2.0)),("inelastic p-p x0.8",dict(Y=0.25,f_inel=0.8)),("inelastic p-p x1.2",dict(Y=0.25,f_inel=1.2)),
           ("neutron in 70% of p-p pion events",dict(Y=0.25,pn_frac=0.70)),("helium break-up x0.7",dict(Y=0.25,f_he=0.7)),("deuteron break-up x2",dict(Y=0.25,f_bu=2.0)),("leading nucleon keeps more",dict(Y=0.25,lead=(0.7,0.95))),("elastic recoil 150 MeV",dict(Y=0.25,rec_mean=150.0))]
    if quick: cases=cases[:2]
    print("deuterons per cosmic-ray PROTON, spectrum (E+980 MeV)^-2.7 (mean 1.4 GeV):  slow neutrons | deuterons made directly and left intact | 3He | 3H | total if every slow neutron is captured by hydrogen")
    for name,kw in cases:
        M=Medium(**kw); yc,avg=averaged(M); out["cases"][name]={"curve":yc.tolist(),"avg":avg}
        print(f"  {name:36s} {avg[0]:.3f} | {avg[1]:.3f} | {avg[2]:.3f} | {avg[3]:.3f} | {avg[0]+avg[1]:.2f}")
        if name in ("hydrogen only","hydrogen and helium (Y = 0.25)"):
            print("       by energy (MeV: neutrons+deuterons per proton): "+", ".join(f"{int(e)}: {a+b:.2f}" for e,(a,b,_,_) in zip(E,yc)))
    M=Medium(Y=0.25); yc,avg=averaged(M,"a"); out["alpha"]={"curve":yc.tolist(),"avg":avg}
    print(f"per cosmic-ray ALPHA (same spectrum per nucleon), hydrogen and helium: {avg[0]:.2f} slow neutrons + {avg[1]:.2f} deuterons = {avg[0]+avg[1]:.2f}; 3He {avg[2]:.2f}, 3H {avg[3]:.2f}")
    json.dump(out,open("deuterium_yield.json","w"),indent=1)
