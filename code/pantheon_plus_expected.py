"""Entry 40: what each curve predicts for the power q measured inside one survey. A curve that bends (the standard one) gives a different q
in each redshift range; a pure power law (Lerner's static law, q = 1/2) gives the same q everywhere."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy import integrate, optimize
A=np.loadtxt(os.path.join(HERE,"pantheon_plus_rows.txt")); z,m,e=A.T; w=1/e**2
lab=np.array(list("".join("".join(l.split()[1:]) for l in open(os.path.join(HERE,"survey_labels.txt")) if l.startswith("c"))))
zg=np.linspace(0,2.4,2401)
def mu_l(om): dc=integrate.cumulative_trapezoid(1/np.sqrt(om*(1+zg)**3+1-om),zg,initial=0); return 5*np.log10((1+z)*np.interp(z,zg,dc))
def qfit(y,s):
    def f(q):
        r=(y-5*np.log10(z*(1+z)**q))[s]; a=np.sum(w[s]*r)/w[s].sum(); return float(np.sum(w[s]*(r-a)**2))
    g=optimize.minimize_scalar(f,bounds=(-3,5),method="bounded",options={"xatol":1e-6}); a=optimize.brentq(lambda x:f(x)-g.fun-1,-3,g.x); b=optimize.brentq(lambda x:f(x)-g.fun-1,g.x,5); return float(g.x),float((b-a)/2)
out={}; rows=[]
print("survey: measured q | q the standard curve predicts there (matter share 0.30 / 0.33) | Lerner 0.5")
for k,nm in (("S","SDSS"),("P","Pan-STARRS Medium Deep"),("D","DES"),("N","SNLS")):
    s=lab==k; q,eq=qfit(m,s); p30=qfit(mu_l(0.30),s)[0]; p33=qfit(mu_l(0.33),s)[0]; zm=float(np.exp(np.sum(w[s]*np.log(z[s]))/w[s].sum()))
    print(f"   {nm:24s} z {z[s].min():.2f}-{z[s].max():.2f} (typical {zm:.2f}): {q:.3f} +- {eq:.3f} | {p30:.3f} / {p33:.3f} | 0.5"); rows.append([nm,float(z[s].min()),float(z[s].max()),zm,q,eq,p30,p33])
out["rows"]=rows; R=np.array([[r[4],r[5],r[6],r[7]] for r in rows])
for nm,col in (("standard, matter share 0.30",2),("standard, matter share 0.33",3)):
    c2=float(np.sum(((R[:,0]-R[:,col])/R[:,1])**2)); print(f"   chi2 of the four against {nm}: {c2:.1f}"); out["chi_"+nm[-4:]]=c2
c2=float(np.sum(((R[:,0]-0.5)/R[:,1])**2)); print(f"   chi2 of the four against Lerner's 0.5: {c2:.1f}"); out["chi_lerner"]=c2
c2=float(np.sum(((R[:,0]-1.0)/R[:,1])**2)); print(f"   chi2 of the four against q = 1 (stretched time): {c2:.1f}"); out["chi_one"]=c2
wm=np.sum(R[:,0]/R[:,1]**2)/np.sum(1/R[:,1]**2); print(f"   weighted mean of the four: {wm:.3f} +- {np.sum(1/R[:,1]**2)**-0.5:.3f}"); out["mean"]=[float(wm),float(np.sum(1/R[:,1]**2)**-0.5)]
# trend of q with the redshift each survey covers
lo=R[:2]; hi=R[2:]; a=np.sum(lo[:,0]/lo[:,1]**2)/np.sum(1/lo[:,1]**2); ea=np.sum(1/lo[:,1]**2)**-0.5; b=np.sum(hi[:,0]/hi[:,1]**2)/np.sum(1/hi[:,1]**2); eb=np.sum(1/hi[:,1]**2)**-0.5
print(f"   the two nearer surveys {a:.3f} +- {ea:.3f}; the two farther {b:.3f} +- {eb:.3f}; difference {a-b:.3f} +- {np.hypot(ea,eb):.3f}; the standard curve predicts {np.mean(lo[:,2])-np.mean(hi[:,2]):.3f}")
out["trend"]=[float(a),float(ea),float(b),float(eb)]
# the camera link: what each curve predicts for Foundation minus Pan-STARRS distant, and what is seen
F=(lab=="F"); P=(lab=="P")&(z>0.1)
def link(y):
    r=m-y; return float(np.sum(w[F]*r[F])/w[F].sum()-np.sum(w[P]*r[P])/w[P].sum())
el=float((1/w[F].sum()+1/w[P].sum())**0.5)
for nm,y in (("standard 0.30",mu_l(0.30)),("standard 0.33",mu_l(0.33)),("Lerner",5*np.log10(z*np.sqrt(1+z))),("ln(1+z) stretched",5*np.log10(np.log(1+z)*(1+z))),("coasting",5*np.log10(z*(1+z/2)))):
    v=link(y); print(f"   Foundation minus Pan-STARRS distant against {nm:18s}: {v:+.3f} +- {el:.3f}"); out["link_"+nm]=[v,el]
json.dump(out,open(os.path.join(HERE,"pantheon_plus_expected.json"),"w"),indent=1)
