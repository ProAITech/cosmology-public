"""Workbench Entry 36: is dark energy changing? DESI's signal rebuilt from published numbers, and taken apart.
Models: flat LCDM (w = -1), and w(a) = w0 + wa (1 - a). Optional: one free factor alpha on the sound horizon (both at the CMB and at the
ruler), which stands for any change to physics before recombination, such as the clumping from primordial magnetic fields.
Significance: the chi2 gain for two extra numbers, turned into a Gaussian sigma the way DESI does."""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
import numpy as np
from scipy import optimize
from scipy.stats import chi2 as _c2, norm as _nm
import cosmo as c
sig=lambda dchi,dof=2: float(_nm.isf(_c2.sf(max(dchi,0),dof)/2)) if dchi>0 else 0.0
def make(parts):
    """parts: list of functions model -> chi2"""
    return lambda m: sum(p(m) for p in parts)
def fit(chi,kind,start=None,free_alpha=False,fix_ob=None):
    """kind: 'lcdm' or 'w0wa'. Returns best chi2 and parameters."""
    def unpack(x):
        ob=fix_ob if fix_ob is not None else x[0]*0.02223; i=0 if fix_ob is not None else 1
        obc=x[i]*0.142; h=x[i+1]*0.68; i+=2
        w0,wa=(-1.0,0.0)
        if kind=='w0wa': w0,wa=x[i],x[i+1]; i+=2
        al=x[i] if free_alpha else 1.0
        return ob,obc,h,w0,wa,al
    def f(x):
        ob,obc,h,w0,wa,al=unpack(x)
        if not (0.01<ob<0.04 and 0.05<obc<0.3 and 0.4<h<1.0 and -3<w0<1 and -3<wa<2 and w0+wa<0 and 0.8<al<1.2): return 1e9
        m=c.Model(ob,obc,h,w0,wa,al)
        if m.ode<=0: return 1e9
        return chi(m)
    best=None
    starts=start or ([[]])
    seeds=[(1.0,1.0),(0.98,1.01),(1.02,0.99)]
    wseeds=[(-1.0,0.0),(-0.7,-1.0),(-0.4,-1.8),(-0.9,-0.4)] if kind=='w0wa' else [None]
    for s in seeds:
        for ws in wseeds:
            x0=([] if fix_ob is not None else [1.0])+[s[0],s[1]]+(list(ws) if ws else [])+([1.0] if free_alpha else [])
            r=optimize.minimize(f,x0,method='Nelder-Mead',options=dict(xatol=1e-5,fatol=1e-5,maxiter=6000,maxfev=6000))
            r=optimize.minimize(f,r.x,method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-6,maxiter=6000,maxfev=6000))
            if best is None or r.fun<best.fun: best=r
    ob,obc,h,w0,wa,al=unpack(best.x); m=c.Model(ob,obc,h,w0,wa,al)
    return dict(chi=float(best.fun),ob=float(ob),obc=float(obc),h=float(h),w0=float(w0),wa=float(wa),alpha=float(al),Om=float((obc+c.OMNU0)/h**2),hrd=float(h*m.rd()),rd=float(m.rd()))
def compare(label,parts,fix_ob=None,free_alpha=False,quiet=False):
    chi=make(parts); L=fit(chi,'lcdm',fix_ob=fix_ob,free_alpha=free_alpha); W=fit(chi,'w0wa',fix_ob=fix_ob,free_alpha=free_alpha); d=L['chi']-W['chi']
    if not quiet: print(f"  {label:62s} LCDM chi2 {L['chi']:8.2f} (Om {L['Om']:.4f}, H0 {100*L['h']:.2f}"+(f", alpha {L['alpha']:.4f}" if free_alpha else "")+f") | w0wa chi2 {W['chi']:8.2f}: w0 {W['w0']:+.2f}, wa {W['wa']:+.2f} | gain {d:5.2f} = {sig(d):.1f} sigma")
    return dict(lcdm=L,w0wa=W,gain=float(d),sigma=sig(d))
if __name__=="__main__":
    out={}
    bao=lambda m: c.chi_bao(m); cmb=lambda m: c.chi_cmb(m); jla=lambda m: c.chi_jla(m); pan=lambda m: c.chi_pan(m)
    print("1. Rebuilding DESI's numbers")
    out["bao"]=compare("ruler alone",[bao],fix_ob=0.02223)
    print(f"      LCDM, ruler alone: matter share {out['bao']['lcdm']['Om']:.4f}, H0 x r_d = {out['bao']['lcdm']['hrd']:.2f}   (DESI: 0.2975 +- 0.0086, 101.54 +- 0.73; gain 4.7 = 1.7 sigma)")
    out["bao_cmb"]=compare("ruler + CMB prior on (theta*, omega_b, omega_bc)",[bao,cmb])
    print("      (DESI: gain 8.0 = 2.4 sigma, w0 = -0.43 +- 0.22, wa = -1.72 +- 0.64; with the full CMB: matter share 0.3027 +- 0.0036, H0 68.17 +- 0.28)")
    print("\n2. Adding supernovae")
    out["bao_cmb_jla"]=compare("ruler + CMB + 740 supernovae (2014, full errors)",[bao,cmb,jla])
    out["bao_cmb_pan"]=compare("ruler + CMB + 1,048 supernovae (2018, statistical errors)",[bao,cmb,pan])
    out["bao_cmb_pan_far"]=compare("   the same, supernovae beyond z = 0.1 only",[bao,cmb,lambda m: c.chi_pan(m,sel=c.ZP>0.1)])
    out["bao_cmb_pan_step"]=compare("   the same, all supernovae, free step at z = 0.1",[bao,cmb,lambda m: c.chi_pan(m,step=0.1)])
    out["bao_cmb_jla_far"]=compare("ruler + CMB + 2014 supernovae beyond z = 0.1 only",[bao,cmb,lambda m: c.chi_jla(m,sel=c.sh.zb>0.1)])
    out["cmb_pan"]=compare("CMB + 2018 supernovae, no ruler",[cmb,pan])
    out["cmb_jla"]=compare("CMB + 2014 supernovae, no ruler",[cmb,jla])
    print("\n3. Which ruler measurement carries it? (ruler + CMB, one tracer left out)")
    out["drop"]={}
    for i,nm in enumerate(c.BAO_NAMES):
        out["drop"][nm]=compare(f"without {nm}",[lambda m,i=i: c.chi_bao(m,drop=i),cmb])
    print("\n4. Which part of the CMB prior carries it?")
    out["cmb_theta_ob"]=compare("ruler + (theta*, omega_b) only: matter density not from the CMB",[bao,lambda m: c.chi_cmb(m,which=(0,1))])
    out["cmb_ob_obc"]=compare("ruler + (omega_b, omega_bc) only: no acoustic angle",[bao,lambda m: c.chi_cmb(m,which=(1,2))])
    print("\n5. A free sound horizon (alpha multiplies it at the CMB and at the ruler alike)")
    out["alpha"]=compare("ruler + CMB, alpha free",[bao,cmb],free_alpha=True)
    out["alpha_pan"]=compare("ruler + CMB + 2018 supernovae, alpha free",[bao,cmb,pan],free_alpha=True)
    out["alpha_jla"]=compare("ruler + CMB + 2014 supernovae, alpha free",[bao,cmb,jla],free_alpha=True)
    a=out["alpha"]["lcdm"]; b=out["bao_cmb"]["lcdm"]; print(f"      LCDM with alpha free: alpha = {a['alpha']:.4f} (sound horizon {100*(a['alpha']-1):+.2f}%), chi2 {a['chi']:.2f} against {b['chi']:.2f} with alpha = 1: gain {b['chi']-a['chi']:.2f} for one number = {sig(b['chi']-a['chi'],1):.1f} sigma; H0 {100*a['h']:.2f} against {100*b['h']:.2f}")
    out["alpha_gain"]=dict(gain=float(b['chi']-a['chi']),sigma=sig(b['chi']-a['chi'],1))
    # error on alpha by profile
    chi=make([bao,cmb])
    def prof(al):
        def f(x):
            ob,obc,h=x[0]*0.02223,x[1]*0.142,x[2]*0.68
            if not (0.01<ob<0.04 and 0.05<obc<0.3 and 0.4<h<1.0): return 1e9
            m=c.Model(ob,obc,h,-1,0,al); return chi(m) if m.ode>0 else 1e9
        r=optimize.minimize(f,[a['ob']/0.02223, a['obc']/0.142, a['h']/0.68],method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-6,maxiter=4000)); return r.fun
    lo=optimize.brentq(lambda v: prof(v)-a['chi']-1,a['alpha']-0.03,a['alpha']); hi=optimize.brentq(lambda v: prof(v)-a['chi']-1,a['alpha'],a['alpha']+0.03)
    print(f"      alpha = {a['alpha']:.4f} +{hi-a['alpha']:.4f} -{a['alpha']-lo:.4f}"); out["alpha_fit"]=[a['alpha'],a['alpha']-lo,hi-a['alpha']]
    json.dump(out,open(os.path.join(HERE,"dark_energy_fit.json"),"w"),indent=1)

