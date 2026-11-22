"""Entry 149: is strong lensing gravity or plasma refraction? SLACS (Bolton et al. 2008, grade-A lenses with SDSS sigma, 58)
and S4TM (Shu et al. 2017, 40). Columns: z_lens, z_source, sigma (km/s), e_sigma, b_SIE (arcsec).
(1) Gravity, no free parameters: an isothermal lens with the measured stellar velocity dispersion has
    theta_E = 4 pi (sigma/c)^2 D_ls/D_s. Distances: flat LCDM, or static Euclidean space with tired light (D proportional to ln(1+z)).
(2) Plasma refraction: a plasma deflects by alpha = (lambda^2 r_e / 2 pi) dN_e/dx, so the electron column that must change
    across the Einstein radius is Delta N_e ~ alpha * 2 pi R_E / (lambda^2 r_e), at the HST F814W wavelength in the lens frame."""
import numpy as np, json
def load(f): return np.array([list(map(float,r.split(','))) for r in open(f).read().strip().split(';')])
c=299792.458; Om=0.31; H0=67.7; Mpc=3.0857e24; re_=2.8179403e-13; arcsec=np.pi/180/3600
def Dc(z):
    zz=np.linspace(0,z,2001); return c/H0*np.trapezoid(1/np.sqrt(Om*(1+zz)**3+1-Om),zz)
res={}
allrows=[]
for name,f in (('SLACS',"slacs.txt"),('S4TM',"s4tm.txt")):
    A=load(f); zl,zs,sig,es,b=A.T
    Dl=np.array([Dc(z) for z in zl]); Ds=np.array([Dc(z) for z in zs])
    rat_L=(Ds-Dl)/Ds                         # flat LCDM: D_ls/D_s with comoving distances
    rat_S=(np.log(1+zs)-np.log(1+zl))/np.log(1+zs)   # static Euclidean, tired-light distances
    thS=4*np.pi*(sig/c)**2/arcsec
    pred_L=thS*rat_L; pred_S=thS*rat_S
    fL=np.sqrt(b/pred_L); fS=np.sqrt(b/pred_S)     # sigma_SIE / sigma_star
    efs=fL*es/sig                                 # error on f from sigma alone
    # slope of log(b * Ds/Dls) against log sigma
    y=np.log10(b/rat_L); x=np.log10(sig); p=np.polyfit(x,y,1)
    # bootstrap slope
    rng=np.random.default_rng(149); bs=[]
    for _ in range(4000):
        k=rng.integers(0,len(x),len(x)); bs.append(np.polyfit(x[k],y[k],1)[0])
    # fraction of scatter explained by sigma errors
    resid=np.log10(b/pred_L)
    # plasma: deflection at the lens, column change across R_E
    alpha=b*arcsec/rat_L                         # radians
    lam=814e-7/(1+zl)                             # cm, F814W in the lens frame
    RE=b*arcsec*Dl/(1+zl)*Mpc                     # proper Einstein radius, cm (angular-diameter distance = Dc/(1+z))
    dN=alpha*2*np.pi*RE/(lam**2*re_)
    lam_r=21.0/(1+zl); dN_r=alpha*2*np.pi*RE/(lam_r**2*re_)
    res[name]=dict(n=int(len(A)),f_LCDM_mean=float(fL.mean()),f_LCDM_sem=float(fL.std(ddof=1)/np.sqrt(len(A))),f_LCDM_rms=float(fL.std(ddof=1)),
                   f_static_mean=float(fS.mean()),f_static_rms=float(fS.std(ddof=1)),med_efrac=float(np.median(efs)),
                   slope=float(p[0]),slope_lo=float(np.percentile(bs,16)),slope_hi=float(np.percentile(bs,84)),
                   resid_rms_dex=float(resid.std(ddof=1)),
                   RE_kpc_med=float(np.median(RE/3.0857e21)),alpha_med_arcsec=float(np.median(alpha/arcsec)),
                   dN_opt_med=float(np.median(dN)),dN_radio_med=float(np.median(dN_r)),
                   sigma_range=[float(sig.min()),float(sig.max())],zs_range=[float(zs.min()),float(zs.max())],b_range=[float(b.min()),float(b.max())])
    for i in range(len(A)): allrows.append([name,float(sig[i]),float(es[i]),float(b[i]),float(pred_L[i]),float(rat_L[i])])
for k,v in res.items(): print(k,{kk:(round(vv,4) if isinstance(vv,float) else vv) for kk,vv in v.items()})
# combined
fl=np.array([np.sqrt(r[3]/r[4]) for r in allrows]); print('combined f mean %.3f +- %.3f rms %.3f n %d'%(fl.mean(),fl.std(ddof=1)/np.sqrt(len(fl)),fl.std(ddof=1),len(fl)))
# Bounds on electron columns: if every baryon in a 1e13 Msun group sat inside R_E (median), electron column:
Msun=1.989e33; mp=1.6726e-24
for nm in res:
    RE=res[nm]['RE_kpc_med']*3.0857e21
    res[nm]['N_all_group']=1e13*Msun/mp/(np.pi*RE**2)
    print(nm,'column if 1e13 Msun of hydrogen sat inside R_E: %.2e cm^-2; needed %.2e (optical), %.2e (21 cm)'%(res[nm]['N_all_group'],res[nm]['dN_opt_med'],res[nm]['dN_radio_med']))
# chromatic factor between HST F555W/F606W and F814W, and optical vs 21 cm
res['chroma']=dict(F814_over_F555=(814/555)**2,F160_over_F814=(1600/814)**2,radio21_over_F814=(21.0/814e-7)**2)
res['combined']=dict(n=int(len(fl)),f_mean=float(fl.mean()),f_sem=float(fl.std(ddof=1)/np.sqrt(len(fl))),f_rms=float(fl.std(ddof=1)))
res['fig']=[[r[1],r[3],round(r[4],3),0 if r[0]=='SLACS' else 1] for r in allrows]
print(res['chroma'])
json.dump(res,open('plens.json','w'))
