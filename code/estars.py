"""Entry 145: are stars powered from outside? Test with Eker et al. (2018) detached eclipsing binaries.
Input: eker_pairs.txt = per pair 'Mp,Ms,dlogL,dlogTeff,e_dlogL,e_dlogTeff' (primary minus secondary), 291 pairs, 0.1-30 Msun.
Built in the browser from VizieR J/MNRAS/479/5491/table1 (logL = 2 log R + 4 log(Teff/5772))."""
import json, numpy as np, sys
rows=np.array([list(map(float,r.split(','))) for r in open(sys.argv[1] if len(sys.argv)>1 else 'eker_pairs.txt').read().split(';')])
Mp,Ms,dL,dT,eL,eT=rows.T
def Lz(M):  # Tout et al. 1996 ZAMS luminosity, Z = 0.02
    return (0.39704170*M**5.5+8.52762600*M**11)/(0.00025546+M**3+5.43288900*M**5+5.56357900*M**7+0.78866060*M**8+0.00586685*M**9.5)
assert abs(Lz(1.0)-0.698)<0.002 and abs(Lz(2.0)-15.9)<0.2
dM=np.log10(Mp/Ms); dLz=np.log10(Lz(Mp))-np.log10(Lz(Ms)); dR=(dL-4*dT)/2
res=dL-dLz
out={}
out['n_pairs']=int(len(rows))
sel=np.abs(dM)>0.02
out['n_massdiff']=int(sel.sum())
out['slope_obs_origin']=float((dM[sel]*dL[sel]).sum()/(dM[sel]**2).sum())
out['slope_zams_origin']=float((dM[sel]*dLz[sel]).sum()/(dM[sel]**2).sum())
med=np.median(res); mad=np.median(np.abs(res-med))
out['res_zams_median']=float(med); out['res_zams_mad']=float(mad)
out['res_zams_gt03']=int((np.abs(res)>0.3).sum())
out['frac_res_pos']=float((res[sel]>0).mean())
out['median_eL']=float(np.median(eL[eL>0]))
# sign agreement: brighter star is the heavier star
s=sel
out['sign_agree']=int((np.sign(dL[s])==np.sign(dM[s])).sum()); out['sign_n']=int(s.sum())
# sign agreement for |dM|>0.1
s2=np.abs(dM)>0.1
out['sign_agree_q']=int((np.sign(dL[s2])==np.sign(dM[s2])).sum()); out['sign_n_q']=int(s2.sum())
# equal-surface-brightness (shared external current) null: dlogT = 0
q=np.abs(dM)>np.log10(1/0.8)
ok=q&(eT>0)
z=dT[ok]/eT[ok]
out['n_q08']=int(q.sum()); out['n_q08_witherr']=int(ok.sum())
out['med_abs_dT_q08']=float(np.median(np.abs(dT[q]))); out['Tratio_med']=float(10**np.median(np.abs(dT[q])))
out['n_gt3sig']=int((np.abs(z)>3).sum()); out['chi2_eqT']=float((z**2).sum())
# hotter star is the heavier one?
out['hotter_is_heavier_q08']=int((np.sign(dT[q])==np.sign(dM[q])).sum())
# regression dL = a dM + b dR  (fusion homology: a~3-5, b<=0; shared external current: a=0, b=2)
X=np.vstack([dM[sel],dR[sel]]).T; y=dL[sel]
coef,*_=np.linalg.lstsq(X,y,rcond=None)
resid=y-X@coef; s2r=resid@resid/(len(y)-2); cov=s2r*np.linalg.inv(X.T@X)
out['reg_a']=float(coef[0]); out['reg_b']=float(coef[1]); out['reg_ea']=float(np.sqrt(cov[0,0])); out['reg_eb']=float(np.sqrt(cov[1,1]))
out['corr_dM_dR']=float(np.corrcoef(dM[sel],dR[sel])[0,1])
# restricted to pairs where both near main sequence, |res|<0.3
m=sel&(np.abs(res)<0.3)
X2=np.vstack([dM[m],dR[m]]).T; c2,*_=np.linalg.lstsq(X2,dL[m],rcond=None)
r2=dL[m]-X2@c2; cv2=(r2@r2/(m.sum()-2))*np.linalg.inv(X2.T@X2)
out['reg2_a']=float(c2[0]); out['reg2_b']=float(c2[1]); out['reg2_ea']=float(np.sqrt(cv2[0,0])); out['reg2_eb']=float(np.sqrt(cv2[1,1])); out['reg2_n']=int(m.sum())
# Same-radius test: pairs with |dR|<0.03 but |dM|>0.05: external model says dL ~ 0
sr=(np.abs(dR)<0.03)&(np.abs(dM)>0.05)
out['sameR_n']=int(sr.sum()); out['sameR_med_dL']=float(np.median(dL[sr]*np.sign(dM[sr]))) if sr.sum() else None
out['sameR_med_dLz']=float(np.median(dLz[sr]*np.sign(dM[sr]))) if sr.sum() else None
for k,v in out.items(): print(k, round(v,4) if isinstance(v,float) else v)
# figure data
fig={'pts':[[round(float(a),4),round(float(b),3),round(float(c),3),round(float(d),3)] for a,b,c,d in zip(dM,dL,dLz,dR)]}
json.dump({'summary':out,'fig':fig},open('estars_results.json','w'))

# ---- Part 2: Teff difference predicted by the fusion ZAMS (Tout et al. 1996 radius fit, Z = 0.02) ----
def Rz(M):
    return (1.71535900*M**2.5+6.59778800*M**6.5+10.08855000*M**11+1.01249500*M**19+0.07490166*M**19.5)/(0.01077422+3.08223400*M**2+17.84778000*M**8.5+M**18.5+0.00022582*M**19.5)
assert abs(Rz(1.0)-0.888)<0.002
dRz=np.log10(Rz(Mp))-np.log10(Rz(Ms)); dTz=(dLz-2*dRz)/4
sg=np.sign(dM); sg[sg==0]=1
X=(dTz*sg)[sel]; Y=(dT*sg)[sel]
out2={}
out2['n']=int(sel.sum())
out2['dTz_med']=float(np.median(X)); out2['dT_med']=float(np.median(Y))
out2['slope_T']=float((X*Y).sum()/(X*X).sum())
out2['med_res_T']=float(np.median(Y-X)); out2['mad_res_T']=float(np.median(np.abs(Y-X-np.median(Y-X))))
out2['heavier_hotter']=int((Y>0).sum())
# shared-current null with error floors
for fl in (0.01,0.02):
    e=np.sqrt(eT**2+fl**2); z=dT[q]/e[q]
    out2['null_chi2_floor%s'%fl]=float((z**2).sum()); out2['null_gt3_floor%s'%fl]=int((np.abs(z)>3).sum())
    zf=(dT[q]-dTz[q])/e[q]
    out2['zams_chi2_floor%s'%fl]=float((zf**2).sum()); out2['zams_gt3_floor%s'%fl]=int((np.abs(zf)>3).sum())
out2['n_q']=int(q.sum())
# mass-binned M-L slopes with bootstrap errors
rng=np.random.default_rng(145); Mt=np.maximum(Mp,Ms); bins=[(0.1,0.8),(0.8,1.5),(1.5,3.0),(3.0,30.0)]; bl=[]
for lo,hi in bins:
    m=np.where((Mt>=lo)&(Mt<hi)&sel)[0]
    so=(dM[m]*dL[m]).sum()/(dM[m]**2).sum(); sz=(dM[m]*dLz[m]).sum()/(dM[m]**2).sum()
    bs=[]
    for _ in range(4000):
        k=rng.choice(m,len(m)); bs.append((dM[k]*dL[k]).sum()/(dM[k]**2).sum())
    lo16,hi84=np.percentile(bs,[16,84])
    bl.append({'bin':[lo,hi],'n':int(len(m)),'obs':float(so),'lo':float(lo16),'hi':float(hi84),'zams':float(sz)})
out2['bins']=bl
for k,v in out2.items(): print(k, v if not isinstance(v,float) else round(v,4))
D=json.load(open('estars_results.json')); D['teff']=out2
D['figT']=[[round(float(a),4),round(float(b),4)] for a,b in zip(X,Y)]
json.dump(D,open('estars_results.json','w'))

# ---- Part 3: the shared-current null slope per mass bin, from measured radii (L proportional to R^2) ----
for b in out2['bins']:
    lo,hi=b['bin']; m=np.where((Mt>=lo)&(Mt<hi)&sel)[0]
    b['null']=float((dM[m]*2*dR[m]).sum()/(dM[m]**2).sum())
D=json.load(open('estars_results.json')); D['teff']['bins']=out2['bins']
FT=np.array(D['figT']); print('figT x range',FT[:,0].min(),FT[:,0].max(),'y range',FT[:,1].min(),FT[:,1].max())
print('null slopes',[round(b['null'],2) for b in out2['bins']])
json.dump(D,open('estars_results.json','w'))
