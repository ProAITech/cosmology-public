/* des_dovekie_fits.js — Workbench Entry 58
   The DES five-year supernovae with the 2025 "Dovekie" calibration (Popovic et al., arXiv:2511.07517), read in a
   browser tab from the collaboration's public table
     github.com/des-science/DES-SN5YR  4_DISTANCES_COVMAT/DES-Dovekie_Metadata.csv   (1,820 rows, SNANA format)
   and fitted there with diagonal errors (MUERR_RENORM). window.SN holds the rows as objects keyed by column name.

   Tests: (1) flat LCDM against the logarithmic law d_L = (c/H0)(1+z) ln(1+z) and the other one-number laws of
   Entry 31, before and after the progenitor-age correction used by Nguyen (arXiv:2608.00806),
   mu -> mu - 0.183 (1 - exp(-2.2 z)); (2) a binned Hubble diagram for the joint fits in age_bias_referee.py;
   (3) the supernovae's own hosts: residual against host colour, mass and star-formation rate at fixed redshift.
   Every law has a free additive constant. Comoving distance from zHD; the (1+z) factor from zHEL, as in the paper. */
window.DV=(function(){
  const ZG=[],N=2400,ZMAX=1.2;for(let i=0;i<=N;i++)ZG.push(i*ZMAX/N);
  function dcTable(om){const t=new Float64Array(N+1);let s=0,prev=1;for(let i=1;i<=N;i++){const z=ZG[i],e=1/Math.sqrt(om*Math.pow(1+z,3)+1-om);s+=0.5*(prev+e)*(ZMAX/N);prev=e;t[i]=s;}return t;}
  function interp(t,z){const x=z/ZMAX*N,i=Math.min(N-1,Math.floor(x));return t[i]+(t[i+1]-t[i])*(x-i);}
  const LAWS={log:z=>Math.log(1+z),lerner:z=>z/Math.sqrt(1+z),cz_stretched:z=>z,tired:z=>Math.log(1+z)/Math.sqrt(1+z),coasting:z=>z*(1+z/2)/(1+z),matter_only:z=>2*(1-1/Math.sqrt(1+z))};  /* comoving-distance part: d_L = (1+z) x this */
  const abc=z=>0.183*(1-Math.exp(-2.2*z));
  function rows(o){o=Object.assign({abc:0,sel:null,err:'MUERR_RENORM'},o||{});return window.SN.filter(r=>!o.sel||o.sel(r)).map(r=>({z:r.zHD,zh:r.zHEL,mu:r.MU-o.abc*abc(r.zHD),w:1/(r[o.err]*r[o.err]),r}));}
  function chi(R,dfn,stepz){/* dfn(z) = comoving distance in c/H0; one free constant, or two with a step at stepz */
    const grp=stepz?[R.filter(p=>p.z<stepz),R.filter(p=>p.z>=stepz)]:[R];let c2=0;const off=[];
    for(const G of grp){let sw=0,sd=0;const d=G.map(p=>p.mu-5*Math.log10((1+p.zh)*dfn(p.z)));G.forEach((p,i)=>{sw+=p.w;sd+=p.w*d[i];});const a=sd/sw;off.push(a);G.forEach((p,i)=>{c2+=p.w*(d[i]-a)**2;});}
    return {chi2:c2,off};}
  function lcdm(R,stepz){const f=om=>{const t=dcTable(om);return chi(R,z=>interp(t,z),stepz).chi2;};let a=0.05,b=0.95;const g=(Math.sqrt(5)-1)/2;let c=b-g*(b-a),d=a+g*(b-a),fc=f(c),fd=f(d);
    for(let i=0;i<60;i++){if(fc<fd){b=d;d=c;fd=fc;c=b-g*(b-a);fc=f(c);}else{a=c;c=d;fc=fd;d=a+g*(b-a);fd=f(d);}}
    const om=(a+b)/2,c0=f(om);const root=(lo,hi)=>{for(let i=0;i<50;i++){const m=(lo+hi)/2;if((f(m)-c0-1)*(f(lo)-c0-1)<=0)hi=m;else lo=m;}return (lo+hi)/2;};
    return {om,chi2:c0,lo:root(0.02,om),hi:root(om,0.98)};}
  function table(o,stepz){const R=rows(o),L=lcdm(R,stepz),out={n:R.length,om:+L.om.toFixed(4),omErr:+((L.hi-L.lo)/2).toFixed(4),lcdm:+L.chi2.toFixed(1)};
    for(const k in LAWS){out[k]=+(chi(R,LAWS[k],stepz).chi2-L.chi2).toFixed(1);}return out;}
  function qfit(o,base){const R=rows(o),f=q=>chi(R,z=>(base==='log'?Math.log(1+z):z)*Math.pow(1+z,q-1)).chi2;let a=-1,b=3;const g=(Math.sqrt(5)-1)/2;let c=b-g*(b-a),d=a+g*(b-a),fc=f(c),fd=f(d);
    for(let i=0;i<60;i++){if(fc<fd){b=d;d=c;fd=fc;c=b-g*(b-a);fc=f(c);}else{a=c;c=d;fc=fd;d=a+g*(b-a);fd=f(d);}}const q=(a+b)/2,c0=f(q);let lo=q,hi=q;while(f(lo)-c0<1)lo-=0.002;while(f(hi)-c0<1)hi+=0.002;return {q:+q.toFixed(3),err:+((hi-lo)/2).toFixed(3),chi2:+c0.toFixed(1)};}
  /* amplitude A of the age correction, fitted together with the model: mu -> mu - A x 0.183 (1 - exp(-2.2 z)) */
  function ampl(o,model,om){const R=rows(o);const t=model==='lcdm'?dcTable(om):null;const dfn=model==='lcdm'?(z=>interp(t,z)):LAWS[model];
    let s1=0,sx=0,sy=0,sxx=0,sxy=0;for(const p of R){const y=p.mu-5*Math.log10((1+p.zh)*dfn(p.z)),x=abc(p.z);s1+=p.w;sx+=p.w*x;sy+=p.w*y;sxx+=p.w*x*x;sxy+=p.w*x*y;}
    const den=s1*sxx-sx*sx;return {A:+((s1*sxy-sx*sy)/den).toFixed(3),err:+Math.sqrt(s1/den).toFixed(3)};}
  /* binned Hubble diagram: residuals from a reference curve averaged in bins, then put back on the curve at the mean redshift */
  function binned(edges,o){const R=rows(o),t=dcTable(0.33),ref=z=>5*Math.log10((1+z)*interp(t,z)),out=[];
    for(let i=0;i<edges.length-1;i++){const G=R.filter(p=>p.z>=edges[i]&&p.z<edges[i+1]);if(!G.length)continue;let sw=0,sz=0,sr=0;for(const p of G){sw+=p.w;sz+=p.w*p.z;sr+=p.w*(p.mu-5*Math.log10((1+p.zh)*interp(t,p.z)));}
      const zb=sz/sw;out.push([+zb.toFixed(5),+(ref(zb)+sr/sw).toFixed(4),+(1/Math.sqrt(sw)).toFixed(4),G.length]);}return out;}
  /* hosts: residual from the best LCDM curve, split at the median of a host property inside each redshift bin */
  function hosts(prop,zedges,o){const R=rows(o),L=lcdm(R),t=dcTable(L.om),c=chi(R,z=>interp(t,z)).off[0];R.forEach(p=>{p.res=p.mu-5*Math.log10((1+p.zh)*interp(t,p.z))-c;p.h=p.r[prop];});
    const out=[];let sw=0,sd=0,sws=0,sds=0;
    for(let i=0;i<zedges.length-1;i++){const G=R.filter(p=>p.z>=zedges[i]&&p.z<zedges[i+1]&&isFinite(p.h)&&p.h>-90);if(G.length<20)continue;const hs=G.map(p=>p.h).sort((a,b)=>a-b),med=hs[hs.length>>1];
      const m=S=>{let w=0,x=0,h=0;for(const p of S){w+=p.w;x+=p.w*p.res;h+=p.w*p.h;}return [x/w,1/Math.sqrt(w),h/w,S.length];};const hi=m(G.filter(p=>p.h>=med)),lo=m(G.filter(p=>p.h<med));
      let s1=0,sx=0,sy=0,sxx=0,sxy=0;for(const p of G){s1+=p.w;sx+=p.w*p.h;sy+=p.w*p.res;sxx+=p.w*p.h*p.h;sxy+=p.w*p.h*p.res;}const den=s1*sxx-sx*sx,sl=(s1*sxy-sx*sy)/den,se=Math.sqrt(s1/den);
      const d=hi[0]-lo[0],e=Math.hypot(hi[1],lo[1]);sw+=1/(e*e);sd+=d/(e*e);sws+=1/(se*se);sds+=sl/(se*se);
      out.push({z:[zedges[i],zedges[i+1]],n:G.length,median:+med.toFixed(3),meanProp:+(sx/s1).toFixed(3),step:+d.toFixed(4),stepErr:+e.toFixed(4),gap:+(hi[2]-lo[2]).toFixed(3),slope:+sl.toFixed(4),slopeErr:+se.toFixed(4)});}
    return {prop,bins:out,step:+(sd/sw).toFixed(4),stepErr:+(1/Math.sqrt(sw)).toFixed(4),slope:+(sds/sws).toFixed(4),slopeErr:+(1/Math.sqrt(sws)).toFixed(4)};}
  return {LAWS,abc,rows,chi,lcdm,table,qfit,ampl,binned,hosts,dcTable,interp};
})();

/* Published run (10-11 November 2026). 1,820 supernovae: 1,623 from the DES camera (z 0.073-1.144) and 197 nearby
   ones from other telescopes (z 0.025-0.087). 101 DES events with PROBCC_BEAMS > 0.5 carry errors about 650 times
   larger and do not count. chi2 of each law minus chi2 of flat LCDM:
                                   LCDM matter share   log law   Lerner z sqrt(1+z)   coasting
     all, as published               0.349 +- 0.011      +208.6      +138.8             +31.9      DV.table()
     all, age correction taken off   0.512 +- 0.014       -12.1       -12.7             +86.0      DV.table({abc:1})
     DES camera only                 0.303 +- 0.016       +98.8       +64.3             +26.8      DV.table({sel:des})
     DES camera only, corrected      0.425 +- 0.020       +12.8        +8.3             +59.7
     free step at z = 0.1            0.305 +- 0.016       +94.7       +61.0
   (the collaboration, with its full error matrix: matter share 0.330 +- 0.015.)
   d_L = ln(1+z)(1+z)^q:  q = 1.130 +- 0.010 (all), 1.134 +- 0.016 (DES camera); after the correction 1.013 and 1.038
   d_L = z (1+z)^q:       q = 0.605 +- 0.010 (all), 0.601 +- 0.016 (DES camera); after the correction, DES camera 0.505
   Correction the logarithmic law needs: 1.10 +- 0.08 of the published one (all), 1.44 +- 0.15 (DES camera).
   Selection corrections undone (MU + biasCor_mu), DES camera only: matter share 0.380, log law +33.8, Lerner +19.4,
   but chi2 rises from 1,451 to 2,714 because the correction also removes scatter; q_log 1.071, q_cz 0.538 +- 0.016.
   All 1,820 with corrections undone: log law +146.9, Lerner +98.7 on a chi2 of 2,992.
   Mean correction: nearby sample -0.059; DES z<0.3 -0.036, 0.3-0.6 -0.064, >0.6 -0.087, >0.8 -0.101 mag.
   Hosts, DES camera, six redshift bins, each split at its median (DV.hosts): redder half minus bluer half
   -0.023 +- 0.008 mag (nearby sample -0.026 +- 0.015; all -0.024 +- 0.007; fixed cut at U-R = 1: -0.030 +- 0.008);
   heavier minus lighter -0.021 +- 0.008; more star-forming minus less +0.032 +- 0.008 (slope 0.0105 +- 0.0028 mag
   per dex). Mean host colour: nearby 1.15; DES z<0.3 0.98, 0.3-0.6 0.83, >0.6 0.89. Mean log sSFR: nearby -10.65;
   DES -9.96 (z 0.07-0.2), -9.57 (0.4-0.5), -10.28 (0.7-0.9). */
