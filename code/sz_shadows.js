/* Workbench Entry 22, Part 7: do the shadows that galaxy clusters cast on the CMB fade with distance?
   A cluster's hot gas leaves a cool spot in the CMB (the Sunyaev-Zeldovich effect) at 98-150 GHz. Catalogues turn the depth of the spot into a mass.
   The same clusters have X-ray masses, which a radio/microwave fog cannot touch. A fog of optical depth tau(z) cuts the spot by exp(-tau),
   so the SZ mass (signal ~ mass^5/3) comes out low by exp(-0.6 tau). Test: log(M_SZ / M_X) against redshift.
   Runs in a browser tab on heasarc.gsfc.nasa.gov; window.tapx is the query helper defined in radio_fog.js. */
window.fire('cl',async()=>{ const D=window.D;
 D.mcxc=(await window.tapx("select name, ra, dec, redshift, lx_500, mass_500, catalog_name, sub_catalog_name from mcxc",5000)).rows;                 // 1,743 X-ray clusters
 D.act=(await window.tapx("select name, ra, dec, snr, y_c, redshift, mass_500c, mass_500c_uncorr, redshift_type from actszclust",9000)).rows;        // ACT DR5, 4,195
 D.psz=(await window.tapx("select name, ra, dec, snr, redshift, mass_sz, y5r500, mcxc_name, cosmology_sample_flag, validation_status from plancksz2",5000)).rows;   // Planck PSZ2, 1,653
 D.spt=(await window.tapx("select name, ra, dec, best_significance, redshift, mass_500c, int_sz_compton from sptszgalcl",3000)).rows;                // SPT-SZ, 677
 D.xcs=(await window.tapx("select name, ra, dec, redshift, kt, kt_neg_err, kt_pos_err, counts from xcs",3000)).rows; return 1; });                    // XMM Cluster Survey temperatures, 503
window.SZ=(function(){ const R=Math.PI/180; const D=window.D;
 const med=a=>{const s=a.slice().sort((x,y)=>x-y);const n=s.length;return n?(n%2?s[(n-1)/2]:(s[n/2-1]+s[n/2])/2):NaN;}; const mean=a=>a.reduce((x,y)=>x+y,0)/a.length; const sd=a=>{const m=mean(a);return Math.sqrt(a.reduce((x,y)=>x+(y-m)*(y-m),0)/(a.length-1));};
 const sep=(a1,d1,a2,d2)=>{const x=Math.sin(d1*R)*Math.sin(d2*R)+Math.cos(d1*R)*Math.cos(d2*R)*Math.cos((a1-a2)*R);return Math.acos(Math.min(1,x))/R*60;};   // arcmin
 // nearest X-ray cluster within rad arcmin and |dz| < 0.03(1+z)
 function match(A,X,rad){ const out=[]; for(const a of A){ if(!(a.z>0)||!(a.M>0))continue; let best=null,bs=rad; for(const x of X){ if(Math.abs(x.dec-a.dec)>0.2||!(x.Mx>0)||Math.abs(x.z-a.z)>0.03*(1+a.z))continue; const s=sep(a.ra,a.dec,x.ra,x.dec); if(s<bs){bs=s;best=x;} }
   if(best)out.push({z:a.z,M:a.M,Mx:best.Mx,r:Math.log10(a.M/best.Mx),snr:a.snr,name:a.name,xn:best.name,sub:best.sub}); } return out; }
 function bins(P,zb){ return zb.slice(0,-1).map((a,i)=>{const B=P.filter(p=>p.z>=a&&p.z<zb[i+1]); const v=B.map(p=>p.r); return [med(B.map(p=>p.z)),B.length,med(v),B.length>2?1.2533*sd(v)/Math.sqrt(B.length):null];}); }
 const E=z=>Math.sqrt(0.3*(1+z)**3+0.7);
 const mcxc=()=>D.mcxc.map(r=>({name:r[0],ra:r[1],dec:r[2],z:r[3],Mx:r[5]/1e14,sub:r[7]}));                                             // mass from X-ray luminosity
 const xcs=()=>D.xcs.filter(r=>r[4]>0&&r[3]>0).map(r=>({name:r[0],ra:r[1],dec:r[2],z:r[3],Mx:3.84*Math.pow(r[4]/5,1.71)/E(r[3]),sub:'XCS'}));   // mass from X-ray temperature (Arnaud et al. 2005)
 const act=()=>D.act.map(r=>({name:r[0],ra:r[1],dec:r[2],snr:r[3],z:r[5],M:r[6],Mu:r[7]}));
 const psz=()=>D.psz.filter(r=>r[4]>0&&r[5]>0).map(r=>({name:r[0],ra:r[1],dec:r[2],snr:r[3],z:r[4],M:r[5]}));
 // Lerner's stated fog: radio ~ D^-0.395 out to 600 Mpc (8% arrives), then optical depth rising linearly to 7.6 at z = 2.2
 const tau=z=>{const d=4282.7*z; return z<=0.14?0.395*Math.log(d):2.526+(7.6-2.526)*(z-0.14)/(2.2-0.14);};
 // fog strength s: r = a + s * (-0.6 tau(z) / ln 10); s = 1 is the stated fog, s = 0 is none
 function fit(P){ const x=P.map(p=>-0.6*tau(p.z)/Math.LN10), y=P.map(p=>p.r); const mx=mean(x),my=mean(y); let sxx=0,sxy=0; for(let i=0;i<x.length;i++){sxx+=(x[i]-mx)**2;sxy+=(x[i]-mx)*(y[i]-my);}
   const s=sxy/sxx; let rss=0; for(let i=0;i<x.length;i++)rss+=(y[i]-my-s*(x[i]-mx))**2; return {n:P.length,s,se:Math.sqrt(rss/(x.length-2)/sxx),scatter:Math.sqrt(rss/(x.length-2))}; }
 return {med,mean,sd,sep,match,bins,mcxc,xcs,act,psz,tau,fit,E}; })();
/* run:  X=SZ.mcxc();  am=SZ.match(SZ.act(),X,2.5).filter(p=>p.z>=0.03);  pm=SZ.match(SZ.psz(),X,5).filter(p=>p.z>=0.03);  ax=SZ.match(SZ.act(),SZ.xcs(),2.5);
         SZ.fit(am); SZ.fit(pm); SZ.fit(pm.filter(p=>/REFLEX|NORAS|BCS/.test(p.sub))); SZ.fit(ax); SZ.bins(am,[0.03,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.9]) */

/* ---------------- Entry 24: shadow depth against X-ray temperature (no distances) ----------------
   Extra downloads: acceptcat (name, ra, dec, redshift, kt), baxgalclus (..., tx), xmmsdssgcs (..., kt_ap), xcs (..., kt), ascaegclus (..., kt).
   For each ACT cluster the nearest temperature within 2.5 arcmin and |dz| < 0.03(1+z) is taken, in the order ACCEPT, XCS, 2XMMi/SDSS, BAX, ASCA.
   y0 = E(z)^2 x (mass_500c_uncorr)^1.08 undoes ACT's conversion from signal to mass and returns the measured central Compton parameter
   (arbitrary units). The 114 matched clusters are in clusters_y_T.txt; the fits are in cluster_density_test.py. */
window.YT_match=function(A,X,rad){ const R=Math.PI/180, E=z=>Math.sqrt(0.3*(1+z)**3+0.7), out=[];
  const sep=(a1,d1,a2,d2)=>{const x=Math.sin(d1*R)*Math.sin(d2*R)+Math.cos(d1*R)*Math.cos(d2*R)*Math.cos((a1-a2)*R);return Math.acos(Math.min(1,x))/R*60;};
  for(const a of A){ let best=null,bs=rad; for(const x of X){ if(Math.abs(x.dec-a.dec)>0.2||Math.abs(x.z-a.z)>0.03*(1+a.z))continue; const s=sep(a.ra,a.dec,x.ra,x.dec); if(s<bs){bs=s;best=x;} }
    if(best){const e=E(a.z); out.push({name:a.name,z:a.z,T:best.T,snr:a.snr,logy0:Math.log10(e*e*Math.pow(a.Mu,1.08))});} } return out; };
