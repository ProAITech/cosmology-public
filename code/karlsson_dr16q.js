/* Workbench Entry 76: is there a Karlsson periodicity, Δlog10(1+z) = 0.089 with peaks at z = 0.06, 0.30, 0.60, 0.96, 1.41, 1.96,
   in the redshifts of 749,742 SDSS DR16 quasars (Lyke et al. 2020, VizieR VII/289/dr16q, column z, QSO = 1)?
   Run in a browser tab on a VizieR mirror:  await KQ.load();  KQ.run();
   Method: histogram of x = log10(1+z) in bins of 0.001; the smooth part is a running mean over exactly one trial period, which removes a
   periodic signal of that period from the trend and leaves it in the residual; the fractional modulation m at that period is
   m_c = 2 Σ resid cos(2π(x - x0)/P) / Σ trend (and m_s with sine). Karlsson's phase is x0 = log10(1.06).
   Survey selection puts its own structure into the redshift distribution, so the modulation at P = 0.089 is compared with the
   modulation at every other trial period between 0.03 and 0.20.                                                                     */
window.KQ=(function(){
const K={P:0.089,x0:Math.log10(1.0605),bw:0.001};
K.load=async function(){const ft=async u=>{for(let k=0;k<3;k++){try{const r=await fetch(u);if(r.ok)return await r.text();}catch(e){}}return 'FAIL';};
 const t=await ft('/viz-bin/asu-tsv?-source=VII/289/dr16q&-out=z,r_z,QSO&-out.max=unlimited&-oc.form=dec');
 const rows=t.split('\n').filter(l=>l.length&&!l.startsWith('#')).slice(3).map(l=>l.split('\t'));
 const Z=[],S=[];for(const r of rows){if(r.length<3||r[2].trim()!=='1')continue;const z=+r[0];if(!(z>0))continue;Z.push(z);S.push(r[1].trim());}
 K.z=Float64Array.from(Z);K.src=S;return Z.length;};
K.hist=function(z,xlo,xhi){const n=Math.round((xhi-xlo)/K.bw),h=new Float64Array(n);for(const v of z){const x=Math.log10(1+v);if(x<xlo||x>=xhi)continue;h[Math.floor((x-xlo)/K.bw)]++;}return h;};
// modulation at period P and reference phase x0, from a histogram h starting at xlo
K.mod=function(h,xlo,P,x0){const w=Math.round(P/K.bw),half=Math.floor(w/2),n=h.length;let sc=0,ss=0,st=0;
 const cs=new Float64Array(n+1);for(let i=0;i<n;i++)cs[i+1]=cs[i]+h[i];
 for(let i=half;i<n-(w-half);i++){const tr=(cs[i-half+w]-cs[i-half])/w,r=h[i]-tr,x=xlo+(i+0.5)*K.bw,ph=2*Math.PI*(x-x0)/P;sc+=r*Math.cos(ph);ss+=r*Math.sin(ph);st+=tr;}
 const mc=2*sc/st,ms=2*ss/st;return {mc,ms,m:Math.hypot(mc,ms),phase:Math.atan2(ms,mc),N:st,sig:Math.sqrt(2/st)};};
K.spectrum=function(h,xlo){const out=[];for(let P=0.030;P<=0.2001;P+=0.001){const r=K.mod(h,xlo,P,K.x0);out.push([+P.toFixed(3),r.m]);}return out;};
K.one=function(z,xlo,xhi){const h=K.hist(z,xlo,xhi),r=K.mod(h,xlo,K.P,K.x0),sp=K.spectrum(h,xlo),others=sp.filter(s=>Math.abs(s[0]-K.P)>0.006).map(s=>s[1]).sort((a,b)=>a-b);
 const rank=others.filter(v=>v>=r.m).length/others.length;return {n:Math.round(r.N),mc:r.mc,ms:r.ms,m:r.m,poisson:r.sig,medOther:others[others.length>>1],p84:others[Math.floor(0.84*others.length)],p95:others[Math.floor(0.95*others.length)],fracOtherAbove:rank,spectrum:sp};};
// the peaks themselves: counts in windows ±0.01 in log(1+z) around each Karlsson peak, against the running-mean expectation
K.peaks=function(z){const xlo=0,xhi=0.75,h=K.hist(z,xlo,xhi),w=Math.round(K.P/K.bw),half=Math.floor(w/2),cs=new Float64Array(h.length+1);for(let i=0;i<h.length;i++)cs[i+1]=cs[i]+h[i];
 return [0.061,0.30,0.60,0.96,1.41,1.96].map(zp=>{const xc=Math.log10(1+zp),i0=Math.floor((xc-0.01-xlo)/K.bw),i1=Math.floor((xc+0.01-xlo)/K.bw);let o=0,e=0;
  for(let i=i0;i<=i1;i++){o+=h[i];const a=Math.max(0,i-half),b=Math.min(h.length,a+w);e+=(cs[b]-cs[a])/(b-a);}return {z:zp,obs:o,exp:Math.round(e),ratio:o/e};});};
K.run=function(){const z=K.z,vi=z.filter((v,i)=>K.src[i]==='VI'),pipe=z.filter((v,i)=>K.src[i]==='PIPE');
 const sets={all_z005_22:[z,Math.log10(1.05),Math.log10(3.2)],all_full:[z,0.0,0.72],vi_z005_22:[vi,Math.log10(1.05),Math.log10(3.2)],pipe_z005_22:[pipe,Math.log10(1.05),Math.log10(3.2)]};
 const out={};for(const [k,[zz,a,b]] of Object.entries(sets)){const r=K.one(zz,a,b);delete r.spectrum;out[k]=r;}
 out.peaks=K.peaks(z);out.spectrum_all=K.one(z,Math.log10(1.05),Math.log10(3.2)).spectrum;return out;};
return K;})();
/* Pairs: quasars around bright nearby galaxies (RC3, VizieR VII/155/rc3: BT <= 13.5, 300 < cz < 7000 km/s), redshift taken in the galaxy's
   frame, (1+zQ)/(1+zG) - 1, as in the original claim.  Null: the same galaxies moved in right ascension by 15 to 334 degrees.
   KQ.pairs(Q, G, mode, raShift) with Q = [[ra, dec, z], ...], G = [{ra, de, cz}, ...]; mode 'arcmin' = within 30', 'kpc' = within 200 kpc. */
window.KQ.pairs=function(Q,G,mode,raShift){const KQ=window.KQ,d2r=Math.PI/180;if(!KQ._grid){KQ._grid=new Map();Q.forEach((q,i)=>{const k=Math.floor(q[0])+','+Math.floor(q[1]);if(!KQ._grid.has(k))KQ._grid.set(k,[]);KQ._grid.get(k).push(i);});}
 const sep=(a1,d1,a2,d2)=>Math.acos(Math.min(1,Math.sin(d1*d2r)*Math.sin(d2*d2r)+Math.cos(d1*d2r)*Math.cos(d2*d2r)*Math.cos((a1-a2)*d2r)))/d2r,out=[];
 for(const g of G){const ra=(g.ra+raShift+360)%360,D=g.cz/70,rad=mode==='kpc'?Math.min(2.0,(0.2/D)/d2r):0.5,cd=Math.max(0.05,Math.cos(g.de*d2r)),rr=Math.ceil(rad/cd)+1;
  for(let b=Math.floor(g.de-rad);b<=Math.floor(g.de+rad);b++)for(let a=Math.floor(ra)-rr;a<=Math.floor(ra)+rr;a++){const L=KQ._grid.get(((a%360+360)%360)+','+b);if(!L)continue;for(const i of L){const q=Q[i],s=sep(ra,g.de,q[0],q[1]);if(s<=rad)out.push([q[2],g.cz/299792.458,s]);}}}return out;};
window.KQ.pairStat=function(P,rel){const KQ=window.KQ,zz=Float64Array.from(P.map(p=>rel?(1+p[0])/(1+p[1])-1:p[0])),lo=Math.log10(1.05),h=KQ.hist(zz,lo,Math.log10(3.2)),r=KQ.mod(h,lo,KQ.P,KQ.x0);return {n:P.length,mc:r.mc,ms:r.ms,sig:r.sig};};
