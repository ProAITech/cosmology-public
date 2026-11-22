/* Workbench Entry 84: is the galaxy distribution a fractal of dimension 2?
   Run in a browser tab on https://skyserver.sdss.org/dr18/SearchTools/sql (same-origin SqlSearch), 11 November 2026.
   Claim tested: matter density falls inversely with the scale on which it is measured, "a fractal with dimension 2",
   up to gigaparsec scales (E. J. Lerner, Asia Times, 11 December 2023).
   Measure: correlation dimension D2(r) = d ln N(<r) / d ln r of counts in spheres around galaxies. A fractal of
   dimension 2 gives D2 = 2 at every r; an even distribution gives 3.                                              */
const SQL = "SELECT TOP 400000 p.ra, p.dec, s.z, p.petroMag_r - p.extinction_r AS r FROM SpecObj s JOIN PhotoObj p ON s.bestObjID = p.objID " +
  "WHERE s.class='GALAXY' AND s.zWarning=0 AND s.z BETWEEN 0.02 AND 0.11 AND p.petroMag_r - p.extinction_r < 17.77 " +
  "AND p.ra BETWEEN 120 AND 250 AND p.dec BETWEEN -3 AND 68";        // 318,207 rows; 318,200 after removing repeats
// Volume-limited sample: 0.02 < z < 0.105, M_r < -20.75 with M = r - 5 log10(dL/10 pc) - 1.0 z (rough K-correction);
// distances either flat LCDM (Om 0.3, H0 70) or Lerner's static law d = cz/H0 with dL = d sqrt(1+z).
// Footprint: 0.5 deg cells (RA x equal-area in sin dec) holding at least two galaxies of the flux-limited sample: 2.072 sr.
// Randoms: 4 x the sample, uniform in the footprint and uniform in volume between the distance limits.
/* Estimator 1 (random-normalised, as WiggleZ): around 3000 galaxies, N(<r) = mean over centres of
   [data within r / randoms within r] x (randoms / data); D2 = 3 + d ln N / d ln r between successive radii.
   Estimator 2 (full spheres, no randoms, no homogeneity assumed): only centres whose whole sphere of radius Rmax lies
   inside the footprint and the distance range (tested at 400 points on and inside the sphere); D2 = d ln <n_data(<r)> / d ln r.
   Rmax 80 Mpc: 800 centres; Rmax 110 Mpc: 184 centres. Errors: jackknife over ten slices of centre right ascension.     */
// ---- the functions as run ----
// GAL = rows [ra, dec, z, r] from SQL above
window.FR = (function(){
  const c=299792.458, H0=70, d2r=Math.PI/180;
  const zt=[],dt=[]; {let s=0,zp=0; for(let i=0;i<=1200;i++){const z=i*0.0001; if(i>0){const E=a=>Math.sqrt(0.3*Math.pow(1+a,3)+0.7); s+=0.0001*0.5*(1/E(zp)+1/E(z));} zt.push(z); dt.push(s*c/H0); zp=z;}}
  const dLCDM=z=>{const i=Math.min(1199,Math.floor(z/0.0001)),f=(z-zt[i])/0.0001;return dt[i]+f*(dt[i+1]-dt[i]);};
  const dLerner=z=>c*z/H0;
  const laws={LCDM:{d:dLCDM, dl:z=>(1+z)*dLCDM(z)}, Lerner:{d:dLerner, dl:z=>Math.sqrt(1+z)*dLerner(z)}};
  const seen=new Set(), gal=[]; for(const g of GAL){const k=g[0].toFixed(5)+':'+g[1].toFixed(5); if(seen.has(k)) continue; seen.add(k); gal.push(g);}
  const cra=0.5, cs=0.5*d2r; const key=(ra,dec)=>Math.floor(ra/cra)+':'+Math.floor(Math.sin(dec*d2r)/cs);
  const occ=new Map(); for(const g of gal){const k=key(g[0],g[1]); occ.set(k,(occ.get(k)||0)+1);}
  const cells=[...occ.entries()].filter(([k,n])=>n>=2).map(([k])=>k.split(':').map(Number)); const cellset=new Set(cells.map(c=>c[0]+':'+c[1]));
  const area=cells.length*(cra*d2r)*cs;
  const xyz=(ra,dec,d)=>[d*Math.cos(dec*d2r)*Math.cos(ra*d2r), d*Math.cos(dec*d2r)*Math.sin(ra*d2r), d*Math.sin(dec*d2r)];
  let seed=12345; const rand=()=>{seed=(seed*16807)%2147483647;return seed/2147483647;};
  const sample=(law,Mlim,zmin,zmax,kz)=>{const L=laws[law]; const pts=[];
    for(const g of gal){ if(!cellset.has(key(g[0],g[1]))) continue; if(g[2]<zmin||g[2]>zmax) continue; const M=g[3]-5*Math.log10(L.dl(g[2])*1e5)-kz*g[2]; if(M>Mlim) continue; pts.push(xyz(g[0],g[1],L.d(g[2])).concat([g[0]]));}
    return {pts,dmin:L.d(zmin),dmax:L.d(zmax)};};
  const randoms=(n,dmin,dmax)=>{const out=[]; while(out.length<n){const c=cells[Math.floor(rand()*cells.length)]; const ra=(c[0]+rand())*cra; const sd=(c[1]+rand())*cs; if(Math.abs(sd)>1) continue; const dec=Math.asin(sd)/d2r; const u=rand(); const d=Math.cbrt(dmin**3+u*(dmax**3-dmin**3)); out.push(xyz(ra,dec,d).concat([ra]));} return out;};
  const inside=(p,Rmax,dmin,dmax)=>{const d=Math.hypot(p[0],p[1],p[2]); if(d-Rmax<dmin||d+Rmax>dmax) return false;
    for(let i=0;i<400;i++){const u=2*((i*0.618034)%1)-1,phi=i*2.399963,s=i<300?1:Math.cbrt((i-300)/100);
      const x=p[0]+Rmax*s*Math.sqrt(1-u*u)*Math.cos(phi),y=p[1]+Rmax*s*Math.sqrt(1-u*u)*Math.sin(phi),z=p[2]+Rmax*s*u;
      const r=Math.hypot(x,y,z); if(r<dmin||r>dmax) return false; const dec=Math.asin(z/r)/d2r; let ra=Math.atan2(y,x)/d2r; if(ra<0) ra+=360;
      if(!cellset.has(Math.floor(ra/cra)+':'+Math.floor(Math.sin(dec*d2r)/cs))) return false;} return true;};
  return {laws,gal,cells,area,sample,randoms,inside};})();
// estimator 1: data and random counts in spheres around data centres (chunked so the page stays responsive)
FR.fastCounts = function (D, Rn, R, ncen, seed0) {
  const cell=20, rmax=R[R.length-1], nc=Math.ceil(rmax/cell);
  const pack=arr=>{const a=new Float64Array(arr.length*3); arr.forEach((p,i)=>{a[3*i]=p[0];a[3*i+1]=p[1];a[3*i+2]=p[2];}); return a;};
  const PD=pack(D), PR=pack(Rn);
  const grid=(P,n)=>{const G=new Map(); for(let i=0;i<n;i++){const k=Math.floor(P[3*i]/cell)*1000003+Math.floor(P[3*i+1]/cell)*1009+Math.floor(P[3*i+2]/cell); let L=G.get(k); if(!L){L=[];G.set(k,L);} L.push(i);} for(const [k,L] of G) G.set(k,Int32Array.from(L)); return G;};
  const GD=grid(PD,D.length), GR=grid(PR,Rn.length);
  let seed=seed0||777; const rnd=()=>{seed=(seed*16807)%2147483647; return seed/2147483647;};
  const cen=[]; const used=new Set(); while(cen.length<Math.min(ncen,D.length)){const i=Math.floor(rnd()*D.length); if(!used.has(i)){used.add(i); cen.push(i);}}
  const R2=R.map(r=>r*r), nR=R.length, job={done:0,total:cen.length,out:[]};
  const one=ci=>{const x=PD[3*ci],y=PD[3*ci+1],z=PD[3*ci+2]; const nd=new Float64Array(nR), nr=new Float64Array(nR); const kx=Math.floor(x/cell),ky=Math.floor(y/cell),kz=Math.floor(z/cell);
    for(let a=-nc;a<=nc;a++)for(let b=-nc;b<=nc;b++)for(let e=-nc;e<=nc;e++){const k=(kx+a)*1000003+(ky+b)*1009+(kz+e);
      const Ld=GD.get(k); if(Ld) for(let t=0;t<Ld.length;t++){const j=Ld[t]; if(j===ci) continue; const dx=PD[3*j]-x,dy=PD[3*j+1]-y,dz=PD[3*j+2]-z,d2=dx*dx+dy*dy+dz*dz; if(d2>R2[nR-1]) continue; for(let q=nR-1;q>=0&&d2<=R2[q];q--) nd[q]++;}
      const Lr=GR.get(k); if(Lr) for(let t=0;t<Lr.length;t++){const j=Lr[t]; const dx=PR[3*j]-x,dy=PR[3*j+1]-y,dz=PR[3*j+2]-z,d2=dx*dx+dy*dy+dz*dz; if(d2>R2[nR-1]) continue; for(let q=nR-1;q>=0&&d2<=R2[q];q--) nr[q]++;}}
    return {ra:D[ci][3], nd:Array.from(nd), nr:Array.from(nr)};};
  const step=()=>{const t0=performance.now(); while(job.done<job.total&&performance.now()-t0<200){job.out.push(one(cen[job.done])); job.done++;} if(job.done<job.total) setTimeout(step,0);};
  setTimeout(step,0); return job;};
FR.d2 = function (res, f, R) {
  const N=list=>R.map((_,k)=>{let s=0,n=0; for(const c of list){ if(c.nr[k]>0){ s+=c.nd[k]/c.nr[k]*f; n++; } } return s/n;});
  const D2=Nr=>R.slice(0,-1).map((r,k)=>[+Math.sqrt(r*R[k+1]).toFixed(1), 3+Math.log(Nr[k+1]/Nr[k])/Math.log(R[k+1]/r)]);
  const all=N(res), d=D2(all);
  const ras=res.map(c=>c.ra).sort((a,b)=>a-b), q=k=>ras[Math.floor(k*ras.length/10)];
  const jk=[]; for(let k=0;k<10;k++){const lo=q(k),hi=k<9?q(k+1):1e9; jk.push(D2(N(res.filter(c=>c.ra<lo||c.ra>=hi))));}
  const err=d.map((_,i)=>{const m=jk.reduce((a,j)=>a+j[i][1],0)/10; return Math.sqrt(9/10*jk.reduce((a,j)=>a+(j[i][1]-m)**2,0));});
  const errN=R.map((_,i)=>{const vals=[]; for(let k=0;k<10;k++){const lo=q(k),hi=k<9?q(k+1):1e9; vals.push(N(res.filter(c=>c.ra<lo||c.ra>=hi))[i]);} const m=vals.reduce((a,b)=>a+b,0)/10; return Math.sqrt(9/10*vals.reduce((a,b)=>a+(b-m)**2,0));});
  return {N:all.map((v,i)=>[R[i],+v.toFixed(4),+errN[i].toFixed(4)]), D2:d.map((x,i)=>[x[0],+x[1].toFixed(3),+err[i].toFixed(3)])};};
// estimator 2: full spheres, data only
FR.fullSphereJob = function (D, R, dmin, dmax, ncenMax, seed0) {
  const Rmax=R[R.length-1]; let seed=seed0||99; const rnd=()=>{seed=(seed*16807)%2147483647; return seed/2147483647;};
  const job={phase:'select',cand:[],tries:0,done:0,per:[],R}; const cell=20, nc=Math.ceil(Rmax/cell); let G=null; const R2=R.map(r=>r*r);
  const step=()=>{const t0=performance.now();
    if(job.phase==='select'){ while(job.tries<60000&&job.cand.length<ncenMax&&performance.now()-t0<200){job.tries++; const i=Math.floor(rnd()*D.length); if(FR.inside(D[i],Rmax,dmin,dmax)) job.cand.push(i);}
      if(job.tries>=60000||job.cand.length>=ncenMax) job.phase='grid'; return setTimeout(step,0);}
    if(job.phase==='grid'){ G=new Map(); D.forEach((p,i)=>{const k=Math.floor(p[0]/cell)*1000003+Math.floor(p[1]/cell)*1009+Math.floor(p[2]/cell); let L=G.get(k); if(!L){L=[];G.set(k,L);} L.push(i);}); job.phase='count'; return setTimeout(step,0);}
    if(job.phase==='count'){ while(job.done<job.cand.length&&performance.now()-t0<200){const ci=job.cand[job.done], p=D[ci], nd=new Float64Array(R.length), kx=Math.floor(p[0]/cell), ky=Math.floor(p[1]/cell), kz=Math.floor(p[2]/cell);
        for(let a=-nc;a<=nc;a++)for(let b=-nc;b<=nc;b++)for(let e=-nc;e<=nc;e++){const L=G.get((kx+a)*1000003+(ky+b)*1009+(kz+e)); if(!L) continue;
          for(let t=0;t<L.length;t++){const j=L[t]; if(j===ci) continue; const q=D[j], d2=(q[0]-p[0])**2+(q[1]-p[1])**2+(q[2]-p[2])**2; if(d2>R2[R.length-1]) continue; for(let k=R.length-1;k>=0&&d2<=R2[k];k--) nd[k]++;}}
        job.per.push({ra:p[3], nd:Array.from(nd)}); job.done++;}
      if(job.done<job.cand.length) return setTimeout(step,0); job.phase='done';}};
  setTimeout(step,0); return job;};
FR.fsResult = function (job) { const R=job.R, per=job.per, n=per.length;
  const meanOf=sub=>R.map((_,i)=>sub.reduce((a,c)=>a+c.nd[i],0)/sub.length);
  const slopes=m=>R.slice(0,-1).map((r,i)=>Math.log(m[i+1]/m[i])/Math.log(R[i+1]/r));
  const mean=meanOf(per), D2=slopes(mean);
  const ras=per.map(c=>c.ra).sort((a,b)=>a-b), q=k=>ras[Math.floor(k*ras.length/10)], jk=[];
  for(let k=0;k<10;k++){const lo=q(k),hi=k<9?q(k+1):1e9; jk.push(slopes(meanOf(per.filter(c=>c.ra<lo||c.ra>=hi))));}
  const err=D2.map((_,i)=>{const mm=jk.reduce((a,j)=>a+j[i],0)/10; return Math.sqrt(0.9*jk.reduce((a,j)=>a+(j[i]-mm)**2,0));});
  return {centres:n, tries:job.tries, mean_counts:mean.map(v=>+v.toFixed(1)), D2:D2.map((v,i)=>[+Math.sqrt(R[i]*R[i+1]).toFixed(1), +v.toFixed(3), +err[i].toFixed(3)])};};
/* Runs: S = FR.sample('LCDM' or 'Lerner', -20.75, 0.02, 0.105, 1.0); Rn = FR.randoms(4*S.pts.length, S.dmin, S.dmax);
   FR.d2(FR.fastCounts(S.pts, Rn, [10,14,20,28,40,50,60,70,80,100,120], 3000, 777 or 778).out, 4, R);
   FR.fsResult(FR.fullSphereJob(S.pts, [10,14,20,28,40,50,60,70,80], S.dmin, S.dmax, 800, 31)) and with
   [10,14,20,28,40,50,60,70,80,95,110], 500 centres sought (184 found), seed 57.                                       */
