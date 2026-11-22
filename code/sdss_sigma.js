// Entry 87: galaxy velocity dispersions of X-ray clusters from SDSS DR18 spectra.
// Run in the browser on skyserver.sdss.org (same-origin SQL search). Cluster list:
// Migkas et al. 2020 (VizieR J/A+A/636/A15), z 0.015-0.15, Dec > -12, b > 20 deg (84 clusters).
// Output per cluster: name, Ngal fetched, r200 used (Mpc), then for each of three selections
// [N members, z_c, sigma, 16th and 84th percentile of 400 bootstraps]:
//   (a) R < r200 (= 1.5 r500, r500 = 1.104 Mpc (T/5 keV)^0.57 / E(z), capped 2.5 Mpc), 1000 km/s gapper
//   (b) R < 1 Mpc, 1000 km/s gapper            (c) as (a) with a 600 km/s gapper (N, sigma only)
const CL={c:299792.458};
CL.E=z=>Math.sqrt(0.3*(1+z)**3+0.7);
CL.DA=z=>{let n=400,s=0;for(let i=0;i<n;i++){const zz=(i+0.5)*z/n;s+=1/CL.E(zz);}return CL.c/70*s*z/n/(1+z);};
CL.sql=async s=>{for(let i=0;i<4;i++){try{const r=await fetch('/dr18/SkyServerWS/SearchTools/SqlSearch?cmd='+encodeURIComponent(s)+'&format=csv');const t=await r.text();if(r.ok&&!/error/i.test(t.slice(0,300)))return t;}catch(e){} await new Promise(r=>setTimeout(r,1500));}return 'FAIL';};
// biweight location (median start) and scale, c = 9 (Beers, Flynn & Gebhardt 1990)
CL.bw=v=>{const n=v.length;const s=[...v].sort((a,b)=>a-b);const med=s[Math.floor(n/2)];const mad=[...v].map(x=>Math.abs(x-med)).sort((a,b)=>a-b)[Math.floor(n/2)]||1;
 let num=0,den=0;for(const x of v){const u=(x-med)/(9*mad);if(Math.abs(u)<1){num+=(x-med)**2*(1-u*u)**4;den+=(1-u*u)*(1-5*u*u);}}
 return {loc:med,scale:Math.sqrt(n)*Math.sqrt(num)/Math.abs(den)};};
// all primary galaxy spectra within 2.5 Mpc and +/-5000 km/s
CL.fetch=async(ra,dec,z0)=>{const da=CL.DA(z0);const th=2.5/da*180/Math.PI;const dz=5000/CL.c*(1+z0);
 const s=`SELECT ra,dec,z FROM SpecObj WHERE class='GALAXY' AND zWarning=0 AND sciencePrimary=1 AND z BETWEEN ${(z0-dz).toFixed(5)} AND ${(z0+dz).toFixed(5)} AND dec BETWEEN ${(dec-th).toFixed(4)} AND ${(dec+th).toFixed(4)} AND ra BETWEEN ${(ra-th/Math.cos(dec*Math.PI/180)).toFixed(4)} AND ${(ra+th/Math.cos(dec*Math.PI/180)).toFixed(4)}`;
 const t=await CL.sql(s); if(t==='FAIL')return null;
 return t.split('\n').filter(l=>/^-?\d/.test(l)).map(l=>l.split(',').map(Number));};
// members: aperture, |v| < 4000 km/s, fixed velocity gapper around the cluster redshift, then 3-sigma biweight clipping
CL.sel=(g,ra,dec,z0,Rap,gap)=>{const da=CL.DA(z0);const d2r=Math.PI/180;
 const R=g.map(([a,d])=>{const cd=Math.sin(dec*d2r)*Math.sin(d*d2r)+Math.cos(dec*d2r)*Math.cos(d*d2r)*Math.cos((a-ra)*d2r);return Math.acos(Math.min(1,cd))*da;});
 let inn=g.filter((r,i)=>R[i]<0.5&&Math.abs(CL.c*(r[2]-z0)/(1+z0))<3000).map(r=>r[2]).sort((a,b)=>a-b);
 let zc=inn.length>=3?inn[Math.floor(inn.length/2)]:z0;
 let idx=g.map((r,i)=>i).filter(i=>R[i]<Rap&&Math.abs(CL.c*(g[i][2]-zc)/(1+zc))<4000);
 if(idx.length<10)return {n:idx.length};
 let v=idx.map(i=>[i,CL.c*(g[i][2]-zc)/(1+zc)]).sort((a,b)=>a[1]-b[1]);
 let k0=0,best=1e9; v.forEach((x,k)=>{if(Math.abs(x[1])<best){best=Math.abs(x[1]);k0=k;}});
 let lo=k0,hi=k0; while(lo>0&&v[lo][1]-v[lo-1][1]<gap)lo--; while(hi<v.length-1&&v[hi+1][1]-v[hi][1]<gap)hi++;
 let mem=v.slice(lo,hi+1).map(x=>x[0]); let sig=0;
 for(let it=0;it<30;it++){const vv=mem.map(i=>CL.c*(g[i][2]-zc)/(1+zc)); if(vv.length<10)return {n:vv.length};
  const b=CL.bw(vv); zc=zc+b.loc*(1+zc)/CL.c; sig=b.scale;
  const nm=mem.filter(i=>Math.abs(CL.c*(g[i][2]-zc)/(1+zc))<3*sig); if(nm.length===mem.length)break; mem=nm;}
 const vv=mem.map(i=>CL.c*(g[i][2]-zc)/(1+zc)); let bs=[];for(let k=0;k<400;k++){const s=vv.map(()=>vv[Math.floor(Math.random()*vv.length)]);bs.push(CL.bw(s).scale);} bs.sort((a,b)=>a-b);
 return {n:mem.length,zc:+zc.toFixed(5),sig:Math.round(sig),lo:Math.round(bs[64]),hi:Math.round(bs[335])};};
// cand: rows [name, z, T, E_T, e_T, LX, RA, Dec] (migkas_cand.txt)
CL.run=async cand=>{const out=[];for(const r of cand){const ra=+r[6],dec=+r[7],z0=+r[1],T=+r[2];
 const g=await CL.fetch(ra,dec,z0); if(!g||g.length<10)continue;
 const r500=1.104*Math.pow(T/5,0.57)/CL.E(z0), r200=1.5*r500;
 const a=CL.sel(g,ra,dec,z0,Math.min(r200,2.5),1000), b=CL.sel(g,ra,dec,z0,1.0,1000), c=CL.sel(g,ra,dec,z0,Math.min(r200,2.5),600);
 out.push([r[0],g.length,+r200.toFixed(2),a.n,a.zc,a.sig,a.lo,a.hi,b.n,b.sig,b.lo,b.hi,c.n,c.sig]);} return out;};
