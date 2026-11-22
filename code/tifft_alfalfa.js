/* Workbench Entry 77: galaxy redshift periodicity (Tifft's 72 km/s; Guthrie & Napier's 37.6 km/s in galactocentric velocities of nearby
   spirals) tested on ALFALFA HI velocities (Haynes et al. 2018, VizieR J/ApJ/861/49/table2, Vhel < 3000 km/s, HI code 1).
   Run in a browser tab on a VizieR mirror:  await TF.load();  TF.run();
   Statistic: m_c = (2/N) Σ cos(2π V/P) is the fractional modulation with peaks at V = n P (the claimed phase); m = (2/N)|Σ e^{2πiV/P}| at any
   phase. Its Poisson error is √(2/N); the spread of m over trial periods 20-100 km/s shows what chance and clustering give.        */
window.TF=(function(){const d2r=Math.PI/180,T={};
 const aG=192.85948*d2r,dG=27.12825*d2r,lN=122.93192*d2r;
 T.gal=(ra,de)=>{ra*=d2r;de*=d2r;const b=Math.asin(Math.sin(de)*Math.sin(dG)+Math.cos(de)*Math.cos(dG)*Math.cos(ra-aG));const y=Math.cos(de)*Math.sin(ra-aG),x=Math.sin(de)*Math.cos(dG)-Math.cos(de)*Math.sin(dG)*Math.cos(ra-aG);return [lN-Math.atan2(y,x),b];};
 T.vgc=(a,U,V,W)=>{const [l,b]=T.gal(a.ra,a.de);return a.v+U*Math.cos(b)*Math.cos(l)+V*Math.cos(b)*Math.sin(l)+W*Math.sin(b);};
 T.ray=(vs,P)=>{let c=0,s=0;for(const v of vs){const p=2*Math.PI*v/P;c+=Math.cos(p);s+=Math.sin(p);}const n=vs.length;return {mc:2*c/n,ms:2*s/n,m:2*Math.hypot(c,s)/n,sig:Math.sqrt(2/n)};};
 T.load=async function(){const ft=async u=>{for(let k=0;k<5;k++){try{const r=await fetch(u);if(r.ok){const t=await r.text();if(!t.includes('STDB_ERR'))return t;}}catch(e){}await new Promise(r=>setTimeout(r,2000));}return 'FAIL';};
  const t=await ft('/viz-bin/asu-tsv?-source=J/ApJ/861/49/table2&-out.add=_RAJ,_DEJ&-out=Vhel,W50,e_W50,SNR,HI&-out.max=unlimited&-oc.form=dec&Vhel=<3000');
  const rows=t.split('\n').filter(l=>l&&l[0]!=='#'),hdr=rows[0].split('\t').map(s=>s.trim());T.A=[];
  for(const l of rows.slice(3)){const r=l.split('\t'),o={};hdr.forEach((h,i)=>o[h]=r[i]);const v=+o.Vhel;if(!isFinite(v))continue;T.A.push({ra:+o._RAJ2000,de:+o._DEJ2000,v,ew:+o.e_W50,snr:+o.SNR,hi:+o.HI});}return T.A.length;};
 T.run=function(){const A=T.A,sel={all:A.filter(a=>a.hi===1),good:A.filter(a=>a.hi===1&&a.snr>=15&&a.ew<=8),best:A.filter(a=>a.hi===1&&a.snr>=30&&a.ew<=4),near:A.filter(a=>a.hi===1&&a.snr>=15&&a.ew<=8&&a.v<2600)};
  const vecs={modern:[11.1,245.0,7.25],iau:[9,232,7],lsr220:[0,220,0]},out={};
  for(const [sk,S] of Object.entries(sel)){out[sk]={n:S.length};for(const [vk,V] of Object.entries(vecs)){const vs=S.map(a=>T.vgc(a,...V)),r1=T.ray(vs,37.6),r2=T.ray(vs,72.0),r3=T.ray(vs,36.0),sp=[];
   for(let P=20;P<=100;P+=0.2)sp.push(T.ray(vs,P).m);sp.sort((a,b)=>a-b);
   out[sk][vk]={P376:[r1.mc,r1.m],P72:[r2.mc,r2.m],P36:[r3.mc,r3.m],sig:r1.sig,spec_med:sp[sp.length>>1],spec_p95:sp[Math.floor(.95*sp.length)],spec_max:sp[sp.length-1],frac_above_376:sp.filter(x=>x>=r1.m).length/sp.length};}}
  return out;};
 return T;})();
