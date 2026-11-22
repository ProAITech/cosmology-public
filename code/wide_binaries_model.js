// Forward model for the Gaia wide-binary test (Plasma Cosmology Workbench, Entry 15).
// Predicts the scaled-velocity distribution per separation bin for Newton or MOND, with
//   (1) an eccentricity distribution f(e) ∝ e^alpha that steepens with separation,
//   (2) hidden companions (extra mass plus a photocentre velocity kick),
//   (3) each pair's own Gaia noise, the catalogue's bound-orbit cut and the same debiasing as the data.
// Run in a browser console on vizier.cfa.harvard.edu after loading the catalogue text into window.__wbT
// (see wide_binaries.js for the query). Results land in window.__wbm.
(()=>{
const d2r=Math.PI/180, G=4*Math.PI*Math.PI, KMS=4.74047, a0=1.2e-10/1.502e-4, SQ2=Math.SQRT2;
const MG=[1.0,2.0,2.6,3.4,4.3,4.6,4.9,5.8,7.0,8.4,9.6,11.6,12.7,14.3], MM=[2.4,1.9,1.6,1.35,1.1,1.0,0.95,0.85,0.70,0.57,0.45,0.23,0.16,0.10];
function mass(Mg){ if(Mg<MG[0]||Mg>MG[MG.length-1]) return NaN; for(let i=0;i<MG.length-1;i++){ if(Mg<=MG[i+1]){const t=(Mg-MG[i])/(MG[i+1]-MG[i]); return Math.exp(Math.log(MM[i])*(1-t)+Math.log(MM[i+1])*t);} } return NaN; }
function basis(ra,de){ra*=d2r;de*=d2r; return [[Math.cos(de)*Math.cos(ra),Math.cos(de)*Math.sin(ra),Math.sin(de)],[-Math.sin(ra),Math.cos(ra),0],[-Math.sin(de)*Math.cos(ra),-Math.sin(de)*Math.sin(ra),Math.cos(de)]];}
function rng(seed){let a=seed>>>0; return ()=>{a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
// ---- per-pair quantities from the catalogue (same steps as wide_binaries.js)
function pairs(text, errCut, absCut){
  const P=[];
  for(const l of text.split('\n')){ if(!l||l[0]==='#'||l[0]==='-') continue; const r=l.split('\t'); if(r.length<26||isNaN(parseFloat(r[0]))) continue; const v=r.map(Number);
    const [ra1,de1,ra2,de2,p1,p2,ep1,ep2,a1,a2,b1,b2,ea1,ea2,eb1,eb2,g1,g2]=v;
    if(Math.abs(p1-p2)>3*Math.hypot(ep1,ep2)) continue;
    const w1=1/ep1**2,w2=1/ep2**2, plx=(p1*w1+p2*w2)/(w1+w2), d=1000/plx; if(d>200) continue;
    const M1=mass(g1+5*Math.log10(plx/100)), M2=mass(g2+5*Math.log10(plx/100)); if(!(M1>0&&M2>0)) continue; const Mt=M1+M2;
    const [n1,A1,D1]=basis(ra1,de1), [n2,A2,D2]=basis(ra2,de2), k=1/plx;
    const v1=[0,1,2].map(i=>k*(a1*A1[i]+b1*D1[i])), v2=[0,1,2].map(i=>k*(a2*A2[i]+b2*D2[i]));
    const nm=[0,1,2].map(i=>(n1[i]+n2[i])/2), nn=Math.hypot(...nm), n=nm.map(x=>x/nn);
    let dv=[0,1,2].map(i=>v2[i]-v1[i]); const dot=dv[0]*n[0]+dv[1]*n[1]+dv[2]*n[2]; dv=dv.map((x,i)=>x-dot*n[i]);
    const th=Math.acos(Math.min(1,n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2])), s=th*206264.806*d, vc=Math.sqrt(G*Mt/s);
    const spm=k*Math.sqrt((ea1**2+ea2**2+eb1**2+eb2**2)/2), sp=(35/KMS)*th/SQ2, st=Math.hypot(spm,sp);
    if(absCut ? st*KMS>absCut : st>errCut*vc) continue;
    const vobs=Math.hypot(...dv), vt=Math.sqrt(Math.max(vobs*vobs-2*st*st,0))/vc;
    P.push({s,vt,vc,sx:st/vc,y:G*Mt/(s*s)/a0,vmax:(2.1/KMS)*Math.sqrt(1000/s)/vc});   // vmax: catalogue's bound-orbit cut, in units of vc
  }
  return P;
}
const EDGES=[200,500,1000,2000,3500,5000,7000,10000,14000,20000,30000], NB=EDGES.length-1;
const binOf=s=>{for(let b=0;b<NB;b++) if(s>=EDGES[b]&&s<EDGES[b+1]) return b; return -1;};
function stats(arrs){ // arrs: array of Float32Array (sorted in place) per bin -> [median,q90,tail>sqrt2]
  return arrs.map(a=>{ a.sort(); const n=a.length; if(!n) return [NaN,NaN,NaN]; const q=p=>{const x=p*(n-1),i=Math.floor(x),f=x-i; return a[i]+(i+1<n?f*(a[i+1]-a[i]):0);};
    let lo=0,hi=n; while(lo<hi){const m=(lo+hi)>>1; if(a[m]>SQ2) hi=m; else lo=m+1;} return [q(0.5),q(0.9),(n-lo)/n]; });
}
// ---- gravity: QUMOND point mass in the Galactic external field (true field 1.8 a0), angle-averaged,
// interpolating between the isolated limit nu(y) and the field-dominated limit nu_e (1 + K_e/3)
const NU={newton:null, rar:y=>1/(1-Math.exp(-Math.sqrt(y))), simple:y=>0.5+Math.sqrt(0.25+1/y), standard:y=>Math.sqrt(0.5+Math.sqrt(0.25+1/(y*y)))};
function yExt(nu){ let a=0.05,b=5; for(let i=0;i<60;i++){const m=(a+b)/2; if(nu(m)*m>1.8) b=m; else a=m;} return (a+b)/2; }
function beta(nu){ if(!nu) return ()=>1; const ye=yExt(nu); return y=>{ const yt=Math.hypot(y,ye), h=1e-4, nt=nu(yt), K=(Math.log(nu(yt*(1+h)))-Math.log(nu(yt*(1-h))))/(2*h); return Math.sqrt(nt*(1+Math.pow(Math.tanh(0.825*ye/y),3.7)*K/3)); }; }
// ---- Newtonian projected v-tilde for one random orbit (a = GM = 1)
function orbitV(e,M,cz,ph){ let E=M; for(let i=0;i<40;i++) E-=(E-e*Math.sin(E)-M)/(1-e*Math.cos(E));
  const q=Math.sqrt(1-e*e), x=Math.cos(E)-e, y=q*Math.sin(E), r=1-e*Math.cos(E), vx=-Math.sin(E)/r, vy=q*Math.cos(E)/r, sz=Math.sqrt(1-cz*cz), nx=sz*Math.cos(ph), ny=sz*Math.sin(ph);
  const dr=x*nx+y*ny, dv=vx*nx+vy*ny; return Math.sqrt(Math.max(vx*vx+vy*vy-dv*dv,0))*Math.pow(Math.max(x*x+y*y-dr*dr,1e-12),0.25); }
// opts: {R, grav:[[name, betaFn]], ALPHA, FC, VK, KAP}; defaults reproduce the four-gravity comparison
function run(tag, errCut, absCut, opts={}){ const R=opts.R||6;
  const out=window.__wbm[tag]={state:'setup'};
  const P=pairs(window.__wbT,errCut,absCut).filter(p=>binOf(p.s)>=0), N=P.length, bins=P.map(p=>binOf(p.s));
  // data vector and errors
  const dat=[...Array(NB)].map(()=>[]); P.forEach((p,i)=>{ if(p.vt<3) dat[bins[i]].push(p.vt); });
  const dArr=dat.map(a=>Float32Array.from(a)), D=stats(dArr), rnd=rng(12345), S=[];
  dArr.forEach((a,b)=>{ const n=a.length, m=[],q=[]; for(let r=0;r<200;r++){ const bs=new Float32Array(n); for(let i=0;i<n;i++) bs[i]=a[Math.floor(rnd()*n)]; const s=stats([bs])[0]; m.push(s[0]); q.push(s[1]); }
    const sd=x=>{const mu=x.reduce((u,v)=>u+v,0)/x.length; return Math.sqrt(x.reduce((u,v)=>u+(v-mu)**2,0)/x.length);}; const p=Math.max(D[b][2],0.004);
    S.push([sd(m),sd(q),Math.sqrt(p*(1-p)/n)]); });
  // common random numbers
  const T=N*R, uE=new Float32Array(T), Mm=new Float32Array(T), cz=new Float32Array(T), ph=new Float32Array(T), uc=new Float32Array(T), zk=new Float32Array(T), pk=new Float32Array(T), um=new Float32Array(T), g1=new Float32Array(T), g2=new Float32Array(T);
  const gauss=()=>Math.sqrt(-2*Math.log(1-rnd()))*Math.cos(2*Math.PI*rnd());
  for(let i=0;i<T;i++){ uE[i]=rnd(); Mm[i]=2*Math.PI*rnd(); cz[i]=2*rnd()-1; ph[i]=2*Math.PI*rnd(); uc[i]=rnd(); zk[i]=gauss(); pk[i]=Math.cos(2*Math.PI*rnd()); um[i]=rnd(); g1[i]=gauss(); g2[i]=gauss(); }
  const ALPHA=opts.ALPHA||[1.0,1.3,1.6,2.0], vN=ALPHA.map(aw=>{ const a=new Float32Array(T); for(let i=0;i<N;i++){ const ls=Math.log10(P[i].s), t=Math.min(1,Math.max(0,(ls-3)/(Math.log10(3000)-3))), al=1+(aw-1)*t;   // thermal below 1,000 AU, alpha_w above 3,000 AU
      for(let r=0;r<R;r++){ const j=i*R+r; a[j]=orbitV(Math.pow(uE[j],1/(al+1)),Mm[j],cz[j],ph[j]); } } return a; });
  const GV=opts.grav||Object.keys(NU).map(g=>[g,beta(NU[g])]), GR=GV.map(g=>g[0]), B=GV.map(g=>Float32Array.from(P,p=>g[1](p.y)));
  const FC=opts.FC||[0,0.05,0.1,0.15,0.2,0.3,0.4], VK=opts.VK||[0.15,0.3,0.6,1.2], KAP=opts.KAP||[0.95,0.96,0.97,0.98,0.99,1.0,1.01,1.02,1.03];
  const grid=[]; GR.forEach((g,gi)=>ALPHA.forEach((al,ai)=>FC.forEach(fc=>(fc?VK:[0]).forEach(vk=>grid.push([gi,ai,fc,vk])))));
  const cnt=new Int32Array(NB); bins.forEach(b=>cnt[b]+=R);
  function evalModel(gi,ai,fc,vk,kap){ const arr=[...cnt].map(c=>new Float32Array(c)), fill=new Int32Array(NB), v0=vN[ai], bt=B[gi];
    for(let i=0;i<N;i++){ const p=P[i], b=bins[i], kk=vk/(p.vc*KMS);
      for(let r=0;r<R;r++){ const j=i*R+r; let v=v0[j]*kap*bt[i];
        if(uc[j]<fc){ v*=Math.sqrt(1+0.4*um[j]); const k=kk*Math.exp(zk[j]); v=Math.sqrt(v*v+k*k+2*v*k*pk[j]); }
        const vo=Math.hypot(v+p.sx*g1[j],p.sx*g2[j]); if(vo>p.vmax+2*p.sx) continue;
        const vt=Math.sqrt(Math.max(vo*vo-2*p.sx*p.sx,0)); if(vt<3) arr[b][fill[b]++]=vt; } }
    return stats(arr.map((a,b)=>a.subarray(0,fill[b]))); }
  const chi=m=>{ const c=[0,0,0]; for(let b=0;b<NB;b++) for(let k=0;k<3;k++) c[k]+=((D[b][k]-m[b][k])**2)/(S[b][k]**2*(1+1/R)); return c; };
  const res=[]; let gi0=0; out.state='running'; out.total=grid.length; out.done=0;
  // background tabs throttle timers, so opts.manual exposes out.step(ms) to be driven from outside instead
  const step=(budget,auto)=>{ const t0=performance.now();
    while(gi0<grid.length && performance.now()-t0<budget){ const [gi,ai,fc,vk]=grid[gi0++]; let best=null;
      for(const kap of KAP){ const m=evalModel(gi,ai,fc,vk,kap), c=chi(m), tot=c[0]+c[1]+c[2]; if(!best||tot<best.tot) best={g:GR[gi],alpha:ALPHA[ai],fc,vk,kap,tot,c,m}; }
      res.push(best); out.done=gi0; }
    if(gi0<grid.length){ if(auto) setTimeout(()=>step(400,true),0); } else {
      const pick=f=>res.filter(f).sort((a,b)=>a.tot-b.tot)[0], r3=x=>+x.toFixed(4), slim=b=>b&&({g:b.g,alpha:b.alpha,fc:b.fc,vk:b.vk,kap:b.kap,chi2:r3(b.tot),parts:b.c.map(r3),med:b.m.map(x=>r3(x[0])),q90:b.m.map(x=>r3(x[1])),tail:b.m.map(x=>r3(x[2]))});
      out.summary={N,R,data:{med:D.map(x=>r3(x[0])),q90:D.map(x=>r3(x[1])),tail:D.map(x=>r3(x[2])),e_med:S.map(x=>r3(x[0])),e_q90:S.map(x=>r3(x[1])),e_tail:S.map(x=>r3(x[2])),n:dArr.map(a=>a.length)},
        best:Object.fromEntries(GR.map(g=>[g,slim(pick(r=>r.g===g))])),
        best_alpha13:Object.fromEntries(GR.map(g=>[g,slim(pick(r=>r.g===g&&r.alpha===1.3))])),
        best_noComp:Object.fromEntries(GR.map(g=>[g,slim(pick(r=>r.g===g&&r.fc===0))])),
        profile_fc:Object.fromEntries(GR.map(g=>[g,FC.map(fc=>{const b=pick(r=>r.g===g&&r.fc===fc); return b?[fc,r3(b.tot),b.alpha,b.vk]:null;})])),
        profile_alpha:Object.fromEntries(GR.map(g=>[g,ALPHA.map(al=>{const b=pick(r=>r.g===g&&r.alpha===al); return b?[al,r3(b.tot),b.fc,b.vk,b.kap]:null;})]))};
      out.state='done'; if(opts.then) setTimeout(opts.then,50); } };
  if(opts.manual) out.step=ms=>{ if(out.state!=='done') step(ms,false); return [out.state,out.done,out.total]; }; else step(400,true);
}
window.__wbm=window.__wbm||{}; window.__wbRun=run; window.__wbFns={beta,NU,orbitV,stats,rng,pairs};
// Runs behind the numbers in Entry 15 (driven in chunks with out.step(ms), because background tabs throttle timers):
//   run('rel25',0.25,null)                      four gravity laws, noise below 25% of the circular speed (39,602 pairs)
//   run('f25',0.25,null,fine())  and  run('f50',null,0.050,fine())
//       boost profile on both samples: gravity = Newton x [1 + mu x (MOND boost - 1)] with the galaxy-calibrated
//       interpolation (mu = 0 Newton, mu = 1 MOND), plus the "standard" and "simple" interpolations
const bR=beta(NU.rar), bS=beta(NU.standard);
const fine=()=>({manual:true,R:4, grav:[0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.7,1.0].map(mu=>['mu'+mu,y=>1+mu*(bR(y)-1)]).concat([['std1',bS],['simple',beta(NU.simple)]]), ALPHA:[1.0,1.3,1.6,2.0], FC:[0.1,0.15,0.2,0.3,0.4], VK:[0.15,0.3,0.6,1.2], KAP:[0.95,0.96,0.97,0.98,0.99,1.0,1.01,1.02]});
window.__wbFine=fine;
// e.g.  __wbRun('f50',null,0.050,__wbFine());  while(__wbm.f50.step(20000)[0]!=='done'){}  __wbm.f50.summary
})();
