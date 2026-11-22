// Filament-current test without the Galactic foreground map (Plasma Cosmology Workbench, Entry 15).
// The published residual RMs (RRM = RM - GRM) subtract a Galactic map that was itself built from these sources,
// so a real small-scale signal could be partly absorbed. Two checks:
//   (1) a cross-spine pair statistic on RAW RMs: differences between sources on opposite sides of the same spine,
//       closer than thetaMax on the sky, so any smooth foreground cancels and no map is used;
//   (2) a transfer test: inject a current into the raw RMs, subtract a foreground rebuilt from the same sources
//       (disc mean, with and without the source itself), and measure how much of the signal survives.
// Run in a browser console on vizier.cfa.harvard.edu with window.__rmT (GAVO lotsspol.rmtable CSV:
// cat_id,ra,dec,rm,rm_err,rrm,grm,grmerr,z_best,b,l) and window.__filT (Tempel+2014 table2 TSV) loaded.
(()=>{
const d2r=Math.PI/180, node=95*d2r, incl=32.5*d2r, AMP18=1.652;      // on-axis rad/m^2 per 1e18 A at n_e = 1e-5 cm^-3
const rm=window.__rmT.trim().split('\n').slice(1).map(l=>l.replace('\r','').split(',')).map(r=>({ra:+r[1],dec:+r[2],rm:+r[3],rrm:+r[5],grm:+r[6],z:r[8]===''?null:+r[8],b:+r[9]})).filter(s=>isFinite(s.rm)&&isFinite(s.rrm));
const pts=[]; for(const l of window.__filT.split('\n')){ if(!l||l[0]==='#'||l[0]==='-') continue; const v=l.trim().split(/\s+/); if(v.length<10||isNaN(+v[0])) continue; pts.push(v.map(Number)); }
const F=new Map(); for(const p of pts){ if(!F.has(p[0]))F.set(p[0],[]); F.get(p[0]).push(p);} const fils=[...F.values()].map(ps=>ps.sort((a,b)=>a[1]-b[1]));
function prep(raShift){
  const src=rm.filter(s=>Math.abs(s.rrm)<20&&Math.abs(s.b)>25).map(s=>{const r=(s.ra+raShift)*d2r,d=s.dec*d2r; const x=Math.cos(r-node)*Math.cos(d), y=Math.sin(r-node)*Math.cos(d), z=Math.sin(d);
    const lam=-Math.asin(x), eta=Math.atan2(z,y)-incl; return {rm:s.rm,rrm:s.rrm,z:s.z,n:[-Math.sin(lam),Math.cos(lam)*Math.cos(eta),Math.cos(lam)*Math.sin(eta)]};});
  const grid=new Map(); src.forEach((s,i)=>{const th=Math.acos(Math.max(-1,Math.min(1,s.n[2])))/d2r, ph=((Math.atan2(s.n[1],s.n[0])/d2r)+360)%360; const k=Math.floor(th)+','+Math.floor(ph); if(!grid.has(k))grid.set(k,[]); grid.get(k).push(i);});
  return {src,grid};
}
function assoc(S,Rh){ const {src,grid}=S, out=[];
  fils.forEach((ps,fi)=>{ const np=ps.length; if(np<3) return;
    const k=Math.min(3,Math.floor(np/2)); let fa=0,fb=0; for(let j=0;j<k;j++){fa+=ps[j][9]; fb+=ps[np-1-j][9];}
    let e=[ps[np-1][2]-ps[0][2],ps[np-1][3]-ps[0][3],ps[np-1][4]-ps[0][4]]; if(fa>fb) e=e.map(v=>-v);
    const best=new Map();
    ps.forEach((p,j)=>{ const P=[p[2],p[3],p[4]], D=p[5], a=(Rh*1.6)/D/d2r; const th=Math.acos(Math.max(-1,Math.min(1,P[2]/D)))/d2r, ph=((Math.atan2(P[1],P[0])/d2r)+360)%360; const pa=Math.min(180,a/Math.max(Math.sin(th*d2r),1e-3));
      for(let ti=Math.floor(th-a);ti<=Math.floor(th+a);ti++) for(let pj=Math.floor(ph-pa);pj<=Math.floor(ph+pa);pj++){ const L=grid.get(ti+','+(((pj%360)+360)%360)); if(!L) continue;
        for(const si of L){ const n=src[si].n, dot=P[0]*n[0]+P[1]*n[1]+P[2]*n[2]; if(dot<=0) continue; const dd=Math.hypot(P[0]-dot*n[0],P[1]-dot*n[1],P[2]-dot*n[2]); if(dd<Rh*1.6){const b=best.get(si); if(!b||dd<b.d) best.set(si,{d:dd,j});} } } });
    for(const [si,b] of best){ const j=b.j; if(j===0||j===np-1) continue; const P=[ps[j][2],ps[j][3],ps[j][4]]; let t=[ps[j][6],ps[j][7],ps[j][8]]; if(t[0]*e[0]+t[1]*e[1]+t[2]*e[2]<0) t=t.map(v=>-v);
      const n=src[si].n, c=[t[1]*n[2]-t[2]*n[1],t[2]*n[0]-t[0]*n[2],t[0]*n[1]-t[1]*n[0]], cn=Math.hypot(...c); if(cn<0.2) continue;
      const p=(P[0]*c[0]+P[1]*c[1]+P[2]*c[2])/cn; if(Math.abs(p)>=Rh) continue;
      const zf=ps[j][5]/2997.9, zs=src[si].z; if(zs!==null&&zs<zf+0.02) continue;
      out.push({fi,si,side:Math.sign(p),g:Math.atan(Math.sqrt(Rh*Rh-p*p)/Math.max(Math.abs(p),1e-6))*2/Math.PI,zf}); } });
  return out;
}
// node-oriented amplitude from single sources (the Entry 13 statistic), for any per-source field x[]
const statB=(A,x,flip)=>{let num=0,den=0; for(const a of A){const o=flip?flip[a.fi]:1; num+=o*a.side*a.g*x[a.si]; den+=a.g*a.g;} return num/den;};
// cross-spine pairs: sources on opposite sides of the same spine, closer than thetaMax degrees
function pairsOf(A,src,thetaMax){ const byF=new Map(); A.forEach(a=>{ if(!byF.has(a.fi)) byF.set(a.fi,[[],[]]); byF.get(a.fi)[a.side>0?0:1].push(a); });
  const cm=Math.cos(thetaMax*d2r), out=[]; for(const [fi,[p,m]] of byF) for(const a of p) for(const b of m){ const u=src[a.si].n, v=src[b.si].n; if(u[0]*v[0]+u[1]*v[1]+u[2]*v[2]>=cm) out.push({fi,i:a.si,j:b.si,w:a.g+b.g}); } return out; }
const statP=(Pr,x,flip)=>{let num=0,den=0; for(const p of Pr){const o=flip?flip[p.fi]:1; num+=o*p.w*(x[p.i]-x[p.j]); den+=p.w*p.w;} return den?num/den:NaN;};
const sd=a=>{const m=a.reduce((u,v)=>u+v,0)/a.length; return Math.sqrt(a.reduce((u,v)=>u+(v-m)**2,0)/a.length);};
function analyse(Rh,thetaMax){
  const S=prep(0), A=assoc(S,Rh), src=S.src, ns=src.length, nf=fils.length, rrm=src.map(s=>s.rrm), raw=src.map(s=>s.rm), Pr=pairsOf(A,src,thetaMax);
  const flips=[]; let seed=7; const rnd=()=>{seed=(seed*1664525+1013904223)>>>0; return seed/4294967296;};
  for(let r=0;r<1000;r++){const f=new Int8Array(nf); for(let i=0;i<nf;i++) f[i]=rnd()<0.5?-1:1; flips.push(f);}
  const nullOf=(st,X,x)=>sd(flips.map(f=>st(X,x,f)));
  const inj=I17=>{const d=new Float64Array(ns); for(const a of A) d[a.si]+=AMP18*(I17/10)*a.side*a.g/((1+a.zf)**2); return d;};
  const res={Rh,thetaMax,crossings:A.length,sources:new Set(A.map(a=>a.si)).size,pairs:Pr.length,pairFilaments:new Set(Pr.map(p=>p.fi)).size,
    B_rrm:statB(A,rrm), B_rrm_sd:nullOf(statB,A,rrm), P_rrm:statP(Pr,rrm), P_rrm_sd:nullOf(statP,Pr,rrm), P_raw:statP(Pr,raw), P_raw_sd:nullOf(statP,Pr,raw)};
  for(const I17 of [1,3]){ const d=inj(I17), rr=rrm.map((v,i)=>v+d[i]), rw=raw.map((v,i)=>v+d[i]);
    res['inj'+I17]={B_rrm:(statB(A,rr)-res.B_rrm), P_raw:(statP(Pr,rw)-res.P_raw), zB:(statB(A,rr)-res.B_rrm)/res.B_rrm_sd, zP:(statP(Pr,rw)-res.P_raw)/res.P_raw_sd}; }
  // transfer through a foreground rebuilt from the same sources: disc mean of raw RM within rad degrees
  const nb=rad=>{const cm=Math.cos(rad*d2r), L=src.map(()=>[]); for(let i=0;i<ns;i++) for(let j=i+1;j<ns;j++){const u=src[i].n,v=src[j].n; if(u[0]*v[0]+u[1]*v[1]+u[2]*v[2]>=cm){L[i].push(j);L[j].push(i);} } return L;};
  res.transfer={}; const d3=inj(3), sig=statB(A,Array.from(d3));
  for(const rad of [0.5,1,2]){ const L=nb(rad), nmean=L.reduce((u,l)=>u+l.length,0)/ns;
    const sub=(x,self)=>x.map((v,i)=>{const l=L[i]; let s=self?v:0,c=self?1:0; for(const j of l){s+=x[j];c++;} return c?v-s/c:v;});
    res.transfer[rad]={neighbours:nmean, leaveOneOut:statB(A,sub(Array.from(d3),false))/sig, includingSelf:statB(A,sub(Array.from(d3),true))/sig}; }
  return {res,src,A,Pr};
}
window.__ff={state:'running'};
setTimeout(()=>{ try{
  const main=analyse(0.7,1.5), out={main:main.res};
  // shifted-sky controls for the pair statistic on raw RMs
  const ctrl=[]; for(const sh of [-50,-38,-26,-14,14,26,38,50]){ const S=prep(sh), A=assoc(S,0.7), Pr=pairsOf(A,S.src,1.5); ctrl.push([statP(Pr,S.src.map(s=>s.rm)),statP(Pr,S.src.map(s=>s.rrm)),statB(A,S.src.map(s=>s.rrm)),Pr.length]); }
  out.controls=ctrl; out.ctrl_P_raw_sd=sd(ctrl.map(c=>c[0])); out.ctrl_P_rrm_sd=sd(ctrl.map(c=>c[1])); out.ctrl_B_rrm_sd=sd(ctrl.map(c=>c[2]));
  out.theta1=analyse(0.7,1.0).res; out.theta3=analyse(0.7,3.0).res;
  window.__ff={state:'done',out};
}catch(e){ window.__ff={state:'error',err:String(e&&e.stack||e)}; } },50);
})();
