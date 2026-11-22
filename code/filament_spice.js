/* Workbench Entry 68: the filament sign-step test of Entries 13 and 15, rerun on the SPICE-RACS first release.
   A straight current along a filament wraps a field around it; behind uniform gas the rotation measure then steps
   from -S to +S across the spine.  S = 1.652 rad/m^2 per 1e18 A at n_e = 1e-5 cm^-3, times the share g of the sightline's
   path inside the filament, divided by (1+z)^2.
   Data, held in browser memory only:
     window.SP    SPICE-RACS DR1 (Thomson et al. 2023), CASDA TAP table AS110.spice_racs_dr1_corrected_cut_v02:
                  components with snr, leakage, channel and Stokes-I fit flags all 0, merged per source_id (weighted mean)
                  -> {ra, dec, rm, e, b}
     window.FILS  Tempel et al. (2014) filament points, VizieR J/MNRAS/438/3465/table2: ID, IDpt, x, y, z, D, dx, dy, dz, fden
                  (SDSS survey Cartesian frame, Mpc/h), grouped by filament and sorted by IDpt.
   Statistics, as in Entry 15 (filament_foreground_check.js):
     B  single sources, side x g x residual RM, oriented toward each filament's denser end (the "node-oriented" current);
        residual = RM minus the median of other sources within R_fg degrees (leave-one-out).
     P  pairs of sources on opposite sides of the same spine within thetaMax degrees, raw RM: any smooth foreground cancels.
     Q  new: the same pairs with no orientation assumed: per filament, (sum of side x g x residual)^2 minus its own
        diagonal, summed; null from random side flips of single sources.  It responds to a current of either direction.
   Errors: 1,000 random flips of filament orientation (B, P) or of single-source sides (Q), and shifted-sky controls.     */
window.FS=(function(){
const d2r=Math.PI/180,node=95*d2r,incl=32.5*d2r,AMP18=1.652,H=2997.9;
let seed=7;const rnd=()=>{seed=(seed*1664525+1013904223)>>>0;return seed/4294967296;};
const nvec=(ra,dec)=>{const r=ra*d2r,d=dec*d2r,x=Math.cos(r-node)*Math.cos(d),y=Math.sin(r-node)*Math.cos(d),z=Math.sin(d),lam=-Math.asin(x),eta=Math.atan2(z,y)-incl;
  return [-Math.sin(lam),Math.cos(lam)*Math.cos(eta),Math.cos(lam)*Math.sin(eta)];};
const sd=a=>{const m=a.reduce((u,v)=>u+v,0)/a.length;return Math.sqrt(a.reduce((u,v)=>u+(v-m)**2,0)/a.length);};
const med=a=>{const s=a.slice().sort((x,y)=>x-y),n=s.length;return n%2?s[n>>1]:(s[n/2-1]+s[n/2])/2;};
function neighbours(src,rad){const cm=Math.cos(rad*d2r),L=src.map(()=>[]);for(let i=0;i<src.length;i++)for(let j=i+1;j<src.length;j++){const u=src[i].n0,v=src[j].n0;if(u[0]*v[0]+u[1]*v[1]+u[2]*v[2]>=cm){L[i].push(j);L[j].push(i);}}return L;}
function resid(x,L,minN){return x.map((v,i)=>L[i].length>=minN?v-med(L[i].map(j=>x[j])):NaN);}
function prep(raShift,opt){
  const src=window.SP.filter(s=>Math.abs(s.b)>25).map(s=>({ra:s.ra,dec:s.dec,rm:s.rm,e:s.e,n0:nvec(s.ra,s.dec),n:nvec(s.ra+raShift,s.dec)}));
  const grid=new Map();src.forEach((s,i)=>{const th=Math.acos(Math.max(-1,Math.min(1,s.n[2])))/d2r,ph=((Math.atan2(s.n[1],s.n[0])/d2r)+360)%360,k=Math.floor(th)+','+Math.floor(ph);if(!grid.has(k))grid.set(k,[]);grid.get(k).push(i);});
  return {src,grid};}
function assoc(S,Rh){const {src,grid}=S,out=[],fils=window.FILS;
  fils.forEach((ps,fi)=>{const np=ps.length;if(np<3)return;
    const k=Math.min(3,Math.floor(np/2));let fa=0,fb=0;for(let j=0;j<k;j++){fa+=ps[j][9];fb+=ps[np-1-j][9];}
    let e=[ps[np-1][2]-ps[0][2],ps[np-1][3]-ps[0][3],ps[np-1][4]-ps[0][4]];if(fa>fb)e=e.map(v=>-v);   // e points to the denser end
    const best=new Map();
    ps.forEach((p,j)=>{const P=[p[2],p[3],p[4]],D=p[5],a=(Rh*1.6)/D/d2r,th=Math.acos(Math.max(-1,Math.min(1,P[2]/D)))/d2r,ph=((Math.atan2(P[1],P[0])/d2r)+360)%360,pa=Math.min(180,a/Math.max(Math.sin(th*d2r),1e-3));
      for(let ti=Math.floor(th-a);ti<=Math.floor(th+a);ti++)for(let pj=Math.floor(ph-pa);pj<=Math.floor(ph+pa);pj++){const L=grid.get(ti+','+(((pj%360)+360)%360));if(!L)continue;
        for(const si of L){const n=src[si].n,dot=P[0]*n[0]+P[1]*n[1]+P[2]*n[2];if(dot<=0)continue;const dd=Math.hypot(P[0]-dot*n[0],P[1]-dot*n[1],P[2]-dot*n[2]);if(dd<Rh*1.6){const b=best.get(si);if(!b||dd<b.d)best.set(si,{d:dd,j});}}}});
    for(const [si,b] of best){const j=b.j;if(j===0||j===np-1)continue;const P=[ps[j][2],ps[j][3],ps[j][4]];let t=[ps[j][6],ps[j][7],ps[j][8]];if(t[0]*e[0]+t[1]*e[1]+t[2]*e[2]<0)t=t.map(v=>-v);
      const n=src[si].n,c=[t[1]*n[2]-t[2]*n[1],t[2]*n[0]-t[0]*n[2],t[0]*n[1]-t[1]*n[0]],cn=Math.hypot(...c);if(cn<0.2)continue;
      const p=(P[0]*c[0]+P[1]*c[1]+P[2]*c[2])/cn;if(Math.abs(p)>=Rh)continue;
      out.push({fi,si,side:Math.sign(p),g:Math.atan(Math.sqrt(Rh*Rh-p*p)/Math.max(Math.abs(p),1e-6))*2/Math.PI,zf:ps[j][5]/H});}});
  return out;}
const statB=(A,x,flip)=>{let num=0,den=0;for(const a of A){if(!isFinite(x[a.si]))continue;const o=flip?flip[a.fi]:1;num+=o*a.side*a.g*x[a.si];den+=a.g*a.g;}return num/den;};
function pairsOf(A,src,thetaMax){const byF=new Map();A.forEach(a=>{if(!byF.has(a.fi))byF.set(a.fi,[[],[]]);byF.get(a.fi)[a.side>0?0:1].push(a);});
  const cm=Math.cos(thetaMax*d2r),out=[];for(const [fi,[p,m]] of byF)for(const a of p)for(const b of m){const u=src[a.si].n,v=src[b.si].n;if(u[0]*v[0]+u[1]*v[1]+u[2]*v[2]>=cm)out.push({fi,i:a.si,j:b.si,w:a.g+b.g});}return out;}
const statP=(Pr,x,flip)=>{let num=0,den=0;for(const p of Pr){const o=flip?flip[p.fi]:1;num+=o*p.w*(x[p.i]-x[p.j]);den+=p.w*p.w;}return den?num/den:NaN;};
// Q: orientation-free.  Per filament s_f = sum side*g*x; Q = sum_f (s_f^2 - sum g^2 x^2) / sum_f sum g^2.  sides: optional per-entry multipliers
function statQ(A,x,sides){const byF=new Map();let den=0;A.forEach((a,k)=>{const v=x[a.si];if(!isFinite(v))return;const s=sides?sides[k]:1;if(!byF.has(a.fi))byF.set(a.fi,[0,0]);const o=byF.get(a.fi);o[0]+=s*a.side*a.g*v;o[1]+=(a.g*v)**2;den+=a.g*a.g;});
  let q=0;for(const [,o] of byF)q+=o[0]*o[0]-o[1];return q/den;}
function analyse(Rh,thetaMax,Rfg,nnull){
  nnull=nnull||1000;const S=prep(0),src=S.src,ns=src.length,A=assoc(S,Rh),L=neighbours(src,Rfg),raw=src.map(s=>s.rm),rrm0=resid(raw,L,4),rrm=rrm0.map(v=>Math.abs(v)<20?v:NaN),Pr=pairsOf(A,src,thetaMax),nf=window.FILS.length;
  seed=7;const flips=[];for(let r=0;r<nnull;r++){const f=new Int8Array(nf);for(let i=0;i<nf;i++)f[i]=rnd()<0.5?-1:1;flips.push(f);}
  const sflips=[];for(let r=0;r<nnull;r++){const f=new Int8Array(A.length);for(let i=0;i<A.length;i++)f[i]=rnd()<0.5?-1:1;sflips.push(f);}
  const res={Rh,thetaMax,Rfg,sources:ns,withResid:rrm.filter(isFinite).length,crossings:A.length,crossSources:new Set(A.map(a=>a.si)).size,filaments:new Set(A.map(a=>a.fi)).size,pairs:Pr.length,pairFilaments:new Set(Pr.map(p=>p.fi)).size,
    medZ:med(A.map(a=>a.zf)),medG:med(A.map(a=>a.g)),
    B:statB(A,rrm),B_sd:sd(flips.map(f=>statB(A,rrm,f))),P:statP(Pr,raw),P_sd:sd(flips.map(f=>statP(Pr,raw,f))),P_rrm:statP(Pr,rrm.map(v=>isFinite(v)?v:0)),
    Q:statQ(A,rrm)};
  const qn=sflips.map(f=>statQ(A,rrm,f));res.Q_mean0=qn.reduce((u,v)=>u+v,0)/qn.length;res.Q_sd=sd(qn);res.Q_p=qn.filter(v=>v>=res.Q).length/qn.length;
  // injection through the whole chain: a current added to the raw RMs, foreground rebuilt
  const inj=(I17,dirs)=>{const d=new Float64Array(ns);for(const a of A)d[a.si]+=(dirs?dirs[a.fi]:1)*AMP18*(I17/10)*a.side*a.g/((1+a.zf)**2);return raw.map((v,i)=>v+d[i]);};
  for(const I17 of [3,10]){const rw=inj(I17),rr=resid(rw,L,4).map((v,i)=>isFinite(rrm[i])?v:NaN);
    const dirs=new Int8Array(nf);for(let i=0;i<nf;i++)dirs[i]=rnd()<0.5?-1:1;const rw2=inj(I17,dirs),rr2=resid(rw2,L,4).map((v,i)=>isFinite(rrm[i])?v:NaN);
    res['inj'+I17]={dB:statB(A,rr)-res.B,dP:statP(Pr,rw)-res.P,dQ_random_dirs:statQ(A,rr2)-res.Q};}
  res.respB=res.inj10.dB/10;res.respP=res.inj10.dP/10;           // per 1e17 A
  res.limB=(Math.max(res.B,0)+2*res.B_sd)/res.respB;res.limP=(Math.max(res.P,0)+2*res.P_sd)/res.respP;
  // Q scales as I^2: limit where the injected excess reaches the observed Q plus 2 sd above the null mean
  const q10=res.inj10.dQ_random_dirs;res.limQ=10*Math.sqrt(Math.max(res.Q-res.Q_mean0+2*res.Q_sd,2*res.Q_sd)/Math.max(q10,1e-9));
  return res;}
function controls(Rh,thetaMax,Rfg){const out=[];for(const sh of [-50,-38,-26,-14,14,26,38,50]){const S=prep(sh),src=S.src,L=neighbours(src,Rfg),raw=src.map(s=>s.rm),rrm=resid(raw,L,4).map(v=>Math.abs(v)<20?v:NaN),A=assoc(S,Rh),Pr=pairsOf(A,src,thetaMax);
  out.push({sh,crossings:A.length,pairs:Pr.length,B:statB(A,rrm),P:statP(Pr,raw),Q:statQ(A,rrm)});}return out;}
return {nvec,prep,assoc,statB,pairsOf,statP,statQ,analyse,controls,neighbours,resid,AMP18};
})();
