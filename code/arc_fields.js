/* Workbench Entry 66: the Giant Arc's own recipe, run on real absorbers from an independent catalogue,
   in the Arc's field and in 53 other volumes of the same sky.

   Catalogue: Seyffert et al. (2013), VizieR J/ApJ/779/161 (SDSS DR7 quasars; 79,595 sightlines, 35,629 Mg II doublets).
   The Giant Arc (Lopez et al. 2022) was found in the Zhu & Menard catalogue (DR7 + DR12) with a different line finder.

   Recipe, as published: slice z = 0.802 +- 0.060; sightlines ("probes") are quasars above the slice; tangent-plane
   coordinates about the cluster PSZ2 G069.39+68.05, scaled to comoving Mpc at the central redshift, flat cosmology with
   Omega_m = 0.27 and H0 = 70; points are linked within 95 Mpc x (504 / N)^(1/3), the paper's own density scaling
   (504 absorbers in 1541 x 1615 x 338 Mpc, as given by Sawala et al. 2025); a group's "significance" is the convex-hull
   test of Clowes et al. (2012), as rebuilt in Entry 52 (hull_significance.py, whose reference table is copied below).

   Null: the absorbers of a field are moved to other sightlines of the same field (weighted by how often a quasar of that
   redshift shows an absorber in the slice), keeping their depths.  A second null also redraws each plate's detection rate.

   To run: open https://vizier.cfa.harvard.edu/viz-bin/VizieR in a browser, paste this file in the console, then
       await GA.fetchTables();  GA.runArc();  GA.runControls(-1) ... GA.runControls(5);  GA.plates();  GA.runDeep();
   Nothing is saved to disk.  Results as published are in arc_fields.json.                                              */
window.GA=(function(){
const OM=0.27,H0=70,CH=299792.458/H0,DZ=0.0005,D2R=Math.PI/180;
const E=z=>Math.sqrt(OM*Math.pow(1+z,3)+1-OM);
const T=[0];{let s=0;for(let i=1;i<=7000;i++){const a=(i-1)*DZ,b=i*DZ;s+=DZ*(1/E(a)+4/E((a+b)/2)+1/E(b))/6;T.push(CH*s);}}
const dc=z=>{const u=z/DZ,i=Math.floor(u);return T[i]+(T[i+1]-T[i])*(u-i);};
const zof=d=>{let lo=0,hi=3.4;for(let k=0;k<50;k++){const m=(lo+hi)/2;if(dc(m)<d)lo=m;else hi=m;}return (lo+hi)/2;};
const sex=(s,ra)=>{const p=s.trim().split(/\s+/);const sg=p[0][0]==='-'?-1:1;const v=Math.abs(+p[0])+(+p[1])/60+(+p[2])/3600;return ra?15*v:sg*v;};
const G={OM,H0,dc,zof,sex};
G.CLUSTER=[215.40494,38.36119];                      // PSZ2 G069.39+68.05 (l = 69.39, b = +68.05) in RA, Dec
G.BOX=[1541,1615,338.2];G.RHO0=504/(1541*1615*338.2);G.L0=95;

/* ---------- data ---------- */
G.fetchTables=async function(){const ft=async u=>{for(let k=0;k<3;k++){try{const r=await fetch(u);if(r.ok)return await r.text();}catch(e){}}throw new Error('fetch failed');};
 const parse=t=>{const L=t.split('\n').filter(l=>l.length&&!l.startsWith('#'));return {hdr:L[0].split('\t'),rows:L.slice(3).map(l=>l.split('\t'))};};
 const T2=parse(await ft('/viz-bin/asu-tsv?-source=J/ApJ/779/161/table2&-out=QSO,zQSO,z2796,f_z2796,W2796,e_W2796,Code&-out.max=unlimited&-oc.form=dec'));
 const T1=parse(await ft('/viz-bin/asu-tsv?-source=J/ApJ/779/161/table1&-out=QSO,RAJ2000,DEJ2000,zQSO,<S/N>,fBAL,dXmax,N2,N,PapI,DR7&-out.max=unlimited&-oc.form=dec'));
 G.D=G.load(T1,T2);return {quasars:G.D.qs.length,absorbers:G.D.ab.length,unmatched:G.D.miss};};
G.load=function(T1,T2){const Q=new Map(),qs=[];
 T1.rows.forEach(r=>{const q={id:r[0].trim(),ra:sex(r[1],1),de:sex(r[2],0),zq:+r[3],abs:[]};Q.set(q.id,q);qs.push(q);});
 let miss=0;const ab=[];
 T2.rows.forEach(r=>{const q=Q.get(r[0].trim());if(!q){miss++;return;}const a={q,z:+r[2],w:+r[4],code:+r[6]};q.abs.push(a);ab.push(a);});
 return {qs,ab,miss};};

/* ---------- geometry, fields, linking ---------- */
G.tp=(ra,de,ra0,de0)=>{const a=(ra-ra0)*D2R,d=de*D2R,d0=de0*D2R;const c=Math.sin(d0)*Math.sin(d)+Math.cos(d0)*Math.cos(d)*Math.cos(a);return [Math.cos(d)*Math.sin(a)/c,(Math.cos(d0)*Math.sin(d)-Math.sin(d0)*Math.cos(d)*Math.cos(a))/c,c];};
G.field=function(D,ra0,de0,zlo,zhi,zc,opt){opt=opt||{};const hx=opt.hx||G.BOX[0]/2,hy=opt.hy||G.BOX[1]/2,dC=dc(zc),pr=[],pts=[];
 for(const q of D.qs){if(q.zq<=zhi)continue;const t=G.tp(q.ra,q.de,ra0,de0);if(t[2]<=0)continue;const x=dC*t[0],y=dC*t[1];if(Math.abs(x)>hx||Math.abs(y)>hy)continue;const p={q,x,y,dz:[],wt:1};
  for(const a of q.abs){if(a.z>=zlo&&a.z<=zhi&&(!opt.wmin||a.w>=opt.wmin)){pts.push([x,y,dc(a.z)-dC,pr.length]);p.dz.push(dc(a.z)-dC);}}pr.push(p);}
 return {pr,pts,hx,hy,dC,zlo,zhi,zc,depth:dc(zhi)-dc(zlo)};};
G.fof=function(pts,link){const n=pts.length,par=new Int32Array(n);for(let i=0;i<n;i++)par[i]=i;const f=i=>{while(par[i]!==i){par[i]=par[par[i]];i=par[i];}return i;};const l2=link*link;
 for(let i=0;i<n;i++){const a=pts[i];for(let j=i+1;j<n;j++){const b=pts[j],dx=a[0]-b[0];if(dx>link||dx<-link)continue;const dy=a[1]-b[1];if(dy>link||dy<-link)continue;const dz=a[2]-b[2];if(dx*dx+dy*dy+dz*dz<=l2){const ra=f(i),rb=f(j);if(ra!==rb)par[ra]=rb;}}}
 const g=new Map();for(let i=0;i<n;i++){const r=f(i);if(!g.has(r))g.set(r,[]);g.get(r).push(i);}return [...g.values()].sort((a,b)=>b.length-a.length);};
G.extent=(pts,idx)=>{let m=0;for(let i=0;i<idx.length;i++)for(let j=i+1;j<idx.length;j++){const a=pts[idx[i]],b=pts[idx[j]],d=Math.hypot(a[0]-b[0],a[1]-b[1]);if(d>m)m=d;}return m;};
G.linkFor=(n,F)=>G.L0*Math.cbrt(G.RHO0/(n/(4*F.hx*F.hy*F.depth)));
G.rng=function(seed){let s=seed>>>0;return ()=>{s=(s+0x6D2B79F5)|0;let t=Math.imul(s^s>>>15,1|s);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};};
G.slice=function(k){const d0=dc(0.742)+k*G.BOX[2];if(k===0)return {k,zlo:0.742,zhi:0.862,zc:0.802};return {k,zlo:zof(d0),zhi:zof(d0+G.BOX[2]),zc:zof(d0+G.BOX[2]/2)};};

/* ---------- the sightline null ---------- */
G.rates=function(D,zlo,zhi){const bw=0.25,nb=16,c=[];for(let i=0;i<nb;i++)c.push([0,0]);   // chance of an absorber in the slice, by quasar redshift
 for(const q of D.qs){if(q.zq<=zhi)continue;const i=Math.min(nb-1,Math.floor((q.zq-zhi)/bw));c[i][0]++;for(const a of q.abs)if(a.z>=zlo&&a.z<=zhi)c[i][1]++;}
 const tot=c.reduce((s,x)=>[s[0]+x[0],s[1]+x[1]],[0,0]),rb=tot[1]/tot[0];return {bw,nb,r:c.map(x=>(x[1]+5*rb)/(x[0]+5)),rb,tot};};
G.wOf=(p,R,zhi)=>R.r[Math.min(R.nb-1,Math.floor((p.q.zq-zhi)/R.bw))];
G.plateOf=q=>q.id.split('-')[1];
G.setW=function(F,R,fac){for(const p of F.pr)p.wt=G.wOf(p,R,F.zhi)*(fac?(fac.get(G.plateOf(p.q))||1):1);};
G.shuffler=function(F){const cum=new Float64Array(F.pr.length);let s=0;for(let i=0;i<F.pr.length;i++){s+=F.pr[i].wt;cum[i]=s;}
 const units=F.pr.filter(p=>p.dz.length).map(p=>p.dz);                                  // a sightline's absorbers move together
 return rand=>{const used=new Set(),pts=[];for(const u of units){let k;do{const t=rand()*s;let lo=0,hi=cum.length-1;while(lo<hi){const m=(lo+hi)>>1;if(cum[m]<t)lo=m+1;else hi=m;}k=lo;}while(used.has(k));used.add(k);const p=F.pr[k];for(const dz of u)pts.push([p.x,p.y,dz,k]);}return pts;};};
G.region=function(F,inside){let sw=0,swA=0,nA=0,n=0;for(const p of F.pr){sw+=p.wt;n+=p.dz.length;if(inside(p.x,p.y)){swA+=p.wt;nA+=p.dz.length;}}return {nA,mu:n*swA/sw,n};};
G.band=(F,h)=>G.region(F,(x,y)=>Math.abs(x)<=500&&Math.abs(y)<=h);
G.pois=function(k,mu){let term=Math.exp(-mu),cdf=0;for(let i=0;i<k;i++){cdf+=term;term*=mu/(i+1);}return Math.max(1-cdf,0);};

/* ---------- the published test (Entry 52's hull_significance.py, ported; checked against scipy to 4 decimals) ---------- */
G.REF=[[5,.562,.548],[6,.596,.455],[7,.617,.400],[8,.645,.366],[10,.674,.302],[12,.727,.278],[15,.765,.244],[20,.816,.187],[25,.866,.162],[30,.893,.140],[36,.922,.127],[44,.950,.106],[52,.972,.097],[60,.994,.082],[72,1.011,.074],[85,1.034,.066],[100,1.045,.057],[120,1.064,.051],[150,1.081,.042],[190,1.094,.036]];
G.refAt=function(n){const R=G.REF,l=Math.log(n);if(n<=R[0][0])return [R[0][1],R[0][2]];for(let i=1;i<R.length;i++){if(n<=R[i][0]){const a=Math.log(R[i-1][0]),b=Math.log(R[i][0]),t=(l-a)/(b-a);return [R[i-1][1]+t*(R[i][1]-R[i-1][1]),R[i-1][2]+t*(R[i][2]-R[i-1][2])];}}const e=R[R.length-1];return [e[1],e[2]];};
G.mstMean=function(P){const n=P.length,d=new Float64Array(n).fill(Infinity),inT=new Uint8Array(n);let cur=0,sum=0;inT[0]=1;for(let s=1;s<n;s++){let best=-1,bd=Infinity;for(let j=0;j<n;j++){if(inT[j])continue;const dx=P[cur][0]-P[j][0],dy=P[cur][1]-P[j][1],dz=P[cur][2]-P[j][2],dd=Math.sqrt(dx*dx+dy*dy+dz*dz);if(dd<d[j])d[j]=dd;if(d[j]<bd){bd=d[j];best=j;}}inT[best]=1;sum+=bd;cur=best;}return sum/(n-1);};
G.hull=function(P){const n=P.length;let cx=0,cy=0,cz=0;for(const p of P){cx+=p[0];cy+=p[1];cz+=p[2];}cx/=n;cy/=n;cz/=n;let V=0,A=0;const edges=new Map();
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(let k=j+1;k<n;k++){const a=P[i],b=P[j],c=P[k];const ux=b[0]-a[0],uy=b[1]-a[1],uz=b[2]-a[2],vx=c[0]-a[0],vy=c[1]-a[1],vz=c[2]-a[2];let nx=uy*vz-uz*vy,ny=uz*vx-ux*vz,nz=ux*vy-uy*vx;const nl=Math.sqrt(nx*nx+ny*ny+nz*nz);if(nl<1e-9)continue;nx/=nl;ny/=nl;nz/=nl;let pos=false,neg=false;
   for(let m=0;m<n;m++){if(m===i||m===j||m===k)continue;const s=(P[m][0]-a[0])*nx+(P[m][1]-a[1])*ny+(P[m][2]-a[2])*nz;if(s>1e-9)pos=true;else if(s<-1e-9)neg=true;if(pos&&neg)break;}if(pos&&neg)continue;
   if(pos){nx=-nx;ny=-ny;nz=-nz;}
   const area=nl/2,h=(a[0]-cx)*nx+(a[1]-cy)*ny+(a[2]-cz)*nz;V+=area*h/3;A+=area;
   for(const e of [[i,j],[j,k],[i,k]]){const key=e[0]*1000+e[1];if(!edges.has(key))edges.set(key,[]);edges.get(key).push([nx,ny,nz]);}}
 let M=0;for(const [key,ns] of edges){if(ns.length!==2)continue;const i=Math.floor(key/1000),j=key%1000,L=Math.hypot(P[i][0]-P[j][0],P[i][1]-P[j][1],P[i][2]-P[j][2]);const c=Math.max(-1,Math.min(1,ns[0][0]*ns[1][0]+ns[0][1]*ns[1][1]+ns[0][2]*ns[1][2]));M+=L*Math.acos(c);}
 return {V,A,M:M/2};};
G.chms=function(P,rho){const n=P.length,e=G.mstMean(P),r=e/2,h=G.hull(P),V=h.V+h.A*r+h.M*r*r+4*Math.PI*r*r*r/3,ref=G.refAt(n),mu=ref[0]*n/rho,sd=ref[1]*mu;return {n,e,V,sig:(mu-V)/sd,over:mu/V-1,hull:h};};
G.bigSig=function(pts,F,link){const g=G.fof(pts,link),P=g[0].map(i=>pts[i]),rho=pts.length/(4*F.hx*F.hy*F.depth);if(P.length<5)return {n:P.length,sig:NaN,over:NaN,e1:0,g};const c=G.chms(P,rho);return {n:c.n,sig:c.sig,over:c.over,e:c.e,V:c.V,e1:G.extent(pts,g[0]),g};};
const qt=(v,f)=>{const s=v.slice().sort((x,y)=>x-y);return s[Math.floor(f*s.length)];};
G.summ=function(obs,N,S,E,n,link){const nsh=N.length;return {n,link,n1:obs.n,e1:obs.e1,sig:obs.sig,over:obs.over,pN:N.filter(x=>x>=obs.n).length/nsh,pS:S.filter(x=>x>=obs.sig).length/nsh,pE:E.filter(x=>x>=obs.e1).length/nsh,medN:qt(N,.5),medS:qt(S,.5),s16:qt(S,.16),s84:qt(S,.84),s95:qt(S,.95),ge45:S.filter(x=>x>=4.5).length/nsh,N,S,E};};
G.calib=function(F,nsh,seed){const n=F.pts.length,link=G.linkFor(n,F),obs=G.bigSig(F.pts,F,link),sh=G.shuffler(F),rand=G.rng(seed),N=[],S=[],E=[];
 for(let k=0;k<nsh;k++){const b=G.bigSig(sh(rand),F,link);N.push(b.n);S.push(b.sig);E.push(b.e1);}return G.summ(obs,N,S,E,n,link);};

/* ---------- plates: does the detection rate vary from plate to plate? ---------- */
G.plateModel=function(D,S,margin){const zl=zof(dc(S.zlo)-margin),zh=zof(dc(S.zhi)+margin),R=G.rates(D,S.zlo,S.zhi);
 const nb=50,bw=0.1,c=[];for(let i=0;i<nb;i++)c.push([0,0]);const out=q=>{let n=0;for(const a of q.abs)if(a.z<zl||a.z>zh)n++;return n;};
 for(const q of D.qs){const i=Math.min(nb-1,Math.floor(q.zq/bw));c[i][0]++;c[i][1]+=out(q);}
 const P=new Map();
 for(const q of D.qs){const pl=G.plateOf(q);if(!P.has(pl))P.set(pl,{pl,nq:0,nOut:0,muOut:0,nIn:0,wIn:0});const o=P.get(pl),i=Math.min(nb-1,Math.floor(q.zq/bw));o.nq++;o.nOut+=out(q);o.muOut+=c[i][1]/Math.max(1,c[i][0]);
  if(q.zq>S.zhi){o.wIn+=R.r[Math.min(R.nb-1,Math.floor((q.zq-S.zhi)/R.bw))];for(const a of q.abs)if(a.z>=S.zlo&&a.z<=S.zhi)o.nIn++;}}
 let tn=0,tw=0;for(const o of P.values()){tn+=o.nIn;tw+=o.wIn;}for(const o of P.values())o.muIn=tn*o.wIn/tw;return {P,R,zl,zh};};
G.kappaOf=function(M){const pls=[...M.P.values()].filter(o=>o.nq>=10&&o.muOut>3);let s=0,mu=0;for(const o of pls){s+=(o.nOut-o.muOut)**2/o.muOut;mu+=o.muOut;}const idx=s/pls.length,mbar=mu/pls.length;return {idx,mbar,kappa:mbar/Math.max(idx-1,0.01)};};
G.factors=function(M,kappa){const f=new Map();for(const o of M.P.values())f.set(o.pl,(o.nOut+kappa)/(o.muOut+kappa));return f;};
G.poisDraw=(lam,rand)=>{const L=Math.exp(-lam);let k=0,p=1;do{k++;p*=rand();}while(p>L);return k-1;};
G.plateDisp=function(D,S,nsim,seed){const M=G.plateModel(D,S,150),pls=[...M.P.values()].filter(o=>o.nq>=10&&o.muIn>0.2);let s=0;for(const o of pls)s+=(o.nIn-o.muIn)**2/o.muIn;const obs=s/pls.length,rand=G.rng(seed),sim=[];
 for(let m=0;m<nsim;m++){let t=0;for(const o of pls){const c=G.poisDraw(o.muIn,rand);t+=(c-o.muIn)**2/o.muIn;}sim.push(t/pls.length);}sim.sort((a,b)=>a-b);
 const p2=[...M.P.values()].filter(o=>o.nq>=10&&o.muOut>3).sort((a,b)=>a.nOut/a.muOut-b.nOut/b.muOut),quint=[];
 for(let i=0;i<5;i++){const sub=p2.slice(Math.floor(i*p2.length/5),Math.floor((i+1)*p2.length/5)),q=[0,0,0,0];for(const o of sub){q[0]+=o.nOut;q[1]+=o.muOut;q[2]+=o.nIn;q[3]+=o.muIn;}quint.push(q);}
 return {k:S.k,plates:pls.length,obs,simMed:sim[nsim>>1],sim95:sim[Math.floor(0.95*nsim)],simMax:sim[nsim-1],quint,kappa:G.kappaOf(M)};};
G.gamma=function(shape,rand){const d=shape-1/3,c=1/Math.sqrt(9*d);for(;;){let x,v;do{const u1=rand()||1e-12,u2=rand();x=Math.sqrt(-2*Math.log(u1))*Math.cos(2*Math.PI*u2);v=1+c*x;}while(v<=0);v=v*v*v;const u=rand();if(u<1-0.0331*x*x*x*x)return d*v;if(Math.log(u)<0.5*x*x+d*(1-v+Math.log(v)))return d*v;}};
G.shufflerRE=function(F,M,kappa,vExtra){const base=F.pr.map(p=>G.wOf(p,M.R,F.zhi)),pl=F.pr.map(p=>G.plateOf(p.q)),names=[...new Set(pl)],idx=new Map(names.map((n,i)=>[n,i])),pi=pl.map(n=>idx.get(n)),po=names.map(n=>M.P.get(n));
 const units=F.pr.filter(p=>p.dz.length).map(p=>p.dz),cum=new Float64Array(F.pr.length),f=new Float64Array(names.length);
 return rand=>{for(let j=0;j<names.length;j++){let v=G.gamma(kappa+po[j].nOut,rand)/(kappa+po[j].muOut);if(vExtra>0)v*=G.gamma(1/vExtra,rand)*vExtra;f[j]=v;}let s=0;for(let i=0;i<cum.length;i++){s+=base[i]*f[pi[i]];cum[i]=s;}
  const used=new Set(),pts=[];for(const u of units){let k;do{const t=rand()*s;let lo=0,hi=cum.length-1;while(lo<hi){const m=(lo+hi)>>1;if(cum[m]<t)lo=m+1;else hi=m;}k=lo;}while(used.has(k));used.add(k);const p=F.pr[k];for(const dz of u)pts.push([p.x,p.y,dz,k]);}return pts;};};
G.calibRE=function(F,M,kappa,vExtra,nsh,seed){const n=F.pts.length,link=G.linkFor(n,F),obs=G.bigSig(F.pts,F,link),sh=G.shufflerRE(F,M,kappa,vExtra),rand=G.rng(seed),N=[],S=[],E=[];
 for(let k=0;k<nsh;k++){const b=G.bigSig(sh(rand),F,link);N.push(b.n);S.push(b.sig);E.push(b.e1);}return G.summ(obs,N,S,E,n,link);};

/* ---------- control fields ---------- */
G.cover=function(F){const c=new Int32Array(16);for(const p of F.pr){const i=Math.min(3,Math.floor((p.x+F.hx)/(F.hx/2))),j=Math.min(3,Math.floor((p.y+F.hy)/(F.hy/2)));c[j*4+i]++;}const s=Array.from(c).sort((a,b)=>a-b);return {min:s[0],med:(s[7]+s[8])/2};};
G.tiles=function(D,S,opt){opt=opt||{};const cand=[];                                     // centres on a 4-degree grid; a field must be fully covered
 for(let de=-4;de<=70;de+=4)for(let ra=104;ra<=268;ra+=4){const F=G.field(D,ra,de,S.zlo,S.zhi,S.zc);if(F.pr.length<(opt.minpr||1600))continue;const c=G.cover(F);if(c.med<(opt.minmed||100)||c.min<0.5*c.med)continue;cand.push({ra,de,F});}
 const picked=[],used=new Set(opt.exclude||[]);                                          // keep fields that share no sightlines
 for(const t of cand){let sh=0;for(const p of t.F.pr)if(used.has(p.q.id))sh++;if(sh>0.01*t.F.pr.length)continue;picked.push(t);for(const p of t.F.pr)used.add(p.q.id);}
 return picked;};

/* ---------- runs ---------- */
G.arcField=opt=>G.field(G.D,G.CLUSTER[0],G.CLUSTER[1],0.742,0.862,0.802,opt);
const short=c=>({n:c.n,link:+c.link.toFixed(1),n1:c.n1,e1:Math.round(c.e1),sig:+c.sig.toFixed(2),over:+c.over.toFixed(2),pN:c.pN,pS:c.pS,pE:c.pE,medN:c.medN,medS:+c.medS.toFixed(2),s16:+c.s16.toFixed(2),s84:+c.s84.toFixed(2),s95:+c.s95.toFixed(2),ge45:c.ge45});
G.relink=function(F,nsh,seed){const link=G.linkFor(F.pts.length,F),sc=link/G.L0;                // does the group hold together at shorter links?
 const brk=pts=>{const g=G.fof(pts,link)[0],P=g.map(i=>pts[i]),a=G.fof(P,65*sc),b=G.fof(P,85*sc);return [g.length,a[0].length/g.length,b[0].length/g.length,(b[0].length+(b[1]?b[1].length:0))/g.length];};
 const real=brk(F.pts),sh=G.shuffler(F),rand=G.rng(seed),B=[];for(let m=0;m<nsh;m++){const r=brk(sh(rand));if(r[0]>=real[0]-6&&r[0]<=real[0]+6)B.push(r);}const col=i=>B.map(r=>r[i]);
 return {real,nsim:B.length,f65:[qt(col(1),.16),qt(col(1),.5),qt(col(1),.84)],f85:[qt(col(2),.16),qt(col(2),.5),qt(col(2),.84)],p65:col(1).filter(x=>x>=real[1]).length/B.length,p85:col(2).filter(x=>x>=real[2]).length/B.length};};
G.runArc=function(){const D=G.D,R=G.rates(D,0.742,0.862),out={};
 const one=(opt,nsh,seed)=>{const F=G.arcField(opt);G.setW(F,R,null);const c=G.calib(F,nsh,seed),g=G.fof(F.pts,c.link)[0],xs=g.map(i=>F.pts[i][0]),ys=g.map(i=>F.pts[i][1]),b1=G.band(F,100),b2=G.band(F,150);
  return Object.assign(short(c),{probes:F.pr.length,box:[Math.min(...xs),Math.max(...xs),Math.min(...ys),Math.max(...ys)].map(Math.round),b100:[b1.nA,+b1.mu.toFixed(1),+G.pois(b1.nA,b1.mu).toFixed(4)],b150:[b2.nA,+b2.mu.toFixed(1),+G.pois(b2.nA,b2.mu).toFixed(4)]});};
 out.base=one({},1000,611);out.swap=one({hx:807.5,hy:770.5},500,802);out.strong=one({wmin:0.6},2000,811);out.w1=one({wmin:1.0},500,804);
 const F=G.arcField({});G.setW(F,R,null);G.F=F;
 const M=G.plateModel(D,G.slice(0),150),kp=G.kappaOf(M);out.platePersist=short(G.calibRE(F,M,kp.kappa,0,1000,613));out.plateFull=short(G.calibRE(F,M,kp.kappa,0.024,1000,614));
 out.relink=G.relink(F,2000,650);
 out.strips=[];for(let y0=-800;y0<800;y0+=100){const r=G.region(F,(x,y)=>Math.abs(x)<=500&&y>=y0&&y<y0+100);out.strips.push([y0,r.nA,+r.mu.toFixed(1)]);}
 out.rot=[];for(let a=0;a<180;a+=15){const c=Math.cos(a*D2R),s=Math.sin(a*D2R),r=G.region(F,(x,y)=>Math.abs(x*c+y*s)<=500&&Math.abs(-x*s+y*c)<=100);out.rot.push([a,r.nA,+r.mu.toFixed(1)]);}
 out.other=[];for(const k of [-2,-1,1,2,3,4]){const Sk=G.slice(k),Rk=G.rates(D,Sk.zlo,Sk.zhi),f=dc(Sk.zc)/dc(0.802),Fk=G.field(D,G.CLUSTER[0],G.CLUSTER[1],Sk.zlo,Sk.zhi,Sk.zc,{hx:770.5*f,hy:807.5*f});G.setW(Fk,Rk,null);
  const b1=G.region(Fk,(x,y)=>Math.abs(x)<=500*f&&Math.abs(y)<=100*f),b2=G.region(Fk,(x,y)=>Math.abs(x)<=500*f&&Math.abs(y)<=150*f);out.other.push([k,b1.nA,+b1.mu.toFixed(1),b2.nA,+b2.mu.toFixed(1)]);}
 return out;};
G.runControls=function(k){const D=G.D,S=G.slice(k),R=G.rates(D,S.zlo,S.zhi),opt=k>=5?{minpr:700,minmed:40}:{};if(k===0)opt.exclude=G.arcField({}).pr.map(p=>p.q.id);
 return G.tiles(D,S,opt).map((t,i)=>{const F=t.F;G.setW(F,R,null);const c=G.calib(F,400,7000+100*(k+2)+i),b1=G.band(F,100),b2=G.band(F,150);
  return [k,t.ra,t.de,c.n,Math.round(c.link),c.n1,Math.round(c.e1),+c.sig.toFixed(2),+c.over.toFixed(2),+c.pN.toFixed(3),+c.pS.toFixed(3),+c.medS.toFixed(2),b1.nA,+b1.mu.toFixed(1),b2.nA,+b2.mu.toFixed(1)];});};
G.runDeep=function(){const D=G.D,S=G.slice(-1),R=G.rates(D,S.zlo,S.zhi),M=G.plateModel(D,S,150),kp=G.kappaOf(M);        // the two strongest control fields, 2,000 shuffles each
 return G.tiles(D,S,{}).slice(0,2).map((t,i)=>{const F=t.F;G.setW(F,R,null);const c=G.calib(F,2000,7100+i),c2=G.calibRE(F,M,kp.kappa,0.024,2000,7200+i),g=G.fof(F.pts,c.link)[0],q=g.map(j=>F.pr[F.pts[j][3]].q);
  return {ra:t.ra,de:t.de,n1:c.n1,e1:Math.round(c.e1),sig:+c.sig.toFixed(2),pN:c.pN,pS:c.pS,s95:+c.s95.toFixed(2),max:+Math.max(...c.S).toFixed(2),pS_plate:c2.pS,sky:[Math.min(...q.map(x=>x.ra)),Math.max(...q.map(x=>x.ra)),Math.min(...q.map(x=>x.de)),Math.max(...q.map(x=>x.de))].map(x=>+x.toFixed(1))};});};
G.plates=function(){return [-1,0,1,2,3,4,5].map(k=>{const r=G.plateDisp(G.D,G.slice(k),300,300+k);return {k,plates:r.plates,index:+r.obs.toFixed(3),simMed:+r.simMed.toFixed(3),sim95:+r.sim95.toFixed(3),simMax:+r.simMax.toFixed(3),quint:r.quint.map(q=>[+(q[0]/q[1]).toFixed(3),q[2],+q[3].toFixed(0)]),kappa:+r.kappa.kappa.toFixed(1)};});};
return G;
})();
