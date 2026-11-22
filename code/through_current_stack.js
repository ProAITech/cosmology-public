/* through_current_stack.js — Workbench Entry 57
   Does the rotation measure of background radio sources behind the halos of nearby inclined galaxies carry the
   pattern that a current running THROUGH each galaxy along its spin axis would give (Peratt's filament picture)?
   Such a current makes rings of field around the axis. Seen from outside, the rotation measure then has one sign
   on one side of the minor axis and the other sign on the other side, the same above and below the disc
   ("two lobes"). Entry 19 tested the four-quadrant pattern of currents running OUT along both poles instead.

   Signal. A line current I along the unit vector n through gas of density n_e(r):
       RM(x,y) = -0.812 (mu0 I / 2 pi) sin(i) x  Int n_e(r) dz / (x^2 + (y cos i - z sin i)^2),   r^2 = x^2+y^2+z^2
   x along the major axis, y along the minor axis (the side the current's projection points to), i the inclination.
   In uniform gas the integral is pi / (|x| sin i): RM = -/+ 0.812 n_e mu0 I / 2 whatever the distance, i.e.
   1.65 rad/m2 for 1e18 A in 1e-5 cm^-3. With halo gas n_e = N0 (r/kpc)^-1.5 (Miller & Bregman 2015,
   N0 = 1.35e-2 cm^-3) the pattern is RM = 710.5 rad/m2 x (b/kpc)^-1.5 x F(theta, i) per 1e18 A, b the impact
   parameter, theta the angle from the major axis; F = 1.748 on the major axis of an edge-on galaxy and
   pi sin^1.5(i) next to the minor axis.

   The catalogue gives no spin sense, so the current's direction is unknown galaxy by galaxy and a plain stack
   cancels. Two sign-free statistics are used instead:
     pairs:    sum over galaxies of sum_{i<j} t_i t_j r_i r_j  /  sum t_i^2 t_j^2      -> estimates I^2
     scatter:  sum t_i^2 (r_i^2 - s^2) / sum t_i^4                                     -> estimates I^2
   (t the template per 1e18 A, r the residual rotation measure, s^2 its mean square over the sight lines used).
   Null distributions: the same sight lines with every galaxy's position angle replaced by a random one.

   Runs in a browser tab. Inputs in window.D as for halo_sign_stack.js (which must be loaded first: HSS.galaxies).
   NVSS residuals are RM minus the median of neighbours within 3 degrees; the template is put through the same
   filter (minus the mean of the neighbours' templates, which is what a median does to a weak signal). */
window.TCS=(function(){
  const DEG=Math.PI/180,num=s=>{const v=parseFloat(s);return isFinite(v)?v:NaN;};
  const KRM=0.812*6.4815e4;      /* rad/m2 per cm^-3 for 1e18 A: 0.812 x (mu0 I / 2 pi in microgauss parsec) */
  let N0=1.35e-2,BCORE=10;       /* halo gas normalisation [cm^-3 at 1 kpc]; template held flat inside BCORE kpc */
  function Ftab(iDeg){const si=Math.sin(iDeg*DEG),ci=Math.cos(iDeg*DEG),tab=new Float64Array(180);
    for(let a=0;a<180;a++){const th=(a+0.5)*0.5*DEG,ct=Math.cos(th),st=Math.sin(th),n=a>=170?40000:4000,h=Math.PI/n;let s=0;
      for(let k=0;k<n;k++){const ph=-Math.PI/2+(k+0.5)*h,cp=Math.cos(ph),sp=Math.sin(ph),q=st*ci*cp-sp*si;s+=Math.pow(cp,1.5)*ct*si/(ct*ct*cp*cp+q*q);}
      tab[a]=s*h;}
    return tab;}
  const tabs={};
  function tabFor(i){const key=Math.max(1,Math.min(90,Math.round(i)));return tabs[key]||(tabs[key]=Ftab(key));}
  function Fval(tab,x,y){const th=Math.atan2(Math.abs(y),Math.abs(x))/DEG;return (x>=0?1:-1)*tab[Math.min(179,Math.floor(th*2))];}
  function tXY(tab,x,y){const b=Math.max(BCORE,Math.hypot(x,y));return KRM*N0*Math.pow(b,-1.5)*Fval(tab,x,y);}
  function rot(PA,e,n){const s=Math.sin(PA*DEG),c=Math.cos(PA*DEG);return [e*s+n*c,e*c-n*s];}   /* east,north -> major,minor */
  function off(g,ra,de){const dd=de-g.de;let da=ra-g.ra;if(da>180)da-=360;if(da<-180)da+=360;const sc=DEG*g.D*1000;return [da*Math.cos(g.de*DEG)*sc,dd*sc];}

  /* NVSS sources away from the Galactic plane, sorted by declination, with neighbour lists */
  function nvss(o){o=Object.assign({rf:3,bmin:20},o||{});
    const S=window.D.nvss.map(r=>({ra:num(r._RAJ2000),de:num(r._DEJ2000),gl:num(r._Glon),gb:num(r._Glat),rm:num(r.RM),e:num(r.e_RM)})).filter(s=>isFinite(s.rm)&&Math.abs(s.gb)>=o.bmin);
    S.sort((a,b)=>a.de-b.de);const de=S.map(s=>s.de);
    const lower=v=>{let lo=0,hi=de.length;while(lo<hi){const m=(lo+hi)>>1;if(de[m]<v)lo=m+1;else hi=m;}return lo;};
    const nb=new Array(S.length);
    for(let k=0;k<S.length;k++){const s=S[k],c=Math.cos(s.de*DEG),v=[];
      for(let j=lower(s.de-o.rf);j<S.length&&de[j]<=s.de+o.rf;j++){if(j===k)continue;let da=Math.abs(S[j].ra-s.ra);if(da>180)da=360-da;da*=c;if(da>o.rf)continue;if(da*da+(de[j]-s.de)**2<o.rf*o.rf)v.push(j);}
      nb[k]=Int32Array.from(v);}
    return {S,nb,de,lower,rf:o.rf};}
  function resid(N,rm){const out=new Float64Array(N.S.length);
    for(let k=0;k<N.S.length;k++){const nb=N.nb[k];if(nb.length<8){out[k]=NaN;continue;}const v=new Array(nb.length);for(let j=0;j<nb.length;j++)v[j]=rm?rm[nb[j]]:N.S[nb[j]].rm;v.sort((a,b)=>a-b);out[k]=(rm?rm[k]:N.S[k].rm)-v[v.length>>1];}
    return out;}

  /* Candidate galaxy-sightline pairs out to bout kpc, with offsets [kpc east, kpc north]; NVSS ones carry their neighbours' offsets */
  function candN(G,N,R,bout){const C=[];G.forEach((g,gi)=>{const tmax=bout/(g.D*1000)/DEG;
    for(let k=N.lower(g.de-tmax);k<N.S.length&&N.de[k]<=g.de+tmax;k++){if(!isFinite(R[k]))continue;const s=N.S[k],[e,n]=off(g,s.ra,s.de),b=Math.hypot(e,n);if(b>bout||b<DEG*g.D*1000/60)continue;
      const nb=N.nb[k],ne=new Float32Array(nb.length),nn=new Float32Array(nb.length);for(let j=0;j<nb.length;j++){const o=off(g,N.S[nb[j]].ra,N.S[nb[j]].de);ne[j]=o[0];nn[j]=o[1];}
      C.push({gi,k,e,n,b,ne,nn});}});return C;}
  function candL(G,bout){const S=window.D.lot.map(r=>({ra:num(r._RAJ2000),de:num(r._DEJ2000),rrm:num(r.RRM),z:num(r.zbest)})).filter(s=>isFinite(s.rrm));const C=[];
    G.forEach((g,gi)=>{const tmax=bout/(g.D*1000)/DEG;S.forEach((s,k)=>{if(Math.abs(s.de-g.de)>tmax)return;const [e,n]=off(g,s.ra,s.de),b=Math.hypot(e,n);if(b>bout||b<DEG*g.D*1000/60)return;if(isFinite(s.z)&&s.z<g.z+0.003)return;C.push({gi,k,e,n,b,r:s.rrm});});});return {C,S};}
  /* Template at a candidate for position angle PA: raw, and after the neighbour filter (NVSS) */
  function tOf(c,g,PA,filt){const tab=tabFor(g.i);const [x,y]=rot(PA,c.e,c.n);let t=tXY(tab,x,y);
    if(filt&&c.ne){let m=0;for(let j=0;j<c.ne.length;j++){const q=rot(PA,c.ne[j],c.nn[j]);m+=tXY(tab,q[0],q[1]);}t-=m/c.ne.length;}
    return {t,x,y};}
  function items(C,G,r,paFn,o){o=Object.assign({bin:30,bout:200,cut:50,filt:true},o||{});const out=[];
    for(const c of C){if(c.b<o.bin||c.b>=o.bout)continue;const rr=r?r[c.k]:c.r;if(!isFinite(rr)||Math.abs(rr)>=o.cut)continue;const g=G[c.gi],q=tOf(c,g,paFn?paFn(g,c.gi):g.PA,o.filt);out.push({gi:c.gi,tp:q.t,x:q.x,y:q.y,r:rr,b:c.b});}
    return out;}
  function stat(I){const by=new Map();let s2=0;for(const it of I)s2+=it.r*it.r;s2/=I.length;
    for(const it of I){let a=by.get(it.gi);if(!a){a={A:0,W:0,Q:0,W4:0,n:0,L:[]};by.set(it.gi,a);}a.A+=it.tp*it.r;a.W+=it.tp*it.tp;a.Q+=it.tp*it.tp*it.r*it.r;a.W4+=it.tp**4;a.n++;a.L.push(it);}
    let P=0,NP=0,V=0,NV=0,S2=0,ng=0,ng2=0,npair=0,so=0,no=0,ss=0,ns=0;
    for(const a of by.values()){P+=(a.A*a.A-a.Q)/2;NP+=(a.W*a.W-a.W4)/2;V+=a.Q;NV+=a.W4;S2+=a.W;ng++;if(a.n>1){ng2++;npair+=a.n*(a.n-1)/2;
      for(let i=0;i<a.n;i++)for(let j=i+1;j<a.n;j++){const p=a.L[i].r*a.L[j].r;if(a.L[i].x*a.L[j].x<0){so+=p;no++;}else{ss+=p;ns++;}}}}
    return {n:I.length,ng,ng2,npair,rms:Math.sqrt(s2),I2pair:P/NP,I2var:(V-s2*S2)/NV,sumT2:S2,nOpp:no,mOpp:so/no,nSame:ns,mSame:ss/ns,split:ss/ns-so/no};}
  function rng(seed){let s=seed||12345;return ()=>{s=(s*16807)%2147483647;return s/2147483647;};}
  function nulls(C,G,r,o,nreal,seed){const rnd=rng(seed),keys=['I2pair','I2var','split'],acc={};keys.forEach(k=>acc[k]=[]);
    for(let q=0;q<nreal;q++){const pa=G.map(()=>rnd()*180),s=stat(items(C,G,r,(g,gi)=>pa[gi],o));keys.forEach(k=>{if(isFinite(s[k]))acc[k].push(s[k]);});}
    const out={};keys.forEach(k=>{const a=acc[k],m=a.reduce((x,y)=>x+y,0)/a.length,sd=Math.sqrt(a.reduce((x,y)=>x+(y-m)**2,0)/a.length),srt=a.slice().sort((x,y)=>x-y);out[k]={mean:m,sd,p95:srt[Math.floor(0.95*srt.length)],n:a.length,all:a};});return out;}
  /* Inject a through-current of I18 x 1e18 A with a random sense per galaxy into the raw NVSS RMs, then re-filter */
  function injectN(G,N,I18,seed){const rnd=rng(seed||777),rm=Float64Array.from(N.S.map(s=>s.rm));
    G.forEach(g=>{const sg=rnd()<0.5?-1:1,tab=tabFor(g.i),tmax=Math.min(40,800/(g.D*1000)/DEG+N.rf+0.3);
      for(let k=N.lower(g.de-tmax);k<N.S.length&&N.de[k]<=g.de+tmax;k++){const s=N.S[k],[e,n]=off(g,s.ra,s.de);if(Math.hypot(e,n)>tmax*DEG*g.D*1000)continue;const [x,y]=rot(g.PA,e,n);rm[k]+=-sg*I18*tXY(tab,x,y);}});
    return resid(N,rm);}
  function set(o){if(o.N0!==undefined)N0=o.N0;if(o.BCORE!==undefined)BCORE=o.BCORE;return {N0,BCORE};}
  return {Ftab,tabFor,Fval,tXY,rot,off,nvss,resid,candN,candL,tOf,items,stat,nulls,injectN,rng,set,KRM};
})();

/* ---- Additions used for the published numbers ---------------------------------------------------------- */
(function(){const TCS=window.TCS,DEG=Math.PI/180;
  /* items with their sky offsets kept, so pairs can be formed explicitly */
  TCS.items2=function(C,G,r,paFn,o){o=Object.assign({bin:30,bout:200,cut:50,filt:true},o||{});const out=[];
    for(const c of C){if(c.b<o.bin||c.b>=o.bout)continue;const rr=r?r[c.k]:c.r;if(!isFinite(rr)||Math.abs(rr)>=o.cut)continue;const g=G[c.gi],q=TCS.tOf(c,g,paFn?paFn(g,c.gi):g.PA,o.filt);out.push({gi:c.gi,tp:q.t,x:q.x,y:q.y,r:rr,b:c.b,e:c.e,n:c.n});}return out;};
  /* pair statistic; minsep [deg] drops pairs of sight lines closer than that on the sky (it made no difference) */
  TCS.pairStat=function(I,G,minsep){const by=new Map();for(const it of I){let a=by.get(it.gi);if(!a){a=[];by.set(it.gi,a);}a.push(it);}
    let P=0,NP=0,np=0;const per=[];for(const [gi,L] of by){const sc=DEG*G[gi].D*1000*minsep;let p=0,d=0;for(let i=0;i<L.length;i++)for(let j=i+1;j<L.length;j++){if(Math.hypot(L[i].e-L[j].e,L[i].n-L[j].n)<sc)continue;const w=L[i].tp*L[j].tp;p+=w*L[i].r*L[j].r;d+=w*w;np++;}P+=p;NP+=d;per.push({gi,n:L.length,p,d});}
    return {I2:P/NP,den:NP,np,per};};
  TCS.pairNull=function(C,G,r,o,minsep,nreal,seed){const rnd=TCS.rng(seed),a=[];for(let k=0;k<nreal;k++){const pa=G.map(()=>rnd()*180);a.push(TCS.pairStat(TCS.items2(C,G,r,(g,gi)=>pa[gi],o),G,minsep).I2);}return a;};
  /* 95% limit: P(statistic <= observed | current I). Random position angles, a random sense per galaxy,
     signal = recovery x filtered template added to the residuals (recovery from injectN: 0.85 for NVSS, 1 for LoTSS). */
  TCS.limit=function(C,G,r,o,rho,Igrid,nreal,seed){const xobs=TCS.pairStat(TCS.items2(C,G,r,null,o),G,0).I2;const rnd=TCS.rng(seed),wide=Object.assign({},o,{cut:1e9}),res=[],draws=[];
    for(let k=0;k<nreal;k++){const pa=G.map(()=>rnd()*180),sg=G.map(()=>rnd()<0.5?-1:1);draws.push({I:TCS.items2(C,G,r,(g,gi)=>pa[gi],wide),sg});}
    for(const I18 of Igrid){let below=0;for(const d of draws){const J=[];for(const it of d.I){const rr=it.r-d.sg[it.gi]*I18*rho*it.tp;if(Math.abs(rr)<o.cut)J.push(Object.assign({},it,{r:rr}));}if(TCS.pairStat(J,G,0).I2<=xobs)below++;}res.push([I18,below/nreal]);}
    return {xobs,p:res};};
  /* One galaxy with a known spin (M31): raw RM = smooth foreground (polynomial of given order in east, north) - A x template */
  TCS.solve=function(A,b){const n=b.length,M=A.map((r,i)=>r.concat([b[i]]));for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(M[r][c])>Math.abs(M[p][c]))p=r;[M[c],M[p]]=[M[p],M[c]];for(let r=0;r<n;r++){if(r===c)continue;const q=M[r][c]/M[c][c];for(let k=c;k<=n;k++)M[r][k]-=q*M[c][k];}}return M.map((r,i)=>r[n]/r[i]);};
  TCS.inv=function(A){const n=A.length,M=A.map((r,i)=>r.concat(r.map((_,j)=>i===j?1:0)));for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(M[r][c])>Math.abs(M[p][c]))p=r;[M[c],M[p]]=[M[p],M[c]];const d=M[c][c];for(let k=0;k<2*n;k++)M[c][k]/=d;for(let r=0;r<n;r++){if(r===c)continue;const q=M[r][c];for(let k=0;k<2*n;k++)M[r][k]-=q*M[c][k];}}return M.map(r=>r.slice(n));};
  TCS.polyFit=function(src,g,PA,o){o=Object.assign({bin:30,bout:150,order:2,clip:3},o||{});const tab=TCS.tabFor(g.i),P=[];
    for(const s of src){const [e,n]=TCS.off(g,s.ra,s.de),b=Math.hypot(e,n);if(b<o.bin||b>=o.bout)continue;const [x,y]=TCS.rot(PA,e,n),u=e/100,v=n/100,row=[];for(let a=0;a<=o.order;a++)for(let c=0;c<=o.order-a;c++)row.push(Math.pow(u,a)*Math.pow(v,c));row.push(-TCS.tXY(tab,x,y));P.push({row,r:s.rm});}
    if(P.length<25)return {n:0,nAll:P.length,A:NaN,err:NaN};
    let keep=P,sol,cov,rms;for(let it=0;it<6;it++){const m=keep[0].row.length,ATA=Array.from({length:m},()=>new Array(m).fill(0)),ATb=new Array(m).fill(0);for(const p of keep)for(let i=0;i<m;i++){ATb[i]+=p.row[i]*p.r;for(let j=0;j<m;j++)ATA[i][j]+=p.row[i]*p.row[j];}
      sol=TCS.solve(ATA,ATb);cov=TCS.inv(ATA);let s2=0;for(const p of keep){p.res=p.r-p.row.reduce((a,v,i)=>a+v*sol[i],0);s2+=p.res*p.res;}rms=Math.sqrt(s2/(keep.length-m));const k2=keep.filter(p=>Math.abs(p.res)<o.clip*rms);if(k2.length===keep.length)break;keep=k2;}
    const m=sol.length;return {n:keep.length,nAll:P.length,A:sol[m-1],err:rms*Math.sqrt(cov[m-1][m-1]),rms};};
  /* Galactic -> equatorial (J2000), for placing the M31 template on empty fields at the same latitude */
  window.gal2eq=function(l,b){const ag=192.85948*DEG,dg=27.12825*DEG,ln=122.93192*DEG;l*=DEG;b*=DEG;const de=Math.asin(Math.sin(b)*Math.sin(dg)+Math.cos(b)*Math.cos(dg)*Math.cos(ln-l));const ra=ag+Math.atan2(Math.cos(b)*Math.sin(ln-l),Math.sin(b)*Math.cos(dg)-Math.cos(b)*Math.sin(dg)*Math.cos(ln-l));return [((ra/DEG)%360+360)%360,de/DEG];};
})();

/* Published run (10 November 2026), currents in units of 1e18 A for N0 = 1.35e-2 ("thin halo"); multiply by 0.367
   for the halo pinned to the density at the Large Magellanic Cloud (N0 = 3.68e-2).
   G = HSS.galaxies()                                   228 galaxies, median distance 17 Mpc
   NV = TCS.nvss(); RR = TCS.resid(NV); CN = TCS.candN(G,NV,RR,300)       28,127 residuals
   LL = TCS.candL(G,300)                                                    LoTSS, catalogue residuals (RRM)
   all galaxies, 30-200 kpc:
     NVSS   o={bin:30,bout:200,cut:50,filt:true}   1,311 sight lines, 30,621 pairs: I2 = 0.29; random angles 0.86 +- 0.73
            TCS.limit(CN,G,RR,o,0.85,grid,300,31337)  -> P = 0.073 at 1.0, 0.027 at 1.25: below 1.1
     LoTSS  o={bin:30,bout:200,cut:30,filt:false}  241 sight lines, 2,344 pairs: I2 = 0.59; random angles 1.06 +- 0.69
            TCS.limit(LL.C,G,null,o,1,grid,300,31337) -> P = 0.073 at 1.25, 0.017 at 1.5: below 1.35
     85% of the NVSS weight sits in eight galaxies 2-3 Mpc away, seven of them dwarfs; 99% of the LoTSS weight in five dwarfs.
   large spirals only (blue magnitude below -19.5, 118 galaxies):
     NVSS (cut 100) I2 = 6.4, random angles median 9.0, P still 0.10 at 5: no limit below 5.  LoTSS I2 = 3.7, median 3.1: below 4.7.
   injection (TCS.injectN, 6 senses each): 1, 2, 3 put in -> sqrt(excess I2) = 0.82, 1.69, 2.62 (NVSS); 1.01, 1.92, 2.81 (LoTSS).
   M31 = {ra:10.6847, de:41.2688, PA:38, i:77, D:0.77}   PA 38 is the receding (north-east) end, so +A is a current along the spin
     NVSS, |b| >= 5 (TCS.nvss({bmin:5})), 3-degree median filter, 30-150 kpc: 472 sight lines, A = -0.71; +-2 injected -> +1.53, -1.52
       (76% returned); the same template on 30 other fields at the same Galactic latitude x 6 angles: scatter 1.00  ->  -0.93 +- 1.3
     LoTSS raw RM, TCS.polyFit(lot,M31,38,{bin:30,bout:150,order:2}): 62 sight lines, A = -0.50 +- 0.85 (formal);
       77 other placements in the survey: scatter 1.0; rotating the template at M31: 1.7; injection returned in full.
     NVSS raw RM with polynomial foregrounds is worse: A = -2.35 (order 2; empty fields scatter 2.65), -0.85 (order 3; 1.57). */
