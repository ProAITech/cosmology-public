/* possum_through.js — Workbench Entry 61
   The through-current test of Entries 57 and 60, galaxy by galaxy, on the first public POSSUM rotation measures
   (Anderson et al. 2024, VizieR J/MNRAS/533/4068: 22,817 measurements in 35 ASKAP fields, about 30 per square degree,
   typical error 1 rad/m2). Runs in a browser tab on the VizieR site, with the catalogues held in memory only:
     D.possum  J/MNRAS/533/4068/table   RAJ2000, DEJ2000, SBID, peakRM, e_peakRM, peakRMc, peakRMfor, polsnr
     D.ungc    J/AJ/145/101/catalog     Karachentsev et al. 2013 (distances, inclinations)
     D.rc3     VII/155/rc3              photometric position angles
     D.lv2, D.lv9  J/MNRAS/478/1611     LVHIS, Koribalski et al. 2018 (table 9: rotation-curve position angles)
     D.phangs  J/ApJ/897/122/table1     Lang et al. 2020
   Model for one galaxy: raw RM = smooth foreground (polynomial in east, north) - A x template(x, y).
   +x points to the receding end of the major axis, so +A is a current along the spin (as in Entry 60).
   Templates per 1e18 A:
     halo   0.812 x 6.4815e4 x N0 x b^-1.5 x F(theta, i)   halo gas n_e = N0 (r/kpc)^-1.5, N0 = 3.68e-2 ("standard halo")
     step   1.65 x sign(x)                                  uniform gas of 1e-5 cm^-3 (the RM does not depend on distance)
   Errors: the same template fitted at other places in the survey where no galaxy is (PT.places), and turned
   through other angles at the galaxy itself (PT.turn). */
window.PT=(function(){
  const DEG=Math.PI/180,KRM=0.812*6.4815e4;let N0=3.68e-2,BCORE=10;
  function Ftab(iDeg){const si=Math.sin(iDeg*DEG),ci=Math.cos(iDeg*DEG),tab=new Float64Array(180);
    for(let a=0;a<180;a++){const th=(a+0.5)*0.5*DEG,ct=Math.cos(th),st=Math.sin(th),n=a>=170?40000:4000,h=Math.PI/n;let s=0;
      for(let k=0;k<n;k++){const ph=-Math.PI/2+(k+0.5)*h,cp=Math.cos(ph),sp=Math.sin(ph),q=st*ci*cp-sp*si;s+=Math.pow(cp,1.5)*ct*si/(ct*ct*cp*cp+q*q);}
      tab[a]=s*h;}
    return tab;}
  const tabs={};
  const tabFor=i=>{const key=Math.max(1,Math.min(90,Math.round(i)));return tabs[key]||(tabs[key]=Ftab(key));};
  const Fval=(tab,x,y)=>{const th=Math.atan2(Math.abs(y),Math.abs(x))/DEG;return (x>=0?1:-1)*tab[Math.min(179,Math.floor(th*2))];};
  const tHalo=(tab,x,y)=>{const b=Math.max(BCORE,Math.hypot(x,y));return KRM*N0*Math.pow(b,-1.5)*Fval(tab,x,y);};
  const tStep=(x)=>1.65*(x>=0?1:-1);
  const rot=(PA,e,n)=>{const s=Math.sin(PA*DEG),c=Math.cos(PA*DEG);return [e*s+n*c,e*c-n*s];};
  /* sources around a point, offsets in kpc east and north at distance D [Mpc] */
  function around(ra,de,D,rmaxKpc){const rd=rmaxKpc/(D*1000)/DEG,sc=DEG*D*1000;return window.psNear(ra,de,rd).map(s=>({k:s.k,e:s.e*sc,n:s.n*sc,rm:window.PS[s.k].rm}));}
  function solve(A,b){const n=b.length,M=A.map((r,i)=>r.concat([b[i]]));for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(M[r][c])>Math.abs(M[p][c]))p=r;[M[c],M[p]]=[M[p],M[c]];for(let r=0;r<n;r++){if(r===c)continue;const q=M[r][c]/M[c][c];for(let k=c;k<=n;k++)M[r][k]-=q*M[c][k];}}return M.map((r,i)=>r[n]/r[i]);}
  function inv(A){const n=A.length,M=A.map((r,i)=>r.concat(r.map((_,j)=>i===j?1:0)));for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(M[r][c])>Math.abs(M[p][c]))p=r;[M[c],M[p]]=[M[p],M[c]];const d=M[c][c];for(let k=0;k<2*n;k++)M[c][k]/=d;for(let r=0;r<n;r++){if(r===c)continue;const q=M[r][c];for(let k=0;k<2*n;k++)M[r][k]-=q*M[c][k];}}return M.map(r=>r.slice(n));}
  /* joint fit. S from around(); o.inj adds a current of that many 1e18 A (halo template) to the data first */
  function fit(S,i,PA,o){o=Object.assign({bin:30,bout:200,order:2,clip:3,tmpl:'halo',inj:0,scale:100},o||{});const tab=tabFor(i),P=[];
    for(const s of S){const b=Math.hypot(s.e,s.n);if(b<o.bin||b>=o.bout)continue;const [x,y]=rot(PA,s.e,s.n),u=s.e/o.scale,v=s.n/o.scale,row=[];
      for(let a=0;a<=o.order;a++)for(let c=0;c<=o.order-a;c++)row.push(Math.pow(u,a)*Math.pow(v,c));
      const t=o.tmpl==='step'?tStep(x):tHalo(tab,x,y);row.push(-t);
      let r=s.rm;if(o.inj)r+=-o.inj*(o.injT==='step'?tStep(x):tHalo(tab,x,y));P.push({row,r,x,y,b,t,k:s.k});}
    const m=P.length?P[0].row.length:0;if(P.length<m+12)return {n:0,nAll:P.length,A:NaN,err:NaN};
    /* first pass: drop wild values against the median */
    const med=P.map(p=>p.r).sort((a,b)=>a-b)[P.length>>1],mad=1.4826*P.map(p=>Math.abs(p.r-med)).sort((a,b)=>a-b)[P.length>>1];
    let keep=P.filter(p=>Math.abs(p.r-med)<6*mad),sol,cov,rms;
    for(let it=0;it<8;it++){const ATA=Array.from({length:m},()=>new Array(m).fill(0)),ATb=new Array(m).fill(0);
      for(const p of keep)for(let a=0;a<m;a++){ATb[a]+=p.row[a]*p.r;for(let c=0;c<m;c++)ATA[a][c]+=p.row[a]*p.row[c];}
      sol=solve(ATA,ATb);cov=inv(ATA);let s2=0;for(const p of keep){p.fg=0;for(let a=0;a<m-1;a++)p.fg+=p.row[a]*sol[a];p.res=p.r-p.fg-p.row[m-1]*sol[m-1];s2+=p.res*p.res;}
      rms=Math.sqrt(s2/(keep.length-m));const k2=keep.filter(p=>Math.abs(p.res)<o.clip*rms);if(k2.length===keep.length)break;keep=k2;}
    return {n:keep.length,nAll:P.length,A:sol[m-1],err:rms*Math.sqrt(cov[m-1][m-1]),rms,pts:o.keep?keep:undefined};}
  /* the template turned through other angles at the galaxy */
  function turn(S,i,o,nang){const a=[];for(let k=0;k<nang;k++){const f=fit(S,i,k*360/nang,o);if(isFinite(f.A))a.push(f.A);}return a;}
  /* the template at other places: centres on a grid over the survey, kept if the ring holds at least frac of the galaxy's count
     and no catalogued galaxy nearer than 12 Mpc lies within the ring */
  function places(g,o,nref,frac,step,nang,seed){let s=seed||4242;const rnd=()=>{s=(s*16807)%2147483647;return s/2147483647;};
    const rd=(o.bout||200)/(g.D*1000)/DEG,out=[];const PSs=window.PS;let dmin=90,dmax=-90;for(const p of PSs){if(p.de<dmin)dmin=p.de;if(p.de>dmax)dmax=p.de;}
    for(let de=dmin;de<=dmax;de+=step*rd){const c=Math.cos(de*DEG);for(let ra=0;ra<360;ra+=step*rd/c){
      if(Math.hypot((ra-g.ra)*Math.cos(g.de*DEG),de-g.de)<1.5*rd)continue;
      const S=around(ra,de,g.D,o.bout||200),nin=S.filter(q=>Math.hypot(q.e,q.n)>=(o.bin||30)).length;if(nin<frac*nref)continue;
      const qd=[0,0,0,0];S.forEach(q=>{qd[(q.e>0?1:0)+(q.n>0?2:0)]++;});if(Math.min(...qd)<0.1*nin)continue;
      for(let k=0;k<nang;k++){const f=fit(S,g.i,rnd()*360,o);if(isFinite(f.A))out.push({ra,de,A:f.A,err:f.err,n:f.n});}}}
    return out;}
  const sd=a=>{const m=a.reduce((x,y)=>x+y,0)/a.length;return Math.sqrt(a.reduce((x,y)=>x+(y-m)**2,0)/a.length);};
  const rob=a=>{const s=a.slice().sort((x,y)=>x-y),m=s[s.length>>1];return 1.4826*a.map(v=>Math.abs(v-m)).sort((x,y)=>x-y)[a.length>>1];};
  function set(o){if(o.N0!==undefined)N0=o.N0;if(o.BCORE!==undefined)BCORE=o.BCORE;return {N0,BCORE};}
  return {Ftab,tabFor,Fval,tHalo,tStep,rot,around,fit,turn,places,sd,rob,set,KRM};
})();

/* ---- Set-up in the tab, and the additions used for the published numbers --------------------------------- */
/* PS: the catalogue with repeat measurements merged (same position within 8 arcsec in two scheduling blocks: 2,213 pairs,
   differences scatter by 2.4 rad/m2), sorted by declination; psNear(ra, dec, radius_deg) returns the sources in a circle.
     raw = D.possum.rows -> {ra, de, sb, rm, e, rmc, fg}; merged with weights 1/e^2  -> 20,442 sources, median error 1.1 rad/m2 */
(function(){const PT=window.PT,DEG=Math.PI/180;
  /* how evenly the sight lines surround a galaxy: smallest share among the four quadrants about its axes */
  PT.cover=function(g,o){const S=PT.around(g.ra,g.de,g.D,o.bout).filter(q=>{const b=Math.hypot(q.e,q.n);return b>=o.bin&&b<o.bout;});const qd=[0,0,0,0];
    S.forEach(q=>{const [x,y]=PT.rot(g.PA,q.e,q.n);qd[(x>0?1:0)+(y>0?2:0)]++;});return {n:S.length,qmin:Math.min(...qd)/Math.max(1,S.length)};};
  /* fit at the galaxy plus the blank-sky set */
  PT.full=function(g,o,frac,step,nang,seed){const S=PT.around(g.ra,g.de,g.D,o.bout*1.02),f=PT.fit(S,g.i,g.PA,o),pl=PT.places(g,o,f.nAll,frac,step,nang,seed);
    return {f,a:pl.map(p=>p.A),npos:new Set(pl.map(p=>p.ra.toFixed(3)+'_'+p.de.toFixed(3))).size};};
  /* a ring's dipole: RM = polynomial over [bin, bout] + (c1 cos phi + c2 sin phi) for sources in [r1, r2]; phi north through east */
  PT.dipole=function(S,o){o=Object.assign({bin:15,bout:150,order:2,r1:15,r2:30,clip:3,scale:100},o||{});const P=[];
    for(const s of S){const b=Math.hypot(s.e,s.n);if(b<o.bin||b>=o.bout)continue;const u=s.e/o.scale,v=s.n/o.scale,row=[];
      for(let a=0;a<=o.order;a++)for(let c=0;c<=o.order-a;c++)row.push(Math.pow(u,a)*Math.pow(v,c));
      const inr=(b>=o.r1&&b<o.r2)?1:0,phi=Math.atan2(s.e,s.n);row.push(inr*Math.cos(phi),inr*Math.sin(phi));P.push({row,r:s.rm,b,inr});}
    const m=P.length?P[0].row.length:0,nr=P.filter(p=>p.inr).length;if(P.length<m+15||nr<8)return {n:0,nr};
    const med=P.map(p=>p.r).sort((a,b)=>a-b)[P.length>>1],mad=1.4826*P.map(p=>Math.abs(p.r-med)).sort((a,b)=>a-b)[P.length>>1];
    let keep=P.filter(p=>Math.abs(p.r-med)<6*mad),sol,rms;
    const solve=(A,b)=>{const n=b.length,M=A.map((r,i)=>r.concat([b[i]]));for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(M[r][c])>Math.abs(M[p][c]))p=r;[M[c],M[p]]=[M[p],M[c]];for(let r=0;r<n;r++){if(r===c)continue;const q=M[r][c]/M[c][c];for(let k=c;k<=n;k++)M[r][k]-=q*M[c][k];}}return M.map((r,i)=>r[n]/r[i]);};
    for(let it=0;it<8;it++){const ATA=Array.from({length:m},()=>new Array(m).fill(0)),ATb=new Array(m).fill(0);
      for(const p of keep)for(let a=0;a<m;a++){ATb[a]+=p.row[a]*p.r;for(let c=0;c<m;c++)ATA[a][c]+=p.row[a]*p.row[c];}
      sol=solve(ATA,ATb);let s2=0;for(const p of keep){p.res=p.r-p.row.reduce((a,v,i)=>a+v*sol[i],0);s2+=p.res*p.res;}
      rms=Math.sqrt(s2/(keep.length-m));const k2=keep.filter(p=>Math.abs(p.res)<o.clip*rms);if(k2.length===keep.length)break;keep=k2;}
    const c1=sol[m-2],c2=sol[m-1];return {n:keep.length,nr:keep.filter(p=>p.inr).length,amp:Math.hypot(c1,c2),dir:((Math.atan2(c2,c1)/DEG)+360)%360,rms};};
  /* the same ring on blank sky: centres on a grid of spacing step x the outer radius, at least 300 sources and no empty quadrant */
  PT.dipNull=function(g,o,step,minr){const rd=(o.bout||150)/(g.D*1000)/DEG,out=[];let dmin=90,dmax=-90;for(const p of window.PS){if(p.de<dmin)dmin=p.de;if(p.de>dmax)dmax=p.de;}
    for(let de=dmin;de<=dmax;de+=step*rd){const c=Math.cos(de*DEG);for(let ra=0;ra<360;ra+=step*rd/c){if(Math.hypot((ra-g.ra)*Math.cos(g.de*DEG),de-g.de)<1.2*rd)continue;
      const S=PT.around(ra,de,g.D,o.bout||150);if(S.length<300)continue;const qd=[0,0,0,0];S.forEach(q=>{qd[(q.e>0?1:0)+(q.n>0?2:0)]++;});if(Math.min(...qd)<0.12*S.length)continue;
      const d=PT.dipole(S,o);if(d.n&&d.nr>=minr)out.push({ra,de,amp:d.amp,dir:d.dir,nr:d.nr});}}
    return out;};
  /* weighted stack and its blank-sky distribution: each galaxy's value replaced by a draw from its own blank-sky set */
  PT.stackNull=function(R,seed,ndraw){const w=R.map(r=>1/(r.rob*r.rob)),W=w.reduce((a,b)=>a+b,0),obs=R.reduce((s,r,i)=>s+w[i]*r.A,0)/W;let s=seed;const rnd=()=>{s=(s*16807)%2147483647;return s/2147483647;};
    const n1=[];let allpos=0;for(let k=0;k<ndraw;k++){let m=0,np=0;for(let i=0;i<R.length;i++){const a=R[i].nulls[Math.floor(rnd()*R[i].nulls.length)];m+=w[i]*a;if(a>0)np++;}n1.push(m/W);if(np===R.length)allpos++;}
    n1.sort((a,b)=>a-b);return {obs,sd:PT.sd(n1),p:n1.filter(v=>v>=obs).length/ndraw,lo95:obs-n1[Math.floor(0.975*ndraw)],hi95:obs-n1[Math.floor(0.025*ndraw)],allpos:allpos/ndraw};};
})();

/* Published run (11 November 2026). Currents in 1e18 A for the standard halo; + = along the spin. o = {bin:30, bout:200, order:2}.
   Catalogue: 22,817 rows, 35 scheduling blocks, 20,442 sources after merging repeats; 31% have a neighbour within 2 arcmin
   (double sources), which is one reason formal errors are too small.
   Galaxies with a rotation-curve position angle (LVHIS table 9: the angle of the receding side; it agrees with THINGS for
   NGC 3621, NGC 5236 and NGC 7793 within 2 deg) or a PHANGS one, and at least 30 sight lines at 30-200 kpc:
     well covered (every quadrant holds at least 10% of the sight lines):
       name         D     i   PA   n     A       formal  blank-sky sd / robust   places   step A (robust)
       NGC 300      2.15  40  290  1212  +0.227  0.077   0.539 / 0.277            44      +2.12 (1.46)
       ESO 154-G023 5.55  80  218   327  +0.127  0.091   0.223 / 0.161           199      +0.79 (1.10)
       NGC 1311     5.20  71   41   319  +0.163  0.078   0.230 / 0.189           151      +0.29 (1.11)
       IC 1959      6.05  79  149   300  +0.110  0.089   0.243 / 0.180           249      +1.21 (1.12)
       IC 1954     15.2   58   61    56  +0.554  0.317   0.575 / 0.359           587      +1.31 (1.81)
       NGC 685     16.0   39   99    36  +1.540  0.897   0.973 / 0.513           702      +5.44 (1.95)     inclination below 40
       NGC 1672    11.9   24  149    77  +0.417  0.341   1.195 / 0.595           403      +0.76 (1.52)     inclination below 40
     poorly covered, not used: ESO 245-G005 +0.78 (smallest quadrant 4%), IC 5152 -0.89 (0%), NGC 1313 +0.2 +- 1.3 (0%);
       NGC 55 and ESO 115-G021 have no usable fit.
   PT.full(g, o, 0.5, step, 6, 1234) with step 0.45, 0.7 or 1.2 for outer radii above 3, above 1.5 and below 1.5 degrees.
   Stacks, weights 1/robust^2, PT.stackNull(R, 1357, 200000):
     five (the cuts of Entry 60):  +0.170, blank-sky sd 0.131, 95% -0.08 to +0.44, P(as large) 0.074, P(all positive) 0.031
     seven:                        +0.217, 0.130, -0.04 to +0.48, 0.040, 0.007      random angles at the true places: sd 0.116
     step pattern: five +1.01 +- 0.70, seven +1.28 +- 0.68 (x 3.3 rad/m2 = receding side minus approaching side, with a minus sign)
   Other settings, seven-galaxy stack: order 1 +0.155, order 3 +0.127, from 50 kpc +0.160, to 300 kpc +0.131, from 15 kpc +0.109,
     100-200 kpc only +0.654 (all seven positive); catalogue's own foreground-corrected RM with a constant +0.117, with a plane +0.179.
     Template turned by 90 deg: +0.28, +0.22, -0.15, +0.44, +1.09, -0.21, -0.86. NGC 300 with the outer disc's angle (332): +0.03.
     Leaving one galaxy out: +0.18 to +0.26. Injected +0.3, -0.3, +1: returned in full.
   Sense unknown (likelihood with either sign, robust errors): seven below 0.30; with six more galaxies on their optical axes
     (NGC 625 -0.21, NGC 1035 +0.51, IC 1914 -0.05, IC 1933 -0.29, NGC 1527 +0.13, NGC 7205 -0.07; robust 0.17 to 0.50) below 0.27.
     Step pattern, seven: below 1.8 (uniform gas of 1e-5 cm^-3).
   NGC 300, ring of 15-30 kpc (35 sources; fit over 15-150 kpc, order 2): dipole 7.79 rad/m2, largest toward position angle 135
     (order 3: 7.91, 128). Sources within 60 deg of the major axis: receding side -6.5 +- 1.2 (12), approaching +6.4 +- 1.9 (13).
     PT.dipNull(N300, {r1:15, r2:30, order:2, bin:15, bout:150}, 0.35, 20): 108 blank rings, median 4.2, 90th percentile 6.8,
     95th 8.7, 7 at or above 7.79. Ring 45-80 kpc: 4.98 toward 293; 23 of 107 blank rings are larger.
   Stuardi et al. 2026 filament (Abell 3651 to Abell 3667): this catalogue has 5 sources within 1 deg of the axis; no test possible. */
