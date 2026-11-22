/* pole_currents.js — Workbench Entry 64
   Entry 61 fitted, galaxy by galaxy, the pattern of a current running THROUGH a galaxy along its spin axis (the same ring sense
   above and below the disc). The Milky Way's own halo shows the other kind: rings of opposite sense above and below the disc,
   which is a current running OUT along both poles (or in along both) and closing through the disc. This script fits both
   patterns at once to the first POSSUM fields (see possum_through.js for the catalogue set-up: window.PS, window.psNear).
     raw RM = smooth foreground (polynomial) - A x t_through(x, y) - Q x t_poles(x, y)
   +x toward the receding end of the major axis, +y a quarter turn on from it (north through east), which is where the spin vector
   points on the sky. With h the height along the spin axis, h = y sin i + z cos i on a sight line (z along the line of sight):
     t_through = K N0 b^-1.5 F(theta, i) sign(x)                       F as in Entry 57
     t_poles   = K N0 b^-1.5 F_odd(theta, i) sign(x) sign(y)           F_odd: the same integral with the integrand times sign(h)
   For an edge-on galaxy the two F are equal. Turning the galaxy the other way up (z -> -z) leaves both unchanged, so the near
   side need not be known. +A: current along the spin. +Q: current out along both poles (ring field with the rotation on the
   side the spin points to, against it on the other), the sense the Milky Way has (Entry 18). Units 1e18 A, standard halo. */
window.PC=(function(){
  const DEG=Math.PI/180,KRM=0.812*6.4815e4;let N0=3.68e-2,BCORE=10;
  function Ftab(iDeg){const si=Math.sin(iDeg*DEG),ci=Math.cos(iDeg*DEG),ev=new Float64Array(180),od=new Float64Array(180);
    for(let a=0;a<180;a++){const th=(a+0.5)*0.5*DEG,ct=Math.cos(th),st=Math.sin(th),n=a>=170?40000:4000,h=Math.PI/n;let s=0,o=0;
      for(let k=0;k<n;k++){const ph=-Math.PI/2+(k+0.5)*h,cp=Math.cos(ph),sp=Math.sin(ph),q=st*ci*cp-sp*si,t=Math.pow(cp,1.5)*ct*si/(ct*ct*cp*cp+q*q);s+=t;o+=t*((st*si*cp+sp*ci)>=0?1:-1);}
      ev[a]=s*h;od[a]=o*h;}
    return {ev,od};}
  const tabs={};
  const tabFor=i=>{const key=Math.max(1,Math.min(90,Math.round(i)));return tabs[key]||(tabs[key]=Ftab(key));};
  const idx=(x,y)=>Math.min(179,Math.floor(Math.atan2(Math.abs(y),Math.abs(x))/DEG*2));
  const amp=(x,y)=>KRM*N0*Math.pow(Math.max(BCORE,Math.hypot(x,y)),-1.5);
  const tThru=(tab,x,y)=>amp(x,y)*(x>=0?1:-1)*tab.ev[idx(x,y)];
  const tPole=(tab,x,y)=>amp(x,y)*(x>=0?1:-1)*(y>=0?1:-1)*tab.od[idx(x,y)];
  const rot=(PA,e,n)=>{const s=Math.sin(PA*DEG),c=Math.cos(PA*DEG);return [e*s+n*c,e*c-n*s];};
  function around(ra,de,D,rmaxKpc){const rd=rmaxKpc/(D*1000)/DEG,sc=DEG*D*1000;return window.psNear(ra,de,rd).map(s=>({k:s.k,e:s.e*sc,n:s.n*sc,rm:window.PS[s.k].rm}));}
  function solve(A,b){const n=b.length,M=A.map((r,i)=>r.concat([b[i]]));for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(M[r][c])>Math.abs(M[p][c]))p=r;[M[c],M[p]]=[M[p],M[c]];for(let r=0;r<n;r++){if(r===c)continue;const q=M[r][c]/M[c][c];for(let k=c;k<=n;k++)M[r][k]-=q*M[c][k];}}return M.map((r,i)=>r[n]/r[i]);}
  function inv(A){const n=A.length,M=A.map((r,i)=>r.concat(r.map((_,j)=>i===j?1:0)));for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(M[r][c])>Math.abs(M[p][c]))p=r;[M[c],M[p]]=[M[p],M[c]];const d=M[c][c];for(let k=0;k<2*n;k++)M[c][k]/=d;for(let r=0;r<n;r++){if(r===c)continue;const q=M[r][c];for(let k=0;k<2*n;k++)M[r][k]-=q*M[c][k];}}return M.map(r=>r.slice(n));}
  /* joint fit of both patterns; o.injA, o.injQ add currents to the data first */
  function fit(S,i,PA,o){o=Object.assign({bin:30,bout:200,order:2,clip:3,scale:100,injA:0,injQ:0},o||{});const tab=tabFor(i),P=[];
    for(const s of S){const b=Math.hypot(s.e,s.n);if(b<o.bin||b>=o.bout)continue;const [x,y]=rot(PA,s.e,s.n),u=s.e/o.scale,v=s.n/o.scale,row=[];
      for(let a=0;a<=o.order;a++)for(let c=0;c<=o.order-a;c++)row.push(Math.pow(u,a)*Math.pow(v,c));
      const t1=tThru(tab,x,y),t2=tPole(tab,x,y);row.push(-t1,-t2);P.push({row,r:s.rm-o.injA*t1-o.injQ*t2,x,y,b});}
    const m=P.length?P[0].row.length:0;if(P.length<m+12)return {n:0,nAll:P.length,A:NaN,Q:NaN};
    const med=P.map(p=>p.r).sort((a,b)=>a-b)[P.length>>1],mad=1.4826*P.map(p=>Math.abs(p.r-med)).sort((a,b)=>a-b)[P.length>>1];
    let keep=P.filter(p=>Math.abs(p.r-med)<6*mad),sol,cov,rms;
    for(let it=0;it<8;it++){const ATA=Array.from({length:m},()=>new Array(m).fill(0)),ATb=new Array(m).fill(0);
      for(const p of keep)for(let a=0;a<m;a++){ATb[a]+=p.row[a]*p.r;for(let c=0;c<m;c++)ATA[a][c]+=p.row[a]*p.row[c];}
      sol=solve(ATA,ATb);cov=inv(ATA);let s2=0;for(const p of keep){p.res=p.r-p.row.reduce((a,v,k)=>a+v*sol[k],0);s2+=p.res*p.res;}
      rms=Math.sqrt(s2/(keep.length-m));const k2=keep.filter(p=>Math.abs(p.res)<o.clip*rms);if(k2.length===keep.length)break;keep=k2;}
    return {n:keep.length,nAll:P.length,A:sol[m-2],Q:sol[m-1],eA:rms*Math.sqrt(cov[m-2][m-2]),eQ:rms*Math.sqrt(cov[m-1][m-1]),rho:cov[m-2][m-1]/Math.sqrt(cov[m-2][m-2]*cov[m-1][m-1]),rms};}
  /* blank sky: the same two patterns at other places and random angles (as PT.places in possum_through.js) */
  function places(g,o,nref,frac,step,nang,seed){let s=seed||4242;const rnd=()=>{s=(s*16807)%2147483647;return s/2147483647;};
    const rd=(o.bout||200)/(g.D*1000)/DEG,out=[];let dmin=90,dmax=-90;for(const p of window.PS){if(p.de<dmin)dmin=p.de;if(p.de>dmax)dmax=p.de;}
    for(let de=dmin;de<=dmax;de+=step*rd){const c=Math.cos(de*DEG);for(let ra=0;ra<360;ra+=step*rd/c){
      if(Math.hypot((ra-g.ra)*Math.cos(g.de*DEG),de-g.de)<1.5*rd)continue;
      const S=around(ra,de,g.D,o.bout||200),nin=S.filter(q=>Math.hypot(q.e,q.n)>=(o.bin||30)).length;if(nin<frac*nref)continue;
      const qd=[0,0,0,0];S.forEach(q=>{qd[(q.e>0?1:0)+(q.n>0?2:0)]++;});if(Math.min(...qd)<0.1*nin)continue;
      for(let k=0;k<nang;k++){const f=fit(S,g.i,rnd()*360,o);if(isFinite(f.Q))out.push({A:f.A,Q:f.Q});}}}
    return out;}
  const sd=a=>{const m=a.reduce((x,y)=>x+y,0)/a.length;return Math.sqrt(a.reduce((x,y)=>x+(y-m)**2,0)/a.length);};
  const rob=a=>{const s=a.slice().sort((x,y)=>x-y),m=s[s.length>>1];return 1.4826*a.map(v=>Math.abs(v-m)).sort((x,y)=>x-y)[a.length>>1];};
  function stackNull(vals,robs,nulls,seed,ndraw){const w=robs.map(r=>1/(r*r)),W=w.reduce((a,b)=>a+b,0),obs=vals.reduce((s,v,i)=>s+w[i]*v,0)/W;let s=seed;const rnd=()=>{s=(s*16807)%2147483647;return s/2147483647;};
    const n1=[];for(let k=0;k<ndraw;k++){let m=0;for(let i=0;i<vals.length;i++)m+=w[i]*nulls[i][Math.floor(rnd()*nulls[i].length)];n1.push(m/W);}
    n1.sort((a,b)=>a-b);return {obs,sd:sd(n1),p_abs:n1.filter(v=>Math.abs(v)>=Math.abs(obs)).length/ndraw,lo95:obs-n1[Math.floor(0.975*ndraw)],hi95:obs-n1[Math.floor(0.025*ndraw)]};}
  return {Ftab,tabFor,tThru,tPole,rot,around,fit,places,sd,rob,stackNull,KRM};
})();

/* Published run (11 November 2026). o = {bin:30, bout:200, order:2}; blank sky as in Entry 61 (PC.places(g, o, n, 0.5, step, 6, seed)).
   The four-quadrant pattern sign(x) sign(y) does not change when the other end of the major axis is taken as +x, so the pole
   test needs only the axis and not the spin sense. Seven galaxies with a known spin (Entry 61) and five on their optical axes:
     name          D     i   PA    n     A (robust)       Q      formal  robust   places x 6
     NGC 300       2.15  40  290  1214  +0.214 (0.298)  +0.160   0.096   0.276     264
     ESO 154-G023  5.55  80  218   331  +0.113 (0.172)  -0.177   0.098   0.161    1194
     NGC 1311      5.20  71   41   317  +0.165 (0.194)  +0.035   0.079   0.171     906
     IC 1959       6.05  79  149   300  +0.095 (0.187)  +0.164   0.088   0.169    1494
     IC 1954      15.2   58   61    56  +0.473 (0.427)  +0.304   0.309   0.461    3522
     NGC 685      16.0   39   99    36  +1.143 (0.663)  +1.427   0.968   0.808    4194
     NGC 1672     11.9   24  149    77  +0.451 (0.665)  +0.190   0.521   1.016    2418
     NGC 1035     14.6   74  150    51                  -0.184   0.192   0.388    3372
     IC 1914      13.6   61   99    51                  +0.708   0.320   0.421    3054
     IC 1933      14.2   60   55    75                  +0.091   0.247   0.454    3096
     NGC 1527     13.7   71   78    63                  +0.210   0.121   0.357    3084
     NGC 7205     21.9   63   73    33                  -0.269   0.299   0.750    4914
     NGC 625 (3.89 Mpc, i 76, PA 92): Q = -2.5 with a formal error of 3.4 against a blank-sky 0.16; the pattern cannot be told
       from the foreground with its sight lines, and it is left out.
   Stacks of Q, weights 1/robust^2, PC.stackNull(..., 200000 draws):
     five with spin (Entry 60's cuts)  +0.029, blank-sky spread 0.142, 95% -0.27 to +0.32
     seven with spin                   +0.047, 0.141, -0.25 to +0.33
     twelve                            +0.067, 0.132, -0.21 to +0.33;  9 of 12 positive;  sum of (Q/robust)^2 = 9.7
     sign ignored (likelihood with either sign): |Q| below 0.24 (five), 0.245 (seven), 0.23 (twelve)
   Through-current A in the same joint fits: five +0.151 +- 0.144, seven +0.177 +- 0.142; all seven still positive.
   Injected Q = 0.5 returns +0.500 in Q and 0.000 in A; injected A = 0.5 the reverse.
   Pattern size per 1e18 A, 45 deg from the axes, inclination 80: 23, 8.2, 3.8, 2.1 rad/m2 at 30, 60, 100, 150 kpc.
   NGC 300, 35 sources at 15-30 kpc, harmonics about the receding axis (290) after a quadratic foreground over 15-150 kpc:
     cos(phi) -7.1 (the two-lobe term), sin(phi) -2.5, sin(2 phi) +1.3 (the four-quadrant term), cos(2 phi) +2.0 rad/m2;
     with the axis at 310: -7.6, +0.1, -0.3, +2.3; at 332: -7.0, +3.0, -1.8, +1.5.
     108 blank rings at random angles: spread of the cos(phi) term 3.3, of the sin(2 phi) term 3.2. */
