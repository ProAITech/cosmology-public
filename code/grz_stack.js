// Entry 150: gravitational-redshift stack around SDSS groups, run in the browser on a VizieR mirror page.
// Tempel et al. 2014 (VizieR J/A+A/566/A1): 82,458 groups and 588,193 galaxies (SDSS DR10 main sample, flux-limited, z < 0.2).
const c=299792.458, H0=67.7, Om=0.31, D2R=Math.PI/180;
const grpTxt=await (await fetch('/viz-bin/asu-tsv?-source=J/A%2BA/566/A1/groups&-out=GroupID,Ngal,RAJ2000,DEJ2000,zcmb,sig.v,Rvir,MNFW&-out.max=200000&-out.form=|')).text();
const GRP=grpTxt.split('\n').filter(x=>x&&x[0]!=='#').slice(3).map(x=>{const p=x.split('|');return {id:+p[0],n:+p[1],ra:+p[2],de:+p[3],z:+p[4],sv:+p[5],rv:+p[6],m:+p[7]}}).filter(g=>isFinite(g.id));   // m: NFW mass, 1e12 Msun
const galTxt=await (await fetch('/viz-bin/asu-tsv?-source=J/A%2BA/566/A1/galaxies&-out=GroupID,Ngal,Rank,zobs,e_z,RAJ2000,DEJ2000,rmag&-out.max=700000&-out.form=|')).text();
const GALS=galTxt.split('\n').filter(x=>x&&x[0]!=='#').slice(3).map(x=>{const p=x.split('|');return {gid:+p[0],ng:+p[1],rk:+p[2],z:+p[3],ez:+p[4],ra:+p[5],de:+p[6],rm:+p[7]}}).filter(o=>isFinite(o.z)&&isFinite(o.ra));
// (rMAG, the k+e-corrected absolute magnitude, was read the same way and used only to build the m - M relation below)
function Dc(z){let s=0,n=200;for(let i=0;i<n;i++){const zz=(i+0.5)/n*z;s+=1/Math.sqrt(Om*(1+zz)**3+1-Om);}return c/H0*s*z/n;}
const bcg=new Map(); for(const g of GALS) if(g.rk===1) bcg.set(g.gid,g);     // central galaxy = brightest member (Rank 1)
// spatial index, 0.5-deg declination strips
const CS=0.5, idx=new Map(), key=(i,j)=>i*10000+j;
for(let n=0;n<GALS.length;n++){const g=GALS[n];const j=Math.floor((g.de+90)/CS);const w=CS/Math.max(0.05,Math.cos(((j+0.5)*CS-90)*D2R));const k=key(j,Math.floor(g.ra/w));if(!idx.has(k))idx.set(k,[]);idx.get(k).push(n);}
// pairs: every galaxy within 6 proper Mpc (projected) and 4,000 km/s of each central of a group with M_NFW >= 5e13 Msun, 0.03 <= z <= 0.2
function pairs(mmin,zmin,zmax){const P=[];
 for(const G of GRP){if(G.m<mmin)continue;const B=bcg.get(G.id);if(!B||B.z<zmin||B.z>zmax)continue;
  const DA=Dc(B.z)/(1+B.z),th=6/DA/D2R,cb=Math.cos(B.de*D2R),sb=Math.sin(B.de*D2R);
  for(let j=Math.floor((B.de-th+90)/CS);j<=Math.floor((B.de+th+90)/CS);j++){const dec=(j+0.5)*CS-90,w=CS/Math.max(0.05,Math.cos(dec*D2R)),nc=Math.ceil(360/w),dra=th/Math.max(0.05,Math.cos(Math.min(89,Math.abs(B.de)+th)*D2R));
   for(let i=Math.floor((B.ra-dra)/w);i<=Math.floor((B.ra+dra)/w);i++){const arr=idx.get(key(j,((i%nc)+nc)%nc));if(!arr)continue;
    for(const n of arr){const g=GALS[n];if(g===B)continue;const dv=c*(g.z-B.z)/(1+B.z);if(Math.abs(dv)>4000)continue;
     const R=Math.acos(Math.min(1,sb*Math.sin(g.de*D2R)+cb*Math.cos(g.de*D2R)*Math.cos((g.ra-B.ra)*D2R)))*DA;if(R>6)continue;P.push([R,dv,G.id,n]);}}}}
 return P;}
const P=pairs(50,0.03,0.2);   // 2,779 groups, 256,242 pairs
// maximum-likelihood fit in each radial bin: two Gaussians with a common centre mu, plus a background linear in velocity, by EM
function fitEM2(V){let mu=0,s1=300,s2=800,w=0.5,b=0.5,beta=0;const N=V.length,r1=new Float64Array(N),r2=new Float64Array(N),rb=new Float64Array(N);
 for(let it=0;it<600;it++){let S1=0,S2=0,SB=0;const c1=(1-b)*w/(Math.sqrt(2*Math.PI)*s1),c2=(1-b)*(1-w)/(Math.sqrt(2*Math.PI)*s2);
  for(let i=0;i<N;i++){const d=V[i]-mu,x=V[i]/4000,p1=c1*Math.exp(-0.5*d*d/(s1*s1)),p2=c2*Math.exp(-0.5*d*d/(s2*s2)),pb=b/8000*Math.max(1e-6,1+beta*x),t=p1+p2+pb;r1[i]=p1/t;r2[i]=p2/t;rb[i]=pb/t;S1+=r1[i];S2+=r2[i];SB+=rb[i];}
  let num=0,den=0;for(let i=0;i<N;i++){const a=r1[i]/(s1*s1)+r2[i]/(s2*s2);num+=a*V[i];den+=a;}const muN=num/den;
  let q1=0,q2=0;for(let i=0;i<N;i++){const d=V[i]-muN;q1+=r1[i]*d*d;q2+=r2[i]*d*d;}
  for(let k=0;k<3;k++){let g1=0,h11=0;for(let i=0;i<N;i++){const x=V[i]/4000,f=Math.max(1e-6,1+beta*x);g1+=rb[i]*x/f;h11-=rb[i]*x*x/(f*f);}beta=Math.max(-0.9,Math.min(0.9,beta-g1/h11));}
  w=S1/(S1+S2);b=SB/N;const conv=Math.abs(muN-mu)<1e-4;mu=muN;s1=Math.max(30,Math.sqrt(q1/S1));s2=Math.max(60,Math.sqrt(q2/S2));if(conv&&it>80)break;}
 return {mu,s1,s2,w,b,beta,N};}
const EDGES=[[0,1.1],[1.1,2.1],[2.1,4.5],[4.5,6.0]], W4=[1,1,3,3];
const runBins=(P)=>EDGES.map(([a,b])=>fitEM2(P.filter(p=>p[0]>=a&&p[0]<b).map(p=>p[1])));
const integ=(R)=>R.reduce((s,r,i)=>s+W4[i]*r.mu,0)/8;
// satellites kept only if brighter than the flux limit r = 17.77 by a margin, judged at the CENTRAL's redshift + 4,000 km/s
// (offAt(z) = median m - M of all galaxies at redshift z), so a satellite's own velocity cannot decide whether it is kept
// const mcut=(zB,margin)=>17.77-(offAt(zB+0.0134)-offAt(zB))-margin;
// errors: jackknife over 20 slices of the centrals' right ascension, each holding 1/20 of the groups
