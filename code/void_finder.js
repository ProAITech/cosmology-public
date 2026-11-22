// Entry 148: projected-void finder run in the browser on the 2M++ catalogue (Lavaux & Hudson 2011), VizieR J/MNRAS/416/2840/catalog.
// Run on a VizieR mirror page (same-origin fetch). Produces the two frozen catalogues voids_slab1.txt and voids_slab2.txt.
const u='/viz-bin/asu-tsv?-source=J/MNRAS/416/2840/catalog&-out=Ksmag,Vcmb,Cln,M0,M1,M2,Ref&-out.add=_RA,_DE&-out.max=100000&-out.form=|';
const t=await (await fetch(u)).text();
const GAL=t.split('\n').filter(x=>x&&x[0]!=='#').slice(3).map(x=>{const p=x.split('|');return {ra:+p[0],de:+p[1],K:+p[2],v:+p[3],cln:+p[4],m0:+p[5],m1:+p[6],m2:+p[7],zoa:/zoa/i.test(p[8])}}).filter(g=>isFinite(g.ra)&&isFinite(g.v)&&g.v>0);
// 72,940 rows, 3,808 Zone-of-Avoidance fakes (excluded); m1 = SDSS region, m2 = 6dFGRS region (K < 12.5); m0 = 2MRS-only region (K < 11.5)
const D2R=Math.PI/180, c=299792.458, H0=67.7, Om=0.31;
function gal(ra,de){const a=ra*D2R,d=de*D2R,aG=192.85948*D2R,dG=27.12825*D2R,lNCP=122.93192*D2R;
 const b=Math.asin(Math.sin(d)*Math.sin(dG)+Math.cos(d)*Math.cos(dG)*Math.cos(a-aG));
 const y=Math.cos(d)*Math.sin(a-aG),x=Math.sin(d)*Math.cos(dG)-Math.cos(d)*Math.sin(dG)*Math.cos(a-aG);
 return [((lNCP-Math.atan2(y,x))/D2R%360+360)%360,b/D2R];}
function Dc(z){let s=0,n=200;for(let i=0;i<n;i++){const zz=(i+0.5)/n*z;s+=1/Math.sqrt(Om*(1+zz)**3+1-Om);}return c/H0*s*z/n;}
const Mlim=(z,K)=>K-5*Math.log10(Dc(z)*(1+z)*1e5)+6*Math.log10(1+z);   // K-correction -6 log10(1+z)
for(const g of GAL){const [l,b]=gal(g.ra,g.de);g.l=l;g.b=b;g.z=g.v/c;g.M=g.K-5*Math.log10(Dc(g.z)*(1+g.z)*1e5)+6*Math.log10(1+g.z);g.deep=(g.m1||g.m2)?1:0;}
// sky map of survey depth: 2 deg in l x equal steps in sin b (16,200 cells); a cell is "deep" if most of its galaxies are
const NL=180,NB=90,cnt=new Array(NL*NB).fill(0),dp=new Array(NL*NB).fill(0);
const cellOf=(l,b)=>Math.min(NB-1,Math.floor((Math.sin(b*D2R)+1)/2*NB))*NL+Math.min(NL-1,Math.floor(l/2));
for(const g of GAL){if(g.zoa)continue;const k=cellOf(g.l,g.b);cnt[k]++;dp[k]+=g.deep;}
const deepMap=cnt.map((n,k)=>n>0?(dp[k]/n>0.5?1:0):-1);   // empty cells filled from neighbours (20 passes, same l-row or b-column)
const inMask=(l,b,needDeep)=>Math.abs(b)>=10&&(!needDeep||deepMap[cellOf(l,b)]===1);
// deep region: 24,062 deg2 (58.3% of the sky, 70.9% of |b| >= 10 deg)
const vec=(l,b)=>[Math.cos(b*D2R)*Math.cos(l*D2R),Math.cos(b*D2R)*Math.sin(l*D2R),Math.sin(b*D2R)];
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;}}
const rnd=mulberry32(145146), Rl=[]; for(let i=0;i<150000;i++){const z=2*rnd()-1,ph=2*Math.PI*rnd();Rl.push([Math.asin(z)/D2R,ph/D2R]);}
const RAND={0:Rl.filter(([b,l])=>inMask(l,b,0)).map(([b,l])=>vec(l,b)),1:Rl.filter(([b,l])=>inMask(l,b,1)).map(([b,l])=>vec(l,b))};
const TR={ s1:GAL.filter(g=>!g.zoa&&g.z>0.002&&g.z<0.03&&g.K<11.5&&g.M<Mlim(0.03,11.5)&&inMask(g.l,g.b,0)).map(g=>vec(g.l,g.b)),   // 6,589 tracers, M_K < -24.09
           s2:GAL.filter(g=>!g.zoa&&g.z>=0.03&&g.z<0.06&&g.K<12.5&&g.deep&&g.M<Mlim(0.06,12.5)&&inMask(g.l,g.b,1)).map(g=>vec(g.l,g.b))};  // 11,075 tracers, M_K < -24.57
// spherical-underdensity finder on the sky: for each of ~9,000 Fibonacci centres, the largest radius (0.5 deg steps to 25 deg)
// with Delta(<theta) <= -0.5, at least 8 tracers expected, and >= 80% of the cap unmasked; then non-overlapping by radius
function findVoids(key,deepFlag,thr,step=0.5,tmax=25,Ncent=9000,minExp=8){
 const T=TR[key],R=RAND[deepFlag],nbar=T.length/R.length,nk=Math.round(tmax/step),cmin=Math.cos(tmax*D2R),NRsky=150000;
 const C=[],ga=Math.PI*(3-Math.sqrt(5));
 for(let i=0;i<Ncent;i++){const z=1-2*(i+0.5)/Ncent,r=Math.sqrt(1-z*z),ph=i*ga,v=[r*Math.cos(ph),r*Math.sin(ph),z],b=Math.asin(z)/D2R,l=((Math.atan2(v[1],v[0])/D2R)+360)%360;if(inMask(l,b,deepFlag))C.push({v,l,b});}
 const out=[];
 for(const c of C){const ht=new Float64Array(nk+1),hr=new Float64Array(nk+1);
  for(const t of T){const d=t[0]*c.v[0]+t[1]*c.v[1]+t[2]*c.v[2];if(d<cmin)continue;const k=Math.floor(Math.acos(Math.min(1,d))/D2R/step);if(k<nk)ht[k]++;}
  for(const r of R){const d=r[0]*c.v[0]+r[1]*c.v[1]+r[2]*c.v[2];if(d<cmin)continue;const k=Math.floor(Math.acos(Math.min(1,d))/D2R/step);if(k<nk)hr[k]++;}
  let st=0,sr=0,best=null;
  for(let k=0;k<nk;k++){st+=ht[k];sr+=hr[k];const th=(k+1)*step,exp=nbar*sr,comp=sr/(NRsky*(1-Math.cos(th*D2R))/2);if(comp<0.8)break;if(exp<minExp)continue;const del=st/exp-1;if(del<=thr)best={th,del,exp,n:st,comp};}
  if(best)out.push({l:c.l,b:c.b,v:c.v,...best});}
 out.sort((a,b)=>b.th-a.th||a.del-b.del);
 const acc=[];for(const o of out){let ok=true;for(const a of acc){const d=Math.acos(Math.min(1,o.v[0]*a.v[0]+o.v[1]*a.v[1]+o.v[2]*a.v[2]))/D2R;if(d<Math.max(o.th,a.th)){ok=false;break;}}if(ok)acc.push(o);}
 return acc;}
const lpcdf=(n,mu)=>{let s=Math.exp(-mu),term=-mu;for(let k=1;k<=n;k++){term+=Math.log(mu)-Math.log(k);s+=Math.exp(term);}return Math.log10(s);};
const S1=findVoids('s1',0,-0.5).filter(x=>lpcdf(x.n,x.exp)<-3);   // 36 voids
const S2=findVoids('s2',1,-0.5).filter(x=>lpcdf(x.n,x.exp)<-3);   // 62 voids
const sep=(a,b)=>Math.acos(Math.min(1,a.v[0]*b.v[0]+a.v[1]*b.v[1]+a.v[2]*b.v[2]))/D2R;
for(const v of S2) v.ov=S1.some(u=>sep(u,v)<0.5*(u.th+v.th));     // 19 overlap a slab-1 inner disk; 43 are on fresh sky
