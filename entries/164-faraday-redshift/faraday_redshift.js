/* Entry 164: does the scatter of Faraday rotation grow with quasar redshift?  Run in the browser, on the VizieR CfA page
   (https://vizier.cfa.harvard.edu/viz-bin/VizieR); all data are held in browser memory and nothing is saved except the results.
   The rule was frozen in the Workbench (Entry 164, "the rule, frozen first") on 14 November 2026, before any analysis.
   Data, loaded from VizieR through /viz-bin/votable:
     window.F164.nv   NVSS rotation measures, Taylor et al. (2009), J/ApJ/702/1230/catalog, RM = -10000..10000 (37,543 sources):
                      RAJ2000, DEJ2000 (degrees), RM and e_RM (rad/m^2).
     Milliquas        VII/294 (the Million Quasars catalogue), columns RAJ2000, DEJ2000, Name, Type, z, fetched in 36 slices of 10 degrees
                      of RA, DEJ2000 > -41, z > 0, one request per Type letter (Q*, A*, B*, K*, N*); only an upper-case first letter is kept.
   Sample: Milliquas objects with Type beginning with Q, matched to the nearest NVSS RM source within 10 arcsec (one quasar per NVSS source,
     the closest), Galactic |b| >= 20 deg, at least 10 other NVSS RM sources within 3 deg.
   Foreground: dRM = RM - median RM of all other NVSS RM sources within 3 deg.
   Statistic: quasars sorted by z, 8 bins of equal count; per bin sigma = 1.4826 x MAD of dRM, error from 500 bootstrap resamples;
     weighted least-squares fit of sigma against the bin's mean log10(1+z); the slope b (rad/m^2 per unit log10(1+z)).
   Calibration: 1,000 shuffles of z among quasars within cells of |b| band (20-35-50-65-90) x tercile of e_RM (12 cells); t = (b - mean)/sd.
   Decision: t >= +3 grows; t <= -3 falls; |t| < 2 flat (95% limit b + 1.645 sd); otherwise unresolved. Checks (README): each Galactic hemisphere,
     the half with the smaller e_RM, all Milliquas classes (Q, A, B, K, N), and a calibration that also splits the sky into four 90-degree
     sectors of Galactic longitude; these use 300 shuffles each (the frozen rule fixed the count for the main calibration only).
   Use:  F164.start()   then poll F164.state (log, stage, result).  F164.state.result holds the results object.                     */
window.F164=(function(){
'use strict';
const d2r=Math.PI/180;
const CFG={matchArcsec:10,fgRadius:3,minNeigh:10,bcut:20,nbins:8,nboot:500,nshuf:1000,nshufCheck:300,seed:20261114,decMin:-41,
  bands:[20,35,50,65,90],classes:['Q','A','B','K','N']};
const state={log:[],stage:'idle',progress:'',result:null,err:null};
const t00=performance.now(), log=m=>state.log.push(((performance.now()-t00)/1000).toFixed(0)+'s '+m);
const sleep=ms=>new Promise(r=>setTimeout(r,ms));
const yieldNow=()=>new Promise(r=>{const c=new MessageChannel();c.port1.onmessage=()=>{c.port1.close();r();};c.port2.postMessage(0);});   // not throttled in a background tab
const rng=seed=>{let a=seed>>>0;return()=>{a=(a+0x6D2B79F5)>>>0;let t=Math.imul(a^(a>>>15),1|a);t=(t+Math.imul(t^(t>>>7),61|t))^t;return((t^(t>>>14))>>>0)/4294967296;};};

async function fetchText(u){for(let i=0;i<4;i++){try{const ac=new AbortController(),to=setTimeout(()=>ac.abort(),240000);const r=await fetch(u,{signal:ac.signal});clearTimeout(to);
  if(!r.ok)throw new Error('HTTP '+r.status);return await r.text();}catch(e){log('retry '+(i+1)+': '+e);await sleep(3000*(i+1));}}throw new Error('fetch failed: '+u.slice(0,140));}
function parseVO(txt){const i0=txt.indexOf('<TABLEDATA>'),i1=txt.indexOf('</TABLEDATA>');if(i0<0)return {fields:[],rows:[]};
  const fields=[...txt.slice(0,i0).matchAll(/<FIELD[^>]*\sname="([^"]+)"/g)].map(m=>m[1]);const rows=[];const trRe=/<TR>([\s\S]*?)<\/TR>/g;let m;const body=txt.slice(i0,i1);
  while((m=trRe.exec(body))){const cells=[];const tdRe=/<TD[^>]*>([^<]*)<\/TD>|<TD\s*\/>/g;let c;while((c=tdRe.exec(m[1])))cells.push(c[1]===undefined?'':c[1].trim());rows.push(cells);}return {fields,rows};}
const num=s=>s===''||s==='-'?NaN:parseFloat(s);
function sexa(s,hours){if(!/\s/.test(s))return parseFloat(s);const p=s.trim().split(/\s+/).map(Number),sg=/^-/.test(s.trim())?-1:1,v=Math.abs(p[0])+p[1]/60+p[2]/3600;return sg*v*(hours?15:1);}

// ---- Galactic coordinates and statistics helpers
const MG=[[-0.0548755604,-0.8734370902,-0.4838350155],[0.4941094279,-0.4448296300,0.7469822445],[-0.8676661490,-0.1980763734,0.4559837762]];
function gal(ra,dec){const a=ra*d2r,d=dec*d2r,v=[Math.cos(d)*Math.cos(a),Math.cos(d)*Math.sin(a),Math.sin(d)];
  const x=MG[0][0]*v[0]+MG[0][1]*v[1]+MG[0][2]*v[2],y=MG[1][0]*v[0]+MG[1][1]*v[1]+MG[1][2]*v[2],z=MG[2][0]*v[0]+MG[2][1]*v[1]+MG[2][2]*v[2];
  return [((Math.atan2(y,x)/d2r)+360)%360,Math.asin(z)/d2r];}
const medOf=a=>{const s=Float64Array.from(a).sort(),n=s.length;return n%2?s[(n-1)/2]:0.5*(s[n/2-1]+s[n/2]);};
function sigMAD(v,n,buf1,buf2){for(let i=0;i<n;i++)buf1[i]=v[i];const s=buf1.subarray(0,n).sort();const m=n%2?s[(n-1)/2]:0.5*(s[n/2-1]+s[n/2]);
  for(let i=0;i<n;i++)buf2[i]=Math.abs(s[i]-m);const d=buf2.subarray(0,n).sort();return 1.4826*(n%2?d[(n-1)/2]:0.5*(d[n/2-1]+d[n/2]));}
function wls(x,y,e){let S=0,Sx=0,Sy=0,Sxx=0,Sxy=0;for(let i=0;i<x.length;i++){const w=1/(e[i]*e[i]);S+=w;Sx+=w*x[i];Sy+=w*y[i];Sxx+=w*x[i]*x[i];Sxy+=w*x[i]*y[i];}
  const D=S*Sxx-Sx*Sx,b=(S*Sxy-Sx*Sy)/D,a=(Sy-b*Sx)/S;let chi2=0;for(let i=0;i<x.length;i++)chi2+=((y[i]-a-b*x[i])/e[i])**2;return {a,b,eb:Math.sqrt(S/D),chi2};}
// bins of equal count along the sort by z; returns {x, sig, err, n} using bootstrap errors
function binsFit(z,d,seed){const n=z.length,idx=Array.from({length:n},(_,i)=>i).sort((i,j)=>z[i]-z[j]),R=rng(seed),K=CFG.nbins,x=[],sg=[],er=[],cn=[];
  const b1=new Float64Array(n),b2=new Float64Array(n),dr=new Float64Array(n),bs=new Float64Array(CFG.nboot);
  for(let k=0;k<K;k++){const lo=Math.round(k*n/K),hi=Math.round((k+1)*n/K),m=hi-lo;let sx=0;const dv=new Float64Array(m);for(let i=0;i<m;i++){const q=idx[lo+i];sx+=Math.log10(1+z[q]);dv[i]=d[q];}
    const s0=sigMAD(dv,m,b1,b2);for(let r=0;r<CFG.nboot;r++){for(let i=0;i<m;i++)dr[i]=dv[(R()*m)|0];bs[r]=sigMAD(dr,m,b1,b2);}
    let mu=0;for(let r=0;r<CFG.nboot;r++)mu+=bs[r];mu/=CFG.nboot;let v=0;for(let r=0;r<CFG.nboot;r++)v+=(bs[r]-mu)**2;
    x.push(sx/m);sg.push(s0);er.push(Math.sqrt(v/(CFG.nboot-1)));cn.push(m);}
  const f=wls(x,sg,er);return {b:f.b,a:f.a,eb:f.eb,chi2:f.chi2,x,sigma:sg,err:er,n:cn};}
function quantile(a,q){const s=Float64Array.from(a).sort();return s[Math.min(s.length-1,Math.max(0,Math.round(q*(s.length-1))))];}

// ---- data
async function loadNVSS(){const t=await fetchText('/viz-bin/votable?-source=J/ApJ/702/1230/catalog&-out.max=unlimited&RM=-10000..10000&-out=RAJ2000&-out=DEJ2000&-out=RM&-out=e_RM&-oc.form=d');
  const {fields,rows}=parseVO(t),ix=n=>fields.indexOf(n),nv=[];
  for(const r of rows){const ra=sexa(r[ix('RAJ2000')],true),dec=sexa(r[ix('DEJ2000')],false),rm=num(r[ix('RM')]),e=num(r[ix('e_RM')]);if([ra,dec,rm,e].every(Number.isFinite))nv.push({ra,dec,rm,e});}
  log('NVSS RM rows '+rows.length+', usable '+nv.length);return nv;}
async function loadMQ(){const cands=[],counts={};const jobs=[];for(let r=0;r<360;r+=10)for(const c of CFG.classes)jobs.push([r,c]);let done=0;
  async function one([r,c]){const t=await fetchText('/viz-bin/votable?-source=VII/294&-out.max=unlimited&RAJ2000='+r+'..'+(r+10)+'&DEJ2000=>'+CFG.decMin+'&Type='+c+'*&z=>0&-out=RAJ2000&-out=DEJ2000&-out=Name&-out=Type&-out=z&-oc.form=d');
    const {fields,rows}=parseVO(t),ix=n=>fields.indexOf(n);let kept=0;
    for(const w of rows){const ty=w[ix('Type')]||'',z=num(w[ix('z')]),ra=num(w[ix('RAJ2000')]),dec=num(w[ix('DEJ2000')]);
      if(ty[0]===c&&z>0&&Number.isFinite(ra)&&Number.isFinite(dec)&&dec>CFG.decMin){cands.push({ra,dec,z,cls:c,type:ty});kept++;counts[c]=(counts[c]||0)+1;}}
    done++;state.progress='Milliquas requests '+done+'/'+jobs.length+', candidates '+cands.length;}
  let p=0;const worker=async()=>{while(p<jobs.length){const j=jobs[p++];await one(j);}};await Promise.all([worker(),worker(),worker()]);
  log('Milliquas candidates '+cands.length+' '+JSON.stringify(counts));return {cands,counts};}

// ---- matching, foreground, sample
function makeGrid(nv){const g=new Map();nv.forEach((s,i)=>{const k=Math.floor(s.ra)*1000+Math.floor(s.dec+90);(g.get(k)||g.set(k,[]).get(k)).push(i);});return g;}
function matchSet(cands,nv,grid,classSet){const best=new Map();const lim=CFG.matchArcsec;
  for(const q of cands){if(!classSet.has(q.cls))continue;let bi=-1,bd=1e9;const c0=Math.floor(q.ra),d0=Math.floor(q.dec+90);
    for(let dr=-1;dr<=1;dr++)for(let dd=-1;dd<=1;dd++){const L=grid.get((((c0+dr)%360+360)%360)*1000+(d0+dd));if(!L)continue;
      for(const i of L){const s=nv[i],dra=((s.ra-q.ra+540)%360-180)*Math.cos(q.dec*d2r)*3600,ddc=(s.dec-q.dec)*3600,sep=Math.hypot(dra,ddc);if(sep<bd){bd=sep;bi=i;}}}
    if(bi>=0&&bd<=lim){const o=best.get(bi);if(!o||bd<o.sep)best.set(bi,{sep:bd,q});}}
  return best;}
function prepNV(nv){const o=nv.map((s,i)=>({i,dec:s.dec})).sort((a,b)=>a.dec-b.dec),n=o.length,dec=new Float64Array(n),X=new Float64Array(n),Y=new Float64Array(n),Z=new Float64Array(n),rm=new Float64Array(n),id=new Int32Array(n);
  o.forEach((e,k)=>{const s=nv[e.i],a=s.ra*d2r,d=s.dec*d2r;dec[k]=s.dec;X[k]=Math.cos(d)*Math.cos(a);Y[k]=Math.cos(d)*Math.sin(a);Z[k]=Math.sin(d);rm[k]=s.rm;id[k]=e.i;});return {dec,X,Y,Z,rm,id,n};}
function foreground(P,ra,dec,self){const a=ra*d2r,d=dec*d2r,x=Math.cos(d)*Math.cos(a),y=Math.cos(d)*Math.sin(a),z=Math.sin(d),cm=Math.cos(CFG.fgRadius*d2r);
  let lo=0,hi=P.n;const t=dec-CFG.fgRadius;while(lo<hi){const m=(lo+hi)>>1;if(P.dec[m]<t)lo=m+1;else hi=m;}const v=[];
  for(let k=lo;k<P.n&&P.dec[k]<=dec+CFG.fgRadius;k++){if(P.id[k]===self)continue;if(x*P.X[k]+y*P.Y[k]+z*P.Z[k]>=cm)v.push(P.rm[k]);}
  return v.length>=CFG.minNeigh?{med:medOf(v),n:v.length}:null;}
function buildSample(best,nv,P){const S={z:[],d:[],e:[],b:[],l:[],type:[]};let nb=0,nn=0;
  for(const [i,m] of best){const s=nv[i],[l,b]=gal(s.ra,s.dec);if(Math.abs(b)<CFG.bcut){nb++;continue;}const f=foreground(P,s.ra,s.dec,i);if(!f){nn++;continue;}
    S.z.push(m.q.z);S.d.push(s.rm-f.med);S.e.push(s.e);S.b.push(b);S.l.push(l);S.type.push(m.q.type);}
  return {S,matched:best.size,droppedLat:nb,droppedNeigh:nn};}

// ---- the statistic, its calibration and the decision
function cellsFor(S,opt){const n=S.z.length,e=S.e.slice().sort((a,b)=>a-b),t1=e[Math.floor(n/3)],t2=e[Math.floor(2*n/3)];
  return S.z.map((_,i)=>{const ab=Math.abs(S.b[i]);let band=0;while(band<3&&ab>=CFG.bands[band+1])band++;const ter=S.e[i]<t1?0:S.e[i]<t2?1:2;let c=band*3+ter;if(opt&&opt.sectors)c=c*4+Math.min(3,Math.floor(S.l[i]/90));return c;});}
async function analyse(S,opt,nshuf,label){const n=S.z.length;const real=binsFit(S.z,S.d,CFG.seed);const cell=cellsFor(S,opt),groups=new Map();cell.forEach((c,i)=>(groups.get(c)||groups.set(c,[]).get(c)).push(i));
  const R=rng(CFG.seed+7),nulls=[];
  for(let s=0;s<nshuf;s++){const z2=S.z.slice();for(const [,ix] of groups){const v=ix.map(i=>S.z[i]);for(let i=v.length-1;i>0;i--){const j=(R()*(i+1))|0;[v[i],v[j]]=[v[j],v[i]];}ix.forEach((q,k)=>z2[q]=v[k]);}
    nulls.push(binsFit(z2,S.d,CFG.seed+1000+s).b);if(s%5===4){state.progress=label+': shuffle '+(s+1)+'/'+nshuf;await yieldNow();}}
  const mu=nulls.reduce((a,b)=>a+b,0)/nshuf,sd=Math.sqrt(nulls.reduce((a,b)=>a+(b-mu)**2,0)/(nshuf-1)),t=(real.b-mu)/sd;
  return {label,n,cells:groups.size,b:real.b,eb_fit:real.eb,chi2:real.chi2,bins:{x:real.x,sigma:real.sigma,err:real.err,n:real.n},
    null:{mean:mu,sd,n:nshuf,q:[0.0013,0.025,0.5,0.975,0.9987].map(q=>quantile(nulls,q)),values:nulls.map(v=>+v.toFixed(2))},
    t,p_upper:(1+nulls.filter(v=>v>=real.b).length)/(1+nshuf),p_lower:(1+nulls.filter(v=>v<=real.b).length)/(1+nshuf),limit95:real.b+1.645*sd,z_median:medOf(S.z),z_range:[Math.min(...S.z),Math.max(...S.z)]};}
const sub=(S,keep)=>{const o={};for(const k in S)o[k]=S[k].filter((_,i)=>keep(i));return o;};
function decide(t){return t>=3?'grows with redshift (t >= +3)':t<=-3?'falls with redshift (t <= -3)':Math.abs(t)<2?'flat, the standard picture (|t| < 2)':'unresolved (2 <= |t| < 3)';}

async function run(){try{
  state.stage='loading NVSS';const nv=await loadNVSS();state.stage='loading Milliquas';const {cands,counts}=await loadMQ();
  state.stage='matching';const grid=makeGrid(nv),P=prepNV(nv);
  const bestQ=matchSet(cands,nv,grid,new Set(['Q'])),bestAll=matchSet(cands,nv,grid,new Set(CFG.classes));
  const sq=buildSample(bestQ,nv,P),sa=buildSample(bestAll,nv,P);
  log('Q matched '+sq.matched+', kept '+sq.S.z.length+' (|b|<20: '+sq.droppedLat+', <10 neighbours: '+sq.droppedNeigh+'); all classes matched '+sa.matched+', kept '+sa.S.z.length);
  const res={generated:new Date().toISOString(),config:CFG,counts:{nvss:nv.length,milliquas_candidates:counts,matched_Q:sq.matched,kept_Q:sq.S.z.length,dropped_latitude_Q:sq.droppedLat,dropped_neighbours_Q:sq.droppedNeigh,matched_all:sa.matched,kept_all:sa.S.z.length}};
  state.stage='primary';res.primary=await analyse(sq.S,null,CFG.nshuf,'primary');res.decision={t:res.primary.t,outcome:decide(res.primary.t)};
  state.stage='checks';res.checks={};const S=sq.S;
  res.checks.north=await analyse(sub(S,i=>S.b[i]>0),null,CFG.nshufCheck,'north');
  res.checks.south=await analyse(sub(S,i=>S.b[i]<0),null,CFG.nshufCheck,'south');
  const e50=medOf(S.e);res.checks.low_eRM=await analyse(sub(S,i=>S.e[i]<=e50),null,CFG.nshufCheck,'low e_RM');
  res.checks.all_classes=await analyse(sa.S,null,CFG.nshufCheck,'all classes');
  res.checks.longitude_sectors=await analyse(S,{sectors:true},CFG.nshufCheck,'longitude sectors');
  const need=Object.values(res.checks).map(c=>c.t),grow=res.primary.t>=3;
  res.decision.checks_t=Object.fromEntries(Object.entries(res.checks).map(([k,c])=>[k,+c.t.toFixed(2)]));
  res.decision.growth_claim=grow?(need.every(t=>t>=1.5)?'claimed: every check reaches t >= 1.5':'unresolved, systematic suspected: not every check reaches t >= 1.5'):'not applicable (primary t < +3)';
  state.result=res;state.stage='done';log('done');
 }catch(e){state.err=String(e&&e.stack||e);state.stage='error';log('ERROR '+state.err);}}
// compact(): the saved results JSON. The raw shuffle slopes are replaced by a histogram (bin width 1 rad/m^2 per unit log10(1+z)), numbers to 5 significant digits.
function compact(res){const o=JSON.parse(JSON.stringify(res,(k,v)=>typeof v==='number'&&Number.isFinite(v)&&!Number.isInteger(v)?+v.toPrecision(5):v));
  const h=c=>{const v=c.null.values,lo=Math.floor(Math.min(...v)),hi=Math.ceil(Math.max(...v)),n=hi-lo,cnt=new Array(n).fill(0);v.forEach(x=>cnt[Math.min(n-1,Math.floor(x-lo))]++);c.null.hist={lo,width:1,counts:cnt};delete c.null.values;};
  h(o.primary);Object.values(o.checks).forEach(h);return o;}
return {state,CFG,compact,start:()=>{if(state.stage==='idle')run();return 'started';}};
})();
