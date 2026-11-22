/* Workbench Entry 22: is intergalactic space a radio fog?
   Lerner's reading of the CMB needs the space between galaxies to absorb or scatter radio waves. Its stated local signature:
   at a fixed infrared luminosity, a galaxy's radio luminosity falls with distance as D^-0.32 (1.4 GHz) to D^-0.42 (150 MHz).
   Tests: (1) IRAS PSCz galaxies (60 micron selected, with redshifts) matched to NVSS (1.4 GHz), WENSS (325 MHz) and LoTSS DR3 (144 MHz);
          (2) optically selected quasars (2QZ, LBQS) matched to NVSS: does the radio-detected share fall with redshift?
   Runs in browser tabs, as a sequence of snippets:
     PART A on heasarc.gsfc.nasa.gov (same-origin TAP queries; results come back as binary VOTables, at most 100,000 rows per query);
     PART B on vo.astron.nl (LoTSS DR3). The PSCz positions are passed from A to B through an upload job on the GAVO TAP service.
   Derived numbers and comparisons with the model are in fog_model.py; the numbers this script returned are in fog_results.json. */

/* ---------------- PART A, step 1: helpers and downloads (heasarc.gsfc.nasa.gov) ---------------- */
window.J={};window.D={};window.fire=(k,f)=>{window.J[k]='running';f().then(r=>{window.J[k]=r}).catch(e=>{window.J[k]='ERR '+e})};
// one ADQL query; decodes the BINARY serialisation of the VOTable that HEASARC returns
window.tapx=async(q,max)=>{const b=new URLSearchParams({REQUEST:'doQuery',LANG:'ADQL',MAXREC:String(max||20000),QUERY:q}); const r=await fetch('/xamin/vo/tap/sync',{method:'POST',body:b}); const t=await r.text();
  const doc=new DOMParser().parseFromString(t,'text/xml'); const info=[...doc.getElementsByTagName('INFO')].map(e=>e.getAttribute('name')+':'+e.getAttribute('value')).join('|');
  const F=[...doc.getElementsByTagName('FIELD')].map(e=>({n:e.getAttribute('name'),t:e.getAttribute('datatype'),a:e.getAttribute('arraysize')})); const st=doc.getElementsByTagName('STREAM')[0]; const rows=[];
  if(st){const bin=atob(st.textContent.replace(/\s+/g,'')); const u=new Uint8Array(bin.length); for(let i=0;i<bin.length;i++)u[i]=bin.charCodeAt(i); const dv=new DataView(u.buffer); let p=0; const dec=new TextDecoder();
    while(p<u.length){const row=[]; for(const f of F){ if(f.t==='char'){let n; if(!f.a||f.a.includes('*')){n=dv.getInt32(p);p+=4;}else n=parseInt(f.a)||1; row.push(dec.decode(u.subarray(p,p+n)).replace(/\0+$/,''));p+=n;}
      else if(f.t==='double'){row.push(dv.getFloat64(p));p+=8;} else if(f.t==='float'){row.push(dv.getFloat32(p));p+=4;} else if(f.t==='int'){row.push(dv.getInt32(p));p+=4;} else if(f.t==='short'){row.push(dv.getInt16(p));p+=2;}
      else if(f.t==='long'){row.push(Number(dv.getBigInt64(p)));p+=8;} else if(f.t==='boolean'||f.t==='unsignedByte'){row.push(dv.getUint8(p));p+=1;} else throw new Error('type '+f.t);} rows.push(row);} }
  return {info,fields:F.map(f=>f.n),rows};};
window.fire('pscz',async()=>{window.D.pscz=await window.tapx("select name, ra, dec, flux_60um, ratio_100_to_60_um_flux, radial_velocity, bii, source_type, redshift_status, optical_class_code from iraspscz",30000); return window.D.pscz.rows.length;});   // 18,411 rows
window.D.nvss=[]; window.fire('nvss',async()=>{const strips=[]; for(let d=-41;d<90;d+=3)strips.push([d,Math.min(d+3,90)]); let idx=0;
  const worker=async()=>{while(idx<strips.length){const [a,b]=strips[idx++]; const o=await window.tapx(`select ra, dec, flux_20_cm from nvss where dec>=${a} and dec<${b}`,400000); for(const r of o.rows)window.D.nvss.push(r);}};
  await Promise.all([worker(),worker(),worker()]); return window.D.nvss.length;});                                                           // 1,773,484 rows
window.fire('q2',async()=>{window.D.q2=await window.tapx("select name, ra, dec, bj_mag, redshift_1, quality_flag_1, broad_type_1 from twodfqsoz",60000); window.D.lbqs=await window.tapx("select name, ra, dec, bjmag, redshift from lbqs",5000); return 1;});
window.D.wenss=[]; window.fire('wenss',async()=>{for(const [a,b] of [[20,34],[34,40],[40,46],[46,53],[53,62],[62,91]]){const o=await window.tapx(`select ra, dec, int_flux_92_cm, flux_92_cm, source_type from wenss where dec>=${a} and dec<${b}`,300000); for(const r of o.rows)window.D.wenss.push(r);} return window.D.wenss.length;});   // 229,420 rows

/* ---------------- PART A, step 2: the galaxy test at 1.4 GHz ----------------
   window.D.pscz.rows = [name, ra, dec, f60 (Jy), f100/f60, cz (km/s), b, source_type, redshift_status, optical_class]
   window.D.nvss      = [ra, dec, flux (mJy)] */
window.FOG=(function(){
  const R=Math.PI/180, H0=70;
  let ix=null;
  function index(){ const a=window.D.nvss.slice().sort((p,q)=>p[1]-q[1]); const n=a.length;
    ix={ra:new Float64Array(n),dec:new Float64Array(n),f:new Float32Array(n)}; for(let i=0;i<n;i++){ix.ra[i]=a[i][0];ix.dec[i]=a[i][1];ix.f[i]=a[i][2];} return n; }
  function lower(d){ let lo=0,hi=ix.dec.length; while(lo<hi){const m=(lo+hi)>>1; if(ix.dec[m]<d)lo=m+1; else hi=m;} return lo; }
  // all NVSS sources within rad arcsec of (ra,dec): returns [sum of flux, nearest separation, flux of nearest, count]
  function cone(ra,dec,rad){ const rd=rad/3600, cd=Math.cos(dec*R); let s=0,best=1e9,bf=0,n=0;
    for(let i=lower(dec-rd); i<ix.dec.length && ix.dec[i]<=dec+rd; i++){ let dr=Math.abs(ix.ra[i]-ra); if(dr>180)dr=360-dr; dr*=cd; if(dr>rd)continue;
      const sep=Math.hypot(dr,ix.dec[i]-dec)*3600; if(sep<=rad){s+=ix.f[i];n++; if(sep<best){best=sep;bf=ix.f[i];}} } return [s,best,bf,n]; }
  // heliocentric -> CMB frame (dipole 369.8 km/s toward RA 167.94, Dec -6.94)
  function vcmb(ra,dec,v){ const a=167.94*R,d=-6.94*R; return v+369.8*(Math.sin(dec*R)*Math.sin(d)+Math.cos(dec*R)*Math.cos(d)*Math.cos(ra*R-a)); }
  /* build the sample. o: {fmin: 60 micron flux limit (Jy), vmin: lowest cz, rkpc: radius (kpc) inside which NVSS components are summed,
     rmin/rmax: limits on that radius (arcsec), shift: declination offset (deg) for the chance-coincidence control} */
  function sample(o){ o=Object.assign({fmin:2,vmin:2000,rkpc:10,rmin:45,rmax:150,shift:0,bmin:10},o||{}); const out=[];
    for(const r of window.D.pscz.rows){ const [name,ra,dec,f60,rat,v,b,st,rs,oc]=r; if(!(v>0)||f60<o.fmin||dec<-39.5||Math.abs(b)<o.bmin)continue; if(oc&&oc[0]!=='g')continue;
      const vc=vcmb(ra,dec,v); if(vc<o.vmin)continue; const D=vc/H0; const rad=Math.min(o.rmax,Math.max(o.rmin,o.rkpc/(D*1000)*206265));
      const c=cone(ra,dec+o.shift,rad); const det=c[3]>0; const S=det?c[0]:2.5;          // mJy; non-detections enter as the NVSS limit
      const f100=rat>0?rat*f60:NaN; const fir=1.26e-14*(2.58*f60+(f100||f60*2));           // W/m^2 (if 100 micron is missing, the typical ratio 2 is used)
      const q=Math.log10(fir/3.75e12)-Math.log10(S*1e-29), q60=Math.log10(f60*1000/S);
      out.push({name,D,z:vc/299792.458,f60,S,det,q,q60,L60:6.014+2*Math.log10(D)+Math.log10(f60),LR:20.08+2*Math.log10(D)+Math.log10(S/1000),near:c[2],sep:c[1],n:c[3],has100:rat>0}); }
    return out; }
  const med=a=>{const s=a.slice().sort((x,y)=>x-y);const n=s.length;return n?(n%2?s[(n-1)/2]:(s[n/2-1]+s[n/2])/2):NaN;};
  const mean=a=>a.reduce((x,y)=>x+y,0)/a.length;
  const sd=a=>{const m=mean(a);return Math.sqrt(a.reduce((x,y)=>x+(y-m)*(y-m),0)/(a.length-1));};
  // ordinary least squares y = X b; returns coefficients and their standard errors
  function ols(X,y){ const n=y.length,k=X[0].length; const A=Array.from({length:k},()=>new Float64Array(k)), g=new Float64Array(k);
    for(let i=0;i<n;i++)for(let a=0;a<k;a++){g[a]+=X[i][a]*y[i];for(let b=0;b<k;b++)A[a][b]+=X[i][a]*X[i][b];}
    const M=A.map((r,i)=>[...r,...Array.from({length:k},(_,j)=>i===j?1:0)]);
    for(let c=0;c<k;c++){let p=c;for(let r=c+1;r<k;r++)if(Math.abs(M[r][c])>Math.abs(M[p][c]))p=r;[M[c],M[p]]=[M[p],M[c]];const d=M[c][c];for(let j=0;j<2*k;j++)M[c][j]/=d;
      for(let r=0;r<k;r++)if(r!==c){const f=M[r][c];for(let j=0;j<2*k;j++)M[r][j]-=f*M[c][j];}}
    const inv=M.map(r=>r.slice(k)); const b=inv.map(r=>r.reduce((s,v,j)=>s+v*g[j],0)); let rss=0; for(let i=0;i<n;i++){let p=0;for(let a=0;a<k;a++)p+=X[i][a]*b[a];rss+=(y[i]-p)**2;}
    const s2=rss/(n-k); return {b,se:inv.map((r,i)=>Math.sqrt(s2*r[i])),rms:Math.sqrt(s2),n}; }
  // (a) Lerner's form: log L_radio = a + b log L_60 + c log D
  function lernerFit(S,key){ return ols(S.map(g=>[1,g.L60-10.5,Math.log10(g.D)-2]),S.map(g=>g[key||'LR'])); }
  // (b) no assumption about how radio scales with infrared: one free level per narrow luminosity bin, one common slope with distance
  function fixedL(S,w,key){ w=w||0.1; const bins={}; for(const g of S){const k=Math.floor(g.L60/w);(bins[k]=bins[k]||[]).push(g);} let sxx=0,sxy=0,n=0,nb=0; const res=[];
    for(const k in bins){const B=bins[k]; if(B.length<8)continue; nb++; const mx=mean(B.map(g=>Math.log10(g.D))), my=mean(B.map(g=>g[key||'q']));
      for(const g of B){const dx=Math.log10(g.D)-mx,dy=g[key||'q']-my; sxx+=dx*dx;sxy+=dx*dy;n++;res.push([dx,dy]);}}
    const c=sxy/sxx; let rss=0; for(const [dx,dy] of res)rss+=(dy-c*dx)**2; return {c,se:Math.sqrt(rss/(n-nb-1)/sxx),n,nb,spread:Math.sqrt(sxx/n)}; }
  // (c) medians of q in distance bins, for a slice in infrared luminosity
  function byDist(S,edges,l0,l1,key){ return edges.slice(0,-1).map((e,i)=>{const B=S.filter(g=>g.D>=e&&g.D<edges[i+1]&&g.L60>=l0&&g.L60<l1); const v=B.map(g=>g[key||'q']);
      return [e,edges[i+1],B.length,+med(v).toFixed(3),B.length>2?+(1.2533*sd(v)/Math.sqrt(B.length)).toFixed(3):null,B.filter(g=>!g.det).length,+med(B.map(g=>g.D)).toFixed(1)];}); }
  const clean=(S,lo,hi)=>S.filter(g=>g.q>(lo??1.64)&&g.q<(hi??3.0));     // drop radio-loud nuclei (q < 1.64) and infrared-excess objects, as Yun et al. (2001) do
  return {index,cone,sample,lernerFit,fixedL,byDist,clean,med,mean,sd,ols};
})();

/* run:  FOG.index(); S=FOG.sample({});  C=FOG.clean(S.filter(g=>g.det));  FOG.fixedL(C,0.1);  FOG.lernerFit(C);  FOG.byDist(S,[28,50,80,125,200,320,700],0,99)
   variants: sample({rkpc:5}), sample({rkpc:20,rmax:240}), sample({vmin:1000,rmax:300}), sample({vmin:4000}), sample({fmin:1.2}), sample({bmin:20}); control: sample({shift:0.3})
   injection check: multiply each radio flux by (D/28.6 Mpc)^-k, drop what falls below 2.5 mJy, and refit. */

/* ---------------- PART A, step 3: quasars ---------------- */
window.FOG.quasars=function(rows,RAD){ RAD=RAD||15; const F=window.FOG; const Q=rows.map(r=>({ra:r[1],dec:r[2],b:r[3],z:r[4]})).filter(g=>g.z>0&&g.dec>-39.5);
  for(const g of Q){const c=F.cone(g.ra,g.dec,RAD); g.S=c[3]?c[2]:0;                       // flux of the nearest NVSS source within RAD arcsec
    let k=0; for(const [dr,dd] of [[0,0.1],[0,-0.1],[0.12,0],[-0.12,0],[0.1,0.1],[-0.1,-0.1]]) if(F.cone(g.ra+dr,g.dec+dd,RAD)[3])k++; g.ch=k/6; }   // chance rate from six offset positions
  return Q; };
window.FOG.byZ=function(Q,zb){ return zb.slice(0,-1).map((a,i)=>{const B=Q.filter(g=>g.z>=a&&g.z<zb[i+1]); const d=B.filter(g=>g.S>0);
  return [window.FOG.med(B.map(g=>g.z)),B.length,d.length,B.reduce((s,g)=>s+g.ch,0),d.filter(g=>g.S>=10).length];}); };
/* run:  Q=FOG.quasars(D.q2.rows.filter(r=>/^QSO/.test(r[6])&&(r[5]===11||r[5]===12)));  FOG.byZ(Q,[0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.5,3.0]);
         L=FOG.quasars(D.lbqs.rows);  FOG.byZ(L,[0.2,0.7,1.2,1.7,2.2,3.4]) */

/* ---------------- PART A, step 4: the same galaxies at 325 MHz (WENSS; declination above 30 degrees) ---------------- */
window.FOG.wenss=function(o){ o=Object.assign({fmin:2,vmin:2000,decmin:30,shift:0},o||{}); const F=window.FOG,R=Math.PI/180;
  const W=window.D.wenss.filter(r=>r[4].trim()!=='C').sort((p,q)=>p[1]-q[1]); const wd=W.map(r=>r[1]);       // 'C' rows are components of multi-component sources
  const lower=d=>{let lo=0,hi=wd.length;while(lo<hi){const m=(lo+hi)>>1;if(wd[m]<d)lo=m+1;else hi=m;}return lo;};
  const cone=(ra,dec,rad)=>{const rd=rad/3600,cd=Math.cos(dec*R);let s=0,n=0;for(let i=lower(dec-rd);i<W.length&&wd[i]<=dec+rd;i++){let dr=Math.abs(W[i][0]-ra);if(dr>180)dr=360-dr;dr*=cd;if(dr>rd)continue;if(Math.hypot(dr,wd[i]-dec)*3600<=rad){s+=W[i][2];n++;}}return [s,n];};
  const P=new Map(window.D.pscz.rows.map(r=>[r[0],r])); const out=[];
  for(const g of F.sample({fmin:o.fmin,vmin:o.vmin})){const p=P.get(g.name); if(p[2]<o.decmin)continue; const rad=Math.min(150,Math.max(60,10/(g.D*1000)*206265)); const c=cone(p[1],p[2]+o.shift,rad); const det=c[1]>0; const S=det?c[0]:18;
    out.push(Object.assign({},g,{q14:g.q,det14:g.det,S325:S,det,q:g.q+Math.log10(g.S/S)}));} return out; };

/* ---------------- handing the galaxy positions to PART B ----------------
   An upload job on the GAVO TAP service (which allows cross-origin requests) holds the list; the second tab reads the job's result. */
window.relay=async function(){ const P=window.D.pscz.rows.filter(r=>r[5]>0&&r[3]>=0.6&&r[2]>=-5&&Math.abs(r[6])>=10&&r[9]&&r[9][0]==='g');     // 7,603 galaxies
  let x='<?xml version="1.0"?><VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3"><RESOURCE><TABLE><FIELD name="id" datatype="int"/><FIELD name="ra" datatype="double" unit="deg"/><FIELD name="dec" datatype="double" unit="deg"/><FIELD name="f60" datatype="double"/><FIELD name="rat" datatype="double"/><FIELD name="v" datatype="double"/><DATA><TABLEDATA>';
  P.forEach((r,i)=>{x+='<TR><TD>'+i+'</TD><TD>'+r[1]+'</TD><TD>'+r[2]+'</TD><TD>'+r[3]+'</TD><TD>'+(r[4]||0)+'</TD><TD>'+r[5]+'</TD></TR>';}); x+='</TABLEDATA></DATA></TABLE></RESOURCE></VOTABLE>';
  const fd=new FormData(); for(const [k,v] of [['REQUEST','doQuery'],['LANG','ADQL'],['FORMAT','csv'],['MAXREC','50000'],['UPLOAD','t1,param:tab'],['QUERY','select id, ra, dec, f60, rat, v from tap_upload.t1']])fd.append(k,v);
  fd.append('tab',new Blob([x],{type:'text/xml'}),'pscz.xml'); const r=await fetch('https://dc.g-vo.org/tap/async',{method:'POST',body:fd}); const job=r.url.split('/').pop();
  await fetch('https://dc.g-vo.org/tap/async/'+job+'/phase',{method:'POST',body:new URLSearchParams({PHASE:'RUN'})}); return job; };

/* ---------------- PART B: the galaxy test at 144 MHz (vo.astron.nl; JOB = the job id returned by relay()) ---------------- */
window.lotss=async function(JOB){ const t=await (await fetch('https://dc.g-vo.org/tap/async/'+JOB+'/results/result')).text(); const P=t.trim().split('\n').slice(1).map(l=>l.split(',').map(Number)); window.D.P=P; window.D.M=[];
  const mk=rows=>{let x='<?xml version="1.0"?><VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3"><RESOURCE><TABLE><FIELD name="id" datatype="int"/><FIELD name="ra" datatype="double" unit="deg"/><FIELD name="dec" datatype="double" unit="deg"/><DATA><TABLEDATA>'; for(const r of rows)x+='<TR><TD>'+r[0]+'</TD><TD>'+r[1]+'</TD><TD>'+r[2]+'</TD></TR>'; return x+'</TABLEDATA></DATA></TABLE></RESOURCE></VOTABLE>';};
  for(let i=0;i<P.length;i+=1500){ const fd=new FormData(); for(const [k,v] of [['REQUEST','doQuery'],['LANG','ADQL'],['FORMAT','csv'],['MAXREC','200000'],['UPLOAD','t1,param:tab'],
      ['QUERY',"SELECT t.id, l.ra, l.dec, l.total_flux, l.peak_flux, l.majax FROM lotss_dr3.main_sources AS l JOIN TAP_UPLOAD.t1 AS t ON 1=CONTAINS(POINT('ICRS', l.ra, l.dec), CIRCLE('ICRS', t.ra, t.dec, 0.0417))"]])fd.append(k,v);
    fd.append('tab',new Blob([mk(P.slice(i,i+1500))],{type:'text/xml'}),'p.xml'); const c=await (await fetch('/tap/sync',{method:'POST',body:fd})).text(); for(const l of c.trim().split('\n').slice(1))window.D.M.push(l.split(',').map(Number)); }
  return window.D.M.length; };                                                                                             // 35,077 LoTSS sources within 150 arcsec of 7,603 positions
window.FOG2=(function(){ const R=Math.PI/180,H0=70;
  const med=a=>{const s=a.slice().sort((x,y)=>x-y);const n=s.length;return n?(n%2?s[(n-1)/2]:(s[n/2-1]+s[n/2])/2):NaN;}; const mean=a=>a.reduce((x,y)=>x+y,0)/a.length; const sd=a=>{const m=mean(a);return Math.sqrt(a.reduce((x,y)=>x+(y-m)*(y-m),0)/(a.length-1));};
  function vcmb(ra,dec,v){ const a=167.94*R,d=-6.94*R; return v+369.8*(Math.sin(dec*R)*Math.sin(d)+Math.cos(dec*R)*Math.cos(d)*Math.cos(ra*R-a)); }
  // slope of q against x with one free level per narrow bin of the other variable: "q against distance at fixed luminosity" and its mirror
  function within(S,xf,binf,w){ const bins={}; for(const g of S){const k=Math.floor(binf(g)/w);(bins[k]=bins[k]||[]).push(g);} let sxx=0,sxy=0,n=0,nb=0; const res=[];
    for(const k in bins){const B=bins[k]; if(B.length<8)continue; nb++; const mx=mean(B.map(xf)), my=mean(B.map(g=>g.q)); for(const g of B){const dx=xf(g)-mx,dy=g.q-my; sxx+=dx*dx;sxy+=dx*dy;n++;res.push([dx,dy]);}}
    const c=sxy/sxx; let rss=0; for(const [dx,dy] of res)rss+=(dy-c*dx)**2; return {c,se:Math.sqrt(rss/(n-nb-1)/sxx),n,nb,spread:Math.sqrt(sxx/n),res}; }
  const fixedL=(S,w)=>within(S,g=>Math.log10(g.D),g=>g.L60,w||0.1), fixedD=(S,w)=>within(S,g=>g.L60,g=>Math.log10(g.D),w||0.1);
  function tab(S,f,edges){ return edges.slice(0,-1).map((e,i)=>{const B=S.filter(g=>f(g)>=e&&f(g)<edges[i+1]); const v=B.map(g=>g.q); return [B.length,med(v),B.length>2?1.2533*sd(v)/Math.sqrt(B.length):null,med(B.map(g=>g.D)),med(B.map(g=>g.L60))];}); }
  // radio flux = all LoTSS sources within 10 kpc of the IRAS position (never less than rmin or more than rmax arcsec)
  function sample(o){ o=Object.assign({fmin:0.6,vmin:2000,rkpc:10,rmin:30,rmax:90},o||{}); const by=new Map(); for(const m of window.D.M){(by.get(m[0])||by.set(m[0],[]).get(m[0])).push(m);} const out=[];
    for(const [id,ra,dec,f60,rat,v] of window.D.P){ if(f60<o.fmin)continue; const ms=by.get(id); if(!ms)continue;               // no LoTSS source within 150 arcsec: outside the survey
      const vc=vcmb(ra,dec,v); if(vc<o.vmin)continue; const D=vc/H0; const rad=Math.min(o.rmax,Math.max(o.rmin,o.rkpc/(D*1000)*206265)); const cd=Math.cos(dec*R);
      let S=0,n=0; for(const m of ms){const sep=Math.hypot((m[1]-ra)*cd,m[2]-dec)*3600; if(sep<=rad){S+=m[3];n++;}} const det=n>0; if(!det)S=1;
      const fir=1.26e-14*(2.58*f60+(rat>0?rat*f60:2*f60)); out.push({id,D,z:vc/299792.458,f60,S,det,q:Math.log10(fir/3.75e12)-Math.log10(S*1e-29),L60:6.014+2*Math.log10(D)+Math.log10(f60),LR:20.08+2*Math.log10(D)+Math.log10(S/1000)}); } return out; }
  const clean=S=>{const det=S.filter(g=>g.det); const m=med(det.map(g=>g.q)); return det.filter(g=>Math.abs(g.q-m)<0.7);};      // drop radio-loud nuclei and mismatches, symmetrically
  return {med,mean,sd,within,fixedL,fixedD,tab,sample,clean}; })();
/* run:  await lotss(JOB);  C=FOG2.clean(FOG2.sample({}));  FOG2.fixedL(C);  FOG2.fixedL(C.filter(g=>g.L60>=11));  FOG2.fixedL(C.filter(g=>g.D>=200));
         FOG2.tab(C.filter(g=>g.D>=100&&g.D<200),g=>g.L60,[10.0,10.3,10.6,10.9,11.2,11.5,12.6]);  FOG2.tab(C.filter(g=>g.L60>=11&&g.L60<11.5),g=>g.D,[60,125,200,320,500,900]) */

/* ---------------- PART C: spectra between 144 MHz and 1.4 GHz for the same galaxies ----------------
   On heasarc: the NVSS flux of each PSCz galaxy (same 10 kpc rule) is added to the position list as column s14, and the join with LoTSS is posted
   straight to the ASTRON service as an upload job (the reply cannot be read from another site, so the job is picked up from the public job list
   in the vo.astron.nl tab, started there and read there). window.D.P = [id, ra, dec, f60, f100/f60, cz, s14]; window.D.M = [id, ra, dec, total_flux]. */
window.SPX=(function(){ const R=Math.PI/180,H0=70, LF=Math.log10(1400/144);
  const med=a=>{const s=a.slice().sort((x,y)=>x-y);const n=s.length;return n?(n%2?s[(n-1)/2]:(s[n/2-1]+s[n/2])/2):NaN;}; const mean=a=>a.reduce((x,y)=>x+y,0)/a.length; const sd=a=>{const m=mean(a);return Math.sqrt(a.reduce((x,y)=>x+(y-m)*(y-m),0)/(a.length-1));};
  function vcmb(ra,dec,v){ const a=167.94*R,d=-6.94*R; return v+369.8*(Math.sin(dec*R)*Math.sin(d)+Math.cos(dec*R)*Math.cos(d)*Math.cos(ra*R-a)); }
  function within(S,xf,binf,w,yk){ const bins={}; for(const g of S){const k=Math.floor(binf(g)/w);(bins[k]=bins[k]||[]).push(g);} let sxx=0,sxy=0,n=0,nb=0; const res=[];
    for(const k in bins){const B=bins[k]; if(B.length<8)continue; nb++; const mx=mean(B.map(xf)), my=mean(B.map(g=>g[yk])); for(const g of B){const dx=xf(g)-mx,dy=g[yk]-my; sxx+=dx*dx;sxy+=dx*dy;n++;res.push([dx,dy]);}}
    const c=sxy/sxx; let rss=0; for(const [dx,dy] of res)rss+=(dy-c*dx)**2; return {c,se:Math.sqrt(rss/(n-nb-1)/sxx),n,nb}; }
  const atFixedL=(S,yk)=>within(S,g=>Math.log10(g.D),g=>g.L60,0.1,yk), atFixedD=(S,yk)=>within(S,g=>g.L60,g=>Math.log10(g.D),0.1,yk);
  function tab(S,f,edges,yk){ return edges.slice(0,-1).map((e,i)=>{const B=S.filter(g=>f(g)>=e&&f(g)<edges[i+1]); const v=B.map(g=>g[yk]); return [med(B.map(f)),B.length,med(v),B.length>2?1.2533*sd(v)/Math.sqrt(B.length):null];}); }
  // galaxies with a flux at both frequencies; alpha > 0 means the flux falls with frequency
  function sample(o){ o=Object.assign({fmin:2,vmin:2000,rkpc:10,rmin:30,rmax:90},o||{}); const by=new Map(); for(const m of window.D.M){(by.get(m[0])||by.set(m[0],[]).get(m[0])).push(m);} const out=[];
    for(const [id,ra,dec,f60,rat,v,s14] of window.D.P){ if(f60<o.fmin)continue; const ms=by.get(id); if(!ms)continue; const vc=vcmb(ra,dec,v); if(vc<o.vmin)continue; const D=vc/H0; const rad=Math.min(o.rmax,Math.max(o.rmin,o.rkpc/(D*1000)*206265)); const cd=Math.cos(dec*R);
      let S=0,n=0; for(const m of ms){const sep=Math.hypot((m[1]-ra)*cd,m[2]-dec)*3600; if(sep<=rad){S+=m[3];n++;}} if(!n||!(s14>0))continue;
      const lf=Math.log10(1.26e-14*(2.58*f60+(rat>0?rat*f60:2*f60))/3.75e12); out.push({id,D,f60,alpha:Math.log10(S/s14)/LF,q144:lf-Math.log10(S*1e-29),q14:lf-Math.log10(s14*1e-29),L60:6.014+2*Math.log10(D)+Math.log10(f60)}); } return out; }
  const clean=S=>{const m=med(S.map(g=>g.alpha)); return S.filter(g=>Math.abs(g.alpha-m)<0.6&&g.q14>1.64&&g.q14<3.0);};
  return {med,mean,sd,within,atFixedL,atFixedD,tab,sample,clean}; })();
/* run:  C=SPX.clean(SPX.sample({fmin:1.2}));  SPX.atFixedD(C.filter(g=>g.L60>=10.4),'alpha');  SPX.atFixedL(C,'alpha');  SPX.atFixedD(C.filter(g=>g.L60>=10.4),'q144');
         SPX.tab(C.filter(g=>g.D>=100&&g.D<200),g=>g.L60,[10.2,10.5,10.8,11.1,11.4,12.6],'alpha') */

/* ---------------- PART D: direction (Entry 22, Part 8) ----------------
   Lerner (2003): absorption should be stronger along the axis of the Local Supercluster (toward Virgo) and weaker at right angles.
   For each galaxy of FOG.clean(FOG.sample({fmin:1.2})): residual r = q - mean q of its 0.1-dex luminosity bin; cosV = cosine of the angle to Virgo
   (RA 187.706, Dec +12.391); sgb = supergalactic latitude (pole at l = 47.37, b = +6.32). Then mean r in bins of acos|cosV| and of |sgb|,
   and the slope of r against |cosV| and against |sin sgb|. Numbers in fog_results.json under "direction". */
