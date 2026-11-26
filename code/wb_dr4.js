// Workbench Entry 114, amendment 1 (26 November 2026, before Gaia DR4): the DR4 rerun of the Entry 15 wide-binary
// forward model, made exact. Load wide_binaries_model.js first (unchanged; it defines __wbRun, __wbFine and __wbFns),
// then this file. Everything runs in browser memory; nothing is written to disk.
//
// 1. parent()   The Entry 15 parent sample, from El-Badry, Rix & Heintz (2021), VizieR J/MNRAS/506/2269/catalog,
//               queried from vizier.cfa.harvard.edu (it allows cross-site reads, so this also works from the Gaia
//               archive's own pages): both stars main sequence (BinType = MSMS), chance-alignment probability R < 0.01,
//               EDR3 RUWE < 1.2 for both stars, parallax of the first star > 5 mas. On 26 November 2026: 49,940 pairs.
//               Through wide_binaries_model.js's own cuts they give 39,607 pairs with noise below 25% of the circular
//               speed and 22,945 with noise below 50 m/s; Entry 15 had 39,602 and 22,942 (its query script,
//               wide_binaries.js, was not kept; this one reproduces its data medians in all ten separation bins).
// 2. dr4()      For each star: its DR4 counterpart (gaiadr4.dr3_neighbourhood, the nearest listed DR4 source), its DR4
//               position, parallax, proper motion, their errors and G (gaiadr4.gaia_source), and a companion flag:
//               an entry in any DR4 non-single-star table listed in NSS below, membership of a resolved pair
//               (nss_resolved_pair) with any source other than its wide companion, or RUWE > 1.2 in gaia_source.
//               Run in a tab on https://gea.esac.esa.int/ (same-origin TAP). checkSchema() first confirms that the
//               tables and columns named here exist and lists any that do not; a substitution is reported with the result.
// 3. layout()   Rows in the 27-column layout wide_binaries_model.js reads, with DR4 values in place of EDR3 values.
// 4. runs       B (decides): flagged pairs removed, companion fraction free as in Entry 15, drawn twice with 16 and 17
//               orbits per pair and decided on the mean of the two chi^2 profiles.
//               A (check, the commitment as first written): all pairs, companion fraction fixed at f_DR4, the share of
//               pairs flagged in that sample. The DR3 rehearsal showed A is biased: f_DR4 is a lower limit, and fixing the
//               fraction there makes the model put the undetected companions' fast tail into gravity.
//               Each on both noise samples. Gravity: Newton x [1 + mu (MOND boost - 1)], galaxy-calibrated interpolation,
//               mu = 0 to 0.3 in steps of 0.025, then 0.4, 0.5, 0.7, 1.0; boost = 16.8 mu per cent (Entry 15).
// 5. readBoost  The chi^2 profile over mu (minimised over everything else) is divided by the reduced chi^2 when that
//               exceeds 1 (30 numbers less the fitted parameters: 4 in A, 5 in B). The boost is the vertex of a parabola
//               through the lowest grid point and its two neighbours; its 68% interval is where the scaled profile rises
//               by 1. Readings, as committed in Entry 114: below 2% favours Newton; 3-4% matches MOND with the sharper
//               ("standard") interpolation; 12% or more (gamma >= 1.25) revives galaxy-calibrated MOND; in between, no
//               reading is named. A reading needs both noise samples of run B to give it, and none is named from a fit
//               whose reduced chi^2 exceeds 5 (Entry 15's fits: 3.0 and 2.0 with 4 orbits per pair; the rehearsal's
//               decision runs 3.9 and 2.6 with 16, which shrink the model-noise allowance in chi^2 by 1.18).
// Rehearsal (DR3, before DR4): mode 'dr3' uses the parent rows' own EDR3 values and stand-in flags that mark close
// resolved companions (ipd_frac_multi_peak > 2 or ipd_gof_harmonic_amplitude > 0.1 for either star). A pipeline check.
(()=>{
const VIZ='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv';
const COLS=['RA_ICRS','DE_ICRS','RA2deg','DE2deg','Plx1','Plx2','e_Plx1','e_Plx2','pmRA1','pmRA2','pmDE1','pmDE2','e_pmRA1','e_pmRA2',
  'e_pmDE1','e_pmDE2','Gmag1','Gmag2','RUWE1','RUWE2','BP-RP1','BP-RP2','sepAU','R','BinType','Source1','Source2',
  'IPDfmp1','IPDfmp2','IPDgofha1','IPDgofha2'];
const MU=[0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.4,0.5,0.7,1.0], BOOST_PER_MU=16.8;
const FC_FREE=[0,0.05,0.1,0.15,0.2,0.3,0.4];
async function parent(){
  const e=encodeURIComponent;
  const q=VIZ+'?-source='+e('J/MNRAS/506/2269/catalog')+'&-out.max=unlimited&-out='+e(COLS.join(','))+'&BinType=MSMS&R='+e('<0.01')+
    '&RUWE1='+e('<1.2')+'&RUWE2='+e('<1.2')+'&Plx1='+e('>5');
  const t=await (await fetch(q)).text();
  const rows=t.split('\n').filter(l=>l&&l[0]!=='#'&&l[0]!=='-'&&!isNaN(parseFloat(l))).map(l=>l.split('\t').map(x=>x.trim()));
  return rows;                                         // rows[i][25], rows[i][26]: the two EDR3 source_ids
}
const standIn=r=>(+r[27]>2||+r[28]>2||+r[29]>0.1||+r[30]>0.1);
// ---- Gaia archive TAP (asynchronous jobs, so no row cap applies), run from a page on gea.esac.esa.int
async function tap(adql){
  const fd=new URLSearchParams({REQUEST:'doQuery',LANG:'ADQL',FORMAT:'csv',PHASE:'RUN',QUERY:adql});
  const r=await fetch('/tap-server/tap/async',{method:'POST',body:fd,redirect:'follow'});
  const job=r.url.replace(/\/$/,'');
  for(let i=0;i<600;i++){ const ph=(await (await fetch(job+'/phase')).text()).trim();
    if(ph==='COMPLETED') break; if(ph==='ERROR'||ph==='ABORTED') throw new Error('TAP '+ph+': '+(await (await fetch(job+'/error')).text()).slice(0,300));
    await new Promise(z=>setTimeout(z,2000)); }
  const csv=await (await fetch(job+'/results/result')).text(), L=csv.trim().split('\n'), H=L[0].split(',');
  return L.slice(1).map(l=>{ const v=l.split(','), o={}; H.forEach((h,i)=>o[h]=v[i]); return o; });
}
const NB_TABLE='gaiadr4.dr3_neighbourhood', NB_DR3='dr3_source_id', NB_DR4='dr4_source_id', NB_DIST='angular_distance';
const SRC='gaiadr4.gaia_source', SRC_COLS=['source_id','ra','dec','parallax','parallax_error','pmra','pmdec','pmra_error','pmdec_error','phot_g_mean_mag','ruwe'];
const NSS=['nss_two_body_orbit','nss_acceleration_astro','nss_non_linear_spectro','nss_multiple_orbits','nss_vim_fl','rvs_two_body_orbit',
  'nss_two_body_orbit_supp','nss_acceleration_astro_supp','nss_multiple_orbits_supp'].map(t=>'gaiadr4.'+t);
const PAIR_TABLES=['gaiadr4.nss_resolved_pair','gaiadr4.nss_resolved_pair_supp'];
async function checkSchema(){
  const want=[[NB_TABLE,[NB_DR3,NB_DR4,NB_DIST]],[SRC,SRC_COLS]].concat(NSS.map(t=>[t,['source_id']])).concat(PAIR_TABLES.map(t=>[t,[]]));
  const rows=await tap("SELECT table_name, column_name FROM TAP_SCHEMA.columns WHERE table_name LIKE 'gaiadr4.%'");
  const have={}; rows.forEach(r=>{ (have[r.table_name]=have[r.table_name]||new Set()).add(r.column_name); });
  const missing=[]; want.forEach(([t,cs])=>{ if(!have[t]) missing.push(t+' (table)'); else cs.forEach(c=>{ if(!have[t].has(c)) missing.push(t+'.'+c); }); });
  const pairCols=Object.fromEntries(PAIR_TABLES.filter(t=>have[t]).map(t=>[t,[...have[t]].filter(c=>/source_id/.test(c))]));
  return {missing, pairCols};
}
const chunks=(a,n)=>{ const o=[]; for(let i=0;i<a.length;i+=n) o.push(a.slice(i,i+n)); return o; };
async function dr4(rows, log=()=>{}){
  const ids3=[...new Set(rows.flatMap(r=>[r[25],r[26]]))], map={}, astro={}, flag=new Set();
  for(const c of chunks(ids3,2000)){                    // nearest listed DR4 counterpart of each DR3 source
    const res=await tap(`SELECT ${NB_DR3}, ${NB_DR4}, ${NB_DIST} FROM ${NB_TABLE} WHERE ${NB_DR3} IN (${c.join(',')})`);
    res.forEach(o=>{ const k=o[NB_DR3], d=+o[NB_DIST]; if(!map[k]||d<map[k].d) map[k]={id:o[NB_DR4],d}; }); log('neighbourhood',Object.keys(map).length); }
  const ids4=[...new Set(Object.values(map).map(m=>m.id))];
  for(const c of chunks(ids4,2000)){
    (await tap(`SELECT ${SRC_COLS.join(', ')} FROM ${SRC} WHERE source_id IN (${c.join(',')})`)).forEach(o=>{ astro[o.source_id]=o; if(+o.ruwe>1.2) flag.add(o.source_id); });
    log('gaia_source',Object.keys(astro).length); }
  for(const t of NSS) for(const c of chunks(ids4,2000)){
    try{ (await tap(`SELECT source_id FROM ${t} WHERE source_id IN (${c.join(',')})`)).forEach(o=>flag.add(o.source_id)); }catch(e){ log('skip '+t,e.message); break; } }
  const partner={}; rows.forEach(r=>{ const a=map[r[25]]&&map[r[25]].id, b=map[r[26]]&&map[r[26]].id; if(a&&b){ partner[a]=b; partner[b]=a; } });
  const sch=await checkSchema();
  for(const [t,cs] of Object.entries(sch.pairCols)) if(cs.length>=2) for(const c of chunks(ids4,1000)){
    const w=cs.map(x=>`${x} IN (${c.join(',')})`).join(' OR ');
    (await tap(`SELECT ${cs.join(', ')} FROM ${t} WHERE ${w}`)).forEach(o=>{ const v=cs.map(x=>o[x]);
      v.forEach(s=>{ if(astro[s]&&v.some(u=>u!==s&&u!==partner[s])) flag.add(s); }); }); }
  return {map, astro, flag, schema:sch};
}
// ---- rows in the model's layout; DR4 values replace columns 0-17 and 25-26 (now DR4 ids); flag in column 31
function layout(rows, D){
  const out=[], lost={noMatch:0,noAstro:0,same:0};
  for(const r of rows){ const m1=D.map[r[25]], m2=D.map[r[26]];
    if(!m1||!m2){ lost.noMatch++; continue; } if(m1.id===m2.id){ lost.same++; continue; }
    const a=D.astro[m1.id], b=D.astro[m2.id]; if(!a||!b||a.pmra===''||b.pmra===''||a.parallax===''||b.parallax===''){ lost.noAstro++; continue; }
    const v=[a.ra,a.dec,b.ra,b.dec,a.parallax,b.parallax,a.parallax_error,b.parallax_error,a.pmra,b.pmra,a.pmdec,b.pmdec,
      a.pmra_error,b.pmra_error,a.pmdec_error,b.pmdec_error,a.phot_g_mean_mag,b.phot_g_mean_mag,a.ruwe,b.ruwe,r[20],r[21],r[22],r[23],r[24],m1.id,m2.id,r[27],r[28],r[29],r[30],
      (D.flag.has(m1.id)||D.flag.has(m2.id))?1:0];
    out.push(v.map(String)); }
  return {rows:out, lost};
}
// ---- the two runs on one noise sample, each in its own worker; fFlag is measured on the pairs that pass the model's cuts
function counted(rows, errCut, absCut){ const F=window.__wbFns, N=[]; rows.forEach(r=>{ const P=F.pairs(r.join('\t'),errCut,absCut); if(P.length&&P[0].s>=200&&P[0].s<30000) N.push(r); }); return N; }
// The grid is split over workers by mu. The model draws its random numbers from a fixed seed, so every worker sees the
// same draws and the merged result is the one a single worker would give. R = draws per pair: 16 for DR4 (Entry 15 used 4).
function launch(tag, rows, errCut, absCut, FC, R, split){
  const text=rows.map(r=>r.join('\t')).join('\n'), n=Math.ceil(MU.length/split);
  for(let i=0;i*n<MU.length;i++) window.__launch({tag:tag+'.'+i, text, errCut, absCut, mu:MU.slice(i*n,(i+1)*n), FC, R}); }
function merged(tag){
  const parts=Object.keys(window.__res).filter(k=>k.startsWith(tag+'.')).map(k=>window.__res[k].summary), best={};
  parts.forEach(p=>Object.assign(best,p.best)); return {best, N:parts.length?parts[0].N:0, data:parts.length?parts[0].data:null, parts:parts.length}; }
// The decision run is drawn twice, with 16 and 17 orbits per pair (independent random numbers), and decided on the mean
// of the two chi^2 profiles: combine(merged(prefix+'B25r16'), merged(prefix+'B25r17')). If the two draws, read alone,
// give different readings, no reading is named.
function startRuns(rows, isFlagged, prefix, {split=6}={}){
  const S={}; for(const [name,eC,aC] of [['25',0.25,null],['50',null,0.050]]){
    const kept=counted(rows,eC,aC), fl=kept.filter(isFlagged), f=fl.length/kept.length, clean=kept.filter(r=>!isFlagged(r));
    S[name]={N:kept.length, flagged:fl.length, fFlag:+f.toFixed(4)};
    launch(prefix+'B'+name+'r16', clean, eC, aC, FC_FREE, 16, split);
    launch(prefix+'B'+name+'r17', clean, eC, aC, FC_FREE, 17, split);
    launch(prefix+'A'+name, kept, eC, aC, [+f.toFixed(4)], 16, split); }
  return S;
}
function combine(a, b){ const best={}; for(const k of Object.keys(a.best)) if(b.best[k]) best[k]={chi2:(a.best[k].chi2+b.best[k].chi2)/2}; return {best}; }
// ---- the boost from a finished run's summary
function readBoost(summary, nPar){
  const prof=MU.map(mu=>{ const b=summary.best['mu'+mu]; return [mu*BOOST_PER_MU, b?b.chi2:NaN]; }).filter(p=>isFinite(p[1]));
  let k=0; prof.forEach((p,i)=>{ if(p[1]<prof[k][1]) k=i; });
  const c2=prof[k][1], dof=30-nPar, s=Math.max(1,c2/dof), D=prof.map(p=>[p[0],(p[1]-c2)/s]);
  let best=D[k][0];
  if(k>0&&k<D.length-1){ const [x0,y0]=D[k-1],[x1,y1]=D[k],[x2,y2]=D[k+1];
    const den=(x0-x1)*(x0-x2)*(x1-x2), A=(x2*(y1-y0)+x1*(y0-y2)+x0*(y2-y1))/den, B=(x2*x2*(y0-y1)+x1*x1*(y2-y0)+x0*x0*(y1-y2))/den;
    if(A>0) best=-B/(2*A); }
  const cross=(dir)=>{ for(let i=k;i>=0&&i<D.length;i+=dir){ const j=i+dir; if(j<0||j>=D.length) return null;
      if(D[j][1]>=1){ const [xa,ya]=D[i],[xb,yb]=D[j]; return xa+(1-ya)*(xb-xa)/(yb-ya); } } return null; };
  const c=cross(-1), lo=k===0||c===null?0:c, hi=cross(1);   // lo 0: consistent with no boost; hi null: beyond the grid's 16.8%
  const reading=b=>c2/dof>5?'none (poor fit)':b<2?'Newton':(b>=3&&b<=4)?'sharp MOND':(b>=12)?'galaxy-calibrated MOND':'between readings';
  return {boost:+best.toFixed(2), lo:lo===null?null:+lo.toFixed(2), hi:hi===null?null:+hi.toFixed(2), chi2min:+c2.toFixed(2), scale:+s.toFixed(3),
    newton_sigma:+Math.sqrt(Math.max(D[0][1],0)).toFixed(2), galaxy_mond_sigma:+Math.sqrt(Math.max(D[D.length-1][1],0)).toFixed(2),
    reading:reading(best), profile:D.map(p=>[+p[0].toFixed(2),+p[1].toFixed(2)])};
}
// ---- load the unchanged model (from the public repository, checked against its sha256) into the page and into workers
const MODEL_URL='https://raw.githubusercontent.com/ProAITech/cosmology-public/main/code/wide_binaries_model.js';
const MODEL_SHA256='8fc27cefe77ac28729e1b7bb98ca9bcc35ee0b907dd200b48268203a601a490e';
async function loadModel(){
  const src=await (await fetch(MODEL_URL)).text();
  const h=[...new Uint8Array(await crypto.subtle.digest('SHA-256',new TextEncoder().encode(src)))].map(b=>b.toString(16).padStart(2,'0')).join('');
  if(h!==MODEL_SHA256) throw new Error('model checksum '+h);
  (0,eval)(src);
  const wsrc='self.window=self;\n'+src+`
self.onmessage=e=>{ const m=e.data; window.__wbT=m.text; const F=window.__wbFns;
  let o=window.__wbFine(); o.manual=false; if(m.FC) o.FC=m.FC; if(m.R) o.R=m.R;
  const bR=F.beta(F.NU.rar); o.grav=m.mu.map(mu=>['mu'+mu,y=>1+mu*(bR(y)-1)]);
  const t0=Date.now(); o.then=()=>postMessage({tag:m.tag,done:true,secs:(Date.now()-t0)/1000,summary:window.__wbm[m.tag].summary});
  window.__wbRun(m.tag,m.errCut,m.absCut,o);
  const iv=setInterval(()=>{const s=window.__wbm[m.tag]; if(!s) return; postMessage({tag:m.tag,progress:[s.state,s.done,s.total]}); if(s.state==='done') clearInterval(iv);},3000); };`;
  window.__res=window.__res||{}; window.__prog=window.__prog||{};
  window.__launch=(msg)=>{ const w=new Worker(URL.createObjectURL(new Blob([wsrc],{type:'text/javascript'})));
    w.onmessage=e=>{ const d=e.data; if(d.done) window.__res[d.tag]=d; else window.__prog[d.tag]=d; };
    w.onerror=e=>{ window.__prog[msg.tag]={err:e.message}; }; w.postMessage(msg); return w; };
  return h;
}
window.__wbDR4={loadModel, parent, standIn, tap, checkSchema, dr4, layout, counted, startRuns, merged, combine, readBoost, MU, FC_FREE, COLS};
})();
// DR3 rehearsal (a tab on vizier.cfa.harvard.edu):
//   await __wbDR4.loadModel();  const rows=await __wbDR4.parent();  __wbDR4.startRuns(rows, __wbDR4.standIn, 'r');
//   then __wbDR4.readBoost(__wbDR4.combine(__wbDR4.merged('rB25r16'), __wbDR4.merged('rB25r17')), 5) (decision),
//   __wbDR4.readBoost(__wbDR4.merged('rA25'), 4) (check), and the same for 50. (The rehearsal of 26 November ran its
//   runs with 4 orbits per pair, one worker each, and the decision runs once more with 16; wb_dr4_rehearsal.json.)
// DR4 (a tab on gea.esac.esa.int, same steps):
//   await __wbDR4.loadModel();  await __wbDR4.checkSchema();  const rows=await __wbDR4.parent();  const D=await __wbDR4.dr4(rows, console.log);
//   const L=__wbDR4.layout(rows, D);  __wbDR4.startRuns(L.rows, r=>r[31]==='1', 'd');
