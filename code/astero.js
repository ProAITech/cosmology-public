// Entry 157 — browser code run on the VizieR CfA page.
// Seismic catalogues: LEGACY (J/ApJ/835/172 table1), APOKASC dwarfs (J/ApJS/233/23 tables 3 and 4), APOKASC-2 giants (J/ApJS/239/32 table5).
// Radii: Berger et al. 2020 Gaia–Kepler catalogue (J/AJ/159/280 table2), matched by KIC.
const get=async(src,cols,max)=>{const u='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source='+encodeURIComponent(src)+'&-out.max='+max+'&-out='+encodeURIComponent(cols)+'&-out.form=|';
  const t=await fetch(u).then(r=>r.text()); const L=t.split('\n').filter(l=>l&&!l.startsWith('#'));
  return L.slice(3).filter(l=>/\d/.test(l)).map(l=>l.split('|').map(x=>x.trim()));};
async function run(){
  const gi=await get('J/ApJS/239/32/table5','KIC,Teff,Numax,Dnu,ES',20000), dw3=await get('J/ApJS/233/23/table3','KIC,numax,Dnu',20000),
        dw4=await get('J/ApJS/233/23/table4','KIC,Teff',20000), leg=await get('J/ApJ/835/172/table1','KIC,numax,Dnu,Teff',2000),
        be=await get('J/AJ/159/280/table2','KIC,Rad,E_Rad,e_Rad,Mass,Teff',300000);
  const BE=new Map(be.map(r=>[r[0],{R:+r[1],M:+r[4]}])), k=1.380649e-23, mp=1.67262e-27, Rs=6.957e8, stars=[], seen=new Set();
  leg.forEach(r=>{stars.push({kic:r[0],dnu:+r[2],T:+r[3],cls:'dwarf'}); seen.add(r[0]);});
  const t4=new Map(dw4.map(r=>[r[0],+r[1]]));
  dw3.forEach(r=>{if(seen.has(r[0])||!t4.get(r[0])) return; stars.push({kic:r[0],dnu:+r[2],T:t4.get(r[0]),cls:'dwarf'}); seen.add(r[0]);});
  gi.forEach(r=>{if(seen.has(r[0])||!(+r[3]>0)||!(+r[1]>0)) return; stars.push({kic:r[0],dnu:+r[3],T:+r[1],cls:r[4]||'giant'}); seen.add(r[0]);});
  const S=[]; stars.forEach(s=>{const b=BE.get(s.kic); if(!b||!(b.R>0)) return; s.R=b.R; s.M=b.M;
    const c=Math.sqrt(5/3*k*s.T/(0.6*mp));            // sound speed of ionised gas at the surface temperature
    s.dcool=c/(2*s.R*Rs)*1e6;                          // ceiling on dnu (muHz) if no layer is hotter than the surface
    s.ratio=s.dnu/s.dcool; s.Tmin=s.T*s.ratio**2;       // lowest peak interior temperature that fits dnu
    S.push(s);});
  return S;   // then: least-squares fits of log dnu on log R (and log M), binned medians, histograms of Tmin
}
