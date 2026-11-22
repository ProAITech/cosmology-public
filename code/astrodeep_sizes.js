/* astrodeep_sizes.js — Workbench Entry 59
   The size-redshift result of Raikov, Tsymbal & Lovyagin (arXiv:2507.19651) re-derived from the catalogue they used:
   ASTRODEEP-JWST (Merlin et al. 2024), VizieR J/A+A/691/A240, seven fields. Runs in a browser tab on the VizieR
   site; tables are fetched as tab-separated text and held in page memory only.
   Size = SExtractor FLUX_RADIUS on the F356W+F444W detection stack (column r50SE), in pixels of 0.03 arcsec.
   It is not corrected for the telescope's resolution: stars in the same tables give the size of a point source. */
(async function(){
  const tsv=async(tab,cols,extra)=>{const u='/viz-bin/asu-tsv?-source='+encodeURIComponent('J/A+A/691/A240/'+tab)+'&-out='+cols+'&-out.max=unlimited'+(extra||'');
    const x=await (await fetch(u)).text();const L=x.split('\n').filter(l=>l.length&&l[0]!=='#');const h=L[0].split('\t');const rows=[];for(let i=3;i<L.length;i++){const c=L[i].split('\t');if(c.length>=h.length)rows.push(c);}return {h,rows};};
  const AD={};
  for(const f of ['ceers','a2744','jadesgn','jadesgs','ngdeep','primerc','primeru']){
    const p=await tsv(f+'p','ID,zspec,zphot,Flag','&zspec='+encodeURIComponent('>0'));
    const o=await tsv(f+'o','ID,r50SE,ClassStarSE,FlagsSE,FautoSE,FF444W,FF356W,FF277W,FF200W,FF150W,FF115W');
    const m=new Map();for(const c of o.rows)m.set(c[0].trim(),c);const spec=[];
    for(const c of p.rows){const q=m.get(c[0].trim());if(q)spec.push({f,id:+c[0],zs:+c[1],zp:+c[2],flag:+c[3],r50:+q[1],cs:+q[2],fauto:+q[4],f444:+q[5]});}
    AD[f]={spec,stars:o.rows.filter(c=>+c[2]>0.95&&+c[5]>1).map(c=>+c[1])};}
  window.AD=AD;window.SPEC=[].concat(...Object.values(AD).map(a=>a.spec));
  const PIX=0.03,FLOOR=3.05,ARC=Math.PI/180/3600,cH=299792.458/70;
  const NZ=4000,ZM=12,dc=new Float64Array(NZ+1);{let s=0,prev=1;for(let i=1;i<=NZ;i++){const z=i*ZM/NZ,e=1/Math.sqrt(0.3*Math.pow(1+z,3)+0.7);s+=0.5*(prev+e)*ZM/NZ;prev=e;dc[i]=s;}}
  const DC=z=>{const x=Math.min(z,ZM)/ZM*NZ,i=Math.min(NZ-1,Math.floor(x));return cH*(dc[i]+(dc[i+1]-dc[i])*(x-i));};
  const dA={standard:z=>DC(z)/(1+z),tired:z=>cH*Math.log(1+z),lerner:z=>cH*z};                    /* Mpc */
  const q=(a,p)=>{const s=a.slice().sort((x,y)=>x-y);return s[Math.min(s.length-1,Math.floor(p*s.length))];};
  const S=window.SPEC.filter(s=>s.flag<80&&s.r50>0&&s.zs>0.02&&s.cs<0.9);                          /* 6,316 galaxies */
  const intr=s=>Math.sqrt(Math.max(0,s.r50*s.r50-FLOOR*FLOOR));
  const edges=[0.02,0.5,1,1.5,2,2.5,3,4,5,6,7,8,10],bins=[];
  for(let i=0;i<edges.length-1;i++){const G=S.filter(s=>s.zs>=edges[i]&&s.zs<edges[i+1]);if(G.length<5)continue;const row={z:q(G.map(s=>s.zs),0.5),n:G.length,angle:q(G.map(s=>2*s.r50*PIX),0.5),near:G.filter(s=>s.r50<1.5*FLOOR).length/G.length};
    for(const k in dA){row[k]=q(G.map(s=>2*s.r50*PIX*ARC*dA[k](s.zs)*1000),0.5);row[k+'_deconvolved']=q(G.map(s=>2*intr(s)*PIX*ARC*dA[k](s.zs)*1000),0.5);}bins.push(row);}
  function slope(G,fn){let n=0,sx=0,sy=0,sxx=0,sxy=0;for(const s of G){const y=fn(s);if(!(y>0))continue;const x=Math.log10(1+s.zs),ly=Math.log10(y);n++;sx+=x;sy+=ly;sxx+=x*x;sxy+=x*ly;}return (n*sxy-sx*sy)/(n*sxx-sx*sx);}
  const fits={};for(const zmin of [0,1,2,3]){const G=S.filter(s=>s.zs>=zmin);fits['z'+zmin]={n:G.length,angle:slope(G,s=>s.r50)};for(const k in dA){fits['z'+zmin][k]=slope(G,s=>s.r50*dA[k](s.zs));fits['z'+zmin][k+'_deconvolved']=slope(G,s=>intr(s)*dA[k](s.zs));}}
  window.R1={n:S.length,bins,fits,starFloorPix:q([].concat(...Object.values(AD).map(a=>a.stars)),0.25)};
  return JSON.stringify(window.R1);
})();
/* Published run (11 November 2026): 20,850 objects with a spectroscopic redshift, 7,069 with flag < 80, 6,316 after
   dropping stars and redshifts below 0.02. Stars brighter than 1 in F444W: r50 = 2.97 / 3.10 / 3.35 pixels at the
   5th / 25th / 50th percentile; 3.05 is used as the point-source size (0.092 arcsec).
   Median r50 by redshift bin: 9.86, 8.50, 7.70, 7.08, 6.77, 5.71, 5.44, 4.88, 4.87, 4.57, 4.43, 3.83 pixels.
   Diameter proportional to (1+z)^a: standard -0.38 (all), -0.77 (z >= 1), -0.91 (z >= 2); tired light +0.67, +0.33, +0.22;
   with the point-source size removed in quadrature, z >= 2: standard -1.09, tired light +0.04.
   The papers report -0.48, -0.87 and -1.00 for the standard model. */
