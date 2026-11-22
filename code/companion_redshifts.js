/* Workbench Entry 78: are companion galaxies redshifted relative to their dominant galaxy (Arp 1994: all 22 major companions in the
   Local Group and M81 group)?  SDSS groups of Tempel et al. (2014), VizieR J/A+A/566/A1/galaxies: GroupID, Ngal, Rank (1 = brightest),
   zobs, e_z, rMAG, pE, pSc.  Run in a browser tab on a VizieR mirror:  await CR.load();  CR.run();
   dv = c (z_i - z_1)/(1 + z_1) for every member i relative to the brightest member; |dv| > 3000 km/s dropped.                        */
window.CR=(function(){const c=299792.458,R={};
 R.load=async function(){const ft=async u=>{for(let k=0;k<5;k++){try{const r=await fetch(u);if(r.ok){const t=await r.text();if(!t.includes('STDB_ERR'))return t;}}catch(e){}await new Promise(r=>setTimeout(r,2000));}return 'FAIL';};
  const t=await ft('/viz-bin/asu-tsv?-source=J/A+A/566/A1/galaxies&-out=GroupID,Ngal,Rank,zobs,e_z,rMAG,pE,pSc&-out.max=unlimited&-oc.form=dec&Ngal=>1');
  R.G=new Map();for(const l of t.split('\n').filter(l=>l&&l[0]!=='#').slice(3)){const r=l.split('\t');if(r.length<8)continue;const g=+r[0],ng=+r[1],z=+r[3];if(!(ng>=2)||!isFinite(z))continue;
   if(!R.G.has(g))R.G.set(g,[]);R.G.get(g).push({rk:+r[2],z,M:+r[5],pE:+r[6],pSc:+r[7]});}return R.G.size;};
 const st=a=>{const n=a.length;if(n<5)return {n};const mu=a.reduce((u,v)=>u+v,0)/n,sd=Math.sqrt(a.reduce((u,v)=>u+(v-mu)**2,0)/n),s=a.slice().sort((x,y)=>x-y);return {n,mean:mu,sem:sd/Math.sqrt(n),sd,median:s[n>>1],fracPos:a.filter(v=>v>0).length/n};};
 R.run=function(){const zb=[0,0.04,0.06,0.08,0.10,0.12,0.15,0.2],sets={all:[],small:[],rich:[],richCore:[],dominant:[],domSpiral:[],domSpiralSmall:[]},byz=zb.slice(1).map(()=>({rich:[],richCore:[],small:[],dsS:[]}));
  for(const [g,mem] of R.G){const p=mem.find(m=>m.rk===1);if(!p||mem.length<2)continue;const s2=mem.find(m=>m.rk===2),gap=s2?s2.M-p.M:0,ng=mem.length;let k=0;while(k<zb.length-2&&p.z>=zb[k+1])k++;
   for(const m of mem){if(m===p)continue;const dv=c*(m.z-p.z)/(1+p.z);if(Math.abs(dv)>3000)continue;sets.all.push(dv);
    if(ng<=4){sets.small.push(dv);byz[k].small.push(dv);}if(ng>=10){sets.rich.push(dv);byz[k].rich.push(dv);if(Math.abs(dv)<300){sets.richCore.push(dv);byz[k].richCore.push(dv);}}
    if(gap>=1.5){sets.dominant.push(dv);if(p.pSc>0.5){sets.domSpiral.push(dv);if(ng<=4){sets.domSpiralSmall.push(dv);byz[k].dsS.push(dv);}}}}}
  const out={};for(const k in sets)out[k]=st(sets[k]);out.byz={zbins:zb,rows:byz.map(b=>({rich:st(b.rich),richCore:st(b.richCore),small:st(b.small),domSpiralSmall:st(b.dsS)}))};return out;};
 return R;})();
