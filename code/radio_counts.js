/* Workbench Entry 83: radio source counts and a universe that does not change.
   Data: NVSS 1.4 GHz catalogue (Condon et al. 1998), VizieR VIII/65/nvss on the CfA mirror, queried 11 November 2026.
   Theorem used: for sources spread evenly and not changing, each luminosity contributes N(>S) proportional to r_max^3,
   where S = L f(r) / (4 pi r^2) with any extra dimming f(r) <= 1 that does not increase with r (redshift energy loss,
   curvature, expansion). Then |d ln r_max / d ln S| <= 1/2, so N(>S) never rises faster than S^-1.5 toward fainter S,
   and N(>S) S^1.5 is non-decreasing in S. The sum over luminosities keeps the property.                          */
// all components brighter than 20 mJy over the NVSS sky (dec > -40), and all components of 2.5-20 mJy within 25 deg of the
// north Galactic pole
const Q_BRIGHT = '/viz-bin/asu-tsv?-source=VIII/65/nvss&-out=_RAJ,_DEJ,S1.4&S1.4=%3E20&-out.max=unlimited&-oc.form=dec';
const Q_CONE = '/viz-bin/asu-tsv?-source=VIII/65/nvss&-c=' + encodeURIComponent('192.859 +27.128') + '&-c.rd=25&-out=_RAJ,_DEJ,S1.4&S1.4=2.5..20&-out.max=unlimited&-oc.form=dec';
const RC = (function(){
  const d2r=Math.PI/180, aG=192.85948*d2r, dG=27.12825*d2r;
  const glat=(ra,de)=>{ra*=d2r;de*=d2r;return Math.asin(Math.sin(de)*Math.sin(dG)+Math.cos(de)*Math.cos(dG)*Math.cos(ra-aG))/d2r;};
  const sep=(a1,d1,a2,d2)=>{const s=Math.sin((d2-d1)*d2r/2)**2+Math.cos(d1*d2r)*Math.cos(d2*d2r)*Math.sin((a2-a1)*d2r/2)**2;return 2*Math.asin(Math.sqrt(s))/d2r;};
  // friends-of-friends: components within `link` degrees are joined and their fluxes summed (NVSS splits extended sources)
  const fof=(rows,link)=>{const cell=Math.max(link,0.05),G=new Map(),key=(a,d)=>Math.floor(d/cell)+':'+Math.floor(a*Math.cos(d*d2r)/cell);
    rows.forEach((r,i)=>{const k=key(r[0],r[1]);if(!G.has(k))G.set(k,[]);G.get(k).push(i);});
    const par=rows.map((_,i)=>i),find=i=>{while(par[i]!==i){par[i]=par[par[i]];i=par[i];}return i;};
    rows.forEach((r,i)=>{const di=Math.floor(r[1]/cell),ai=Math.floor(r[0]*Math.cos(r[1]*d2r)/cell);
      for(let dd=-1;dd<=1;dd++)for(let da=-2;da<=2;da++){const L=G.get((di+dd)+':'+(ai+da));if(!L)continue;for(const j of L){if(j<=i)continue;const s=rows[j];if(Math.abs(s[1]-r[1])>link)continue;if(sep(r[0],r[1],s[0],s[1])<=link){const a=find(i),b=find(j);if(a!==b)par[a]=b;}}}});
    const M=new Map();rows.forEach((r,i)=>{const p=find(i);if(!M.has(p))M.set(p,[0,0,0,0]);const m=M.get(p);m[0]+=r[0]*r[2];m[1]+=r[1]*r[2];m[2]+=r[2];m[3]++;});
    return [...M.values()].map(m=>[m[0]/m[2],m[1]/m[2],m[2],m[3]]);};
  const area=(bcut,decmin)=>{let A=0;const n=720,m=360;for(let i=0;i<m;i++){const de=-90+(i+0.5)*180/m;if(de<decmin)continue;const w=Math.cos(de*d2r)*(Math.PI/m)*(2*Math.PI/n);for(let j=0;j<n;j++){const ra=(j+0.5)*360/n;if(Math.abs(glat(ra,de))>bcut)A+=w;}}return A;};
  return {glat,sep,fof,area};})();
/* Steps reported in Entry 83:
   sky = bright rows with dec > -40 and |b| > 10 deg (area 8.665 sr); merged with links of 1', 2' (main) and 3'
   cone = faint rows plus bright rows within 25 deg of the north Galactic pole (0.589 sr), unmerged
   integral counts N(>S) S^1.5 (S in Jy, per sr) at edges 5 mJy x 10^(0.2 k); differential S^2.5 dN/dS in the same bins
   slope of N(>S) between flux pairs: ln(N(>b)/N(>a)) / ln(a/b), Poisson error sqrt(1/Na + 1/Nb)/ln(a/b),
   for the merged and unmerged sky and for the northern and southern Galactic hemispheres separately                */
