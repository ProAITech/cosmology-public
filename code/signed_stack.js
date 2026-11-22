/* signed_stack.js — Workbench Entry 60
   The through-current test of Entry 57 with each galaxy's spin sense added, so that the pattern can be stacked
   with its sign. Runs in a browser tab after through_current_stack.js (which provides TCS) with the same two
   rotation-measure catalogues in window.D (nvss, lot), plus three kinematic catalogues:
     D.things  THINGS, Walter et al. 2008      J/AJ/136/2563/table1     Name, Dist, BMAG, i, PA
     D.phangs  PHANGS, Lang et al. 2020        J/ApJ/897/122/table1     ID, Dist, PA, i
     D.ghasp   GHASP, Epinat et al. 2008       J/MNRAS/390/466/galaxies Name, Dist, BMAG, i, iK, Vmax, PA-FP, n_PA-FP, f_PA-FP
   Position angles are those of the receding half of the major axis, north through east. THINGS and PHANGS list
   them over 0-360. The GHASP table lists 0-180 with a flag 'a' on those to which 180 must be added; that reading
   is the one that agrees with THINGS and PHANGS for all six galaxies in common. A fourth survey (Gomez-Lopez et al.
   2019) disagreed with GHASP for both galaxies in common and is not used.
   Sign: with +x toward the receding side, a current along the spin gives RM = -I x template. Halo gas: the
   standard halo of Entry 57, TCS.set({N0: 3.68e-2}). */
(function(){
  const num=v=>parseFloat(v),DEG=Math.PI/180;
  const galb=(ra,de)=>{const ag=192.85948*DEG,dg=27.12825*DEG;ra*=DEG;de*=DEG;return Math.asin(Math.sin(de)*Math.sin(dg)+Math.cos(de)*Math.cos(dg)*Math.cos(ra-ag))/DEG;};
  const L=[];
  for(const t of D.things)L.push({name:t.Name.replace(/\s+/g,' ').trim(),src:'THINGS',ra:num(t._RAJ2000),de:num(t._DEJ2000),D:num(t.Dist),i:num(t.i),PA:num(t.PA),big:num(t.BMAG)<=-19});
  for(const p of D.phangs)L.push({name:p.ID.trim(),src:'PHANGS',ra:num(p._RAJ2000),de:num(p._DEJ2000),D:num(p.Dist),i:num(p.i),PA:((num(p.PA)%360)+360)%360,big:true});
  for(const g of D.ghasp){const pa=num(g['PA-FP']);if(!isFinite(pa)||g['f_PA-FP']==='*')continue;const ik=num(g.iK);
    L.push({name:(g.NGC?('NGC '+g.NGC.trim()):g.Name.replace(/\s+/g,' ').trim()),src:'GHASP',ra:num(g._RAJ2000),de:num(g._DEJ2000),D:num(g.Dist),i:isFinite(ik)?ik:num(g.i),PA:(pa+(g['n_PA-FP']==='a'?180:0))%360,big:num(g.Vmax)>=100||num(g.BMAG)<=-19});}
  const U=[];for(const g of L){if(!isFinite(g.ra)||!isFinite(g.D)||!isFinite(g.i)||!isFinite(g.PA))continue;if(U.some(u=>Math.hypot((u.ra-g.ra)*Math.cos(g.de*DEG),u.de-g.de)<0.05))continue;g.b=galb(g.ra,g.de);g.z=g.D*70/299792.458;U.push(g);}
  window.KIN=U;window.GS=U.filter(g=>g.i>=40&&Math.abs(g.b)>=20&&g.D>=1.5&&g.D<=80);window.GB=window.GS.filter(g=>g.big);
  TCS.signed=function(I,flip){let sx=0,sy=0,s2=0;for(const it of I){const s=flip?flip[it.gi]:1;sx+=it.tp*it.tp;sy+=-s*it.tp*it.r;s2+=it.r*it.r;}return {n:I.length,A:sy/sx,formal:Math.sqrt(s2/I.length/sx)};};
  /* null: every galaxy's sense flipped at random, which leaves axes, sight lines and weights as they are */
  TCS.flipNull=function(I,G,nflip,seed){const rnd=TCS.rng(seed||2718),a=[];for(let k=0;k<nflip;k++){const s=G.map(()=>rnd()<0.5?-1:1);a.push(TCS.signed(I,s).A);}const m=a.reduce((x,y)=>x+y,0)/a.length;return Math.sqrt(a.reduce((x,y)=>x+(y-m)**2,0)/a.length);};
})();
/* Published run (11 November 2026). 249 galaxies with a kinematic position angle; 141 with inclination >= 40 deg,
   |b| >= 20 deg and 1.5-80 Mpc; 111 of them large (rotation >= 100 km/s, or blue magnitude <= -19, or in PHANGS):
   16 THINGS, 28 PHANGS, 67 GHASP; median distance 16 Mpc, median inclination 58 deg.
   Conventions checked on shared galaxies: PHANGS-THINGS NGC 2903, 3351, 3521, 3621, 3627 agree within 7 deg
   (NGC 4826 differs by 175 deg: its inner and outer gas counter-rotate); GHASP-THINGS NGC 925, Ho II, NGC 5055,
   NGC 6946 agree within 12 deg; GHASP-PHANGS NGC 2775, NGC 3596 within 32 deg.
   Large spirals, currents in 1e18 A, standard halo, positive = along the spin (error from 2,000 random flips):
     NVSS  30-200 kpc: 405 sight lines, 70 galaxies: +0.122 +- 0.214      15-200: +0.081 +- 0.118
           30-100: -0.033 +- 0.216     30-300: +0.041 +- 0.218
     LoTSS 30-200 kpc:  71 sight lines, 19 galaxies: +0.379 +- 0.267      30-300: +0.317 +- 0.293
   Injection through the NVSS filter: +0.5 -> +0.44, +1 -> +0.90, -1 -> -0.98 (about 0.9 returned); LoTSS in full.
   All 141: NVSS +0.073 +- 0.125, LoTSS +0.153 +- 0.110 (30-200 kpc).
   Sight lines by ring (15-30, 30-60, 60-100, 100-200 kpc): NVSS 10, 34, 66, 305; LoTSS 0, 3, 7, 61.
   Per-galaxy amplitudes, chi2 about zero for galaxies with two or more sight lines: NVSS 56.8 for 55 (random
   angles: median 63.7); LoTSS 29.3 for 13 (median 24.4, 95th percentile 30.7).
   Heaviest galaxies (NVSS): NGC 3031 +0.13 +- 0.48 (15% of the weight), NGC 4449 -0.72 +- 0.63, NGC 6503 +1.38 +- 0.65,
   NGC 4736 +0.50 +- 0.69, NGC 6689 -1.28 +- 0.79, NGC 7741 +0.89 +- 0.87, NGC 2403 +2.01 +- 0.95. */
