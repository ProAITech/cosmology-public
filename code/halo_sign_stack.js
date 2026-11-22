/* halo_sign_stack.js — Workbench Entry 19
   Does the rotation measure of background radio sources seen through the halos of nearby inclined galaxies
   carry the four-quadrant sign pattern that outward axial currents would give (Myserlis & Contopoulos 2021)?
   With the galaxy's major axis horizontal on the sky (north up, east left, rotated without mirroring),
   outward current along both poles means RM > 0 in the upper-left and lower-right quadrants and RM < 0 in the
   other two:  RM = -A sin(2 theta), A > 0, theta measured from the major axis.

   Runs in a browser tab. Inputs in window.D, parsed from VizieR tab-separated output into arrays of row objects:
     D.lot   LoTSS-DR2 RM grid (O'Sullivan+ 2023; J/MNRAS/519/5723/rmcat): _RAJ2000,_DEJ2000,RM,e_RM,RRM,GRM,zbest
             RRM = RM minus the Galactic foreground of Hutschenreuter+ 2022, as supplied in the catalogue
     D.nvss  NVSS RM catalogue (Taylor+ 2009; J/ApJ/702/1230/catalog): _RAJ2000,_DEJ2000,_Glat,RM,e_RM
     D.rc3   RC3 (de Vaucouleurs+ 1991; VII/155/rc3): _RAJ2000,_DEJ2000,name,T,D25,R25,PA,BT,VGSR,V3K,GLAT
     D.nbg   Tully, Nearby Galaxies Catalog (VII/145/catalog): _RAJ2000,_DEJ2000,R [Mpc]
*/
window.HSS=(function(){
  const DEG=Math.PI/180;
  const num=s=>{const v=parseFloat(s);return isFinite(v)?v:NaN;};
  const incl=(logR25)=>{const q=Math.pow(10,-logR25),q0=0.2;const c2=(q*q-q0*q0)/(1-q0*q0);return c2<=0?90:Math.acos(Math.sqrt(c2))/DEG;};

  /* Disc galaxies, bright (B_T <= BTmax), inclined, with a measured position angle, out of the Milky Way's plane. */
  function galaxies(o){
    o=Object.assign({BTmax:12.5,imin:55,Dmin:2,Dmax:110,bmin:20,Tmin:0},o||{});
    const nb=window.D.nbg.map(r=>({ra:num(r._RAJ2000),de:num(r._DEJ2000),R:num(r.R)})).filter(r=>r.R>0);
    const out=[];
    for(const r of window.D.rc3){
      const ra=num(r._RAJ2000),de=num(r._DEJ2000),PA=num(r.PA),lr=num(r.R25),BT=num(r.BT),T=num(r.T),gb=num(r.GLAT);
      if(!isFinite(ra)||!isFinite(de)||!isFinite(PA)||!isFinite(lr)||!isFinite(BT)||!isFinite(T)) continue;
      if(BT>o.BTmax||T<o.Tmin||Math.abs(gb)<o.bmin) continue;
      const i=incl(lr); if(i<o.imin) continue;
      let D=NaN,src='';
      for(const k of nb){ if(Math.abs(k.de-de)<0.04&&Math.abs(k.ra-ra)*Math.cos(de*DEG)<0.04){D=k.R;src='Tully';break;} }
      const v=isFinite(num(r.V3K))?num(r.V3K):num(r.VGSR);
      if(!(D>0)&&v>300){D=v/70;src='Hubble';}
      if(!(D>=o.Dmin&&D<=o.Dmax)) continue;
      out.push({name:(r.name||r.altname||('PGC'+r.PGC)).replace(/\s+/g,' ').trim(),ra,de,PA,i,D,BT,T,src,z:isFinite(v)?v/299792.458:D*70/299792.458});
    }
    return out;
  }
  /* Background sightlines. LoTSS: residual RM from the catalogue. NVSS: RM minus the median of neighbours within rf degrees. */
  function sightlines(which,o){
    o=Object.assign({rf:3,bmin:20},o||{});
    if(which==='lotss') return window.D.lot.map(r=>({ra:num(r._RAJ2000),de:num(r._DEJ2000),rrm:num(r.RRM),rm:num(r.RM),e:num(r.e_RM),z:num(r.zbest)})).filter(s=>isFinite(s.rrm));
    const S=window.D.nvss.map(r=>({ra:num(r._RAJ2000),de:num(r._DEJ2000),gb:num(r._Glat),rm:num(r.RM),e:num(r.e_RM),z:NaN})).filter(s=>isFinite(s.rm)&&Math.abs(s.gb)>=o.bmin);
    S.sort((a,b)=>a.de-b.de); const de=S.map(s=>s.de);
    const lower=v=>{let lo=0,hi=de.length;while(lo<hi){const m=(lo+hi)>>1;if(de[m]<v)lo=m+1;else hi=m;}return lo;};
    for(let k=0;k<S.length;k++){const s=S[k],c=Math.cos(s.de*DEG),v=[];
      for(let j=lower(s.de-o.rf);j<S.length&&de[j]<=s.de+o.rf;j++){ if(j===k) continue; let da=Math.abs(S[j].ra-s.ra); if(da>180) da=360-da; da*=c; if(da>o.rf) continue; if(da*da+(de[j]-s.de)**2<o.rf*o.rf) v.push(S[j].rm); }
      if(v.length<8){s.rrm=NaN;continue;} v.sort((a,b)=>a-b); s.rrm=s.rm-v[v.length>>1]; }
    return S.filter(s=>isFinite(s.rrm));
  }
  /* All galaxy-sightline pairs with impact parameter below bmax [kpc]. theta is measured from the major axis
     in the frame described above; s2 = sin(2 theta). paFn lets the controls substitute another position angle. */
  function pairs(G,S,bmax,paFn){
    const out=[];
    G.forEach((g,gi)=>{
      const tmax=bmax/(g.D*1000)/DEG, c=Math.cos(g.de*DEG), PA=(paFn?paFn(g,gi):g.PA)*DEG, sp=Math.sin(PA), cp=Math.cos(PA);
      for(let k=0;k<S.length;k++){const s=S[k];
        const dd=s.de-g.de; if(Math.abs(dd)>tmax) continue; let da=s.ra-g.ra; if(da>180)da-=360; if(da<-180)da+=360; const xi=da*c;
        const th2=xi*xi+dd*dd; if(th2>tmax*tmax||th2<(1/60)**2) continue;            // skip the innermost arcminute: the galaxy's own emission
        if(isFinite(s.z)&&s.z<g.z+0.003) continue;                                    // must lie behind the galaxy
        const x=xi*sp+dd*cp, y=xi*cp-dd*sp;                                          // major-axis and minor-axis coordinates
        out.push({gi,k,b:Math.sqrt(th2)*DEG*g.D*1000,x,y,s2:2*x*y/(x*x+y*y),rrm:s.rrm,e:s.e});
      }
    });
    return out;
  }
  /* Amplitude of RRM = -A sin(2 theta): least squares, each sightline weighted equally (the scatter between
     sightlines, not the measurement error, sets the noise). sgn: the sign-only version, mean of -RRM*sign(xy). */
  function fit(P){
    let sxx=0,sxy=0,ss=0,n=0,sg=0; for(const p of P){sxx+=p.s2*p.s2;sxy+=-p.rrm*p.s2;n++;sg+=-p.rrm*Math.sign(p.s2);ss+=p.rrm*p.rrm;}
    const A=sxy/sxx; let res=0; for(const p of P) res+=(p.rrm+A*p.s2)**2;
    return {n,A,err:Math.sqrt(res/Math.max(1,n-1)/sxx),sgn:sg/n,rms:Math.sqrt(ss/n)};
  }
  /* Null distribution: the same sightlines with each galaxy's position angle replaced by a random one. */
  function nullA(G,S,bmax,nreal,seed){
    let s=seed||12345; const rnd=()=>{s=(s*16807)%2147483647;return s/2147483647;}; const A=[];
    for(let r=0;r<nreal;r++){const pa=G.map(()=>rnd()*180); const f=fit(pairs(G,S,bmax,(g,gi)=>pa[gi])); if(isFinite(f.A)) A.push(f.A);}
    const m=A.reduce((a,b)=>a+b,0)/A.length; return {mean:m,sd:Math.sqrt(A.reduce((a,b)=>a+(b-m)**2,0)/A.length),n:A.length,A};
  }
  /* Faster path used for the published numbers: find candidate pairs once (they do not depend on the position
     angle), then re-project them for the true angles, for mirrored angles (180 - PA), or for random ones. */
  function cands(G,S,bmax){const out=[];G.forEach((g,gi)=>{const tmax=bmax/(g.D*1000)/DEG,c=Math.cos(g.de*DEG);
    for(let k=0;k<S.length;k++){const s=S[k],dd=s.de-g.de;if(Math.abs(dd)>tmax)continue;let da=s.ra-g.ra;if(da>180)da-=360;if(da<-180)da+=360;
      const xi=da*c,th2=xi*xi+dd*dd;if(th2>tmax*tmax||th2<(1/60)**2)continue;if(isFinite(s.z)&&s.z<g.z+0.003)continue;
      out.push({gi,k,xi,dd,b:Math.sqrt(th2)*DEG*g.D*1000,rrm:s.rrm});}});return out;}
  function project(C,G,paFn){return C.map(c=>{const PA=(paFn?paFn(G[c.gi],c.gi):G[c.gi].PA)*DEG,x=c.xi*Math.sin(PA)+c.dd*Math.cos(PA),y=c.xi*Math.cos(PA)-c.dd*Math.sin(PA);
    return {gi:c.gi,k:c.k,b:c.b,x,y,s2:2*x*y/(x*x+y*y),rrm:c.rrm};});}
  function nullC(C,G,nreal,seed){let s=seed||12345;const rnd=()=>{s=(s*16807)%2147483647;return s/2147483647;};const A=[];
    for(let r=0;r<nreal;r++){const pa=G.map(()=>rnd()*180);const f=fit(project(C,G,(g,gi)=>pa[gi]));if(isFinite(f.A))A.push(f.A);}
    const m=A.reduce((a,b)=>a+b,0)/A.length;return {mean:m,sd:Math.sqrt(A.reduce((a,b)=>a+(b-m)**2,0)/A.length),n:A.length};}
  return {galaxies,sightlines,pairs,fit,nullA,cands,project,nullC,incl};
})();

/* Published run (9 November 2026):
   G = HSS.galaxies()                        228 disc galaxies: B_T <= 12.5, inclination >= 55 deg, 2-110 Mpc, |b| >= 20 deg
   LoTSS: S = HSS.sightlines('lotss');  C = HSS.cands(G,S,200) keeping |RRM| < 30 rad/m2
   NVSS:  S = HSS.sightlines('nvss');   C = HSS.cands(G,S,200) keeping |RRM| < 50 rad/m2
   amplitude: HSS.fit(HSS.project(C.filter(c => c.b < bmax), G)).A ; error = larger of the formal one and the
   scatter over 400 sets of random position angles (HSS.nullC). Mirror control: HSS.project(C, G, g => 180 - g.PA). */
