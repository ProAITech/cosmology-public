/* gas_vs_stars.js — Workbench Entry 20
   If magnetic forces, which act on gas but not on stars, supplied part of the centripetal force that keeps the
   Galaxy's gas on its orbit, stars would leave the gas behind as soon as they formed: they would stream outward,
   lag the gas, and scatter. Compare gas (masers in star-forming regions, VLBI parallaxes: Reid+ 2019,
   VizieR J/ApJ/885/131/table1) with star clusters of known age (Hunt & Reffert 2023, J/A+A/673/A114/clusters).
   Runs in a browser tab; window.D.ms and window.D.cl are arrays of row objects parsed from VizieR output. */
window.GVS=(function(){
  const DEG=Math.PI/180, K=4.74047;                       // km/s per (mas/yr x kpc)
  const R0=8.15, U0=11.1, VTOT=248.0, W0=7.25;            // Sun: distance from the centre [kpc]; velocity toward the centre, along rotation, toward the north pole [km/s]
  const USTD=10.27, VSTD=15.32, WSTD=7.74;                // the solar motion built into published V_LSR values
  const AG=[[-0.0548755604162154,-0.8734370902348850,-0.4838350155487132],
            [ 0.4941094278755837,-0.4448296299600112, 0.7469822444972189],
            [-0.8676661490190047,-0.1980763734312015, 0.4559837761750669]];   // ICRS -> Galactic
  const num=s=>{const v=parseFloat(s);return isFinite(v)?v:NaN;};
  const rot=v=>AG.map(r=>r[0]*v[0]+r[1]*v[1]+r[2]*v[2]);
  /* Position and velocity in the Galaxy's frame from sky position, distance [kpc], proper motion [mas/yr] and
     radial velocity [km/s]. vlsr=true: the radial velocity is a published V_LSR. Returns R, z [kpc] and
     vR (outward +), vphi (along rotation +), vz [km/s]. */
  function kin(ra,de,d,pmra,pmde,vr,vlsr){
    const a=ra*DEG,dl=de*DEG,ca=Math.cos(a),sa=Math.sin(a),cd=Math.cos(dl),sd=Math.sin(dl);
    const r=[cd*ca,cd*sa,sd],p=[-sa,ca,0],q=[-sd*ca,-sd*sa,cd],g=rot(r);
    if(vlsr) vr=vr-(USTD*g[0]+VSTD*g[1]+WSTD*g[2]);        // back to the heliocentric radial velocity
    const v=rot([0,1,2].map(i=>vr*r[i]+K*d*(pmra*p[i]+pmde*q[i])));
    const X=d*g[0]-R0,Y=d*g[1],Z=d*g[2],VX=v[0]+U0,VY=v[1]+VTOT,VZ=v[2]+W0,R=Math.hypot(X,Y);
    return {R,z:Z,vR:(X*VX+Y*VY)/R,vphi:-(X*VY-Y*VX)/R,vz:VZ};
  }
  function masers(o){
    o=Object.assign({snr:5,Rmin:5,Rmax:14},o||{});const out=[];
    for(const r of window.D.ms){const plx=num(r.plx),e=num(r.e_plx);if(!(plx>0)||plx/e<o.snr)continue;
      const k=kin(num(r._RAJ2000),num(r._DEJ2000),1/plx,num(r.pmE),num(r.pmN),num(r.VLSR),true);
      if(k.R<o.Rmin||k.R>o.Rmax||!isFinite(k.vphi))continue;out.push(Object.assign(k,{name:r.Name}));}
    return out;
  }
  function clusters(o){
    o=Object.assign({cst:5,cmd:0.5,nrv:3,erv:5,zmax:0.5,Rmin:5,Rmax:14},o||{});const out=[];
    for(const r of window.D.cl){ if(r.Type!=='o')continue;
      if(!(num(r.CST)>=o.cst)||!(num(r.CMDCl50)>=o.cmd)||!(num(r.n_RV)>=o.nrv)||!(num(r.e_RV)<=o.erv))continue;
      const d=num(r.dist50)/1000,age=Math.pow(10,num(r.logAge50))/1e6;if(!(d>0)||!isFinite(age))continue;
      const k=kin(num(r.RA_ICRS),num(r.DE_ICRS),d,num(r.pmRA),num(r.pmDE),num(r.RV),false);
      if(k.R<o.Rmin||k.R>o.Rmax||Math.abs(k.z)>o.zmax||!isFinite(k.vphi))continue;out.push(Object.assign(k,{name:r.Name,age,d}));}
    return out;
  }
  const med=v=>{const s=Float64Array.from(v).sort();const n=s.length;return n?(n%2?s[n>>1]:(s[n/2-1]+s[n/2])/2):NaN;};
  const rob=v=>{const m=med(v),s=1.4826*med(v.map(x=>Math.abs(x-m)));return {n:v.length,med:m,sig:s,err:1.2533*s/Math.sqrt(Math.max(1,v.length))};};
  /* Straight-line fit of the gas rotation speed against radius (least squares, 3-sigma clipped once). */
  function gasCurve(M){
    let pts=M.map(m=>[m.R,m.vphi]);
    for(let it=0;it<2;it++){const n=pts.length,sx=pts.reduce((a,p)=>a+p[0],0),sy=pts.reduce((a,p)=>a+p[1],0),sxx=pts.reduce((a,p)=>a+p[0]*p[0],0),sxy=pts.reduce((a,p)=>a+p[0]*p[1],0);
      const b=(n*sxy-sx*sy)/(n*sxx-sx*sx),a=(sy-b*sx)/n,res=pts.map(p=>p[1]-a-b*p[0]),s=1.4826*med(res.map(Math.abs));
      if(it===0) pts=pts.filter((p,i)=>Math.abs(res[i])<3*s); else return {a,b,sig:s,n,at:R=>a+b*R};}
  }
  /* Radial streaming, lag behind the gas, and scatter, in bins of age [Myr]. */
  function byAge(C,curve,edges){
    return edges.slice(0,-1).map((lo,i)=>{const S=C.filter(c=>c.age>=lo&&c.age<edges[i+1]);
      return {lo,hi:edges[i+1],n:S.length,vR:rob(S.map(c=>c.vR)),dphi:rob(S.map(c=>c.vphi-curve.at(c.R))),R:med(S.map(c=>c.R)),age:med(S.map(c=>c.age))};});
  }
  return {kin,masers,clusters,gasCurve,byAge,rob,med,R0,VTOT};
})();

/* Published run (9 November 2026):
   M = GVS.masers()                 129 masers with parallax better than 5 sigma, 5-14 kpc from the centre
   curve = GVS.gasCurve(M)          straight-line gas rotation curve: 232.1 km/s at 8.15 kpc, slope -1.16 km/s/kpc
   C = GVS.clusters()               1,563 open clusters with at least 3 member radial velocities, |z| < 0.5 kpc
   GVS.byAge(C, curve, [0,20,50,120,300,1000,20000])   -> radial streaming, lead over the gas, scatter per age bin
   Predicted values for a share f of the gas's centripetal force that stars do not feel, per cluster of age t at radius R:
     Vg = curve.at(R); Vc = Vg*sqrt(1-f); dv = Vg-Vc; kappa = sqrt(2)*Vc/R
     vR = sqrt(2)*dv*sin(kappa*t);  vphi - Vg = dv*(cos(kappa*t) - 1)        (ages smeared by 0.2 dex, medians per bin)
   Cepheids: VizieR J/A+A/674/A37/table2 (3,306 stars; Dist in kpc, logAge), positions matched within 1 arcsec to
   gaia.dr3lite at dc.g-vo.org for pmra, pmdec and radial_velocity (all 3,306 source ids agree; 2,123 have a
   radial velocity). Kept: radial-velocity error < 10 km/s, RUWE < 1.4, distance-modulus error <= 0.25, |z| < 1.5 kpc:
   1,901 stars, then GVS.kin(...) and the same statistics in bins of radius. */
