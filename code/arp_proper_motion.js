/* arp_proper_motion.js — Workbench Entry 17
   Test of the claim that quasars are ejected from nearby galaxies (Arp; Burbidge; Bell).
   If a quasar left a galaxy a time t ago and now sits an angle theta from it, it must drift
   outward on the sky at mu = theta / t, whatever the galaxy's distance. Gaia measures mu.

   Runs in a browser tab. Inputs, loaded into window.D beforehand:
     D.mq   Million Quasars catalogue v8 (VizieR VII/294): ra, de [deg], z, ty (type), rz (redshift source)
     D.g    Gaia DR3 (GAVO gaia.dr3lite), nearest source within 0.6 arcsec of each Milliquas row:
            pmra, pmde, ea, ed [mas/yr], plx, eplx [mas], G, ruwe, aps, sep [arcsec], nm (number of matches)
     D.nbg  Tully, Nearby Galaxies Catalog (VizieR VII/145): Name, _RAJ2000, _DEJ2000, R [Mpc], MB
     D.ungc Karachentsev+ 2013, Updated Nearby Galaxy Catalog (VizieR J/AJ/145/101): Name, _RAJ2000, _DEJ2000, Dist [Mpc]
*/
window.ARP=(function(){
  const DEG=Math.PI/180, MASRAD=206264806.2, KV=4.740; // km/s per (mas/yr x kpc)
  const BINS=[0,50,100,200,500,1000];                  // projected distance from the galaxy, kpc
  const unit=(ra,de)=>{const a=ra*DEG,d=de*DEG;return [Math.cos(d)*Math.cos(a),Math.cos(d)*Math.sin(a),Math.sin(d)];};
  const NGP=unit(192.85948,27.12825);
  const glat=(ra,de)=>{const u=unit(ra,de);return Math.asin(u[0]*NGP[0]+u[1]*NGP[1]+u[2]*NGP[2])/DEG;};
  const sepDeg=(a,b)=>Math.acos(Math.min(1,a[0]*b[0]+a[1]*b[1]+a[2]*b[2]))/DEG;
  const mad=(v)=>{const s=Float64Array.from(v).sort();const m=s[s.length>>1];const d=s.map(x=>Math.abs(x-m)).sort();return 1.4826*d[d.length>>1];};

  /* Parent galaxies: Tully's catalogue, luminous (M_B <= MBmax), away from the Milky Way's plane,
     beyond the Local Group. Distance from Karachentsev where the galaxy is in that catalogue. */
  function galaxies(o){
    o=Object.assign({MBmax:-19.0,Dmin:2.0,Dmax:40,bmin:20},o||{});
    const U=window.D.ungc.map(r=>({u:unit(+r._RAJ2000,+r._DEJ2000),D:+r.Dist})).filter(r=>isFinite(r.D));
    const out=[];
    for(const r of window.D.nbg){
      const ra=+r._RAJ2000,de=+r._DEJ2000,R0=+r.R;let MB=+r.MB,D=R0;
      if(!isFinite(ra)||!isFinite(de)||!isFinite(MB)||!(D>0)||r.MB==='') continue;
      const u=unit(ra,de);let src='Tully';
      for(const k of U){ if(Math.abs(k.u[2]-u[2])<0.002&&sepDeg(k.u,u)<0.1){D=k.D;src='Karachentsev';break;} }
      MB-=5*Math.log10(D/R0);   // luminosity rescaled to the adopted distance
      const b=glat(ra,de);
      if(MB>o.MBmax||D<o.Dmin||D>o.Dmax||Math.abs(b)<o.bmin) continue;
      out.push({name:r.Name.replace(/\s+/g,' '),ra,de,D,MB,b,u,src});
    }
    return out;
  }

  /* Quasars: confirmed broad-line quasars (type contains Q) with a spectroscopic redshift >= zmin and one
     clean Gaia match. No cut on proper motion in the primary sample. Errors scaled by the measured
     scatter of the whole sample in each magnitude band. */
  function quasars(o){
    o=Object.assign({zmin:0.5,ruwe:1.4,sep:0.5,plxsig:5,clip:0,src:null},o||{});
    const m=window.D.mq,g=window.D.g,N=m.ra.length,bad={GAIA3:1,MQ:1,'':1};
    const keep=[];
    for(let i=0;i<N;i++){
      if(m.ty[i].indexOf('Q')<0||!(m.z[i]>=o.zmin)||bad[m.rz[i]]) continue;
      if(o.src&&!o.src.test(m.rz[i])) continue;
      if(g.nm[i]!==1||!(g.sep[i]<o.sep)||!(g.aps[i]===31||g.aps[i]===95)) continue;
      if(!isFinite(g.pmra[i])||!isFinite(g.pmde[i])||!(g.ea[i]>0)||!(g.ed[i]>0)||!(g.ruwe[i]<o.ruwe)) continue;
      if(Math.abs((g.plx[i]+0.017)/g.eplx[i])>o.plxsig) continue;
      keep.push(i);
    }
    // error scale per G band from the all-sky scatter of pm/error (true quasar motions are ~0)
    const edges=[0,17,18,19,19.5,20,20.5,30],kfac=[];
    for(let b=0;b<edges.length-1;b++){const v=[];for(const i of keep){const G=g.G[i];if(G>=edges[b]&&G<edges[b+1]){v.push(g.pmra[i]/g.ea[i]);v.push(g.pmde[i]/g.ed[i]);}}kfac.push(v.length>50?Math.max(1,mad(v)):1);}
    const kOf=G=>{for(let b=0;b<edges.length-1;b++) if(G<edges[b+1]) return kfac[b];return kfac[kfac.length-1];};
    let rows=keep.map(i=>{const k=kOf(g.G[i]);return {i,de:m.de[i],ra:m.ra[i],k};});
    let nclip=0;
    if(o.clip>0) rows=rows.filter(r=>{const i=r.i,c=(g.pmra[i]/(r.k*g.ea[i]))**2+(g.pmde[i]/(r.k*g.ed[i]))**2;if(c>o.clip*o.clip){nclip++;return false;}return true;});
    rows.sort((a,b)=>a.de-b.de);
    const n=rows.length,Q={n,kfac,edges,nclip,idx:new Int32Array(n),de:new Float64Array(n),x:new Float64Array(n),y:new Float64Array(n),z:new Float64Array(n),
      sa:new Float64Array(n),ca:new Float64Array(n),sd:new Float64Array(n),cd:new Float64Array(n),pa:new Float64Array(n),pd:new Float64Array(n),ea:new Float64Array(n),ed:new Float64Array(n)};
    rows.forEach((r,j)=>{const i=r.i,a=r.ra*DEG,d=r.de*DEG;Q.idx[j]=i;Q.de[j]=r.de;Q.sa[j]=Math.sin(a);Q.ca[j]=Math.cos(a);Q.sd[j]=Math.sin(d);Q.cd[j]=Math.cos(d);
      Q.x[j]=Q.cd[j]*Q.ca[j];Q.y[j]=Q.cd[j]*Q.sa[j];Q.z[j]=Q.sd[j];Q.pa[j]=g.pmra[i];Q.pd[j]=g.pmde[i];Q.ea[j]=r.k*g.ea[i];Q.ed[j]=r.k*g.ed[i];});
    return Q;
  }
  const lower=(a,v)=>{let lo=0,hi=a.length;while(lo<hi){const m=(lo+hi)>>1;if(a[m]<v)lo=m+1;else hi=m;}return lo;};

  /* For every galaxy-quasar pair inside Rmax: outward (radial) and sideways (tangential) motion.
     each(g, j, theta[rad], rp[kpc], mu_r, sig_r, mu_t, sig_t) */
  function pairs(G,Q,each,Rmax){
    Rmax=Rmax||BINS[BINS.length-1];
    for(let gi=0;gi<G.length;gi++){
      const g=G[gi],Dk=g.D*1000,tmax=Math.min(25*DEG,Math.asin(Math.min(1,Rmax/Dk))),cmax=Math.cos(tmax),u=g.u;
      const lo=lower(Q.de,g.de-tmax/DEG),hi=lower(Q.de,g.de+tmax/DEG);
      for(let j=lo;j<hi;j++){
        const dot=u[0]*Q.x[j]+u[1]*Q.y[j]+u[2]*Q.z[j]; if(dot<cmax||dot>=1) continue;
        const th=Math.acos(dot),s=Math.sin(th); if(s<1e-7) continue;
        // unit vector at the quasar pointing away from the galaxy, in (east, north) components
        const ox=-(u[0]-dot*Q.x[j])/s,oy=-(u[1]-dot*Q.y[j])/s,oz=-(u[2]-dot*Q.z[j])/s;
        const oE=-Q.sa[j]*ox+Q.ca[j]*oy, oN=-Q.sd[j]*Q.ca[j]*ox-Q.sd[j]*Q.sa[j]*oy+Q.cd[j]*oz;
        const mr=Q.pa[j]*oE+Q.pd[j]*oN, mt=-Q.pa[j]*oN+Q.pd[j]*oE;
        const sr=Math.sqrt((Q.ea[j]*oE)**2+(Q.ed[j]*oN)**2), st=Math.sqrt((Q.ea[j]*oN)**2+(Q.ed[j]*oE)**2);
        each(gi,j,th,Dk*s,mr,sr,mt,st);
      }
    }
  }

  /* Stack. Per galaxy and distance bin, inverse-variance sums of:
     r: outward motion [mas/yr]   t: sideways motion [mas/yr] (must be zero in every picture)
     v: outward speed if the quasar is at the galaxy's distance [km/s]
     i: 1/travel-time = mu_r/theta [1/yr] */
  function stack(G,Q){
    const nb=BINS.length-1,A=G.map(()=>Array.from({length:nb},()=>new Float64Array(9)));
    pairs(G,Q,(gi,j,th,rp,mr,sr,mt,st)=>{
      let b=0;while(b<nb&&rp>=BINS[b+1])b++; if(b>=nb) return;
      const a=A[gi][b],Dk=G[gi].D*1000,wr=1/(sr*sr),wt=1/(st*st),sv=KV*Dk*sr,wv=1/(sv*sv),tm=th*MASRAD,wi=tm*tm/(sr*sr);
      a[0]++;a[1]+=wr;a[2]+=wr*mr;a[3]+=wt;a[4]+=wt*mt;a[5]+=wv;a[6]+=wv*KV*Dk*mr;a[7]+=wi;a[8]+=wi*mr/tm;
    });
    return A;
  }
  /* Combine galaxies; errors both formal and by leaving one galaxy out at a time. */
  function summarise(A,sel,b0,b1){
    const idx=[];for(let g=0;g<A.length;g++) if(!sel||sel(g)) idx.push(g);
    const tot=new Float64Array(9),per=idx.map(g=>{const s=new Float64Array(9);for(let b=b0;b<=b1;b++)for(let k=0;k<9;k++)s[k]+=A[g][b][k];for(let k=0;k<9;k++)tot[k]+=s[k];return s;});
    const est=(t,w,x)=>t[w]>0?t[x]/t[w]:NaN, res={pairs:tot[0],gals:per.filter(s=>s[0]>0).length};
    [['r',1,2],['t',3,4],['v',5,6],['i',7,8]].forEach(([nm,w,x])=>{
      const m=est(tot,w,x);let v=0,n=0;
      for(const s of per){ if(s[w]<=0) continue; const mi=(tot[x]-s[x])/(tot[w]-s[w]); if(isFinite(mi)){v+=(mi-m)**2;n++;} }
      res[nm]={mean:m,formal:1/Math.sqrt(tot[w]),jack:n>1?Math.sqrt((n-1)/n*v):NaN};
    });
    return res;
  }
  /* The same galaxies moved in Galactic longitude: a sky with no parent galaxies at these spots. */
  function shifted(G,dl){
    const X=unit(266.40499,-28.93617),Z=NGP,Y=[Z[1]*X[2]-Z[2]*X[1],Z[2]*X[0]-Z[0]*X[2],Z[0]*X[1]-Z[1]*X[0]];
    return G.map(g=>{const u=g.u,gx=u[0]*X[0]+u[1]*X[1]+u[2]*X[2],gy=u[0]*Y[0]+u[1]*Y[1]+u[2]*Y[2],gz=u[0]*Z[0]+u[1]*Z[1]+u[2]*Z[2];
      const l=Math.atan2(gy,gx)+dl*DEG,c=Math.hypot(gx,gy),nx=c*Math.cos(l),ny=c*Math.sin(l);
      const v=[nx*X[0]+ny*Y[0]+gz*Z[0],nx*X[1]+ny*Y[1]+gz*Z[1],nx*X[2]+ny*Y[2]+gz*Z[2]];
      let ra=Math.atan2(v[1],v[0])/DEG;if(ra<0)ra+=360;
      return Object.assign({},g,{u:v,ra,de:Math.asin(v[2])/DEG});});
  }
  /* Every catalogued quasar with a Gaia motion within rad arcmin of a position (for Arp's own cases). */
  function field(ra,de,Dmpc,radArcmin,zmin){
    const m=window.D.mq,g=window.D.g,u=unit(ra,de),out=[],cm=Math.cos(radArcmin/60*DEG);
    for(let i=0;i<m.ra.length;i++){
      if(Math.abs(m.de[i]-de)>radArcmin/60) continue;
      const q=unit(m.ra[i],m.de[i]),dot=u[0]*q[0]+u[1]*q[1]+u[2]*q[2]; if(dot<cm) continue;
      if(!(m.z[i]>=(zmin||0))) continue;
      const th=Math.acos(Math.min(1,dot)),s=Math.sin(th)||1e-12,a=m.ra[i]*DEG,d=m.de[i]*DEG;
      const ox=-(u[0]-dot*q[0])/s,oy=-(u[1]-dot*q[1])/s,oz=-(u[2]-dot*q[2])/s;
      const oE=-Math.sin(a)*ox+Math.cos(a)*oy,oN=-Math.sin(d)*Math.cos(a)*ox-Math.sin(d)*Math.sin(a)*oy+Math.cos(d)*oz;
      const ok=isFinite(g.pmra[i])&&g.sep[i]<0.6;
      out.push({i,z:+m.z[i].toFixed(3),ty:m.ty[i],sepArcmin:+(th/DEG*60).toFixed(2),pa:+((Math.atan2(oE,oN)/DEG+360)%360).toFixed(0),G:ok?+g.G[i].toFixed(2):null,
        mr:ok?g.pmra[i]*oE+g.pmde[i]*oN:null,sr:ok?Math.sqrt((g.ea[i]*oE)**2+(g.ed[i]*oN)**2):null,mt:ok?-g.pmra[i]*oN+g.pmde[i]*oE:null,ruwe:ok?+g.ruwe[i].toFixed(2):null,
        pmra:ok?g.pmra[i]:null,pmde:ok?g.pmde[i]:null,ea:ok?g.ea[i]:null,ed:ok?g.ed[i]:null,Dk:Dmpc*1000,th});
    }
    return out.sort((a,b)=>a.sepArcmin-b.sepArcmin);
  }
  /* Second, independent sample: SDSS DR16Q (VizieR VII/289) with the Gaia DR2 motions listed in that
     catalogue. Its "error" columns are inverse variances. D.dr16: ra, de, z, pa, ea, pd, ed, plx, eplx, G. */
  function quasarsDR16(o){
    o=Object.assign({zmin:0.5,plxsig:5,clip:0},o||{});const d=window.D.dr16,keep=[];
    for(let i=0;i<d.ra.length;i++){
      if(!(d.z[i]>=o.zmin)||!(d.ea[i]>0)||!(d.ed[i]>0)||!isFinite(d.pa[i])||!isFinite(d.pd[i])) continue;
      if(d.eplx[i]>0&&Math.abs(d.plx[i]*Math.sqrt(d.eplx[i]))>o.plxsig) continue;
      keep.push(i);
    }
    const edges=[0,18,19,19.5,20,20.5,30],kfac=[];
    for(let b=0;b<edges.length-1;b++){const v=[];for(const i of keep){const G=d.G[i];if(G>=edges[b]&&G<edges[b+1]){v.push(d.pa[i]*Math.sqrt(d.ea[i]));v.push(d.pd[i]*Math.sqrt(d.ed[i]));}}kfac.push(v.length>50?Math.max(1,mad(v)):1);}
    const kOf=G=>{for(let b=0;b<edges.length-1;b++) if(G<edges[b+1]) return kfac[b];return kfac[kfac.length-1];};
    let rows=keep.map(i=>({i,de:d.de[i],ra:d.ra[i],k:kOf(d.G[i])}));
    if(o.clip>0) rows=rows.filter(r=>{const i=r.i;return (d.pa[i]*Math.sqrt(d.ea[i])/r.k)**2+(d.pd[i]*Math.sqrt(d.ed[i])/r.k)**2<=o.clip*o.clip;});
    rows.sort((a,b)=>a.de-b.de);
    const n=rows.length,Q={n,kfac,edges,idx:new Int32Array(n),de:new Float64Array(n),x:new Float64Array(n),y:new Float64Array(n),z:new Float64Array(n),
      sa:new Float64Array(n),ca:new Float64Array(n),sd:new Float64Array(n),cd:new Float64Array(n),pa:new Float64Array(n),pd:new Float64Array(n),ea:new Float64Array(n),ed:new Float64Array(n)};
    rows.forEach((r,j)=>{const i=r.i,a=r.ra*DEG,dd=r.de*DEG;Q.idx[j]=i;Q.de[j]=r.de;Q.sa[j]=Math.sin(a);Q.ca[j]=Math.cos(a);Q.sd[j]=Math.sin(dd);Q.cd[j]=Math.cos(dd);
      Q.x[j]=Q.cd[j]*Q.ca[j];Q.y[j]=Q.cd[j]*Q.sa[j];Q.z[j]=Q.sd[j];Q.pa[j]=d.pa[i];Q.pd[j]=d.pd[i];Q.ea[j]=r.k/Math.sqrt(d.ea[i]);Q.ed[j]=r.k/Math.sqrt(d.ed[i]);});
    return Q;
  }

  /* The numbers quoted in Entry 17. Error = largest of formal, leave-one-galaxy-out, and the scatter over
     19 copies of the galaxy list shifted in Galactic longitude. Units: micro-arcsec/yr, km/s, 1e-9 per yr. */
  function run(){
    const G=galaxies(),Q=quasars({clip:5}),S=stack(G,Q),sets=[[0,0],[1,1],[2,2],[3,3],[4,4],[0,2],[0,3]];
    const shifts=[-165,-150,-135,-120,-105,-90,-75,-60,-45,45,60,75,90,105,120,135,150,165,180];
    const ctl=sets.map(()=>({r:[],t:[],v:[],i:[]}));
    for(const dl of shifts){const Sx=stack(shifted(G,dl),Q);sets.forEach(([a,b],k)=>{const o=summarise(Sx,null,a,b);for(const q of ['r','t','v','i'])ctl[k][q].push(o[q].mean);});}
    const rms=v=>Math.sqrt(v.reduce((s,x)=>s+x*x,0)/v.length),sc={r:1000,t:1000,v:1,i:1e9};
    return {galaxies:G.length,quasars:Q.n,removedAbove5sigma:Q.nclip,rows:sets.map(([a,b],k)=>{const o=summarise(S,null,a,b),row={kpc:[BINS[a],BINS[b+1]],pairs:o.pairs,galaxies:o.gals};
      for(const q of ['r','t','v','i']) row[q]={mean:o[q].mean*sc[q],err:Math.max(o[q].formal,o[q].jack,rms(ctl[k][q]))*sc[q]};return row;})};
  }
  return {BINS,galaxies,quasars,quasarsDR16,pairs,stack,summarise,shifted,field,run,unit,glat,KV,MASRAD};
})();

/* How the inputs were loaded (9 November 2026), from a tab open on vizier.cfa.harvard.edu:
   D.mq    /viz-bin/asu-tsv?-source=VII/294/catalog&-out=RAJ2000,DEJ2000,Type,z,rz&-out.max=unlimited&-oc.form=dec
   D.nbg   /viz-bin/asu-tsv?-source=VII/145/catalog&-out=Name,_RAJ2000,_DEJ2000,R,MB,D25,RVel,HubCode,GLAT&-out.max=unlimited&-oc.form=dec
   D.ungc  /viz-bin/asu-tsv?-source=J/AJ/145/101/catalog&-out=Name,_RAJ2000,_DEJ2000,Dist,BMag,KLum,a26,TT&-out.max=unlimited&-oc.form=dec
   D.dr16  /viz-bin/asu-tsv?-source=VII/289/dr16q&-out=RAJ2000,DEJ2000,z,pmRA,e_pmRA,pmDE,e_pmDE,Plx,e_Plx,Gmag&-out.max=unlimited&-oc.form=dec
   D.g     Milliquas positions uploaded in blocks of 25,000 to https://dc.g-vo.org/tap/sync (TAP_UPLOAD) with
           SELECT t.id,g.pmra,g.pmdec,g.pmra_error,g.pmdec_error,g.parallax,g.parallax_error,g.phot_g_mean_mag,g.ruwe,
                  g.astrometric_params_solved,g.astrometric_excess_noise,DISTANCE(t.ra,t.dec,g.ra,g.dec)*3600 AS sep
           FROM TAP_UPLOAD.t1 AS t JOIN gaia.dr3lite AS g ON DISTANCE(t.ra,t.dec,g.ra,g.dec)<0.6/3600
           keeping the nearest source and counting matches (nm); 661,767 of 1,021,800 rows matched.
   Then: ARP.run() for the stack; ARP.field(ra,dec,Mpc,arcmin,zmin) for single fields;
   variants by passing options to ARP.quasars / ARP.galaxies or a galaxy filter to ARP.summarise. */
