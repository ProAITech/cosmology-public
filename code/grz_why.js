// Entry 153 — browser code run on the VizieR CfA page, on top of Entry 150's grz_stack.js
// (GALS, BCGM, PB.P = [R_Mpc, dv_kms, central_gid, satellite_index], fitEMw, offAt, mcut).
// Colours: rest-frame g - r from J/A+A/566/A1/galaxies (gMAG - rMAG), matched by position.

// 20 sky jackknife zones: quantiles of central RA
const ras=[...new Set(PB.P.map(p=>p[2]))].map(id=>BCGM.get(id).ra).sort((a,b)=>a-b);
const QRA=Array.from({length:21},(_,k)=>ras[Math.min(ras.length-1,Math.floor(k*ras.length/20))]); QRA[20]=361;
const inZ=(p,k)=>{const ra=BCGM.get(p[2]).ra; return ra>=QRA[k]&&ra<QRA[k+1];};
const sd=a=>{const m=a.reduce((s,x)=>s+x,0)/a.length; return Math.sqrt((a.length-1)/a.length*a.reduce((s,x)=>s+(x-m)**2,0));};

// colour table
async function loadColours(){
  const u='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2BA/566/A1/galaxies&-out.max=700000&-out=RAJ2000,DEJ2000,gMAG,rMAG,pE&-out.form=|';
  const t=await fetch(u).then(r=>r.text()); const m=new Map();
  t.split('\n').forEach(l=>{if(!l||l.startsWith('#')||!/^\s*\d/.test(l))return; const s=l.split('|');
    m.set(Math.round(+s[0]*1e5)+'_'+Math.round(+s[1]*1e5),[(+s[2])-(+s[3]),+s[4]]);});
  window.GR=new Float32Array(GALS.length);
  GALS.forEach((g,i)=>{const v=m.get(Math.round(g.ra*1e5)+'_'+Math.round(g.de*1e5)); GR[i]=v?v[0]:NaN;});
}
const red=p=>GR[p[3]]>0.68-0.03*(GALS[p[3]].M+20);

// complete sample, luminosity judged at the central's distance
const C=PB.P.filter(p=>p[0]<4.5&&GALS[p[3]].rm<mcut(BCGM.get(p[2]).z,0.3));
const Mc=p=>GALS[p[3]].rm-offAt(BCGM.get(p[2]).z);

// statistics
const asym=Q=>{let Np=0,Nm=0; for(const p of Q){const v=p[1]; if(v>2000&&v<4000)Np++; else if(v<-2000&&v>-4000)Nm++;} return (Np-Nm)/(Np+Nm);};
const winMean=(Q,V)=>{let Np=0,Nm=0,n=0,sv=0; for(const p of Q){const v=p[1]; if(v>2000&&v<4000)Np++; else if(v<-2000&&v>-4000)Nm++; else if(Math.abs(v)<V){n++; sv+=v;}}
  const a=(Np+Nm)/4000, sl=(Np-Nm)/1.2e7; return (sv-sl*2*V**3/3)/(n-2*V*a);};   // linear background from the wings
const winMedian=(Q,V)=>{const a=Q.map(p=>p[1]).filter(v=>Math.abs(v)<V).sort((x,y)=>x-y); return a[Math.floor(a.length/2)];};
const fit=Q=>fitEMw(Q.map(p=>p[1]),null,true).mu;                                  // Entry 150's two-Gaussian fit

// regression of centre on wing asymmetry, 4 luminosity quarters x 2 redshift halves, repeated per jackknife zone
async function regression(){
  const s=C.map(Mc).sort((a,b)=>a-b); const q=[-99,s[s.length>>2],s[s.length>>1],s[(3*s.length)>>2],99];
  const subs=[]; for(let i=0;i<4;i++) for(const [zl,zh] of [[0,0.08],[0.08,0.3]])
    subs.push(C.filter(p=>{const m=Mc(p),z=BCGM.get(p[2]).z; return m>=q[i]&&m<q[i+1]&&z>=zl&&z<zh;}));
  const res=subs.map(Q=>({A:asym(Q),mu:fit(Q),jA:[],jmu:[]}));
  for(let k=0;k<20;k++){ subs.forEach((Q,i)=>{const Qk=Q.filter(p=>!inZ(p,k)); res[i].jA.push(asym(Qk)); res[i].jmu.push(fit(Qk));}); await new Promise(r=>setTimeout(r,0)); }
  const wls=(x,y,w)=>{let S=0,Sx=0,Sy=0,Sxx=0,Sxy=0; for(let i=0;i<x.length;i++){S+=w[i];Sx+=w[i]*x[i];Sy+=w[i]*y[i];Sxx+=w[i]*x[i]*x[i];Sxy+=w[i]*x[i]*y[i];}
    const D=S*Sxx-Sx*Sx; return [(Sxx*Sy-Sx*Sxy)/D,(S*Sxy-Sx*Sy)/D];};
  const w=res.map(r=>1/sd(r.jmu)**2); const F=wls(res.map(r=>r.A),res.map(r=>r.mu),w);
  const J=[]; for(let k=0;k<20;k++) J.push(wls(res.map(r=>r.jA[k]),res.map(r=>r.jmu[k]),w));
  return {c0:F[0],e0:sd(J.map(j=>j[0])),slope:F[1],es:sd(J.map(j=>j[1]))};
}

// isolated pairs: Tempel groups with exactly two members
function isolatedPairs(){
  const by=new Map(); GALS.forEach((g,i)=>{if(g.gid>0&&g.ng===2){if(!by.has(g.gid))by.set(g.gid,[]); by.get(g.gid).push(i);}});
  const P=[]; by.forEach(a=>{if(a.length!==2)return; let [i,j]=a; if(GALS[i].rk>GALS[j].rk)[i,j]=[j,i];
    const dv=299792.458*(GALS[j].z-GALS[i].z)/(1+GALS[i].z); if(Math.abs(dv)<1000) P.push({dv,rmF:GALS[j].rm,redF:GR[j]>0.68-0.03*(GALS[j].M+20)});});
  return P;
}
