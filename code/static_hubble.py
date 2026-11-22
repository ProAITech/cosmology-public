"""Workbench Entry 31: the static universe's Hubble law against the supernovae.
Lerner, Falomo & Scarpa (2014) and Lerner (2018) use a static, flat (Euclidean) universe in which distance is d = cz/H0 at every redshift, and light
loses energy on the way but nothing else, so a supernova's apparent brightness follows  m - M = 5 log(cz/H0) + 2.5 log(1+z) + 25.
They state that this fits the supernova data without dark energy. This script tests that with the binned supernova distances of the
Joint Light-curve Analysis (Betoule et al. 2014, 740 supernovae; their tables of 31 control points and the full covariance matrix).
Each model gets a free additive constant (the supernova luminosity and H0 are not known separately)."""
import numpy as np, json
from scipy import integrate, optimize
zb=np.array([0.010,0.012,0.014,0.016,0.019,0.023,0.026,0.031,0.037,0.043,0.051,0.060,0.070,0.082,0.097,0.114,0.134,0.158,0.186,0.218,0.257,0.302,0.355,0.418,0.491,0.578,0.679,0.799,0.940,1.105,1.300])
mub=np.array([32.9538,33.8790,33.8421,34.1185,34.5934,34.9390,35.2520,35.7485,36.0697,36.4345,36.6511,37.1580,37.4301,37.9566,38.2532,38.6128,39.0678,39.3414,39.7921,40.1565,40.5649,40.9052,41.4214,41.7909,42.2314,42.6170,43.0527,43.5041,43.9725,44.5140,44.8218])
tri="""21282 -10840 1918 451 946 614 785 686 581 233 881 133 475 295 277 282 412 293 337 278 219 297 156 235 133 179 -25 -106 0 137 168
28155 -2217 1702 74 322 380 273 424 487 266 303 406 468 447 398 464 403 455 468 417 444 351 399 83 167 -86 15 -2 76 243
6162 -1593 1463 419 715 580 664 465 613 268 570 376 405 352 456 340 412 355 317 341 242 289 119 152 -69 -33 -44 37 209
5235 -722 776 588 591 583 403 651 212 555 353 355 323 442 319 372 337 288 343 210 272 92 167 -48 -29 -21 50 229
7303 -508 1026 514 596 315 621 247 493 320 375 290 383 286 350 300 269 313 198 251 99 126 18 46 13 10 203
3150 -249 800 431 358 414 173 514 231 248 221 293 187 245 198 175 231 126 210 103 170 51 66 -8 -51 308
3729 -88 730 321 592 188 546 316 342 290 389 267 341 285 252 301 189 242 122 159 35 72 30 28 255
3222 -143 568 421 203 491 257 280 240 301 221 275 227 210 249 148 220 123 160 43 69 27 7 253
3225 -508 774 156 502 273 323 276 370 260 316 273 231 273 171 226 111 154 0 29 19 23 206
5646 -1735 691 295 362 316 305 370 280 346 313 276 310 217 274 131 175 38 118 78 48 303
8630 -1642 944 152 253 184 274 202 254 233 196 237 156 207 27 115 -32 7 -15 0 176
3855 -754 502 225 278 294 274 285 253 239 255 173 229 181 177 93 124 132 108 227
4340 -634 660 240 411 256 326 276 235 290 184 256 135 222 90 152 67 17 318
2986 -514 479 340 363 377 362 315 343 265 311 144 198 17 62 86 147 226
3592 -134 606 333 422 374 333 349 267 300 157 184 9 71 85 136 202
1401 22 431 343 349 302 322 245 284 171 186 70 70 93 142 202
1491 141 506 386 356 394 278 306 188 212 79 71 106 145 240
1203 200 435 331 379 281 311 184 209 49 51 110 197 181
1032 258 408 398 305 330 197 223 78 79 113 174 225
1086 232 453 298 328 120 189 -48 22 42 142 204
1006 151 329 282 169 195 58 80 95 192 188
1541 124 400 199 261 150 166 202 251 251
1127 72 227 222 93 118 93 171 161
1723 -105 406 -3 180 190 198 247
1550 144 946 502 647 437 215
1292 187 524 393 387 284
3941 587 1657 641 346
2980 360 1124 305
4465 -1891 713
23902 -1826
19169"""
C=np.zeros((31,31))
for i,line in enumerate(tri.split("\n")):
    v=[float(x) for x in line.split()]; assert len(v)==31-i,(i,len(v)); C[i,i:]=v; C[i:,i]=v
C*=1e-6; assert np.all(np.linalg.eigvalsh(C)>0); Ci=np.linalg.inv(C); one=np.ones(31)
def chi2(mu_model):
    d=mub-mu_model; a=(one@Ci@d)/(one@Ci@one); r=d-a; return float(r@Ci@r),float(a),r
mu=lambda dl: 5*np.log10(dl)
E=lambda z,om: np.sqrt(om*(1+z)**3+1-om)
dl_lcdm=lambda om: np.array([(1+z)*integrate.quad(lambda x:1/E(x,om),0,z)[0] for z in zb])
fit=optimize.minimize_scalar(lambda om: chi2(mu(dl_lcdm(om)))[0],bounds=(0.05,0.9),method="bounded"); om=fit.x
# error on Omega_m from delta chi2 = 1
f=lambda o: chi2(mu(dl_lcdm(o)))[0]-fit.fun-1; om_lo=optimize.brentq(f,0.05,om); om_hi=optimize.brentq(f,om,0.9)
models={"flat LCDM, matter density fitted":dl_lcdm(om),
        "static, Lerner: d = cz/H0, energy loss only":zb*np.sqrt(1+zb),
        "static, tired light: d = (c/H0) ln(1+z), energy loss only":np.log(1+zb)*np.sqrt(1+zb),
        "static, d = cz/H0, photons also stretched in time":zb*(1+zb),
        "static, d = (c/H0) ln(1+z), photons also stretched in time":np.log(1+zb)*(1+zb),
        "empty expanding universe (coasting)":zb*(1+zb/2),
        "expanding, matter only (no dark energy)":2*(1+zb)*(1-1/np.sqrt(1+zb))}
out={"om":[float(om),float(om_lo),float(om_hi)],"models":{},"z":zb.tolist(),"sig":np.sqrt(np.diag(C)).tolist()}
print(f"31 binned supernova distances, z = {zb[0]} to {zb[-1]}; errors {np.sqrt(np.diag(C)).min():.3f} to {np.sqrt(np.diag(C)).max():.3f} mag")
print(f"flat LCDM: matter density {om:.3f} (+{om_hi-om:.3f}, -{om-om_lo:.3f}); Betoule et al. quote 0.295 +- 0.034\n")
base=chi2(mu(models["flat LCDM, matter density fitted"]))
print(f"{'model':62s} chi2 (31 points)   delta chi2   model minus LCDM in mag at z = 0.1, 0.5, 1.0, 1.3 (matched at low z)")
ref=mu(models["flat LCDM, matter density fitted"])
for name,dl in models.items():
    c,a,r=chi2(mu(dl)); d=mu(dl)-ref; d=d-d[0]; pick=[float(np.interp(x,zb,d)) for x in (0.1,0.5,1.0,1.3)]
    nfree=2 if "LCDM" in name else 1; from scipy.stats import chi2 as ch; p=ch.sf(c,31-nfree)
    print(f"{name:62s} {c:7.1f}          {c-base[0]:+7.1f}     "+"  ".join(f"{x:+.3f}" for x in pick)+f"    p = {p:.1e}")
    out["models"][name]=dict(chi2=c,dchi2=c-base[0],offset=a,resid=r.tolist(),diff=pick,p=float(p))
# a general static law: d_L = z (1+z)^q
g=optimize.minimize_scalar(lambda q: chi2(mu(zb*(1+zb)**q))[0],bounds=(0,2),method="bounded"); fq=lambda q: chi2(mu(zb*(1+zb)**q))[0]-g.fun-1
qlo=optimize.brentq(fq,0,g.x); qhi=optimize.brentq(fq,g.x,2); print(f"\nd = cz/H0 with brightness falling as an extra (1+z)^(2q): best q = {g.x:.3f} (+{qhi-g.x:.3f}, -{g.x-qlo:.3f}), chi2 = {g.fun:.1f}; energy loss only is q = 0.5, energy loss and stretched time q = 1")
out["q"]=[float(g.x),float(qlo),float(qhi),float(g.fun)]
g2=optimize.minimize_scalar(lambda q: chi2(mu(np.log(1+zb)*(1+zb)**q))[0],bounds=(0,2),method="bounded"); fq2=lambda q: chi2(mu(np.log(1+zb)*(1+zb)**q))[0]-g2.fun-1
print(f"d = (c/H0) ln(1+z): best q = {g2.x:.3f} (+{optimize.brentq(fq2,g2.x,2)-g2.x:.3f}, -{g2.x-optimize.brentq(fq2,0,g2.x):.3f}), chi2 = {g2.fun:.1f}")
out["q_log"]=[float(g2.x),float(optimize.brentq(fq2,0,g2.x)),float(optimize.brentq(fq2,g2.x,2)),float(g2.fun)]
# residuals for the figure: data minus LCDM, and each model minus LCDM, each at its own best constant
out["fig"]={"data":[[float(z),float(r),float(s)] for z,r,s in zip(zb,base[2],np.sqrt(np.diag(C)))]}
zz=np.geomspace(0.01,1.3,60); ref_c=5*np.log10(np.array([(1+z)*integrate.quad(lambda x:1/E(x,om),0,z)[0] for z in zz]))
for key,fn in (("lerner",lambda z: z*np.sqrt(1+z)),("log",lambda z: np.log(1+z)*np.sqrt(1+z)),("lerner_td",lambda z: z*(1+z)),("log_td",lambda z: np.log(1+z)*(1+z)),("coast",lambda z: z*(1+z/2))):
    name=[n for n in models if np.allclose(models[n],fn(zb))][0]; a=out["models"][name]["offset"]
    out["fig"][key]=[[float(z),float(v)] for z,v in zip(zz,5*np.log10(fn(zz))+a-(ref_c+base[1]))]
json.dump(out,open("static_hubble.json","w"),indent=1)
