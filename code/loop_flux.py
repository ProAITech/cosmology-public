"""Workbench Entry 54, check: how much net magnetic flux does a tangled field put through a large loop?
A field made by causal processes in the early universe is tangled on a small scale lam.  Because field lines close, the flux through
a disc of radius R >> lam does not grow as the area: it is set by the rim.  Expectation: flux_rms ~ B lam^2 (2 pi R / lam)^(1/2).
A field with equal power on all scales (as inflation could make) is coherent across the disc: flux ~ B R^2.
Here: random divergence-free fields on a grid, flux through discs of growing radius."""
import os, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); rng=np.random.default_rng(54); N=160
k1=np.fft.fftfreq(N)*2*np.pi; kx,ky,kz=np.meshgrid(k1,k1,k1,indexing='ij'); k=np.sqrt(kx**2+ky**2+kz**2); k[0,0,0]=1.0
def field(PB):
    """Gaussian random B with power spectrum PB(k): draw a vector potential with spectrum PB/k^2 and take its curl."""
    amp=np.sqrt(PB(k))/k; amp[0,0,0]=0; A=[np.fft.fftn(rng.standard_normal((N,N,N)))*amp for _ in range(3)]
    Bz=np.fft.ifftn(1j*(kx*A[1]-ky*A[0])).real; Bx=np.fft.ifftn(1j*(ky*A[2]-kz*A[1])).real; By=np.fft.ifftn(1j*(kz*A[0]-kx*A[2])).real
    return Bx,By,Bz
def disc_flux(Bz,R):
    x=np.fft.fftfreq(N)*N; X,Y=np.meshgrid(x,x,indexing='ij'); ker=(X**2+Y**2<=R*R).astype(float); fk=np.fft.fft2(ker)
    F=np.fft.ifft2(np.fft.fft2(Bz,axes=(0,1))*fk[:,:,None],axes=(0,1)).real          # flux through a disc centred on every point of every plane
    return F.std()
out={}; radii=[4,6,8,12,16,24,32]
kc=2*np.pi/8.0                                              # energy spectrum k^4 exp(-2k^2/kc^2) peaks at k = kc: lam = 8 cells
res=[]
for trial in range(3):
    Bx,By,Bz=field(lambda q: q**2*np.exp(-2*q**2/kc**2)); Brms=np.sqrt((Bx**2+By**2+Bz**2).mean()); res.append([disc_flux(Bz,R)/Brms for R in radii])
res=np.mean(res,axis=0); lam=8.0; sl=np.polyfit(np.log(radii[2:]),np.log(res[2:]),1)[0]; c=np.mean(res[2:]/(lam**1.5*np.sqrt(np.array(radii[2:]))))
print(f"tangled field (spectrum rising as k^2, coherence lam = {lam:.0f} cells):")
for R,f in zip(radii,res): print(f"   R = {R:2d} ({R/lam:3.1f} lam): flux_rms / B_rms = {f:7.1f} cells^2   = {f/(np.pi*R*R):.4f} of a uniform field's flux")
print(f"   slope of log flux on log R beyond one coherence length: {sl:.2f} (rim scaling predicts 0.5; area scaling 2)")
print(f"   flux_rms = {c:.2f} x B_rms lam^1.5 R^0.5")
out['tangled']=dict(radii=radii,flux=[float(x) for x in res],slope=float(sl),c=float(c),lam=lam)
kmin=2*np.pi/N; res=[]
for trial in range(3):
    Bx,By,Bz=field(lambda q: np.where(q>=kmin,q**-3.0,0.0)*np.exp(-(q/2.0)**2)); Brms=np.sqrt((Bx**2+By**2+Bz**2).mean()); res.append([disc_flux(Bz,R)/Brms for R in radii])
res=np.mean(res,axis=0); sl2=np.polyfit(np.log(radii),np.log(res),1)[0]
print(f"\nfield with equal power per octave (spectrum k^-3):")
for R,f in zip(radii,res): print(f"   R = {R:2d}: flux_rms / B_rms = {f:8.1f} cells^2 = {f/(np.pi*R*R):.3f} of a uniform field of the full strength")
print(f"   slope: {sl2:.2f}")
out['flat']=dict(radii=radii,flux=[float(x) for x in res],slope=float(sl2),frac=[float(f/(np.pi*R*R)) for R,f in zip(radii,res)])
json.dump(out,open(os.path.join(HERE,'loop_flux.json'),'w'),indent=1)
