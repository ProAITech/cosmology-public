"""Workbench Entry 18: the Milky Way's halo circuit. What drives ~1e17 A, and is it fed from outside?
Inputs are measured quantities from the literature plus our own medians from the NVSS rotation-measure
catalogue (Taylor, Stil & Sunstrum 2009, VizieR J/ApJ/702/1230; 37,543 sources), computed in a browser tab."""
import json, numpy as np
MU0=4e-7*np.pi; KPC=3.0857e19; MYR=3.156e13; GYR=1e3*MYR; KMSKPC=1e3/KPC   # s^-1 per (km/s/kpc)
uG=1e-10
out={}

# ---------- 1. NVSS rotation measures: medians at 20 < |b| < 50 deg in 30-deg longitude bins [l0, north, err, n, south, err, n]
lon=[[0,10.1,0.9,745,-10.7,1.3,693],[30,40.0,1.3,863,-8.3,1.4,757],[60,20.5,0.8,952,-27.2,1.4,919],[90,2.5,1.2,939,-50.8,1.4,979],
     [120,-15.6,0.8,894,-49.8,1.4,1002],[150,1.3,0.7,995,0.4,1.1,912],[180,15.4,0.7,949,12.4,1.2,820],[210,14.1,1.0,992,16.5,1.1,810],
     [240,-6.1,1.1,847,24.7,3.0,140],[270,-17.4,1.2,799,None,None,0],[300,-21.4,1.1,849,None,None,0],[330,-7.2,0.9,779,8.9,6.0,12]]
quad={"Q1 north":[18.3,0.5,3186],"Q1 south":[-11.9,0.7,2924],"Q4 north":[-9.8,0.5,3050],"Q4 south":[11.5,4.3,14]}
pole={"north |b|>70":[2.6,0.4,1753],"south |b|>70":[6.1,0.6,1235],"north |b|>80":[3.2,0.8,452],"south |b|>80":[6.2,1.1,337]}
L=np.array([r[0]+15 for r in lon if r[4] is not None]); N=np.array([r[1] for r in lon if r[4] is not None]); S=np.array([r[4] for r in lon if r[4] is not None])
eN=np.array([r[2] for r in lon if r[4] is not None]); eS=np.array([r[5] for r in lon if r[4] is not None])
s=np.sin(np.radians(L)); odd=(N-S)/2; even=(N+S)/2; e=np.sqrt(eN**2+eS**2)/2
A_odd=(odd*s).sum()/(s*s).sum(); A_even=(even*s).sum()/(s*s).sum(); eA=np.sqrt((e*e*s*s).sum())/(s*s).sum()
print(f"NVSS: part that reverses across the plane, fitted as A sin(l): A = {A_odd:+.1f} +- {eA:.1f} rad/m2 (statistical); part that does not reverse: {A_even:+.1f}")
DM_perp=25.0; b=35.0
DM=DM_perp/np.sin(np.radians(b)); Bphi_nvss=A_odd/(0.812*DM*np.cos(np.radians(b)))
print(f"  -> toroid field {Bphi_nvss:.2f} uG for an electron column of {DM_perp:.0f} pc/cm3 toward the pole (path at |b|=35: {DM:.0f}); Xu & Han: 0.73 +- 0.12")
Bz_S=pole["south |b|>70"][0]/(0.812*DM_perp); eBz_S=pole["south |b|>70"][1]/(0.812*DM_perp); Bz_N=-pole["north |b|>70"][0]/(0.812*DM_perp)
print(f"  vertical field from the polar caps: south {Bz_S:+.2f} +- {eBz_S:.2f} uG (pointing north), north {Bz_N:+.2f} uG (pointing south); Taylor+2009: +0.30, -0.14; Mao+2010 north: 0.0")
out["nvss"]={"lon":lon,"quad":quad,"pole":pole,"A_odd":A_odd,"A_odd_err":eA,"A_even":A_even,"Bphi":Bphi_nvss,"Bz_south":Bz_S,"Bz_north":Bz_N}

# ---------- 2. The circuit, from Ampere's law
R=8*KPC; Bphi=0.73*uG; Bz=0.30*uG; V=220e3
I=2*np.pi*R*Bphi/MU0; K=2*Bphi/MU0
print(f"\nCurrent through each half of the halo inside 8 kpc: {I:.2e} A; opposite senses north and south -> radial current in the disk {2*I:.2e} A ({K:.2e} A per metre of circumference)")
R1,R2=2*KPC,15*KPC
emf=V*Bz*(R2-R1); P=emf*2*I
print(f"EMF of the rotating disk across its vertical field (v x B from 2 to 15 kpc): {emf:.2e} V; power into the circuit {P:.2e} W")
Lmw=1e37; Psn=1e44*2/(100*3.156e7)
print(f"  = {P/Lmw:.1e} of the Galaxy's luminosity, {100*P/Psn:.1f}% of its supernova power ({Psn:.1e} W)")
out["circuit"]={"I_half":I,"I_disk":2*I,"K":K,"emf":emf,"P":P,"P_over_L":P/Lmw,"P_over_SN":P/Psn}

# ---------- 3. Winding: dB_phi/dt = B_R * R dOmega/dR + B_z * dv_phi/dz
lag=15.0; elag=4.0; rad=220/8.0          # km/s/kpc: halo lag (Marasco & Fraternali 2011), V/R for a flat rotation curve
ratio=Bphi/Bz; eratio=ratio*np.hypot(0.12/0.73,0.03/0.30)
t_lag=ratio/(lag*KMSKPC)/MYR; t_lag_lo=ratio/((lag+elag)*KMSKPC)/MYR; t_lag_hi=ratio/((lag-elag)*KMSKPC)/MYR
t_both=ratio/((lag+rad)*KMSKPC)/MYR
Prot=2*np.pi*8*KPC/V/MYR
print(f"\nToroid / vertical field = {ratio:.2f} +- {eratio:.2f}")
print(f"Winding time to make it: {t_lag:.0f} Myr (range {t_lag_lo:.0f}-{t_lag_hi:.0f}) from the halo lag alone; {t_both:.0f} Myr if radial shear acts on an equal radial field; Unger & Farrar fit: 54.7 +- 1.1 Myr")
print(f"  one Galactic rotation at the Sun = {Prot:.0f} Myr -> {t_both/Prot:.2f} to {t_lag/Prot:.2f} of a turn")
tau=1.6*KPC/20e3/MYR
print(f"Turnover time of halo gas (1.6 kpc at 20 km/s): {tau:.0f} Myr -> predicted toroid/vertical = {lag*KMSKPC*tau*MYR:.1f} (lag only) to {(lag+rad)*KMSKPC*tau*MYR:.1f} (both shears); measured {ratio:.1f}; Unger & Farrar {3.2/0.98:.1f}")
B10=Bz/uG*np.array([lag,lag+rad])*KMSKPC*10*GYR
print(f"If the winding had run for 10 Gyr with nothing lost: {B10[0]:.0f}-{B10[1]:.0f} uG, against 0.73 measured")
out["winding"]={"ratio":ratio,"ratio_err":eratio,"t_lag":t_lag,"t_lag_range":[t_lag_lo,t_lag_hi],"t_both":t_both,"P_rot":Prot,"tau":tau,"pred":[lag*KMSKPC*tau*MYR,(lag+rad)*KMSKPC*tau*MYR],"B10":list(B10)}

# ---------- 4. Energy and torque
vol=2*np.pi*(R2**2-R1**2)*3*KPC; U=Bphi**2/(2*MU0)*vol
Pshear=Bphi*Bz*(lag+rad)*KMSKPC/MU0*vol
print(f"\nEnergy stored in the toroids: {U:.1e} J; work done by shear {Pshear:.1e} W -> refill time {U/Pshear/MYR:.0f} Myr")
stress=Bphi*Bz/MU0; T=2*stress*(2*np.pi/3)*(R2**3-R1**3)
Lgas=1e10*1.989e30*R*V; Lhalo=3.2e8*1.989e30*R*(V-20e3)
print(f"Magnetic torque between disk and halo: {T:.1e} N m; time to brake the gas disk {Lgas/T/GYR:.0f} Gyr; to spin up the HI halo {Lhalo/T/GYR:.1f} Gyr = {Lhalo/T/MYR/tau:.0f} turnover times")
out["energy"]={"U":U,"P_shear":Pshear,"refill_Myr":U/Pshear/MYR,"torque":T,"t_brake_disk_Gyr":Lgas/T/GYR,"t_spin_halo_Gyr":Lhalo/T/GYR}

# ---------- 5. A current threading the Galaxy pole to pole would give the same sense north and south
BN,BS,eB=3.3,-3.1,0.3
evenpart=(BN+BS)/2; oddpart=(BN-BS)/2; ee=eB/np.sqrt(2)
frac=(evenpart+2*ee)/oddpart
print(f"\nUnger & Farrar halo toroids: north {BN:+.1f}, south {BS:+.1f} uG -> same-sense part {evenpart:+.2f} +- {ee:.2f}, i.e. below {100*frac:.0f}% of the reversing part (95%)")
print(f"  -> any current threading the Galaxy along its axis is below about {frac*I:.1e} A; Peratt's galactic currents: 1e17-1e19 A")
out["through"]={"even":evenpart,"err":ee,"frac":frac,"I_max":frac*I}

# ---------- 6. Where does the vertical flux come from?
flux=Bz*np.pi*R2**2; comp=(1e-3/2.5e-7)**(2/3)
print(f"\nVertical flux through the disk inside 15 kpc: {flux:.1e} Wb. Intergalactic field needed if only compressed (x{comp:.0f}): {Bz/uG*1e3/comp:.1f} nG; CMB clumping limit 0.045 nG -> a further factor {Bz/uG*1e3/comp/0.045:.0f} ({np.log(Bz/uG*1e3/comp/0.045):.1f} e-folds) from a dynamo or from flux carried in by accretion")
out["seed"]={"flux_Wb":flux,"compression":comp,"B_igm_nG":Bz/uG*1e3/comp,"extra_factor":Bz/uG*1e3/comp/0.045}

# ---------- 7. Sign test for the source of the vertical field (Part 5)
# Outward axial current <=> vertical field against the spin. Dynamo or primordial seed: 50% of galaxies. Cosmic Battery: all.
capS=[6.3,3.3,2.4,9.4,6.2,6.9]; capN=[4.8,2.6,-0.6,6.0,-2.9,3.7]   # NVSS polar-cap medians in six 60-deg sectors of longitude, rad/m2
print(f"\nPolar caps by sector: south {min(capS):+.1f} to {max(capS):+.1f} (all positive), north {min(capN):+.1f} to {max(capN):+.1f}")
print("Milky Way: spin points south (clockwise seen from north), vertical field points north -> against the spin -> outward current (as measured)")
print("CHANG-ES stack (Myserlis & Contopoulos 2021): correlation 0.70, p ~ 5e-3 (~3 sigma), 1,931 RMs from 24 galaxies; jets (Christodoulou+ 2016): ~3 sigma")
out["sign_test"]={"capS":capS,"capN":capN}
json.dump(out,open("halo_results.json","w"),indent=1)
