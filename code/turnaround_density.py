"""Workbench Entry 35, Part 4: the Alfven-Klein cloud at the turnaround density its 1962 paper states (about 1e-2 atoms per cm3).
A cloud of mass M at that density has radius R; its Schwarzschild radius is 2GM/c^2. Inside it, a contracting cloud cannot turn round."""
import math, json, os
G=6.674e-11; c=2.998e8; mp=1.6735e-27; ly=9.461e15; Msun=1.989e30
rho=1e4*mp                                             # 1e-2 atoms per cm3
H0=68*1000/3.0857e22; rho_b=0.0224/0.68**2*3*H0**2/(8*math.pi*G)
out={"rho_b_per_cm3":rho_b/mp/1e6,"rows":[]}
cases=(("ordinary matter inside the Hubble distance",4/3*math.pi*(c/H0)**3*rho_b),("ordinary matter in view, to 46.5 billion light years",4/3*math.pi*(46.5e9*ly)**3*rho_b),("the same with an equal mass of antimatter",2*4/3*math.pi*(46.5e9*ly)**3*rho_b))
for nm,M in cases:
    R=(3*M/(4*math.pi*rho))**(1/3); Rs=2*G*M/c**2
    print(f"{nm:52s} M = {M/Msun:.2e} Msun; radius {R/ly/1e9:.2f} Gly; Schwarzschild radius {Rs/ly/1e9:.2f} Gly; ratio {R/Rs:.3f}; GM/Rc^2 = {G*M/(R*c*c):.1f}")
    out["rows"].append([nm,M/Msun,R/ly/1e9,Rs/ly/1e9,R/Rs,G*M/(R*c*c)])
Mmax=c**3/(2*G)**1.5*(3/(4*math.pi*rho))**0.5; print(f"largest mass outside its Schwarzschild radius at this density: {Mmax/Msun:.2e} Msun")
M=cases[1][1]; nmax=3*c**6/(32*math.pi*G**3*M**2)/mp/1e6; print(f"highest density at which the matter in view is outside its Schwarzschild radius: {nmax:.1e} per cm3 ({1e-2/nmax:.0f} times below the paper's); mean ordinary-matter density today {rho_b/mp/1e6:.1e} per cm3")
out["Mmax"]=Mmax/Msun; out["nmax"]=nmax
json.dump(out,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"turnaround_density.json"),"w"),indent=1)
