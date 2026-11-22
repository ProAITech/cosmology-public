"""Workbench Entry 26: El Gordo. How much of the claimed conflict with LCDM sits in the collision speed?
The claim (Asencio, Banik & Kroupa 2021, 2023): two clusters this massive should not be colliding this fast this early; more than 5 sigma
unless the infall speed is below 2,300 km/s, while simulations that match the X-ray picture start the clusters at 2,500 km/s
at a separation of twice the sum of their radii. Three questions that need only Newtonian gravity:
  1. What speed does the pair's own gravity give at that separation?
  2. What does starting at 2,500 km/s instead imply?
  3. How different are the two cases at closest approach, where the X-ray picture is made?"""
import numpy as np, json
from scipy.optimize import brentq
G=4.30091e-9                      # Mpc (km/s)^2 / Msun
kmsMpc_to_perGyr=1.0227e-3        # 1 km/s/Mpc = 1.0227e-3 /Gyr
H0,Om=70.0,0.3
E=lambda z: np.sqrt(Om*(1+z)**3+1-Om)
def age(z,N=200000):              # Gyr
    a=np.linspace(1e-6,1/(1+z),N); return float(np.trapezoid(1/(a*np.sqrt(Om/a**3+1-Om)),a)/(H0*kmsMpc_to_perGyr))
rho_c=lambda z: 2.775e11*(H0/100)**2*E(z)**2            # Msun/Mpc^3
r200=lambda M,z: (3*M/(4*np.pi*200*rho_c(z)))**(1/3)    # Mpc
z=0.87; t_obs=age(z); out={"age_obs":t_obs}
print(f"age of the universe at z = {z}: {t_obs:.2f} Gyr")
cases={"Valdarnini (2024) simulation: 1.6e15 + 6.9e14":(1.6e15,6.9e14),"Kim et al. (2021) lensing total 2.13e15, ratio 1.5":(2.13e15*0.6,2.13e15*0.4),"Zhang et al. (2015) / Asencio et al. (2021): 3e15, ratio 3.6":(3e15*3.6/4.6,3e15/4.6)}
v_sim, v_peri_sim = 2500.0, 4000.0
def timing(M,t_peri):
    """two point masses that parted at t = 0 and meet again at t_peri (Gyr): largest separation they ever had (Mpc)"""
    T=t_peri*kmsMpc_to_perGyr; a=(G*M*(T/(2*np.pi))**2)**(1/3); return 2*a
def speed_on_orbit(M,rmax,r):       # km/s at separation r on the way back in
    return float(np.sqrt(max(2*G*M*(1/r-1/rmax),0)))
def timing_mass(d,v,t_peri):
    """mass for which two bodies that parted at t = 0 are d apart, closing at v, and meet at t_peri"""
    T=t_peri*kmsMpc_to_perGyr
    f=lambda lm: G*10**lm*(2/d-1/((G*10**lm*(T/(2*np.pi))**2)**(1/3)))-v*v
    return 10**brentq(f,14,18)
rows=[]
t_peri=t_obs-0.24
for name,(m1,m2) in cases.items():
    M=m1+m2; touch=r200(m1,z)+r200(m2,z); d=2*touch; vpar=np.sqrt(2*G*M/d)
    rmax=timing(M,t_peri); v_touch=speed_on_orbit(M,rmax,touch)
    vinf=np.sqrt(max(v_sim**2-vpar**2,0)); Mneed=v_sim**2*d/(2*G); Mt=timing_mass(d,v_sim,t_peri)
    print(f"\n{name}\n  radii {r200(m1,z):.2f} + {r200(m2,z):.2f} Mpc; simulations start at twice their sum, {d:.2f} Mpc")
    print(f"  timing argument (parted at the beginning, meet at {t_peri:.2f} Gyr): never more than {rmax:.2f} Mpc apart; closing at {v_touch:.0f} km/s when their edges touch ({touch:.2f} Mpc)")
    print(f"  at {d:.2f} Mpc the most their own gravity can give (falling from rest at infinity) is {vpar:.0f} km/s")
    print(f"  starting at {v_sim:.0f} km/s there: {v_sim/vpar:.2f} x that speed, {(v_sim/vpar)**2:.2f} x the energy; like arriving from far away at {vinf:.0f} km/s; needs {Mneed:.2e} Msun to be bound ({Mneed/M:.1f} x), or {Mt:.2e} Msun by the timing argument ({Mt/M:.1f} x)")
    rows.append({"name":name,"M":M,"d":d,"touch":touch,"rmax":rmax,"v_touch":v_touch,"v_parabolic":vpar,"v_inf":vinf,"M_needed":Mneed,"M_timing":Mt})
out["cases"]=rows
# 3. closest approach: energy conservation with the simulation's own two numbers (2,500 km/s at the start, about 4,000 km/s at closest approach)
gain=v_peri_sim**2-v_sim**2
print(f"\nspeed gained falling in, from the simulation: sqrt({v_peri_sim:.0f}^2 - {v_sim:.0f}^2) = {np.sqrt(gain):.0f} km/s (in quadrature)")
tab=[]
for v0 in (2500,2300,2000,1800,1500,1000,0):
    vp=np.sqrt(v0**2+gain); tab.append([v0,float(vp),float(vp/v_peri_sim)]); print(f"  start at {v0:5d} km/s -> closest approach at {vp:.0f} km/s ({100*(vp/v_peri_sim-1):+.0f}%); shock heating ~ v^2: {100*((vp/v_peri_sim)**2-1):+.0f}%")
out["peri"]=tab
for lab,v0 in (("falling from rest at infinity",rows[1]["v_parabolic"]),("Asencio et al.'s 5 sigma line",2300.0),("Ng et al. (2015)-based estimate, 3.5 sigma",1200.0),("Kim et al. nominal, 2.7 sigma",450.0)):
    vp=np.sqrt(v0**2+gain); print(f"  {lab}: {v0:.0f} km/s at the start -> {vp:.0f} km/s at closest approach, {100*(1-vp/v_peri_sim):.0f}% below the simulation's 4,000")
    out.setdefault("named",[]).append([lab,float(v0),float(vp)])
# check: Ng et al. (2015) quote 2,400 (+400/-200) km/s at closest approach; their method lets the clusters fall from rest at touching separation
touch=rows[0]["touch"]; d=rows[0]["d"]; M=rows[0]["M"]; g2=gain-2*G*M*(1/touch-1/d)
print(f"\nfalling from rest at touching separation ({touch:.2f} Mpc) with the same potential: {np.sqrt(g2):.0f} km/s at closest approach (Ng et al.: 2,400 +400/-200)")
out["peri_from_touch"]=float(np.sqrt(g2))
json.dump(out,open("infall.json","w"),indent=1)
