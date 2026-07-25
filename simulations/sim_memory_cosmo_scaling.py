# ==============================================================================
# FSG v6.1 -- Section "Progress Toward a Covariant Completion" (the AeST route)
# ------------------------------------------------------------------------------
# AeST reproduces the CMB because its scalar carries an energy density scaling
# dust-like (~a^-3) in matter domination. Does FSG's OWN memory field U=box^-1 R
# have that scaling, so the dark sector comes from FSG's memory rather than being
# added by hand? We integrate U on the FLRW background and examine rho_U~(dU/dt)^2.
#   Homogeneous box:  U'' + (3 + H'/H) U' = -R/H^2 ,  R/H^2 = 12 + 6 H'/H  (N=ln a).
# Consistency: U(a=1) must reproduce the U0 ~ -16 that fixes a0 in sim_friedmann_fsg.
# ==============================================================================
import numpy as np
from scipy.integrate import solve_ivp

Om=0.315; Or=9.0e-5; OL=1-Om-Or
def E2(N): a=np.exp(N); return Or*a**-4+Om*a**-3+OL
def dlnH(N):
    a=np.exp(N); return 0.5*(-4*Or*a**-4-3*Om*a**-3)/E2(N)
def rhs(N,Y):
    U,Up=Y; hh=dlnH(N)
    return [Up,-(3+hh)*Up-(12+6*hh)]

N=np.linspace(-15,0,6000)
s=solve_ivp(rhs,[N[0],N[-1]],[0,0],t_eval=N,rtol=1e-9,atol=1e-12)
U,Up=s.y; a=np.exp(N); rhoU=E2(N)*Up**2          # rho_U ~ (dU/dt)^2 = H^2 U'^2

def slope(av):
    i=np.argmin(np.abs(a-av)); return np.gradient(np.log(np.abs(rhoU)),N)[i]

print("="*70)
print(" FSG memory field U=box^-1 R on FLRW -- cosmological energy scaling")
print("="*70)
print(f"  U(a=1) = {U[-1]:.1f}   (fixes a0 in sim_friedmann_fsg: U0 ~ -16, consistent)")
print(f"  U'(matter era) = {Up[np.argmin(np.abs(a-0.05))]:.2f}  (analytic steady state -2.00)")
print("-"*70)
for lbl,av in [("radiation",1e-5),("matter (CMB epoch)",1e-3),("matter (late)",0.1),("today",1.0)]:
    s_=slope(av)
    tag=("~a^-3 (dark-matter-like)" if -3.4<s_<-2.3 else
         "~const (dark-energy-like)" if s_>-0.5 else "transition")
    print(f"  {lbl:20s} a={av:.0e}  dln(rho_U)/dln a = {s_:+.2f}  {tag}")
print("-"*70)
print("  => rho_U ~ a^-3 through the CMB epoch (dark-matter-like) and -> const today")
print("     (dark-energy-like): one memory field, both dark sectors, same U0 as a0.")
print("  NECESSARY not sufficient: the amplitude (Omega-equivalent, peak heights)")
print("  requires the perturbed non-local calculation. At quasi-static level the void")
print("  ISW sign comes out wrong (Cold Spot) -- see Section on honest boundaries.")
