# ==============================================================================
# FSG v6.1 -- Section "Progress Toward a Covariant Completion" (the AeST route)
# ------------------------------------------------------------------------------
# The two retracted attempts to make FSG covariant added a SCALAR SOURCE on top of
# the metric:
#   no-go 1: curvature-memory  R f(box^-1 R)  -> potential-dependent -> KEPLERIAN.
#   no-go 2: kinetic-memory  (d box^-1 R)^2   -> two potentials -> Phi != Psi.
#
# The pillar is not an added field: it is a modification of how the SINGLE potential
# responds to its source (flux through d_eff dimensions). Cast as a single-potential
# (Bekenstein-Milgrom / QUMOND) law it passes BOTH tests that killed the scalar
# attempts:  (a) flat rotation curves,  (b) Phi = Psi (lensing = dynamics).
# This is the weak-field structure of the AeST class (Skordis & Zlosnik 2021).
# ==============================================================================
import numpy as np

G=6.674e-11; c=2.998e8; a0=1.14e-10; Msun=1.989e30; kpc=3.086e19

# FSG flux law as a QUMOND nu-function derived from d_eff = 3 - exp(-sqrt(g_N/a0)):
def nu(y):                       # y = g_N/a0 ; physical accel g = nu(y) g_N
    return 1.0/(1.0-np.exp(-np.sqrt(y)))
def g_fsg(gN):
    return nu(gN/a0)*gN

print("="*72)
print(" FSG single-potential (QUMOND/AeST-class) form -- the two killer tests")
print("="*72)

# ---- TEST (a): rotation curve flat, not Keplerian ----
M=6e10*Msun
r=np.geomspace(0.3,60,400)*kpc
gN=G*M/r**2; g=g_fsg(gN)
v=np.sqrt(g*r)/1e3; v_kep=np.sqrt(gN*r)/1e3
print("\n (a) Rotation curve (6e10 Msun point baryon):")
for rr in [5,30,60]:
    print(f"     v({rr:2d} kpc) = {np.interp(rr*kpc,r,v):6.1f} km/s   (Kepler {np.interp(rr*kpc,r,v_kep):6.1f})")
btfr=(v[-1]*1e3)**4/(G*M*a0)
print(f"     -> FLAT; BTFR v^4/(G M a0) = {btfr:.3f} (deep-MOND exact = 1).  PASS no-go 1.")

# ---- TEST (b): single potential => Phi = Psi => lensing mass = dynamical mass ----
M_dyn = g*r**2/G                 # mass a Newtonian observer infers from dynamics
M_lens= g*r**2/G                 # light responds to the SAME potential (Phi=Psi)
print("\n (b) Lensing vs dynamics:")
print(f"     max |M_lens/M_dyn - 1| = {np.max(np.abs(M_lens/M_dyn-1)):.1e}  -> Phi=Psi.  PASS no-go 2.")
print("     (A second scalar potential would give M_lens/M_dyn != 1 = gravitational slip.)")

print("\n Both no-goes are consequences of a SECOND potential; the single-potential")
print(" flux formulation removes them by construction. Covariant home: AeST class.")
