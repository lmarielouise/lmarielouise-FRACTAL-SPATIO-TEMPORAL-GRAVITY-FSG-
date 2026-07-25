# ==============================================================================
# FSG v6.1 -- Section "Progress Toward a Covariant Completion" (the AeST route)
# ------------------------------------------------------------------------------
# Reconstruct the AQUAL / AeST free function F(X) implied by FSG's d_eff, and check
# it lies in the AeST family (Skordis & Zlosnik 2021): correct deep-MOND and
# Newtonian limits.
#   AQUAL:  L = -(a0^2/8piG) F(X),  X=|grad Phi|^2/a0^2,  mu(x)=F'(x^2).
#   FSG (QUMOND):  g = nu(g_N/a0) g_N,  nu(y)=1/(1-exp(-sqrt y)),  from d_eff.
#   Convert:  x = nu(y) y = g/a0 ,  mu(x)=1/nu(y)=g_N/g ;  F(X)=int_0^X mu(sqrt X')dX'.
# ==============================================================================
import numpy as np
from scipy.integrate import cumulative_trapezoid

y=np.geomspace(1e-8,1e8,200000)          # y = g_N/a0
nu=1.0/(1.0-np.exp(-np.sqrt(y)))
x=nu*y                                    # = g/a0
mu=1.0/nu                                 # = mu(x)  (AQUAL interpolation)
X=x**2
o=np.argsort(X); Xs=X[o]; mus=mu[o]
F=cumulative_trapezoid(mus,Xs,initial=0.0)

print("="*70)
print(" FSG d_eff -> AeST free function F(X)   (X = g^2/a0^2)")
print("="*70)
print(f"{'X':>10}{'mu':>12}{'F(X)':>14}{'(2/3)X^1.5':>14}{'X':>10}")
for Xq in [1e-4,1e-2,1.0,1e2]:
    print(f"{Xq:>10.0e}{np.interp(Xq,Xs,mus):>12.4f}{np.interp(Xq,Xs,F):>14.4e}"
          f"{(2/3)*Xq**1.5:>14.4e}{Xq:>10.1f}")

iL=np.argmin(np.abs(Xs-1e-5)); iH=np.argmin(np.abs(Xs-1e5))
print("-"*70)
print(f" deep-MOND  X->0 :  F/((2/3)X^3/2) = {F[iL]/((2/3)*Xs[iL]**1.5):.4f}  (AeST needs 1)")
print(f" Newtonian  X->oo:  F'(X)=mu       = {mus[iH]:.4f}  (AeST needs 1)")
print(" => F_FSG interpolates between (2/3)X^3/2 (AeST MOND regime) and X (canonical);")
print("    it lies in the AeST free-function family, with the exponential shape as the")
print("    geometric fingerprint of d_eff. NECESSARY, not sufficient: CMB amplitudes")
print("    still require the perturbed non-local (Boltzmann) calculation.")
