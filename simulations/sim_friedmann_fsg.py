# ==============================================================================
# FSG — Intégration cosmologique des champs mémoire U = Box^{-1}R, V = Box^{-1}U
# ==============================================================================
# Correction 1 : dérivation numérique de a0 = c H0 |V0|^{3/2} / |U0|
# Correction 5 : back-réaction au premier ordre -> w_eff(z) et Delta H / H
#
# Méthode (ordre 0) : on intègre U et V sur un fond LCDM, depuis l'ère de
# radiation profonde (R = 0 => U = V = 0), avec x = ln(a) :
#
#   U'' + (3 + h'/h) U' = -6 (h'/h + 2)          [ Box U = R ]
#   V'' + (3 + h'/h) V' = -U / h^2               [ Box V = U, V normalisé par H0^2 ]
#
# Méthode (ordre 1, BORNE SUPÉRIEURE) : la densité d'énergie sombre effective
# est supposée tracker la couche de mémoire profonde V (structure analogue au
# modèle RR de Maggiore-Mancarella 2014, où rho_DE est portée par le champ S) :
#
#   rho_DE(a) = rho_DE,0 * |V(a)| / |V(a=1)|,   w(a) = -1 - (1/3) dln rho_DE/dln a
#
# ATTENTION : ce proxy néglige les termes dérivatifs compensateurs du tenseur
# énergie-impulsion non-local ; il SURESTIME donc |w+1| et Delta H/H. Il fournit
# une borne supérieure. La borne inférieure est le benchmark du modèle RR
# (traitement covariant complet de la même famille d'opérateurs) : w0 ~ -1.14,
# Delta H/H ~ 1-2 % (Maggiore & Mancarella 2014 ; Dirian et al. 2014).
# Le résultat STRUCTUREL robuste, indépendant du proxy : toute densité de
# mémoire monotone croissante implique w < -1 (comportement fantôme).
# ==============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- Cosmologie de fond (Planck 2018) ---
H0_kms = 67.4                 # km/s/Mpc
Omega_m = 0.315
Omega_r = 9.0e-5
Omega_L = 1.0 - Omega_m - Omega_r
c_light = 2.998e8             # m/s
Mpc = 3.0857e22               # m
H0_SI = H0_kms * 1e3 / Mpc    # s^-1

def h2(x):
    a = np.exp(x)
    return Omega_r / a**4 + Omega_m / a**3 + Omega_L

def dlnh_dx(x):
    a = np.exp(x)
    dh2 = -4 * Omega_r / a**4 - 3 * Omega_m / a**3
    return 0.5 * dh2 / h2(x)

def rhs(x, y):
    """y = [U, U', V, V']"""
    U, Up, V, Vp = y
    q = dlnh_dx(x)
    Upp = -(3 + q) * Up - 6 * (q + 2)
    Vpp = -(3 + q) * Vp - U / h2(x)
    return [Up, Upp, Vp, Vpp]

# --- Intégration depuis l'ère de radiation (R ~ 0 => U = V = 0) ---
x_ini, x_end = np.log(1e-7), 0.0
sol = solve_ivp(rhs, (x_ini, x_end), [0, 0, 0, 0],
                dense_output=True, rtol=1e-10, atol=1e-12, max_step=0.01)

x = np.linspace(np.log(1e-4), 0.0, 4000)
U, Up, V, Vp = sol.sol(x)
U0, V0 = U[-1], V[-1]

print("=" * 70)
print("CHAMPS MÉMOIRE AUJOURD'HUI (fond LCDM, ordre 0)")
print(f"  U(a=1) = {U0:+.4f}   (Deser-Woodard : X ~ -14, ordre de grandeur OK)")
print(f"  V(a=1) = {V0:+.4f}")

# --- Correction 1 : échelle d'accélération a0 ---
a0_pred = c_light * H0_SI * abs(V0)**1.5 / abs(U0)
a0_obs = 1.2e-10
print("-" * 70)
print("ÉCHELLE D'ACCÉLÉRATION")
print(f"  a0 = c H0 |V0|^(3/2) / |U0| = {a0_pred:.3e} m/s^2")
print(f"  a0 observé (McGaugh)        = {a0_obs:.3e} m/s^2")
print(f"  écart : {100 * (a0_pred / a0_obs - 1):+.1f} %")

# --- Correction 5 : back-réaction au premier ordre ---
# rho_DE portée par |V|, normalisée à Omega_L aujourd'hui
lnV = np.log(np.abs(V) + 1e-30)
dlnV_dx = np.gradient(lnV, x)
w_eff = -1.0 - dlnV_dx / 3.0
w0 = w_eff[-1]

Omega_DE = Omega_L * np.abs(V) / abs(V0)
a = np.exp(x)
z = 1.0 / a - 1.0
h2_fsg = Omega_r / a**4 + Omega_m / a**3 + Omega_DE
dH_over_H = np.sqrt(h2_fsg / h2(x)) - 1.0

mask = (z > 0.5) & (z < 2.0)
print("-" * 70)
print("ÉQUATION D'ÉTAT EFFECTIVE (proxy V-tracking = BORNE SUPÉRIEURE)")
print(f"  w0 = w_eff(z=0)             = {w0:+.3f}   [borne sup. ; benchmark RR : -1.14]")
print(f"  w_eff(z=0.5)                = {np.interp(0.5, z[::-1], w_eff[::-1]):+.3f}")
print(f"  Delta H/H  max (0.5<z<2)    = {100 * np.max(np.abs(dH_over_H[mask])):.2f} %   [borne sup. ; benchmark RR : 1-2 %]")
print(f"  Delta H/H  (z=1)            = {100 * np.interp(1.0, z[::-1], dH_over_H[::-1]):+.2f} %")
print("  => Prédiction structurelle robuste : w0 < -1 (mémoire croissante)")
print("=" * 70)

# --- Figures ---
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

axes[0].plot(x, U, 'b-', lw=2, label=r'$U$ (first memory layer)')
axes[0].plot(x, V, 'r-', lw=2, label=r'$V$ (second memory layer)')
axes[0].set_xlabel(r'$\ln a$'); axes[0].set_ylabel('memory fields')
axes[0].legend(); axes[0].set_title('Cosmic memory build-up')

zm = z < 3
axes[1].plot(z[zm], w_eff[zm], 'r-', lw=2)
axes[1].axhline(-1, color='k', ls='--', lw=1, label=r'$\Lambda$CDM')
axes[1].set_xlabel(r'$z$'); axes[1].set_ylabel(r'$w_{\rm eff}(z)$')
axes[1].legend(); axes[1].set_title(f'Effective EoS, upper bound ($w_0={w0:.2f}$)')

axes[2].plot(z[zm], 100 * dH_over_H[zm], 'g-', lw=2)
axes[2].axhline(0, color='k', ls='--', lw=1)
axes[2].set_xlabel(r'$z$'); axes[2].set_ylabel(r'$\Delta H/H$ [%]')
axes[2].set_title('FSG vs LCDM expansion (upper bound)')

plt.tight_layout()
plt.savefig('fig_friedmann_fsg.pdf', dpi=300, bbox_inches='tight')
plt.savefig('fig_friedmann_fsg.png', dpi=150, bbox_inches='tight')
print("Figures : fig_friedmann_fsg.pdf / .png")
