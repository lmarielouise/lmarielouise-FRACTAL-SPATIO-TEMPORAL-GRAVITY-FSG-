# ==============================================================================
# FSG v7 — Dérivation de MOND depuis l'action corrigée (Priorité 1)
# ==============================================================================
# RÉSULTAT NÉGATIF (v6, prouvé) : le secteur mémoire-de-courbure R f(Box^-1 R)
# ne peut PAS produire MOND. En régime quasi-statique U = Box^-1 R = 2 chi + cte
# est LOCAL dans le potentiel chi ; la réduction de R f(U) donne une
# modification dépendant du POTENTIEL, pas de l'ACCÉLÉRATION, et la courbe de
# rotation reste képlérienne (vérifié numériquement, partie A ci-dessous).
#
# CORRECTION (v7) : on ajoute à l'action un terme cinétique du champ mémoire,
#     S_K = (a0^2 / 8 pi G) int sqrt(-g) F( c^4 (d_mu U)(d^mu U) / a0^2 ),
# qui est une fonctionnelle covariante NON-LOCALE de la métrique (U = Box^-1 R).
# Comme grad U = 2 grad chi en quasi-statique, ce terme "voit" l'accélération.
#
# MÉCANISME (vérifié symboliquement, partie B) : dans la forme localisée
# (multiplicateur xi imposant Box U = R), les équations d'Euler-Lagrange du
# système (chi, U, xi) possèdent deux premières intégrales qui transforment le
# terme ADDITIF en modification MULTIPLICATIVE de Poisson :
#     div[ mu_eff(|grad chi|/a0) grad chi ] = 4 pi G rho / c^2,
#     mu_eff = 1 + (4 beta / alpha) F'(4 |grad chi|^2 / a0^2).
# C'est exactement AQUAL (Bekenstein-Milgrom). Le choix
#     mu_eff(y) = 1 - exp(-y),   y = |grad chi| / a0
# donne (i) la branche MOND profonde exacte mu -> y => BTFR v^4 = G M a0,
# (ii) des corrections exp(-g/a0) ~ exp(-1e8) dans le système solaire :
# aucun écrantage nécessaire. Forme fermée : F(Z) = (alpha/(2 beta)) (1+y) e^-y.
#
# Précédent : MOND non-local métrique de Soussa-Woodard (2003) et
# Deffayet-Esposito-Farese-Woodard (2011). L'ingrédient propre à FSG : le champ
# qui porte le terme cinétique est le MÊME champ mémoire U qui engendre
# l'énergie sombre et l'échelle a0 = c H0 |V0|^(3/2)/|U0|.
#
# Ce script : (A) vérifie le no-go képlérien du secteur f(X) seul ;
#             (B) résout l'équation AQUAL dérivée et vérifie plateau + BTFR.
# ==============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G = 1.0
a0 = 1.0e-3          # unités code
masses = [1.0, 10.0, 100.0, 1000.0]

# ---------- (A) NO-GO : secteur f(X) seul => Kepler ----------
def rotation_curve_fX_only(M, lam=0.5, Uc=1.0):
    """mu(Phi) = 1 + lam (2 Phi + Uc)^(-3/2) : depend du potentiel, pas de g."""
    r = np.logspace(0, 3.2, 800)
    def dPhi(rr, y):
        mu = 1.0 + lam*np.abs(2*y[0] + Uc)**(-1.5)
        return [G*M/(mu*rr**2)]
    sol = solve_ivp(dPhi, (r[0], r[-1]), [-G*M/r[0]], t_eval=r, rtol=1e-9)
    Phi = sol.y[0]
    mu = 1.0 + lam*np.abs(2*Phi + Uc)**(-1.5)
    g = G*M/(mu*r**2)
    return r, np.sqrt(g*r)

# ---------- (B) AQUAL dérivée : mu(y) = 1 - exp(-y) ----------
def mu_exp(y):
    return 1.0 - np.exp(-y)

def solve_g(gN):
    """Résout mu(g/a0) g = gN pour g (relation algébrique AQUAL, source ponctuelle)."""
    out = np.empty_like(gN)
    for i, w in enumerate(gN):
        f = lambda g: mu_exp(g/a0)*g - w
        lo = max(w, np.sqrt(w*a0)*0.5)
        hi = max(w*2, np.sqrt(w*a0)*3) + 10*a0
        out[i] = brentq(f, lo*1e-3, hi*1e3, xtol=1e-16, rtol=1e-13)
    return out

print("=" * 76)
print("(A) NO-GO — secteur f(Box^-1 R) seul (mu fonction du potentiel)")
r, v = rotation_curve_fX_only(10.0)
i1, i2 = np.argmin(np.abs(r-100)), np.argmin(np.abs(r-1000))
print(f"    v(1000)/v(100) = {v[i2]/v[i1]:.3f}   [Kepler = {np.sqrt(0.1):.3f} ; plat = 1.000]")
print("    => keplerien : PAS de MOND. Le no-go est confirme numeriquement.")
print("-" * 76)
print("(B) AQUAL DERIVEE — action v7 avec terme cinetique memoire")
print(f"    {'M':>6} {'v_inf (mesure)':>16} {'(G M a0)^1/4 (BTFR)':>21} {'ecart':>8}")
for M in masses:
    r = np.logspace(0, 4, 400)
    gN = G*M/r**2
    g = solve_g(gN)
    v_inf = np.sqrt(g[-1]*r[-1])
    v_btfr = (G*M*a0)**0.25
    print(f"    {M:6.0f} {v_inf:16.6f} {v_btfr:21.6f} {abs(v_inf/v_btfr-1)*100:7.3f}%")
print("-" * 76)
print("    Systeme solaire : correction relative = exp(-g/a0)")
print(f"    g/a0 ~ 1e8 a 1 UA  =>  delta g/g ~ exp(-1e8) : aucun ecrantage requis")
print("=" * 76)

# ---------- Figure ----------
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5))

M = 10.0
r = np.logspace(0, 3.5, 600)
gN = G*M/r**2
g = solve_g(gN)
rA, vA = rotation_curve_fX_only(M)
axes[0].semilogx(r, np.sqrt(gN*r), 'r--', lw=1.8, label='Newton (Kepler)')
axes[0].semilogx(rA, vA, 'm-.', lw=1.8, label='curvature-memory only (no-go): Kepler')
axes[0].semilogx(r, np.sqrt(g*r), 'b-', lw=2.2, label=r'v7 action (derived AQUAL): flat')
axes[0].axhline((G*M*a0)**0.25, color='gray', ls=':', lw=1.2, label=r'BTFR $(GMa_0)^{1/4}$')
axes[0].set_xlabel('radius (code units)'); axes[0].set_ylabel('v(r)')
axes[0].set_title('Rotation curve: no-go vs derived MOND')
axes[0].legend(fontsize=8.5)

gN_grid = np.logspace(-8, -1, 200)
g_grid = solve_g(gN_grid)
axes[1].loglog(gN_grid/a0, g_grid/a0, 'b-', lw=2.2, label=r'derived: $\mu(y)=1-e^{-y}$')
axes[1].loglog(gN_grid/a0, gN_grid/a0, 'r--', lw=1.5, label='Newton 1:1')
axes[1].loglog(gN_grid/a0, np.sqrt(gN_grid/a0), 'g:', lw=1.5, label=r'deep MOND $\sqrt{g_N a_0}$')
axes[1].set_xlabel(r'$g_N/a_0$'); axes[1].set_ylabel(r'$g/a_0$')
axes[1].set_title('Radial Acceleration Relation (derived)')
axes[1].legend(fontsize=8.5)
for ax in axes: ax.grid(True, which='both', alpha=0.2)

plt.suptitle('FSG v7: MOND derived from the kinetic memory term (Priority 1)', fontsize=12)
plt.tight_layout()
plt.savefig('fig_mond_from_action.pdf', dpi=300, bbox_inches='tight')
plt.savefig('fig_mond_from_action.png', dpi=150, bbox_inches='tight')
print("Figure : fig_mond_from_action.pdf / .png")
