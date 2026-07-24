# ==============================================================================
# FSG — Régime de champ fort : trous noirs et séparation d'échelles
# ==============================================================================
# CORRECTION v6/v7 : la version précédente de ce script appliquait la force
# MOND (~1/r) PARTOUT, y compris à l'horizon. C'est physiquement à l'envers.
#
# Physique correcte :
#  * Dans le vide, R = 0 (Schwarzschild) => X = Box^-1 R n'est sourcé que par
#    l'histoire cosmique => f(X) ~ const à l'échelle du trou noir => l'action
#    redevient la RG. Schwarzschild est donc solution exacte de FSG.
#  * La correction fractale ne devient O(1) qu'au RAYON MOND
#        r_M = sqrt(G M / a0),
#    où l'accélération propre du trou noir tombe à a0. Pour tout trou noir,
#    r_M depasse l'horizon de 8 a 11 ordres de grandeur : sur toute la zone
#    de champ fort, FSG = RG. Aucune modification des observables EHT/LIGO.
#  * FSG ne resout PAS la singularite (modification IR, pas UV) : pres de r=0
#    la courbure diverge, X -> grand, f -> 0, RG pure.
#
# Ce script calcule r_s et r_M pour plusieurs masses, la loi d'accélération
# g(r) interpolée (fonction RAR de McGaugh), et produit une figure log-log
# montrant la séparation d'échelles.
# ==============================================================================

import numpy as np
import matplotlib.pyplot as plt

# --- Constantes (SI) ---
G = 6.674e-11          # m^3 kg^-1 s^-2
c = 2.998e8            # m/s
M_sun = 1.989e30       # kg
a0 = 1.2e-10           # m/s^2 (échelle d'accélération FSG/MOND)
AU = 1.496e11          # m
pc = 3.086e16          # m

def r_schwarzschild(M):
    return 2 * G * M / c**2

def r_mond(M):
    """Rayon où g_Newton = a0."""
    return np.sqrt(G * M / a0)

def g_newton(r, M):
    return G * M / r**2

def g_fsg(r, M):
    """Accélération FSG via la fonction d'interpolation RAR (McGaugh 2016) :
       g = g_N / (1 - exp(-sqrt(g_N/a0))).
       g_N >> a0 -> g = g_N (Newton/RG) ; g_N << a0 -> g = sqrt(g_N a0)."""
    gN = g_newton(r, M)
    return gN / (1.0 - np.exp(-np.sqrt(gN / a0)))

# --- Séparation d'échelles pour un échantillon de masses ---
cases = [
    ("Soleil",              1.0),
    ("Trou noir stellaire", 10.0),
    ("Trou noir intermed.", 1.0e3),
    ("Sgr A* (SMBH)",       4.0e6),
    ("M87* (SMBH)",         6.5e9),
]

print("=" * 78)
print("SÉPARATION D'ÉCHELLES FSG EN CHAMP FORT")
print(f"{'Objet':<22}{'r_s':>13}{'r_M':>15}{'r_M/r_s':>14}")
print("-" * 78)
for name, m in cases:
    M = m * M_sun
    rs = r_schwarzschild(M)
    rM = r_mond(M)
    # affichage lisible
    rs_str = f"{rs/1e3:.2e} km" if rs < AU else f"{rs/AU:.2e} AU"
    rM_str = f"{rM/pc:.2e} pc"
    print(f"{name:<22}{rs_str:>13}{rM_str:>15}{rM/rs:>14.1e}")
print("=" * 78)
print("Interprétation : sur toute la zone r_s < r < r_M (8-11 ordres de")
print("grandeur), FSG coïncide avec la RG. La modification n'apparaît qu'au")
print("rayon MOND, bien au-delà de toute observable de champ fort.")
print("=" * 78)

# --- Figure : g(r) pour un trou noir stellaire ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

for ax, (name, m) in zip(axes, [("Trou noir stellaire (10 M_sun)", 10.0),
                                 ("SMBH type M87* (6.5e9 M_sun)", 6.5e9)]):
    M = m * M_sun
    rs = r_schwarzschild(M)
    rM = r_mond(M)
    r = np.logspace(np.log10(rs), np.log10(30 * rM), 2000)

    ax.loglog(r / pc, g_newton(r, M), 'r--', lw=1.8, label=r'Newton / GR ($GM/r^2$)')
    ax.loglog(r / pc, g_fsg(r, M), 'b-', lw=2.2, label='FSG (interpolated)')
    ax.axhline(a0, color='gray', ls=':', lw=1.2, label=r'$a_0$')
    ax.axvline(rs / pc, color='k', lw=1.5, label=r'Horizon $r_s$')
    ax.axvline(rM / pc, color='green', ls='-.', lw=1.5, label=r'MOND radius $r_M$')

    ax.set_xlabel('radius [pc]')
    ax.set_ylabel(r'acceleration $g$ [m/s$^2$]')
    ax.set_title(f"{name}\n" + rf"$r_M/r_s \approx {rM/rs:.0e}$")
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, which='both', alpha=0.2)

plt.suptitle("FSG in the strong-field regime: GR is recovered from the horizon "
             "up to the MOND radius", fontsize=12)
plt.tight_layout()
plt.savefig('fig_black_hole_scales.pdf', dpi=300, bbox_inches='tight')
plt.savefig('fig_black_hole_scales.png', dpi=150, bbox_inches='tight')
print("Figure : fig_black_hole_scales.pdf / .png")
