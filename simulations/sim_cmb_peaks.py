# ==============================================================================
# FSG — Positions des pics acoustiques CMB (calcul honnête, sans code Boltzmann)
# ==============================================================================
# Correction 4 : ce qui EST calculable sans code Boltzmann complet, c'est la
# POSITION des pics, via l'échelle angulaire acoustique :
#
#   l_n ~ n * pi * chi(z_dec) / r_s(z_dec)
#
# où r_s est l'horizon sonore à la recombinaison et chi la distance comobile.
# La modification FSG de l'énergie sombre tardive (w0 < -1) ne change pas r_s
# (physique pré-recombinaison inchangée) mais déplace chi(z_dec).
# On compare LCDM (w = -1) et FSG (w0 = -1.1 constant, cas conservateur).
#
# Les AMPLITUDES des pics, elles, requièrent la hiérarchie de Boltzmann
# complète (couplage Thomson photon-baryon, potentiels métriques, neutrinos) :
# hors de portée de ce script, voir Section "Future Work" du papier.
# ==============================================================================

import numpy as np
from scipy.integrate import quad

# --- Paramètres (Planck 2018) ---
H0 = 67.4                    # km/s/Mpc
h = H0 / 100.0
Omega_b = 0.0493
Omega_m = 0.315
Omega_g = 2.47e-5 / h**2     # photons
N_eff = 3.046
Omega_r = Omega_g * (1 + 0.2271 * N_eff)
Omega_L = 1.0 - Omega_m - Omega_r
c = 2.998e5                  # km/s
z_dec = 1089.9

def E(z, w0=-1.0):
    """H(z)/H0 ; l'énergie sombre a une EoS constante w0."""
    return np.sqrt(Omega_r * (1 + z)**4 + Omega_m * (1 + z)**3
                   + Omega_L * (1 + z)**(3 * (1 + w0)))

def sound_speed(z):
    """Vitesse du son du fluide photon-baryon (en unités de c)."""
    R_b = (3.0 * Omega_b) / (4.0 * Omega_g) / (1 + z)
    return 1.0 / np.sqrt(3.0 * (1.0 + R_b))

def sound_horizon(w0=-1.0):
    """r_s(z_dec) en Mpc (comobile)."""
    integrand = lambda z: sound_speed(z) / E(z, w0)
    val, _ = quad(integrand, z_dec, 1e7, limit=400)
    return (c / H0) * val

def comoving_distance(w0=-1.0):
    """chi(z_dec) en Mpc (comobile)."""
    integrand = lambda z: 1.0 / E(z, w0)
    val, _ = quad(integrand, 0, z_dec, limit=400)
    return (c / H0) * val

print("=" * 72)
print("POSITIONS DES PICS ACOUSTIQUES : LCDM vs FSG (calcul géométrique)")
print("=" * 72)

results = {}
for name, w0 in [("LCDM  (w = -1.0)", -1.0), ("FSG   (w0 = -1.1)", -1.1)]:
    rs = sound_horizon(w0)
    chi = comoving_distance(w0)
    l1 = np.pi * chi / rs
    results[name] = (rs, chi, l1)
    print(f"{name} :  r_s = {rs:7.2f} Mpc   chi = {chi:8.1f} Mpc   l_1 ~ {l1:6.1f}")

l1_lcdm = results["LCDM  (w = -1.0)"][2]
l1_fsg = results["FSG   (w0 = -1.1)"][2]
shift = 100 * (l1_fsg / l1_lcdm - 1)

print("-" * 72)
print(f"Shift FSG de la position des pics : {shift:+.2f} %")
print(f"Precision Planck sur theta_star   : ~0.03 %")
print("-" * 72)
print("NOTES :")
print(" * l_1 calcule ici (~302) > 220 observe : l'ecart vient de la phase")
print("   acoustique (~0.27 pi) negligee dans l'approx l = pi*chi/r_s ;")
print("   il s'annule dans le RAPPORT FSG/LCDM, qui est la quantite testee.")
print(" * r_s est identique dans les deux modeles : la modification FSG est")
print("   infrarouge/tardive et ne touche pas la physique pre-recombinaison.")
print(" * Le shift de +0.8% A H0 FIXE est absorbe par la degenerescence")
print("   geometrique (w, H0) du CMB : un w0 < -1 se compense par un H0 plus")
print("   eleve a theta* constant (direction favorable a la tension de Hubble).")
print("   La geometrie du CMB seule n'exclut donc pas FSG ; le vrai test est")
print("   l'AMPLITUDE des pics, qui requiert un code Boltzmann non-local.")
print("=" * 72)
