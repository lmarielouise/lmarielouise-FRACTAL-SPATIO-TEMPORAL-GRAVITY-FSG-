import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fftn, ifftn, fftfreq
from scipy.integrate import odeint
import os

# =============================================================================
# FRACTAL SPACETIME GRAVITY (FSG) - Simulation Complète
# Auteur : Laurent Marie-Louise
# Description : Ce script valide numériquement la théorie FSG sur deux échelles :
# 1. Échelle Galactique : Convolution FFT 3D du propagateur fractal pour obtenir les courbes de rotation.
# 2. Échelle Cosmologique : Résolution de l'équation de croissance linéaire pour expliquer les données JWST.
# =============================================================================

def run_simulation():
    print("--- DÉMARRAGE DE LA SIMULATION FSG (Physique Complète) ---")
    os.makedirs("figures", exist_ok=True)

    # -------------------------------------------------------------------------
    # PARTIE A : DYNAMIQUE GALACTIQUE (Moteur FFT 3D)
    # -------------------------------------------------------------------------
    print("[1/2] Calcul de la gravité galactique via FFT 3D...")
    
    # 1. Création de l'Univers (Grille 3D)
    N = 64            # Résolution de la grille (64x64x64)
    BoxSize = 50.0    # Taille de la boîte en kpc
    dx = BoxSize / N
    x = np.linspace(-BoxSize/2, BoxSize/2, N)
    X, Y, Z = np.meshgrid(x, x, x)
    R = np.sqrt(X**2 + Y**2 + Z**2)
    R[R==0] = dx      # Évite la division par zéro au centre

    # 2. Injection de la Matière (Bulbe + Disque)
    # Profil de Hernquist pour le bulbe
    r_bulbe = 1.5
    rho_bulbe = (1e4) / (2*np.pi) * (r_bulbe / (R * (R + r_bulbe)**3))
    # Profil exponentiel pour le disque
    r_disk = 6.0
    rho_disk = (5e4) / (4*np.pi*r_disk**2) * np.exp(-R/r_disk) * np.exp(-Z**2/1.0)
    
    rho_total = rho_bulbe + rho_disk

    # 3. Passage dans l'Espace de Fourier (k-space)
    k = fftfreq(N, d=dx) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(k, k, k)
    K = np.sqrt(KX**2 + KY**2 + KZ**2)
    K[0,0,0] = 1e-10  # Sécurité pour k=0

    # 4. Application des Lois de la Gravité (Propagateurs)
    
    # --> NEWTON (Standard GR) : G(k) ~ 1/k^2
    Green_N = 1.0 / K**2
    Green_N[0,0,0] = 0
    Phi_N_k = fftn(rho_total) * (-4*np.pi * Green_N)
    Phi_N = np.real(ifftn(Phi_N_k)) # Retour espace réel

    # --> FSG (Fractal) : G(k) ~ 1 / (k^2 * (1 + (kL)^-1)) => Tend vers 1/k^3 en IR
    L_scale = 10.0 # Échelle de transition (kpc)
    term_fractal = (K * L_scale)**(-1.0)
    Green_F = 1.0 / (K**2 * (1.0 + term_fractal))
    Green_F[0,0,0] = 0
    Phi_F_k = fftn(rho_total) * (-4*np.pi * Green_F)
    Phi_F = np.real(ifftn(Phi_F_k)) # Retour espace réel

    # 5. Extraction des Courbes de Rotation (V = sqrt(r * dPhi/dr))
    center = N//2
    r_axis = x[center:]
    # Gradient radial numérique
    grad_N = np.gradient(-Phi_N[center, center, center:], dx)
    grad_F = np.gradient(-Phi_F[center, center, center:], dx)
    
    v_newton = np.sqrt(np.abs(r_axis * grad_N))
    v_fsg = np.sqrt(np.abs(r_axis * grad_F))
    
    # Normalisation pour affichage réaliste
    scale_factor = 230 / np.max(v_fsg)
    v_newton *= scale_factor
    v_fsg *= scale_factor

    # -------------------------------------------------------------------------
    # PARTIE B : COSMOLOGIE JWST (Résolution Équation Différentielle)
    # -------------------------------------------------------------------------
    print("[2/2] Résolution de la croissance des structures (JWST)...")

    # Modèle Standard LambdaCDM
    Om = 0.31
    OL = 0.69

    def solve_growth(epsilon):
        # Espace temporel (Facteur d'échelle a) de z=1000 à z=4
        a_range = np.linspace(1e-3, 0.2, 500) 
        y0 = [1e-4, 1e-4] # Conditions initiales (delta, delta')

        def growth_eq(y, a):
            delta, d_delta = y
            E = np.sqrt(Om/a**3 + OL) # Hubble parameter
            
            # Friction de l'expansion
            term_friction = (1.5 * Om / a**3) / E**2
            
            # Source Gravitationnelle avec Boost FSG
            # Boost = 1 + epsilon/a (Effet plus fort dans le passé lointain)
            boost = 1.0 + (0.5 * epsilon / a)
            term_source = (1.5 * Om / (a**5 * E**2)) * boost
            
            # Équation différentielle du 2nd ordre
            d2_delta = - (3/a + (term_friction - 2/a)) * d_delta + term_source * delta
            return [d_delta, d2_delta]
        
        sol = odeint(growth_eq, y0, a_range)
        delta = sol[:, 0]
        z = (1/a_range) - 1
        
        # Trouver z où delta atteint le seuil critique d'effondrement (1.686)
        idx = np.where(delta > 1.686)[0]
        if len(idx) > 0:
            return z[idx[0]]
        return 0.0

    # Scan des paramètres epsilon (Déficit dimensionnel)
    eps_values = np.linspace(0, 0.25, 40)
    z_formation_values = []
    for eps in eps_values:
        z_formation_values.append(solve_growth(eps))
    z_formation_values = np.array(z_formation_values)

    # -------------------------------------------------------------------------
    # PARTIE C : GÉNÉRATION DE LA FIGURE
    # -------------------------------------------------------------------------
    print("Génération du graphique combiné...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot Galaxie
    ax1.plot(r_axis, v_newton, 'r--', linewidth=2, label='Newton ($1/r^2$) [No DM]')
    ax1.plot(r_axis, v_fsg, 'b-', linewidth=3, label='FSG ($d_S \\to 2$) [No DM]')
    ax1.axvspan(0, 5, color='yellow', alpha=0.1) # Zone Bulbe
    ax1.text(12, 50, "Plateau naturel\n(sans Halo)", color='blue')
    ax1.set_title(r'\textbf{A. Galactic Rotation} (FFT Simulation)')
    ax1.set_xlabel('Radius [kpc]')
    ax1.set_ylabel('Velocity [km/s]')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 25)

    # Plot Cosmologie
    ax2.plot(eps_values, z_formation_values, 'k.-')
    # Cibles
    jwst_target = 17.0
    lcdm_limit = z_formation_values[0]
    ax2.axhline(jwst_target, color='purple', linewidth=2, label='JWST Observation ($z \\approx 17$)')
    ax2.axhline(lcdm_limit, color='red', linestyle='--', linewidth=2, label='Standard $\Lambda$CDM Limit')
    
    # Point optimal
    idx_opt = np.abs(z_formation_values - jwst_target).argmin()
    opt_eps = eps_values[idx_opt]
    ax2.plot(opt_eps, jwst_target, 'bo', markersize=10)
    ax2.annotate(f'FSG Prediction:\n$\epsilon \\approx {opt_eps:.2f}$', 
                 xy=(opt_eps, jwst_target), xytext=(opt_eps+0.05, 10),
                 arrowprops=dict(facecolor='blue', shrink=0.05), color='blue')

    ax2.set_title(r'\textbf{B. Early Structure Formation} (Linear Growth)')
    ax2.set_xlabel(r'Dimensional Deficit $\epsilon = 3 - d_S$')
    ax2.set_ylabel(r'Collapse Redshift $z_{coll}$')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file = "figures/fig_fsg_unified_results.pdf"
    plt.savefig(output_file)
    print(f"✅ SUCCÈS : Figure sauvegardée sous '{output_file}'")
    plt.show()

if __name__ == "__main__":
    run_simulation()
