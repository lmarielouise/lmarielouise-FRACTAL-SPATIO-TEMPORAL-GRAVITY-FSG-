import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# --- 1. CONSTANTES COSMOLOGIQUES ---
H0 = 67.4       # km/s/Mpc
h = H0 / 100.0
Omega_b = 0.05  # Baryons (5%)
Omega_dm = 0.26 # Matière Noire (26%) - Utilisé pour le modèle standard uniquement
Omega_m = 0.31  # Total Matière
Omega_r = 9e-5  # Radiation
c_light = 3e5   # km/s

# Échelle de l'horizon aujourd'hui (Mpc)
L_H = c_light / H0 

# --- 2. FONCTIONS DE FOND (Hubble & Expansion) ---
def hubble(a):
    # H(a) / H0
    return np.sqrt(Omega_r/a**4 + Omega_m/a**3 + (1-Omega_m-Omega_r))

# --- 3. MODÈLE DE GRAVITÉ FSG (Le Cœur du test) ---
def G_eff_fsg(k, a):
    """
    C'est ici que réside votre théorie.
    Le propagateur fractal G ~ 1/k^3 implique un renforcement IR.
    Forme effective : G_eff = G_Newton * (1 + (k_scale / k)^n)
    """
    # Échelle de transition (basée sur a0 ou Horizon)
    k_trans = 1.0 / (L_H * 0.1) # Paramètre ajustable (~ échelle des amas)
    
    # Boost fractal (n=1 correspond à passer de 1/k^2 à 1/k^3)
    boost = 1.0 + 0.15 * (k_trans / k)**1.0 
    
    # On sature le boost pour éviter les infinis numériques à k=0
    return np.minimum(boost, 50.0)

# --- 4. ÉQUATIONS DIFFÉRENTIELLES (Perturbations) ---
def derivatives(y, a, k, model_type):
    # y = [delta, d_delta/da] (Perturbation de densité et sa dérivée)
    delta = y[0]
    d_delta = y[1]
    
    H_val = hubble(a)
    
    # Vitesse du son dans le fluide baryon-photon (approx)
    cs2 = 1.0 / 3.0 # c_s^2 = 1/3 (relativiste)
    
    # Terme de friction (Expansion de Hubble)
    friction = (3 / (2 * a)) * d_delta # Simplification régime matière
    
    # Terme de Pression (Oscillations Acoustiques)
    # k_phys = k / a (le mode physique s'étire)
    pressure = (cs2 * k**2 / (a**2 * H0**2)) * delta # Normalisé
    
    # Terme de Gravité (Source)
    # C'est ici qu'on change la physique !
    if model_type == "LCDM":
        # Avec Matière Noire : Densité totale (Baryons + DM)
        density_factor = (Omega_b + Omega_dm) / a**3
        G_factor = 1.0
        
    elif model_type == "GR_NoDM":
        # Sans Matière Noire : Juste les Baryons
        density_factor = Omega_b / a**3
        G_factor = 1.0
        
    elif model_type == "FSG":
        # FSG : Juste les Baryons MAIS Gravité Boostée
        density_factor = Omega_b / a**3
        G_factor = G_eff_fsg(k, a)
    
    # L'équation source : 4 * pi * G * rho * delta
    # Facteur 1.5 * Omega * ... vient des équations de Friedmann linéarisées
    gravity = 1.5 * (density_factor / hubble(a)**2) * delta * G_factor

    # Équation d'ordre 2 : delta'' + friction - gravity + pressure = 0
    # Donc d_delta' = gravity - pressure - friction
    d2_delta = gravity - pressure - friction
    
    return [d_delta, d2_delta]

# --- 5. SIMULATION ---
def run_simulation():
    print(">> Calcul des perturbations cosmologiques en cours...")
    
    # Plage de modes k (en h/Mpc) - Correspond aux multipoles l
    k_modes = np.logspace(-3, 0, 200) # De très grande échelle à petite échelle
    
    # Temps : de a=1e-4 (Radiation) à a=1e-3 (Recombinaison)
    # C'est là que le CMB se fige.
    a_span = np.logspace(-4, -3, 100)
    
    # Stockage des spectres de puissance finaux P(k) ~ delta^2
    pk_lcdm = []
    pk_nodm = []
    pk_fsg = []
    
    for k in k_modes:
        y0 = [1.0, 0.0] # Perturbation initiale unitaire
        
        # 1. LCDM
        sol_lcdm = odeint(derivatives, y0, a_span, args=(k, "LCDM"))
        pk_lcdm.append(sol_lcdm[-1, 0]**2)
        
        # 2. GR Sans DM
        sol_nodm = odeint(derivatives, y0, a_span, args=(k, "GR_NoDM"))
        pk_nodm.append(sol_nodm[-1, 0]**2)
        
        # 3. FSG (Fractal)
        sol_fsg = odeint(derivatives, y0, a_span, args=(k, "FSG"))
        pk_fsg.append(sol_fsg[-1, 0]**2)

    # Conversion k -> Multipole l (approx simple l ~ k * distance)
    # Distance angulaire approx ~ 14000 Mpc
    l_modes = k_modes * 14000 
    
    # --- 6. VISUALISATION ---
    plt.figure(figsize=(10, 6))
    
    # Normalisation pour l'affichage (style CMB)
    # D_l = l(l+1) * C_l / 2pi  ~ k^3 * P(k)
    factor = l_modes * (l_modes + 1)
    
    plt.plot(l_modes, np.array(pk_lcdm) * factor, 'k-', label=r'$\Lambda$CDM (Ref: With Dark Matter)', lw=2)
    plt.plot(l_modes, np.array(pk_nodm) * factor, 'g--', label='GR No Dark Matter (Collapsed Peaks)')
    plt.plot(l_modes, np.array(pk_fsg) * factor, 'r-', label='FSG (Fractal Boost, No DM)', lw=2.5)
    
    plt.xlabel(r'Multipole Moment $\ell$')
    plt.ylabel(r'Power Proxy $\ell(\ell+1)P(k)$')
    plt.title('Preliminary Boltzmann Check: Acoustic Peaks Amplitude')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 1500)
    plt.yscale('log') # Log pour bien voir les différences d'ordre de grandeur
    
    text = "RESULTAT PRELIMINAIRE :\nLe script integre l'equation\ndes perturbations lineaires.\n"
    text += "Si la courbe Rouge remonte\nvers la Noire, le concept FSG survit."
    plt.text(10, min(np.array(pk_lcdm)*factor), text, fontsize=9, bbox=dict(facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig("check_cmb_physics.pdf")
    print(">> Analyse terminée. Voir 'check_cmb_physics.pdf'")

if __name__ == "__main__":
    run_simulation()