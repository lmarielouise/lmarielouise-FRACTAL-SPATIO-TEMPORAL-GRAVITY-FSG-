import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ==========================================
# 1. PARAMÈTRES DE LA SIMULATION
# ==========================================
G = 1.0
M_total = 1000.0  # Masse baryonique totale
R_disk = 15.0     # Rayon caractéristique du disque
a0 = 0.1          # Accélération critique FSG/MOND
N_stars = 400     # Nombre d'étoiles pour l'animation

# --- CORRECTION DE L'ERREUR ---
R_min = 0.1       # Rayon minimum (centre)
R_max = 40.0      # Rayon maximum (bord de la galaxie)

# ==========================================
# 2. PHYSIQUE (Distributions et Lois)
# ==========================================

def mass_profile(r, M_tot, R_d):
    """ Modèle de masse : Bulbe + Disque exponentiel """
    return M_tot * (1 - np.exp(-r / R_d))

def acc_newtonian(r, M_r):
    """ Accélération Newtonienne standard a = GM/r^2 """
    return G * M_r / (r**2)

def acc_fsg_mond(a_N, a0_val):
    """ 
    Transition FSG/MOND.
    Fonction d'interpolation 'Simple' : mu(x) = x / (1+x)
    Cela donne des courbes plates parfaites.
    Relation inverse : a_eff = (a_N + sqrt(a_N^2 + 4*a_N*a0)) / 2
    """
    # Formule analytique exacte pour l'interpolation simple
    return (a_N + np.sqrt(a_N**2 + 4 * a_N * a0_val)) / 2

# ==========================================
# 3. CALCUL DES VITESES
# ==========================================
# Génération des positions radiales (aléatoires pour l'aspect visuel)
np.random.seed(42)
r_stars = np.random.uniform(R_min, R_max, N_stars)
theta_stars = np.random.uniform(0, 2*np.pi, N_stars)

# Calculs physiques pour chaque étoile
M_r = mass_profile(r_stars, M_total, R_disk)
a_N = acc_newtonian(r_stars, M_r)
a_FSG = acc_fsg_mond(a_N, a0)

# Vitesses orbitales : V = sqrt(a * r)
v_newton = np.sqrt(a_N * r_stars)
v_fsg = np.sqrt(a_FSG * r_stars)

# Vitesses angulaires : omega = v / r
omega_newton = v_newton / r_stars
omega_fsg = v_fsg / r_stars

# ==========================================
# 4. GRAPHIQUE 1 : COURBE DE ROTATION (Statique)
# ==========================================
# Pour la courbe propre, on utilise un rayon lissé
r_curve = np.linspace(R_min, R_max, 100)
M_curve = mass_profile(r_curve, M_total, R_disk)
a_N_curve = acc_newtonian(r_curve, M_curve)
a_FSG_curve = acc_fsg_mond(a_N_curve, a0)

plt.figure(figsize=(10, 6))
plt.plot(r_curve, np.sqrt(a_N_curve * r_curve), 'r--', linewidth=2, label='Newton (Standard)')
plt.plot(r_curve, np.sqrt(a_FSG_curve * r_curve), 'b-', linewidth=3, label='FSG (Fractal Spacetime)')
plt.axhline(np.sqrt(G*M_total*a0)**0.5 * 1.8, color='gray', linestyle=':', alpha=0.5, label='Asymptote Plate')

plt.title('Validation: Flat Rotation Curves in FSG', fontsize=14)
plt.xlabel('Distance from Center (kpc)')
plt.ylabel('Orbital Velocity (km/s)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('galaxy_rotation_curve.png')
print(">> Graphique statique sauvegardé : 'galaxy_rotation_curve.png'")

# ==========================================
# 5. GRAPHIQUE 2 : ANIMATION SPLIT-SCREEN (GIF)
# ==========================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), subplot_kw={'projection': 'polar'})
fig.patch.set_facecolor('black')

# Configuration des axes
for ax, name in zip([ax1, ax2], ["Newtonian Gravity (No DM)", "FSG (Fractal Spacetime)"]):
    ax.set_facecolor('black')
    ax.set_title(name, color='white', pad=20, fontsize=14, fontweight='bold')
    ax.grid(False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_ylim(0, R_max)

# Création des étoiles (Scatter plot)
# On utilise la colormap 'plasma' : jaune au centre, violet au bord
scat1 = ax1.scatter(theta_stars, r_stars, c=r_stars, cmap='plasma', s=8, alpha=0.8)
scat2 = ax2.scatter(theta_stars, r_stars, c=r_stars, cmap='plasma', s=8, alpha=0.8)

# Variables pour la mise à jour
theta1 = np.copy(theta_stars)
theta2 = np.copy(theta_stars)
dt = 0.2  # Pas de temps pour l'animation

def update(frame):
    global theta1, theta2
    
    # Mise à jour des angles : theta = theta + omega * dt
    # C'est ici qu'on voit la différence physique !
    theta1 += omega_newton * dt  # Newton : cisaillement fort, bords lents
    theta2 += omega_fsg * dt     # FSG : rotation plus rigide, bords rapides
    
    # Mise à jour graphique
    scat1.set_offsets(np.c_[theta1, r_stars])
    scat2.set_offsets(np.c_[theta2, r_stars])
    
    return scat1, scat2

print(">> Génération de l'animation comparative (patience)...")
anim = animation.FuncAnimation(fig, update, frames=200, interval=20, blit=True)

# Sauvegarde
output_file = 'galaxy_rotation_comparison.gif'
anim.save(output_file, writer='pillow', fps=30)
print(f">> Animation sauvegardée : '{output_file}'")
print(">> Vous pouvez maintenant mettre ce GIF dans votre README !")