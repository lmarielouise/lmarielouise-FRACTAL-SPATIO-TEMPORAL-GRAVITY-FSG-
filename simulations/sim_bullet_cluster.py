import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# GÉNÉRATEUR FIGURE : BULLET CLUSTER (PDF)
# ==========================================

# Paramètres FDTD
L = 100.0
Nx = 400
dx = L / Nx
Nt = 601 
c_grav = 1.0
dt = 0.5 * dx / c_grav

v_gas = 0.6 * c_grav
t_collision = 300

x = np.linspace(0, L, Nx)
phi = np.zeros(Nx)
phi_prev = np.zeros(Nx)
source = np.zeros(Nx)
m_eff = 0.1 

# Captures
snapshots = []
times_to_capture = [100, 250, 400, 600]

# Boucle de Simulation
for n in range(Nt):
    # Dynamique du Gaz
    if n < t_collision:
        pos_gas = 20.0 + v_gas * (n * dt)
    else:
        pos_gas = 20.0 + v_gas * (t_collision * dt) # Arrêt brutal

    source = np.exp(-((x - pos_gas)**2) / 4.0)
    
    # Laplacien
    laplacian = np.zeros(Nx)
    laplacian[1:-1] = (phi[2:] - 2*phi[1:-1] + phi[:-2]) / dx**2
    
    # Équation d'Onde FSG
    accel = (c_grav**2 * laplacian) - (m_eff**2 * phi) + source
    
    if n == 0:
        phi_next = phi
    else:
        phi_next = 2*phi - phi_prev + (dt**2 * accel)
    
    phi_prev = phi.copy()
    phi = phi_next.copy()
    
    # Conditions limites (Zéro aux bords)
    phi[0] = 0; phi[-1] = 0

    if n in times_to_capture:
        snapshots.append((phi.copy(), source.copy(), n))

# Création du Graphique Vectoriel
plt.figure(figsize=(12, 8)) # Taille standard A4/Article

for i, (p, s, t) in enumerate(snapshots):
    plt.subplot(4, 1, i+1)
    
    p_norm = p / np.max(p) if np.max(p) > 0 else p
    
    # Calcul décalage
    peak_gas = x[np.argmax(s)]
    peak_grav = x[np.argmax(p)]
    shift = peak_grav - peak_gas
    
    status = "PRE-COLLISION" if t < t_collision else "POST-COLLISION"
    
    # Esthétique "Papier Scientifique"
    plt.plot(x, s, 'r-', linewidth=2, label='Gas (Baryons)' if i==0 else "")
    plt.plot(x, p_norm, 'b-', linewidth=2, label='FSG Potential (Gravity)' if i==0 else "")
    plt.fill_between(x, p_norm, color='blue', alpha=0.1)
    
    plt.text(1, 0.8, f"t={t} ({status})", fontsize=10, fontweight='bold')
    plt.text(1, 0.5, f"Offset: {shift:.1f}", fontsize=10, color='black')

    plt.ylim(-0.2, 1.3)
    plt.yticks([]) # On enlève les axes Y pour la clarté
    plt.grid(True, alpha=0.2)
    
    if i == 0:
        plt.legend(loc='upper right', frameon=True)
    
    if i < 3:
        plt.xticks([]) # On enlève les axes X sauf pour le dernier
    else:
        plt.xlabel("Spatial Coordinate (Arbitrary Units)")

plt.tight_layout()

# SAUVEGARDE EN PDF (Vectoriel)
filename = "fig_bullet_simulation.pdf"
plt.savefig(filename, format='pdf', bbox_inches='tight')
print(f"Fichier '{filename}' généré avec succès.")
plt.close()