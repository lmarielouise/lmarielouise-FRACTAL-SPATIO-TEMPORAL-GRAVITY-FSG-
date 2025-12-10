import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

# ==========================================
# 1. PARAMÈTRES PHYSIQUES (Unités Astronomiques & Années)
# ==========================================
G = 4 * np.pi**2  # Gravité ajustée pour UA et Années

# Données [Distance (UA), Vitesse (UA/an), Couleur, Taille, Nom]
planets_data = [
    [0.39, 10.0, 'gray', 10, 'Mercure'],
    [0.72, 7.3, 'orange', 15, 'Vénus'],
    [1.00, 6.28, 'blue', 15, 'Terre'],
    [1.52, 5.0, 'red', 12, 'Mars'],
    [5.20, 2.7, 'brown', 40, 'Jupiter'],
    [9.54, 2.0, 'gold', 35, 'Saturne']
]
n_planets = len(planets_data)

# Initialisation
pos = np.zeros((n_planets, 2))
vel = np.zeros((n_planets, 2))

# On place les planètes sur l'axe X avec leur vitesse tangentielle sur Y
for i, p in enumerate(planets_data):
    pos[i, 0] = p[0]  # Distance r
    vel[i, 1] = p[1]  # Vitesse v

# ==========================================
# 2. PRÉPARATION GRAPHIQUE
# ==========================================
fig, ax = plt.subplots(figsize=(10, 10))
fig.patch.set_facecolor('black') # Fond noir (fenêtre)
ax.set_facecolor('black')        # Fond noir (graphique)

# Zoom pour bien voir jusqu'à Saturne
ax.set_xlim(-11, 11)
ax.set_ylim(-11, 11)
ax.set_aspect('equal')
ax.axis('off') # On cache les axes moches

# Le Soleil (Statique au centre)
ax.plot(0, 0, 'o', color='yellow', markersize=15, zorder=10)
ax.text(0, -1.5, "Soleil", color='yellow', ha='center', fontsize=8)

# Les Planètes (Scatter plot)
sizes = [p[3]*2 for p in planets_data]
colors = [p[2] for p in planets_data]
scat = ax.scatter(pos[:, 0], pos[:, 1], s=sizes, c=colors, zorder=10)

# Les Traînées (Trails)
trails_x = [[] for _ in range(n_planets)]
trails_y = [[] for _ in range(n_planets)]
lines = [ax.plot([], [], color=p[2], lw=1, alpha=0.6)[0] for p in planets_data]

# Titre
title = ax.text(0.5, 1.05, "FSG Solar System Stability Check (Vainshtein Regime)", 
                transform=ax.transAxes, ha='center', color='white', fontsize=14)

# ==========================================
# 3. MOTEUR PHYSIQUE
# ==========================================
dt = 0.005 # ~1.8 jour par image

def update(frame):
    global pos, vel
    
    # Affichage dans la console pour vérifier que ça ne plante pas
    if frame % 10 == 0:
        sys.stdout.write(f"\rCalcul Frame : {frame}")
        sys.stdout.flush()

    # Calcul Gravité (Newton pur car écrantage total)
    for i in range(n_planets):
        r_vec = pos[i]
        r = np.sqrt(r_vec[0]**2 + r_vec[1]**2)
        
        # F = G*M*m / r^2  -> a = G*M / r^2 (M_soleil = 1)
        acc_mag = G / (r**2)
        
        # Vecteur unitaire vers le centre (-pos / r)
        acc_vec = -pos[i] / r * acc_mag
        
        # Euler semi-implicite (simple et stable pour orbites circulaires)
        vel[i] += acc_vec * dt
        pos[i] += vel[i] * dt
        
        # Mise à jour des traînées
        trails_x[i].append(pos[i, 0])
        trails_y[i].append(pos[i, 1])
        
        # On garde les 100 dernières positions pour la queue
        if len(trails_x[i]) > 100:
            trails_x[i].pop(0)
            trails_y[i].pop(0)
            
        lines[i].set_data(trails_x[i], trails_y[i])

    # Mise à jour des points
    scat.set_offsets(pos)
    return scat, *lines

# ==========================================
# 4. LANCEMENT
# ==========================================
print(">> Génération de l'animation... (Fermez la fenêtre pour arrêter)")

# blit=False est MOINS performant mais PLUS compatible (évite l'écran figé)
ani = FuncAnimation(fig, update, frames=None, interval=20, blit=False)

# Si vous voulez sauvegarder en GIF directement (décommentez les lignes ci-dessous)
# print("\n>> Sauvegarde en cours (solar_system.gif)...")
# ani.save('solar_system.gif', writer='pillow', fps=30)
# print(">> Sauvegarde terminée !")

plt.show()