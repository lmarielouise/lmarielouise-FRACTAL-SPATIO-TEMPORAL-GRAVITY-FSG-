import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

G = 1.0
M = 10.0
Event_Horizon = 2 * G * M  
N_sondes = 20
angles = np.linspace(0, 2*np.pi, N_sondes)
radius_init = 15.0

pos = np.zeros((N_sondes, 2))
pos[:, 0] = radius_init * np.cos(angles)
pos[:, 1] = radius_init * np.sin(angles)
vel = np.zeros((N_sondes, 2))
# Vitesse tangentielle pour orbiter
vel[:, 0] = -0.6 * np.sin(angles) - 0.2 * np.cos(angles) 
vel[:, 1] = 0.6 * np.cos(angles) - 0.2 * np.sin(angles)

fig, ax = plt.subplots(figsize=(8, 8))
plt.style.use('dark_background')
circle = plt.Circle((0, 0), Event_Horizon, color='black', zorder=10)
horizon = plt.Circle((0, 0), Event_Horizon, color='white', fill=False, linestyle='--')
ax.add_artist(circle)
ax.add_artist(horizon)
ax.set_xlim(-20, 20)
ax.set_ylim(-20, 20)
ax.set_aspect('equal')
ax.set_title("Trajectoires Trou Noir (Potentiel FSG)", color='cyan')

scat = ax.scatter(pos[:, 0], pos[:, 1], s=10, c='cyan')

def update(frame):
    global pos, vel
    dt = 0.1
    r = np.sqrt(pos[:, 0]**2 + pos[:, 1]**2)
    r[r < 0.1] = 0.1 
    
    # FSG: Force en 1/r (Log potentiel), pas de singularite violente
    acc_mag = G * M / (r**1.0)
    
    ax_grav = -acc_mag * (pos[:, 0] / r)
    ay_grav = -acc_mag * (pos[:, 1] / r)
    
    vel[:, 0] += ax_grav * dt
    vel[:, 1] += ay_grav * dt
    pos[:, 0] += vel[:, 0] * dt
    pos[:, 1] += vel[:, 1] * dt
    
    scat.set_offsets(pos)
    return scat,

print("Simulation Trou Noir en cours...")
anim = FuncAnimation(fig, update, frames=300, interval=20, blit=True)
plt.show()
