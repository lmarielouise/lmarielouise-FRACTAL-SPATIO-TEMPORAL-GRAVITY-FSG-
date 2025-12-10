import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# MODE SELECTION
MODE_FSG = True # Mettre False pour Newton

N = 100          
G = 0.5          
Softening = 1.0  
BoxSize = 30.0   
Duration = 200   

pos = np.random.rand(N, 2) * BoxSize - (BoxSize/2)
vel = np.zeros((N, 2)) 
mass = np.ones(N)      

fig, ax = plt.subplots(figsize=(8, 8))
plt.style.use('dark_background')
ax.set_xlim(-BoxSize, BoxSize)
ax.set_ylim(-BoxSize, BoxSize)
ax.axis('off')

particles, = ax.plot([], [], '.', color='cyan', markersize=4, alpha=0.7)
title_mode = "FSG (Fractal 2D Longue Portee)" if MODE_FSG else "Newton (Standard 3D)"
ax.set_title(f"Formation des Structures\nMode: {title_mode}", color='white')

def compute_accelerations(pos):
    acc = np.zeros_like(pos)
    for i in range(N):
        for j in range(N):
            if i != j:
                dx = pos[j,0] - pos[i,0]
                dy = pos[j,1] - pos[i,1]
                r = np.sqrt(dx**2 + dy**2 + Softening**2)
                
                if MODE_FSG:
                    force_mag = G * mass[j] / (r**1.0) # Force en 1/r (Fractal)
                else:
                    force_mag = G * mass[j] / (r**2.0) # Force en 1/r^2 (Newton)
                
                acc[i,0] += force_mag * (dx/r)
                acc[i,1] += force_mag * (dy/r)
    return acc

def update(frame):
    global pos, vel
    acc = compute_accelerations(pos)
    vel += acc * 0.05
    pos += vel * 0.05
    particles.set_data(pos[:,0], pos[:,1])
    return particles,

print(f"Simulation Toile Cosmique en cours ({title_mode})...")
anim = FuncAnimation(fig, update, frames=Duration, interval=30, blit=True)
plt.show()
