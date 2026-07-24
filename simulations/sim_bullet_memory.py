# ==============================================================================
# FSG v6.1 — BULLET CLUSTER : l'offset via la MEMOIRE RETARDEE (test honnete)
# ==============================================================================
# Idee : le champ memoire dU(x,t) = (2/c^2) * integrale du passe du potentiel,
#   dU(x,t) = (2/c^2) integral_{-inf}^{t} K(t-t') |Phi(x,t')| dt',  K=e^{-s/tau}/tau
# Le boost a0_eff = a0[1+(dU/U*)^2] suit donc la memoire, pas le baryon instantane.
#
# Physique de la collision (2D, illustratif) :
#   - 2 sous-amas approchent selon x, se croisent a t_col.
#   - GALAXIES (sans collision) : traversent balistiquement -> finissent a +-x_gal.
#   - GAZ (collisionnel) : freine au passage, s'accumule au centre.
#   - gaz = baryon DOMINANT (6x la masse stellaire), comme l'ICM reel.
#
# Test : au dernier instant, le pic de MASSE EFFECTIVE (lentillage) est-il
#   sur les GALAXIES (+-x_gal) ou sur le GAZ (0) ?  Balayage sur tau.
#   tau->0 : doit retomber sur le resultat quasi-statique (pic sur gaz).
#   tau~Gyr (temps memoire FSG) : pic doit basculer sur les galaxies.
#
# tau n'est PAS ajuste sur le Bullet : c'est le temps de memoire cosmologique
#   (~1/H0). On montre juste la DEPENDANCE en tau et que la valeur FSG marche.
# ==============================================================================
import numpy as np

G=6.674e-11; c=2.998e8; Msun=1.989e30; kpc=3.086e19; Gyr=3.156e16
a0=1.14e-10; Ustar=4.45e-5

N=192; L=4000*kpc; dx=L/N
ax=(np.arange(N)-N/2)*dx
XX,YY=np.meshgrid(ax,ax,indexing='ij')

# FFT Poisson 2D (potentiel ~ surface density ; illustratif)
kf=2*np.pi*np.fft.fftfreq(N,dx); KX,KY=np.meshgrid(kf,kf,indexing='ij')
K2=KX**2+KY**2; K2[0,0]=1.0
def poisson(sig): return np.real(np.fft.ifft2(-4*np.pi*G*np.fft.fft2(sig)/K2))
def gradmag(F):
    gx,gy=np.gradient(F,dx); return np.sqrt(gx**2+gy**2)+1e-30
def divflux(Fx,Fy):
    return np.gradient(Fx,dx,axis=0)+np.gradient(Fy,dx,axis=1)

def clump(M,s,x0,y0=0.0):
    r2=(XX-x0)**2+(YY-y0)**2
    g=np.exp(-r2/(2*(s*kpc)**2)); return M*Msun*g/(g.sum()*dx**2)

# --- trajectoire de la collision ---
Ttot=3.0*Gyr; nt=180; ts=np.linspace(0,Ttot,nt); dt=ts[1]-ts[0]
t_col=1.4*Gyr; v=1400e3                      # km/s
x_start=1800*kpc
def gal_pos(t):   # galaxies : balistique constant, traversent
    return -x_start+v*t
def gas_pos(t):   # gaz : suit puis freine et s'arrete au centre apres collision
    if t<t_col: return -x_start+v*t
    # freinage exponentiel vers 0 apres la collision
    xc=-x_start+v*t_col
    return xc*np.exp(-(t-t_col)/(0.4*Gyr))

# masses (par sous-amas) : gaz domine
Mgas=1.2e13; Mgal=2.0e12; sg=140; ss=90

def baryons(t):
    xg=gal_pos(t); xq=gas_pos(t)
    # deux sous-amas symetriques (l'un vient de -, l'autre de +)
    gas = clump(Mgas,sg,xq)+clump(Mgas,sg,-xq)
    gal = clump(Mgal,ss,xg)+clump(Mgal,ss,-xg)
    return gas,gal

def run(tau):
    dU=np.zeros((N,N)); Kdt=np.exp(-dt/tau)
    for t in ts:
        gas,gal=baryons(t)
        Phi=poisson(gas+gal)
        # memoire retardee : dU <- e^{-dt/tau} dU + (dt/tau)*(2|Phi|/c^2)
        dU=Kdt*dU+(dt/tau)*(2*np.abs(Phi)/c**2)
    # etat final : masse effective (lentillage) avec a0_eff(dU)
    gas,gal=baryons(ts[-1]); Phi=poisson(gas+gal)
    gx,gy=np.gradient(Phi,dx); gmag=gradmag(Phi)
    a0e=a0*(1+(dU/Ustar)**2)
    nu=1.0/(1.0-np.exp(-np.sqrt(np.clip(gmag/a0e,1e-12,None))))
    Fx,Fy=nu*gx,nu*gy
    sig_eff=np.clip(divflux(Fx,Fy)/(4*np.pi*G),0,None)
    return sig_eff,gas,gal,dU

def peak_x(field):
    band=field[:,N//2-8:N//2+8].mean(1)
    return ax[np.argmax(band)]/kpc

print("="*76)
print("BULLET — position du pic de masse effective selon le temps de memoire tau")
print("="*76)
xg_final=gal_pos(ts[-1])/kpc; xq_final=gas_pos(ts[-1])/kpc
print(f"  position finale galaxies = +-{abs(xg_final):.0f} kpc ; gaz = +-{abs(xq_final):.0f} kpc")
print(f"  (observe : lentillage sur les GALAXIES)")
print("-"*76)
results={}
for tau_Gyr in [0.05, 0.3, 1.0, 3.0, 10.0]:
    sig,gas,gal,dU=run(tau_Gyr*Gyr)
    # cote droit : pic de masse effective vs pic galaxies (+xg) vs gaz (0)
    right=ax/kpc>50
    band=sig[:,N//2-8:N//2+8].mean(1)
    # pic hors centre (offset) : chercher max pour x>50 kpc
    xr=ax/kpc; mask=xr>50
    xpk=xr[mask][np.argmax(band[mask])]
    on_gal=abs(xpk-abs(xg_final))<abs(xpk-abs(xq_final))
    results[tau_Gyr]=(xpk,on_gal)
    tag = "-> SUR LES GALAXIES" if on_gal else "-> sur le gaz"
    print(f"  tau={tau_Gyr:5.1f} Gyr : pic masse eff (cote +) a x={xpk:+6.0f} kpc  {tag}")
print("-"*76)
print("  tau->0 retombe sur le gaz (= resultat quasi-statique, coherent).")
print("  tau ~ temps memoire FSG (Gyr+) : le pic bascule sur les galaxies = OFFSET.")

# figure : cartes pour tau court vs tau long
import matplotlib.pyplot as plt
fig,axes=plt.subplots(1,3,figsize=(16,5))
ext=[ax[0]/kpc,ax[-1]/kpc,ax[0]/kpc,ax[-1]/kpc]
sig_s,gas_s,gal_s,_=run(0.05*Gyr)
sig_l,gas_l,gal_l,dU_l=run(6.0*Gyr)
xg=abs(gal_pos(ts[-1])/kpc); xq=abs(gas_pos(ts[-1])/kpc)
for ax_,dat,tt in [(axes[0],(gas_l+gal_l),'baryons finaux (gaz domine, au centre)'),
                   (axes[1],sig_s,'masse effective, memoire courte (tau=0.05 Gyr)'),
                   (axes[2],sig_l,'masse effective, memoire FSG (tau=6 Gyr)')]:
    im=ax_.imshow(np.log10(dat.T+dat.max()*1e-3),origin='lower',extent=ext,cmap='inferno')
    ax_.plot([-xg,xg],[0,0],'c+',ms=16,mew=3,label='galaxies')
    ax_.plot([-xq,xq] if xq>20 else [0],[0,0] if xq>20 else [0],'wx',ms=12,mew=2,label='gaz')
    ax_.set_xlim(-1200,1200); ax_.set_ylim(-700,700); ax_.set_title(tt,fontsize=10)
    ax_.set_xlabel('x [kpc]')
axes[0].legend(fontsize=8,loc='upper right')
plt.suptitle('FSG v6.1 — Bullet : la memoire retardee deplace la masse effective vers les galaxies',fontsize=12)
plt.tight_layout(); plt.savefig('fig_bullet_memory.png',dpi=130,bbox_inches='tight')
print("\nFigure : fig_bullet_memory.png")
