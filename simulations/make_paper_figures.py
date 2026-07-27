# ==============================================================================
# FSG v6.1 -- regenerate the paper figures that are not produced by the physics
# scripts, directly from the FSG relations. HONEST content only: exact relations,
# real memory-field integration, or clearly-labelled schematics/model curves.
# No fabricated data points. Output to the repo root (found by \graphicspath).
# ==============================================================================
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os
OUT=os.path.join(os.path.dirname(__file__),"..")
G=6.674e-11;c=2.998e8;a0=1.14e-10;Msun=1.989e30;kpc=3.086e19
def nu(y): return 1.0/(1.0-np.exp(-np.sqrt(y)))
def gfsg(gN): return nu(gN/a0)*gN
def save(fig,name): fig.savefig(os.path.join(OUT,name),bbox_inches='tight'); plt.close(fig); print("  ",name)

# 1) BTFR -- exact relation V^4 = G M a0 (no data points)
f,ax=plt.subplots(figsize=(5.2,4)); M=np.logspace(7,12,80)
ax.loglog(M,(G*M*Msun*a0)**0.25/1e3,'b-',lw=2.2,label=r'FSG (exact): $V_{flat}^4=GMa_0$')
ax.set_xlabel(r'$M_{bary}\,[M_\odot]$');ax.set_ylabel(r'$V_{flat}$ [km/s]')
ax.set_title('Baryonic Tully–Fisher (exact FSG relation)');ax.legend();save(f,'fig_btfr.pdf')

# 2) RAR -- FSG relation g_obs(g_bar) (theory)
f,ax=plt.subplots(figsize=(5.2,4)); gN=np.logspace(-13,-8,300)
ax.loglog(gN,gfsg(gN),'b-',lw=2.2,label='FSG'); ax.loglog(gN,gN,'k:',lw=1,label='Newton (1:1)')
ax.axvline(a0,color='r',ls='--',lw=1,label=r'$a_0$')
ax.set_xlabel(r'$g_{bar}\,[m/s^2]$');ax.set_ylabel(r'$g_{obs}\,[m/s^2]$')
ax.set_title('Radial Acceleration Relation (FSG law)');ax.legend();save(f,'fig_rar.pdf')

# 3) Newton vs FSG rotation curve (theory)
f,ax=plt.subplots(figsize=(5.2,4)); Mg=6e10*Msun; r=np.linspace(0.3,40,300)*kpc; gN=G*Mg/r**2
ax.plot(r/kpc,np.sqrt(gN*r)/1e3,'r--',lw=2,label='Newton (baryons only)')
ax.plot(r/kpc,np.sqrt(gfsg(gN)*r)/1e3,'b-',lw=2.2,label='FSG (no dark matter)')
ax.set_xlabel('r [kpc]');ax.set_ylabel('v [km/s]');ax.set_title(r'Rotation curve ($M_b=6\times10^{10}M_\odot$)')
ax.legend();ax.set_ylim(0,260);save(f,'fig_newton_vs_fractal.pdf')

# 9) NGC 6503-like -- FSG MODEL only, explicitly not a data fit
f,ax=plt.subplots(figsize=(5.2,4)); Mg=1.0e10*Msun; r=np.linspace(0.3,22,200)*kpc; gN=G*Mg/r**2
ax.plot(r/kpc,np.sqrt(gfsg(gN)*r)/1e3,'b-',lw=2.2,label='FSG model')
ax.plot(r/kpc,np.sqrt(gN*r)/1e3,'r--',lw=1.5,label='Newton (baryons)')
ax.set_xlabel('r [kpc]');ax.set_ylabel('v [km/s]')
ax.set_title(r'NGC 6503-mass FSG model ($M_b\!\sim\!10^{10}M_\odot$; not a fit)')
ax.legend();ax.set_ylim(0,150);save(f,'fig_ngc6503.pdf')

# 5) spectral dimension d_S: 2 (IR) -> 4 (UV), FSG prediction (schematic)
f,ax=plt.subplots(figsize=(5.2,4)); k=np.logspace(-4,3,400); dS=2+2/(1+(1.0/k)**1.2)
ax.semilogx(k,dS,'b-',lw=2.2);ax.axhline(2,color='gray',ls=':');ax.axhline(4,color='gray',ls=':')
ax.set_ylim(1.8,4.2);ax.set_xlabel('k (probe scale, arb. units)');ax.set_ylabel(r'$d_S$')
ax.set_title('Spectral dimension: 2 (IR) $\\to$ 4 (UV) [FSG, schematic]');save(f,'fig_dimension_spectral.pdf')

# 4) three regimes (schematic)
f,ax=plt.subplots(figsize=(6,3.6))
for x0,x1,d,lab,col in [(1e-4,1e-2,2,'IR fractal\n$d_S\\to2$','#4a90d9'),
                        (1e-2,1,3,'transition\n$d_S\\approx3$','#7ac74f'),
                        (1,1e2,4,'Newton/GR\n$d_S=4$','#e0813a')]:
    ax.fill_between([x0,x1],[d-.07]*2,[d+.07]*2,color=col);ax.text(np.sqrt(x0*x1),d+0.25,lab,ha='center',fontsize=8)
ax.set_xscale('log');ax.set_ylim(1.5,4.7);ax.set_xlabel('scale (large $\\to$ small)');ax.set_ylabel(r'$d_S$')
ax.set_title('The three regimes of FSG (schematic)');save(f,'fig_three_scales.pdf')

# 6) propagator (schematic: IR enhancement)
f,ax=plt.subplots(figsize=(5.2,4)); k=np.logspace(-3,2,300)
ax.loglog(k,1/k**2,'k:',lw=1.5,label=r'Newton $\sim1/k^2$')
ax.loglog(k,1/k**2*(1+0.05/k),'b-',lw=2.2,label='FSG (IR-enhanced, schematic)')
ax.set_xlabel('k');ax.set_ylabel('G(k)');ax.set_title('Gravitational propagator (schematic)');ax.legend();save(f,'fig_propagator.pdf')

# 7,8) w(z), Delta H/H from the REAL memory-field integration (same as sim_friedmann_fsg)
Om,Or=0.315,9.0e-5;OL=1-Om-Or
def E2(N):a=np.exp(N);return Or*a**-4+Om*a**-3+OL
def dlnH(N):a=np.exp(N);return 0.5*(-4*Or*a**-4-3*Om*a**-3)/E2(N)
def rhs(N,Y):
    U,Up,V,Vp=Y;q=dlnH(N);return[Up,-(3+q)*Up-6*(q+2),Vp,-(3+q)*Vp-U/E2(N)]
s=solve_ivp(rhs,[np.log(1e-7),0],[0,0,0,0],dense_output=True,rtol=1e-9)
N=np.linspace(np.log(1e-4),0,2000);U,Up,V,Vp=s.sol(N);a=np.exp(N);z=1/a-1
w=-1-np.gradient(np.log(np.abs(V)+1e-30),N)/3
ODE=OL*np.abs(V)/abs(V[-1]);dHH=np.sqrt((Or*a**-4+Om*a**-3+ODE)/E2(N))-1;m=z<3
f,ax=plt.subplots(figsize=(5.2,4));ax.plot(z[m],w[m],'r-',lw=2);ax.axhline(-1,color='k',ls='--',label=r'$\Lambda$CDM')
ax.set_xlabel('z');ax.set_ylabel('w(z)');ax.set_title('Effective EoS from cosmic memory (upper bound)')
ax.legend();ax.invert_xaxis();save(f,'fig_wz.pdf')
f,ax=plt.subplots(figsize=(5.2,4));ax.plot(z[m],100*dHH[m],'g-',lw=2);ax.axhline(0,color='k',ls='--')
ax.set_xlabel('z');ax.set_ylabel(r'$\Delta H/H$ [%]');ax.set_title(r'Expansion history vs $\Lambda$CDM (upper bound)')
ax.invert_xaxis();save(f,'fig_Hz.pdf')

# 10) summary: two REAL FSG panels (rotation curve + RAR)
f,(A,B)=plt.subplots(1,2,figsize=(9,3.8))
Mg=6e10*Msun;r=np.linspace(0.3,40,300)*kpc;gN=G*Mg/r**2
A.plot(r/kpc,np.sqrt(gfsg(gN)*r)/1e3,'b-',lw=2,label='FSG');A.plot(r/kpc,np.sqrt(gN*r)/1e3,'r--',lw=1.5,label='Newton')
A.set_xlabel('r [kpc]');A.set_ylabel('v [km/s]');A.set_title('A. Flat rotation curve');A.set_ylim(0,260);A.legend(fontsize=8)
gg=np.logspace(-13,-8,300);B.loglog(gg,gfsg(gg),'b-',lw=2);B.loglog(gg,gg,'k:',lw=1);B.axvline(a0,color='r',ls='--')
B.set_xlabel(r'$g_{bar}$');B.set_ylabel(r'$g_{obs}$');B.set_title('B. Radial acceleration relation')
save(f,'fig_fsg_simulation_results.pdf')
print("done -- honest figures only")
