# ==============================================================================
# FSG v6 (corrected) — Cluster mass budget on the full X-COP sample (12 clusters)
# ==============================================================================
# Data: X-COP high-level products (Eckert/Ghirardini/Ettori 2019).
#   *_hydro_mass.fits : M_tot(<r) [Msun], RADIUS [kpc], header R500 [kpc], M500 [1e14 Msun]
#   *_fgas_profile    : MGAS(<r) [Msun], RADIUS [R/R500]
#   *_mstar.fits      : MSTAR(<r) [Msun], RADIUS [R/R500]   (absent for 4 clusters -> gas only)
#
# IMPORTANT (units, audited): fgas/mstar RADIUS is in R/R500 (from FITS TUNIT),
# converted to kpc via R500 from the header. hydro RADIUS is in kpc.
#
# Test: deficit g_obs/g_pillar(g_bar) [the MOND cluster problem], then closure by
#   a0_eff = a0 [1 + (dU/U*)^n], dU = 2|Phi|/c^2 (Phi from g_obs), single (U*, n).
# ==============================================================================
import numpy as np, glob
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
from scipy.optimize import minimize, minimize_scalar

G=6.674e-11; c=2.998e8; Msun=1.989e30; kpc=3.086e19; a0=1.14e-10
CL=["A85","A644","A1644","A1795","A2029","A2142","A2255","A2319","A3158","A3266","RXC1825","ZW1215"]

def col(fn,h,n):
    with fits.open(fn) as x: return np.array(x[h].data[n],float)
def hdr(fn,k):
    with fits.open(fn) as x: return float(x[1].header[k])
def g_pillar(gN,a): return gN/(1.0-np.exp(-np.sqrt(gN/a)))

def load(cl):
    hm=glob.glob(f"{cl}/*hydro_mass.fits")[0]
    R500=hdr(hm,"R500")*kpc
    r=col(hm,1,"RADIUS")*kpc; Mt=col(hm,1,"M_NFW")*Msun
    fg=glob.glob(f"{cl}/*fgas_profile.fits")[0]
    rg=col(fg,1,"RADIUS")*R500; Mg=col(fg,1,"MGAS")*Msun     # R/R500 -> m
    ms=glob.glob(f"{cl}/*mstar.fits")
    if ms:
        rs=col(ms[0],1,"RADIUS")*R500; Ms=col(ms[0],1,"MSTAR")*Msun
    else:
        rs=rg; Ms=np.zeros_like(Mg)
    Mb=interp1d(rg,Mg,bounds_error=False,fill_value=(Mg[0],Mg[-1]))(r) \
      +interp1d(rs,Ms,bounds_error=False,fill_value=(Ms[0],Ms[-1]))(r)
    go=G*Mt/r**2; gb=G*Mb/r**2
    Phi=-(cumulative_trapezoid(go[::-1],-r[::-1],initial=0)[::-1]+G*Mt[-1]/r[-1])
    return dict(r=r,go=go,gb=gb,dU=2*np.abs(Phi)/c**2,R500=R500)

data={cl:load(cl) for cl in CL}
def mask(d): return (d['r']>80*kpc)&(d['r']<1500*kpc)

# global fit (U*, n)
def resid(p):
    U=10**p[0]; n=p[1]; t=[]
    for cl in CL:
        d=data[cl]; m=mask(d)
        t.append(np.log10(d['go'][m]/g_pillar(d['gb'][m],a0*(1+(d['dU'][m]/U)**n))))
    return np.concatenate(t)
res=minimize(lambda p:np.sqrt(np.mean(resid(p)**2)),[np.log10(4.5e-5),2.0],method='Nelder-Mead')
Ufit,nfit=10**res.x[0],res.x[1]
r2=minimize_scalar(lambda lu:np.sqrt(np.mean(resid([lu,2.0])**2)),bounds=(-5.5,-3.5),method='bounded')
Ustar=10**r2.x

print("="*72)
print("X-COP 12 clusters — memory modulation of the MOND cluster deficit")
print("="*72)
print(f"  free fit : U* = {Ufit:.2e}, n = {nfit:.2f}, RMS = {res.fun:.3f} dex")
print(f"  n=2 fix  : U* = {Ustar:.2e}, RMS = {r2.fun:.3f} dex")
print(f"  vs CLASH lensing calibration (4.9e-5): {abs(Ustar/4.9e-5-1)*100:.0f}%")
sstd=[]; smod=[]
for cl in CL:
    d=data[cl]; m=mask(d)
    sstd.append(np.sqrt(np.mean(np.log10(d['go'][m]/g_pillar(d['gb'][m],a0))**2)))
    a0e=a0*(1+(d['dU'][m]/Ustar)**2)
    smod.append(np.sqrt(np.mean(np.log10(d['go'][m]/g_pillar(d['gb'][m],a0e))**2)))
print(f"  deficit (mean) : {np.mean(sstd):.3f} dex  ->  after modulation : {np.mean(smod):.3f} dex")

# figure
import matplotlib.pyplot as plt
fig,axes=plt.subplots(1,3,figsize=(16.5,5)); cols=plt.cm.turbo(np.linspace(0.05,0.95,12))
ax=axes[0]; gg=np.geomspace(1e-12,1e-8,200)
ax.loglog(gg,g_pillar(gg,a0),'k--',lw=1.5,label='galactic pillar'); ax.loglog(gg,gg,'k:',lw=1)
for cl,cc in zip(CL,cols): d=data[cl]; ax.loglog(d['gb'],d['go'],'o',color=cc,ms=2.5,alpha=0.5)
ax.set_xlabel(r'$g_{bar}$ [m/s$^2$]'); ax.set_ylabel(r'$g_{obs}$ [m/s$^2$]'); ax.set_title('Cluster RAR (X-COP)'); ax.legend(fontsize=8)

ax=axes[1]; rg=np.geomspace(80,1500,40)*kpc; ss=[]; sm=[]
for cl in CL:
    d=data[cl]
    go=interp1d(d['r'],d['go'],bounds_error=False)(rg); gb=interp1d(d['r'],d['gb'],bounds_error=False)(rg)
    dU=interp1d(d['r'],d['dU'],bounds_error=False)(rg)
    ss.append(go/g_pillar(gb,a0)); sm.append(go/g_pillar(gb,a0*(1+(dU/Ustar)**2)))
ss=np.array(ss); sm=np.array(sm)
ax.plot(rg/kpc,np.nanmedian(ss,0),'m-.',lw=2,label='standard (deficit)')
ax.fill_between(rg/kpc,np.nanpercentile(sm,16,0),np.nanpercentile(sm,84,0),color='b',alpha=0.2)
ax.plot(rg/kpc,np.nanmedian(sm,0),'b-',lw=2.5,label='memory modulation')
ax.axhline(1,color='k',ls=':',lw=1); ax.set_xscale('log'); ax.set_ylim(0,4)
ax.set_xlabel('r [kpc]'); ax.set_ylabel(r'$g_{obs}/g_{pred}$ (median $\pm1\sigma$)')
ax.set_title('Deficit closed on 12 clusters'); ax.legend(fontsize=9)

ax=axes[2]
for cl,cc in zip(CL,cols):
    d=data[cl]; m=mask(d)
    ax.loglog(d['dU'][m],(d['go'][m]/g_pillar(d['gb'][m],a0))**2,'o',color=cc,ms=2,alpha=0.35)
dU_g=np.geomspace(3e-6,4e-4,100); ax.loglog(dU_g,1+(dU_g/Ustar)**2,'k-',lw=2.5,label=f'$a_0^{{eff}}/a_0$, U*={Ustar:.1e}')
ax.set_xlabel(r'$\delta U=2|\Phi|/c^2$'); ax.set_ylabel(r'required boost'); ax.set_title('Boost vs memory depth'); ax.legend(fontsize=9)

plt.suptitle('FSG v6 (corrected): memory modulation on the full X-COP cluster sample',fontsize=12)
plt.tight_layout()
plt.savefig('fig_xcop_clusters.pdf',dpi=300,bbox_inches='tight')
plt.savefig('fig_xcop_clusters.png',dpi=140,bbox_inches='tight')
print("Figure : fig_xcop_clusters.pdf / .png")
