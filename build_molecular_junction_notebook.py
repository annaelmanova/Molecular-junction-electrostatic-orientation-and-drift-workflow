import json, base64
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "Molecular_Junction_Workflow.ipynb"

def md(s): return {"cell_type":"markdown","metadata":{},"source":s.splitlines(True)}
def code(s): return {"cell_type":"code","execution_count":None,"metadata":{},"outputs":[],"source":s.splitlines(True)}

cells=[]
cells.append(md(r'''# Molecular junction tensor-aligned workflow

This is the **single source notebook** for the complete calculation and plot sequence. It does not require project helper scripts. Edit only the user-input cell, then run all cells from top to bottom.

Sequence: electrostatic potential → electric field and derivatives → anisotropic orientation/tensor analysis → isotropic reference → anisotropic force/drift → model comparison → full drift atlas → CSV/PNG/PDF export.
'''))

cells.append(code(r'''# ========================= USER INPUTS: EDIT THIS CELL =========================
from pathlib import Path
import numpy as np

CFG = {
    # User-facing system names (used in layer strips, status text, and reports)
    "molecule_name": "Y6",
    "left_layer_name": "ZnO",
    "interlayer_name": "SAM",
    "active_layer_name": "PM6:Y6",

    # Output
    "output_dir": Path.cwd() / "Molecular_Junction_Workflow_outputs",
    "save_png": True, "save_pdf": True, "dpi": 300,

    # Device/electrostatic parameters
    "active_layer_nm": 100.0,
    "plot_zmax_nm": 50.0,
    "z_step_nm": 0.10,
    "V_bi_V": 0.85,
    "delta_phi_ZnO_V": 0.30,
    "lambda_ZnO_nm": 30.0,
    "lambda_SAM_nm": 30.0,
    "Vapp_sweep_V": [1.0, 0.0, -1.0, -5.0],
    "SAM_sweep_V": [-2.0, -0.7, 0.7, 2.0],
    "lambda_sweep_nm": [10.0, 20.0, 30.0],
    "lambda_SAM_sweep_nm": [10.0, 20.0, 30.0],

    # Transport/orientational parameters
    "temperature_K": 418.15,
    "orientation_temperature_K": 418.15,
    "diffusion_m2_s": [1e-19, 1e-18, 1e-17],
    "times_min": [10.0, 60.0],
    "orientation_samples": 8192,
    "boltzmann_chunk": 128,

    # Geometry-aligned Y6 ground-state dipole and polarizability tensor
    "mu_D": np.array([11.4904, -0.2662, -0.9078], dtype=float),
    "alpha_A3": np.array([[331.972, -24.940, -27.115],
                           [-24.940, 199.188, 6.388],
                           [-27.115, 6.388, 106.596]], dtype=float),
    # Matching quadrupole was not supplied; replace this zero tensor if available.
    "Q_DA": np.zeros((3,3), dtype=float),

    # Plot selections (must be members of the sweeps above)
    "representative_Vapp_V": [0.0, -5.0],
    "derivative_SAM_V": [-0.7, 0.7],
}

print("System:", f'{CFG["left_layer_name"]}/{CFG["interlayer_name"]}/{CFG["active_layer_name"]}',
      "| molecule:", CFG["molecule_name"])
print("Output:", CFG["output_dir"])
print("Cases:", len(CFG["Vapp_sweep_V"])*len(CFG["SAM_sweep_V"]),
      "electrostatic combinations;")
print("transport combinations:", len(CFG["diffusion_m2_s"])*len(CFG["times_min"]))
'''))

cells.append(md('''## 1. Imports, validation, constants, and publication style\n'''))
cells.append(code(r'''import math, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import AutoMinorLocator

required=["molecule_name","left_layer_name","interlayer_name","active_layer_name",
          "Vapp_sweep_V","SAM_sweep_V","diffusion_m2_s","times_min","mu_D","alpha_A3"]
for key in required:
    if key not in CFG: raise ValueError(f"Missing CFG[{key!r}]")
for key in ["molecule_name","left_layer_name","interlayer_name","active_layer_name"]:
    if not isinstance(CFG[key],str) or not CFG[key].strip(): raise ValueError(f"CFG[{key!r}] must be a non-empty name")
if np.asarray(CFG["alpha_A3"]).shape != (3,3): raise ValueError("alpha_A3 must be 3x3")
if not np.allclose(CFG["alpha_A3"],np.asarray(CFG["alpha_A3"]).T): raise ValueError("alpha_A3 must be symmetric")
if np.asarray(CFG["mu_D"]).shape != (3,): raise ValueError("mu_D must have three Cartesian components")

q=1.602176634e-19; kB=1.380649e-23; eps0=8.8541878128e-12
D2CM=3.33564e-30; NM=1e-9
mu=np.asarray(CFG["mu_D"])*D2CM
alpha=4*np.pi*eps0*np.asarray(CFG["alpha_A3"])*1e-30
Q=np.asarray(CFG["Q_DA"])*D2CM*1e-10
alpha_iso=float(np.trace(alpha)/3); mu_mag=float(np.linalg.norm(mu)); Q_iso=float(np.trace(Q)/3)

OUT=Path(CFG["output_dir"]); FIG=OUT/"figures"; PDF_FIG=OUT/"UPDATED_PUBLICATION_PDF_PLOTS"; DAT=OUT/"data"
for p in (OUT,FIG,PDF_FIG,DAT): p.mkdir(parents=True,exist_ok=True)
z_nm=np.arange(0,CFG["plot_zmax_nm"]+CFG["z_step_nm"]*.5,CFG["z_step_nm"]); z=z_nm*NM
d=CFG["active_layer_nm"]*NM

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':7,'axes.labelsize':7,'axes.titlesize':7.4,'legend.fontsize':6.3,
 'xtick.labelsize':6.3,'ytick.labelsize':6.3,'pdf.fonttype':42,'ps.fonttype':42,
 'axes.linewidth':.7,'lines.linewidth':1.25,'figure.facecolor':'white','axes.facecolor':'white',
 'savefig.facecolor':'white','savefig.edgecolor':'white','savefig.transparent':False})
COL_SAM={s:c for s,c in zip(CFG["SAM_sweep_V"],['#1F78B4','#202020','#4D4D4D','#D95F02','#7570B3','#E7298A'])}
COL_V={v:c for v,c in zip(CFG["Vapp_sweep_V"],['#D55E00','#202020','#009E73','#0072B2','#CC79A7','#56B4E9'])}

def device_strip(ax):
    tr=ax.get_xaxis_transform(); y0=-.54; h=.14; xmax=CFG["plot_zmax_nm"]
    for x,w,l,fc,tc in [(-10,4.8,CFG['left_layer_name'],'#dbe8f3','black'),(-5.2,5.2,CFG['interlayer_name'],'#eef3f7','black'),(0,xmax,CFG['active_layer_name'],'#2f73b3','white')]:
        ax.add_patch(Rectangle((x,y0),w,h,transform=tr,facecolor=fc,edgecolor='.25',lw=.8,clip_on=False))
        ax.text(x+w/2,y0+h/2,l,transform=tr,ha='center',va='center',fontsize=5.2,color=tc,clip_on=False)
def style(ax,xlabel=False):
    ax.set_xlim(0,CFG["plot_zmax_nm"]); ax.grid(True,lw=.55,alpha=.35,color='.75');
    ax.xaxis.set_minor_locator(AutoMinorLocator(2)); ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(direction='out',top=False,right=False,width=.7,length=5)
    for s in ['top','right']: ax.spines[s].set_visible(False)
    ax.axvspan(0,2,color='#c9d8e4',alpha=.75,zorder=-10); device_strip(ax)
    if xlabel: ax.set_xlabel('Position z (nm)'); ax.xaxis.set_label_coords(.5,-.76)
def letters(axs):
    for i,a in enumerate(np.asarray(axs).ravel()): a.text(-.15,1.04,f'({chr(97+i)})',transform=a.transAxes,fontweight='bold')
def save(fig,name):
    if CFG["save_png"]: fig.savefig(FIG/f'{name}.png',dpi=CFG["dpi"],bbox_inches='tight')
    if CFG["save_pdf"]:
        fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight')
        fig.savefig(PDF_FIG/f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)
'''))

cells.append(md(r'''## 2. Junction electrostatics: exponential and parabolic closures

For the exponential interfacial model,
\[\phi(z)=-\frac{V_{bi}-V_{app}}d z+\Delta\phi_{ZnO}e^{-z/\lambda_{ZnO}}+\Delta\phi_{SAM}e^{-z/\lambda_{SAM}},\quad E=-\frac{d\phi}{dz}.\]
For comparison, a parabolic depletion closure replaces each exponential by
\[p(z;\lambda)=(1-z/\lambda)^2\ (0\le z<\lambda),\qquad p=0\ (z\ge\lambda).\]
A uniform-field model has no interfacial gradient, an unscreened Coulomb form is singular at the interface, and a self-consistent Poisson/drift-diffusion model requires charge-density, dielectric, injection, and boundary-condition data that are not supplied. Exponential screening and the parabolic depletion reference are therefore transparent analytically differentiable limiting models.
'''))
cells.append(code(r'''def profile(Vapp,SAM,lambda_ZnO_nm=None,lambda_SAM_nm=None):
    lz=(CFG["lambda_ZnO_nm"] if lambda_ZnO_nm is None else lambda_ZnO_nm)*NM
    ls=(CFG["lambda_SAM_nm"] if lambda_SAM_nm is None else lambda_SAM_nm)*NM
    ez=np.exp(-z/lz); es=np.exp(-z/ls)
    phi=-(CFG["V_bi_V"]-Vapp)*z/d+CFG["delta_phi_ZnO_V"]*ez+SAM*es
    E=(CFG["V_bi_V"]-Vapp)/d+CFG["delta_phi_ZnO_V"]/lz*ez+SAM/ls*es
    dE=-CFG["delta_phi_ZnO_V"]/lz**2*ez-SAM/ls**2*es
    d2E=CFG["delta_phi_ZnO_V"]/lz**3*ez+SAM/ls**3*es
    return phi,E,dE,d2E

def profile_parabolic(Vapp,SAM,lambda_ZnO_nm=None,lambda_SAM_nm=None):
    lz=(CFG["lambda_ZnO_nm"] if lambda_ZnO_nm is None else lambda_ZnO_nm)*NM
    ls=(CFG["lambda_SAM_nm"] if lambda_SAM_nm is None else lambda_SAM_nm)*NM
    def term(A,L):
        inside=z<L; u=np.maximum(1-z/L,0)
        return A*u*u,np.where(inside,2*A*u/L,0.0),np.where(inside,-2*A/L**2,0.0),np.zeros_like(z)
    pz,ez,gz,hz=term(CFG["delta_phi_ZnO_V"],lz); ps,es,gs,hs=term(SAM,ls)
    return -(CFG["V_bi_V"]-Vapp)*z/d+pz+ps,(CFG["V_bi_V"]-Vapp)/d+ez+es,gz+gs,hz+hs

profiles={(V,S):profile(V,S) for V in CFG["Vapp_sweep_V"] for S in CFG["SAM_sweep_V"]}
profiles_parabolic={(V,S):profile_parabolic(V,S) for V in CFG["Vapp_sweep_V"] for S in CFG["SAM_sweep_V"]}
rows=[]
for (V,S),(ph,E,G,H) in profiles.items():
    rows.extend(dict(Vapp_V=V,SAM_V=S,z_nm=zz,phi_V=p,E_V_m=e,dE_V_m2=g,d2E_V_m3=h)
                for zz,p,e,g,h in zip(z_nm,ph,E,G,H))
electro=pd.DataFrame(rows); electro.to_csv(DAT/'electrostatic_full_sweep.csv',index=False)
electro.head()
'''))

cells.append(code(r'''# Figure 01: potential and field for exponential and parabolic closures
from matplotlib.lines import Line2D
Vgroups=[[1.0,-1.0],[0.0,-5.0]]
fig,ax=plt.subplots(4,2,figsize=(176/25.4,190/25.4),sharex=True)
sam_styles=['-','--','-.',':']
for c,(model,source) in enumerate([('Exponential',profiles),('Parabolic',profiles_parabolic)]):
 for g,Vs in enumerate(Vgroups):
  for V in Vs:
   for j,S in enumerate(CFG["SAM_sweep_V"]):
    ph,E,_,_=source[(V,S)]; kw=dict(color=COL_V[V],ls=sam_styles[j],lw=1.05)
    ax[2*g,c].plot(z_nm,ph-ph[0],**kw); ax[2*g+1,c].plot(z_nm,E/1e5,**kw)
  ax[2*g,c].set_title(model+rf' potential, $V_{{app}}={Vs[0]:g}, {Vs[1]:g}$ V')
  ax[2*g+1,c].set_title(model+rf' field, $V_{{app}}={Vs[0]:g}, {Vs[1]:g}$ V')
  style(ax[2*g,c]); style(ax[2*g+1,c],g==1)
for g in range(2): ax[2*g,0].set_ylabel(r'$\phi(z)-\phi(0)$ (V)'); ax[2*g+1,0].set_ylabel(r'$E$ (kV cm$^{-1}$)')
letters(ax)
vhandles=[Line2D([0],[0],color=COL_V[V],lw=1.6,label=rf'{V:g} V') for V in CFG["Vapp_sweep_V"]]
shandles=[Line2D([0],[0],color='.15',ls=sam_styles[j],lw=1.6,label=rf'{S:+g} V') for j,S in enumerate(CFG["SAM_sweep_V"])]
fig.legend(handles=vhandles,title=r'Color: $V_{app}$',loc='center left',bbox_to_anchor=(.82,.68),frameon=False)
fig.legend(handles=shandles,title=r'Line style: $\Delta\phi_{SAM}$',loc='center left',bbox_to_anchor=(.82,.31),frameon=False)
fig.tight_layout(rect=(.035,.065,.81,.99)); fig.subplots_adjust(hspace=1.45,wspace=.34); save(fig,'01_potential_and_field'); plt.show()
'''))

cells.append(code(r'''# Figure 01b: one column per left-layer decay length; remaining sweeps overlaid
from matplotlib.lines import Line2D
Vgroups=[[1.0,-1.0],[0.0,-5.0]]
fig,ax=plt.subplots(4,len(CFG["lambda_sweep_nm"]),figsize=(190/25.4,190/25.4),sharex=True,squeeze=False)
sam_styles=['-','--','-.',':']
for i,L in enumerate(CFG["lambda_sweep_nm"]):
 for g,Vs in enumerate(Vgroups):
  for V in Vs:
   for j,S in enumerate(CFG["SAM_sweep_V"]):
    ph,E,_,_=profile(V,S,L); kw=dict(color=COL_V[V],ls=sam_styles[j],lw=1.05)
    ax[2*g,i].plot(z_nm,ph-ph[0],**kw); ax[2*g+1,i].plot(z_nm,E/1e5,**kw)
  style(ax[2*g,i]); style(ax[2*g+1,i],g==1)
  ax[2*g,i].set_title(rf'$\lambda_{{ZnO}}={L:g}$ nm; $V_{{app}}={Vs[0]:g},{Vs[1]:g}$ V')
for g in range(2): ax[2*g,0].set_ylabel(r'$\phi(z)-\phi(0)$ (V)'); ax[2*g+1,0].set_ylabel(r'$E$ (kV cm$^{-1}$)')
letters(ax)
vhandles=[Line2D([0],[0],color=COL_V[V],lw=1.6,label=rf'{V:g} V') for V in CFG["Vapp_sweep_V"]]
shandles=[Line2D([0],[0],color='.15',ls=sam_styles[j],lw=1.6,label=rf'{S:+g} V') for j,S in enumerate(CFG["SAM_sweep_V"])]
fig.legend(handles=vhandles,title=r'Color: $V_{app}$',loc='center left',bbox_to_anchor=(.82,.68),frameon=False)
fig.legend(handles=shandles,title=r'Line style: $\Delta\phi_{SAM}$',loc='center left',bbox_to_anchor=(.82,.31),frameon=False)
fig.tight_layout(rect=(.035,.065,.81,.99)); fig.subplots_adjust(hspace=1.45,wspace=.32); save(fig,'01b_voltage_lambda_sweeps'); plt.show()
'''))

cells.append(code(r'''# Figure 01c: one column per interlayer decay length; voltage and interlayer-step sweeps overlaid
from matplotlib.lines import Line2D
Vgroups=[[1.0,-1.0],[0.0,-5.0]]
fig,ax=plt.subplots(4,len(CFG["lambda_SAM_sweep_nm"]),figsize=(190/25.4,190/25.4),sharex=True,squeeze=False)
sam_styles=['-','--','-.',':']
for i,L in enumerate(CFG["lambda_SAM_sweep_nm"]):
 for g,Vs in enumerate(Vgroups):
  for V in Vs:
   for j,S in enumerate(CFG["SAM_sweep_V"]):
    ph,E,_,_=profile(V,S,CFG["lambda_ZnO_nm"],L); kw=dict(color=COL_V[V],ls=sam_styles[j],lw=1.05)
    ax[2*g,i].plot(z_nm,ph-ph[0],**kw); ax[2*g+1,i].plot(z_nm,E/1e5,**kw)
  style(ax[2*g,i]); style(ax[2*g+1,i],g==1)
  ax[2*g,i].set_title(rf'$\lambda_{{SAM}}={L:g}$ nm; $V_{{app}}={Vs[0]:g},{Vs[1]:g}$ V')
for g in range(2): ax[2*g,0].set_ylabel(r'$\phi(z)-\phi(0)$ (V)'); ax[2*g+1,0].set_ylabel(r'$E$ (kV cm$^{-1}$)')
letters(ax)
vhandles=[Line2D([0],[0],color=COL_V[V],lw=1.6,label=rf'{V:g} V') for V in CFG["Vapp_sweep_V"]]
shandles=[Line2D([0],[0],color='.15',ls=sam_styles[j],lw=1.6,label=rf'{S:+g} V') for j,S in enumerate(CFG["SAM_sweep_V"])]
fig.legend(handles=vhandles,title=r'Color: $V_{app}$',loc='center left',bbox_to_anchor=(.82,.68),frameon=False)
fig.legend(handles=shandles,title=r'Line style: $\Delta\phi_{SAM}$',loc='center left',bbox_to_anchor=(.82,.31),frameon=False)
fig.tight_layout(rect=(.035,.065,.81,.99)); fig.subplots_adjust(hspace=1.45,wspace=.32); save(fig,'01c_voltage_lambda_SAM_sweeps'); plt.show()
'''))

cells.append(code(r'''# Figure 02: field and derivatives for both electrostatic closures
fig,ax=plt.subplots(2,3,figsize=(190/25.4,110/25.4),sharex=True)
V=CFG["representative_Vapp_V"][0]
for r,(model,source) in enumerate([('Exponential',profiles),('Parabolic',profiles_parabolic)]):
    for S in CFG["derivative_SAM_V"]:
        _,E,G,H=source[(V,S)]
        col=COL_SAM[S]; ax[r,0].plot(z_nm,E/1e5,color=col); ax[r,1].plot(z_nm,G/1e15,color=col); ax[r,2].plot(z_nm,H/1e24,color=col)
    for c in range(3): style(ax[r,c],r==1); ax[r,c].axhline(0,color='.3',lw=.6)
    ax[r,0].set_ylabel(model+'\n'+r'$E$ (kV cm$^{-1}$)')
ax[0,0].set_title(r'$E$'); ax[0,1].set_title(r'$dE/dz$'); ax[0,2].set_title(r'$d^2E/dz^2$')
for r in range(2): ax[r,1].set_ylabel(r'$dE/dz$ ($10^{15}$ V m$^{-2}$)'); ax[r,2].set_ylabel(r'$d^2E/dz^2$ ($10^{24}$ V m$^{-3}$)')
letters(ax); fig.tight_layout(rect=(.045,.10,.99,.95)); fig.subplots_adjust(hspace=1.22,wspace=.72); save(fig,'02_field_first_second_derivatives'); plt.show()
'''))

cells.append(md(r'''## 3. Tensor-aligned anisotropic orientation calculation

For each position and electrostatic case, the notebook evaluates
\[U(\mathbf n,z)=-E\,\mathbf n\!\cdot\!\boldsymbol\mu-\tfrac12E^2\mathbf n^T\boldsymbol\alpha\mathbf n-\tfrac16E'\mathbf n^T\mathbf Q\mathbf n\]
on a user-controlled Fibonacci sphere. It reports the hard minimum and the finite-temperature Boltzmann average. The dipole, polarizability tensor, and optional quadrupole are always expressed in the same Cartesian frame.
'''))
cells.append(code(r'''def fibonacci(n):
    i=np.arange(n,dtype=float); golden=(1+5**.5)/2; y=1-(2*i+1)/n
    r=np.sqrt(np.maximum(0,1-y*y)); t=2*np.pi*i/golden
    return np.column_stack((r*np.cos(t),y,r*np.sin(t)))
nvec=fibonacci(int(CFG["orientation_samples"])); mup=nvec@mu
alp=np.einsum('ni,ij,nj->n',nvec,alpha,nvec); qproj=np.einsum('ni,ij,nj->n',nvec,Q,nvec)

hard_rows=[]; boltz_rows=[]; iso_rows=[]
for V in CFG["Vapp_sweep_V"]:
  for S in CFG["SAM_sweep_V"]:
    ph,E,G,H=profiles[(V,S)]
    for zz,p,e,g,h in zip(z_nm,ph,E,G,H):
      U=-e*mup-.5*e*e*alp-(1/6)*g*qproj
      j=int(np.argmin(U)); n=nvec[j]; me=mup[j]; ae=alp[j]; qe=qproj[j]
      Fmu=me*g; Fa=ae*e*g; Fq=(1/6)*qe*h; Ft=Fmu+Fa+Fq
      base=dict(Vapp_V=V,SAM_V=S,z_nm=zz,phi_V=p,E_V_m=e,dE_V_m2=g,d2E_V_m3=h)
      hr=base|dict(model='anisotropic_hard',nx=n[0],ny=n[1],nz=n[2],theta_deg=np.degrees(np.arccos(np.clip(n[2],-1,1))),phi_deg=np.degrees(np.arctan2(n[1],n[0])),mu_eff_D=me/D2CM,alpha_eff_A3=ae/(4*np.pi*eps0*1e-30),F_mu_N=Fmu,F_alpha_N=Fa,F_Q_N=Fq,F_total_N=Ft,U_J=U[j])
      us=(U-U.min())/(kB*CFG["orientation_temperature_K"]); w=np.exp(-np.clip(us,0,700)); w/=w.sum()
      nb=w@nvec; meb=w@mup; aeb=w@alp; qeb=w@qproj; Fmb=meb*g; Fab=aeb*e*g; Fqb=(1/6)*qeb*h
      br=base|dict(model='anisotropic_boltzmann',nx=nb[0],ny=nb[1],nz=nb[2],theta_deg=np.degrees(np.arccos(np.clip(nb[2]/max(np.linalg.norm(nb),1e-30),-1,1))),phi_deg=np.degrees(np.arctan2(nb[1],nb[0])),mu_eff_D=meb/D2CM,alpha_eff_A3=aeb/(4*np.pi*eps0*1e-30),F_mu_N=Fmb,F_alpha_N=Fab,F_Q_N=Fqb,F_total_N=Fmb+Fab+Fqb,U_J=np.nan)
      Fmi=mu_mag*g; Fai=alpha_iso*e*g; Fqi=(1/6)*Q_iso*h
      ir=base|dict(model='isotropic',nx=np.nan,ny=np.nan,nz=np.nan,theta_deg=np.nan,phi_deg=np.nan,mu_eff_D=mu_mag/D2CM,alpha_eff_A3=alpha_iso/(4*np.pi*eps0*1e-30),F_mu_N=Fmi,F_alpha_N=Fai,F_Q_N=Fqi,F_total_N=Fmi+Fai+Fqi,U_J=np.nan)
      for row in (hr,br,ir):
        for D in CFG["diffusion_m2_s"]:
          for tm in CFG["times_min"]:
            row[f'L_nm_D{D:.0e}_t{tm:g}min']=D*row['F_total_N']/(kB*CFG["temperature_K"])*tm*60*1e9
      hard_rows.append(hr); boltz_rows.append(br); iso_rows.append(ir)
hard=pd.DataFrame(hard_rows); boltz=pd.DataFrame(boltz_rows); iso=pd.DataFrame(iso_rows)
hard.to_csv(DAT/'anisotropic_hard_full_sweep.csv',index=False); boltz.to_csv(DAT/'anisotropic_boltzmann_full_sweep.csv',index=False); iso.to_csv(DAT/'isotropic_full_sweep.csv',index=False)
print(len(hard),'rows per model')
'''))

cells.append(code(r'''# Figure 03: anisotropic orientation and effective tensor projections
fig,ax=plt.subplots(4,2,figsize=(176/25.4,190/25.4),sharex=True)
spec=[('theta_deg',r'$\theta$ (deg)'),('phi_deg',r'$\varphi$ (deg)'),('mu_eff_D',r'$\mu_{eff}$ (D)'),('alpha_eff_A3',r'$\alpha_{eff}$ ($\AA^3$)')]
for c,V in enumerate(CFG["representative_Vapp_V"]):
  for r,(col,yl) in enumerate(spec):
    for S in CFG["SAM_sweep_V"]:
      g=hard[np.isclose(hard.Vapp_V,V)&np.isclose(hard.SAM_V,S)]; ax[r,c].plot(g.z_nm,g[col],color=COL_SAM[S])
    style(ax[r,c],r==3); ax[r,c].set_title(rf'$V_{{app}}={V:g}$ V');
    if c==0: ax[r,c].set_ylabel(yl)
letters(ax); fig.tight_layout(rect=(.035,.085,.82,.99)); fig.subplots_adjust(hspace=1.52,wspace=.34); save(fig,'03_anisotropic_orientation_tensor'); plt.show()
'''))

cells.append(code(r'''# Figure 04: isotropic force reference and drift
fig,ax=plt.subplots(2,2,figsize=(176/25.4,112/25.4),sharex=True)
for c,V in enumerate(CFG["representative_Vapp_V"]):
  S=CFG["derivative_SAM_V"][-1]; g=iso[np.isclose(iso.Vapp_V,V)&np.isclose(iso.SAM_V,S)]
  ax[0,c].plot(g.z_nm,g.F_mu_N*1e15,label='dipolar'); ax[0,c].plot(g.z_nm,g.F_alpha_N*1e15,label='polarizability'); ax[0,c].plot(g.z_nm,g.F_Q_N*1e15,label='quadrupolar'); ax[0,c].plot(g.z_nm,g.F_total_N*1e15,'k--',label='total')
  for D in CFG["diffusion_m2_s"]:
    col=f'L_nm_D{D:.0e}_t{CFG["times_min"][0]:g}min'; ax[1,c].plot(g.z_nm,g[col],label=rf'$D={D:.0e}$ m$^2$/s')
  for r in range(2): style(ax[r,c],r==1); ax[r,c].set_title(rf'Isotropic, $V_{{app}}={V:g}$ V')
ax[0,0].set_ylabel('Force (fN)'); ax[1,0].set_ylabel(r'$L_{drift}$ (nm)'); letters(ax)
fig.legend(*ax[0,0].get_legend_handles_labels(),loc='center left',bbox_to_anchor=(.84,.53),frameon=False)
fig.tight_layout(rect=(.035,.08,.82,.98)); fig.subplots_adjust(hspace=1.22,wspace=.34); save(fig,'04_isotropic_force_drift'); plt.show()
'''))

cells.append(code(r'''# Figure 05: anisotropic hard/Boltzmann force and drift comparison
fig,ax=plt.subplots(2,2,figsize=(176/25.4,112/25.4),sharex=True)
for c,V in enumerate(CFG["representative_Vapp_V"]):
  S=CFG["derivative_SAM_V"][-1]
  for df,col,lab in [(hard,'#0072B2','hard minimum'),(boltz,'#009E73','Boltzmann')]:
    g=df[np.isclose(df.Vapp_V,V)&np.isclose(df.SAM_V,S)]; ax[0,c].plot(g.z_nm,g.F_total_N*1e15,color=col,label=lab)
    D=CFG["diffusion_m2_s"][1]; tm=CFG["times_min"][0]; ax[1,c].plot(g.z_nm,g[f'L_nm_D{D:.0e}_t{tm:g}min'],color=col)
  for r in range(2): style(ax[r,c],r==1); ax[r,c].axhline(0,color='.3',lw=.6); ax[r,c].set_title(rf'Anisotropic, $V_{{app}}={V:g}$ V')
ax[0,0].set_ylabel('Total force (fN)'); ax[1,0].set_ylabel(r'$L_{drift}$ (nm)'); letters(ax)
fig.legend(*ax[0,0].get_legend_handles_labels(),loc='center left',bbox_to_anchor=(.84,.53),frameon=False)
fig.tight_layout(rect=(.035,.08,.82,.98)); fig.subplots_adjust(hspace=1.22,wspace=.34); save(fig,'05_anisotropic_force_drift'); plt.show()
'''))

cells.append(code(r'''# Figure 06: full user-controlled Ldrift atlas, all D values on each plot
for tm in CFG["times_min"]:
  fig,ax=plt.subplots(len(CFG["SAM_sweep_V"]),len(CFG["representative_Vapp_V"]),figsize=(176/25.4,190/25.4),sharex=True,squeeze=False)
  for r,S in enumerate(CFG["SAM_sweep_V"]):
    for c,V in enumerate(CFG["representative_Vapp_V"]):
      a=ax[r,c]
      for df,color,lab in [(iso,'#0072B2','isotropic')]:
        g=df[np.isclose(df.Vapp_V,V)&np.isclose(df.SAM_V,S)]
        for i,D in enumerate(CFG["diffusion_m2_s"]): a.plot(g.z_nm,g[f'L_nm_D{D:.0e}_t{tm:g}min'],color=color,ls=['-','--','-.',':'][i%4],label=lab if i==0 else None)
      style(a,r==len(CFG["SAM_sweep_V"])-1); a.axhline(0,color='.3',lw=.6); a.set_title(rf'$V_{{app}}={V:g}$ V, $\Delta\phi_{{SAM}}={S:+g}$ V')
      if c==0:a.set_ylabel(rf'$L_{{drift}}$ (nm), {tm:g} min')
  letters(ax); fig.tight_layout(rect=(.035,.085,.82,.99)); fig.subplots_adjust(hspace=1.52,wspace=.34); save(fig,f'06_Ldrift_atlas_{tm:g}min'); plt.show()
'''))

cells.append(code(r'''# Figure 07: anisotropic hard-minimum force-component atlas
fig,ax=plt.subplots(len(CFG["SAM_sweep_V"]),len(CFG["representative_Vapp_V"]),figsize=(176/25.4,190/25.4),sharex=True,squeeze=False)
for r,S in enumerate(CFG["SAM_sweep_V"]):
 for c,V in enumerate(CFG["representative_Vapp_V"]):
  a=ax[r,c];g=hard[np.isclose(hard.Vapp_V,V)&np.isclose(hard.SAM_V,S)]
  a.plot(g.z_nm,g.F_mu_N*1e15,color='#0072B2',label='dipolar');a.plot(g.z_nm,g.F_alpha_N*1e15,color='#009E73',label='polarizability');a.plot(g.z_nm,g.F_Q_N*1e15,color='#D55E00',label='quadrupolar');a.plot(g.z_nm,g.F_total_N*1e15,color='.1',ls='--',label='total')
  a.axhline(0,color='.3',lw=.6);style(a,r==len(CFG["SAM_sweep_V"])-1);a.set_title(rf'$V_{{app}}={V:g}$ V, $\Delta\phi_{{SAM}}={S:+g}$ V')
  if c==0:a.set_ylabel('Force (fN)')
letters(ax);fig.legend(*ax[0,0].get_legend_handles_labels(),loc='center left',bbox_to_anchor=(.84,.53),frameon=False,title='Anisotropic force');fig.tight_layout(rect=(.035,.085,.82,.99));fig.subplots_adjust(hspace=1.52,wspace=.34);save(fig,'07_anisotropic_force_components');plt.show()
'''))

cells.append(code(r'''# Figure 08: isotropic versus hard-minimum versus Boltzmann drift
D=max(CFG["diffusion_m2_s"]);tm=max(CFG["times_min"])
fig,ax=plt.subplots(len(CFG["SAM_sweep_V"]),len(CFG["representative_Vapp_V"]),figsize=(176/25.4,190/25.4),sharex=True,squeeze=False)
for r,S in enumerate(CFG["SAM_sweep_V"]):
 for c,V in enumerate(CFG["representative_Vapp_V"]):
  a=ax[r,c]
  for df,color,lab in [(iso,'#D55E00','isotropic'),(hard,'#0072B2','hard minimum'),(boltz,'#009E73','Boltzmann')]:
   g=df[np.isclose(df.Vapp_V,V)&np.isclose(df.SAM_V,S)];a.plot(g.z_nm,g[f'L_nm_D{D:.0e}_t{tm:g}min'],color=color,label=lab)
  a.axhline(0,color='.3',lw=.6);style(a,r==len(CFG["SAM_sweep_V"])-1);a.set_title(rf'$V_{{app}}={V:g}$ V, $\Delta\phi_{{SAM}}={S:+g}$ V')
  if c==0:a.set_ylabel(r'$L_{drift}$ (nm)')
letters(ax);fig.legend(*ax[0,0].get_legend_handles_labels(),loc='center left',bbox_to_anchor=(.84,.53),frameon=False,title=rf'$D={D:.0e}$ m$^2$/s; $t={tm:g}$ min');fig.tight_layout(rect=(.035,.085,.82,.99));fig.subplots_adjust(hspace=1.52,wspace=.34);save(fig,'08_isotropic_anisotropic_comparison');plt.show()
'''))

cells.append(code(r'''# Figure 09: energy-minimization and orientational-closure diagnostics
fig,ax=plt.subplots(3,2,figsize=(176/25.4,145/25.4),sharex=True)
for c,V in enumerate(CFG["representative_Vapp_V"]):
 for S in CFG["SAM_sweep_V"]:
  gh=hard[np.isclose(hard.Vapp_V,V)&np.isclose(hard.SAM_V,S)];gb=boltz[np.isclose(boltz.Vapp_V,V)&np.isclose(boltz.SAM_V,S)]
  ax[0,c].plot(gh.z_nm,gh.U_J/1e-21,color=COL_SAM[S]);ax[1,c].plot(gh.z_nm,gh.theta_deg,color=COL_SAM[S]);ax[2,c].plot(gb.z_nm,np.sqrt(gb.nx**2+gb.ny**2+gb.nz**2),color=COL_SAM[S])
 for r in range(3):style(ax[r,c],r==2);ax[r,c].set_title(rf'$V_{{app}}={V:g}$ V')
ax[0,0].set_ylabel(r'$U_{min}$ ($10^{-21}$ J)');ax[1,0].set_ylabel(r'$\theta_{min}$ (deg)');ax[2,0].set_ylabel(r'$|\langle n\rangle|$');letters(ax)
fig.tight_layout(rect=(.035,.08,.99,.99));fig.subplots_adjust(hspace=1.34,wspace=.45);save(fig,'09_orientation_minimization_diagnostics');plt.show()
'''))

cells.append(md(r'''## 7. Paper-sequence diagnostics and fixed-case calculations

The paper figures below use the fixed diagnostic case $V_{app}=-5$ V and $\Delta\phi_{SAM}=-2$ V where requested. Isotropic results are shown before orientation minimization. Every anisotropic result uses the exponential electrostatic closure only.
'''))

cells.append(code(r'''# Figure 10: exponential/parabolic potential and field components, fixed case
V0,S0=-5.0,-2.0
fig,ax=plt.subplots(2,2,figsize=(176/25.4,112/25.4),sharex=True)
for c,(model,fn) in enumerate([('Exponential',profile),('Parabolic',profile_parabolic)]):
 ph,E,G,H=fn(V0,S0); base=-(CFG['V_bi_V']-V0)*z/d; Eb=np.full_like(z,(CFG['V_bi_V']-V0)/d)
 if model=='Exponential':
  pz=CFG['delta_phi_ZnO_V']*np.exp(-z/(CFG['lambda_ZnO_nm']*NM)); ps=S0*np.exp(-z/(CFG['lambda_SAM_nm']*NM))
  Ez=CFG['delta_phi_ZnO_V']/(CFG['lambda_ZnO_nm']*NM)*np.exp(-z/(CFG['lambda_ZnO_nm']*NM)); Es=S0/(CFG['lambda_SAM_nm']*NM)*np.exp(-z/(CFG['lambda_SAM_nm']*NM))
 else:
  pz,Ez,_,_=profile_parabolic(CFG['V_bi_V'],0); ps,Es,_,_=profile_parabolic(CFG['V_bi_V'],S0); ps-=pz; Es-=Ez
 for y,lab,col,ls in [(base,'background','#202020','--'),(pz,f"{CFG['left_layer_name']} step",'#0072B2','-'),(ps,f"{CFG['interlayer_name']} step",'#D55E00','-'),(ph,'total','#009E73','-')]: ax[0,c].plot(z_nm,y-y[0],color=col,ls=ls,label=lab)
 for y,lab,col,ls in [(Eb,'background','#202020','--'),(Ez,f"{CFG['left_layer_name']} step",'#0072B2','-'),(Es,f"{CFG['interlayer_name']} step",'#D55E00','-'),(E,'total','#009E73','-')]: ax[1,c].plot(z_nm,y/1e5,color=col,ls=ls,label=lab)
 ax[0,c].set_title(model+' potential components'); ax[1,c].set_title(model+' field components'); style(ax[0,c]); style(ax[1,c],True)
ax[0,0].set_ylabel(r'$\Delta\phi$ (V)'); ax[1,0].set_ylabel(r'$E$ (kV cm$^{-1}$)'); letters(ax)
fig.legend(*ax[0,0].get_legend_handles_labels(),loc='center left',bbox_to_anchor=(.82,.52),frameon=False);fig.tight_layout(rect=(.035,.08,.81,.98));fig.subplots_adjust(hspace=1.22,wspace=.34);save(fig,'10_fixed_case_model_components');plt.show()
'''))

cells.append(code(r'''# Figure 11: isotropic analytical derivatives, exponential fixed case
ph,E,G,H=profile(-5.0,-2.0); fig,ax=plt.subplots(2,2,figsize=(176/25.4,112/25.4),sharex=True)
for i,(a,y,title,yl) in enumerate(zip(ax.ravel(),[ph-ph[0],E/1e5,G/1e15,H/1e24],[r'$\phi(z)-\phi(0)$',r'$E$',r'$dE/dz$',r'$d^2E/dz^2$'],['V',r'kV cm$^{-1}$',r'$10^{15}$ V m$^{-2}$',r'$10^{24}$ V m$^{-3}$'])): a.plot(z_nm,y,color='#0072B2');a.set_title(title);a.set_ylabel(yl);style(a,i>=2)
letters(ax);fig.tight_layout(rect=(.035,.08,.99,.98));fig.subplots_adjust(hspace=1.22,wspace=.36);save(fig,'11_isotropic_analytical_derivatives');plt.show()
'''))

cells.append(code(r'''# Figure 12: isotropic force components and total, exponential fixed case
g=iso[np.isclose(iso.Vapp_V,-5)&np.isclose(iso.SAM_V,-2)];fig,ax=plt.subplots(figsize=(65/25.4,55/25.4))
for col,lab,color,ls in [('F_mu_N','dipolar','#0072B2','-'),('F_alpha_N','polarizability','#009E73','-'),('F_Q_N','quadrupolar','#D55E00','-'),('F_total_N','total','#202020','--')]:ax.plot(g.z_nm,g[col]*1e15,label=lab,color=color,ls=ls)
ax.axhline(0,color='.3',lw=.6);ax.set_ylabel('Force (fN)');style(ax,True);fig.text(.5,.96,'blue dipole | green polarizability | orange quadrupole | black dashed total',ha='center',va='top',fontsize=3.8);fig.tight_layout(rect=(.08,.16,.99,.82));save(fig,'12_isotropic_fixed_force_components');plt.show()
'''))

cells.append(code(r'''# Figure 13: isotropic total-force sweep, six lines
fig,ax=plt.subplots(figsize=(65/25.4,55/25.4)); styles=['-','--','-.']
for V in [-1.0,-5.0]:
 for j,S in enumerate([-2.0,-.7,.7]):
  g=iso[np.isclose(iso.Vapp_V,V)&np.isclose(iso.SAM_V,S)];ax.plot(g.z_nm,g.F_total_N*1e15,color=COL_V[V],ls=styles[j],label=rf'$V_{{app}}={V:g}$ V, $\Delta\phi_{{SAM}}={S:+g}$ V')
ax.axhline(0,color='.3',lw=.6);ax.set_ylabel('Total force (fN)');style(ax,True)
fig.text(.5,.97,'green: V=-1 V | blue: V=-5 V',ha='center',va='top',fontsize=4.2);fig.text(.5,.90,'solid: SAM=-2 V | dashed: -0.7 V | dash-dot: +0.7 V',ha='center',va='top',fontsize=3.8);fig.tight_layout(rect=(.08,.16,.99,.78));save(fig,'13_isotropic_total_force_sweep');plt.show()
'''))

cells.append(code(r'''# Figures 14-16: three distinct exponential orientation minima, forces, and drift
V0,S0=-5.0,-2.0; ph,E,G,H=profile(V0,S0); states=[]
for zz,e,g,h in zip(z_nm,E,G,H):
 U=-e*mup-.5*e*e*alp-(1/6)*g*qproj; order=np.argsort(U); chosen=[]
 for idx in order:
  if all(np.degrees(np.arccos(np.clip(abs(nvec[idx]@nvec[k]),-1,1)))>12 for k in chosen): chosen.append(int(idx))
  if len(chosen)==3: break
 for rank,idx in enumerate(chosen,1):
  n=nvec[idx]; me=mup[idx]; ae=alp[idx]; qe=qproj[idx]; Fm=me*g; Fa=ae*e*g; Fq=(1/6)*qe*h
  states.append(dict(z_nm=zz,state=rank,U_J=U[idx],theta_deg=np.degrees(np.arccos(np.clip(n[2],-1,1))),phi_deg=np.degrees(np.arctan2(n[1],n[0])),nx=n[0],ny=n[1],nz=n[2],F_mu_N=Fm,F_alpha_N=Fa,F_Q_N=Fq,F_total_N=Fm+Fa+Fq))
states=pd.DataFrame(states);states.to_csv(DAT/'three_orientation_minima_exponential.csv',index=False)
cols=['#0072B2','#D55E00','#009E73'];fig,ax=plt.subplots(3,1,figsize=(65/25.4,120/25.4),sharex=True)
for s,col in zip([1,2,3],cols):q3=states[states.state==s];ax[0].plot(q3.z_nm,q3.U_J/1e-21,color=col,label=f'minimum {s}');ax[1].plot(q3.z_nm,q3.theta_deg,color=col);ax[2].plot(q3.z_nm,q3.phi_deg,color=col)
for i,a in enumerate(ax):style(a,i==2);a.axhline(0,color='.3',lw=.5)
ax[0].set_ylabel(r'$U$ ($10^{-21}$ J)');ax[1].set_ylabel(r'$\theta$ (deg)');ax[2].set_ylabel(r'$\varphi$ (deg)');ax[0].legend(frameon=False);letters(ax);fig.tight_layout(rect=(.05,.08,.99,.98));fig.subplots_adjust(hspace=1.18);save(fig,'14_three_orientation_energy_minima');plt.show()
fig,ax=plt.subplots(3,1,figsize=(65/25.4,125/25.4),sharex=True)
for s,col in zip([1,2,3],cols):
 q3=states[states.state==s]
 for key,lab,ls in [('F_mu_N','dipolar','-'),('F_alpha_N','polarizability','--'),('F_Q_N','quadrupolar','-.')]:ax[s-1].plot(q3.z_nm,q3[key]*1e15,color=col,ls=ls,label=lab)
 ax[s-1].plot(q3.z_nm,q3.F_total_N*1e15,color='.15',lw=1.3,label='total');ax[s-1].set_title(f'Orientation minimum {s}');ax[s-1].set_ylabel('Force (fN)');style(ax[s-1],s==3)
fig.text(.5,.99,'solid dipolar | dashed polarizability | dash-dot quadrupolar | black total',ha='center',va='top',fontsize=3.8);letters(ax);fig.tight_layout(rect=(.05,.08,.99,.95));fig.subplots_adjust(hspace=1.2);save(fig,'15_three_orientation_force_components');plt.show()
D=max(CFG['diffusion_m2_s']);tm=max(CFG['times_min']);fig,ax=plt.subplots(figsize=(65/25.4,55/25.4))
for s,col in zip([1,2,3],cols):q3=states[states.state==s];L=D*q3.F_total_N/(kB*CFG['temperature_K'])*(tm*60)*1e9;ax.plot(q3.z_nm,L,color=col,label=f'orientation {s}')
ax.axhline(0,color='.3',lw=.6);ax.set_ylabel(r'$L_{drift}$ (nm)');style(ax,True);fig.text(.5,.96,'blue orientation 1 | orange orientation 2 | green orientation 3',ha='center',va='top',fontsize=4);fig.tight_layout(rect=(.08,.16,.99,.82));save(fig,'16_three_orientation_Ldrift');plt.show()
'''))

cells.append(code(r'''# Figure 17: exponential anisotropic hard-minimum Ldrift, full sweeps, <=6 lines/panel
Vgroups=[[1.0,-1.0],[0.0,-5.0]]; dstyles=['-','--','-.',':']
for tm in CFG['times_min']:
 fig,ax=plt.subplots(len(CFG['SAM_sweep_V']),2,figsize=(176/25.4,190/25.4),sharex=True,squeeze=False)
 for r,S in enumerate(CFG['SAM_sweep_V']):
  for c,Vs in enumerate(Vgroups):
   a=ax[r,c]
   for V in Vs:
    g=hard[np.isclose(hard.Vapp_V,V)&np.isclose(hard.SAM_V,S)]
    for j,D in enumerate(CFG['diffusion_m2_s']): a.plot(g.z_nm,g[f'L_nm_D{D:.0e}_t{tm:g}min'],color=COL_V[V],ls=dstyles[j],label=rf'$V={V:g}$, $D={D:.0e}$')
   a.axhline(0,color='.3',lw=.6);style(a,r==len(CFG['SAM_sweep_V'])-1);a.set_title(rf'$\Delta\phi_{{SAM}}={S:+g}$ V; $V_{{app}}={Vs[0]:g},{Vs[1]:g}$ V')
   if c==0:a.set_ylabel(r'$L_{drift}$ (nm)')
 letters(ax);fig.legend(*ax[0,0].get_legend_handles_labels(),loc='center left',bbox_to_anchor=(.82,.52),frameon=False,title=rf'Exponential anisotropic; {tm:g} min');fig.tight_layout(rect=(.035,.065,.81,.99));fig.subplots_adjust(hspace=1.5,wspace=.34);save(fig,f'17_anisotropic_Ldrift_full_sweep_{tm:g}min');plt.show()
'''))

cells.append(code(r'''# Final numerical summary and reproducibility manifest
summary=[]
for name,df in [('isotropic',iso),('anisotropic_hard',hard),('anisotropic_boltzmann',boltz)]:
  for (V,S),g in df.groupby(['Vapp_V','SAM_V']):
    row={'model':name,'Vapp_V':V,'SAM_V':S}
    for D in CFG["diffusion_m2_s"]:
      for tm in CFG["times_min"]:
        col=f'L_nm_D{D:.0e}_t{tm:g}min'; a=g[col].to_numpy(); row[f'max_abs_{col}']=float(np.max(np.abs(a))); row[f'rms_{col}']=float(np.sqrt(np.mean(a*a)))
    summary.append(row)
summary=pd.DataFrame(summary); summary.to_csv(DAT/'final_drift_summary.csv',index=False)
manifest={k:(v.tolist() if isinstance(v,np.ndarray) else str(v) if isinstance(v,Path) else v) for k,v in CFG.items()}
(OUT/'user_inputs.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('Complete. Outputs:',OUT.resolve()); display(summary.head())
'''))

# Embed every static/non-numerical manuscript asset so the notebook remains a
# single-file archive. Running this final cell restores the exact workflow,
# band-bending, molecular, PyMOL, and legacy comparison panels used by the DOCX.
asset_names = ['Figure_01_workflow_or_cover.png','Figure_02_band_bending_scheme.png','Figure_09a_PyMOL_unrotated.png','Figure_09b_PyMOL_before.png','Figure_09c_PyMOL_after.png']
embedded={name:base64.b64encode((ROOT/'assets'/name).read_bytes()).decode('ascii') for name in asset_names}
payload=json.dumps(embedded,separators=(',',':'))
cells.append(md('''## 8. Static schematics and molecular render views\n\nThe workflow and band-bending scheme are embedded for reproducible restoration. The packaged molecular views belong specifically to the default Y6 example; for another molecule they are not used or relabeled. Every numerical curve is regenerated from the current GUI inputs by the calculation cells above.\n'''))
cells.append(code("""import base64\nMANUSCRIPT_ASSETS = json.loads(r'''"""+payload+"""''')\nfor filename,encoded in MANUSCRIPT_ASSETS.items():\n    (FIG/filename).write_bytes(base64.b64decode(encoded))\nprint(f'Restored {len(MANUSCRIPT_ASSETS)} exact manuscript assets to {FIG}')\n"""))

cells.append(code(r'''# Figure 14b: packaged views are valid only for the default Y6 geometry
if CFG['molecule_name'].strip().casefold() == 'y6':
 views=[('Figure_09a_PyMOL_unrotated.png','Unrotated'),('Figure_09b_PyMOL_before.png','Minimum orientation 1'),('Figure_09c_PyMOL_after.png','Minimum orientation 2')]
 fig,ax=plt.subplots(1,3,figsize=(190/25.4,68/25.4))
 for i,(name,title) in enumerate(views):
  ax[i].imshow(plt.imread(FIG/name));ax[i].set_title(f"{CFG['molecule_name']}: {title}");ax[i].axis('off')
 letters(ax);fig.tight_layout(rect=(.01,.01,.99,.98));save(fig,'14b_positioned_Y6_views');plt.show()
else:
 print(f"Skipping packaged Y6 views for custom molecule {CFG['molecule_name']!r}. Supply molecule-specific rendered views before report generation.")
'''))

nb={"cells":cells,"metadata":{"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"},"language_info":{"name":"python","version":"3"}},"nbformat":4,"nbformat_minor":5}
OUT.write_text(json.dumps(nb,indent=1),encoding='utf-8')
print(OUT)
