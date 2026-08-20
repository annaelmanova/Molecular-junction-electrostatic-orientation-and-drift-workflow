"""GUI front end for Molecular_Junction_Workflow.ipynb."""
from __future__ import annotations
import ast, json, os, re, subprocess, sys, threading, traceback, webbrowser
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

ROOT=Path(__file__).resolve().parent
TEMPLATE=ROOT/'Molecular_Junction_Workflow.ipynb'
OUTROOT=ROOT/'Molecular_Junction_GUI_runs'

DEFAULTS={
 'molecule_name':'Y6','left_layer_name':'ZnO','interlayer_name':'SAM','active_layer_name':'PM6:Y6',
 'active_layer_nm':'100.0','plot_zmax_nm':'50.0','z_step_nm':'0.10','V_bi_V':'0.85',
 'delta_phi_ZnO_V':'0.30','lambda_ZnO_nm':'30.0','lambda_SAM_nm':'30.0',
 'Vapp_sweep_V':'1, 0, -1, -5','SAM_sweep_V':'-2, -0.7, 0.7, 2',
 'lambda_sweep_nm':'10, 20, 30','lambda_SAM_sweep_nm':'10, 20, 30','temperature_K':'418.15','orientation_temperature_K':'418.15',
 'diffusion_m2_s':'1e-19, 1e-18, 1e-17','times_min':'10, 60','orientation_samples':'8192',
 'mu_D':'11.4904, -0.2662, -0.9078',
 'alpha_A3':'331.972, -24.940, -27.115; -24.940, 199.188, 6.388; -27.115, 6.388, 106.596',
 'Q_DA':'0,0,0; 0,0,0; 0,0,0','representative_Vapp_V':'0, -5','derivative_SAM_V':'-0.7, 0.7'}

HELP="""Molecular junction tensor-workflow GUI - complete instructions

1. QC INPUT
Paste a Gaussian, ORCA, Q-Chem, or plain-text output containing a Cartesian dipole and polarizability tensor. Click Parse QC text. The parser recognizes common dipole-component lines and 3x3 tensor blocks. Always inspect the parsed values; program formats differ. The dipole, polarizability tensor, quadrupole, and XYZ geometry must share the same Cartesian frame.

If automatic parsing is not possible, enter the three dipole components and the symmetric 3x3 polarizability tensor manually. Dipole units are Debye; polarizability units are cubic angstrom. The quadrupole is optional and defaults to zero because mixing tensors from different calculations or frames is invalid.

2. DEVICE AND SWEEPS
The molecule and all three layer labels are user editable. Defaults reproduce Y6/ZnO/SAM/PM6:Y6. Custom names propagate to notebook layer strips and reports; packaged molecular renders remain Y6-specific and are skipped for another molecule.

Comma-separated fields define sweeps. Applied voltage and the interlayer potential step accept signed values. Diffusion coefficients are entered in square metres per second (1e-17 m2/s equals 1e-13 cm2/s). Times are minutes. The active-layer thickness sets the uniform background field; plot_zmax controls only the displayed/calculated interfacial range.

3. NUMERICS
z_step_nm controls spatial resolution. orientation_samples controls the Fibonacci-sphere orientation grid. 8192 is the recommended normal setting; 2048 is useful for previews and 16384 for convergence checks. Runtime scales approximately linearly with both orientation_samples and the number of z points and electrostatic cases.

4. VALIDATE AND CREATE
Click Validate all inputs. Errors must be fixed before notebook generation. Create configured notebook writes a new self-contained IPYNB and a JSON manifest into a timestamped run folder. The original V15 template is never overwritten.

5. RUN
Run calculation executes the configured notebook code headlessly with the available Codex plotting environment. In Jupyter, open the generated notebook and use Run All to obtain inline outputs. PNG, PDF, CSV, and JSON files are written below the selected run folder.

6. PHYSICAL INTERPRETATION
Potential determines orientation through U(n,z). Translation requires an inhomogeneous field: the permanent-dipole force scales with dE/dz, the induced-dipole force with E*dE/dz, and the optional quadrupole force with d2E/dz2. Hard-minimum branch exchange near E=0 is a degeneracy; the finite-temperature Boltzmann result is the smooth equilibrium closure. Drift lengths are local Einstein-relation estimates, not full morphology trajectories.

7. REPRODUCIBILITY
Archive the configured notebook, user_inputs.json, original QC output, XYZ geometry, and generated CSV tables together. Do not combine a dipole/tensor from one coordinate frame with a differently oriented geometry.
"""

def nums(text): return [float(x.strip()) for x in re.split(r'[,\s]+',text.strip()) if x.strip()]
def matrix(text):
 rows=[nums(r) for r in text.strip().split(';') if r.strip()]
 if len(rows)!=3 or any(len(r)!=3 for r in rows): raise ValueError('matrix must contain three semicolon-separated rows of three numbers')
 return rows

def parse_qc(text):
 out={}
 # Common x/y/z dipole component formats.
 pats=[r'X\s*=\s*([-+\d.Ee]+)\s+Y\s*=\s*([-+\d.Ee]+)\s+Z\s*=\s*([-+\d.Ee]+)',
       r'Dipole[^\n]*\n[^\n]*?([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)']
 for p in pats:
  m=re.search(p,text,re.I)
  if m: out['mu_D']=[float(m.group(i)) for i in range(1,4)]; break
 # Explicit GUI-friendly block or common "polarizability tensor" block.
 m=re.search(r'(?:polarizability(?:\s+tensor)?|alpha)\s*[:=]?\s*\n((?:[^\n]*\n){3})',text,re.I)
 if m:
  rows=[]
  for line in m.group(1).splitlines():
   v=re.findall(r'[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?',line)
   if len(v)>=3: rows.append([float(x) for x in v[-3:]])
  if len(rows)==3: out['alpha_A3']=rows
 return out

class App(tk.Tk):
 def __init__(self):
  super().__init__(); self.title('Molecular junction tensor workflow'); self.geometry('1040x760'); self.minsize(900,650)
  self.vars={k:tk.StringVar(value=v) for k,v in DEFAULTS.items()}; self.run_dir=None
  self._style(); self._build()
 def _style(self):
  s=ttk.Style(self); s.theme_use('clam'); s.configure('TButton',padding=7); s.configure('TNotebook.Tab',padding=(14,8)); s.configure('Header.TLabel',font=('Segoe UI',15,'bold'),foreground='#17365D')
 def _build(self):
  top=ttk.Frame(self,padding=12); top.pack(fill='x'); ttk.Label(top,text='Molecular interfacial electrostatics and tensor drift',style='Header.TLabel').pack(side='left')
  nb=ttk.Notebook(self); nb.pack(fill='both',expand=True,padx=12,pady=(0,8)); self.nb=nb
  self._form_tab(nb,'Device & sweeps',[('molecule_name','Molecule name'),('left_layer_name','Left layer name'),('interlayer_name','Interlayer name'),('active_layer_name','Active layer name'),('active_layer_nm','Active layer thickness (nm)'),('plot_zmax_nm','Calculation/plot range (nm)'),('z_step_nm','Spatial step (nm)'),('V_bi_V','Built-in voltage (V)'),('delta_phi_ZnO_V','Left-layer potential step (V)'),('lambda_ZnO_nm','Left-layer decay length (nm)'),('lambda_SAM_nm','Interlayer decay length (nm)'),('Vapp_sweep_V','Applied-voltage sweep (V)'),('SAM_sweep_V','Interlayer-step sweep (V)'),('lambda_sweep_nm','Left-layer decay-length sweep (nm)'),('lambda_SAM_sweep_nm','Interlayer decay-length sweep (nm)')])
  self._form_tab(nb,'Transport & numerics',[('temperature_K','Transport temperature (K)'),('orientation_temperature_K','Orientation temperature (K)'),('diffusion_m2_s','Diffusion coefficients (m²/s)'),('times_min','Drift times (min)'),('orientation_samples','Orientation samples'),('representative_Vapp_V','Representative voltages (V)'),('derivative_SAM_V','Derivative-panel SAM steps (V)')])
  self._qc_tab(nb); self._instructions(nb); self._run_tab(nb)
 def _form_tab(self,nb,title,fields):
  f=ttk.Frame(nb,padding=18); nb.add(f,text=title)
  for i,(key,label) in enumerate(fields):
   ttk.Label(f,text=label).grid(row=i,column=0,sticky='w',pady=6,padx=(0,18)); ttk.Entry(f,textvariable=self.vars[key],width=55).grid(row=i,column=1,sticky='ew',pady=6)
  f.columnconfigure(1,weight=1)
 def _qc_tab(self,nb):
  f=ttk.Frame(nb,padding=14); nb.add(f,text='QC input')
  bar=ttk.Frame(f); bar.pack(fill='x'); ttk.Button(bar,text='Load QC file',command=self.load_qc).pack(side='left'); ttk.Button(bar,text='Parse QC text',command=self.do_parse).pack(side='left',padx=8)
  self.qc=tk.Text(f,height=15,wrap='none',font=('Consolas',9)); self.qc.pack(fill='both',expand=True,pady=8)
  grid=ttk.Frame(f); grid.pack(fill='x')
  for i,(k,l) in enumerate([('mu_D','Dipole x,y,z (D)'),('alpha_A3','Polarizability rows (Å³)'),('Q_DA','Quadrupole rows (D Å)')]):
   ttk.Label(grid,text=l).grid(row=i,column=0,sticky='w',pady=5); ttk.Entry(grid,textvariable=self.vars[k],width=85).grid(row=i,column=1,sticky='ew',padx=10,pady=5)
  grid.columnconfigure(1,weight=1)
 def _instructions(self,nb):
  f=ttk.Frame(nb,padding=12); nb.add(f,text='Instructions'); t=tk.Text(f,wrap='word',font=('Segoe UI',10)); t.insert('1.0',HELP); t.configure(state='disabled'); t.pack(fill='both',expand=True)
 def _run_tab(self,nb):
  f=ttk.Frame(nb,padding=14); nb.add(f,text='Validate & run'); b=ttk.Frame(f); b.pack(fill='x')
  ttk.Button(b,text='Validate all inputs',command=self.validate_dialog).pack(side='left'); ttk.Button(b,text='Create configured notebook',command=self.create_notebook).pack(side='left',padx=8); ttk.Button(b,text='Run calculation',command=self.run).pack(side='left'); ttk.Button(b,text='Open output folder',command=self.open_folder).pack(side='left',padx=8)
  self.log=tk.Text(f,wrap='word',font=('Consolas',9)); self.log.pack(fill='both',expand=True,pady=12)
 def load_qc(self):
  p=filedialog.askopenfilename(title='Select QC text output',filetypes=[('QC/text files','*.out *.log *.txt'),('All files','*.*')])
  if p: self.qc.delete('1.0','end'); self.qc.insert('1.0',Path(p).read_text(errors='replace'))
 def do_parse(self):
  got=parse_qc(self.qc.get('1.0','end'))
  if 'mu_D' in got:self.vars['mu_D'].set(', '.join(map(str,got['mu_D'])))
  if 'alpha_A3' in got:self.vars['alpha_A3'].set('; '.join(', '.join(map(str,r)) for r in got['alpha_A3']))
  messagebox.showinfo('QC parser',f"Parsed: {', '.join(got) if got else 'nothing recognized; use manual fields'}")
 def values(self):
  v={'save_png':True,'save_pdf':True,'dpi':300,'boltzmann_chunk':128}; scalar=['active_layer_nm','plot_zmax_nm','z_step_nm','V_bi_V','delta_phi_ZnO_V','lambda_ZnO_nm','lambda_SAM_nm','temperature_K','orientation_temperature_K']
  for k in ['molecule_name','left_layer_name','interlayer_name','active_layer_name']:
   v[k]=self.vars[k].get().strip()
   if not v[k]:raise ValueError(f'{k.replace("_"," ")} cannot be empty')
  for k in scalar:v[k]=float(self.vars[k].get())
  for k in ['Vapp_sweep_V','SAM_sweep_V','lambda_sweep_nm','lambda_SAM_sweep_nm','diffusion_m2_s','times_min','representative_Vapp_V','derivative_SAM_V']:v[k]=nums(self.vars[k].get())
  v['orientation_samples']=int(self.vars['orientation_samples'].get()); v['mu_D']=nums(self.vars['mu_D'].get()); v['alpha_A3']=matrix(self.vars['alpha_A3'].get()); v['Q_DA']=matrix(self.vars['Q_DA'].get())
  if len(v['mu_D'])!=3:raise ValueError('dipole must have exactly 3 components')
  if v['z_step_nm']<=0 or v['plot_zmax_nm']<=0 or v['active_layer_nm']<v['plot_zmax_nm']:raise ValueError('Require 0 < z step, 0 < plot range <= active-layer thickness')
  if any(x<=0 for x in v['diffusion_m2_s']+v['times_min']):raise ValueError('Diffusion coefficients and times must be positive')
  if v['orientation_samples']<128:raise ValueError('Use at least 128 orientation samples')
  return v
 def validate_dialog(self):
  try:v=self.values(); messagebox.showinfo('Validation','All inputs are valid.\nEstimated cases: '+str(len(v['Vapp_sweep_V'])*len(v['SAM_sweep_V'])))
  except Exception as e:messagebox.showerror('Validation error',str(e))
 def create_notebook(self):
  try:
   v=self.values(); nb=json.loads(TEMPLATE.read_text(encoding='utf-8')); stamp=__import__('datetime').datetime.now().strftime('%Y%m%d_%H%M%S'); self.run_dir=OUTROOT/stamp; self.run_dir.mkdir(parents=True)
   # Lists are intentionally retained in the generated configuration; the
   # scientific cells validate and convert them with np.asarray.
   cfg="from pathlib import Path\nimport numpy as np\nCFG = "+repr(v)+"\nCFG['output_dir']=Path(r'"+str(self.run_dir/'outputs')+"')\n"
   nb['cells'][1]['source']=cfg.splitlines(True); p=self.run_dir/'Configured_Molecular_Junction_Workflow.ipynb'; p.write_text(json.dumps(nb,indent=1),encoding='utf-8'); (self.run_dir/'user_inputs.json').write_text(json.dumps(v,indent=2),encoding='utf-8'); (self.run_dir/'QC_input.txt').write_text(self.qc.get('1.0','end'),encoding='utf-8')
   self.log.insert('end',f'Created {p}\n'); return p
  except Exception as e:messagebox.showerror('Create error',str(e)); self.log.insert('end',traceback.format_exc()+'\n')
 def run(self):
  p=self.create_notebook()
  if not p:return
  self.log.insert('end','Running...\n'); threading.Thread(target=self._run_worker,args=(p,),daemon=True).start()
 def _run_worker(self,p):
  runner=ROOT/'.local_plot_deps'; py=Path(os.environ.get('MOLECULAR_WORKFLOW_PYTHON',sys.executable))
  cmd=[str(py),'-c',f"import sys,json;sys.path.insert(0,r'{runner}');import matplotlib;matplotlib.use('Agg');nb=json.load(open(r'{p}',encoding='utf-8'));ns={{'display':print}};[exec(compile(''.join(c['source']),f'cell_{{i}}','exec'),ns) for i,c in enumerate(nb['cells']) if c['cell_type']=='code']"]
  cp=subprocess.run(cmd,cwd=self.run_dir,text=True,capture_output=True)
  extra=''
  if cp.returncode==0:
   report=subprocess.run([str(py),str(ROOT/'build_molecular_junction_report.py'),str(self.run_dir)],cwd=self.run_dir,text=True,capture_output=True)
   extra+=report.stdout+report.stderr
   if report.returncode==0:
    docx=self.run_dir/'outputs'/'Molecular_Junction_Report.docx'; pdf=self.run_dir/'outputs'/'Molecular_Junction_Report.pdf'
    ps=f"$ErrorActionPreference='Stop';$w=New-Object -ComObject Word.Application;$w.Visible=$false;try{{$d=$w.Documents.Open('{docx}',$false,$true);$d.ExportAsFixedFormat('{pdf}',17);$d.Close($false)}}finally{{$w.Quit()}}"
    exp=subprocess.run(['powershell.exe','-NoProfile','-Command',ps],text=True,capture_output=True)
    extra+=exp.stdout+exp.stderr
    if exp.returncode==0:extra+=f'Created {docx}\nCreated {pdf}\n'
  self.after(0,lambda:self.log.insert('end',cp.stdout+cp.stderr+extra+f'\nExit code {cp.returncode}\n'))
 def open_folder(self):
  p=self.run_dir or OUTROOT; p.mkdir(parents=True,exist_ok=True); subprocess.Popen(['explorer.exe',str(p)])

if __name__=='__main__': App().mainloop()
