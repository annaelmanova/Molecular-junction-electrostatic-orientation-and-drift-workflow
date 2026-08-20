import json, sys
from pathlib import Path

if len(sys.argv)<2:
    raise SystemExit('Usage: python repair_molecular_junction_run.py <run-directory>')
root=Path(__file__).resolve().parent
run=Path(sys.argv[1]).resolve()
template=root/'Molecular_Junction_Workflow.ipynb'
values=json.loads((run/'user_inputs.json').read_text(encoding='utf-8'))
values.update(save_png=True,save_pdf=True,dpi=300,boltzmann_chunk=128)
nb=json.loads(template.read_text(encoding='utf-8'))
cfg=("from pathlib import Path\nimport numpy as np\nCFG = "+repr(values)+
     "\nCFG['output_dir']=Path(r'"+str(run/'outputs')+"')\n")
nb['cells'][1]['source']=cfg.splitlines(True)
(run/'Configured_Molecular_Junction_Workflow.ipynb').write_text(json.dumps(nb,indent=1),encoding='utf-8')
print(run/'Configured_Molecular_Junction_Workflow.ipynb')
