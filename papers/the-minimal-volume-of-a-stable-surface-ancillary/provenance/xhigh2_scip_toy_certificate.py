"""Check rational input and VIPR emission on 1/2*x=1/4, x integral.
Run through [job runner]. This verifies backend capability, not the RR model.
"""
from pathlib import Path
import json,subprocess
R=Path(__file__).resolve().parent
exe=json.loads((R/'xhigh2_scip_bundle_probe.json').read_text())['binary']
model=R/'xhigh2_scip_toy.mps'
model.write_text('''NAME TOY
ROWS
 N OBJ
 E HALF
COLUMNS
 MARK0 'MARKER' 'INTORG'
 X OBJ 0
 X HALF 1/2
 MARK1 'MARKER' 'INTEND'
RHS
 RHS1 HALF 1/4
BOUNDS
 LO BND X 0
 UP BND X 1
ENDATA
''')
commands=['set exact enable TRUE','set certificate filename '+str(R/'xhigh2_scip_toy.vipr'),
          'set presolving emphasis off','set separating emphasis off','read '+str(model),
          'optimize','quit']
args=[exe]
for c in commands:args+=['-c',c]
p=subprocess.run(args,capture_output=True,text=True,timeout=45,cwd=R)
(R/'xhigh2_scip_toy.log').write_text(p.stdout+p.stderr)
print(p.stdout,p.stderr,flush=True)
cert=R/'xhigh2_scip_toy.vipr'
if cert.exists():print('CERTIFICATE\n'+cert.read_text(),flush=True)
