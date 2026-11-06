"""Audit checker execution and actual rounding path using a valid and bad toy.
This addresses a source-level concern in the upstream bool return handling.
Run through [job runner].
"""
from pathlib import Path
import subprocess,json
R=Path(__file__).resolve().parent
source=R/'xhigh2_vipr/vipr-99eec197a49b2bd9ad33c1f4c7afabe93caa6511/code/viprchk.cpp'
s=source.read_text();at=s.index('bool Constraint::round()')
print('ACTUAL ROUND SOURCE\n'+s[at:at+1000],flush=True)
exe=R/'xhigh2_vipr/build/viprchk';toy=R/'xhigh2_scip_toy.vipr'
bad=R/'xhigh2_scip_toy_bad.vipr'
bad.write_text(toy.read_text().replace('ACT_R5 L 0 1 0 1','ACT_R5 L -1 1 0 1'))
for path in (toy,bad):
    p=subprocess.run([str(exe),str(path)],capture_output=True,text=True,timeout=20)
    print('TEST',path.name,'returncode',p.returncode,p.stdout,p.stderr,flush=True)
