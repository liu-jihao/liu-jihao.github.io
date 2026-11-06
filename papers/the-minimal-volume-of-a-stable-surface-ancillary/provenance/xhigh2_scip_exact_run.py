"""One bounded exact-SCIP proof-logging attempt on the complete degree45 model.
Run through [job runner]. No verdict is trusted without certificate audit.
"""
from pathlib import Path
import json,subprocess,os,time
R=Path(__file__).resolve().parent
exe=json.loads((R/'xhigh2_scip_bundle_probe.json').read_text())['binary']
model=R/'xhigh2_exact_prefix_45.mps';assert model.exists()
cert=R/'xhigh2_scip_exact_45.vipr';log=R/'xhigh2_scip_exact_45.log'
commands=['set exact enable TRUE','set certificate filename '+str(cert),
          'set presolving emphasis off','set limits time 180','set limits memory 1400',
          'set parallel maxnthreads 1','read '+str(model),'optimize','quit']
args=[exe]
for command in commands:args.extend(['-c',command])
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
start=time.monotonic()
with log.open('w') as output:
    proc=subprocess.run(args,stdout=output,stderr=subprocess.STDOUT,timeout=210,cwd=R,env=env)
lines=log.read_text().splitlines()
print('\n'.join(lines[-45:]),flush=True)
print('returncode',proc.returncode,'elapsed',time.monotonic()-start,
      'certificate_bytes',cert.stat().st_size if cert.exists() else 0,flush=True)
