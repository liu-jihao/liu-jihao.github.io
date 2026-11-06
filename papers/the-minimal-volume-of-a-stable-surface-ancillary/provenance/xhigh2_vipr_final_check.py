"""Official proof trimming then independent exact GMP certificate checking.
Run through [job runner] default tier.
"""
from pathlib import Path
import subprocess,time,json,hashlib,sys
R=Path(__file__).resolve().parent;build=R/'xhigh2_vipr/build'
cert=R/'xhigh2_scip_exact_45.complete.vipr';start=time.monotonic()
if len(sys.argv)<2 or sys.argv[1]!='checkonly':
    with (R/'xhigh2_vipr_official_trim.log').open('w') as log:
        p=subprocess.run([str(build/'viprttn'),str(cert)],stdout=log,stderr=subprocess.STDOUT,timeout=200)
    print('trim_returncode',p.returncode,flush=True)
# Upstream viprttn at the pinned commit initializes rs=-1 and never updates
# it, even on success. Only the independent viprchk result is authoritative.
small=Path(str(cert)+'.opt');print('trimmed_bytes',small.stat().st_size,flush=True)
with (R/'xhigh2_vipr_final_check.log').open('w') as log:
    p=subprocess.run([str(build/'viprchk'),str(small)],stdout=log,stderr=subprocess.STDOUT,timeout=340)
log=(R/'xhigh2_vipr_final_check.log').read_text()
digest=hashlib.sha256()
with small.open('rb') as src:
    for block in iter(lambda:src.read(1<<20),b''):digest.update(block)
result={'returncode':p.returncode,'verified':p.returncode==0 and 'Successfully verified infeasibility.' in log,
        'certificate':small.name,'bytes':small.stat().st_size,'sha256':digest.hexdigest(),
        'seconds':time.monotonic()-start}
(R/'xhigh2_vipr_final_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(log[-3000:],flush=True);print(result,flush=True)
