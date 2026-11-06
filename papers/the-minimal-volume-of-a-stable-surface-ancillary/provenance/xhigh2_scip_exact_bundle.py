"""Inspect an isolated official SCIP bundle for exact proof generation.
Run through [job runner]; no system installation or host configuration writes.
"""
from pathlib import Path
import urllib.request,tarfile,hashlib,subprocess,json,os
R=Path(__file__).resolve().parent
name='scipoptsuite-10.1.0-glibc2_28-amd64.tgz'
url='https://github.com/scipopt/scip/releases/download/v10.1.0/'+name
expected='330a884428af4d644370cdf36ba405b879a9e17341c5ad2646917ba8972767eb'
archive=R/name;root=R/'xhigh2_scip_exact_bundle';root.mkdir(exist_ok=True)
if not archive.exists():
    with urllib.request.urlopen(url,timeout=45) as src,archive.open('wb') as dst:
        size=0
        while True:
            block=src.read(1<<20)
            if not block:break
            size+=len(block);assert size<=110_000_000;dst.write(block)
digest_obj=hashlib.sha256()
with archive.open('rb') as src:
    for block in iter(lambda:src.read(1<<20),b''):digest_obj.update(block)
digest=digest_obj.hexdigest()
assert digest==expected
with tarfile.open(archive) as tf:
    members=tf.getmembers();assert sum(z.size for z in members)<1_000_000_000
    for member in members:
        assert (root/member.name).resolve().is_relative_to(root.resolve())
        if member.issym():
            assert (root/member.name).parent.joinpath(member.linkname).resolve().is_relative_to(root.resolve())
        if member.islnk():
            assert (root/member.linkname).resolve().is_relative_to(root.resolve())
    tf.extractall(root)
bins=list(root.rglob('bin/scip'));assert len(bins)==1
exe=bins[0]
commands=['set exact enable TRUE','set certificate filename '+str(R/'xhigh2_scip_probe.vipr'),'quit']
args=[str(exe)]
for command in commands:args.extend(['-c',command])
proc=subprocess.run(args,capture_output=True,text=True,timeout=45,cwd=R)
out={'archive':name,'sha256':digest,'binary':str(exe),'returncode':proc.returncode,
     'stdout':proc.stdout,'stderr':proc.stderr}
(R/'xhigh2_scip_bundle_probe.json').write_text(json.dumps(out,indent=2)+'\n')
print(proc.stdout,proc.stderr,flush=True)
