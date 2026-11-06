"""Build pinned official VIPR checker in computation, through [job runner]."""
from pathlib import Path
import urllib.request,tarfile,subprocess,json,time,shutil
R=Path(__file__).resolve().parent
sha='99eec197a49b2bd9ad33c1f4c7afabe93caa6511'
archive=R/'xhigh2_vipr_source.tgz';root=R/'xhigh2_vipr';root.mkdir(exist_ok=True)
if not archive.exists():
    with urllib.request.urlopen('https://codeload.github.com/scipopt/vipr/tar.gz/'+sha,timeout=30) as src:
        data=src.read(20_000_001);assert len(data)<=20_000_000;archive.write_bytes(data)
with tarfile.open(archive) as tf:
    for m in tf.getmembers():
        assert (root/m.name).resolve().is_relative_to(root.resolve())
        assert not m.issym() and not m.islnk()
    tf.extractall(root)
source=root/('vipr-'+sha)/'code';build=root/'build';build.mkdir(exist_ok=True)
print({name:shutil.which(name) for name in ('gcc','g++','clang++','cmake','make')},flush=True)
if shutil.which('cmake') is None:
    import sys
    subprocess.run([sys.executable,'-m','pip','install','--target',str(R/'xhigh2_build_deps'),'cmake'],check=True,timeout=90)
    cmake=str(R/'xhigh2_build_deps/cmake/data/bin/cmake')
else:cmake=shutil.which('cmake')
log=R/'xhigh2_vipr_build.log'
with log.open('w') as out:
    p=subprocess.run([cmake,'-S',str(source),'-B',str(build),'-DVIPRCOMP=off','-DCMAKE_BUILD_TYPE=Release'],stdout=out,stderr=subprocess.STDOUT,timeout=60)
    assert p.returncode==0
    p=subprocess.run([cmake,'--build',str(build),'-j1'],stdout=out,stderr=subprocess.STDOUT,timeout=430)
print(log.read_text()[-3500:],flush=True)
print('build_returncode',p.returncode,flush=True)
