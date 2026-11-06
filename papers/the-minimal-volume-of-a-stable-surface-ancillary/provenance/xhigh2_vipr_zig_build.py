"""Isolated compiler and upstream checker build, no system installation."""
from pathlib import Path
import subprocess,sys,urllib.request,tarfile
R=Path(__file__).resolve().parent;deps=R/'xhigh2_build_deps'
if not (deps/'ziglang/zig').exists():
    subprocess.run([sys.executable,'-m','pip','install','--target',str(deps),'ziglang==0.13.0'],check=True,timeout=100)
zig=deps/'ziglang/zig';assert zig.exists()
archive=R/'xhigh2_gmp_dev.deb'
url='https://deb.debian.org/debian/pool/main/g/gmp/libgmp-dev_6.2.1+dfsg1-1.1_amd64.deb'
if not archive.exists():urllib.request.urlretrieve(url,archive)
dev=R/'xhigh2_gmp_dev';dev.mkdir(exist_ok=True)
subprocess.run(['ar','x',str(archive)],cwd=dev,check=True,timeout=20)
data=list(dev.glob('data.tar.*'));assert len(data)==1
with tarfile.open(data[0]) as tf:
    for m in tf.getmembers():
        assert (dev/m.name).resolve().is_relative_to(dev.resolve())
    tf.extractall(dev)
source=R/'xhigh2_vipr/vipr-99eec197a49b2bd9ad33c1f4c7afabe93caa6511/code'
build=R/'xhigh2_vipr/build';build.mkdir(exist_ok=True)
(build/'CMakeConfig.hpp').write_text('#define VIPR_VERSION_MAJOR 1\n#define VIPR_VERSION_MINOR 1\n')
lib=R/'xhigh2_scip_exact_bundle/scipoptsuite-10.1.0/lib64'
# libgmpxx in the official bundle targets libstdc++; Zig uses libc++.
# Supply only the two text stream operators required by the unchanged checker.
# Values in this certificate are integer/fraction tokens. Invalid tokens fail.
shim=build/'gmp_stream_shim.cpp'
shim.write_text('''#include <iostream>
#include <string>
#include <vector>
#include <gmp.h>
std::istream& operator>>(std::istream& in, mpq_ptr q) {
  std::string s; if(!(in>>s)) return in;
  if(mpq_set_str(q,s.c_str(),10)!=0 || mpz_sgn(mpq_denref(q))==0) {
    in.setstate(std::ios::failbit); return in;
  }
  mpq_canonicalize(q); return in;
}
std::ostream& operator<<(std::ostream& out, mpq_srcptr q) {
  std::vector<char> s(mpz_sizeinbase(mpq_numref(q),10)+mpz_sizeinbase(mpq_denref(q),10)+4);
  mpq_get_str(s.data(),10,q); return out<<s.data();
}
''')
args=[str(zig),'c++','-O2','-DNDEBUG','-std=c++14','-I'+str(build),
      '-I'+str(dev/'usr/include'),'-I'+str(dev/'usr/include/x86_64-linux-gnu'),
      str(source/'viprchk.cpp'),str(shim),str(lib/'libgmp.so.10'),'-Wl,-rpath,'+str(lib),
      '-o',str(build/'viprchk')]
with (R/'xhigh2_vipr_zig_build.log').open('w') as out:
    p=subprocess.run(args,stdout=out,stderr=subprocess.STDOUT,timeout=350,
                     env={**__import__('os').environ,'ZIG_GLOBAL_CACHE_DIR':str(R/'xhigh2_zig_cache')})
print((R/'xhigh2_vipr_zig_build.log').read_text()[-4000:],flush=True)
print('returncode',p.returncode,flush=True)
