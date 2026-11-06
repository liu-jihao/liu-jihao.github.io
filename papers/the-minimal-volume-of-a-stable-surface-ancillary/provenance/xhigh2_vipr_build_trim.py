from pathlib import Path
import subprocess,os
R=Path(__file__).resolve().parent;build=R/'xhigh2_vipr/build'
source=R/'xhigh2_vipr/vipr-99eec197a49b2bd9ad33c1f4c7afabe93caa6511/code'
p=subprocess.run([str(R/'xhigh2_build_deps/ziglang/zig'),'c++','-O2','-DNDEBUG','-std=c++14',
                  '-I'+str(build),str(source/'viprttn.cpp'),'-o',str(build/'viprttn')],
                 capture_output=True,text=True,timeout=180,env={**os.environ,'ZIG_GLOBAL_CACHE_DIR':str(R/'xhigh2_zig_cache')})
print('build',p.returncode,p.stdout,p.stderr,flush=True)
p=subprocess.run([str(build/'viprchk'),str(R/'xhigh2_scip_toy.vipr')],capture_output=True,text=True,timeout=30)
print('toy checker',p.returncode,p.stdout,p.stderr,flush=True)
(R/'xhigh2_vipr_toy_check.log').write_text(p.stdout+p.stderr)
