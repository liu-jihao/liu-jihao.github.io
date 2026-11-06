"""Check all branch domains and backward references in the final certificate.
Also confirm the compiled checker source equals its pinned upstream file.
Run through [job runner].
"""
from pathlib import Path
from fractions import Fraction as F
import urllib.request,hashlib,json,time
R=Path(__file__).resolve().parent;start=time.monotonic()
sha='99eec197a49b2bd9ad33c1f4c7afabe93caa6511'
source=R/('xhigh2_vipr/vipr-'+sha+'/code/viprchk.cpp')
with urllib.request.urlopen('https://raw.githubusercontent.com/scipopt/vipr/'+sha+'/code/viprchk.cpp',timeout=30) as src:official=src.read()
assert official==source.read_bytes()
path=R/'xhigh2_scip_exact_45.complete.vipr.opt';kinds={};assumptions=0
with path.open('rb') as f:
    while True:
        line=f.readline();assert line
        if line.startswith(b'VAR '):nv=int(line.split()[1])
        if line.startswith(b'CON '):nc=int(line.split()[1])
        if line.startswith(b'DER '):nd=int(line.split()[1]);break
    for i in range(nd):
        line=f.readline();lhs,right=line.split(b'{',1);reason=right.split(b'}',1)[0].split()
        kind=reason[0].decode();kinds[kind]=kinds.get(kind,0)+1
        if kind in ('lin','rnd'):
            n=int(reason[1]);assert len(reason)==2+2*n
            refs=[int(reason[2+2*j]) for j in range(n)]
            assert len(refs)==len(set(refs)) and all(0<=j<nc+i for j in refs)
        elif kind=='uns':assert len(reason)==5 and all(0<=int(j)<nc+i for j in reason[1:])
        elif kind=='asm':
            t=lhs.split();assert t[1] in (b'L',b'G') and F(t[2].decode()).denominator==1
            nn=int(t[3]);assert len(t)==4+2*nn
            assert all(0<=int(t[4+2*j])<nv and F(t[5+2*j].decode()).denominator==1 for j in range(nn))
            assumptions+=1
        else:raise ValueError(kind)
out={'checker_source_matches_pinned_upstream':True,'source_sha256':hashlib.sha256(official).hexdigest(),
     'all_branch_assumptions_integral':True,'all_references_backward_and_distinct':True,
     'assumptions':assumptions,'derivation_counts':kinds,'seconds':time.monotonic()-start}
(R/'xhigh2_vipr_structure_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(out,flush=True)
