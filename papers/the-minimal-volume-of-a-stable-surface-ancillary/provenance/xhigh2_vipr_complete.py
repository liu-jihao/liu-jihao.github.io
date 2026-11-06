"""Complete weak linear derivations by explicit rational variable bounds.
No theorem is trusted here: a separate upstream checker validates output.
Memory use follows last-use indices; all source coefficients are unchanged.
Run only with [job runner].
"""
from pathlib import Path
import sys,json,time,subprocess
R=Path(__file__).resolve().parent
try:import gmpy2
except ImportError:
    dep=R/'xhigh2_gmpy2_deps'
    if not dep.exists():subprocess.run([sys.executable,'-m','pip','install','--target',str(dep),'gmpy2'],check=True,timeout=90)
    sys.path.insert(0,str(dep));import gmpy2
F=gmpy2.mpq
source=R/'xhigh2_scip_exact_45.trim.vipr';out=R/'xhigh2_scip_exact_45.complete.vipr'
store={};globalbounds={};expires={};completed=0;start=time.monotonic()
def con(tokens):
    name,sense,beta,count=tokens[:4];sense={b'G':1,b'L':-1,b'E':0}[sense]
    if count==b'OBJ':coef=dict(objective)
    else:
        nn=int(count);assert len(tokens)==4+2*nn
        coef={int(tokens[4+2*k]):F(tokens[5+2*k].decode()) for k in range(nn)}
        assert len(coef)==nn
        coef={j:a for j,a in coef.items() if a}
    return sense,F(beta.decode()),coef
def updatebounds(index,data):
    sense,beta,coef=data
    if len(coef)!=1:return
    j,a=next(iter(coef.items()))
    if not sense:return
    typ=1 if sense*a>0 else -1;b=beta/a
    key=(j,typ)
    if key not in globalbounds or typ*b>typ*globalbounds[key][0]:globalbounds[key]=(b,index)
with source.open('rb') as f,out.open('wb') as dst:
    line=f.readline();assert line.strip()==b'VER 1.0';dst.write(line)
    # All headers emitted by this solver have one variable/integer index per line.
    line=f.readline();dst.write(line);assert line.startswith(b'VAR ');nv=int(line.split()[1])
    variables=[]
    for _ in range(nv):line=f.readline();variables.append(line.strip().decode());dst.write(line)
    line=f.readline();dst.write(line);assert line.startswith(b'INT ');ni=int(line.split()[1]);ints=[]
    while len(ints)<ni:
        line=f.readline();ints.extend(map(int,line.split()));dst.write(line)
    assert len(ints)==ni
    line=f.readline();dst.write(line);assert line.startswith(b'OBJ ')
    line=f.readline();dst.write(line);t=line.split();nn=int(t[0]);assert len(t)==1+2*nn
    objective={int(t[1+2*k]):F(t[2+2*k].decode()) for k in range(nn)}
    line=f.readline();dst.write(line);assert line.startswith(b'CON ');nc=int(line.split()[1])
    for i in range(nc):
        line=f.readline();dst.write(line);data=con(line.split());store[i]=data;updatebounds(i,data)
    line=f.readline();dst.write(line);assert line.strip()==b'RTP infeas'
    line=f.readline();dst.write(line);assert line.strip()==b'SOL 0'
    line=f.readline();dst.write(line);assert line.startswith(b'DER ');nd=int(line.split()[1])
    for no in range(nd):
        i=nc+no;line=f.readline();lhs,reasonpart=line.split(b'{',1)
        reason,suffix=reasonpart.rsplit(b'}',1);t=reason.split();tail=suffix.split()
        globalflag=b'global' in tail;last=int(tail[0]);data=con(lhs.split());store[i]=data
        if len(t)>1 and t[1]==b'weak':
            assert t[0]==b'lin' and t[2]==b'{';close=t.index(b'}');nb=int(t[3]);assert close==4+4*nb
            bounds=dict(globalbounds)
            for k in range(nb):
                typ,j,ref,val=t[4+4*k:8+4*k];j=int(j);ref=int(ref);val=F(val.decode())
                s,b,a=store[ref];assert len(a)==1 and j in a
                sign=1 if typ==b'L' else -1;assert sign*s*a[j]>0 and b/a[j]==val
                bounds[(j,sign)]=(val,ref)
            p=close+1;n=int(t[p]);assert len(t)==p+1+2*n
            multipliers={int(t[p+1+2*k]):F(t[p+2+2*k].decode()) for k in range(n)}
            sense,rhs,target=data;assert sense
            coeff={};beta=F(0)
            for ref,mult in multipliers.items():
                s,b,a=store[ref];assert sense*s*mult>=0
                beta+=mult*b
                for j,v in a.items():coeff[j]=coeff.get(j,F(0))+mult*v
            for j in set(coeff)|set(target):
                diff=target.get(j,F(0))-coeff.get(j,F(0))
                if not diff:continue
                typ=1 if sense*diff>0 else -1
                bound,ref=bounds[(j,typ)];s,b,a=store[ref];mult=diff/a[j]
                assert sense*s*mult>=0
                multipliers[ref]=multipliers.get(ref,F(0))+mult;beta+=mult*b
            assert sense*(beta-rhs)>=0,(i,str(beta),str(rhs))
            multipliers={j:a for j,a in multipliers.items() if a}
            reason=b'lin '+str(len(multipliers)).encode()+b' '+b' '.join(
                str(j).encode()+b' '+str(a).encode() for j,a in sorted(multipliers.items()))
            completed+=1
        elif len(t)>1 and t[1]==b'incomplete':raise RuntimeError('LP completion required')
        dst.write(lhs+b'{ '+reason.strip()+b' } -1\n')
        if globalflag:updatebounds(i,data)
        if last>=0:expires.setdefault(last,[]).append(i)
        protected={ref for val,ref in globalbounds.values()}
        for j in expires.pop(i,[]):
            if j>=nc and j not in protected:store.pop(j,None)
        if no%20000==0:print('completed rows',no,'weak',completed,'resident',len(store),'seconds',round(time.monotonic()-start,1),flush=True)
report={'weak_completed':completed,'derivations':nd,'output_bytes':out.stat().st_size,
        'seconds':time.monotonic()-start,'independently_checked':False}
(R/'xhigh2_vipr_completion.json').write_text(json.dumps(report)+'\n');print(report,flush=True)
