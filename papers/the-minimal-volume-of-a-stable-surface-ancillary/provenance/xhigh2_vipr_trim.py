"""Keep only derivation ancestors of the final contradiction, no validation.
Produce standard VIPR with remapped references and exact last-use indices.
Run only through [job runner]. A separate checker must validate the result.
"""
from pathlib import Path
from array import array
import json,time
R=Path(__file__).resolve().parent;path=R/'xhigh2_scip_exact_45.vipr'
start=time.monotonic();offsets=array('Q');dep_offsets=array('Q',[0]);deps=array('I');kinds={};globals_=[]
with path.open('rb') as f:
    header=[];nc=None
    while True:
        line=f.readline();assert line
        if line.startswith(b'CON '):nc=int(line.split()[1])
        if line.startswith(b'DER '):
            count=int(line.split()[1]);break
        header.append(line)
    assert nc is not None
    for i in range(count):
        offsets.append(f.tell());line=f.readline();assert line
        reason=line.split(b'{',1)[1].rsplit(b'}',1)[0].split();kind=reason[0]
        if line.rstrip().endswith(b'global'):globals_.append(i)
        kinds[kind.decode()]=kinds.get(kind.decode(),0)+1
        if kind in (b'lin',b'rnd'):
            if reason[1]==b'weak':
                kinds['weak']=kinds.get('weak',0)+1
                close=reason.index(b'}');nlocal=int(reason[3]);refs=[]
                assert close==4+4*nlocal
                refs.extend(int(reason[6+4*j]) for j in range(nlocal))
                p=close+1;n=int(reason[p]);assert len(reason)==p+1+2*n
                refs.extend(int(reason[p+1+2*j]) for j in range(n))
            elif reason[1]==b'incomplete':
                kinds['incomplete']=kinds.get('incomplete',0)+1;refs=list(map(int,reason[2:]))
            else:
                n=int(reason[1]);assert len(reason)==2+2*n
                refs=[int(reason[2+2*j]) for j in range(n)]
        elif kind==b'uns':refs=list(map(int,reason[1:]));assert len(refs)==4
        elif kind in (b'asm',b'sol'):refs=[]
        else:raise ValueError(reason[:10])
        assert all(j<nc+i for j in refs)
        deps.extend(refs);dep_offsets.append(len(deps))
        if i%100000==0:print('indexed',i,'seconds',round(time.monotonic()-start,1),flush=True)
    keep=bytearray(count);stack=[count-1]+globals_
    while stack:
        i=stack.pop()
        if keep[i]:continue
        keep[i]=1
        stack.extend(j-nc for j in deps[dep_offsets[i]:dep_offsets[i+1]] if j>=nc)
    mapping=array('i',[-1])*count;selected=[]
    for i,yes in enumerate(keep):
        if yes:mapping[i]=nc+len(selected);selected.append(i)
    last=array('i',[-1])*(nc+len(selected))
    def remap(j):return j if j<nc else mapping[j-nc]
    for i in selected:
        for j in deps[dep_offsets[i]:dep_offsets[i+1]]:last[remap(j)]=mapping[i]
    out=R/'xhigh2_scip_exact_45.trim.vipr'
    with out.open('wb') as dst:
        dst.writelines(header);dst.write(f'DER {len(selected)}\n'.encode())
        for i in selected:
            f.seek(offsets[i]);line=f.readline();left,rest=line.split(b'{',1)
            reason=rest.rsplit(b'}',1)[0].split();kind=reason[0]
            if kind in (b'lin',b'rnd'):
                if reason[1]==b'weak':
                    close=reason.index(b'}');nlocal=int(reason[3])
                    for k in range(6,close,4):reason[k]=str(remap(int(reason[k]))).encode()
                    for k in range(close+2,len(reason),2):reason[k]=str(remap(int(reason[k]))).encode()
                elif reason[1]==b'incomplete':reason=reason[:2]+[str(remap(int(z))).encode() for z in reason[2:]]
                else:
                    for k in range(2,len(reason),2):reason[k]=str(remap(int(reason[k]))).encode()
            elif kind==b'uns':reason=[kind]+[str(remap(int(z))).encode() for z in reason[1:]]
            suffix=b' global' if line.rstrip().endswith(b'global') else b''
            dst.write(left+b'{ '+b' '.join(reason)+b' } '+str(last[mapping[i]]).encode()+suffix+b'\n')
report={'source_derivations':count,'kept_derivations':len(selected),'original_constraints':nc,
        'source_bytes':path.stat().st_size,'trimmed_bytes':out.stat().st_size,'kinds':kinds,
        'seconds':time.monotonic()-start,'validated':False}
(R/'xhigh2_vipr_trim.json').write_text(json.dumps(report,indent=2)+'\n')
print(report,flush=True)
