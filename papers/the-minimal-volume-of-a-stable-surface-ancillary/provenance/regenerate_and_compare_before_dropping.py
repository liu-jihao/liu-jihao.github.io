"""Regenerates, with exact rational arithmetic, the 1368 local types of the finite list (running
through all parameter tuples and applying the two filters), their invariants and all 2903
constraints of the integer system, and compares them with the local rows and with the constraint
section of the proof file, as multisets of exact normalized rational coefficient vectors.
Usage: python3 regenerate_and_compare.py local_rows.json large_mld_proof.vipr
"""
from pathlib import Path
import sys
from fractions import Fraction as F
from math import gcd,lcm
from collections import Counter
import json,time,hashlib
ROWS,CERT=Path(sys.argv[1]),Path(sys.argv[2]);start=time.monotonic()
W=(139754,122195,98574,75044,51760,31572,15936,7272)
def hj(r,q):
    out=[]
    while q:e=(r+q-1)//q;out.append(e);r,q=q,e*q-r
    return out
def data(label):
    if label[0]=='cyclic':
        r,q=label[1:];weights=hj(r,q);av=[];left,cur=F(1),F(q+1,r)
        for e in weights:av.append(cur);left,cur=cur,e*cur-left
        assert cur==1;edges=[(j,j+1) for j in range(len(weights)-1)];G=F(r)
    else:
        e,arms=label[1:];sigma=sum(F(1,r) for r,q in arms)-1;eta=e-sum(F(q,r) for r,q in arms)
        a0=sigma/eta;G=4*eta/sigma**2;weights=[e];av=[a0];edges=[]
        for r,q in arms:
            left,cur=a0,(q*a0+1)/r;prior=0
            for e in hj(r,q):
                j=len(weights);weights.append(e);av.append(cur);edges.append((prior,j))
                left,cur=cur,e*cur-left;prior=j
            assert cur==1
    if min(av)<F(5,46):return None
    assert G.denominator==1
    k=[e-2 for e in weights];b=[1-a for a in av]
    for i in range(len(weights)):
        assert weights[i]*b[i]-sum(b[v] if u==i else b[u] for u,v in edges if u==i or v==i)==k[i]
    den=lcm(*(a.denominator for a in av));ai=[int(a*den) for a in av]
    gamma=len(weights)-sum(v*z for v,z in zip(k,b))
    def c(m):
        ff=[(-m*a)%den for a in ai]
        return F(den*sum(v*z for v,z in zip(k,ff))-sum(e*z*z for e,z in zip(weights,ff))+2*sum(ff[i]*ff[j] for i,j in edges),2*den*den)
    cc=[c(m) for m in range(2,46)]
    cost=57820*(gamma+1-1/G)+sum(w*cc[i] for i,w in enumerate(W))
    if cost>=151733:return None
    assert cost>=91*(len(weights)+1)
    return {'label':label,'gamma':gamma,'group':G,'ell':len(weights),'c':cc,'cost':cost}
labels=[]
for r in range(2,305):
    for q in range(1,r):
        if gcd(r,q)==1 and q<=pow(q,-1,r):labels.append(['cyclic',r,q])
for e in range(3,7):
    for r in range(2,8//(e-2)+1):
        for q in range(1,r):
            if gcd(r,q)==1 and (e-1)*r-q<=9:labels.append(['star',e,[(2,1),(2,1),(r,q)]])
for d in range(2,10):
    for t in range(1,d):
        if gcd(d,t)!=1:continue
        for u in range(1,4):labels.append(['star',2,[(2,1),(2,1),(u*d+t,(u-1)*d+t)]])
for r in (3,4,5):
    for e in ((2,3) if r==3 else (2,)):
        for q in (1,2):
            for t in range(1,r):
                if gcd(r,t)==1 and (r!=3 or q<=t):labels.append(['star',e,[(2,1),(3,q),(r,t)]])
rows=[]
for label in labels:
    z=data(label)
    if z:rows.append(z)
n=len(rows);nc=n+28
old=json.loads(ROWS.read_text());assert len(old)==n
for a,b in zip(rows,old):
    assert json.dumps(a['label'])==json.dumps(b['label'])
    assert a['gamma']==F(b['gamma']) and a['group']==F(b['group']) and a['ell']==b['ell']
    assert a['cost']==F(b['certificate_cost']) and a['c']==list(map(F,b['c'][:44]))
print('complete regenerated catalogue',n,'seconds',round(time.monotonic()-start,1),flush=True)
cons=[]
def add(s,a,b):cons.append((s,{j:F(v) for j,v in a.items() if v},F(b)))
s=[[c-m*(m-1)//2*z['c'][0] for m,c in zip(range(2,46),z['c'])] for z in rows]
assert all(a.denominator==1 for row in s for a in row)
add('E',{j:a[1] for j,a in enumerate(s)},2)
for m in range(4,46):
    d={j:a[m-2]-2*a[m-3]+a[m-4] for j,a in enumerate(s)}
    for k,v in ((m,-1),(m-1,2),(m-2,-1)):
        if k>=18:d[n+k-18]=v
    add('E',d,1)
add('L',{j:z['c'][0] for j,z in enumerate(rows)},-1-F(1,84235998))
add('L',{j:-z['c'][0] for j,z in enumerate(rows)},1+F(1,48983))
add('L',{j:z['gamma']-z['c'][0] for j,z in enumerate(rows)},10)
add('L',{j:z['gamma']+1-1/z['group']-F(4,3)*z['c'][0] for j,z in enumerate(rows)},F(40,3))
add('L',{j:z['ell']+1 for j,z in enumerate(rows)},1667)
add('L',{j:z['cost'] for j,z in enumerate(rows)},151733)
for a in range(18,46):
    for k in range(2,45//a+1):add('L',{n+a-18:1,n+k*a-18:-1},0)
    for b in range(a+1,46):
        if gcd(a,b)<=17 and lcm(a,b)<=156:add('L',{n+a-18:1,n+b-18:1},1)
for a in range(18,46):
    for b in range(a,46-a):
        d={n+a-18:1,n+b-18:1,n+a+b-18:-1}
        if a==b:d[n+a-18]=2
        add('L',d,1)
base=len(cons)
ub=[min(833,(151733*z['cost'].denominator-1)//z['cost'].numerator,1667//(z['ell']+1)) for z in rows]+[1]*28
for j,u in enumerate(ub):add('G',{j:1},0);add('L',{j:1},u)
def key(s,a,b):
    # Normalize by a positive rational scale only. Equality may also change sign.
    first=a[min(a)] if a else b
    scale=abs(first) if first else F(1)
    if s=='E' and first<0:scale=-scale
    return s,tuple((j,str(v/scale)) for j,v in sorted(a.items())),str(b/scale)
expected=Counter(key(*z) for z in cons)
cert=CERT
with cert.open() as f:
    assert f.readline().strip()=='VER 1.0'
    assert f.readline().strip()==f'VAR {nc}'
    var=[int(f.readline().strip().removeprefix('t_X')) for _ in range(nc)]
    assert sorted(var)==list(range(nc))
    assert f.readline().strip()==f'INT {nc}'
    iv=[]
    while len(iv)<nc:iv.extend(map(int,f.readline().split()))
    assert sorted(iv)==list(range(nc))
    assert f.readline().startswith('OBJ ');f.readline()
    info=f.readline().split();assert info[0]=='CON';count=int(info[1]);assert count==len(cons)
    actual=Counter()
    for _ in range(count):
        t=f.readline().split();sense=t[1];b=F(t[2]);nn=int(t[3]);assert len(t)==4+2*nn
        a={var[int(t[4+2*k])]:F(t[5+2*k]) for k in range(nn)};assert len(a)==nn
        a={j:v for j,v in a.items() if v}
        actual[key(sense,a,b)]+=1
assert actual==expected,{'missing':len(expected-actual),'extra':len(actual-expected)}
out={'complete_local_rows':n,'variables':nc,'integer_variables':nc,'base_constraints':base,
     'bounds':2*nc,'certificate_constraints':count,'exact_model_match':True,
     'seconds':time.monotonic()-start}
print(json.dumps(out,indent=2),flush=True)
