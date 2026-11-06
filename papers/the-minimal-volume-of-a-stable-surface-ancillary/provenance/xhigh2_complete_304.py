"""Direct complete local enumeration from accepted cost bounds.
Run only with [job runner]; enumerate cyclic r<=304 and finite star parameters.
Keep canonical and RR-zero rows. No Pareto or packet deletion is performed.
"""
from fractions import Fraction as F
from math import gcd,lcm
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
W=(139754,122195,98574,75044,51760,31572,15936,7272)
def hj(r,q):
    z=[]
    while q:
        e=(r+q-1)//q;z.append(e);r,q=q,e*q-r
    return z
def row(label,es,av,edges,G):
    if min(av)<F(5,46): return None
    d=lcm(*(a.denominator for a in av)); aa=[int(a*d) for a in av]
    k=[e-2 for e in es]
    gamma=len(es)-sum(F(kk)*(1-a) for kk,a in zip(k,av))
    def corr(m):
        ff=[(-m*a)%d for a in aa]
        energy=sum(e*f*f for e,f in zip(es,ff))-2*sum(ff[i]*ff[j] for i,j in edges)
        return F(d*sum(kk*f for kk,f in zip(k,ff))-energy,2*d*d)
    cc=[corr(m) for m in range(2,10)]
    C=57820*(gamma+1-1/G)+sum(w*c for w,c in zip(W,cc))
    assert C>=91*(len(es)+1),(label,C)
    if C>=151733:return None
    cc.extend(corr(m) for m in range(10,222))
    assert all((c-m*(m-1)//2*cc[0]).denominator==1 for m,c in zip(range(2,222),cc))
    return dict(label=label,weights=es,discrepancies=list(map(str,av)),edges=edges,
                group=str(G),gamma=str(gamma),ell=len(es),c=list(map(str,cc)),certificate_cost=str(C))
def star(e0,arms):
    eta=e0-sum(F(q,r) for r,q in arms);sigma=sum(F(1,r) for r,q in arms)-1
    a0=sigma/eta
    es=[e0];av=[a0];edges=[]
    for r,q in arms:
        old,cur,prior=a0,(q*a0+1)/r,0
        for e in hj(r,q):
            index=len(es);es.append(e);av.append(cur);edges.append((prior,index))
            old,cur,prior=cur,e*cur-old,index
        assert cur==1
    G=4*eta/(sigma*sigma);assert G.denominator==1
    return row(['star',e0,arms],es,av,edges,G)
rows=[];cyclic=0;starcount=0
for r in range(2,305):
    for q in range(1,r):
        if gcd(r,q)!=1 or q>pow(q,-1,r):continue
        es=hj(r,q);old,cur=F(1),F(q+1,r);av=[]
        for e in es:av.append(cur);old,cur=cur,e*cur-old
        assert cur==1
        z=row(['cyclic',r,q],es,av,[(j,j+1) for j in range(len(es)-1)],F(r))
        if z:rows.append(z);cyclic+=1
    if r%75==0:print('cyclic r',r,'accepted rows',cyclic,flush=True)
# D-I: complete fork>=3 parameter list.
for e in range(3,7):
    for r in range(2,8//(e-2)+1):
        for q in range(1,r):
            if gcd(r,q)!=1 or (e-1)*r-q>9:continue
            z=star(e,[(2,1),(2,1),(r,q)])
            if z:rows.append(z);starcount+=1
# D-II: d>=2, u<=3; d=1 gives canonical D and cost>151733.
for d in range(2,10):
    for v in range(1,d):
        if gcd(d,v)!=1:continue
        for u in range(1,4):
            r=u*d+v;q=r-d
            z=star(2,[(2,1),(2,1),(r,q)])
            if z:rows.append(z);starcount+=1
# E: full finite parameter ranges, including canonical cases before the cost filter.
for r in (3,4,5):
    for e in ((2,3) if r==3 else (2,)):
        for q2 in (1,2):
            for q3 in range(1,r):
                if gcd(r,q3)!=1 or (r==3 and q2>q3):continue
                z=star(e,[(2,1),(3,q2),(r,q3)])
                if z:rows.append(z);starcount+=1
old=json.loads((ROOT/'xhigh2_sparse_catalog_full.json').read_text())
def key(r):
    if r['label'][0]=='cyclic':return ('cyclic',*min(tuple(r['weights']),tuple(reversed(r['weights']))))
    return json.dumps(r['label'],separators=(',',':'))
oldkeys={key(r) for r in old};newkeys={key(r) for r in rows}
summary=dict(rows=len(rows),cyclic=cyclic,noncyclic=starcount,
             canonical=sum(all(e==2 for e in r['weights']) for r in rows),
             zero_correction=sum(all(F(c)==0 for c in r['c']) for r in rows),
             old_missing=[r['label'] for r in rows if key(r) not in oldkeys],
             old_extra=[r['label'] for r in old if key(r) not in newkeys])
(ROOT/'xhigh2_complete_304.json').write_text(json.dumps(rows,separators=(',',':'))+'\n')
(ROOT/'xhigh2_complete_304_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
