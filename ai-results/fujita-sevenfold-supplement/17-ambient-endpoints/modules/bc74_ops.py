"""Independent reader of the complete written bc74 rational certificate.

Run only with danus compute DEFAULT. Exact Fraction arithmetic; no producer
program, optimizer, floating-point maximum, or cached numerical output is used.
Polynomial utilities follow sub11's earlier independent audit, with the new
old-place divided capacity reconstructed explicitly from equations (18)-(20).
"""
from .common import Q as F, record_lines, require
from pathlib import Path
from collections import Counter
import hashlib
import json
import re

TARGET = 'bc74f9bc7c0fc216'
EPS, R, U, V, T = F(1,10**6), F(9,10), F(26,25), F(27,25), F(29,25)
H0 = F(263,2000)
ROOT = (F(263,1800),F(268501,180000),F(7,30),F(2,5))
QM, QP = F(16020581569,90194313216), F(8010290789,45097156608)
SCALE = 10**9

def trim(p):
    p = list(p)
    while len(p)>1 and not p[-1]: p.pop()
    return tuple(p)
def add(*ps):
    out = [F(0)]*max(map(len,ps))
    for p in ps:
        for i,c in enumerate(p): out[i] += c
    return trim(out)
def scale(p,c): return trim([c*x for x in p])
def sub(p,q): return add(p,scale(q,-1))
def mul(p,q):
    out = [F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): out[i+j] += a*b
    return trim(out)
def power(p,n):
    out = (F(1),)
    for _ in range(n): out = mul(out,p)
    return out
def value(p,t): return sum(c*t**i for i,c in enumerate(p))
def integral(p,l,z): return 7*sum(c*(z**(i+1)-l**(i+1))/F(i+1) for i,c in enumerate(p))
ZERO, TP = (F(0),), (F(0),F(1))
def line(h,a,r): return (h-a*r,a)
def root(p): return -p[0]/p[1] if len(p)==2 and p[1] else None
def positive(p,x): return p if value(p,x)>0 else ZERO
def nonnegative(p,l,z):
    assert len(p)<=2
    return value(p,l)>=0 and value(p,z)>=0
def profiles(family,box,x):
    al,ah,ll,lh = box
    hk = F(2017,10000) if family=='Q' else F(129,625)
    h1 = max(hk,H0+al*(U-R)-EPS)
    y0 = line(H0,ah if x<R else al,R)
    y1 = line(h1,(lh+EPS-h1)/(V-U),U) if x<U else scale(TP,h1/U)
    fixed = line(ll,12,V) if x<V else line(ll,25*(ll-F(1,5)),V)
    return y0,y1,fixed
def breaks(family,box,l,z):
    points = sorted({l,z,*[x for x in [R,F(99,100),U,V] if l<x<z]})
    more = set(points)
    for a,b in zip(points,points[1:]):
        for p in profiles(family,box,(a+b)/2):
            rt = root(p)
            if rt is not None and a<rt<b: more.add(rt)
    points = sorted(more)
    thresholds = [F(0),F(1,7),F(1,6),F(2,11),QM,QP,F(1,4),F(1,3),F(1,2),F(1)]
    for a,b in zip(points,points[1:]):
        y0,y1,fp = [positive(p,(a+b)/2) for p in profiles(family,box,(a+b)/2)]
        d = sub(TP,fp)
        # The divided image uses y0, not the incident y1.
        for base,y in [(TP,y0),(TP,y1),(d,y0)]:
            candidates = [root(base)]
            candidates += [root(sub(y,scale(base,c))) for c in thresholds]
            candidates += [root(sub(base,scale(y,j))) for j in [1,2,3,4]]
            for rt in candidates:
                if rt is not None and a<rt<b: more.add(rt)
    return sorted(more)
def polynomial(family,box,label,l,z):
    mid = (l+z)/2
    y0,y1,fp = [positive(p,mid) for p in profiles(family,box,mid)]
    if label.startswith('B'):
        c = F(label[1:])
        baseline = F(1) if mid<R else F(61,100) if mid<F(99,100) else F(11,25) if mid<U else F(17,50)
        assert c>=baseline,(label,l,z,baseline)
        return scale(power(TP,6),c)
    prefix,name = label.split(':')
    assert prefix in ['0','1','D']
    y = y1 if prefix=='1' else y0
    d = sub(TP,fp) if prefix=='D' else TP
    assert nonnegative(d,l,z) and value(d,mid)>0,(label,l,z)
    assert nonnegative(y,l,z)
    def ge(c): return nonnegative(sub(y,scale(d,c)),l,z)
    def le(c): return nonnegative(sub(scale(d,c),y),l,z)
    q = 3 if prefix=='1' else 5
    if name in ['Nlo','Nhi','C']:
        assert prefix!='1' or family=='M'
        if name=='Nlo':
            assert le(F(1,7)),(label,l,z)
            return sub(power(d,6),scale(mul(power(d,6-q),power(y,q)),F(7*q,q+1)**q))
        assert ge(F(1,7)) and le(F(1,2)),(label,l,z)
        if name=='C': return scale(power(sub(d,y),6),F(7,6)**6)
        inner = add(scale(d,6*q-1),scale(y,7))
        return sub(power(d,6),scale(mul(power(d,6-q),power(inner,q)),F(1,6*(q+1))**q))
    assert prefix=='1',(family,label)
    if name=='Q':
        assert family=='Q' and le(QM),(label,l,z)
        return mul(power(sub(d,scale(y,2)),5),add(d,scale(y,10)))
    if name=='A':
        assert family=='Q' and ge(QP) and le(1),(label,l,z)
        return power(sub(d,y),6)
    assert family=='M'
    if name=='H':
        assert ge(F(1,6)),(label,l,z)
        return scale(power(d,6),F(256,729))
    if name=='Hlip':
        assert le(F(1,6)),(label,l,z)
        return sub(scale(power(d,6),F(256,729)+3),scale(mul(power(d,5),y),18))
    assert name=='R' and ge(F(2,11)),(label,l,z)
    return add(*[scale(power(positive(sub(d,scale(y,j)),mid),6),w)
                 for j,w in enumerate([1,1,-1],1)])

