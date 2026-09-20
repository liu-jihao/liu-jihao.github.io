"""Independent exact audit of the three written smooth27 source certificates.

Run ONLY with danus compute (DEFAULT). No producer code or numeric cache is used.
All outputs stay in this computation directory; source records are read-only.
"""
from .common import Q as F, record_lines, require
from pathlib import Path
from math import comb
import hashlib
import json
import re

TARGETS = ['55cf5ca7a5348bdb', 'e1931f5132e0ca4d', 'd9363b99f849506e']
EPS, r0, u0, v0, T0 = F(1, 10**6), F(9,10), F(26,25), F(27,25), F(29,25)
QM, QP = F(16020581569,90194313216), F(8010290789,45097156608)
SCALE = 10**9

def trim(p):
    p = list(p)
    while len(p)>1 and p[-1] == 0:
        p.pop()
    return tuple(p)

def add(*ps):
    ans = [F(0)]*max(map(len,ps))
    for p in ps:
        for j,x in enumerate(p): ans[j] += x
    return trim(ans)

def mul(p,q):
    ans = [F(0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q): ans[i+j] += x*y
    return trim(ans)

def scale(p,c): return trim([x*c for x in p])
def sub(p,q): return add(p,scale(q,-1))
def power(p,n):
    ans = (F(1),)
    for _ in range(n): ans = mul(ans,p)
    return ans
def value(p,t): return sum(x*t**i for i,x in enumerate(p))
def integral(p,l,z): return 7*sum(x*(z**(j+1)-l**(j+1))/F(j+1) for j,x in enumerate(p))
Z, ONE, tpoly = (F(0),), (F(1),), (F(0),F(1))

def line(h,a,r): return (h-a*r,a)
def zero_of(p):
    return -p[0]/p[1] if len(p)==2 and p[1] else None
def positive_poly(p,t): return p if value(p,t)>0 else Z

def profile_pieces(kind,family,box,x):
    if kind=='d936':
        al,ah,ll,lh = box
        h0 = F(263,2000)
        hk = F(2017,10000) if family=='Q' else F(129,625)
        h1 = max(hk,h0+al*(u0-r0)-EPS)
        y0 = line(h0, ah if x<r0 else al, r0)
        slope = F(5)
    else:
        ll,lh = box
        h0 = F(131,1000)
        h1 = F(2777,10000) if kind=='e193' else (F(37,200) if family=='Q' else F(41,200))
        y0 = line(h0,(lh+EPS-h0)/(v0-r0),r0) if x<r0 else scale(tpoly,h0/r0)
        slope = F(12)
    y1 = line(h1,(lh+EPS-h1)/(v0-u0),u0) if x<u0 else scale(tpoly,h1/u0)
    fixed = line(ll,slope,v0) if x<v0 or kind=='d936' else line(ll,25*(ll-F(1,5)),v0)
    return y0,y1,fixed

def all_breaks(kind,family,box,l,z):
    points = sorted({l,z,*[p for p in [F(0),r0,F(99,100),u0,v0,T0] if l<p<z]})
    more = set(points)
    for a,b in zip(points,points[1:]):
        pieces = profile_pieces(kind,family,box,(a+b)/2)
        for p in pieces:
            root=zero_of(p)
            if root is not None and a<root<b: more.add(root)
    points = sorted(more)
    ratios = [F(0),F(1,7),F(1,6),F(2,11),QM,QP,F(1,4),F(1,3),F(1,2),F(1)]
    for a,b in zip(points,points[1:]):
        y0,y1,fp = [positive_poly(p,(a+b)/2) for p in profile_pieces(kind,family,box,(a+b)/2)]
        d = sub(tpoly,fp)
        for base,y in [(tpoly,y0),(tpoly,y1),(d,y1)]:
            roots = [zero_of(base)]+[zero_of(sub(y,scale(base,c))) for c in ratios]
            roots += [zero_of(sub(base,scale(y,j))) for j in [1,2,3,4]]
            for root in roots:
                if root is not None and a<root<b: more.add(root)
    return sorted(more)

def nonnegative(p,a,b):
    assert len(p)<=2, p
    return value(p,a)>=0 and value(p,b)>=0

def branch_poly(kind,family,box,label,a,b):
    mid=(a+b)/2
    y0,y1,fp = [positive_poly(p,mid) for p in profile_pieces(kind,family,box,mid)]
    if label.startswith('B'):
        c=F(label[1:])
        baseline = F(1) if mid<r0 else F(61,100) if mid<F(99,100) else F(11,25) if mid<u0 else F(17,50)
        assert c>=baseline,(label,a,b,baseline)
        return scale(power(tpoly,6),c)
    if label=='P': return power(positive_poly(sub(tpoly,fp),mid),6)
    prefix,name=label.split(':')
    y = y0 if prefix=='0' else y1
    d = sub(tpoly,fp) if prefix=='D' else tpoly
    assert nonnegative(d,a,b) and value(d,mid)>0,(kind,label,a,b)
    assert nonnegative(y,a,b)
    def ge(c): return nonnegative(sub(y,scale(d,c)),a,b)
    def le(c): return nonnegative(sub(scale(d,c),y),a,b)
    q = 3 if kind=='d936' else 5
    if name=='Nlo':
        assert le(F(1,7)),(kind,label,a,b)
        return sub(power(d,6),scale(mul(power(d,6-q),power(y,q)),F(7*q,q+1)**q))
    if name=='Nhi':
        assert ge(F(1,7)) and le(F(1,2)),(kind,label,a,b)
        inner=add(scale(d,6*q-1),scale(y,7))
        return sub(power(d,6),scale(mul(power(d,6-q),power(inner,q)),1/F(6*(q+1))**q))
    if name=='Q':
        assert family=='Q' and le(QM),(kind,label,a,b)
        return mul(power(sub(d,scale(y,2)),5),add(d,scale(y,10)))
    if name=='A':
        assert family=='Q' and ge(QP),(kind,label,a,b)
        return power(positive_poly(sub(d,y),mid),6)
    if name=='H':
        assert (prefix=='0' and kind=='d936') or family=='M'
        assert ge(F(1,6)),(kind,label,a,b)
        return scale(power(d,6),F(256,729))
    if name=='Hlip':
        assert (prefix=='0' and kind=='d936') or family=='M'
        assert le(F(1,6)),(kind,label,a,b)
        return sub(scale(power(d,6),F(256,729)+3),scale(mul(power(d,5),y),18))
    assert family=='M' and name in ['R','S'],(kind,label)
    if name=='R': assert ge(F(2,11)),(kind,label,a,b)
    weights = [1,1,-1] if name=='R' else [1,3,-5,2]
    return add(*[scale(power(positive_poly(sub(d,scale(y,j)),mid),6),w) for j,w in enumerate(weights,1)])

