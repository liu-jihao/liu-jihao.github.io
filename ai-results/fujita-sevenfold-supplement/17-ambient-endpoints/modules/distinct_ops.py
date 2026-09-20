"""Independent exact replay of the native f4 certificate. Run with danus compute.

Only native proof text is input. No producer program or prior checker is used.
"""
from .common import Q as R, record_lines, require
from pathlib import Path
from math import comb
import hashlib
import json
import re

ZERO = (R(0),)
ONE = (R(1),)
T = (R(0), R(1))

def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)

def plus(*ps):
    out = [R(0)] * max(map(len, ps))
    for p in ps:
        for i, x in enumerate(p):
            out[i] += x
    return trim(out)

def scale(p, a):
    return trim([a*x for x in p])

def minus(p, q):
    return plus(p, scale(q, -1))

def times(p, q):
    out = [R(0)] * (len(p)+len(q)-1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i+j] += x*y
    return trim(out)

def power(p, n):
    out = ONE
    for _ in range(n):
        out = times(out, p)
    return out

def value(p, x):
    ans = R(0)
    for a in reversed(p):
        ans = ans*x+a
    return ans

def integral(p, l, z):
    return sum(a*(z**(j+1)-l**(j+1))/(j+1) for j, a in enumerate(p))

def line(h, anchor, slope):
    return (h-slope*anchor, slope)

def positive(p, m):
    return p if value(p, m) > 0 else ZERO

def root(p):
    assert len(p) <= 2
    return -p[0]/p[1] if len(p) == 2 and p[1] else None

def splits(l, z, polys):
    result = {l, z}
    for p in polys:
        x = root(p)
        if x is not None and l < x < z:
            result.add(x)
    return sorted(result)

def nonnegative(p, l, z, depth=0):
    """Exact Bernstein certificate; subdivision is proof, not sampling."""
    n = len(p)-1
    affine = (l, z-l)
    pp = ZERO
    for j, a in enumerate(p):
        pp = plus(pp, scale(power(affine, j), a))
    pp = list(pp)+[R(0)]*(n+1-len(pp))
    bb = [sum(pp[j]*R(comb(k,j),comb(n,j)) for j in range(k+1))
          for k in range(n+1)]
    if min(bb) >= 0:
        return 1
    assert depth < 40 and value(p,(l+z)/2) >= 0
    return nonnegative(p,l,(l+z)/2,depth+1)+nonnegative(p,(l+z)/2,z,depth+1)

def divide(p, q):
    out = [R(0)]*max(1,len(p)-len(q)+1)
    while p != ZERO and len(p) >= len(q):
        k, a = len(p)-len(q), p[-1]/q[-1]
        out[k] += a
        p = minus(p,(R(0),)*k+scale(q,a))
    return trim(out), p

def sturm(p):
    seq = [p, trim([j*p[j] for j in range(1,len(p))])]
    while seq[-1] != ZERO:
        rem = divide(seq[-2], seq[-1])[1]
        if rem == ZERO:
            break
        seq.append(scale(rem,-1))
    return seq

def variations(seq,x):
    signs = [1 if value(p,x)>0 else -1 for p in seq if value(p,x)]
    return sum(a!=b for a,b in zip(signs,signs[1:]))

r,u,v = R(9,10),R(26,25),R(27,25)
cstar,h0 = R(23801,100000),R(263,2000)
heights = {'NQ':R(1857,10000),'NE':R(189,1000),'IQ':R(2017,10000),'IM':R(129,625)}
cm,cp = R(16020581569,90194313216),R(8010290789,45097156608)
rawF = line(R(7,30),v,R(148,125))
raw0 = line(h0,r,(cstar-h0)/(v-r))

def branches(family,l,z):
    h1 = heights[family]
    raw1 = line(h1,u,(cstar-h1)/(v-u))
    first = sorted(set(splits(l,z,[rawF,raw0,raw1])+[x for x in [r,R(99,100),u] if l<x<z]))
    out = []
    for a,b in zip(first,first[1:]):
        m = (a+b)/2
        F = positive(rawF,m)
        H0 = positive(raw0,m) if m<r else scale(T,h0/r)
        H1 = positive(raw1,m) if m<u else scale(T,h1/u)
        D = minus(T,F)
        exprs = [F,H0,H1,D,minus(F,H1),minus(T,H1)]
        for Z,Y in [(T,H0),(T,H1),(D,H1)]:
            exprs += [minus(Y,scale(Z,c)) for c in [R(0),R(1,7),R(1,6),cm,cp,R(1,2),R(1)]]
        for j in [0,1,2,3]:
            for U in [F,H1,T]:
                exprs.append(minus(minus(T,U),scale(H1,j)))
        finer = splits(a,b,exprs)
        for aa,bb in zip(finer,finer[1:]):
            out.append((aa,bb,F,H0,H1,D))
    return out

def polynomial(family,label,l,z,F,H0,H1,D):
    m = (l+z)/2
    if label=='B1':
        return power(T,6)
    if label=='P':
        assert min(value(D,l),value(D,z)) >= 0
        return power(D,6)
    if label=='J':
        assert family in ['IQ','IM']
        # Independently integrate the e-slices at theta=0.
        # D in the native J notation is the undivided degree T here.
        U = T if value(T,m)<value(H1,m) else H1
        M = F if value(F,m)>value(H1,m) else H1
        ans = power(positive(minus(T,M),m),6)
        if value(U,m)>value(F,m):
            coeffs = [(2,1)] if family=='IQ' else [(2,3),(3,-2)]
            for j,c in coeffs:
                left = power(positive(minus(minus(T,F),scale(H1,j)),m),6)
                right = power(positive(minus(minus(T,U),scale(H1,j)),m),6)
                ans = plus(ans,scale(minus(left,right),c))
        return ans
    prefix,kind = label.split(':')
    if prefix=='0':
        Z,Y,q = T,H0,5
    elif prefix=='1':
        Z,Y,q = T,H1,3
    else:
        assert prefix=='D' and family in ['NQ','NE']
        Z,Y,q = D,H1,3
    assert value(Z,m)>0
    ratios = [value(Y,x)/value(Z,x) for x in [l,z] if value(Z,x)>0]
    lo,hi = min(ratios),max(ratios)
    if kind=='Nlo':
        assert 0<=lo<=hi<=R(1,7),(family,label,l,z,lo,hi)
        return minus(power(Z,6),scale(times(power(Z,6-q),power(Y,q)),R(7*q,q+1)**q))
    if kind=='Nhi':
        assert R(1,7)<=lo<=hi<=1,(family,label,l,z,lo,hi)
        base = plus(scale(Z,6*q-1),scale(Y,7))
        return minus(power(Z,6),scale(times(power(Z,6-q),power(base,q)),R(1,6*(q+1))**q))
    if kind=='C':
        assert R(1,7)<=lo<=hi<=1
        return scale(power(minus(Z,Y),6),R(7,6)**6)
    if kind=='Q':
        assert family in ['NQ','IQ'] and 0<=lo<=hi<=cm,(family,label,l,z,lo,hi)
        return times(power(minus(Z,scale(Y,2)),5),plus(Z,scale(Y,10)))
    if kind=='A':
        assert family in ['NQ','IQ'] and cp<=lo<=hi<=1,(family,label,l,z,lo,hi)
        return power(minus(Z,Y),6)
    assert kind=='H' and family=='IM' and prefix=='1' and lo>=R(1,6)
    return scale(power(T,6),R(256,729))

