"""Independent exact audit of the printed b7b52fa8c778247e certificate.

Run only with: danus compute computation/sub11_same_prime_containment_audit.py
Reads theorem text as data; imports no author's checker or result file.
"""
from .common import Q as F, require
from functools import lru_cache
from pathlib import Path
from math import comb
import ast
import hashlib
import json
import re
import time

R, H, U, T, Q = F(9,10), F(131,1000), F(27,25), F(29,25), F(13,5)
HP, AP_MIN, CUT, EPS = F(3003,10000), F(8021,5600), F(71,400), F(1,5000)
HEIGHT = {'generic': F(26,125), 'four6': F(59,250), 'four7': F(17,80)}
KNOTS = sorted({F(0), R, U, T, Q} | {F(j,50) for j in range(1,130)})
CHECKS, FAILURES = [], []

def check(name, condition, **data):
    record = {'name': name, 'passed': bool(condition), **data}
    CHECKS.append(record)
    require(condition,name)

def positive(x):
    return max(F(0), x)

@lru_cache(maxsize=250000)
def power_integral(c, d, left, right):
    """7 times integral of (c*t+d)_+^6, including a zero slope."""
    if not c:
        return 7 * positive(d)**6 * (right-left)
    return (positive(c*right+d)**7-positive(c*left+d)**7)/c

def balanced_integral(c, d, u, v, left, right):
    """Independent binomial expansion of 7*(c*t+d)^5*(u*t+v)."""
    total = F(0)
    for j in range(6):
        coefficient = 7*comb(5,j)*c**j*d**(5-j)
        total += coefficient*(u*(right**(j+2)-left**(j+2))/F(j+2)
                              +v*(right**(j+1)-left**(j+1))/F(j+1))
    return total

def inner_knots(left, right, candidates):
    return sorted({left, right} | {x for x in candidates if left < x < right})

@lru_cache(maxsize=250000)
def kernel_integral(kind, anchor, height, slope, prime_slope, left, right):
    intercept = slope*anchor-height
    dp = None if prime_slope is None else prime_slope*T-HP
    pieces = inner_knots(left, right, [] if dp is None else [dp/prime_slope])
    total = F(0)
    for l, v in zip(pieces,pieces[1:]):
        mid = (l+v)/2
        if dp is None or prime_slope*mid <= dp:
            cd, dd = F(1), F(0)
        else:
            cd, dd = 1-prime_slope, dp
        candidates = [] if not cd else [-dd/cd]
        for threshold in (F(0), CUT, F(2,11)):
            denominator = slope-threshold*cd
            if denominator:
                candidates.append((intercept+threshold*dd)/denominator)
        subpieces = inner_knots(l,v,candidates)
        for a,b in zip(subpieces,subpieces[1:]):
            m = (a+b)/2
            degree, retained = cd*m+dd, slope*m-intercept
            if degree <= 0:
                continue
            if retained <= 0:
                total += power_integral(cd,dd,a,b)
                continue
            ratio = retained/degree
            def ej(j):
                return power_integral(cd-j*slope,dd+j*intercept,a,b)
            if kind == 'generic':
                if ratio < CUT:
                    # Split rules ensure the balanced root is positive here.
                    assert cd-2*slope == 0 or min((cd-2*slope)*a+dd+2*intercept,
                                                (cd-2*slope)*b+dd+2*intercept) >= 0
                    total += balanced_integral(cd-2*slope, dd+2*intercept,
                                               cd+10*slope, dd-10*intercept,a,b)
                else:
                    total += ej(1)+EPS*power_integral(cd,dd,a,b)
            elif kind == 'four6':
                total += (ej(1)+3*ej(2)-5*ej(3)+2*ej(4)
                          if ratio < F(2,11) else ej(1)+ej(2)-ej(3))
            elif kind == 'four7':
                total += 6*ej(2)-8*ej(3)+3*ej(4)
            else:
                raise ValueError(kind)
    assert total >= 0
    return total

@lru_cache(maxsize=50000)
def prime_integral(slope,left,right):
    dp = slope*T-HP
    knots = inner_knots(left,right,[dp/slope])
    total = F(0)
    for a,b in zip(knots,knots[1:]):
        total += (b**7-a**7 if slope*(a+b)/2 <= dp
                  else power_integral(1-slope,dp,a,b))
    return total

def box_integrals(mode,kind,box):
    incoming, outgoing = F(0), F(0)
    selections = [0]*6
    for left,right in zip(KNOTS,KNOTS[1:]):
        a0 = box[0][1 if right <= R else 0]
        a1 = box[1][1 if right <= U else 0]
        ap = box[2][1 if right <= T else 0]
        values = [right**7-left**7,
                  kernel_integral('generic',R,H,a0,None,left,right),
                  kernel_integral(kind,U,HEIGHT[kind],a1,None,left,right),
                  prime_integral(ap,left,right)]
        values.append(kernel_integral('generic',R,H,a0,ap,left,right)
                      if mode == 'old' else
                      kernel_integral(kind,U,HEIGHT[kind],a1,ap,left,right))
        if left >= T:
            values.append(power_integral(F(-52,73),F(125,73),left,right))
        best = min(range(len(values)),key=values.__getitem__)
        selections[best] += 1
        if right <= T:
            incoming += values[best]
        else:
            assert left >= T
            outgoing += values[best]
    return incoming,outgoing,selections

def volume(box):
    product = F(1)
    for left,right in box:
        product *= right-left
    return product

def leaf_boxes(tree,box,path=''):
    if isinstance(tree,int):
        yield path,tree,box
        return
    assert isinstance(tree,tuple) and len(tree)==2
    widths = [(v-l)/(1+(v+l)/2) for l,v in box]
    axis = max(range(3),key=widths.__getitem__)
    left,right = box[axis]
    mid = (left+right)/2
    assert left < mid < right
    low,high = list(box),list(box)
    low[axis],high[axis] = (left,mid),(mid,right)
    assert volume(low)+volume(high)==volume(box)
    yield from leaf_boxes(tree[0],tuple(low),path+'0')
    yield from leaf_boxes(tree[1],tuple(high),path+'1')

def scalar_checks():
    finite = [('old',F(5,38),H),('A',F(10559,44675),HEIGHT['four6']),
              ('B',F(15392,72425),HEIGHT['four7']),('C',F(13637,65475),HEIGHT['generic']),
              ('P',F(225,749),HP)]
    for name,price,height in finite:
        check(f'{name}: uniform height slack',price-height>F(1,100000),gap=price-height)
    check('fixed big beta',F(1,500)-F(7,10000)*F(7,6)**6>0,
          lower=F(1,500)-F(7,10000)*F(7,6)**6)
    check('prime slope floor',(HP-F(13,70))/(T-U)==AP_MIN)
    check('Q ratio maximum',2*F(10,11)**5<F(4,3))
    ha,hb,hc = HEIGHT['four6'],HEIGHT['four7'],HEIGHT['generic']
    tails = [('old slope >=2',R**7+F(4,3)*(R-H)**7,F(691,1000)),
             ('27 slope >=4 A',F(19,25)+(U-ha)**7/3+(U-2*ha)**7/7,F(867,1000)),
             ('27 slope >=4 B',F(19,25)+6*(U-2*hb)**7/7,F(805,1000)),
             ('27 slope >=4 C',F(19,25)+F(4,3)*(U-hc)**7/3,F(931,1000)),
             ('prime slope >=60 old',F(499,500)+F(4,3)*(T-HP-H*T/R)**7/59,F(9999,10000)),
             ('prime slope >=60 new C',F(499,500)+F(4,3)*(T-HP-hc*T/U)**7/59,F(9999,10000)),
             ('prime slope >=60 new A',F(499,500)+((T-HP-ha*T/U)**7+(T-HP-2*ha*T/U)**7)/59,F(9999,10000)),
             ('prime slope >=60 new B',F(499,500)+6*(T-HP-2*hb*T/U)**7/59,F(9999,10000))]
    for name,value,bound in tails:
        check(name,value<bound,value=value,bound=bound,slack=bound-value)
    check('A sharp branch at u',ha/U>F(2,11))
    check('A sharp branch at T after division',(ha*T/U)/(T-HP)>F(2,11))
    for name,height in [('old',H*T/R),('A',ha*T/U),('2A',2*ha*T/U),('2B',2*hb*T/U),('C',hc*T/U)]:
        check(f'{name}: positive prime-tail root',T-HP-height>0,root=T-HP-height)
    for name,a,p,bcap,expected in [('smooth fivefold',F(1),F(143,125),F(5),F(632,125)),
                                  ('quadruple e6 fivefold',F(33,25),F(76,125),F(2),F(439,190))]:
        check(f'{name}: positive cost denominator',a*U-1>0)
        forced=((8-p)*U-7)/(a*U-1)
        check(f'{name}: exact deficit obstruction',forced==expected and forced>bcap,forced=forced,cap=bcap)
    check('quintuple e6 fivefold deficit',F(7,5)>1)
    check('paid endpoint error margin',F(12949,1625)+F(1,1000)<8-F(1,100),
          remaining=8-F(12949,1625)-F(1,1000))

