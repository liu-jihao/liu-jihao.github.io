#!/usr/bin/env python3
"""Independent exact audit of c5e4d763b03998ce. Run only via danus compute.

The printed trees are read as data. No worker computation is imported.
All decisions and comparisons use fractions; decimal values are display only.
"""
import ast
import hashlib
import json
import re
from .common import Q as F, require
from functools import lru_cache
from math import comb
from pathlib import Path


r, h, u, T, Q, H = F(9,10), F(131,1000), F(27,25), F(29,25), F(13,5), F(257,500)
HEIGHT = {"generic": F(26,125), "four6": F(59,250), "four7": F(17,80)}
KIND = {"generic": "Q", "four6": "A", "four7": "B"}
CUTS = (F(0), F(71,400), F(2,11), F(1,2), F(1), F(1,10), F(8,25))
BASE = {F(0), r, u, T} | {F(j,50) for j in range(1,58)}
checks, failures = [], []

def check(label, condition, **details):
    record = {"label": label, "passed": bool(condition), **details}
    checks.append(record)
    require(condition,label)

def positive(x):
    return max(F(0), x)

def lin_power_integral(a, b, power, q, z):
    """Integral of (a*t+b)_+**power, clipping at its actual zero."""
    if a == 0:
        return positive(b)**power * (z-q)
    zero = -b/a
    if a > 0:
        q = max(q, zero)
    else:
        z = min(z, zero)
    if z <= q:
        return F(0)
    return ((a*z+b)**(power+1) - (a*q+b)**(power+1)) / (a*(power+1))

def D(j, a, d, q, z):
    return 7 * lin_power_integral(1-j*a, j*d, 6, q, z)

def balanced_primitive(a, d, t):
    C, DD = 1-2*a, 2*d
    if C == 0:
        return 21*DD**5*t**2-35*DD**6*t
    return ((6-5*C)*(C*t+DD)**7-7*DD*(C*t+DD)**6)/C**2

def balanced_expanded(a, d, q, z):
    """Independent binomial integration of 7*((1-2a)t+2d)^5*((1+10a)t-10d)."""
    C, DD, E, FF = 1-2*a, 2*d, 1+10*a, -10*d
    result = F(0)
    for p in range(6):
        coeff = 7*comb(5,p)*C**p*DD**(5-p)
        result += coeff*(E*(z**(p+2)-q**(p+2))/F(p+2) + FF*(z**(p+1)-q**(p+1))/F(p+1))
    return result

@lru_cache(maxsize=100000)
def incoming(kind, a, anchor, height, q, z):
    d = a*anchor-height
    mid = (q+z)/2
    c = positive(a-d/mid)
    whole = z**7-q**7
    if c == 0:
        return whole
    if kind == "Q":
        if c < F(71,400):
            return balanced_primitive(a,d,z)-balanced_primitive(a,d,q)
        return D(1,a,d,q,z) + whole/5000
    if kind == "A":
        if c < F(2,11):
            return D(1,a,d,q,z)+3*D(2,a,d,q,z)-5*D(3,a,d,q,z)+2*D(4,a,d,q,z)
        return D(1,a,d,q,z)+D(2,a,d,q,z)-D(3,a,d,q,z)
    if kind == "B":
        return 6*D(2,a,d,q,z)-8*D(3,a,d,q,z)+3*D(4,a,d,q,z)
    if kind == "third":
        if c <= F(1,10):
            return whole
        return (2 if c < F(8,25) else 1)*D(1,a,d,q,z)
    raise ValueError(kind)

def knots_for(box, height):
    (l0,v0),(l1,v1),(l2,v2) = box
    knots = set(BASE)
    for a,anchor,hh in ((l0,r,h),(v0,r,h),(l1,u,height),(v1,u,height),(v2,T,H)):
        d = a*anchor-hh
        for c in CUTS:
            if a != c:
                crossing = d/(a-c)
                if 0 < crossing < T:
                    knots.add(crossing)
    return sorted(knots)

@lru_cache(maxsize=None)
def outgoing(A):
    d = A*T-H
    a1,b1,a2,b2 = 1-A,d,-F(52,73),F(125,73)
    knots = {T,Q}
    if a1 != a2:
        crossing = (b2-b1)/(a1-a2)
        if T < crossing < Q:
            knots.add(crossing)
    knots = sorted(knots)
    result = F(0)
    for q,z in zip(knots,knots[1:]):
        mid=(q+z)/2
        a,b = (a1,b1) if a1*mid+b1 <= a2*mid+b2 else (a2,b2)
        result += 7*lin_power_integral(a,b,6,q,z)
    return result

def enclosure(regime, box):
    (l0,v0),(l1,v1),(l2,v2)=box
    hh = HEIGHT[regime]
    knots = knots_for(box,hh)
    total=F(0)
    for q,z in zip(knots,knots[1:]):
        a0=v0 if z<=r else l0
        a1=v1 if z<=u else l1
        vals=(z**7-q**7, incoming("Q",a0,r,h,q,z),
              incoming(KIND[regime],a1,u,hh,q,z),incoming("third",v2,T,H,q,z))
        assert all(v>=0 for v in vals), (regime,box,q,z,vals)
        total += min(vals)
    tail=outgoing(l2)
    return min(total,F(499,500))+tail,total,tail,len(knots)-1

def leaves(node, box, path=""):
    if isinstance(node,int):
        yield node,box,path
        return
    assert isinstance(node,tuple) and len(node)==2
    widths=[(v-l)/(1+(v+l)/2) for l,v in box]
    axis=max(range(3),key=lambda j:(widths[j],-j))
    l,v=box[axis]
    mid=(l+v)/2
    left,right=list(box),list(box)
    left[axis]=(l,mid)
    right[axis]=(mid,v)
    yield from leaves(node[0],tuple(left),path+"0")
    yield from leaves(node[1],tuple(right),path+"1")

def frac_record(value):
    return {"exact":str(value),"decimal":float(value)}

def scalar_checks():
    # Finite original prices, unbounded tails, and the genuinely available third height.
    margins={
        "old_profile_height":F(5,38)-h,
        "A_profile_height":F(10559,44675)-HEIGHT["four6"],
        "B_profile_height":F(15392,72425)-HEIGHT["four7"],
        "C_profile_height":F(13637,65475)-HEIGHT["generic"],
        "third_nonlinear_height":F(87116,338285)-H/2,
    }
    for label,value in margins.items():
        threshold=F(1,100000) if label in ("old_profile_height","third_nonlinear_height") else F(1,50000)
        check(label,value>threshold,value=str(value),required=str(threshold))
    singular={2:F(225,874),3:F(9829,38050),4:F(577,2235)}
    for m,height in singular.items():
        margin=m*height-H
        check(f"singular_mass_m{m}",margin>F(4,100000),margin=str(margin))
        check(f"singular_mass_after_uniform_error_m{m}",m*(height-F(1,100000))>H)
    tails={
        "a0_ge_2":(r**7+F(4,3)*(r-h)**7,F(691,1000)),
        "a1_ge_4_A":(F(19,25)+(u-HEIGHT["four6"])**7/3+(u-2*HEIGHT["four6"])**7/7,F(867,1000)),
        "a1_ge_4_B":(F(19,25)+6*(u-2*HEIGHT["four7"])**7/7,F(805,1000)),
        "a1_ge_4_C":(F(19,25)+F(4,3)*(u-HEIGHT["generic"])**7/3,F(931,1000)),
        "A_ge_30":(F(499,500)+(T-H)**7/29,F(9997,10000)),
    }
    for name,(value,bound) in tails.items():
        check(name,value<bound,value=str(value),claimed_bound=str(bound))
    check("Q_balanced_ratio_max",2*F(10,11)**5<F(4,3),value=str(2*F(10,11)**5))
    check("big_interval_margin",F(1,500)-F(7,10000)*F(7,6)**6>0,
          margin=str(F(1,500)-F(7,10000)*F(7,6)**6))
    check("third_entire_outgoing_factor_one",H/T>F(8,25),argument=str(H/T))
    check("A_entire_outgoing_sharp_branch",HEIGHT["four6"]/u>F(2,11))
    c0=F(71,400)
    branch_gap=(1-2*c0)**5*(1+10*c0)-(1-c0)**6
    check("Q_rational_branch_gap",0<branch_gap<F(1,5000),gap=str(branch_gap))
    derivative_upper=-120*c0*F(13,21)**4+6*F(329,400)**5
    check("Q_branch_gap_decreases_through_4_over_21",derivative_upper<0,value=str(derivative_upper))
    check("factor_one_below_Q_branch_crossing",(1-2*F(4,21))**5*(1+10*F(4,21))<(1-F(4,21))**6)
    check("scaled_K2_upper_artificial_jump",2*F(9,10)**6>1)
    # Formula (17), independently expanded, including its exceptional C=0 branch.
    for a,d,q,z in [(F(1,2),F(1,8),F(1,3),F(2,5)),(F(1,3),F(1,7),F(4,5),F(6,7)),
                    (F(2),F(3,2),F(4,5),F(5,6)),(F(131,900),F(0),F(1,10),F(3,5))]:
        value=balanced_primitive(a,d,z)-balanced_primitive(a,d,q)
        check("balanced_primitive_"+str(a),value==balanced_expanded(a,d,q,z))
    # Printed nongeneration-to-height rows, recomputed from original full prices.
    row_specs=[
        ("div1",F(5,4),F(187,250),F(17654,5625)),
        ("div2",F(5,4),F(499,500),F(14029,5625)),
        ("div3",F(1201,1000),F(223,200),F(24665,9829)),
        ("div4",F(63,50),F(121,125),F(7232,2885)),
        ("five_double",F(1149,1000),F(247,250),F(28348,8321)),
        ("five3",F(623,500),F(39,40),F(28725,11134)),
        ("five4",F(33,25),F(77,100),F(3467,1328)),
        ("five5",F(69,50),F(673,1000),F(37483,15020)),
    ]
    for name,a,p,beta in row_specs:
        actual=((8-p)*T-7)/(a*T-1)
        inverse=T/(7-actual)
        check("full_price_"+name,actual==beta,beta=str(actual),inverse=str(inverse))
    check("nonlinear_common_height",min(T/(7-b) for _,_,_,b in row_specs[-3:])==F(87116,338285))
    endpoint_max=F(12949,1625)
    check("endpoint_error_gap",endpoint_max+F(1,1000)<8-F(1,100),
          exact_remaining_gap=str(8-endpoint_max-F(1,1000)))
    return {"margins":{k:frac_record(v) for k,v in margins.items()},
            "tails":{k:frac_record(v[0]) for k,v in tails.items()}}

