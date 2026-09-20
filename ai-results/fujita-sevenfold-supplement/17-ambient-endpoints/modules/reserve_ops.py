"""Independent exact audit; execute only with danus compute (DEFAULT).

Polynomial expansion integrates the kernels independently of the printed
antiderivatives. Source files are read only. No floats or quadrature are used.
"""
from .common import Q as F, require
from math import comb
from pathlib import Path
import ast
import hashlib
import json
import re

def check(name,truth,evidence=None):
    require(truth,name)

def mul(p, q):
    z = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            z[i+j] += a*b
    return z


def power(p, n):
    z = [F(1)]
    for _ in range(n):
        z = mul(z, p)
    return z


def integral(p, l, v):
    return sum(F(7)*x*(v**(i+1)-l**(i+1))/F(i+1)
               for i, x in enumerate(p))


def linear_integral(slope, intercept, l, v):
    """Integral of 7 times the positive sixth power of an affine root."""
    points = {l, v}
    if slope and l < -intercept/slope < v:
        points.add(-intercept/slope)
    total = F(0)
    points = sorted(points)
    for x, y in zip(points, points[1:]):
        if slope*(x+y)/2 + intercept > 0:
            total += integral(power([intercept, slope], 6), x, y)
    return total


r0, h0, u, r, q = map(F, ['9/10', '131/1000', '27/25', '29/25', '13/5'])
cstar, eps = F(71, 400), F(1, 5000)
heights = {'A': F(59, 250), 'B': F(17, 80), 'C': F(26, 125)}
branch_gap = (1-2*cstar)**5*(1+10*cstar)-(1-cstar)**6
check('Q rational majorant branch gap',
      0 < branch_gap == F(320*129**5*111-329**6, 400**6) < eps)
check('Q rational majorant derivative bound',
      -120*cstar*F(13,21)**4+6*F(329,400)**5 < 0)


def cap_integral(kind, anchor, height, slope, l, v):
    d = slope*anchor-height
    switches = [F(0), F(1), F(1, 2), F(1, 3), F(1, 4)]
    switches += [cstar] if kind in ('Q', 'C') else [F(2, 11)]
    points = {l, v}
    for c in switches:
        if slope != c and l < d/(slope-c) < v:
            points.add(d/(slope-c))
    points = sorted(points)
    total = F(0)
    for x, y in zip(points, points[1:]):
        mid = (x+y)/2
        height_mid = max(F(0), slope*mid-d)
        c = height_mid/mid
        if not height_mid:
            total += y**7-x**7
            continue
        if kind in ('Q', 'C'):
            if c <= cstar:
                total += integral(mul(power([2*d, 1-2*slope], 5),
                                      [-10*d, 1+10*slope]), x, y)
            else:
                total += linear_integral(1-slope, d, x, y)
                total += eps*(y**7-x**7)
        else:
            terms = ([(1, 1), (2, 3), (3, -5), (4, 2)] if c < F(2, 11)
                     else [(1, 1), (2, 1), (3, -1)]) if kind == 'A' else [(2, 6), (3, -8), (4, 3)]
            total += sum(k*linear_integral(1-j*slope, j*d, x, y) for j, k in terms)
    return total


def box_integral(kind, box):
    l0, v0, l1, v1 = box
    H = heights[kind]
    points = {F(0), r0, u, r}
    if kind != 'C':
        points.update(F(j, 50) for j in range(1, 58))
    total = F(0)
    points = sorted(points)
    for l, v in zip(points, points[1:]):
        mid = (l+v)/2
        a = v0 if mid < r0 else l0
        b = v1 if mid < u else l1
        if kind == 'C':
            sub = {l, v}
            if a != b:
                crossing = (a*r0-h0-b*u+H)/(a-b)
                if l < crossing < v:
                    sub.add(crossing)
            sub = sorted(sub)
            for x, y in zip(sub, sub[1:]):
                z = (x+y)/2
                if h0+a*(z-r0) >= H+b*(z-u):
                    total += cap_integral('Q', r0, h0, a, x, y)
                else:
                    total += cap_integral('Q', u, H, b, x, y)
        else:
            total += min(cap_integral('Q', r0, h0, a, l, v),
                         cap_integral(kind, u, H, b, l, v))
    return total


def leaves(tree, box):
    if isinstance(tree, int):
        yield box, tree
        return
    assert isinstance(tree, tuple) and len(tree) == 2
    l0, v0, l1, v1 = box
    if v0-l0 >= v1-l1:
        middle = (l0+v0)/2
        halves = ((l0, middle, l1, v1), (middle, v0, l1, v1))
    else:
        middle = (l1+v1)/2
        halves = ((l0, v0, l1, middle), (l0, v0, middle, v1))
    yield from leaves(tree[0], halves[0])
    yield from leaves(tree[1], halves[1])


