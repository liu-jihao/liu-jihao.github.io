"""Independent reader/checker of the written f08 rational certificates.

No producer code, numerical cache, floating point, or optimizer is imported.
Execute only through danus compute DEFAULT. All evidence stays here.
"""
from pathlib import Path
from .common import Q
from math import comb
from itertools import product
from collections import Counter
import hashlib
import json
import re

R, U, Z, T = map(Q, ("9/10", "26/25", "27/25", "29/25"))
H0, EPS = Q("131/1000"), Q("1/100000")
QLO, QHI = Q("16020581569/90194313216"), Q("8010290789/45097156608")
DEN = 10**10
ZERO = (Q(0),) * 7
checks = Counter()

def check(condition, what, details=None):
    checks[what] += 1
    if not condition:
        raise AssertionError((what, details))

def linear_power(a, b, n, offset=0, coefficient=Q(1)):
    ans = list(ZERO)
    for j in range(n + 1):
        ans[j + offset] = coefficient * comb(n, j) * a**j * b**(n-j)
    return tuple(ans)

def add(*terms):
    return tuple(sum(x) for x in zip(*terms))

def scale(k, p):
    return tuple(k*x for x in p)

def multiply(p, q):
    a = list(ZERO)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            if x and y:
                check(i+j <= 6, "polynomial_degree")
                a[i+j] += x*y
    return tuple(a)

def monomial(k, n=6):
    return tuple(Q(k) if j == n else Q(0) for j in range(7))

def baseline(t):
    if t < R: return Q(1)
    if t < Q("99/100"): return Q("61/100")
    if t < U: return Q("11/25")
    return Q("17/50")

def poly_for_label(label, midpoint, specs, slopes, heights):
    if label.startswith("B"):
        b = Q(label[1:])
        check(b >= baseline(midpoint), "baseline_domain", (label, str(midpoint)))
        return monomial(7*b)
    index, name = label.split(":")
    i = int(index)
    kind, base, anchor, cutoff = specs[i]
    a = slopes[2*i+1] if midpoint < anchor else slopes[2*i]
    h = heights[i]
    d = a*anchor-h
    c = max(Q(0), a-d/midpoint)
    allowed = {
        "G": {"N5lo", "N5hi", "C", "E1"},
        "Q": {"Qbal", "E1"},
        "M": {"Mspan", "Mref", "N3lo", "N3hi", "C", "E1", "H3", "H3lip"},
        "D": {"Dlo", "Dhi"},
    }
    check(name in allowed[kind], "capacity_family", (kind, name))
    check(c > 0, "positive_profile", (label, str(c), str(midpoint)))

    def pj(j):
        return linear_power(1-j*a, j*d, 6, coefficient=Q(7)) if 1-j*c > 0 else ZERO

    if name.startswith("N"):
        q = int(name[1])
        if name.endswith("lo"):
            check(c <= Q(1,7), "low_generic_domain")
            return add(monomial(7), linear_power(a, -d, q, 6-q, -7*Q(7*q,q+1)**q))
        check(Q(1,7) < c < Q(1,2), "high_generic_domain")
        return add(monomial(7), linear_power(6*q-1+7*a, -7*d, q, 6-q, -Q(7,(6*(q+1))**q)))
    if name == "C":
        check(Q(1,7) < c < Q(1,2), "C_domain")
        return scale(Q(7,6)**6, pj(1))
    if name == "E1":
        check(c >= (QHI if kind == "Q" else Q(1,2)), "E1_domain")
        return pj(1)
    if name == "Qbal":
        check(c <= QLO, "Qbal_domain")
        return scale(7, multiply(linear_power(1-2*a, 2*d, 5), linear_power(1+10*a, -10*d, 1)))
    if name == "Dlo": return add(scale(6,pj(2)),scale(-8,pj(3)),scale(3,pj(4)))
    if name == "Dhi":
        check(c >= Q(1,6), "Dhi_domain")
        return monomial(7*Q(138259,900000))
    if name == "Mspan": return add(pj(1),scale(3,pj(2)),scale(-5,pj(3)),scale(2,pj(4)))
    if name == "Mref":
        check(c >= Q(2,11), "Mref_domain")
        return add(pj(1),pj(2),scale(-1,pj(3)))
    if name == "H3":
        check(c >= Q(1,6), "H3_domain")
        return monomial(7*Q(256,729))
    if name == "H3lip":
        check(c < Q(1,6), "H3lip_domain")
        return add(monomial(7*(Q(256,729)+3-18*a)),monomial(126*d,5))
    raise AssertionError(label)

def integrate_cell(left, right, label, specs, slopes, heights):
    check(0 <= left < right <= T, "cell_endpoints", (str(left), str(right)))
    cuts = {left, right}
    cuts.update(t for t in (R, Q("99/100"), U, Z) if left < t < right)
    if not label.startswith("B"):
        i = int(label.split(":")[0])
        anchor = specs[i][2]
        if left < anchor < right: cuts.add(anchor)
        for a in slopes[2*i:2*i+2]:
            d = a*anchor-heights[i]
            for c in (Q(0), Q(1,7), Q(1,6), Q(2,11), QLO, QHI, Q(1,4), Q(1,3), Q(1,2), Q(1)):
                if a != c:
                    t = d/(a-c)
                    if left < t < right: cuts.add(t)
    total = Q(0)
    polynomials = set()
    cuts = sorted(cuts)
    for l, r in zip(cuts, cuts[1:]):
        p = poly_for_label(label, (l+r)/2, specs, slopes, heights)
        polynomials.add(p)
        area = sum(v*(r**(j+1)-l**(j+1))/(j+1) for j,v in enumerate(p))
        check(area >= 0, "nonnegative_area")
        total += area
    check(len(polynomials) == 1, "one_polynomial_per_printed_cell", (label,str(left),str(right)))
    return total

def derive_heights(specs, slopes):
    heights = []
    for i, (_, base, anchor, _) in enumerate(specs):
        heights.append(max([base] + [heights[j] + slopes[2*j]*(anchor-specs[j][2])-EPS for j in range(i)]))
    return heights

def initial_leaves(grids):
    leaves = {}
    n = len(grids)
    for indices in product(*(range(1,len(grids[i])) for i in range(n))):
        tag = ".".join(map(str,indices))
        leaves[tag] = tuple(x for i,j in enumerate(indices) for x in (grids[i][j-1],grids[i][j]))
    return leaves

def apply_splits(grids, splits):
    leaves = initial_leaves(grids)
    for tag, coordinate, midpoint in splits:
        check(tag in leaves, "split_parent_present", tag)
        box = leaves.pop(tag)
        check(0 <= coordinate < len(grids), "split_coordinate")
        l, r = box[2*coordinate:2*coordinate+2]
        check(midpoint == (l+r)/2, "exact_bisection", tag)
        child0, child1 = list(box), list(box)
        child0[2*coordinate+1] = midpoint
        child1[2*coordinate] = midpoint
        for suffix, child in (("0",child0),("1",child1)):
            check(tag+suffix not in leaves, "split_child_new")
            leaves[tag+suffix] = tuple(child)
    return leaves

