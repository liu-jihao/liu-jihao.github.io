"""Independent, exact check of the written certificate in 29467d401ea8b3df.

Finite question, framed before execution: do all 14 printed rational partitions,
all 16 enlarged beta intervals, all B-cell Bernstein coefficient enclosures,
all S-cell integrals, and all cell/row integer enclosures follow from the native
proof's displayed polynomials? Run ONLY with `danus compute` in DEFAULT tier.
No producer program or producer arithmetic output is imported. Only the frozen
native Markdown certificate is an input. All generated material stays here.
The output is audit evidence, not a formal proof or a fact submission.
"""

from .common import Q, require, record_lines
from functools import lru_cache
from hashlib import sha256
from math import comb
from pathlib import Path
import json
import re
import time

R, U, HEIGHT, HP = Q(9, 10), Q(26, 25), Q(263, 2000), Q(1, 5)
ZERO, ONE, T = (Q(0),), (Q(1),), (Q(0), Q(1))
BETAS = [(Q(3, 5) + Q(j, 40), Q(3, 5) + Q(j + 1, 40)) for j in range(16)]
ZLEVELS = [Q(1), Q(1, 2), Q(1, 6), Q(0)]


def trim(a):
    a = list(a)
    while len(a) > 1 and not a[-1]:
        a.pop()
    return tuple(a)


def plus(a, b):
    return trim(tuple((a[i] if i < len(a) else 0) +
                      (b[i] if i < len(b) else 0)
                      for i in range(max(len(a), len(b)))))


def scale(a, b):
    return trim(tuple(v * b for v in a))


def minus(a, b):
    return plus(a, scale(b, -1))


def times(a, b):
    ans = [Q(0)] * (len(a) + len(b) - 1)
    for i, v in enumerate(a):
        for j, w in enumerate(b):
            ans[i + j] += v * w
    return trim(ans)


@lru_cache(maxsize=None)
def power(a, n):
    if n == 0:
        return ONE
    return times(a, power(a, n - 1))


def at(a, t):
    val = Q(0)
    for v in reversed(a):
        val = val * t + v
    return val


def roots(a, l, v):
    assert len(a) <= 2
    if len(a) < 2 or not a[1]:
        return []
    root = -a[0] / a[1]
    return [root] if l < root < v else []


def line(slope, anchor, height):
    return (height - slope * anchor, slope)


def ceilq(a):
    return -((-a.numerator) // a.denominator)


def integral(a, l, v):
    return sum((c * (v ** (m + 1) - l ** (m + 1)) / (m + 1)
                for m, c in enumerate(a)), Q(0))


def bplus(a, b):
    ans = dict(a)
    for monomial, coeff in b.items():
        ans[monomial] = ans.get(monomial, Q(0)) + coeff
    return {k: v for k, v in ans.items() if v}


def bscale(a, q):
    return {k: v * q for k, v in a.items() if v * q}


def btimes(a, b):
    ans = {}
    for (i, j), v in a.items():
        for (k, m), w in b.items():
            monomial = (i + k, j + m)
            ans[monomial] = ans.get(monomial, Q(0)) + v * w
    return {k: v for k, v in ans.items() if v}


def bpower(a, n):
    ans = {(0, 0): Q(1)}
    for _ in range(n):
        ans = btimes(ans, a)
    return ans


@lru_cache(maxsize=None)
def antiderivative(kind, h, B, C):
    # Independent expansion in global t and slice coordinate e.
    x = {(1, 0): Q(1), (0, 1): Q(-1)}
    y = {(i, 0): coeff / C for i, coeff in enumerate(h) if coeff}
    y[(0, 1)] = -B / C
    xy = bplus(x, bscale(y, -1))
    if kind == 'axis':
        poly = bpower(xy, 5)
    elif kind == 'convex':
        poly = bscale(bpower(xy, 5), Q(6, 5) ** 5)
    elif kind == 'centroid':
        argument = bscale(bplus(bscale(x, 19), bscale(y, 6)), Q(1, 25))
        poly = bplus(bpower(x, 5), bscale(btimes(x, bpower(argument, 4)), -1))
    elif kind == 'low':
        poly = bplus(bpower(x, 5),
                     bscale(btimes(x, bpower(y, 4)), -Q(24, 5) ** 4))
    else:
        raise ValueError(kind)
    return tuple((i, j + 1, coeff / (j + 1))
                 for (i, j), coeff in sorted(poly.items()))


@lru_cache(maxsize=None)
def endpoint_integral(kind, h, B, C, endpoint):
    ans = ZERO
    for i, j, coeff in antiderivative(kind, h, B, C):
        term = (Q(0),) * i + scale(power(endpoint, j), coeff)
        ans = plus(ans, term)
    return ans


def slice_integral(kind, h, B, C, lower, upper):
    return minus(endpoint_integral(kind, h, B, C, upper),
                 endpoint_integral(kind, h, B, C, lower))


def bernstein(poly, l, v):
    # Binomial substitution to z, followed by t^k = sum B_i^6 comb(i,k)/comb(6,k).
    assert len(poly) <= 7
    zpower = [sum((poly[m] * comb(m, k) * l ** (m - k)
                   for m in range(k, len(poly))), Q(0)) * (v - l) ** k
              for k in range(7)]
    result = [sum((zpower[k] * Q(comb(i, k), comb(6, k))
                   for k in range(i + 1)), Q(0)) for i in range(7)]
    # A separate integral identity checks the conversion, exactly.
    assert (v - l) * sum(result) == 7 * integral(poly, l, v)
    assert result[0] == at(poly, l) and result[6] == at(poly, v)
    return result


def sign_uniform(a, l, v):
    lo, hi, mid = at(a, l), at(a, v), at(a, (l + v) / 2)
    assert (lo >= 0 and hi >= 0) if mid > 0 else (
        (lo <= 0 and hi <= 0) if mid < 0 else lo == hi == 0)


def profiles(row, l, v):
    mid = (l + v) / 2
    h = line(row['a1'] if mid < R else row['a0'], R, HEIGHT)
    f = line(Q(3), U, HP) if row['case'] == 'prime3' else ZERO
    sign_uniform(h, l, v)
    sign_uniform(f, l, v)
    return h if at(h, mid) > 0 else ZERO, f if at(f, mid) > 0 else ZERO


def edge(h, B, C, z):
    assert B - C * z > 0
    return scale(minus(h, scale(T, C * z)), 1 / (B - C * z))


def partition(row):
    initial = [Q(0), R, Q(99, 100), U]
    positive = set(initial)
    for l, v in zip(initial, initial[1:]):
        h = line(row['a1'] if (l + v) / 2 < R else row['a0'], R, HEIGHT)
        f = line(Q(3), U, HP) if row['case'] == 'prime3' else ZERO
        positive.update(roots(h, l, v))
        positive.update(roots(f, l, v))
    positive = sorted(positive)
    branch = set(positive)
    for l, v in zip(positive, positive[1:]):
        h, f = profiles(row, l, v)
        branch.update(roots(minus(f, T), l, v))
        for lowbeta, B in BETAS:
            C = 1 - lowbeta
            branch.update(roots(minus(h, scale(T, B)), l, v))
            for z in ZLEVELS:
                branch.update(roots(minus(f, edge(h, B, C, z)), l, v))
    branch = sorted(branch)
    refined = [Q(0)]
    for l, v in zip(branch, branch[1:]):
        n = max(1, ceilq(40 * (v - l)))
        refined.extend(l + (v - l) * Q(k, n) for k in range(1, n + 1))
    return refined, branch


def interval_polynomial(h, f, B, C, l, v):
    mid = (l + v) / 2
    for comparison in [h, f, minus(f, T), minus(h, scale(T, B))]:
        sign_uniform(comparison, l, v)
    metadata = {'B': str(B), 'C': str(C)}
    if at(h, mid) == 0:
        return power(minus(T, f), 6), '-', dict(metadata, ordinary=True)
    edges = [edge(h, B, C, z) for z in ZLEVELS]
    for ep in edges:
        sign_uniform(minus(f, ep), l, v)
    metadata['edges_z_1_half_sixth_zero'] = [[str(q) for q in ep] for ep in edges]
    if at(f, mid) >= mid or at(h, mid) >= B * mid:
        return ZERO, 'X', dict(metadata, empty=True)
    for left, right in zip(edges, edges[1:] + [T]):
        sign_uniform(minus(right, left), l, v)
        assert at(left, mid) < at(right, mid)
    e1, ehalf, esixth, e0 = edges

    def larger(a, b):
        sign_uniform(minus(a, b), l, v)
        return a if at(a, mid) >= at(b, mid) else b

    lower_tail = larger(f, e0)
    poly = power(minus(T, lower_tail), 6)
    slices = []
    choices = []
    for kind, left, right in [('axis', e1, ehalf), ('middle', ehalf, esixth), ('low', esixth, e0)]:
        lower = larger(f, left)
        sign_uniform(minus(right, lower), l, v)
        if at(lower, mid) >= at(right, mid):
            continue
        if kind == 'middle':
            center = scale(plus(lower, right), Q(1, 2))
            halves = [(lower, center), (center, right)]
            for aa, bb in halves:
                centroid = slice_integral('centroid', h, B, C, aa, bb)
                convex = slice_integral('convex', h, B, C, aa, bb)
                take_convex = at(convex, mid) < at(centroid, mid)
                selected = convex if take_convex else centroid
                choices.append(int(take_convex))
                slices.append({'kind': 'convex' if take_convex else 'centroid',
                               'lower': [str(q) for q in aa], 'upper': [str(q) for q in bb],
                               'midpoint_centroid': str(at(centroid, mid)),
                               'midpoint_convex': str(at(convex, mid))})
                poly = plus(poly, scale(selected, 6))
        else:
            slices.append({'kind': kind, 'lower': [str(q) for q in lower],
                           'upper': [str(q) for q in right]})
            poly = plus(poly, scale(slice_integral(kind, h, B, C, lower, right), 6))
    assert len(choices) in (0, 2)
    code = str(2 * choices[0] + choices[1]) if choices else '-'
    metadata.update({'tail_lower': [str(q) for q in lower_tail], 'slices': slices})
    return poly, code, metadata


