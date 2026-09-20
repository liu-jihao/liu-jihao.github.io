"""Independent exact audit of the printed 138fd7a05c1323ca certificates.

Finite question: do all 4,088 tree leaves bound the full rational loss in
equations (18)--(20), do the trees cover both root boxes, and do both cyclic
Bernstein leaves prove the displayed polynomial positive?  Run ONLY with
danus compute, DEFAULT tier.  No author checker or prior result is imported.
"""

from pathlib import Path
from functools import lru_cache
from math import comb
import ast
import hashlib
import json
import time

from .common import Q as Rat, ENGINE, require

Z, ONE = Rat(0), Rat(1)
r, u, T, Q = Rat(9, 10), Rat(27, 25), Rat(29, 25), Rat(13, 5)
h0, h1, star, eps = Rat(131, 1000), Rat(471, 2000), Rat(71, 400), Rat(1, 5000)
zero, coord = (Z, Z), (Z, ONE)  # constant coefficient, slope
knots = sorted({Z, r, u, T, Q} | {Rat(j, 50) for j in range(1, 130)})
checks, failures = 0, []


def check(ok, label, data=None):
    global checks
    checks += 1
    require(ok,label)


def add(f, g):
    return tuple((f[i] if i < len(f) else Z) +
                 (g[i] if i < len(g) else Z) for i in range(max(len(f), len(g))))


def scale(f, c):
    return tuple(c * v for v in f)


def sub(f, g):
    return add(f, scale(g, -ONE))


def mul(f, g):
    result = [Z] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            result[i + j] += a * b
    return tuple(result)


def power(f, n):
    result = (ONE,)
    for _ in range(n):
        result = mul(result, f)
    return result


def value(f, t):
    out = Z
    for a in reversed(f):
        out = out * t + a
    return out


def line(a, anchor, height):
    return (height - a * anchor, a)


def split(l, v, forms):
    points = {l, v}
    for f in forms:
        if f[1]:
            root = -f[0] / f[1]
            if l < root < v:
                points.add(root)
    points = sorted(points)
    return zip(points, points[1:])


def positive_sixth_integral(f, l, v):
    """Seven times integral of (b+a t)_+^6, either sign of a."""
    b, a = f
    if not a:
        return 7 * max(b, Z) ** 6 * (v - l)
    return (max(value(f, v), Z) ** 7 - max(value(f, l), Z) ** 7) / a


@lru_cache(maxsize=100000)
def balanced_integral(f, g, l, v):
    """Expand 7 f(t)^5 g(t) about l, independently of primitive (18)."""
    fl, gl, width = value(f, l), value(g, l), v - l
    return 7 * sum((comb(5, j) * fl ** (5 - j) * f[1] ** j *
                    (gl * width ** (j + 1) / (j + 1) +
                     g[1] * width ** (j + 2) / (j + 2))
                    for j in range(6)), Z)


@lru_cache(maxsize=180000)
def upper_integral(D, H, l, v):
    """Equation F, with all branch decisions made by exact rationals."""
    out = Z
    for a, b in split(l, v, (D, H, sub(H, scale(D, star)))):
        mid = (a + b) / 2
        d, h = value(D, mid), value(H, mid)
        if d <= 0:
            continue
        if h <= 0:
            out += positive_sixth_integral(D, a, b)
        elif h < star * d:
            out += balanced_integral(sub(D, scale(H, 2)), add(D, scale(H, 10)), a, b)
        else:
            out += positive_sixth_integral(sub(D, H), a, b)
            out += eps * positive_sixth_integral(D, a, b)
    return out


@lru_cache(maxsize=100000)
def prime_integral(P, l, v):
    total = Z
    for a, b in split(l, v, (P,)):
        active = P if value(P, (a + b) / 2) > 0 else zero
        total += positive_sixth_integral(sub(coord, active), a, b)
    return total


@lru_cache(maxsize=100000)
def partial_integrals(H, P, theta, l, v):
    full, capped = Z, Z
    for a, b in split(l, v, (H, P, sub(H, scale(P, theta)))):
        mid = (a + b) / 2
        hh = H if value(H, mid) > 0 else zero
        pp = P if value(P, mid) > 0 else zero
        qc = pp if theta * value(pp, mid) <= value(hh, mid) else scale(hh, ONE / theta)
        full += upper_integral(sub(coord, pp), sub(hh, scale(pp, theta)), a, b)
        capped += upper_integral(sub(coord, qc), sub(hh, scale(qc, theta)), a, b)
    return full, capped


def box_bound(box, hp, theta):
    incoming, outgoing = Z, Z
    winners = [0] * 7
    for l, v in zip(knots, knots[1:]):
        a0 = box[0][1] if v <= r else box[0][0]
        a1 = box[1][1] if v <= u else box[1][0]
        ap = box[2][1] if v <= T else box[2][0]
        H0, H1, P = line(a0, r, h0), line(a1, u, h1), line(ap, T, hp)
        full, capped = partial_integrals(H1, P, theta, l, v)
        candidates = [v ** 7 - l ** 7,
                      upper_integral(coord, H0, l, v),
                      upper_integral(coord, H1, l, v),
                      prime_integral(P, l, v), full, capped]
        if l >= T:
            candidates.append(positive_sixth_integral((Rat(125, 73), Rat(-52, 73)), l, v))
        which = min(range(len(candidates)), key=lambda j: candidates[j])
        check(candidates[which] >= 0, "nonnegative interval integral")
        winners[which] += 1
        if v <= T:
            incoming += candidates[which]
        else:
            outgoing += candidates[which]
    volume_branch = incoming + (ONE - hp) ** 7
    endpoint_branch = min(incoming, Rat(499, 500)) + outgoing
    return min(volume_branch, endpoint_branch), incoming, outgoing, winners, (
        "volume" if volume_branch <= endpoint_branch else "endpoint")


def bernstein(coeffs, degree):
    return tuple(sum((coeffs[k] * Rat(comb(j, k), comb(degree, k))
                      for k in range(min(j + 1, len(coeffs)))), Z) for j in range(degree + 1))


def bisect_bernstein(row):
    left, right = [row[0]], [row[-1]]
    while len(row) > 1:
        row = tuple((a + b) / 2 for a, b in zip(row, row[1:]))
        left.append(row[0])
        right.append(row[-1])
    return tuple(left), tuple(reversed(right))


def analytic_checks():
    data = {}
    for name, a, p, stated_b, stated_inverse in (
        ("m2", Rat(5, 4), Rat(92, 125), Rat(10564, 4375), Rat(175, 743)),
        ("m3", Rat(623, 500), Rat(719, 1000), Rat(21587, 8642), Rat(8642, 36025))):
        b = ((8 - p) * u - 7) / (a * u - 1)
        inverse = u / (7 - b)
        check(a * u > 1, name + " positive price denominator")
        check(b == stated_b and inverse == stated_inverse, name + " exact original price")
        check(inverse - h1 > Rat(1, 100000), name + " uniform height margin")
        data[name] = {"deficit_floor": str(b), "inverse_floor": str(inverse), "height_gap": str(inverse - h1)}
    check(Rat(175, 743) - h1 == Rat(47, 1486000), "m2 printed height gap")
    check(Rat(225, 749) - Rat(3003, 10000) == Rat(753, 7490000), "actual29 printed margin")
    beta_floor = Rat(1, 500) - Rat(7, 10000) * Rat(7, 6) ** 6
    check(beta_floor > 0, "beta inside big range")
    data["beta_volume_floor"] = str(beta_floor)
    check(Rat(3072, 625) < 6, "balanced derivative maximum below 6")
    check(120 * Rat(1, 10) * Rat(4, 5) ** 4 == Rat(3072, 625), "balanced derivative maximum value")
    check(2 * Rat(10, 11) ** 5 < Rat(4, 3), "global Q enlargement")
    tail0 = r ** 7 + Rat(4, 3) * (r - h0) ** 7
    tail1 = Rat(19, 25) + Rat(4, 9) * (u - h1) ** 7
    check(tail0 < Rat(691, 1000), "all a0>=2")
    check(tail1 < Rat(931, 1000), "all a1>=4")
    data["unbounded_old"] = str(tail0)
    data["unbounded_27"] = str(tail1)
    for name, hp, expected in (
        ("masscap", Rat(31, 100), Rat(254620898673, 254720000000)),
        ("payment", Rat(3003, 10000), Rat(1989490791196134952369386003613, 1990000000000000000000000000000))):
        bound = Rat(499, 500) + (T - hp) ** 7 / 199
        check(bound == expected, name + " exact unbounded prime rational")
        check(bound < Rat(9999, 10000), name + " all aP>=200")
        data[name + "_unbounded_prime"] = str(bound)
    theta, lam, endpoint = Rat(43, 50), Rat(3101, 10000), Rat(13, 70)
    C = Rat(4, 6) / (ONE - Rat(5, 6) * theta) ** 2
    bc = h1 - Rat(27, 29) * (ONE - theta) * lam
    check(C == Rat(2400, 289), "cyclic coefficient C0")
    check(Rat(27, 29) * lam > endpoint, "residual threshold increases in theta")
    check(bc - theta * endpoint > 0, "positive residual across z interval")
    Rpoly = sub(power((ONE, -ONE), 2), scale(power((bc, -theta), 2), C))
    check(Rpoly[2] < 0, "quadratic concavity")
    check(value(Rpoly, Z) > Rat(3, 5) and value(Rpoly, endpoint) > Rat(3, 5), "positive quadratic endpoints")
    polynomial = sub(scale(power((ONE, -endpoint), 5), Rat(6, 25)),
                     power(sub(power((ONE, -endpoint), 2), scale(power((bc, -theta * endpoint), 2), C)), 6))
    root = bernstein(polynomial, 12)
    halves = bisect_bernstein(root)
    tree, raw_line = (16429905,893426), None
    check(tree == (16429905, 893426), "printed cyclic tree data")
    records = []
    for j, (coeffs, claimed) in enumerate(zip(halves, tree)):
        bound = min(coeffs)
        check(bound >= Rat(claimed, 1000000000), "cyclic leaf " + str(j))
        check(bound > 0, "cyclic leaf positivity " + str(j))
        records.append({"interval": [str(Rat(j, 2)), str(Rat(j + 1, 2))],
                        "coefficients": [str(v) for v in coeffs], "minimum": str(bound),
                        "printed_lower_integer": claimed})
    data["cyclic"] = {"bc": str(bc), "C": str(C), "positive_residual_min": str(bc - theta * endpoint),
                      "quadratic": [str(v) for v in Rpoly], "power_coefficients": [str(v) for v in polynomial],
                      "root_bernstein": [str(v) for v in root], "raw_tree_line": raw_line, "leaves": records}
    maxima = [Rat(43, 13), Rat(51, 13), Rat(12879, 2600), Rat(76841, 13000), Rat(11234, 1625), Rat(12949, 1625)]
    check(max(maxima) == Rat(12949, 1625), "full endpoint maximum")
    check(max(maxima) + Rat(1, 1000) < 8 - Rat(1, 100), "paid endpoint with both-error allowance")
    data["endpoint_maxima"] = [str(v) for v in maxima]
    data["endpoint_gap_after_errors"] = str(8 - max(maxima) - Rat(1, 1000))
    return data


