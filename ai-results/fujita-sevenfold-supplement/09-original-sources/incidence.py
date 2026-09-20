"""Exact independent 938 incidence integrals, adapted from sub11 audit.
The unrelated 718 family is outside this package and is not read.
"""
from fractions import Fraction as Q
from collections import Counter

R, U, B = Q(9, 10), Q(26, 25), Q(241, 100)
ZERO, T = (Q(0), Q(0)), (Q(1), Q(0))  # slope, intercept
checks = Counter()


def check(truth, category):
    checks[category] += 1
    if not truth:
        raise AssertionError(category)


def at(line, t):
    return line[0] * t + line[1]


def linear_combination(*terms):
    return tuple(sum(k * f[i] for k, f in terms) for i in (0, 1))


def anchored(slope, anchor, height):
    return (slope, height - slope * anchor)


def split(l, u, lines):
    points = {l, u}
    for slope, intercept in lines:
        if slope:
            root = -intercept / slope
            if l < root < u:
                points.add(root)
    points = sorted(points)
    return list(zip(points, points[1:]))


def sign(line, l, u):
    a, b = at(line, l), at(line, u)
    check(a * b >= 0, 'affine_sign_constant')
    return 1 if a + b > 0 else -1 if a + b < 0 else 0


def interpolation_weights():
    # Integrate the seven Lagrange basis polynomials on [0,1].
    nodes = [Q(j, 6) for j in range(7)]
    weights = []
    for j, x in enumerate(nodes):
        polynomial = [Q(1)]
        for m, y in enumerate(nodes):
            if m == j:
                continue
            updated = [Q(0)] * (len(polynomial) + 1)
            for i, coefficient in enumerate(polynomial):
                updated[i] -= y * coefficient / (x - y)
                updated[i + 1] += coefficient / (x - y)
            polynomial = updated
        weights.append(sum(v / (i + 1) for i, v in enumerate(polynomial)))
    for power in range(7):
        check(sum(w * x**power for w, x in zip(weights, nodes)) == Q(1, power + 1),
              'exact_interpolation_moments')
    return nodes, weights


NODES, WEIGHTS = interpolation_weights()


def integrate(function, l, u):
    return (u - l) * sum(w * function(l + x * (u - l))
                         for x, w in zip(NODES, WEIGHTS))


def rate(t):
    return (Q(1) if t < R else Q(61, 100) if t < Q(99, 100)
            else Q(11, 25) if t < U else Q(17, 50) if t < Q(29, 25)
            else Q(1))


def ceiling_million(x):
    return -((-x.numerator * 10**6) // x.denominator)


def certificate(a0, a1, d0, d1, end, height, theta):
    total, records = Q(0), []
    anchors = sorted({q for q in (Q(0), R, Q(99, 100), U, Q(29, 25), end)
                      if 0 <= q <= end})
    for left, right in zip(anchors, anchors[1:]):
        hline = anchored(a1 if right <= R else a0, R, height)
        pline = anchored(d1 if right <= U else d0, U, Q(1, 5))
        for l1, u1 in split(left, right, (hline, pline)):
            h = hline if sign(hline, l1, u1) > 0 else ZERO
            p = pline if sign(pline, l1, u1) > 0 else ZERO
            qstar = linear_combination((7 / (7 * theta - 1), h),
                                      (-1 / (7 * theta - 1), T))
            qminus = linear_combination((1, qstar), (-1, p))
            for l2, u2 in split(l1, u1, (qstar, qminus)):
                q = (ZERO if sign(qstar, l2, u2) <= 0
                     else p if sign(qminus, l2, u2) >= 0 else qstar)
                for exponent in (ZERO, p, q):
                    check(min(at(exponent, l2), at(exponent, u2)) >= 0,
                          'legal_nonnegative_partial_exponent')
                    check(min(at(linear_combination((1, p), (-1, exponent)), l2),
                              at(linear_combination((1, p), (-1, exponent)), u2)) >= 0,
                          'partial_exponent_below_same_prime')
                pairs = [(linear_combination((1, T), (-1, exponent)),
                          linear_combination((1, h), (-theta, exponent)))
                         for exponent in (ZERO, p, q)]
                boundaries = []
                for x, y in pairs:
                    boundaries.extend((x, y, linear_combination((1, y), (-1, x)),
                                       linear_combination((7, y), (-1, x)),
                                       linear_combination((2, y), (-1, x))))
                for l3, u3 in split(l2, u2, boundaries):
                    shell_rate = rate((l3 + u3) / 2)
                    bounds = [('shell', lambda t, s=shell_rate: s * t**6)]
                    for j, (x, y) in enumerate(pairs):
                        syx = sign(linear_combination((1, y), (-1, x)), l3, u3)
                        sx, sy = sign(x, l3, u3), sign(y, l3, u3)
                        s7 = sign(linear_combination((7, y), (-1, x)), l3, u3)
                        s2 = sign(linear_combination((2, y), (-1, x)), l3, u3)
                        if sx <= 0 or syx >= 0:
                            bounds.append((f'q{j}:zero', lambda t: Q(0)))
                        elif sy <= 0:
                            bounds.append((f'q{j}:ordinary', lambda t, x=x: at(x, t)**6))
                        else:
                            if s7 <= 0:
                                bounds.append((f'q{j}:lower_chord', lambda t, x=x, y=y:
                                               at(x, t)**6 - Q(35, 6)**5 * at(x, t) * at(y, t)**5))
                            else:
                                bounds.append((f'q{j}:upper_chord', lambda t, x=x, y=y:
                                               at(x, t)**6 - at(x, t) * ((29 * at(x, t) + 7 * at(y, t)) / 36)**5))
                                bounds.append((f'q{j}:convex_order', lambda t, x=x, y=y:
                                               Q(7, 6)**6 * (at(x, t) - at(y, t))**6))
                            if s2 >= 0:
                                bounds.append((f'q{j}:axis', lambda t, x=x, y=y:
                                               (at(x, t) - at(y, t))**6))
                    midpoint = (l3 + u3) / 2
                    for lo, hi in ((l3, midpoint), (midpoint, u3)):
                        values = [(label, integrate(fun, lo, hi)) for label, fun in bounds]
                        check(all(value >= 0 for _, value in values), 'nonnegative_integrals')
                        winner, value = min(values, key=lambda v: v[1])
                        total += 7 * value
                        records.append({'interval': [str(lo), str(hi)], 'winner': winner,
                                        'all_integrals': [[label, str(v)] for label, v in values]})
    check(records[0]['interval'][0] == '0' and records[-1]['interval'][1] == str(end),
          'entire_order_interval')
    for before, after in zip(records, records[1:]):
        check(before['interval'][1] == after['interval'][0], 'adjacent_order_cells')
    return total, records


def coverage(incoming, outgoing, amin):
    check(incoming[0][0] == amin and incoming[-1][1] == 3, 'incoming_endpoints')
    for i, row in enumerate(incoming):
        check(row[0] < row[1], 'incoming_positive_width')
        if i:
            check(incoming[i - 1][1] == row[0], 'incoming_no_gap_or_overlap')
    # The full arrangement of all rectangle edges, rather than an area sum.
    aa = sorted({v for row in outgoing for v in row[:2]})
    dd = sorted({v for row in outgoing for v in row[2:4]})
    check(aa[0] == amin and aa[-1] == 3 and dd[0] == B and dd[-1] == 5,
          'outgoing_domain_endpoints')
    for left, right in zip(aa, aa[1:]):
        for low, high in zip(dd, dd[1:]):
            owners = [row for row in outgoing
                      if row[0] <= left < right <= row[1] and row[2] <= low < high <= row[3]]
            check(len(owners) == 1, 'outgoing_exact_atomic_rectangle_cover')

