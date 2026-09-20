#!/usr/bin/env python3
"""Independent exact checker for the rank-two finite certificate.

Adapted from the separately written sub27_e13_mesh audit, not the producer.
Run: python3 check.py --data data --report report.json
All mathematical conditions use explicit exceptions, including under python -O.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from pathlib import Path
import json
import sys


class CertificateError(ValueError):
    """A malformed payload or a failed mathematical condition."""


def require(condition, message):
    if not condition:
        raise CertificateError(message)


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f'duplicate JSON key: {key}')
        result[key] = value
    return result


def reject_number(value):
    raise CertificateError(f'nonintegral JSON number is forbidden: {value}')


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'),
                      object_pairs_hook=unique_object,
                      parse_float=reject_number, parse_constant=reject_number)


def exact_keys(row, keys, where):
    require(type(row) is dict, f'{where}: expected an object')
    require(set(row) == set(keys), f'{where}: missing or extra fields')


def read_payload(data_dir):
    expected = {'state_rows.json', 'mesh_rows.json'}
    require(data_dir.is_dir(), 'data directory is missing')
    require({p.name for p in data_dir.iterdir()} == expected,
            'data directory must contain exactly state_rows.json and mesh_rows.json')
    for name in expected:
        require((data_dir / name).is_file(), f'{name}: expected a regular file')
    states = read_json(data_dir / 'state_rows.json')
    mesh = read_json(data_dir / 'mesh_rows.json')
    require(type(states) is list and len(states) == 28, 'states: expected exactly 28 records')
    require(type(mesh) is list and len(mesh) == 235, 'mesh: expected exactly 235 records')
    state_keys = {'state', 'T', 'intrinsic_cap', 'root_cap', 'radius',
                  'lambda_price', 'theta_price', 'theta_kind', 'descendants'}
    mesh_keys = {'j', 'gamma', 'Y', 'expressions', 'L', 'E',
                 'maximizers', 'rounded_prefix_reserve'}
    previous = None
    for i, row in enumerate(states):
        exact_keys(row, state_keys, f'states[{i}]')
        v = row['state']
        require(type(v) is list and len(v) == 3 and all(type(x) is int for x in v),
                f'states[{i}]: state must be a triple of integers')
        require(previous is None or previous < v, 'states: duplicate or unordered state')
        previous = v
    for i, row in enumerate(mesh):
        exact_keys(row, mesh_keys, f'mesh[{i}]')
        require(type(row['j']) is int and row['j'] == 1731 + i,
                f'mesh[{i}]: missing, duplicate, unordered or out-of-range index')
        require(type(row['E']) is int and row['E'] >= 0,
                f'mesh[{i}]: E must be a nonnegative integer')
    return states, mesh


def compare(actual, expected, where):
    """Consume every supplied field; rational strings must equal reduced values."""
    require(type(actual) is type(expected), f'{where}: wrong JSON type')
    if type(expected) is dict:
        require(set(actual) == set(expected), f'{where}: missing or extra fields/candidates')
        for key in expected:
            compare(actual[key], expected[key], f'{where}.{key}')
    elif type(expected) is list:
        require(len(actual) == len(expected), f'{where}: wrong list length')
        for i, (left, right) in enumerate(zip(actual, expected)):
            compare(left, right, f'{where}[{i}]')
    else:
        require(actual == expected, f'{where}: supplied witness differs from exact value')


def validate(data_dir):
    state_payload, mesh_payload = read_payload(data_dir)
    B, r, lam, theta, step = Q(18, 5), Q(25, 27), Q(1731, 2000), Q(983, 1000), Q(1, 2000)
    zero = Q(0)
    radius_lists = {
        1: ['1'],
        2: ['1', '283/200', '1733/1000', '2', '2237/1000'],
        3: ['1', '63/50', '1443/1000', '397/250', '171/100', '909/500', '1913/1000', '2', '2081/1000'],
        4: ['1', '119/100', '1317/1000', '283/200', '187/125', '783/500', '1627/1000'],
    }
    qv, moment, cap, rootcap = {}, {}, {}, {}
    
    def first_moment(c, m):
        if c == 0:
            require(m == 1, 'arithmetic condition: m == 1')
            return 0
        remaining, j, total = m, 0, 0
        while remaining:
            take = min(remaining, comb(c + j - 1, j))
            total += j * take
            remaining -= take
            j += 1
        return total
    
    # Exhaustive bounds come from the positive-deficit first-moment inequality.
    for d in range(1, 5):
        for e in range(d, 7):
            if d == 1 and e != 1:
                continue
            c = e - d
            ell = d // 2 + 1
            max_m = (2 * ell * comb(c + ell - 1, ell - 1) - 1) // (2 * ell - d) if c else 1
            for m in range(c + 1 if c else 1, max_m + 1):
                T = first_moment(c, m)
                intrinsic = d - Q(2 * T, m)
                C = min(intrinsic, Q(8 - m, 2)) if d == 4 else intrinsic
                if C <= 0:
                    continue
                v = (d, e, m)
                moment[v], cap[v], rootcap[v] = T, intrinsic, C
                qv[v] = Q(radius_lists[d][m - 1])
                require(qv[v] ** d >= m, 'arithmetic condition: qv[v] ** d >= m')
    
    states = sorted(qv)
    require(len(states) == 28, 'arithmetic condition: len(states) == 28')
    expected = {(1, 1, 1)} | {(2, m + 1, m) for m in range(1, 6)}
    expected |= {(3, 3, 1), (4, 4, 1)}
    for e, ms in [(4, range(2, 4)), (5, range(3, 7)), (6, range(4, 10))]:
        expected |= {(3, e, m) for m in ms}
    for e, ms in [(5, range(2, 5)), (6, range(3, 8))]:
        expected |= {(4, e, m) for m in ms}
    require(set(states) == expected, 'arithmetic condition: set(states) == expected')
    edges = {}
    for v in states:
        d, e, m = v
        edges[v] = []
        for w in states:
            f, h, mw = w
            if f >= d or h > e:
                continue
            C = cap[w]
            if f == d - 1:
                C = min(C, Q(d * m - 2 * moment[v] - mw, m))
            if C > 0:
                edges[v].append((w, C))
    
    @lru_cache(None)
    def J(v, k, C):
        require(0 <= C <= cap[v], 'arithmetic condition: 0 <= C <= cap[v]')
        q = qv[v]
        return max(q - k, zero) * C + max(
            [zero] + [J(w, max(k, q), min(C, bound)) for w, bound in edges[v]]
        )
    
    def kappa(gamma):
        return (7 - B) * r + B / gamma
    
    def price(v, gamma, C=None):
        if C is None:
            C = rootcap[v]
        return kappa(gamma) + J(v, 1 / gamma, min(B, C))
    
    quartic_child = (3, 6, 4)
    cubic = (4, 6, 3)
    exceptional = {cubic, (4, 6, 4), (4, 6, 5)}
    lambda_prices = {v: price(v, lam, min(rootcap[v], 1) if v == quartic_child else rootcap[v]) for v in states}
    paid_lambda = {v: p for v, p in lambda_prices.items() if v not in exceptional}
    require(max(paid_lambda.values()) == Q(155528897, 19473750) < Q(799, 100), 'arithmetic condition: max(paid_lambda.values()) == Q(155528897, 19473750) < Q(799, 100)')
    theta_prices = {v: price(v, theta) for v in states if v != cubic}
    require(max(theta_prices.values()) == Q(211925069, 26541000) < Q(7998, 1000), 'arithmetic condition: max(theta_prices.values()) == Q(211925069, 26541000) < Q(7998, 1000)')
    require(kappa(theta) < Q(7998, 1000), 'arithmetic condition: kappa(theta) < Q(7998, 1000)')  # Direct point outcome.
    
    q = {3: Q(1317, 1000), 4: Q(283, 200), 5: Q(187, 125)}
    p = {4: Q(379, 1000), 5: Q(741, 2500)}
    for m in (4, 5):
        v = (4, 6, m)
        require(J(v, q[m], cap[v]) == p[m], 'arithmetic condition: J(v, q[m], cap[v]) == p[m]')
    
    threshold = {3: Q(751, 100), 4: 8 - Q(1, 500) - p[4], 5: 8 - Q(1, 500) - p[5]}
    A = {m: q[m] - r for m in q}
    C = {m: threshold[m] - 7 * r for m in q}
    D = {m: A[m] * B - C[m] for m in q}
    require([A[m] for m in q] == [Q(10559, 27000), Q(2641, 5400), Q(1924, 3375)], 'arithmetic condition: [A[m] for m in q] == [Q(10559, 27000), Q(2641, 5400), Q(1924, 3375)]')
    require([C[m] for m in q] == [Q(2777, 2700), Q(30713, 27000), Q(41179, 33750)], 'arithmetic condition: [C[m] for m in q] == [Q(2777, 2700), Q(30713, 27000), Q(41179, 33750)]')
    require([D[m] for m in q] == [Q(12803, 33750), Q(673, 1080), Q(5617, 6750)], 'arithmetic condition: [D[m] for m in q] == [Q(12803, 33750), Q(673, 1080), Q(5617, 6750)]')
    
    def R(m, gamma):
        return (q[m] * gamma - 1) / D[m]
    
    def Y(m, gamma):
        return B - gamma / R(m, gamma)
    
    require(C[3] > Q(7, 3) * A[3], 'arithmetic condition: C[3] > Q(7, 3) * A[3]')
    require(C[4] > 2 * A[4], 'arithmetic condition: C[4] > 2 * A[4]')
    require(C[5] > Q(3, 2) * A[5], 'arithmetic condition: C[5] > Q(3, 2) * A[5]')
    for m, factor in [(3, Q(19, 25)), (4, Q(4, 5)), (5, Q(21, 25))]:
        require(factor * R(m, lam) > Q(1, 4), 'arithmetic condition: factor * R(m, lam) > Q(1, 4)')
    require(Y(4, lam) > 1 and Y(5, lam) > 1, 'arithmetic condition: Y(4, lam) > 1 and Y(5, lam) > 1')
    require(Y(3, theta) > Q(7, 3) and Y(4, theta) > 2 and Y(5, theta) > Q(3, 2), 'arithmetic condition: Y(3, theta) > Q(7, 3) and Y(4, theta) > 2 and (Y(5, theta) > Q(3, 2))')
    cubic_theta = kappa(theta) + (q[3] - 1 / theta) * Q(7, 3)
    require(cubic_theta == Q(199315493, 26541000) < Q(751, 100), 'arithmetic condition: cubic_theta == Q(199315493, 26541000) < Q(751, 100)')
    
    # The derivative conditions used to cover continuous germ parameters.
    require(4 * Q(4, 5) * R(3, lam) - lam > step, 'arithmetic condition: 4 * Q(4, 5) * R(3, lam) - lam > step')
    require(4 * Q(21, 25) * R(4, lam) - lam > step, 'arithmetic condition: 4 * Q(21, 25) * R(4, lam) - lam > step')
    require(Q(5, 2) * R(3, lam) - lam > step, 'arithmetic condition: Q(5, 2) * R(3, lam) - lam > step')
    require(4 * Q(4, 5) * q[3] / D[3] > 1, 'arithmetic condition: 4 * Q(4, 5) * q[3] / D[3] > 1')
    require(4 * Q(21, 25) * q[4] / D[4] > 1, 'arithmetic condition: 4 * Q(21, 25) * q[4] / D[4] > 1')
    require(Q(5, 2) * q[3] / D[3] > 1, 'arithmetic condition: Q(5, 2) * q[3] / D[3] > 1')
    x = Q(1, 3)
    require(110 - 440*x + Q(745, 2)*x*x - 75*x*x*x == Q(35, 18), 'arithmetic condition: 110 - 440 * x + Q(745, 2) * x * x - 75 * x * x * x == Q(35, 18)')
    require(-440 + 745*x < 0, 'arithmetic condition: -440 + 745 * x < 0')
    require(4*B - 12 > 0, 'arithmetic condition: 4 * B - 12 > 0')
    require(5*(B + 1) - 10*Q(5, 2) + Q(5, 2)**2 == Q(17, 4), 'arithmetic condition: 5 * (B + 1) - 10 * Q(5, 2) + Q(5, 2) ** 2 == Q(17, 4)')
    
    # Exact complete price interfaces, including all cubic descendants.
    a, affine_epsilon = Q(5000, 3997), Q(49267, 107919)
    direct_affine = kappa(lam) + (a - 1/lam) * Q(3, 2)
    require(direct_affine == Q(463961965, 62269263), 'arithmetic condition: direct_affine == Q(463961965, 62269263)')
    require(7 + affine_epsilon - direct_affine == Q(349935, 62269263) > Q(1, 1000), 'arithmetic condition: 7 + affine_epsilon - direct_affine == Q(349935, 62269263) > Q(1, 1000)')
    cubic_corrections = {w: J(w, q[3], bound) for w, bound in edges[cubic] if w != quartic_child}
    require(max(cubic_corrections.values()) == Q(451, 1000) <= Q(977, 2000), 'arithmetic condition: max(cubic_corrections.values()) == Q(451, 1000) <= Q(977, 2000)')
    require(J(quartic_child, q[3], Q(1)) == Q(477, 1000), 'arithmetic condition: J(quartic_child, q[3], Q(1)) == Q(477, 1000)')
    require(affine_epsilon - Q(51, 100) - a + q[3] == Q(1357633, 107919000) > Q(1, 2000), 'arithmetic condition: affine_epsilon - Q(51, 100) - a + q[3] == Q(1357633, 107919000) > Q(1, 2000)')
    require(1 - Q(51, 100) - Q(977, 2000) - Q(1, 2000) > Q(1, 400000), 'arithmetic condition: 1 - Q(51, 100) - Q(977, 2000) - Q(1, 2000) > Q(1, 400000)')
    require(Q(1, 500) - Q(1, 1000) > Q(1, 400000), 'arithmetic condition: Q(1, 500) - Q(1, 1000) > Q(1, 400000)')
    require(1 - r == Q(2, 27), 'arithmetic condition: 1 - r == Q(2, 27)')
    
    # The supplied indices were checked before evaluating any arithmetic.
    certificate = {row['j']: row['E'] for row in mesh_payload}
    require(sum(certificate.values()) == 28148081, 'mesh: sum of delivered E_j')
    
    def efunc(x):
        return (8 - 5*x) / (10 - 5*x)
    
    def ffunc(x):
        return (1 + x) * efunc(x)
    
    def loss_expressions(gamma, upper):
        def Phi(z):
            return max(upper - z, zero)**5 - max(gamma - z, zero)**5
        rr, yy = {m: R(m, gamma) for m in q}, {m: Y(m, gamma) for m in q}
        for m in q:
            require(D[m] > 0 and q[m]*gamma - 1 > 0, 'arithmetic condition: D[m] > 0 and q[m] * gamma - 1 > 0')
            require(rr[m] * (B - yy[m]) == gamma, 'arithmetic condition: rr[m] * (B - yy[m]) == gamma')
            require(yy[m] == (threshold[m] - kappa(gamma)) / (q[m] - 1/gamma), 'arithmetic condition: yy[m] == (threshold[m] - kappa(gamma)) / (q[m] - 1 / gamma)')
        require(upper <= 4*Q(4, 5)*rr[3], 'arithmetic condition: upper <= 4 * Q(4, 5) * rr[3]')
        require(upper <= 4*Q(21, 25)*rr[4], 'arithmetic condition: upper <= 4 * Q(21, 25) * rr[4]')
        require(upper <= Q(5, 2)*rr[3], 'arithmetic condition: upper <= Q(5, 2) * rr[3]')
        expr = {}
        if yy[3] < Q(7, 3):
            xstar = max(zero, yy[3] - 2)
            require(0 <= xstar < Q(1, 3), 'arithmetic condition: 0 <= xstar < Q(1, 3)')
            expr['smooth_cubic_extreme'] = Phi(rr[3]*ffunc(xstar))/(1-xstar) + (1-2*xstar)/(1-xstar)*Phi(Q(38, 25)*rr[3])
            expr['smooth_cubic_balanced'] = 2*Phi(Q(57, 50)*rr[3])
        if yy[3] < 2:
            expr['other_cubic'] = Phi(Q(4, 5)*rr[3]) + Phi(Q(8, 5)*rr[3])
        if yy[4] < Q(3, 2):
            expr['distinct_quartic'] = Phi(Q(21, 25)*rr[4]) + Phi(Q(63, 25)*rr[4])
        if yy[4] < 2:
            expr['cartier_quartic'] = 2*Phi(Q(8, 5)*rr[4])
            expr['power_quartic'] = 2*Phi(gamma*yy[4]/(4-yy[4]))
        if yy[5] < Q(3, 2):
            expr['cartier_quintic'] = 2*Phi(Q(42, 25)*rr[5])
            expr['power_quintic'] = 2*Phi(gamma*(1+yy[5])/(4-yy[5]))
        require(all(val >= 0 for val in expr.values()), 'arithmetic condition: all((val >= 0 for val in expr.values()))')
        return expr, yy
    
    def ceilq(x):
        return (x.numerator + x.denominator - 1) // x.denominator
    
    M0 = 1 - 2*lam**5
    require(M0 == Q(458771731542349, 16000000000000000), 'arithmetic condition: M0 == Q(458771731542349, 16000000000000000)')
    rows, prefix_spent, ceiling_mismatches, active_ranges, maximizers = [], 0, [], {}, []
    for j in range(1731, 1966):
        gamma, upper = Q(j, 2000), Q(j+1, 2000)
        expr, yy = loss_expressions(gamma, upper)
        L = max([zero] + list(expr.values()))
        E = certificate[j]
        for name, value in expr.items():
            require(value * 10**9 <= E, 'arithmetic condition: value * 10 ** 9 <= E')
            active_ranges.setdefault(name, []).append(j)
        require(L * 10**9 <= E, 'arithmetic condition: L * 10 ** 9 <= E')
        require(E - 1 < L * 10**9 < E, f'mesh[{j}]: strict ceiling bracket')
        if ceilq(L * 10**9) != E:
            ceiling_mismatches.append(j)
        prefix_spent += E
        reserve = M0 - Q(prefix_spent, 10**9)
        require(reserve > Q(1, 2000), 'arithmetic condition: reserve > Q(1, 2000)')
        require(reserve - Q(1, 4000) > Q(1, 4000), 'arithmetic condition: reserve - Q(1, 4000) > Q(1, 4000)')
        winning = sorted(name for name, value in expr.items() if value == L)
        maximizers.append((j, winning))
        rows.append({'j': j, 'gamma': str(gamma), 'Y': {str(m): str(yy[m]) for m in q},
                     'expressions': {name: str(value) for name, value in expr.items()},
                     'L': str(L), 'E': E, 'maximizers': winning, 'rounded_prefix_reserve': str(reserve)})
    require(not ceiling_mismatches, 'arithmetic condition: not ceiling_mismatches')
    final_reserve = M0 - Q(prefix_spent, 10**9)
    require(final_reserve == Q(8402435542349, 16000000000000000), 'arithmetic condition: final_reserve == Q(8402435542349, 16000000000000000)')
    require(final_reserve - Q(1, 2000) == Q(402435542349, 16000000000000000), 'arithmetic condition: final_reserve - Q(1, 2000) == Q(402435542349, 16000000000000000)')
    
    def group_winners(items):
        groups = []
        for j, names in items:
            if groups and groups[-1]['names'] == names and groups[-1]['end'] + 1 == j:
                groups[-1]['end'] = j
            else:
                groups.append({'start': j, 'end': j, 'names': names})
        return groups
    
    # All continuous extinction thresholds: solve Y_m(gamma)=cap exactly.
    extinction = {}
    for m, caps in [(3, [Q(2), Q(7, 3)]), (4, [Q(3, 2), Q(2)]), (5, [Q(3, 2)])]:
        for bound in caps:
            g = (B-bound)/(q[m]*(B-bound)-D[m])
            require(Y(m, g) == bound, 'arithmetic condition: Y(m, g) == bound')
            extinction[f'm{m}_cap{bound}'] = str(g)
    
    summary = {
        'all_checks_pass': True,
        'scope': 'Finite source/whole-price arithmetic; local geometric premises imported only at read scopes.',
        'state_count': len(states),
        'allowed_state_edges': sum(len(v) for v in edges.values()),
        'mesh_points': 236, 'mesh_intervals': len(rows),
        'all_printed_E_are_exact_ceilings': not ceiling_mismatches,
        'sum_E': prefix_spent,
        'initial_mass': str(M0), 'final_rounded_reserve': str(final_reserve),
        'margin_over_1_over_2000': str(final_reserve-Q(1, 2000)),
        'lambda_paid_max': str(max(paid_lambda.values())),
        'lambda_paid_max_state': [list(v) for v, val in paid_lambda.items() if val == max(paid_lambda.values())],
        'theta_noncubic_max': str(max(theta_prices.values())),
        'theta_noncubic_max_state': [list(v) for v, val in theta_prices.items() if val == max(theta_prices.values())],
        'theta_noncubic_gap_to_8': str(8-max(theta_prices.values())),
        'theta_cubic_affine_max': str(cubic_theta),
        'theta_cubic_entry_slack': str(Q(751, 100)-cubic_theta),
        'direct_quartic_affine_max': str(direct_affine),
        'direct_quartic_entry_slack': str(7+affine_epsilon-direct_affine),
        'Y_lambda': {str(m): str(Y(m, lam)) for m in q},
        'Y_theta': {str(m): str(Y(m, theta)) for m in q},
        'minimum_fixed_coefficients_at_lambda': {str(m): str(factor*R(m,lam)) for m, factor in [(3,Q(19,25)),(4,Q(4,5)),(5,Q(21,25))]},
        'active_family_index_ranges': {name: [min(js), max(js)] for name, js in active_ranges.items()},
        'maximizer_index_ranges': group_winners(maximizers),
        'continuous_extinction_gammas': extinction,
        'cubic_nonquartic_max_correction': str(max(cubic_corrections.values())),
    }
    
    state_rows = []
    for v in states:
        state_rows.append({'state': list(v), 'T': moment[v], 'intrinsic_cap': str(cap[v]), 'root_cap': str(rootcap[v]),
                           'radius': str(qv[v]), 'lambda_price': str(lambda_prices[v]),
                           'theta_price': str(theta_prices[v]) if v != cubic else str(cubic_theta),
                           'theta_kind': 'cubic affine endpoint' if v == cubic else 'complete scalar tail',
                           'descendants': [{'state': list(w), 'cap': str(bound)} for w, bound in edges[v]]})
    
    compare(mesh_payload, rows, 'mesh')
    compare(state_payload, state_rows, 'states')
    summary['schema'] = 'fujita-rank-two-check-v1'
    summary['status'] = 'pass'
    summary['families'] = {
        '20-rank-two-feedback/states': {
            'records': len(state_rows),
            'successor_records': sum(len(row['descendants']) for row in state_rows),
            'lambda_prices': len(state_rows),
            'theta_prices': len(state_rows),
            'interface': 'Complete positive-cap census, 170 successor caps, moments, radii, full recurrence, both whole-price columns; theta cubic is affine.'
        },
        '20-rank-two-feedback/mesh': {
            'records': len(rows), 'endpoints': len(rows) + 1,
            'first_index': 1731, 'last_index': 1965,
            'candidate_formulas': 8,
            'enabled_candidate_records': sum(len(row['expressions']) for row in rows),
            'strict_ceiling_brackets': len(rows),
            'prefix_reserves': len(rows),
            'interface': 'Exact active-key set, every candidate value, maximum and all ties, E-1 < 10^9 L < E, every reserve and terminal mass.'
        }
    }
    summary['rounding'] = 'Ordinary ceiling is the rule; all 235 supplied maxima are nonintegral after scaling, so floor-plus-one agrees here.'
    summary['geometric_scope'] = 'Conditional on the full symbolic/principal estimates, continuous reductions and complete-source injection proved in the article; no geometric proof is executed.'
    summary['dependencies'] = 'Python >=3.9 standard library; only the two declared JSON inputs.'
    return summary



def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args(argv)
    report_path = args.report.resolve()
    data_path = args.data.resolve()
    if report_path == Path(__file__).resolve() or data_path == report_path or data_path in report_path.parents:
        parser.error('--report must not overwrite the checker or its data')
    try:
        report = validate(data_path)
        exit_code = 0
    except (CertificateError, OSError, ValueError, TypeError, KeyError, IndexError, ZeroDivisionError) as error:
        report = {'schema': 'fujita-rank-two-check-v1', 'status': 'fail',
                  'all_checks_pass': False, 'error': str(error)}
        exit_code = 1
    try:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    except OSError as error:
        print(f'cannot write report: {error}', file=sys.stderr)
        return 1
    print(json.dumps({'status': report['status'],
                      'report': str(args.report),
                      **({'error': report['error']} if exit_code else {'families': report['families']})},
                     sort_keys=True))
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
