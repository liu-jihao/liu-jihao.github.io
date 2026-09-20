#!/usr/bin/env python3
"""Independent exact checker for the portable intrinsic-double certificates.

Python >= 3.9, standard library only. No candidate producer is imported.
Adapted from the archived independent 3672 and 18f checkers; see manifest.
On the Danus host execute only through danus compute DEFAULT.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from math import comb
from pathlib import Path
import re
import sys


class Invalid(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise Invalid(message)


def rational(text):
    require(re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?", text),
            "malformed rational: " + repr(text))
    value = Q(text)
    require(str(value) == text, "noncanonical rational: " + text)
    return value


def integer(text):
    require(re.fullmatch(r"0|[1-9][0-9]*", text),
            "malformed nonnegative integer: " + repr(text))
    return int(text)


SCHEMAS = {
    "states.tsv": ("d e m moment cap radius J_50_49 P_lambda P_theta", 28, 3),
    "transitions.tsv": ("id lo hi shift loss remaining_mass", 8, 1),
    "ambient-prefix.tsv": ("id lo hi numerator", 19, 1),
    "ambient-boxes.tsv": ("id a_lo a_hi d_lo d_hi numerator", 237, 1),
    "rank-two-prices.tsv": ("d e m lambda_b3 theta_b3 lambda_b18_5 J_eta", 28, 3),
    "rank-two-powers.tsv": ("M xi shift", 7, 1),
    "identities.tsv": ("name value", None, 1),
}


def read_table(path, schema):
    header, count, key_size = schema
    raw = path.read_bytes()
    require(len(raw) < 2_000_000, f"{path.name}: oversized payload")
    require(raw.endswith(b"\n") and re.fullmatch(rb"[\t\n\x20-\x7e]*", raw),
            f"{path.name}: expected ASCII TSV with LF-only records")
    lines = raw[:-1].decode("ascii").split("\n")
    fields = header.split()
    require(lines[0] == "\t".join(fields), f"{path.name}: header mismatch")
    require(count is None or len(lines) == count + 1,
            f"{path.name}: expected {count} rows, received {len(lines)-1}")
    rows, keys = [], set()
    for index, line in enumerate(lines[1:], 1):
        parts = line.split("\t")
        require(len(parts) == len(fields), f"{path.name}:{index}: field count")
        values = []
        for field, text in zip(fields, parts):
            if field == "name":
                require(re.fullmatch(r"[a-z][a-z0-9_]*", text),
                        f"{path.name}:{index}: invalid name")
                values.append(text)
            elif field in {"id", "d", "e", "m", "moment", "numerator", "M"}:
                values.append(integer(text))
            else:
                values.append(rational(text))
        key = tuple(values[:key_size])
        require(key not in keys, f"{path.name}: duplicate key {key}")
        keys.add(key)
        rows.append(tuple(values))
    if fields[0] == "id":
        require([row[0] for row in rows] == list(range(1, len(rows)+1)),
                f"{path.name}: missing, extra, or reordered indices")
    return rows


def moment(c, m):
    if c == 0:
        require(m == 1, "codimension-zero state")
        return 0
    left, degree, total = m, 0, 0
    while left:
        take = min(left, comb(c+degree-1, degree))
        total += degree*take
        left -= take
        degree += 1
    return total


RADII = {
    1: ["1"],
    2: ["1", "283/200", "1733/1000", "2", "2237/1000"],
    3: ["1", "63/50", "1443/1000", "397/250", "171/100",
        "909/500", "1913/1000", "2", "2081/1000"],
    4: ["1", "119/100", "1317/1000", "283/200", "187/125",
        "783/500", "1627/1000"],
}


def state_space():
    result = {}
    for d in range(1, 5):
        for e in range(d, 7):
            if d == 1 and e != 1:
                continue
            c = e-d
            # T(c,m) >= L(m-H), H=# monomials of degree < L.
            # Positive cap implies m < L H/(L-d/2); no trial cutoff.
            L = d//2+1
            H = comb(c+L-1, L-1) if c else 1
            limit = int(Q(L*H)/(L-Q(d, 2)))+1
            for m in ([1] if not c else range(c+1, limit+1)):
                mu = moment(c, m)
                cap = Q(d)-Q(2*mu, m)
                first = min(cap, Q(8-m, 2)) if d == 4 else cap
                if cap <= 0 or first <= 0:
                    continue
                require(m <= len(RADII[d]), "missing radius")
                radius = Q(RADII[d][m-1])
                require(radius**d >= m, "radius below the true radius")
                result[d, e, m] = (mu, cap, first, radius)
    require(len(result) == 28, "complete state census")
    return result


STATES = state_space()
r0, lam, theta, bs = Q(25, 27), Q(1731, 2000), Q(21, 20), Q(18, 5)
eta = Q(999, 1000)
cubic, quartic3, quartic4, quintic4 = (4, 6, 3), (3, 6, 4), (4, 6, 4), (4, 6, 5)
q3, q4, q5 = Q(1317, 1000), Q(283, 200), Q(187, 125)
qa, qe = Q(5000, 3997), Q(49267, 107919)


def children(state):
    d, e, m = state
    mu = STATES[state][0]
    for child, (_, cap, _, _) in STATES.items():
        if child[0] >= d or child[1] > e:
            continue
        if child[0] == d-1:
            cap = min(cap, Q(d*m-2*mu-child[2], m))
        if cap > 0:
            yield child, cap


@lru_cache(None)
def J(state, h, cap):
    require(cap >= 0, "negative scalar domain")
    radius = STATES[state][3]
    return max(radius-h, 0)*cap + max(
        [Q(0)] + [J(child, max(h, radius), min(cap, child_cap))
                   for child, child_cap in children(state)])


def price(state, order, bmax, quartic_cap=None):
    cap = min(bmax, STATES[state][2])
    if state == quartic3 and quartic_cap is not None:
        cap = min(cap, quartic_cap)
    return (7-bmax)*r0+bmax/order+J(state, 1/order, cap)


def check_states(rows):
    require({row[:3] for row in rows} == set(STATES), "states: complete index set")
    for row in rows:
        state = row[:3]
        mu, _, cap, radius = STATES[state]
        first = price(state, lam, bs, Q(1))
        # The cubic theta entry is an affine price, not a scalar tail.
        last = (3*r0+4/theta+(q3-1/theta)*Q(5, 2)
                if state == cubic else price(state, theta, Q(4)))
        expected = (mu, cap, radius, J(state, Q(50, 49), cap), first, last)
        require(row[3:] == expected, f"states: altered witness at {state}")
        if state not in {cubic, quartic4, quintic4}:
            require(first < Q(799, 100), f"lambda payment {state}")
        require(last < (7+Q(51, 100) if state == cubic else Q(799, 100)),
                f"theta payment {state}")
    return {"rows": 28, "checked_numeric_columns": 6,
            "theta_cubic_interpretation": "whole affine price",
            "scalar_domain": "all 28 positive-cap states and all dimension skips"}


def scalar_identities():
    values = {}
    values["eta_noncubic_max"] = max(
        J(v, 1/eta, x[2]) for v, x in STATES.items() if v != cubic)
    require(values["eta_noncubic_max"] <= Q(1209, 1000), "eta grouped bound")
    values["eta_grouped_price"] = 3*r0+4/eta+Q(1209, 1000)
    values["eta_post_error_margin"] = 8-values["eta_grouped_price"]-Q(1, 1000)
    require(values["eta_post_error_margin"] > Q(1, 400000), "cutting error margin")
    values["cubic_correction_max"] = max(
        J(w, q3, cap) for w, cap in children(cubic) if w != quartic3)
    require(values["cubic_correction_max"] <= Q(977, 2000), "cubic correction")
    values["cubic_quartic_low_correction"] = J(quartic3, q3, Q(1))
    values["cubic_affine_entry_slack"] = qe-Q(51, 100)-qa+q3
    require(values["cubic_affine_entry_slack"] > Q(1, 2000), "affine slack")
    values["quartic_intercept"] = J(quartic4, q4, STATES[quartic4][2])
    values["quintic_intercept"] = J(quintic4, q5, STATES[quintic4][2])
    values["direct_quartic_price"] = (7-bs)*r0+bs/lam+(qa-1/lam)*Q(3, 2)
    require(values["direct_quartic_price"] < 7+qe, "direct quartic payment")
    return values


def shift_factors(values):
    slopes = {3: q3, 4: q4, 5: q5}
    constants = {3: 7+Q(51, 100)-7*r0,
                 4: 8-Q(1, 500)-7*r0-values["quartic_intercept"],
                 5: 8-Q(1, 500)-7*r0-values["quintic_intercept"]}
    factors = {}
    for m, slope in slopes.items():
        a, c = slope-r0, constants[m]
        require(c > (4-Q(m, 2))*a > 0 and a*bs-c > 0,
                f"shift denominator and monotonicity for {m}")
        factors[m] = (bs-4+Q(m, 2))/(a*bs-c)
        values[f"shift_factor_{m}"] = factors[m]
    for m in (4, 5):
        for order in (lam, theta):
            require(factors[m]*(slopes[m]*order-1) >= factors[3]*(q3*order-1),
                    "affine shift domination at endpoints")
    return factors, constants


def check_transitions(rows, values, factors, constants):
    mesh = [lam]+[Q(j, 50) for j in (44, 45, 46, 47, 48, 49, 50)]+[theta]
    mass = 1-2*lam**5
    values["initial_mass"] = mass
    for row, lo, hi in zip(rows, mesh, mesh[1:]):
        shift = factors[3]*(q3*lo-1)
        loss = 2*(max(hi-shift, 0)**5-max(lo-shift, 0)**5)
        require(loss >= 0, "negative intrinsic interval")
        mass -= loss
        require(row[1:] == (lo, hi, shift, loss, mass),
                f"transitions: altered witness {row[0]}")
        require(mass > Q(1, 1000), "cumulative source reserve")
    values["terminal_mesh_mass"] = mass
    a2 = Q(3, 2)/(4*(q3-r0)-constants[3])
    z2 = a2*(q3*eta-1)
    require(lam < eta < theta and z2 < eta, "degree-two interval domain")
    values["degree_two_shift"] = z2
    values["degree_two_mass"] = 2-2*eta**5-2*((theta-z2)**5-(eta-z2)**5)
    require(values["degree_two_mass"] > Q(9, 1000), "degree-two source mass")
    values["degree_three_mass"] = 3-2*theta**5
    require(values["degree_three_mass"] > 0, "degree-three kernel")
    return {"rows": 8, "endpoint_count": 9,
            "minimum_partial_reserve": str(mass), "reserved_error": "1/1000"}


def check_rank_two(rows, powers, values, factors, constants):
    require({row[:3] for row in rows} == set(STATES), "rank-two price coverage")
    rho = Q(19, 20)  # Companion source; different from theta=21/20.
    for row in rows:
        state = row[:3]
        expected = (price(state, lam, Q(3), Q(1)), price(state, rho, Q(3)),
                    price(state, lam, bs, Q(1)), J(state, 1/eta, STATES[state][2]))
        require(row[3:] == expected, f"rank-two prices: altered witness {state}")
    values["rank_two_terminal_price"] = max(
        price(v, rho, Q(3)) for v in STATES if v != cubic)
    values["rank_two_initial_price"] = max(
        price(v, lam, Q(3), Q(1)) for v in STATES if v not in {cubic, quartic4})
    require(max(values["rank_two_terminal_price"], values["rank_two_initial_price"])
            < Q(7998, 1000), "rank-two scalar prices")
    k = 4*r0+3/lam
    values["rank_two_terminal_affine"] = 4*r0+3/rho+(q3-1/rho)*Q(7, 3)
    values["rank_two_small_deficit_affine"] = k+(q3-1/lam)*Q(21, 10)
    require(max(values["rank_two_terminal_affine"], values["rank_two_small_deficit_affine"])
            < 7+Q(51, 100), "rank-two affine prices")
    values["rank_two_distinct_factor_quartic"] = k+J(quartic4, 1/lam, Q(3, 2))
    require(values["rank_two_distinct_factor_quartic"] < Q(7998, 1000),
            "distinct-factor quartic payment")
    z4, z5 = factors[4]*(q4*lam-1), Q(4, 5)*factors[5]*(q5*lam-1)
    require(0 < z4 <= z5 < lam < eta, "full principal interval arguments")
    values["factorial_quartic_shift"], values["factorial_quintic_shift"] = z4, z5
    values["factorial_mass"] = values["initial_mass"]-2*((eta-z4)**5-(lam-z4)**5)
    require(values["factorial_mass"] > Q(1, 250), "factorial source mass")
    a3, a4 = q3-r0, q4-r0
    ratio = (q3*lam-1)/(3*a3-constants[3])
    values["rank_two_cubic_ratio"] = ratio
    require([row[0] for row in powers] == list(range(3, 10)), "power index coverage")
    for M, xi, shift in powers:
        KM = 2-Q(1, M)
        require(constants[3] > (4-KM)*a3 > 0, "cyclic coefficient monotonicity")
        require(xi == Q(M-1, 2*M-1)*ratio and shift == xi*(1+Q(1, M)),
                f"smooth-power altered witness M={M}")
        require(xi > Q(3, 8) and shift > Q(1, 2), "full character shifts")
    # Polynomial identity on the entire closed real interval [3,9].
    # -7M^2+80M-153=(9-M)(7M-17), denominators M(2M-1)>0.
    require((-153, 80, -7) == (-9*17, 9*7+17, -7),
            "smooth-power interval factorization")
    values["rank_two_cubic_uniform_shift"] = ratio*Q(80, 153)
    require(values["rank_two_cubic_uniform_shift"] > Q(1, 2), "uniform cubic shift")
    # For every integer M>=10, 2+1/M<=21/10; the affine price increases
    # with the deficit. This checks the cutoff, not an infinite sampling.
    require(q3-1/lam > 0 and 2+Q(1, 10) == Q(21, 10), "all-M cutoff")
    xi4 = Q(1, 2)*(q4*lam-1)/(3*a4-constants[4])
    y4 = (8-Q(1, 500)-k-values["quartic_intercept"])/(q4-1/lam)
    values["rank_two_quartic_coefficient"], values["rank_two_quartic_deficit"] = xi4, y4
    values["rank_two_quartic_shift"] = xi4*y4
    require(xi4 > Q(1, 3) and xi4*y4 > Q(1, 2), "quartic full shift")
    values["rank_two_mass"] = values["initial_mass"]-2*((rho-Q(1, 2))**5-(lam-Q(1, 2))**5)
    require(values["rank_two_mass"] > Q(1, 250), "rank-two companion source mass")
    return {"price_rows": 28, "price_columns": 4, "exceptional_power_rows": 7,
            "companion_order": "19/20", "companion_deficit_cap": "3",
            "all_power_reduction": "M>=10 paid by y<=2+1/M<=21/10; [3,9] factored inequality",
            "geometry_scope": "the article and README retain the analytic family hypotheses"}


zero = (Q(0),)*7


def poly(a, b, power=6, shift=0):
    return tuple(Q(0) if j < shift or j > shift+power else
                 Q(comb(power, j-shift))*a**(j-shift)*b**(power-j+shift)
                 for j in range(7))


def linear_combination(*terms):
    return tuple(sum(c*p[j] for c, p in terms) for j in range(7))


mon = poly(Q(1), Q(0))
r, v, T = Q(9, 10), Q(27, 25), Q(29, 25)
h, eps, H = Q(263, 2000), Q(1, 1000000), Q(5399983, 17000000)


def primitive_difference(p, lo, hi):
    return sum(p[j]*(hi**(j+1)-lo**(j+1))/Q(j+1) for j in range(7))


def candidates(which, line, t):
    a, b = line
    f = a*t+b
    if f <= 0:
        return [linear_combination((7, mon))]
    c = f/t
    if which == "old":
        if c >= 1:
            return [zero]
        ps = []
        if c <= Q(1, 7):
            ps.append(linear_combination((7, mon), (-7*Q(21, 4)**3, poly(a, b, 3, 3))))
        else:
            ps.extend([linear_combination((7, mon), (-Q(7, 24**3), poly(17+7*a, 7*b, 3, 3))),
                       linear_combination((7*Q(7, 6)**6, poly(1-a, -b)))])
        if c >= Q(1, 2):
            ps.append(linear_combination((7, poly(1-a, -b))))
        ps.append(linear_combination((7*(Q(256, 729)+3), mon), (-126, poly(a, b, 1, 5)))
                  if c <= Q(1, 6) else linear_combination((7*Q(256, 729), mon)))
        return ps
    es = [poly(1-j*a, -j*b) if t-j*f > 0 else zero for j in (1, 2, 3)]
    ps = [linear_combination((7, es[0]), (7, es[1]), (-7, es[2]))]
    if c >= Q(4, 21):
        ps.append(linear_combination((7, es[0])))
    return ps


@lru_cache(None)
def integrate(lo, hi, old, new):
    require(0 <= lo <= hi <= T, "integration domain")
    if lo == hi:
        return Q(0), 0
    knots = {lo, hi}
    for anchor in (r, Q(99, 100), Q(26, 25), v, T):
        if lo < anchor < hi:
            knots.add(anchor)
    for line in (old, new):
        if line is None:
            continue
        a, b = line
        for c in map(Q, ("0", "1/7", "1/6", "4/21", "1/3", "1/2", "1")):
            if a != c and lo < -b/(a-c) < hi:
                knots.add(-b/(a-c))
    for j in range(int(lo*100), int(hi*100)+1):
        if lo < Q(j, 100) < hi:
            knots.add(Q(j, 100))
    knots = sorted(knots)
    answer = Q(0)
    for left, right in zip(knots, knots[1:]):
        mid = (left+right)/2
        baseline = Q(1) if mid < r else Q(61, 100) if mid < Q(99, 100) else Q(11, 25) if mid < Q(26, 25) else Q(17, 50)
        ps = [linear_combination((7*baseline, mon))]
        if old is not None:
            ps += candidates("old", old, mid)
        if new is not None:
            ps += candidates("new", new, mid)
        value = min(primitive_difference(p, left, right) for p in ps)
        require(value >= 0, "negative candidate minimum")
        answer += value
    return answer, len(knots)-1


def numerator(value):
    return (10**9*value.numerator)//value.denominator+1


def check_ambient(prefix, boxes, values):
    ordered = sorted(row[1:] for row in prefix)
    require(ordered[0][0] == h/r and ordered[-1][1] == 3, "prefix outer endpoints")
    require(all(a[1] == b[0] for a, b in zip(ordered, ordered[1:])), "prefix coverage")
    receipts, cell_count = [], 0
    for index, lo, hi, expected in prefix:
        require(lo < hi, "empty or reversed prefix record")
        parts = [integrate(Q(0), r, (hi, h-hi*r), None),
                 integrate(r, v, (lo, h-lo*r), None)]
        value = sum(p[0] for p in parts)
        cell_count += sum(p[1] for p in parts)
        require(numerator(value) == expected and expected < 891000000,
                f"prefix {index}: altered bound or failed margin")
        receipts.append({"kind": "prefix", "id": index, "numerator": expected})
    area = Q(0)
    for position, (index, al, au, dl, du, expected) in enumerate(boxes):
        require(h/r <= al < au <= 3 and H/v <= dl < du <= 20, f"box {index}: domain")
        area += (au-al)*(du-dl)
        for _, bl, bu, el, eu, _ in boxes[:position]:
            require(min(au, bu) <= max(al, bl) or min(du, eu) <= max(dl, el),
                    f"box {index}: overlapping interiors")
        linked = max(H, h+(v-r)*al-eps)
        dd = max(dl, linked/v)
        require(dd <= du, f"box {index}: source table unexpectedly infeasible")
        oldpre, oldpost = (au, h-au*r), (al, h-al*r)
        newpre, newpost = (du, linked-du*v), (dd, linked-dd*v)
        parts = [integrate(Q(0), r, oldpre, newpre),
                 integrate(r, v, oldpost, newpre),
                 integrate(v, T, oldpost, newpost)]
        value = sum(p[0] for p in parts)
        cell_count += sum(p[1] for p in parts)
        require(numerator(value) == expected, f"box {index}: altered rational witness")
        require(expected <= 899892923, f"box {index}: failed uniform margin")
        receipts.append({"kind": "box", "id": index, "numerator": expected})
        if index % 50 == 0:
            print(f"Checked {index}/237 ambient boxes.", flush=True)
    require(area == (3-h/r)*(20-H/v), "closed rectangle coverage area")
    # Independent coverage by elementary vertical strips, including endpoint
    # fibers. Closed intervals cover all boundary segments and singleton ties.
    knots = sorted({row[j] for row in boxes for j in (1, 2)})
    probes = sorted(set(knots + [(a+b)/2 for a, b in zip(knots, knots[1:])]))
    for a in probes:
        intervals = sorted((dl, du) for _, al, au, dl, du, _ in boxes if al <= a <= au)
        covered = H/v
        for lo, hi in intervals:
            require(lo <= covered, "uncovered vertical interval")
            covered = max(covered, hi)
        require(covered == 20, "incomplete vertical fiber")
    values["ambient_area"] = area
    values["ambient_max_numerator"] = Q(max(row[-1] for row in boxes))
    require(values["ambient_max_numerator"] == 899892923, "ambient maximum identity")
    tail = r**7+integrate(r, T, (Q(3), h-3*r), None)[0]
    values["ambient_a_tail_numerator"] = Q(numerator(tail))
    require(v > H and 20 > 1 and H/v > Q(4, 21), "unbounded d-tail applicability")
    values["ambient_d_tail_numerator"] = Q(numerator(Q(891, 1000)+(v-H)**7/19))
    require(max(values["ambient_a_tail_numerator"], values["ambient_d_tail_numerator"])
            < 899892923, "unbounded total margin")
    h2 = Q(341, 1000)
    kappa, z0 = T-h2, 25*(T-h2)/21
    q0 = 1-2*h2/T
    z = Q(125)/(52+73*q0)
    require(0 < kappa < z0 < T < z < Q(125, 52) < Q(13, 5), "terminal interval endpoints")
    require(3*kappa-2*z0 > 0 and 2*kappa-T > 0 and q0 > 0, "terminal positive parts")
    values["ambient_incoming"] = 3*kappa**7-2*(2*kappa-z0)**7+(3*kappa-2*z0)**7-(2*kappa-T)**7
    values["ambient_outgoing"] = q0**6*(Q(125, 52)*z**6-T**7)
    require(values["ambient_incoming"] < Q(2, 3) and values["ambient_outgoing"] < Q(7, 50),
            "terminal small-slope margins")
    values["ambient_large_slope_loss"] = (T-2*h2)**7
    values["ambient_divisorial_bound"] = (1-h2)**7
    values["ambient_new_height"] = Q(29, 85)-Q(56, 27)*eps
    require(values["ambient_large_slope_loss"] < Q(1, 100) and
            values["ambient_divisorial_bound"] < Q(3, 50) and
            values["ambient_new_height"] > h2, "terminal threshold inequalities")
    require(numerator(Q(1, 2)) == 500000001, "strict rounding at an integral scaled value")
    require(integrate(r, r, (Q(0), Q(0)), None) == (0, 0), "zero integration interval")
    for c in map(Q, ("0", "1/7", "1/6", "4/21", "1/3", "1/2", "1")):
        integrate(Q(0), r, (c, Q(0)), (c, Q(0)))
    return {"prefix_intervals": 19, "closed_rectangles": 237,
            "integration_cells_reconstructed": cell_count, "coverage_fibers": len(probes),
            "maximum_numerator": int(values["ambient_max_numerator"]), "denominator": 10**9,
            "margin_below_nine_tenths": str(Q(9, 10)-values["ambient_max_numerator"]/10**9),
            "unbounded_ranges": ["a0>=3", "a0<=3,a1>=20"],
            "records": receipts}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data"))
    parser.add_argument("--report", type=Path, default=Path("report.json"))
    args = parser.parse_args(argv)
    report = {"schema": "intrinsic-double-check-report-v1", "success": False, "families": {}}
    try:
        require(args.data.is_dir(), "missing data directory")
        require({p.name for p in args.data.iterdir()} == set(SCHEMAS), "unexpected or missing data file")
        tables = {name: read_table(args.data/name, schema) for name, schema in SCHEMAS.items()}
        report["families"]["18-intrinsic-double/scalar"] = check_states(tables["states.tsv"])
        values = scalar_identities()
        factors, constants = shift_factors(values)
        report["families"]["18-intrinsic-double/transitions"] = check_transitions(
            tables["transitions.tsv"], values, factors, constants)
        report["families"]["18-intrinsic-double/rank-two"] = check_rank_two(
            tables["rank-two-prices.tsv"], tables["rank-two-powers.tsv"], values, factors, constants)
        report["families"]["18-intrinsic-double/ambient"] = check_ambient(
            tables["ambient-prefix.tsv"], tables["ambient-boxes.tsv"], values)
        supplied = dict(tables["identities.tsv"])
        require(set(supplied) == set(values), "identities: missing or extra named record")
        for name, value in values.items():
            require(supplied[name] == value, "identities: altered witness " + name)
        report["families"]["18-intrinsic-double/identities"] = {
            "count": len(values), "values": {k: str(v) for k, v in sorted(values.items())}}
        report["files"] = {name: hashlib.sha256((args.data/name).read_bytes()).hexdigest()
                           for name in SCHEMAS}
        report["success"] = True
        report["scope"] = "Exact arithmetic and complete declared data coverage; geometric implications require the article."
    except (Invalid, OSError, UnicodeError, ValueError, ZeroDivisionError, IndexError) as error:
        report["error"] = str(error)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps({"success": report["success"], "error": report.get("error"),
                      "family_ids": list(report["families"])}, sort_keys=True), flush=True)
    return 0 if report["success"] else 1


if __name__ == "__main__":
    sys.exit(main())
