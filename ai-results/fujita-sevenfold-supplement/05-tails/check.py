#!/usr/bin/env python3
"""Portable exact checker for 05-tails/ambient71 and double-fivefold28.

Python >= 3.9, standard library only. No source conversion, candidate search,
private imports, live graph, network, or floating-point arithmetic is used.
The affine test adapts volume_two_early_rigidity_independent_audit.py; the
support/flag evaluators adapt sub11_control1920_written_replay.py. See README.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
from math import comb
from pathlib import Path
import re
import sys


ZERO = Q(0)
LINE_COUNTS = (
    1, 1, 1, 1, 1, 1, 1, 3, 4, 1, 4, 3, 1, 1, 4, 2, 1, 1, 1, 1,
    3, 3, 4, 1, 1, 1, 1, 4, 4, 3, 3, 7, 4, 5, 5, 3, 3, 3, 5, 5,
    4, 5, 3, 3, 4, 4, 3, 3, 4, 4, 7, 6, 6, 4, 7, 6, 6, 7, 6, 5,
    6, 5, 3, 2, 1, 7, 7, 7, 6, 6, 4,
)
COMMON = {"degree_floor": "1", "root_grid": 1000, "curve": "smooth-only"}
PARAMETERS = {
    "ambient71": dict(COMMON, ambient_dimension=7, cap_rule="intrinsic"),
    "double-fivefold28": dict(COMMON, parent_state=[5, 6, 2], parent_moment=1,
                             support_slope="50/49", cap_rule="first-successor-only"),
}


class CertificateError(ValueError):
    pass


def require(ok, message):
    if not ok:
        raise CertificateError(message)


def keys(obj, expected, where):
    require(type(obj) is dict and set(obj) == set(expected), where + ": wrong fields")


def integer(value, where, low=0, high=10000):
    require(type(value) is int and low <= value <= high, where + ": invalid integer")
    return value


def rational(value, where):
    require(type(value) is str and len(value) <= 128 and
            re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?", value) is not None,
            where + ": malformed rational")
    result = Q(value)
    require(str(result) == value, where + ": rational must be reduced and canonical")
    return result


def no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def no_float(value):
    raise CertificateError("non-exact JSON number: " + value)


def read_data(path, family):
    obj = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=no_duplicates,
                     parse_float=no_float, parse_constant=no_float)
    keys(obj, ("schema_version", "family", "parameters", "rows"), family)
    require(type(obj["schema_version"]) is int and obj["schema_version"] == 1,
            family + ": unsupported schema")
    require(obj["family"] == "05-tails/" + family, family + ": wrong family")
    require(json.dumps(obj["parameters"], sort_keys=True) ==
            json.dumps(PARAMETERS[family], sort_keys=True), family + ": wrong parameters")
    require(type(obj["rows"]) is list, family + ": rows must be an array")
    result = []
    for rid, raw in enumerate(obj["rows"], 1):
        where = family + " row " + str(rid)
        extra = ("cap", "lines") if family == "ambient71" else (
            "ambient_id", "intrinsic_cap", "first_cap", "support")
        keys(raw, ("id", "state", "moment", "radius", *extra), where)
        require(integer(raw["id"], where + " id", 1) == rid, where + ": unordered/duplicate id")
        state = raw["state"]
        require(type(state) is list and len(state) == 3, where + ": malformed state")
        state = tuple(integer(x, where + " state", 1) for x in state)
        row = {"id": rid, "state": state, "moment": integer(raw["moment"], where + " moment"),
               "radius": rational(raw["radius"], where + " radius")}
        if family == "ambient71":
            row["cap"] = rational(raw["cap"], where + " cap")
            require(rid <= len(LINE_COUNTS) and type(raw["lines"]) is list and
                    len(raw["lines"]) == LINE_COUNTS[rid - 1], where + ": affine coverage")
            pairs = []
            for j, pair in enumerate(raw["lines"], 1):
                require(type(pair) is list and len(pair) == 2, where + ": malformed affine pair")
                a, p = (rational(x, where + " line " + str(j)) for x in pair)
                require(a >= 1 and p >= 0, where + ": affine coefficient out of range")
                pairs.append((a, p))
            require(len(set(pairs)) == len(pairs) and pairs == sorted(pairs),
                    where + ": duplicate or unordered affine pairs")
            row["lines"] = tuple(pairs)
        else:
            row["ambient_id"] = integer(raw["ambient_id"], where + " ambient_id", 1, 71)
            for field in ("intrinsic_cap", "first_cap", "support"):
                row[field] = rational(raw[field], where + " " + field)
        result.append(row)
    return result


def moment(c, m):
    """Fill each graded monomial capacity; independent of delivered moments."""
    if c == 0:
        require(m == 1, "zero codimension with nonsmooth multiplicity")
        return 0
    remaining, value, degree = m, 0, 0
    while remaining:
        capacity = 1 if degree == 0 else comb(c + degree - 1, degree)
        take = min(remaining, capacity)
        value += degree * take
        remaining -= take
        degree += 1
    return value


def enumerate_states(n, max_dimension):
    states, bounds = {}, []
    for d in range(1, max_dimension + 1):
        for e in range(d, n + 1):
            if d == 1 and e != 1:
                continue
            c = e - d
            ell = d // 2 + 1
            top = 1 if c == 0 else (2 * ell * comb(c + ell - 1, ell - 1) - 1) // (2 * ell - d)
            bounds.append({"d": d, "e": e, "inclusive_m_bound": top})
            for m in range(1 if c == 0 else c + 1, top + 1):
                t = moment(c, m)
                beta = Q(d) - Q(2 * t, m)
                if beta > 0:
                    states[(d, e, m)] = (t, beta)
    return states, bounds


def radius(d, m):
    """Smallest j/1000 with j**d >= m*1000**d; integer binary search."""
    low, high, target = 0, 1000 * m, m * 1000 ** d
    while high - low > 1:
        mid = (low + high) // 2
        if mid ** d >= target:
            high = mid
        else:
            low = mid
    return Q(high, 1000)


def edge(v, w, primitive):
    if w[0] >= v[0] or w[1] > v[1]:
        return ZERO
    cap = primitive[w][1]
    if w[0] == v[0] - 1:
        cap = min(cap, primitive[v][1] - Q(w[2], v[2]))
    return max(ZERO, cap)


def points(row, cap):
    require(cap >= 0, "negative affine comparison domain")
    result = {ZERO, cap}
    for (a, p), (b, q) in combinations(row["lines"], 2):
        if a != b:
            y = (q - p) / (a - b)
            if 0 < y < cap:
                result.add(y)
    return sorted(result)


def envelope(row, y):
    return min(a * y + p for a, p in row["lines"])


def ambient_check(rows, report):
    primitive, bounds = enumerate_states(7, 6)
    require([r["state"] for r in rows] == list(primitive) and len(rows) == 71,
            "ambient71: complete state coverage failed")
    for row in rows:
        v = row["state"]
        require((row["moment"], row["cap"]) == primitive[v], "ambient71: moment/cap mismatch " + str(v))
        require(row["radius"] == radius(v[0], v[2]), "ambient71: radius mismatch " + str(v))
    allowed = {row["id"]: [(nxt, min(row["cap"], edge(row["state"], nxt["state"], primitive)))
                           for nxt in rows if edge(row["state"], nxt["state"], primitive) > 0]
               for row in rows}
    tests, successor_tests, evaluations = [], 0, 0
    for row in rows:
        for j, (a, p) in enumerate(row["lines"], 1):
            extra, witness = ZERO, None
            for nxt, cap in allowed[row["id"]]:
                candidates = points(nxt, cap)
                value, y = max((envelope(nxt, y) - max(a, row["radius"]) * y, y)
                               for y in candidates)
                successor_tests += 1
                evaluations += len(candidates)
                if value > extra:
                    extra, witness = value, {"row": nxt["id"], "cap": str(cap), "y": str(y)}
            rhs = max(row["radius"] - a, ZERO) * row["cap"] + extra
            require(p >= rhs, "ambient71: failed affine inequality at row %s line %s: %s < %s" %
                    (row["id"], j, p, rhs))
            tests.append({"row": row["id"], "line": j, "slope": str(a), "intercept": str(p),
                          "rhs": str(rhs), "slack": str(p - rhs), "positive_successor_witness": witness})
    result = {"status": "PASS", "states": len(rows), "states_by_dimension": dict(Counter(r["state"][0] for r in rows)),
              "enumeration_bounds": bounds, "moments_checked": len(rows), "minimal_radii_checked": len(rows),
              "affine_lines": len(tests), "allowed_edges": sum(map(len, allowed.values())),
              "line_successor_tests": successor_tests, "endpoint_intersection_evaluations": evaluations,
              "minimum_affine_slack": str(min(Q(t["slack"]) for t in tests)), "line_checks": tests}
    report["families"]["05-tails/ambient71"] = result
    result["consequences"] = consequences(rows)
    return primitive


def consequences(rows):
    """Check all printed census/payment summaries, not just the line witnesses."""
    byid = {r["id"]: r for r in rows}
    records = []

    def failure_max(selected, u, truncate=None):
        return max(((7 - y) / u + envelope(r, y), r["id"])
                   for r in selected for y in points(r, r["cap"] if truncate is None else min(r["cap"], truncate)))

    def record(name, observed, expected, row=None):
        require(observed == Q(expected), name + ": printed value mismatch")
        require(observed < 8, name + ": missing strict payment margin")
        records.append({"test": name, "maximum": str(observed), "gap_to_8": str(8 - observed), "row": row})

    for dims, cap, value, rid in [((1, 2, 3), None, "20773/2600", 21),
                                  ((4,), Q(3, 2), "41553/5200", 39),
                                  ((5, 6), Q(7, 5), "38951/4875", 60)]:
        observed, winner = failure_max([r for r in rows if r["state"][0] in dims], Q(26, 25), cap)
        require(winner == rid, "order 26/25: maximizing row mismatch")
        record("order 26/25 dimensions " + str(dims), observed, value, rid)
    exceptional = (51, 55, 56, 57)
    for d, value, rid in [(1, "179/29", 1), (2, "95919/14500", 4), (3, "43109/5800", 21),
                          (4, "92763/11600", 39), (5, "115473/14500", 58)]:
        observed, winner = failure_max([r for r in rows if r["state"][0] == d and r["id"] not in exceptional], Q(29, 25))
        require(winner == rid, "order 29/25: maximizing row mismatch")
        record("order 29/25 dimension " + str(d), observed, value, rid)
    pairs = (("1149/1000", "247/250"), ("623/500", "39/40"),
             ("33/25", "77/100"), ("69/50", "673/1000"))
    for rid, pair, value in zip(exceptional, pairs, ("33107/4250", "135909/17000", "13041/1700", "127437/17000")):
        row = byid[rid]
        a, p = row["lines"][0]
        require((a, p) == tuple(map(Q, pair)) and a * Q(29, 25) > 1, "exceptional deficit pair/denominator mismatch")
        record("order 34/25 exceptional row " + str(rid), (7 - row["cap"]) / Q(34, 25) + a * row["cap"] + p, value, rid)
    divisor_records = []
    for rid, expected in zip(range(66, 72), ("73/125", "1267/1000", "761/500", "447/250", "1957/1000", "1501/375")):
        row = byid[rid]
        a, p = row["lines"][0]
        m, denominator = row["state"][2], 7 * a + p - 8
        require(denominator == Q(expected) and denominator > 0, "divisor denominator mismatch")
        value, slope = m * (a - 1) / denominator, m * a / denominator - Q(125, 73)
        require(value >= 0 and slope >= 0, "divisor fixed-mass inequality failed")
        divisor_records.append({"m": m, "D": str(denominator), "slack_at_u_1": str(value), "slope_slack": str(slope)})
    for d, expected in enumerate(("43/13", "51/13", "12879/2600", "76841/13000", "11234/1625", "12949/1625"), 1):
        observed, rid = failure_max([r for r in rows if r["state"][0] == d], Q(13, 5))
        record("order 13/5 dimension " + str(d), observed, expected, rid)
    final_margin = 8 - Q(1, 100) - Q(12949, 1625) - Q(1, 1000)
    require(final_margin > 0, "high-order point: cutting-error budget exhausted")
    return {"payment_maxima": records, "divisor_fixed_mass": divisor_records,
            "weighted_cutting_error_budget": "1/1000", "gap_after_error_to_8_minus_1_over_100": str(final_margin)}


def double_check(rows, ambient, report):
    primitive, bounds = enumerate_states(6, 4)
    first_cap = {v: min(beta, Q(8 - v[2], 2)) if v[0] == 4 else beta
                 for v, (_, beta) in primitive.items()}
    expected = [v for v in primitive if first_cap[v] > 0]
    require([r["state"] for r in rows] == expected and len(rows) == 28,
            "double-fivefold28: complete first-successor coverage failed")
    bystate = {r["state"]: r for r in rows}
    ambient_by_state = {r["state"]: r for r in ambient}
    for row in rows:
        v = row["state"]
        require((row["moment"], row["intrinsic_cap"]) == primitive[v], "double-fivefold28: moment/intrinsic cap mismatch")
        require(row["first_cap"] == first_cap[v], "double-fivefold28: inherited first cap mismatch " + str(v))
        require(row["radius"] == radius(v[0], v[2]), "double-fivefold28: radius mismatch")
        require(row["ambient_id"] == ambient_by_state[v]["id"], "double-fivefold28: ambient mapping mismatch")
    children = {v: [(w, edge(v, w, primitive)) for w in primitive if edge(v, w, primitive) > 0] for v in expected}
    require(all(w in bystate for choices in children.values() for w, cap in choices), "restricted domain is not closed under descendants")

    @lru_cache(None)
    def support(v, a, h):
        require(0 <= h <= primitive[v][1], "support argument outside intrinsic domain")
        r = bystate[v]["radius"]
        return max(r - a, ZERO) * h + max([ZERO] + [support(w, max(a, r), min(h, cap)) for w, cap in children[v]])

    @lru_cache(None)
    def flags(v):
        return ((v,),) + tuple((v,) + rest for w, cap in children[v] for rest in flags(w))

    def flag_value(flag, a, h):
        caps = [h]
        for v, w in zip(flag, flag[1:]):
            caps.append(min(caps[-1], edge(v, w, primitive)))
        grid = sorted({ZERO, *caps})
        coeffs = [bystate[flag[0]]["radius"] - a]
        coeffs += [bystate[w]["radius"] - bystate[v]["radius"] for v, w in zip(flag, flag[1:])]
        return max(sum((q * y for q, y in zip(coeffs, ys)), ZERO)
                   for ys in product(grid, repeat=len(flag))
                   if all(0 <= y <= cap for y, cap in zip(ys, caps)) and
                   all(y >= z for y, z in zip(ys, ys[1:])))

    checks = []
    for row in rows:
        v, h, a = row["state"], row["first_cap"], Q(50, 49)
        observed = support(v, a, h)
        require(observed == row["support"], "double-fivefold28: support witness mismatch at " + str(v))
        expanded = max(flag_value(flag, a, h) for flag in flags(v))
        require(expanded == observed, "double-fivefold28: full flag evaluation disagrees at " + str(v))
        require(support(v, a, ZERO) == 0, "zero-interval support is not zero")
        checks.append({"state": list(v), "first_cap": str(h), "support": str(observed),
                       "expanded_flag_support": str(expanded), "flags": len(flags(v))})
    maximum = max(r["support"] for r in rows)
    require(maximum == Q(384791, 294000) < Q(1309, 1000) < Q(131, 100), "double-fivefold28: maximum/margins failed")
    require([r["state"] for r in rows if r["support"] == maximum] == [(4, 6, 3)], "double-fivefold28: maximizing state mismatch")
    curve = support((1, 1, 1), Q(20, 21), Q(1))
    require(curve == Q(1, 21), "curve correction below slope one omitted")
    repairs = []
    for m, a, cap in ((8, Q(2), Q(1, 4)), (9, Q(9, 4), Q(1, 9))):
        v = (3, 6, m)
        require(primitive[v][1] == cap and support(v, a, cap) == 0, "repaired threefold tail failed")
        surfaces = [w[2] for w, h in children[v] if w[0] == 2]
        require(surfaces == ([1] if m == 8 else []), "repaired threefold surface coverage failed")
        repairs.append({"state": list(v), "moment": primitive[v][0], "cap": str(cap),
                        "cost_slope": str(a), "support": "0", "surface_multiplicities": surfaces})
    report["families"]["05-tails/double-fivefold28"] = {
        "status": "PASS", "states": len(rows), "states_by_dimension": dict(Counter(v[0] for v in expected)),
        "enumeration_bounds": bounds, "intrinsic_states_before_parent_cap": len(primitive),
        "excluded_by_first_cap": [list(v) for v in primitive if first_cap[v] <= 0],
        "allowed_descendant_edges": sum(map(len, children.values())), "moments_checked": len(rows),
        "minimal_radii_checked": len(rows), "support_values_checked": len(checks),
        "flags_checked": sum(len(flags(v)) for v in expected), "zero_intervals_checked": len(rows),
        "support_checks": checks, "maximum": str(maximum), "gap_to_1309_over_1000": str(Q(1309, 1000) - maximum),
        "gap_to_131_over_100": str(Q(131, 100) - maximum), "subunit_curve_support_at_h_1": str(curve),
        "embedding_six_repairs": repairs,
        "embedding_seven_rows": [{"state": r["state"], "moment": r["moment"], "cap": str(r["cap"])}
                                  for r in ambient if r["state"][:2] == (3, 7) and r["state"][2] >= 8]}


def validate(data):
    require(sys.version_info >= (3, 9), "Python 3.9 or later is required")
    expected_files = {name + ".json" for name in PARAMETERS}
    require(data.is_dir() and {p.name for p in data.iterdir()} == expected_files,
            "data directory must contain exactly ambient71.json and double-fivefold28.json")
    require(all((data / name).is_file() and not (data / name).is_symlink() for name in expected_files),
            "data payloads must be ordinary files")
    report = {"status": "PASS", "checker_schema": 1, "python_version": sys.version.split()[0],
              "arithmetic": "exact integer and fractions.Fraction", "families": {},
              "data_sha256": {name: hashlib.sha256((data / name).read_bytes()).hexdigest() for name in sorted(expected_files)}}
    rows71 = read_data(data / "ambient71.json", "ambient71")
    rows28 = read_data(data / "double-fivefold28.json", "double-fivefold28")
    ambient_check(rows71, report)
    double_check(rows28, rows71, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    try:
        require(args.report.resolve().parent != args.data.resolve(), "report must be outside data directory")
        result = validate(args.data)
    except (ValueError, OSError, TypeError, KeyError, OverflowError) as exc:
        result = {"status": "FAIL", "error": str(exc), "families": {}}
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "families": list(result["families"]),
                      **({"error": result["error"]} if "error" in result else {})}))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
