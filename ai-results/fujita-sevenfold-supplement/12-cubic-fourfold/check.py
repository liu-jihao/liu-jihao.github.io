#!/usr/bin/env python3
"""Portable exact checker; Python >=3.10, standard library only.

Adapted from the independent sub9 state/path audit and the cubic-affine
part of the corrected-3672 audit. No producer, archived proof, live store,
optimizer, network, or host path is used. The delivered data are inputs.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
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


def keys(obj, expected, where):
    require(type(obj) is dict and set(obj) == set(expected),
            f"{where}: missing or extra fields")


def rational(value):
    require(type(value) is str and re.fullmatch(r"-?(?:0|[1-9]\d*)(?:/[1-9]\d*)?", value),
            f"malformed rational: {value!r}")
    result = Q(value)
    require(str(result) == value, f"noncanonical rational: {value!r}")
    return result


def integer(value):
    require(type(value) is int, f"noninteger: {value!r}")
    return value


def state_id(value):
    require(type(value) is list and len(value) == 3, "invalid state tuple")
    result = tuple(integer(x) for x in value)
    require(all(x > 0 for x in result), "nonpositive state coordinate")
    return result


def unique_pairs(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON field: {key}")
        result[key] = value
    return result


def read_json(path):
    def invalid_number(value):
        raise Invalid(f"floating or nonfinite JSON number: {value}")
    return json.loads(path.read_text(encoding="utf-8"),
                      object_pairs_hook=unique_pairs,
                      parse_float=invalid_number, parse_constant=invalid_number)


def records(rows, expected_fields):
    require(type(rows) is list, "rows must be a list")
    result = {}
    for row in rows:
        keys(row, expected_fields, "row")
        v = state_id(row["state"])
        require(v not in result, f"duplicate state: {v}")
        result[v] = row
    return result


def moment(c, m):
    if c == 0:
        require(m == 1, "regular state has multiplicity one")
        return 0
    remaining, degree, answer = m, 0, 0
    while remaining:
        capacity = comb(c + degree - 1, degree)
        take = min(remaining, capacity)
        answer += degree * take
        remaining -= take
        degree += 1
    return answer


def intrinsic(v):
    d, e, m = v
    return Q(d) - Q(2 * moment(e - d, m), m)


def first_cap(v):
    return min(intrinsic(v), Q(8 - v[2], 3)) if v[0] == 3 else intrinsic(v)


def census():
    positive, excluded = [(1, 1, 1)], []
    for d in (2, 3):
        for e in range(d, 7):
            if d == e:
                positive.append((d, e, 1))
                continue
            m = e - d + 1
            while intrinsic((d, e, m)) > 0:
                positive.append((d, e, m))
                m += 1
            excluded.append((d, e, m))
    return positive, excluded


def transition(v, w):
    require(w[0] < v[0] and w[1] <= v[1], "incompatible transition")
    cap = intrinsic(w)
    if w[0] == v[0] - 1:
        cap = min(cap, Q(v[0] * v[2] - 2 * moment(v[1] - v[0], v[2]) - w[2], v[2]))
    return cap


def solve(matrix, right):
    """Unique solution of a square rational system, or None if singular."""
    n = len(right)
    rows = [list(a) + [b] for a, b in zip(matrix, right)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j]), None)
        if pivot is None:
            return None
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = rows[j][j]
        rows[j] = [x / scale for x in rows[j]]
        for i in range(n):
            if i != j:
                scale = rows[i][j]
                rows[i] = [x - scale * y for x, y in zip(rows[i], rows[j])]
    return tuple(row[-1] for row in rows)


def inequalities(caps):
    # y_i >= 0; y_i <= caps[i]; and y_i <= y_(i-1).
    n = len(caps)
    result = []
    for i, cap in enumerate(caps):
        unit = tuple(Q(int(j == i)) for j in range(n))
        result.extend([(tuple(-x for x in unit), Q(0)), (unit, cap)])
        if i:
            row = tuple(Q(int(j == i) - int(j == i - 1)) for j in range(n))
            result.append((row, Q(0)))
    return result


@lru_cache(None)
def vertices(caps):
    n = len(caps)
    constraints = inequalities(caps)
    result = set()
    for selected in combinations(constraints, n):
        candidate = solve([a for a, _ in selected], [b for _, b in selected])
        if candidate is not None and all(
                sum(a * x for a, x in zip(row, candidate)) <= b
                for row, b in constraints):
            result.add(candidate)
    return tuple(sorted(result))


class Engine:
    def __init__(self, states, radii):
        self.states = tuple(states)
        self.radii = radii
        self.domains = {}
        self.requests = {}

    @lru_cache(None)
    def paths(self, v):
        # Enumerate every dimension/embedding-compatible path, including
        # empty or zero-cap paths. Every prefix is a direct-point ending.
        result = [(v,)]
        for w in self.states:
            if w[0] < v[0] and w[1] <= v[1]:
                result.extend((v,) + path for path in self.paths(w))
        return tuple(result)

    def surplus(self, v, kappa, a):
        require(kappa >= 1 and a <= intrinsic(v), "surplus outside declared domain")
        candidates = []
        for path in self.paths(v):
            caps = (a,) + tuple(transition(x, y) for x, y in zip(path, path[1:]))
            vv = vertices(caps)
            domain_key = (path, caps)
            self.domains[domain_key] = vv
            slopes = (self.radii[path[0]] - kappa,) + tuple(
                self.radii[y] - self.radii[x] for x, y in zip(path, path[1:]))
            candidates.extend((sum(c * z for c, z in zip(slopes, point)), path, point)
                              for point in vv)
        answer = max((row[0] for row in candidates), default=None)
        request_key = (v, kappa, a)
        self.requests[request_key] = {
            "state": list(v), "kappa": str(kappa), "cap": str(a),
            "maximum": None if answer is None else str(answer),
            "maximizers": [{"path": [list(w) for w in path],
                            "point": list(map(str, point))}
                           for val, path, point in candidates if val == answer]}
        return answer

    def price(self, v, gamma, b, sigma=Q(0)):
        value = self.surplus(v, 1 / gamma, min(b, first_cap(v)))
        require(value is not None, "price requested for empty first-root domain")
        return Q(25, 27) * (7 - b - sigma) + b / gamma + value


def table_view(rows, kind):
    header = (r" v&(3/4,2)&(3/4,8/3)&(41/50,8/3)\\ \hline"
              if kind == "integral-prices" else r" v&\lambda&\theta_0&\theta_1\\ \hline")
    lines = [r"\begin{array}{c|r|r|r}", header]
    for row in rows:
        lines.append(" (" + ",".join(map(str, row["state"])) + ")&" +
                     "&".join(map(str, row["ceilings"])) + r"\\")
    return "\n".join(lines + [r"\end{array}", ""])


def affine_view(rows):
    lines = [r"\begin{array}{c|c@{\qquad}c|c}",
             r"v&J_v(q,C_v)&v&J_v(q,C_v)\\ \hline"]
    for j in range(0, len(rows), 2):
        terms = []
        for row in rows[j:j+2]:
            terms.extend(["(" + ",".join(map(str, row["state"])) + ")", row["surplus"]])
        lines.append("&".join(terms) + r"\\")
    return "\n".join(lines + [r"\end{array}", ""])


def check(data):
    required_files = {"states.json", "integral-prices.json", "reduced-prices.json",
                      "affine.json", "claims.json", "integral-prices.tex",
                      "reduced-prices.tex", "affine.tex"}
    require(data.is_dir(), "missing data directory")
    require({p.name for p in data.iterdir()} == required_files,
            "missing or unexpected data file")
    require(all((data / name).is_file() for name in required_files), "nonfile in data")
    positive, excluded = census()
    state_data = read_json(data / "states.json")
    keys(state_data, ["schema", "family", "rows", "first_excluded", "embedding_seven"], "states")
    require(state_data["schema"] == 1 and type(state_data["schema"]) is int and
            state_data["family"] == "12-cubic-fourfold/states", "wrong states schema")
    rows = records(state_data["rows"], ["state", "moment", "intrinsic_cap", "root_cap", "radius"])
    require(list(rows) == positive and len(rows) == 19, "incomplete or reordered state census")
    radii = {}
    for v, row in rows.items():
        require(integer(row["moment"]) == moment(v[1] - v[0], v[2]), f"wrong moment: {v}")
        require(rational(row["intrinsic_cap"]) == intrinsic(v), f"wrong intrinsic cap: {v}")
        require(rational(row["root_cap"]) == first_cap(v), f"wrong parent cap: {v}")
        radii[v] = rational(row["radius"])
        require(radii[v] >= 1 and radii[v] ** v[0] >= v[2], f"invalid radius: {v}")
        require((1000*radii[v]).denominator == 1 and
                (radii[v]-Q(1, 1000))**v[0] < v[2], f"radius differs from the stated upward grid: {v}")
    sentinels = records(state_data["first_excluded"], ["state", "moment", "cap"])
    require(list(sentinels) == excluded, "incomplete first-excluded census")
    for v, row in sentinels.items():
        require(integer(row["moment"]) == moment(v[1] - v[0], v[2]) and
                rational(row["cap"]) == intrinsic(v) <= 0, "bad exclusion witness")
    e7 = records(state_data["embedding_seven"], ["state", "moment", "cap"])
    require(list(e7) == [(3, 7, 8), (3, 7, 9)], "wrong embedding-seven comparison")
    for v, row in e7.items():
        require(integer(row["moment"]) == moment(4, v[2]) and
                rational(row["cap"]) == intrinsic(v), "bad embedding-seven cap")
    engine = Engine(positive, radii)
    admissible = [v for v in positive if first_cap(v) > 0]
    require(len(admissible) == 17, "wrong positive-root count")
    params = {
        "integral-prices": [("3/4", "2"), ("3/4", "8/3"), ("41/50", "8/3")],
        "reduced-prices": [("3/4", "5/2"), ("381/500", "5/2"), ("39/50", "5/2")]}
    report = {"status": "pass", "scope": "Exact finite arithmetic, conditional on the article's geometric reductions.",
              "families": {}, "arithmetic": "fractions.Fraction; no floating point"}
    computed_tables = {}
    for family, settings in params.items():
        payload = read_json(data / (family + ".json"))
        keys(payload, ["schema", "family", "columns", "denominator", "rows", "maxima"], family)
        require(type(payload["schema"]) is int and payload["schema"] == 1 and
                payload["family"] == "12-cubic-fourfold/" + family, "wrong table identity")
        expected_columns = [{"gamma": g, "b": b, "sigma": "0", "n": 7} for g, b in settings]
        require(payload["columns"] == expected_columns, "wrong parameter mapping")
        require(integer(payload["denominator"]) == 10000, "wrong rounding denominator")
        table = records(payload["rows"], ["state", "costs", "ceilings"])
        require(list(table) == admissible, f"{family}: missing/extra/reordered state")
        values = []
        for v, row in table.items():
            require(type(row["costs"]) is list and len(row["costs"]) == 3 and
                    type(row["ceilings"]) is list and len(row["ceilings"]) == 3, "wrong column count")
            exact = [engine.price(v, Q(g), Q(b)) for g, b in settings]
            require([rational(x) for x in row["costs"]] == exact, f"{family}: wrong costs at {v}")
            for x, number in zip(exact, row["ceilings"]):
                number = integer(number)
                require(Q(number - 1) < 10000 * x <= number, f"{family}: invalid ceiling at {v}")
            values.append(exact)
        maxima = [max(row[j] for row in values) for j in range(3)]
        require(type(payload["maxima"]) is list and
                [rational(x) for x in payload["maxima"]] == maxima, "wrong table maxima")
        require((data / (family + ".tex")).read_text() == table_view(payload["rows"], family),
                f"{family}: rendered TeX differs from canonical data")
        unpaid = [[list(v) for v, row in zip(admissible, values) if row[j] >= Q(7999, 1000)]
                  for j in range(3)]
        expected_unpaid = ([[], [[3, 6, 4], [3, 6, 5]], []] if family == "integral-prices"
                           else [[[3, 6, 4]], [[3, 6, 4]], []])
        require(unpaid == expected_unpaid, "wrong complete changed-center classification")
        computed_tables[family] = maxima
        report["families"]["12-cubic-fourfold/" + family] = {
            "rows": len(table), "entries": 3 * len(table), "maxima": list(map(str, maxima)),
            "unpaid_states_at_7999_over_1000": unpaid}

    affine = read_json(data / "affine.json")
    keys(affine, ["schema", "family", "q", "intercept", "quartic_slope",
                  "quartic_intercept", "rows", "quartic_low"], "affine")
    require(type(affine["schema"]) is int and affine["schema"] == 1 and
            affine["family"] == "12-cubic-fourfold/affine", "wrong affine identity")
    q, intercept, qa, qe = map(rational, [affine["q"], affine["intercept"],
                                         affine["quartic_slope"], affine["quartic_intercept"]])
    require((q, intercept, qa, qe) == (Q(1317, 1000), Q(51, 100), Q(5000, 3997), Q(49267, 107919)),
            "wrong shared affine interface")
    require(q**4 >= 3 and qa < q, "bad affine slope comparison")
    quartic = (3, 6, 4)
    correction = records(affine["rows"], ["state", "surplus"])
    require(list(correction) == [v for v in admissible if v != quartic], "affine coverage")
    for v, row in correction.items():
        require(rational(row["surplus"]) == engine.surplus(v, q, first_cap(v)),
                f"wrong affine correction: {v}")
    require(rational(affine["quartic_low"]) == engine.surplus(quartic, q, Q(1)),
            "wrong quartic split")
    require((data / "affine.tex").read_text() == affine_view(affine["rows"]), "affine rendered mismatch")

    claims = read_json(data / "claims.json")
    keys(claims, ["schema", "family", "values"], "claims")
    require(type(claims["schema"]) is int and claims["schema"] == 1 and
            claims["family"] == "12-cubic-fourfold/claims", "wrong claims identity")
    result = {}
    result["affine_nonquartic_max"] = max(rational(r["surplus"]) for r in correction.values())
    result["affine_entry_margin"] = qe - intercept - qa + q
    result["affine_nonquartic_post_error_margin"] = 1 - intercept - Q(977, 2000) - Q(1, 2000)
    result["affine_low_quartic_post_error_margin"] = 1 - intercept - rational(affine["quartic_low"]) - Q(1, 2000)
    require(result["affine_nonquartic_max"] < Q(977, 2000), "affine correction reserve")
    require(result["affine_entry_margin"] > Q(1, 2000), "affine entry does not absorb errors")
    require(min(result[k] for k in ("affine_nonquartic_post_error_margin",
                                   "affine_low_quartic_post_error_margin")) > Q(1, 400000),
            "affine endpoint margin")
    require(Q(1, 200000) > Q(1, 400000), "imported quartic endpoint margin")
    for k, slope, cap, tail, target_order in [
            (4, Q(397, 250), Q(4, 3), Q(103, 500), Q(15471, 36775)),
            (5, Q(171, 100), Q(1), Q(29, 500), Q(45765, 102278))]:
        actual = engine.surplus((3, 6, k), slope, cap)
        require(actual == tail, "incomplete upper tail")
        A = slope - Q(25, 27)
        D = 8 - Q(1, 1000) - Q(175, 27) - tail
        require(D > Q(8-k, 3) * A > 0 and Q(8, 3)*A-D > 0, "forced-order denominator")
        eta = (slope * Q(3, 4) - 1)/(Q(8, 3)*A-D)
        require(eta == target_order and eta > Q(21, 50), "forced prime order")
        result[f"tail_{k}"] = actual
        result[f"A_{k}"], result[f"D_{k}"], result[f"eta_{k}"] = A, D, eta
    for k, slope in [(8, Q(2)), (9, Q(9, 4))]:
        require(engine.surplus((3, 6, k), slope, intrinsic((3, 6, k))) == 0,
                "multiplicity-eight/nine full tail")
    for v in positive:
        require(engine.surplus(v, Q(1), Q(0)) == 0, "zero-deficit endpoint")
    require(engine.surplus((3, 6, 8), Q(2), first_cap((3, 6, 8))) == 0,
            "multiplicity-eight point domain")
    require(engine.surplus((3, 6, 9), Q(9, 4), first_cap((3, 6, 9))) is None,
            "multiplicity-nine empty first-root domain")
    require(vertices((Q(-1, 3),)) == (), "empty interval handling")
    require(vertices((Q(0), Q(1))) == ((Q(0), Q(0)),), "point handling")
    require(vertices((Q(1), Q(0))) == ((Q(0), Q(0)), (Q(1), Q(0))), "segment handling")
    # All b>=2 surpluses used for the slack reduction have stabilized.
    for v in admissible:
        for gamma in (Q(3, 4), Q(381, 500), Q(39, 50)):
            require(engine.surplus(v, 1/gamma, min(Q(2), first_cap(v))) ==
                    engine.surplus(v, 1/gamma, min(Q(8, 3), first_cap(v))),
                    "surplus stabilization")
            require(1/gamma > Q(25, 27) and
                    (1/gamma-Q(25, 27))/2-Q(25, 27) < 0, "parameter monotonicity")
    result["low_b_price"] = max(engine.price(v, Q(3, 4), Q(19, 8)) for v in admissible)
    result["large_sigma_price"] = max(engine.price(v, Q(3, 4), Q(8, 3), Q(1, 8)) for v in admissible)
    result["later_low_deficit_price"] = Q(25, 6) + Q(5, 2)/Q(381, 500) + engine.surplus(
        quartic, Q(500, 381), Q(5, 4))
    paid = [result[k] for k in ("low_b_price", "large_sigma_price", "later_low_deficit_price")]
    paid += [computed_tables["reduced-prices"][-1], computed_tables["integral-prices"][0],
             computed_tables["integral-prices"][-1]]
    require(all(x < Q(7999, 1000) and x + Q(1, 10000) < 8-Q(1, 2000) for x in paid),
            "complete paid endpoint fails after cutting errors")
    result["dimension_margin_slope"] = 1-Q(25, 27)
    lam, theta, theta0, theta1 = Q(3, 4), Q(41, 50), Q(381, 500), Q(39, 50)
    result["initial_mass"] = 1-3*lam**4
    result["integral_loss"] = 3*((theta-Q(14, 25))**4-(lam-Q(14, 25))**4)
    result["integral_mass"] = result["initial_mass"]-result["integral_loss"]
    result["first_reduced_mass"] = 1-lam**4-2*theta0**4
    result["final_loss"] = 3*((theta1-Q(1, 2))**4-(theta0-Q(1, 2))**4)
    result["final_mass"] = result["first_reduced_mass"]-result["final_loss"]
    result["first_keta"] = Q(9, 4)-3/Q(19, 8)
    result["first_eta"] = result["first_keta"]/4
    result["second_eta"] = theta0*Q(25, 36)
    require(result["initial_mass"] > 0 and result["integral_mass"] > Q(1, 25) and
            result["first_reduced_mass"] > Q(9, 1000) and result["final_mass"] > Q(1, 250),
            "source mass bound")
    require(result["first_keta"] > theta0 and result["first_eta"] > Q(1, 5) and
            result["second_eta"] > Q(1, 2), "component-shift bound")
    result["final_mass_after_intrinsic_error"] = result["final_mass"]-Q(1, 500)
    require(result["final_mass_after_intrinsic_error"] > 0, "intrinsic error consumes source")
    keys(claims["values"], result, "claim values")
    require({k: rational(v) for k, v in claims["values"].items()} == result,
            "delivered claim differs from recomputed rational identity")
    report["families"]["12-cubic-fourfold/states"] = {
        "positive_states": len(positive), "first_excluded": len(excluded), "positive_roots": len(admissible),
        "embedding_seven_comparisons": len(e7), "multiplicities_eight_nine_full_tails_checked": True}
    report["families"]["12-cubic-fourfold/affine"] = {
        "nonquartic_rows": len(correction), "quartic_low_deficit_cases": 1,
        "high_deficit_domain": "1<z<=4/3", "entry_margin": str(result["affine_entry_margin"])}
    report["families"]["12-cubic-fourfold/claims"] = {"identities": len(result),
        "values": {k: str(v) for k, v in result.items()},
        "paid_prices": list(map(str, paid)),
        "post_error_gaps_to_endpoint": [str(8-Q(1, 2000)-x-Q(1, 10000)) for x in paid]}
    domains = [{"path": [list(v) for v in path], "caps": list(map(str, caps)),
                "vertices": [list(map(str, v)) for v in vv]}
               for (path, caps), vv in sorted(engine.domains.items())]
    report["coverage"] = {"distinct_paths": len({path for path, _ in engine.domains}),
        "path_domains": len(domains), "vertex_records": sum(len(d["vertices"]) for d in domains),
        "empty_domains": sum(not d["vertices"] for d in domains),
        "point_domains": sum(len(d["vertices"]) == 1 for d in domains),
        "segment_domains": sum(len(d["vertices"]) == 2 for d in domains),
        "surplus_requests": len(engine.requests),
        "all_paths_and_vertices": domains,
        "all_maximizers": [engine.requests[k] for k in sorted(engine.requests)]}
    require(report["coverage"]["empty_domains"] > 0 and report["coverage"]["segment_domains"] > 0 and
            report["coverage"]["point_domains"] > 0,
            "degenerate coverage absent")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    try:
        report, code = check(args.data), 0
    except (Invalid, OSError, ValueError, TypeError, KeyError, IndexError, ZeroDivisionError) as error:
        report, code = {"status": "fail", "error": str(error)}, 1
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    brief = {k: v for k, v in report.items() if k != "coverage"}
    if "coverage" in report:
        brief["coverage"] = {k: v for k, v in report["coverage"].items()
                             if not k.startswith("all_")}
    print(json.dumps(brief, sort_keys=True))
    return code


if __name__ == "__main__":
    sys.exit(main())
