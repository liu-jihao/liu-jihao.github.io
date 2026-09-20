#!/usr/bin/env python3
"""Independent exact validator for the delivered Section 16 certificate.

Adapted from the full-flag evaluator in sub11_control1920_written_replay.py
and the independent state/Bellman check in sub7_factorial_double5_cold_audit.py.
No producer, historical receipt, source prose, live graph, or network is read.
All failures use explicit exceptions and remain active under python -O.
"""
import sys
sys.dont_write_bytecode = True
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement
import json
from math import comb
from pathlib import Path
import re
from contract import SCHEMA, PARAMETERS, GROUPS, SQUARE_GROUPS, group_for, state_id, query_domain
from render import views


class Invalid(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise Invalid(message)


def keys(obj, expected, name):
    require(type(obj) is dict and set(obj) == set(expected), name + ": wrong keys")


def integer(value, name):
    require(type(value) is int, name + ": expected an integer")
    return value


def rational(value, name):
    require(type(value) is str and re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?", value),
            name + ": malformed rational")
    return F(value)


def no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def read_json(path):
    return json.loads(path.read_text(), object_pairs_hook=no_duplicates,
                      parse_constant=lambda s: (_ for _ in ()).throw(Invalid("nonfinite JSON: " + s)))


def indexed(rows, field, expected, name):
    require(type(rows) is list, name + ": expected a list")
    result = {}
    for row in rows:
        require(type(row) is dict and field in row, name + ": malformed record")
        ident = row[field]
        require(type(ident) in (str, int) and ident not in result, name + ": duplicate/invalid record")
        result[ident] = row
    require(set(result) == set(expected), name + ": missing or extra coverage records")
    return result


def moment(c, m):
    """Closed cumulative-binomial formula, independent of producer packing."""
    if c == 0:
        require(m == 1, "embedding equals dimension requires multiplicity one")
        return 0
    q = 0
    while comb(c+q, q) < m:
        q += 1
    return q*m - (comb(c+q, q-1) if q else 0)


def upward_root(m, d):
    low, high = 0, 1000*m
    while high-low > 1:
        middle = (high+low)//2
        if middle**d >= m*1000**d:
            high = middle
        else:
            low = middle
    return F(high, 1000)


def enumerate_states():
    states, exclusions = {}, []
    for d in range(1, 5):
        for e in range(d, 7):
            if d == 1 and e != 1:
                continue
            m = e-d+1
            while True:
                t = moment(e-d, m)
                beta = F(d)-F(2*t, m)
                cap = min(beta, F(8-m, 2)) if d == 4 else beta
                if cap <= 0:
                    exclusions.append({"state": state_id((d, e, m)), "beta": str(beta), "cap": str(cap)})
                    break
                states[d, e, m] = {"moment": t, "beta": beta, "cap": cap, "q": upward_root(m, d)}
                if e == d:
                    break
                m += 1
    return states, exclusions


def expected_label(name, square=False):
    if name == "d<=2":
        return r"\(d\le2\)"
    if name == "d3e6":
        return r"\((3,6,m),\ 4\le m\le9\)"
    if name == "d3e6m4":
        return r"\((3,6,4)\)"
    if name == "d3e6m5to9":
        return r"\((3,6,m),\ 5\le m\le9\)"
    if name.startswith("d4e6m"):
        m = int(name[-1])
        if square:
            return (r"\((4,6,3),\ C_\nu=7/3\)" if m == 3 else rf"\((4,6,{m})\)")
        caps = {3: "5/2", 4: "2", 5: "3/2", 6: "1", 7: "1/2"}
        return rf"\((4,6,{m}),\ C_\nu={caps[m]}\)"
    return rf"\(d={name[1]},e={name[3]}\)"


def validate(data_dir):
    files = {"tables.json", "flags.json", "queries.json", "arithmetic.json"}
    require({p.name for p in data_dir.iterdir()} == files, "data directory has missing/undeclared files")
    tables = read_json(data_dir / "tables.json")
    keys(tables, ["schema", "parameters", "states", "generic", "square", "replacements"], "tables")
    require(tables["schema"] == SCHEMA, "wrong schema")
    require(tables["parameters"] == PARAMETERS, "wrong article parameters")
    states, exclusions = enumerate_states()
    names = {state_id(v): v for v in states}
    rows = indexed(tables["states"], "id", names, "states")
    for name, v in names.items():
        row = rows[name]
        keys(row, ["id", "d", "e", "m", "moment", "beta", "cap", "q", "terminal_j"], name)
        require(tuple(integer(row[k], name) for k in ("d", "e", "m")) == v, name + ": wrong state")
        require(integer(row["moment"], name) == states[v]["moment"], name + ": wrong moment")
        for k in ("beta", "cap", "q"):
            require(rational(row[k], name) == states[v][k], name + ": wrong " + k)
        rational(row["terminal_j"], name)
    require(len(states) == 28, "state census is not 28")
    edges = {}
    for v, r in states.items():
        edges[v] = []
        for w, s in states.items():
            if w[0] >= v[0] or w[1] > v[1]:
                continue
            cap = s["beta"]
            if w[0] == v[0]-1:
                cap = min(cap, F(v[0]*v[2]-2*r["moment"]-w[2], v[2]))
            if cap > 0:
                edges[v].append((w, cap))
    edge_count = sum(map(len, edges.values()))
    require(edge_count == 170, "edge census is not 170")

    @lru_cache(None)
    def all_flags(v):
        return ((v,),) + tuple((v,) + tail for w, _ in edges[v] for tail in all_flags(w))

    flags_payload = read_json(data_dir / "flags.json")
    keys(flags_payload, ["schema", "flags"], "flags")
    require(flags_payload["schema"] == SCHEMA and type(flags_payload["flags"]) is list, "flag schema")
    expected_flags = {tuple(map(state_id, flag)) for v in states for flag in all_flags(v)}
    actual_flags = set()
    for flag in flags_payload["flags"]:
        require(type(flag) is list and all(type(s) is str and s in names for s in flag), "invalid flag")
        key = tuple(flag)
        require(key not in actual_flags, "duplicate flag")
        actual_flags.add(key)
    require(actual_flags == expected_flags and len(actual_flags) == 739, "missing/extra flag coverage")

    @lru_cache(None)
    def flag_grid(flag, cap):
        """Every vertex, plus harmless feasible grid points; no floating LP."""
        caps = [cap]
        for v, w in zip(flag, flag[1:]):
            caps.append(min(caps[-1], dict(edges[v])[w]))
        levels = sorted({F(0), *caps})
        points = []
        for ascending in combinations_with_replacement(levels, len(flag)):
            yy = tuple(reversed(ascending))
            if all(y <= bound for y, bound in zip(yy, caps)):
                cost = states[flag[0]]["q"]*yy[0]
                cost += sum(((states[w]["q"]-states[v]["q"])*y
                             for v, w, y in zip(flag, flag[1:], yy[1:])), F(0))
                points.append((yy, cost))
        return tuple(points)

    @lru_cache(None)
    def vertex_J(v, baseline, cap):
        return max(cost-baseline*point[0]
                   for flag in all_flags(v) for point, cost in flag_grid(flag, cap))

    # A third, cap-level Bellman implementation checks the historical factorial
    # algorithm independently of the producer's support-function recurrence.
    base_levels = {F(0)}
    for v, row in states.items():
        base_levels.update((row["beta"], row["cap"]))
        base_levels.update(cap for _, cap in edges[v])

    @lru_cache(None)
    def bellman_tail(v, y):
        q = states[v]["q"]
        best = q*y
        for w, cap in edges[v]:
            upper = min(y, cap)
            levels = {level for level in base_levels if level <= upper} | {upper}
            for z in levels:
                best = max(best, q*(y-z)+bellman_tail(w, z))
        return best

    queries_payload = read_json(data_dir / "queries.json")
    keys(queries_payload, ["schema", "queries"], "queries")
    require(queries_payload["schema"] == SCHEMA, "query schema")
    domain = query_domain(states)
    queries = indexed(queries_payload["queries"], "id", [r[0] for r in domain], "queries")
    query_values, query_prices, counts = {}, {}, {"empty": 0, "point": 0, "positive_cap": 0}
    for name, v, baseline, cap, gamma, bmax in domain:
        r = queries[name]
        keys(r, ["id", "state", "baseline", "cap", "value", "flag", "point", "price"], name)
        require(r["state"] == state_id(v) and rational(r["baseline"], name) == baseline
                and rational(r["cap"], name) == cap, name + ": changed query domain")
        if cap < 0:
            require(all(r[k] is None for k in ("value", "flag", "point", "price")), name + ": empty domain")
            counts["empty"] += 1
            continue
        counts["point" if cap == 0 else "positive_cap"] += 1
        value = rational(r["value"], name)
        exact = vertex_J(v, baseline, cap)
        require(value == exact, name + ": altered scalar witness")
        require(type(r["flag"]) is list and tuple(r["flag"]) in expected_flags, name + ": invalid witness flag")
        flag = tuple(names[s] for s in r["flag"])
        require(flag[0] == v and type(r["point"]) is list and len(r["point"]) == len(flag), name + ": witness shape")
        point = tuple(rational(s, name) for s in r["point"])
        require(any(point == yy and cost-baseline*point[0] == value
                    for yy, cost in flag_grid(flag, cap)), name + ": invalid attaining point")
        levels = {level for level in base_levels if level <= cap} | {cap}
        bellman = max(bellman_tail(v, y)-baseline*y for y in levels)
        require(value == bellman, name + ": Bellman mismatch")
        query_values[name] = value
        if gamma is None:
            require(r["price"] is None, name + ": unexpected source price")
        else:
            price = (7-bmax)*F(25, 29)+bmax/gamma+value
            require(rational(r["price"], name) == price, name + ": wrong whole-boundary price")
            query_prices[name] = price
    for name in names:
        require(F(rows[name]["terminal_j"]) == query_values["terminal/"+name], name + ": wrong printed terminal J")
    terminal_max = max(query_values["terminal/"+name] for name in names)
    require(terminal_max < F(131, 100), "terminal J strict bound")

    grouped_receipts = {}
    for family, ids, column_names, denom in [
        ("generic", GROUPS, ["generic-lambda", "generic-mu"], 10000),
        ("square", SQUARE_GROUPS, ["square"], 100000),
    ]:
        table = indexed(tables[family], "id", ids, family)
        covered = set()
        receipts = []
        for ident in ids:
            row = table[ident]
            keys(row, ["id", "label", "numerators"], family + "/" + ident)
            require(row["label"] == expected_label(ident, family == "square"), "changed group label " + ident)
            require(type(row["numerators"]) is list and len(row["numerators"]) == len(column_names), "group columns")
            members = [v for v in states if group_for(v, family == "square") == ident]
            require(bool(members), "empty required group")
            covered.update(members)
            for column, numerator in zip(column_names, row["numerators"]):
                integer(numerator, ident)
                value = max(query_prices[column+"/"+state_id(v)] for v in members)
                require(value < F(numerator, denom), ident + ": failed strict ceiling")
                require(numerator == (value*denom).numerator//(value*denom).denominator+1,
                        ident + ": ceiling must be strict floor-plus-one")
                receipts.append({"group": ident, "column": column, "maximum": str(value),
                                 "strict_ceiling": f"{numerator}/{denom}"})
        require(covered == set(states), family + ": incomplete group union")
        grouped_receipts[family] = receipts
    replacements = indexed(tables["replacements"], "m", range(3, 7), "replacements")
    for m, cap in [(3, F(2)), (4, F(2)), (5, F(1)), (6, F(1))]:
        row = replacements[m]
        keys(row, ["m", "cap", "numerators"], "replacement")
        require(rational(row["cap"], "replacement") == cap and len(row["numerators"]) == 2, "replacement cap")
        for column, numerator in zip(["lambda", "mu"], row["numerators"]):
            integer(numerator, "replacement")
            value = query_prices[f"distinct-{column}/4-6-{m}"]
            require(value < F(numerator, 10000), "replacement strict ceiling")
            require(numerator == (value*10000).numerator//(value*10000).denominator+1, "replacement rounding")

    # Exact arithmetic identities below are derived from the paper formulas;
    # their asserted values are delivered records, not hardcoded answers.
    lam, mu, theta = (F(PARAMETERS[k]) for k in ("lambda", "mu", "theta"))
    def price(scenario, v):
        return query_prices[scenario+"/"+state_id(v)]
    def maximum(scenario, predicate=lambda v: True):
        return max(price(scenario, v) for v in states if predicate(v))
    def phi(eta, z, a, b):
        return max(F(0), b-eta*z)**5-max(F(0), a-eta*z)**5
    def kstar(exponent, eta, a, b):
        return F(exponent, exponent-1)*phi(eta, 1+F(1, exponent), a, b)+F(exponent-2, exponent-1)*phi(eta, F(2), a, b)
    computed = {
        "initial_mass": 1-2*lam**5, "square_mass": 1-2*F(87, 100)**5,
        "terminal_gap": 8-(F(75, 29)+F(200, 49)+F(131, 100)),
        "generic_lambda_cubic": price("generic-lambda", (4, 6, 3)),
        "generic_mu_cubic": price("generic-mu", (4, 6, 3)),
        "generic_lambda_quartic": price("generic-lambda", (4, 6, 4)),
        "generic_mu_quartic": price("generic-mu", (4, 6, 4)),
        "generic_lambda_quintic": price("generic-lambda", (4, 6, 5)),
        "generic_lambda_quartic_threefold": price("generic-lambda", (3, 6, 4)),
        "low_b_price": maximum("low-b"), "quartic_low_price": query_prices["quartic-low"],
        "quartic_retained_price": F(75, 29)+F(8, 3)/lam+F(4, 3),
        "square_maximum": maximum("square"), "square_quartic": price("square", (4, 6, 4)),
        "square_quintic": price("square", (4, 6, 5)),
        "square_threefolds": maximum("square", lambda v: v[:2] == (3, 6)),
        "square_gap": 8-maximum("square"),
        "factorial_quintic_price": price("factorial-lambda", (4, 6, 5)),
        "factorial_shell_loss": 2*(mu-F(7, 10))**5,
        "noncartier_quartic_price": query_prices["noncartier-quartic"],
        "distinct_cubic_lambda_price": price("distinct-lambda", (4, 6, 3)),
        "singular_septic_lambda_price": price("singular-lambda", (4, 6, 7)),
        "smooth_large_mu_price": maximum("smooth-7-mu"),
        "smooth_theta_price": maximum("smooth-3-theta"),
    }
    order_conditions = []
    for m in (3, 4, 5):
        q = states[4, 6, m]["q"]
        p = query_values[f"affine/{m}"]
        computed[f"affine_{m}"] = p
        A = q-F(25, 29)
        for label, eps, K in [("generic", F(1, 500), F(m, 2)),
                              ("power", F(0), F(2) if m == 3 else F(m, 2))]:
            C = 8-eps-F(175, 29)-p
            require(C > (4-K)*A > 0 and 4*A > C and q*lam > 1, "fixed-order denominator/monotonicity")
            value = (q*lam-1)/(4*A-C)
            computed[f"{label}_order_{m}"] = value
            lower = F(33, 100) if label == "generic" or m == 3 else (F(7, 20) if m == 4 else F(17, 50))
            require(value > lower > F(1, 10), "full-source fixed order and dilution")
            order_conditions.append({"family": label, "m": m, "lower": str(lower), "excess": str(value-lower)})
            if m == 3 and label == "generic":
                computed["mu_order_3"] = (q*mu-1)/(4*A-C)
                require(computed["mu_order_3"] > F(49, 100), "mu order")
    computed["power_order_5_excess"] = computed["power_order_5"]-F(17, 50)
    c0 = computed["initial_mass"]
    computed["integral_first_loss"] = 2*phi(F(1), F(99, 200), lam, mu)
    computed["integral_second_loss"] = 2*phi(F(1), F(147, 200), mu, F(49, 50))
    computed["integral_reserve"] = c0-computed["integral_first_loss"]-computed["integral_second_loss"]
    computed["cubic_reserve"] = c0-phi(F(33, 100), F(1), lam, mu)-phi(F(33, 100), F(2), lam, mu)
    computed["cartier_reserve"] = c0-2*phi(F(7, 20), F(2), lam, mu)
    computed["singular_quartic_reserve"] = c0-phi(F(7, 20), F(1), lam, mu)-phi(F(7, 20), F(3), lam, mu)
    computed["singular_quintic_reserve"] = c0-phi(F(17, 50), F(1), lam, mu)-phi(F(17, 50), F(4), lam, mu)
    computed["smooth_first_loss"] = kstar(6, F(33, 100), lam, mu)
    computed["smooth_second_loss"] = kstar(6, F(49, 100), mu, theta)
    computed["smooth_reserve"] = c0-computed["smooth_first_loss"]-computed["smooth_second_loss"]
    arithmetic = read_json(data_dir / "arithmetic.json")
    keys(arithmetic, ["schema", "claims"], "arithmetic")
    require(arithmetic["schema"] == SCHEMA, "arithmetic schema")
    keys(arithmetic["claims"], computed, "arithmetic claims")
    for name, value in computed.items():
        require(rational(arithmetic["claims"][name], name) == value, name + ": altered rational identity")

    margin, error = F(1, 2000), F(1, 10000)
    gaps = {}
    paid_exceptions = {(3, 6, 4), (4, 6, 3), (4, 6, 4), (4, 6, 5)}
    payments = [
        ("generic-lambda", lambda v: v not in paid_exceptions),
        ("generic-mu", lambda v: v != (4, 6, 3)),
        ("square", lambda v: True), ("low-b", lambda v: True),
        ("factorial-mu", lambda v: v != (4, 6, 3)),
        ("factorial-lambda", lambda v: v not in {(3, 6, 4), (4, 6, 3), (4, 6, 4)}),
        ("distinct-mu", lambda v: v != (4, 6, 7)),
        ("distinct-lambda", lambda v: v not in {(3, 6, 4), (4, 6, 3), (4, 6, 4), (4, 6, 7)}),
        ("singular-mu", lambda v: True),
        ("singular-lambda", lambda v: v not in paid_exceptions),
        ("smooth-7-mu", lambda v: True), ("smooth-3-theta", lambda v: True),
    ]
    for scenario, allowed in payments:
        gaps[scenario] = 8-maximum(scenario, allowed)
        require(gaps[scenario] > margin+error, scenario + ": insufficient cutting margin")
    for scenario in ("factorial-lambda", "distinct-lambda", "singular-lambda"):
        allowed = lambda v: v != (3, 6, 4) and not (scenario == "factorial-lambda" and v == (4, 6, 3))
        gap = 9-F(25, 29)-maximum(scenario, allowed)
        require(gap > F(3, 100), "n>=8 price increment")
        gaps[scenario+"-n8"] = gap
    for key in ("quartic_low_price", "noncartier_quartic_price"):
        require(8-computed[key] > margin+error, key + ": insufficient strict margin")
    require(F(1, 100)+F(1, 1000) < computed["terminal_gap"], "terminal error")
    require(computed["quartic_retained_price"] < 7+F(1, 1000), "quartic retained price")
    require(F(1, 1000)+F(3009, 7994) < F(49267, 107919), "quartic affine bridge")
    require(F(1, 500) > margin+error, "fixed order epsilon")
    require(F(1, 50) < c0 and computed["factorial_shell_loss"] < F(1, 50), "factorial mass")
    require(c0-computed["integral_first_loss"] > F(1, 100), "first integral reserve")
    require(computed["integral_second_loss"] < F(1, 500), "second integral loss")
    for name, lower in [("integral_reserve", F(1, 100)), ("cubic_reserve", F(3, 10000)),
                        ("cartier_reserve", F(1, 100)), ("singular_quartic_reserve", F(4, 1000)),
                        ("singular_quintic_reserve", F(3, 1000)), ("smooth_reserve", F(1, 500))]:
        require(computed[name] > lower > F(PARAMETERS["intrinsic_reserve"]), "source reserve: " + name)
    require(computed["square_mass"] > 0 and 2-2*F(49, 50)**5 > 0, "original Hilbert source positivity")
    require(1-F(25, 29) == F(4, 29), "ambient dimension monotonicity")

    # These exact checks support the analytic all-exponent proof in the article;
    # this finite list is expressly NOT a proof of geometric or infinite scope.
    characters, endpoints = [], []
    for exponent in range(3, 8):
        for k in range(1, exponent):
            value = 3-F(k*(exponent-k), exponent)
            require(value <= 2+F(1, exponent), "character cap")
            characters.append({"exponent": exponent, "character": k, "cap": str(value)})
        if exponent <= 6:
            for eta, a, b in [(F(33, 100), lam, mu), (F(49, 100), mu, theta)]:
                val = kstar(exponent, eta, a, b)
                require(val <= kstar(6, eta, a, b) and b < 4*eta, "analytic endpoint conditions")
                require(F(exponent, exponent-1)*(1+F(1, exponent))+F(exponent-2, exponent-1)*2 == 3,
                        "balanced convexity weights")
                endpoints.append({"exponent": exponent, "eta": str(eta), "interval": [str(a), str(b)], "loss": str(val)})
    for m, intrinsic, ambient in [(8, F(1, 4), F(1, 2)), (9, F(1, 9), F(1, 3))]:
        require(states[3, 6, m]["beta"] == intrinsic, "embedding-six high-multiplicity cap")
        require(F(3)-F(2*moment(4, m), m) == ambient, "separate embedding-seven cap")
        v = (3, 6, m)
        require(vertex_J(v, states[v]["q"], intrinsic) == 0, "complete high-multiplicity tail")
        require(states[v]["q"] <= (F(2) if m == 8 else F(9, 4)), "high-multiplicity whole-tail slope")

    view_dir = data_dir.parent / "views"
    rendered = views(tables, arithmetic)
    require({p.name for p in view_dir.iterdir()} == set(rendered), "missing/extra rendered views")
    for name, text in rendered.items():
        require((view_dir / name).read_text() == text, "rendered view differs from canonical data: " + name)
    return {
        "status": "PASS", "schema": SCHEMA,
        "families": {
            "16-double-twentynine/states": {"states": len(states), "edges": edge_count, "flags": len(actual_flags),
                "terminal_entries": len(rows), "maximum_terminal_J": str(terminal_max), "first_exclusions": exclusions},
            "16-double-twentynine/prices": {"queries": len(queries), "domain_kinds": counts,
                "flag_grids": flag_grid.cache_info().currsize,
                "feasible_grid_points": sum(len(flag_grid(flag, cap)) for flag, cap in grid_keys(states, domain, all_flags)),
                "generic_entries": 26, "square_entries": 12, "replacement_entries": 8,
                "grouped": grouped_receipts, "cutting_gaps": {k: str(v) for k, v in gaps.items()}},
            "16-double-twentynine/arithmetic": {"identities": len(computed), "values": {k: str(v) for k, v in computed.items()},
                "fixed_order_conditions": order_conditions, "finite_character_checks": characters,
                "finite_endpoint_checks": endpoints},
        },
        "rendered_views": sorted(rendered),
        "scope": "Exact finite arithmetic conditional on the article's geometric, analytic all-parameter, and uniform-error proofs; not formal or whole-graph verification.",
    }


def grid_keys(states, domain, all_flags):
    return {(flag, cap) for _, v, _, cap, _, _ in domain if cap >= 0 for flag in all_flags(v)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = validate(args.data.resolve())
        code = 0
    except (Invalid, ValueError, TypeError, KeyError, IndexError, OSError) as exc:
        result = {"status": "FAIL", "error": str(exc)}
        code = 1
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "report": str(args.report),
                      "error": result.get("error")}), flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
