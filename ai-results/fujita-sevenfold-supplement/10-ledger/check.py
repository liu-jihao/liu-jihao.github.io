#!/usr/bin/env python3
"""Independent exact checker for the complete 10-ledger certificate.

Python >= 3.9; standard library only. No generator, provenance file, live fact
store, network, private path or cached acceptance is an input. See README.md
and the article's certificate checking lemma for the mathematical implication.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path
import re
import sys

R, U, Q, H0 = F(9, 10), F(27, 25), F(13, 5), F(131, 1000)
DEN = 10**9
TARGET = F(9999, 10000)
CONSTANTS = {"old_anchor": "9/10", "new_anchor": "27/25", "endpoint": "13/5", "old_height": "131/1000", "rounding_denominator": DEN, "bounded_target": "9999/10000", "comparison_error": "1/1000000"}
LINE_COUNTS = [1,1,1,1,1,1,1,3,4,1,4,3,1,1,4,2,1,1,1,1,3,3,4,1,1,1,1,4,4,3,3,7,4,5,5,3,3,3,5,5,4,5,3,3,4,4,3,3,4,4,7,6,6,4,7,6,6,7,6,5,6,5,3,2,1,7,7,7,6,6,4]
SPECS = {
    "P2": {"height": F(24, 125), "grid1": [F(8, 45), F(1), F(2), F(4), F(5)], "leaves": 44, "splits": 24, "cells": 340, "maximum": 999872274, "maximum_tags": ["2.4"], "tail": F(4191719943087, 476837158203125)},
    "T": {"height": F(103, 500), "grid1": [F(103, 540), F(1), F(2), F(3), F(4)], "leaves": 82, "splits": 62, "cells": 717, "maximum": 999778839, "maximum_tags": ["2.400"], "tail": F(3622557586593623, 427246093750000000)},
}


class CertificateError(ValueError):
    pass


class Audit:
    def __init__(self):
        self.context = "input"
        self.checks = Counter()
        self.families = {}

    def require(self, condition, message):
        self.checks[message] += 1
        if not condition:
            raise CertificateError(self.context + ": " + message)

    def keys(self, value, expected):
        self.require(type(value) is dict and set(value) == set(expected), "exact object fields " + ",".join(expected))

    def array(self, value, count=None):
        self.require(type(value) is list and (count is None or len(value) == count), "array shape/count")
        return value

    def integer(self, value, low=0, high=10**12):
        self.require(type(value) is int and low <= value <= high, "integer type/range")
        return value

    def rational(self, value):
        self.require(type(value) is str and len(value) <= 250 and re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?", value) is not None, "rational syntax")
        result = F(value)
        self.require(str(result) == value, "canonical reduced rational")
        return result

    def tag(self, value):
        self.require(type(value) is str and re.fullmatch(r"[1-5]\.[1-4][01]{0,30}", value) is not None, "tag type/range")
        return value


def load_json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise CertificateError("duplicate JSON field: " + key)
            result[key] = value
        return result

    def reject_float(value):
        raise CertificateError("floating/nonfinite JSON number: " + value)

    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique,
                      parse_float=reject_float, parse_constant=reject_float)


def minimum_moment(codimension, multiplicity):
    if codimension == 0:
        return 0
    result, degree = 0, 1
    while comb(codimension + degree - 1, codimension) < multiplicity:
        result += multiplicity - comb(codimension + degree - 1, codimension)
        degree += 1
    return result


def envelope(row, x):
    return min(a*x+p for a, p in row["lines"])


def endpoints(row, cap):
    # This also defines the point case cap=0; parallel lines add no breakpoint.
    result = {F(0), cap}
    for i, (a, p) in enumerate(row["lines"]):
        for aa, pp in row["lines"][i+1:]:
            if a != aa:
                x = (pp-p)/(a-aa)
                if 0 < x < cap:
                    result.add(x)
    return sorted(result)


def audit_scalar(raw, audit):
    audit.context = "10-ledger/scalar"
    audit.keys(raw, ("schema", "rows"))
    audit.require(raw["schema"] == "10-ledger/scalar/v1", "scalar schema")
    rows = []
    for row in audit.array(raw["rows"], 71):
        audit.keys(row, ("id", "d", "e", "m", "moment", "cap", "radius", "lines"))
        item = {key: audit.integer(row[key], 0, 1000) for key in ("id", "d", "e", "m", "moment")}
        audit.require(1 <= item["id"] <= 71, "scalar row id range")
        item.update(cap=audit.rational(row["cap"]), radius=audit.rational(row["radius"]))
        item["lines"] = [tuple(audit.rational(x) for x in audit.array(pair, 2)) for pair in audit.array(row["lines"], LINE_COUNTS[item["id"]-1])]
        audit.require(len(set(item["lines"])) == len(item["lines"]), "no duplicate affine witnesses")
        audit.require(all(a > 0 and p >= 0 for a, p in item["lines"]), "affine sign conditions")
        rows.append(item)
    audit.require([r["id"] for r in rows] == list(range(1, 72)), "all scalar ids once in order")
    expected = {(1, 1, 1)}
    rejected = []
    for d in range(2, 7):
        for e in range(d, 8):
            c = e-d
            if c == 0:
                expected.add((d, e, 1))
                continue
            m = c+1
            while 2*minimum_moment(c, m) < d*m:
                expected.add((d, e, m))
                m += 1
            rejected.append({"d": d, "e": e, "first_excluded_m": m, "moment": minimum_moment(c, m)})
    audit.require([(r["d"], r["e"], r["m"]) for r in rows] == sorted(expected), "exhaustive scalar state enumeration")
    records = []
    for row in rows:
        audit.context = "scalar row " + str(row["id"])
        d, e, m = row["d"], row["e"], row["m"]
        audit.require(row["moment"] == minimum_moment(e-d, m), "exact Hilbert moment")
        audit.require(row["cap"] == d-2*F(row["moment"], m) > 0, "exact positive deficit cap")
        radius = row["radius"]
        audit.require(radius > 0 and (1000*radius).denominator == 1 and radius**d >= m and (radius-F(1,1000))**d < m, "least upward thousandth radius")
        for a, p in row["lines"]:
            extra, witness, admissible = F(0), None, 0
            for successor in rows:
                if successor["d"] >= d or successor["e"] > e:
                    continue
                cap = min(row["cap"], successor["cap"])
                if successor["d"] == d-1:
                    cap = min(cap, row["cap"]-F(successor["m"], m))
                if cap <= 0:
                    continue
                admissible += 1
                for y in endpoints(successor, cap):
                    surplus = envelope(successor, y)-max(a, radius)*y
                    if surplus > extra:
                        extra = surplus
                        witness = {"successor": successor["id"], "deficit": str(y)}
            rhs = max(radius-a, F(0))*row["cap"]+extra
            audit.require(p >= rhs, "complete successor recurrence")
            records.append({"row": row["id"], "line": [str(a), str(p)], "rhs": str(rhs), "slack": str(p-rhs), "positive_cap_successors": admissible, "largest_positive_surplus_at": witness})
    audit.families["10-ledger/scalar"] = {"states": len(rows), "states_by_dimension": [sum(r["d"] == d for r in rows) for d in range(1,7)], "affine_witnesses": len(records), "first_excluded_multiplicities": rejected, "recurrence_checks": records, "interface": "all degree-floor-one descendants, including adjacent and skipped dimensions and the direct point successor"}


def make_integrator(audit):
    # Rational interpolation, adapted from the written audit, not its producer.
    nodes = [F(j, 6) for j in range(7)]
    weights = []
    for j, x in enumerate(nodes):
        polynomial = [F(1)]
        for k, y in enumerate(nodes):
            if k == j:
                continue
            following = [F(0)]*(len(polynomial)+1)
            for i, coefficient in enumerate(polynomial):
                following[i] -= y*coefficient/(x-y)
                following[i+1] += coefficient/(x-y)
            polynomial = following
        weights.append(sum(coefficient/F(i+1) for i, coefficient in enumerate(polynomial)))
    for degree in range(7):
        audit.require(sum(w*x**degree for x, w in zip(nodes, weights)) == F(1, degree+1), "exact integration moment")

    def integrate(function, left, right):
        return (right-left)*sum(w*function(left+x*(right-left)) for x, w in zip(nodes, weights))
    return integrate


def cell_function(name, left, right, bounds, kind, height, old_tail, audit):
    audit.require(type(name) is str, "polynomial name type")
    if name in ("B1", "B61/100", "B11/25", "B17/50"):
        rate = F(name[1:])
        # For t>=29/25 the high baseline decreases from (1617/2117)^6.
        # 17/50 is a conservative bound there. All other branches are constant.
        cuts = sorted({left, right} | {t for t in (R, F(99,100), F(26,25), F(29,25)) if left < t < right})
        for l, r in zip(cuts, cuts[1:]):
            t = (l+r)/2
            baseline = F(1) if t < R else F(61,100) if t < F(99,100) else F(11,25) if t < F(26,25) else F(17,50)
            audit.require(rate >= baseline, "constant baseline dominates on whole subinterval")
        return lambda t: 7*rate*t**6
    if name == "Bhigh":
        audit.require(F(29,25) <= left < right <= F(125,52), "high baseline branch")
        return lambda t: 7*((125-52*t)/73)**6
    if name == "Bzero":
        audit.require(left >= F(125,52), "zero baseline branch")
        return lambda t: F(0)
    allowed = {"0:N5lo", "0:N5hi", "0:C", "0:E1", "1:E2", "1:Tlow"}
    audit.require(name in allowed and (kind == "T" or name != "1:Tlow"), "known polynomial and family")
    index = int(name[0])
    v, h = (R, H0) if index == 0 else (U, height)
    if old_tail:
        audit.require(index == 0 and left >= R, "old tail uses only old profile after anchor")
        slope = F(16)
    else:
        audit.require(not left < v < right, "selected slope does not change inside cell")
        slope = bounds[2*index+1] if right <= v else bounds[2*index]
    offset = slope*v-h
    audit.require(offset >= 0 and slope*left-offset >= 0 and slope*right-offset > 0, "unclipped profile positive on open cell")
    cmin = slope-offset/left if left else slope
    cmax = slope-offset/right
    audit.require(cmin <= cmax, "normalized profile monotone on cell")
    if name == "0:N5lo":
        audit.require(0 <= cmin <= cmax <= F(1,7), "N5lo branch")
        return lambda t: 7*t**6-7*F(35,6)**5*t*(slope*t-offset)**5
    if name == "0:N5hi":
        audit.require(F(1,7) <= cmin <= cmax <= F(1,2), "N5hi branch")
        return lambda t: 7*t**6-7/F(36)**5*t*((29+7*slope)*t-7*offset)**5
    if name == "0:C":
        audit.require(F(1,7) <= cmin <= cmax <= F(1,2), "convex branch")
        needed = [1]
    elif name == "0:E1":
        audit.require(cmin >= F(1,2), "old E1 branch")
        needed = [1]
    elif name == "1:E2":
        audit.require(cmin >= (F(4,25) if kind == "T" else 0), "new E2 branch")
        needed = [2]
    else:
        audit.require(0 <= cmin <= cmax <= F(4,25), "Tlow branch")
        needed = [2, 3]
    signs = {}
    for j in needed:
        lvalue = (1-j*slope)*left+j*offset
        rvalue = (1-j*slope)*right+j*offset
        audit.require(not min(lvalue,rvalue) < 0 < max(lvalue,rvalue), "no positive-part zero in cell interior")
        signs[j] = lvalue+rvalue > 0

    def p(j, t):
        return 7*((1-j*slope)*t+j*offset)**6 if signs[j] else F(0)
    if name == "0:C":
        return lambda t: F(7,6)**6*p(1,t)
    if name == "1:Tlow":
        return lambda t: 3*p(2,t)-2*p(3,t)
    return lambda t: p(needed[0], t)


def audit_cells(cells, bounds, kind, height, printed_total, integrate, audit, old_tail=False):
    audit.array(cells)
    audit.require(0 < len(cells) <= 100, "nonempty bounded cell count")
    previous, exact_sum, upper_sum = F(0), F(0), 0
    records = []
    base_context = audit.context
    for i, cell in enumerate(cells):
        audit.context = base_context + " cell " + str(i+1)
        audit.keys(cell, ("l", "r", "name", "upper"))
        left, right = audit.rational(cell["l"]), audit.rational(cell["r"])
        numerator = audit.integer(cell["upper"], 1, DEN)
        audit.require(left == previous and 0 <= left < right <= Q, "complete ordered cell partition without gap/overlap")
        fn = cell_function(cell["name"], left, right, bounds, kind, height, old_tail, audit)
        exact = integrate(fn, left, right)
        audit.require(0 <= exact < F(numerator, DEN), "strict exact cell upper bound")
        audit.require(numerator == (DEN*exact).numerator//(DEN*exact).denominator+1, "floor-plus-one upward rounding")
        records.append({"l": str(left), "r": str(right), "name": cell["name"], "exact": str(exact), "upper": numerator, "strict_gap": str(F(numerator,DEN)-exact)})
        previous, exact_sum, upper_sum = right, exact_sum+exact, upper_sum+numerator
    audit.context = base_context
    audit.require(previous == Q, "cell partition reaches full order endpoint")
    audit.require(upper_sum == printed_total and exact_sum < F(printed_total,DEN), "all cell bounds consumed in strict rectangle total")
    return {"exact_total": str(exact_sum), "upper": printed_total, "rounding_gap": str(F(printed_total,DEN)-exact_sum), "cells": records}


def audit_profile_family(family, kind, integrate, audit):
    spec = SPECS[kind]
    audit.context = "10-ledger/" + kind
    audit.keys(family, ("id", "height", "grid0", "grid1", "splits", "rectangles", "new_tail"))
    audit.require(family["id"] == kind, "family identity")
    height = audit.rational(family["height"])
    audit.require(height == spec["height"], "actual source-specific height")
    grid0 = [audit.rational(x) for x in audit.array(family["grid0"], 6)]
    grid1 = [audit.rational(x) for x in audit.array(family["grid1"], 5)]
    audit.require(grid0 == [H0/R,F(1),F(2),F(4),F(8),F(16)] and grid1 == spec["grid1"] and grid1[0] == height/U, "entire prescribed slope grids")
    audit.require(all(a < b for grid in (grid0,grid1) for a,b in zip(grid,grid[1:])), "strictly increasing grids")
    active = {f"{i+1}.{j+1}": [a,b,c,d] for i,(a,b) in enumerate(zip(grid0,grid0[1:])) for j,(c,d) in enumerate(zip(grid1,grid1[1:]))}
    for split in audit.array(family["splits"], spec["splits"]):
        audit.keys(split, ("tag", "coordinate", "midpoint"))
        tag, coordinate = audit.tag(split["tag"]), audit.integer(split["coordinate"], 0, 1)
        midpoint = audit.rational(split["midpoint"])
        audit.require(tag in active, "split existing active leaf exactly once")
        bounds = active.pop(tag)
        i = 2*coordinate
        audit.require(bounds[i] < midpoint < bounds[i+1] and midpoint == (bounds[i]+bounds[i+1])/2, "nondegenerate exact bisection")
        low, high = bounds.copy(), bounds.copy()
        low[i+1], high[i] = midpoint, midpoint
        audit.require(tag+"0" not in active and tag+"1" not in active, "unique split children")
        active[tag+"0"], active[tag+"1"] = low, high
    rectangles = audit.array(family["rectangles"], spec["leaves"])
    tags, checked, total_cells = [], {}, 0
    for rectangle in rectangles:
        audit.keys(rectangle, ("tag", "bounds", "upper", "cells"))
        tag = audit.tag(rectangle["tag"])
        audit.context = "10-ledger/" + kind + " rectangle " + tag
        audit.require(tag not in checked and tag in active, "each terminal rectangle exactly once")
        bounds = [audit.rational(x) for x in audit.array(rectangle["bounds"], 4)]
        audit.require(bounds == active[tag], "all terminal bounds agree with split tree")
        upper = audit.integer(rectangle["upper"], 1, DEN)
        audit.require(F(upper,DEN) < TARGET, "strict bounded-domain target")
        checked[tag] = audit_cells(rectangle["cells"], bounds, kind, height, upper, integrate, audit)
        total_cells += len(rectangle["cells"])
        tags.append(tag)
    audit.require(tags == sorted(active), "complete ordered terminal tag set")
    audit.require(total_cells == spec["cells"], "all expected bounded integration cells")
    maximum = max(item["upper"] for item in checked.values())
    maximum_tags = [tag for tag in tags if checked[tag]["upper"] == maximum]
    audit.require(maximum == spec["maximum"] and maximum_tags == spec["maximum_tags"], "exact reported worst rectangle")
    audit.families["10-ledger/" + kind] = {"root": [str(grid0[0]),str(grid0[-1]),str(grid1[0]),str(grid1[-1])], "initial_rectangles": 20, "splits": len(family["splits"]), "leaves": len(checked), "cells": total_cells, "maximum_upper_numerator": maximum, "maximum_tags": maximum_tags, "target": str(TARGET), "margin_below_target": str(TARGET-F(maximum,DEN)), "volume_lower_margin": str(1-F(maximum,DEN)), "rectangles": checked}
    audit.context = "10-ledger/" + kind + "-new-tail"
    tail = family["new_tail"]
    audit.keys(tail, ("slope_min", "denominator", "loss", "target", "volume_floor"))
    a, denominator, loss, target, floor = [audit.rational(tail[k]) for k in ("slope_min", "denominator", "loss", "target", "volume_floor")]
    audit.require(a == grid1[-1] and denominator == 2*a-1 > 0, "entire new-slope tail and positive denominator")
    audit.require(a*U-height >= 0 and height/U > F(4,25) and U-2*height > 0, "new-tail E2 branch on every order after anchor")
    zero = 2*(a*U-height)/(2*a-1)
    audit.require(zero > U, "tail positive segment followed by identically zero ray")
    exact = integrate(lambda t: 7*((1-2*a)*t+2*(a*U-height))**6, U, zero)
    audit.require(exact == (U-2*height)**7/denominator == loss == spec["tail"], "actual prose tail fraction and exact antiderivative")
    audit.require(target == F(1,100) and floor == F(1,32) and loss < target < floor, "strict tail loss and imported volume margins")
    audit.families["10-ledger/" + kind + "-new-tail"] = {"slope_domain": ">="+str(a), "normalized_height_lower": str(height/U), "branch_margin": str(height/U-F(4,25)), "denominator_lower": str(denominator), "zero_at_minimum_slope": str(zero), "loss_upper": str(loss), "margin_below_one_hundredth": str(target-loss), "volume_lower_margin": str(floor-loss), "interface": "all a1 in the unbounded ray, all t>=27/25; denominator increases with a1"}


def audit_profiles(raw, audit):
    audit.context = "profiles"
    audit.keys(raw, ("schema", "constants", "families", "old_tail"))
    audit.require(raw["schema"] == "10-ledger/profiles/v1", "profile schema")
    audit.keys(raw["constants"], tuple(CONSTANTS))
    audit.require(raw["constants"] == CONSTANTS and type(raw["constants"]["rounding_denominator"]) is int, "all fixed profile constants")
    audit.require(((125-52*F(29,25))/(73*F(29,25)))**6 < F(17,50), "constant baseline controls decreasing high branch")
    prices = [("old",H0,F(5,38)), ("P2",SPECS["P2"]["height"],F(27,140)), ("T",SPECS["T"]["height"],F(124,601))]
    margin_records = []
    for name, height, available in prices:
        error = F(raw["constants"]["comparison_error"])
        audit.require(height+error < available, "source-specific height absorbs uniform comparison error")
        margin_records.append({"profile": name, "available": str(available), "height": str(height), "error": str(error), "strict_margin": str(available-height-error)})
    audit.require(R/(7-F(4,25)) == F(5,38) and U/(7-F(7,5)) == F(27,140) and F(124,601) < F(26,125), "actual old/divisor/fivefold source-price arithmetic")
    audit.families["10-ledger/source-margins"] = {"comparisons": margin_records}
    integrate = make_integrator(audit)
    for family, kind in zip(audit.array(raw["families"], 2), ("P2", "T")):
        audit_profile_family(family, kind, integrate, audit)
    audit.context = "10-ledger/old-tail"
    tail = raw["old_tail"]
    audit.keys(tail, ("slope_min", "upper", "cells"))
    audit.require(audit.rational(tail["slope_min"]) == 16, "old unbounded slope ray")
    upper = audit.integer(tail["upper"], 1, DEN)
    audit.require(upper == 502032677 and F(upper,DEN) < 1, "exact old-tail total and strict margin")
    audit.array(tail["cells"], 5)
    checked = audit_cells(tail["cells"], None, "P2", SPECS["P2"]["height"], upper, integrate, audit, old_tail=True)
    checked.update(slope_domain=">=16", volume_lower_margin=str(1-F(upper,DEN)), interface="all a1; rate one before 9/10 and old slope sixteen thereafter")
    audit.families["10-ledger/old-tail"] = checked


def validate(data, audit):
    audit.require(data.is_dir(), "data directory exists")
    entries = list(data.iterdir())
    audit.require({p.name for p in entries} == {"scalar.json", "profiles.json"} and all(p.is_file() and not p.is_symlink() for p in entries), "exact input file set")
    audit_scalar(load_json(data / "scalar.json"), audit)
    audit_profiles(load_json(data / "profiles.json"), audit)
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(entries)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data"))
    parser.add_argument("--report", type=Path, default=Path("report.json"))
    args = parser.parse_args()
    audit = Audit()
    report = {"checker": "10-ledger/v1", "status": "FAIL"}
    try:
        report["data_sha256"] = validate(args.data, audit)
        report["status"] = "PASS"
    except Exception as error:
        report["error"] = type(error).__name__ + ": " + str(error)
    report["checks"] = dict(sorted(audit.checks.items()))
    report["families"] = audit.families
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary = {"status": report["status"], "checks": sum(audit.checks.values()), "families": list(audit.families)}
    if "error" in report:
        summary["error"] = report["error"]
    print(json.dumps(summary, sort_keys=True))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
