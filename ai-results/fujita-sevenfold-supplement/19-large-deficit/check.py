#!/usr/bin/env python3
"""Exact validator for 19-large-deficit/paired-profile-v1.

Python 3.9+ standard library only. No generator, private state, or stored
verdict is imported. On the Danus host execute only through danus compute.
The polynomial arithmetic is adapted from check_paired_certificate.py;
README.md and manifest.json identify its preserved original source.
"""

import argparse
import ast
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path
import re
import sys

FAMILY = "19-large-deficit/paired-profile-v1"
DATA_NAME = "paired-profile.txt"
r, v, T, R = Q(9, 10), Q(27, 25), Q(29, 25), Q(13, 5)
eps, h, H = Q(1, 10**6), Q(263, 2000), Q(27, 85) - Q(1, 10**6)
cm, cp = Q(1455087, 8192000), Q(116407, 655360)
DEN, LIMIT, MAX_TOTAL = 10**12, 999900000000, 999894402477
ZERO = Q(0)
RATIONAL = r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?"
ROW = re.compile(r"(" + RATIONAL + r"),(" + RATIONAL
                 + r"),([A-Za-z0-9/]+),([0-9]+)")
BOX = re.compile(
    r"Box ([0-9]+): a=\[(" + RATIONAL + r"),(" + RATIONAL
    + r")\], d=\[(" + RATIONAL + r"),(" + RATIONAL
    + r")\], eta=(" + RATIONAL + r"), d_lower=(" + RATIONAL + r")\."
)


def demand(condition, message):
    if not condition:
        raise ValueError(message)


def rational(value):
    demand(isinstance(value, str) and re.fullmatch(RATIONAL, value),
           "invalid exact rational: " + repr(value))
    return Q(value)


class Lines:
    """A total parser: every nonempty line must be consumed in its place."""

    def __init__(self, raw):
        demand(b"\r" not in raw and b"\x00" not in raw and raw.endswith(b"\n"),
               "data must be an LF-terminated UTF-8 text stream")
        self.lines = [(i, s) for i, s in enumerate(
            raw.decode("utf-8").splitlines(), 1) if s]
        self.index = 0

    def peek(self):
        return self.lines[self.index][1] if self.index < len(self.lines) else ""

    def take(self, expected=None, pattern=None):
        demand(self.index < len(self.lines), "unexpected end of certificate")
        number, text = self.lines[self.index]
        self.index += 1
        if expected is not None:
            demand(text == expected, f"line {number}: expected {expected!r}")
            return text
        if pattern is not None:
            match = re.fullmatch(pattern, text)
            demand(match is not None, f"line {number}: malformed record {text!r}")
            return match
        return text

    def end(self):
        demand(self.index == len(self.lines), "unrecognized or duplicate trailing record")


def read_dictionary(lines, name, count):
    match = lines.take(pattern=re.escape(name) + r": (.+)")
    entries = ast.literal_eval(match[1])
    demand(type(entries) is list and len(entries) == count,
           name + ": wrong number of entries")
    result = []
    for entry in entries:
        demand(type(entry) is tuple and len(entry) == 2,
               name + ": expected an intercept/slope pair")
        result.append(tuple(rational(x) for x in entry))
    return result


def read_rows(lines):
    rows = []
    while ROW.fullmatch(lines.peek()):
        match = lines.take(pattern=ROW)
        rows.append((rational(match[1]), rational(match[2]), match[3], int(match[4])))
    demand(rows, "missing integral rows")
    return rows


def parse(raw):
    lines = Lines(raw)
    lines.take("# Exact rational cell certificate")
    lines.take("Every cell integral upper numerator has denominator 10^12.")
    match = lines.take(pattern=r"Balanced cutoff (" + RATIONAL
                       + r"); endpoint cutoff (" + RATIONAL + r")\.")
    cutoffs = (rational(match[1]), rational(match[2]))
    match = lines.take(pattern=r"Exact branch differences: (" + RATIONAL
                       + r"), (" + RATIONAL + r")\.")
    differences = (rational(match[1]), rational(match[2]))
    chord = read_dictionary(lines, "Chord tangents (intercept,slope)", 6)
    call = read_dictionary(lines, "Call tangents (intercept,slope)", 9)
    boxes = {}
    order = []
    while lines.peek().startswith("Box "):
        match = lines.take(pattern=BOX)
        ident = int(match[1])
        demand(ident not in boxes, f"duplicate box ID {ident}")
        fields = tuple(rational(x) for x in match.groups()[1:])
        floor = lines.peek() == "Use incoming 31/32."
        if floor:
            lines.take("Use incoming 31/32.")
            incoming, in_sum = [], None
        else:
            lines.take("Incoming cells l,u,density,upper numerator:")
            incoming = read_rows(lines)
            match = lines.take(pattern=r"Incoming numerator sum ([0-9]+)\.")
            in_sum = int(match[1])
        lines.take("Outgoing cells l,u,density,upper numerator:")
        outgoing = read_rows(lines)
        match = lines.take(pattern=r"Outgoing numerator sum ([0-9]+); reserve numerator "
                           r"<=([0-9]+); total <=([0-9]+)<999900000000\.")
        boxes[ident] = {
            "bounds": fields[:4], "eta": fields[4], "dlo": fields[5],
            "floor": floor, "in": incoming, "out": outgoing,
            "in_sum": in_sum, "out_sum": int(match[1]),
            "reserve": int(match[2]), "total": int(match[3])
        }
        order.append(ident)
    demand(order == list(range(301)), "missing, out-of-range, or reordered box IDs")
    match = lines.take(pattern=r"Largest written total ([0-9]+)/10\^12\.")
    largest = int(match[1])
    lines.take("Unbounded a>=3: use r^7+1/16 and the following old-slope3 outgoing cells:")
    tail = read_rows(lines)
    match = lines.take(pattern=r"Total numerator ([0-9]+)<999900000000\.")
    a_total = int(match[1])
    match = lines.take(pattern=r"Unbounded d>=20 exact bound: (" + RATIONAL
                       + r")<(" + RATIONAL + r")\.")
    d_exact, d_upper = rational(match[1]), rational(match[2])
    lines.take("# Exact coverage tree")
    match = lines.take(pattern=r"Root rectangle: \[(" + RATIONAL + r"),(" + RATIONAL
                       + r")\] x\[(" + RATIONAL + r"),(" + RATIONAL + r")\]\.")
    root = tuple(rational(x) for x in match.groups())
    lines.take("Node [0,left,right] bisects the a interval at its midpoint; "
               "node [1,left,right] bisects the d interval.")
    lines.take("A quoted numeric leaf is the Box number in the full cell witness; "
               "E denotes the exact infeasible-box test.")
    match = lines.take(pattern=r"Leaf count ([0-9]+)\. Largest integer leaf ([0-9]+)\.")
    leaves, coarse_max = int(match[1]), int(match[2])
    match = lines.take(pattern=r"Full covered area (" + RATIONAL + r")\.")
    area = rational(match[1])
    match = lines.take(pattern=r"Unbounded a>=3 bound (" + RATIONAL
                       + r"); d>=20 bound (" + RATIONAL + r")\.")
    coarse_a, coarse_d = rational(match[1]), rational(match[2])
    tree = json.loads(lines.take())
    lines.end()
    return dict(boxes=boxes, cutoffs=cutoffs, differences=differences,
                chord=chord, call=call, largest=largest, tail=tail,
                a_total=a_total, d_exact=d_exact, d_upper=d_upper,
                root=root, leaves=leaves, coarse_max=coarse_max,
                area=area, coarse_a=coarse_a, coarse_d=coarse_d, tree=tree)


def add(*polys):
    result = [ZERO] * max(map(len, polys))
    for poly in polys:
        for j, value in enumerate(poly):
            result[j] += value
    return result


def scale(poly, scalar):
    return [scalar * value for value in poly]


def mul(left, right):
    result = [ZERO] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i + j] += x * y
    return result


def power(A, B, exponent):
    return [Q(comb(exponent, j)) * A**j * B**(exponent-j)
            for j in range(exponent + 1)]


def integrate(poly, lo, hi):
    return sum((value * (hi**(j+1) - lo**(j+1)) / (j+1)
                for j, value in enumerate(poly)), ZERO)


def positive_power(A, B, exponent, lo, hi):
    left, right = A*lo+B, A*hi+B
    demand(not (left < 0 < right or right < 0 < left),
           "unsplit positive-part zero")
    # A zero endpoint is included by continuity; a zero affine function
    # contributes the zero polynomial.
    return power(A, B, exponent) if left+right > 0 else [ZERO] * (exponent+1)


def baseline(t):
    if t < r:
        return Q(1)
    if t < Q(99, 100):
        return Q(61, 100)
    if t < Q(26, 25):
        return Q(11, 25)
    return Q(17, 50)


def polynomial(label, A, B, lo, hi, chord, call):
    t6 = [ZERO] * 6 + [Q(7)]
    if label in ("B0", "B1", "B61/100", "B11/25", "B17/50"):
        beta = Q(1) if label == "B0" else Q(label[1:])
        cuts = sorted({lo, hi, *[x for x in (r, Q(99, 100), Q(26, 25))
                                 if lo < x < hi]})
        demand(all(beta >= baseline((left+right)/2)
                   for left, right in zip(cuts, cuts[1:])),
               "baseline density fails on a subinterval")
        return scale(t6, beta)
    demand(lo > 0, "a nonbaseline row cannot start at zero")
    demand(B <= 0, "profile ratio must be nondecreasing")
    cmin, cmax = A+B/lo, A+B/hi

    def domain(low, high=None):
        demand(cmin >= low and (high is None or cmax <= high),
               f"{label} outside its valid ratio interval: {cmin}, {cmax}")

    def pj(j):
        return scale(positive_power(1-j*A, -j*B, 6, lo, hi), Q(7))

    if label == "Ol":
        domain(0, Q(1, 7))
        return add(t6, scale([ZERO]+power(A, B, 5), -7*Q(35, 6)**5))
    if label == "Oh":
        domain(Q(1, 7), Q(1, 2))
        return add(t6, scale([ZERO]+power(29+7*A, 7*B, 5), -Q(7, 36**5)))
    if label == "Oc":
        domain(Q(1, 7), Q(1))
        return scale(pj(1), Q(7, 6)**6)
    if label == "Oe":
        domain(Q(1, 3))
        return pj(1)
    if re.fullmatch(r"Oq[0-5]", label):
        domain(Q(1, 6), Q(1, 4))
        C, D = chord[int(label[2:])]
        demand(min(C+D*cmin, C+D*cmax) > 0, "nonpositive chord tangent")
        return add(t6, scale([ZERO]+power(C+D*A, D*B, 5), -Q(7)))
    if re.fullmatch(r"Ok[0-8]", label):
        domain(Q(1, 6), Q(8, 25))
        C, D = call[int(label[2:])]
        demand(min(C+D*cmin, C+D*cmax) > 0, "nonpositive call tangent")
        return scale([ZERO]+power(C+D*A, D*B, 5), 7*Q(7, 6)**6)
    if label == "Nb":
        domain(0, cm)
        return scale(mul(power(1-2*A, -2*B, 5), [10*B, 1+10*A]), Q(7))
    if label == "Ne":
        domain(cp)
        return pj(1)
    if label == "Ns":
        domain(0)
        return add(pj(1), pj(2), scale(pj(3), -1))
    if label == "Z":
        domain(1)
        return [ZERO]
    raise ValueError("unknown density label: " + label)


def clipped_dimension(bounds):
    """Intersect the closed box with both linear feasibility inequalities."""
    al, au, dl, du = bounds
    poly = [(al, dl), (au, dl), (au, du), (al, du)]
    for value in (lambda p: v*p[1]-H,
                  lambda p: v*p[1]-(h+(v-r)*p[0]-eps)):
        result = []
        if not poly:
            return "empty"
        previous, fp = poly[-1], value(poly[-1])
        for current in poly:
            fc = value(current)
            if (fp >= 0) != (fc >= 0):
                ratio = fp/(fp-fc)
                result.append(tuple(x+ratio*(y-x)
                                    for x, y in zip(previous, current)))
            if fc >= 0:
                result.append(current)
            previous, fp = current, fc
        poly = list(dict.fromkeys(result))
    if not poly:
        return "empty"
    if len(poly) == 1:
        return "point"
    x0, y0 = poly[0]
    x1, y1 = poly[1]
    return "segment" if all((x1-x0)*(y-y0) == (y1-y0)*(x-x0)
                           for x, y in poly[2:]) else "area"


def validate(data):
    chord, call, boxes = data["chord"], data["call"], data["boxes"]
    demand(data["cutoffs"] == (cm, cp), "branch cutoff mismatch")
    aq = [Q(1, 6), Q(9, 50), Q(1, 5), Q(11, 50), Q(6, 25), Q(1, 4)]
    ak = [Q(1, 6), Q(9, 50), Q(1, 5), Q(11, 50), Q(6, 25),
          Q(13, 50), Q(7, 25), Q(3, 10), Q(8, 25)]
    demand(chord == [((29-28*c)/(36-42*c)-c*210/(36-42*c)**2,
                      210/(36-42*c)**2) for c in aq], "chord dictionary mismatch")
    demand(call == [((1-2*c)/(1-c)+c/(1-c)**2, -1/(1-c)**2)
                     for c in ak], "call dictionary mismatch")
    differences = tuple((1-2*c)**5*(1+10*c)-(1-c)**6 for c in (cm, cp))
    demand(data["differences"] == differences and differences[0] > 0 > differences[1],
           "incorrect branch crossing witnesses")
    demand(Q(5, 6)**6 < 1-Q(44, 51)**5, "chord endpoint comparison")
    demand(Q(7, 6)**6*Q(9, 17)**5 > Q(17, 25)**6, "call endpoint comparison")
    demand(Q(7, 6)**6 < Q(2**17, 6**6), "exact endpoint at one third")

    density_counts = Counter()
    clips = Counter()
    row_margins = []
    box_margins = []
    floor_count = 0
    finite_count = 0
    zero_integrals = 0
    zero_reserves = 0

    def rows_checked(rows, start, end, box=None):
        nonlocal zero_integrals
        demand(rows and rows[0][0] == start and rows[-1][1] == end,
               "row partition has missing endpoints")
        demand(all(left[1] == right[0] for left, right in zip(rows, rows[1:])),
               "row partition has a gap, overlap, duplicate, or missing row")
        total = 0
        for lo, hi, label, upper in rows:
            demand(start <= lo < hi <= end, "empty, reversed, or out-of-range row")
            if box is not None:
                demand(not lo < r < hi and not lo < v < hi,
                       "finite row crosses a profile anchor")
            if box is None:
                A, B = Q(3), h-3*r
                demand(label.startswith("O") or label == "Z", "non-old density in a tail")
            elif label.startswith("O") or label == "Z":
                A = box["bounds"][1] if hi <= r else box["bounds"][0]
                B = h-A*r
            else:
                A = box["bounds"][3] if hi <= v else box["dlo"]
                B = box["eta"]-A*v
            poly = polynomial(label, A, B, lo, hi, chord, call)
            exact = integrate(poly, lo, hi)
            demand(0 <= exact < Q(upper, DEN),
                   f"failed strict row integral: [{lo},{hi}], {label}")
            density_counts[label] += 1
            zero_integrals += exact == 0
            row_margins.append(Q(upper, DEN)-exact)
            total += upper
        return total

    for ident in range(301):
        box = boxes[ident]
        al, au, dl, du = box["bounds"]
        demand(al < au and dl < du, f"box {ident}: degenerate or reversed rectangle")
        eta = max(H, h+(v-r)*al-eps)
        dlo = max(dl, eta/v)
        demand(box["eta"] == eta and box["dlo"] == dlo,
               f"box {ident}: wrong linked lower height or slope")
        demand(dlo <= du, f"box {ident}: empty feasible domain cannot be discarded")
        shape = clipped_dimension(box["bounds"])
        demand(shape != "empty", f"box {ident}: no feasible point")
        clips[shape] += 1
        if box["floor"]:
            floor_count += 1
            demand(not box["in"] and box["in_sum"] is None, "stray capped incoming data")
            incoming = 31*DEN//32
        else:
            raw_incoming = rows_checked(box["in"], ZERO, v, box)
            demand(raw_incoming == box["in_sum"], f"box {ident}: incoming sum mismatch")
            incoming = min(raw_incoming, 31*DEN//32)
        outgoing = rows_checked(box["out"], v, T, box)
        demand(outgoing == box["out_sum"], f"box {ident}: outgoing sum mismatch")
        M = eta+dlo*(T-v)-eps
        reserve = max(ZERO, 1-M)**7
        zero_reserves += reserve == 0
        demand(reserve <= Q(box["reserve"], DEN), f"box {ident}: reserve underestimated")
        total = incoming+outgoing+box["reserve"]
        demand(total == box["total"] and total < LIMIT, f"box {ident}: failed total")
        box_margins.append(Q(LIMIT-total, DEN))
        finite_count += len(box["in"])+len(box["out"])
    demand(finite_count == 4258 and floor_count == 7, "finite record or incoming-cap count")
    largest = max(box["total"] for box in boxes.values())
    demand(data["largest"] == largest == MAX_TOTAL, "largest finite total mismatch")

    root = (h/r, Q(3), H/v, Q(20))
    demand(data["root"] == root, "root rectangle mismatch")
    leaves = []
    nodes = 0
    depth_max = 0

    def walk(node, bounds, depth=0):
        nonlocal nodes, depth_max
        demand(depth <= 64, "invalid subdivision depth")
        depth_max = max(depth_max, depth)
        if type(node) is str:
            demand(re.fullmatch(r"0|[1-9][0-9]*", node), "unknown or empty leaf")
            ident = int(node)
            demand(ident in boxes and ident not in leaves, "missing or duplicate tree leaf")
            demand(boxes[ident]["bounds"] == bounds, f"leaf {ident}: rectangle mismatch")
            leaves.append(ident)
            return
        demand(type(node) is list and len(node) == 3 and type(node[0]) is int
               and node[0] in (0, 1), "invalid binary subdivision node")
        nodes += 1
        axis = 2*node[0]
        midpoint = (bounds[axis]+bounds[axis+1])/2
        low, high = list(bounds), list(bounds)
        low[axis+1], high[axis] = midpoint, midpoint
        walk(node[1], tuple(low), depth+1)
        walk(node[2], tuple(high), depth+1)

    walk(data["tree"], root)
    demand(data["leaves"] == len(leaves) == 301 and nodes == 300
           and set(leaves) == set(boxes), "tree does not cover every box exactly once")
    area = sum(((b["bounds"][1]-b["bounds"][0])*(b["bounds"][3]-b["bounds"][2])
                for b in boxes.values()), ZERO)
    demand(area == data["area"] == (root[1]-root[0])*(root[3]-root[2])
           == Q(1858566687329, 33048000000), "root/leaf area mismatch")
    demand(data["coarse_max"] == 999894403
           and Q(largest, DEN) < Q(data["coarse_max"], 10**9) < Q(LIMIT, DEN),
           "coarse finite rounding failed")
    # Exercise the closed feasibility convention without deleting boundary sets.
    degenerate_examples = [
        ((h/r, h/r, H/v, H/v), "point"),
        ((h/r, h/r, H/v, H/v+1), "segment"),
        ((h/r, h/r, H/v-1, H/v-1), "empty")]
    for bounds, expected in degenerate_examples:
        demand(clipped_dimension(bounds) == expected, "closed clipping convention failed")

    demand(len(data["tail"]) == 11, "unbounded-a row count")
    a_rows = rows_checked(data["tail"], r, T)
    a_bound = r**7+Q(a_rows, DEN)+Q(1, 16)
    demand(a_bound < Q(data["a_total"], DEN) == Q(667049244496, DEN)
           < Q(largest, DEN), "unbounded-a total failed")
    demand(Q(659, 1000)**7 < Q(1, 16), "unbounded-a reserve failed")
    d_bound = Q(31, 32)+(v-H)**7/19
    demand(d_bound == data["d_exact"] < data["d_upper"] == Q(976627, 10**6)
           < Q(largest, DEN), "unbounded-d total failed")
    demand(H/v > cp and H+20*(T-v)-eps > 1, "unbounded-d branch or reserve failed")
    demand(a_bound < data["coarse_a"] == Q(667049245, 10**9)
           and d_bound < data["coarse_d"] == Q(976626676, 10**9), "tail rounding failed")
    demand(len(row_margins) == 4269, "complete integral record count")
    demand(Q(5, 38)-h == Q(3, 38000) > eps, "initial finite comparison slack")
    height_floor = T/v*H-eps
    demand(height_floor-Q(341, 1000) == Q(3, 17000)-Q(56, 27)*eps > 0,
           "linked-height finite comparison slack")
    reserve_margin = 1-Q(largest, DEN)
    demand(reserve_margin > Q(1, 10000), "strict source reserve failed")

    h2 = Q(341, 1000)
    kappa, z0 = T-h2, 25*(T-h2)/21
    demand(kappa == Q(819, 1000) and z0 == Q(39, 40), "small-slope breakpoints")
    demand(0 < kappa < z0 < T and 3*kappa-2*z0 > 0 and 2*kappa-T > 0,
           "incoming positive-part ranges")
    demand((z0-kappa)/z0 == Q(4, 25) and h2/T > Q(4, 25),
           "spanning-capacity switching threshold")
    I0 = (3*kappa**7-2*(2*kappa-z0)**7+(3*kappa-2*z0)**7-(2*kappa-T)**7)
    incoming_direct = (
        kappa**7
        + integrate(scale(add(scale(power(-1, 2*kappa, 6), 3),
                              scale(power(-2, 3*kappa, 6), -2)), 7), kappa, z0)
        + integrate(scale(power(-1, 2*kappa, 6), 7), z0, T))
    demand(I0 == incoming_direct == Q(315891331947925570507, 5*10**20)
           < Q(2, 3), "small-slope incoming exact integral")
    q0, z = 1-2*h2/T, 125/(52+73*(1-2*h2/T))
    demand(q0 == Q(239, 580) and T < z < Q(125, 52) < R,
           "small-slope outgoing crossing")
    I1 = q0**6*(Q(125, 52)*z**6-T**7)
    outgoing_direct = q0**6*(z**7-T**7)+(125-52*z)**7/(52*73**6)
    demand(125-52*z == 73*q0*z and I1 == outgoing_direct
           == Q(7860756587864860626248390838442045303558341047,
                59119162918550643946713187123828125000000000000)
           < Q(7, 50), "small-slope outgoing exact integral")
    demand(1-Q(2, 3)-Q(7, 50) == Q(29, 150) and h2 > T-1,
           "terminal reserve and large-slope comparison")
    return {
        "id": FAMILY, "boxes": 301, "finite_rows": finite_count,
        "unbounded_a_rows": 11, "total_integral_records": len(row_margins),
        "incoming_caps": floor_count, "tree_leaves": len(leaves),
        "tree_internal_nodes": nodes, "tree_max_depth": depth_max,
        "root": list(map(str, root)), "root_area": str(area),
        "closed_feasible_dimensions": dict(clips),
        "degenerate_clipping_checks": ["empty", "segment", "point"],
        "density_counts": dict(sorted(density_counts.items())),
        "zero_integrals": zero_integrals, "zero_reserves": zero_reserves,
        "chord_pairs": len(chord), "call_pairs": len(call),
        "branch_differences": list(map(str, differences)),
        "max_total": str(Q(largest, DEN)),
        "worst_boxes": [i for i, b in boxes.items() if b["total"] == largest],
        "min_row_rounding_gap": str(min(row_margins)),
        "min_box_gap_to_9999_10000": str(min(box_margins)),
        "source_reserve": str(reserve_margin),
        "excess_over_1_10000": str(reserve_margin-Q(1, 10000)),
        "unbounded_a_bound": str(a_bound), "unbounded_delta_bound": str(d_bound),
        "linked_height_floor": str(height_floor),
        "small_slope_incoming": str(I0), "small_slope_outgoing": str(I1),
        "tested_interface": [
            "strict total parsing, exact expected indices and counts",
            "all row domains and positive parts, rational integrals and upper rounding",
            "incoming/outgoing partitions, same-height reserve, all box sums",
            "closed root bisection and complete unique leaves, feasible boundary strata",
            "both unbounded slope ranges and terminal source inequalities"]
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="directory containing paired-profile.txt")
    parser.add_argument("--report", type=Path, required=True, help="JSON result, including failures")
    args = parser.parse_args(argv)
    report = {"schema_version": 1, "checker": FAMILY, "arithmetic": "exact rational"}
    try:
        demand(args.data.is_dir(), "data directory does not exist")
        demand({p.name for p in args.data.iterdir()} == {DATA_NAME},
               "data directory has missing or unexpected files")
        raw = (args.data/DATA_NAME).read_bytes()
        family = validate(parse(raw))
        report.update(status="pass", families=[family],
                      data_sha256=hashlib.sha256(raw).hexdigest(),
                      scope="Certificate arithmetic and coverage; geometric applicability is proved in the article.")
        code = 0
    except (ValueError, TypeError, KeyError, IndexError, OSError, SyntaxError,
            ZeroDivisionError, RecursionError) as error:
        report.update(status="fail", error=f"{type(error).__name__}: {error}")
        code = 1
    try:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    except OSError as error:
        print("cannot write report: " + str(error), file=sys.stderr)
        return 1
    print(json.dumps({"status": report["status"], "family": FAMILY,
                      "report": str(args.report), **({"error": report["error"]} if code else {})}))
    return code


if __name__ == "__main__":
    sys.exit(main())
