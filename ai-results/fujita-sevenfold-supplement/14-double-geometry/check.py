"""Portable exact checker for the Section 14 rational certificates.

Python 3.10+, standard library only. No generator or private state is loaded.
Usage: python3 check.py --data data --report report.json
On the Danus host, stage and execute only through danus compute DEFAULT.
Acceptance checks finite arithmetic conditional on the article's geometric
reductions; it is not a formal or whole-paper proof verification.
"""

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import json
import re
import argparse
import sys


R, U, Z, T = map(F, ("9/10", "26/25", "27/25", "29/25"))
H0, H2, EPS = map(F, ("131/1000", "1609/5000", "1/100000"))
QLO, QHI = map(F, ("16020581569/90194313216", "8010290789/45097156608"))
C4 = F(138259, 900000)
SCALE, TARGET = 10**10, 9990000000
ANCHORS = (R, U, Z)
HMIN = {"Q": F(37, 200), "M": F(41, 200), "D": F(189, 1000)}
CUT_C = tuple(map(F, ("0", "1/7", "4/25", "1/6", "2/11", "1/4", "1/3", "1/2", "1"))) + (QLO, QHI)
ZERO = (F(0),) * 7
T6 = ZERO[:6] + (F(1),)
COUNTS = Counter()
CELL_RECORDS = []
EXPECTED = {"Q": (22, 172, 1445, 9986634132),
            "M": (38, 188, 1827, 9989635195),
            "D": (0, 150, 1442, 9796511857)}
TAIL_AFFINE_COUNTS = (1,1,1,1,1,1,1,3,4,1,4,3,1,1,4,2,1,1,1,1,3,3,4,1,1,1,1,
                      4,4,3,3,7,4,5,5,3,3,3,5,5,4,5,3,3,4,4,3,3,4,4,7,6,6,4,
                      7,6,6,7,6,5,6,5,3,2,1,7,7,7,6,6,4)


def need(condition, message):
    if not condition:
        raise AssertionError(message)


def rational(token):
    need(re.fullmatch(r"(?:0|-?[1-9][0-9]*)(?:/[1-9][0-9]*)?", token) is not None,
         ("invalid rational token", token))
    value = F(token)
    need(str(value) == token, ("noncanonical rational", token))
    return value


def natural(token):
    need(re.fullmatch(r"0|[1-9][0-9]*", token) is not None, ("invalid integer", token))
    return int(token)


def add(*polys):
    return tuple(sum(p[i] for p in polys) for i in range(7))


def scale(poly, q):
    return tuple(q * c for c in poly)


def linear_power(intercept, slope, exponent):
    return tuple(F(comb(exponent, i)) * intercept ** (exponent-i) * slope**i
                 if i <= exponent else F(0) for i in range(7))


def mul(p, q):
    return tuple(sum(p[j] * q[i-j] for j in range(i+1)) for i in range(7))


def integrate(p, left, right):
    return sum(coef * (right**(i+1) - left**(i+1)) / (i+1)
               for i, coef in enumerate(p))


def baseline(t):
    return F(1) if t < R else F(61, 100) if t < F(99, 100) else F(11, 25) if t < U else F(17, 50)


def height_pair(family, box):
    h1 = max(HMIN[family], H0 + (U-R)*box[0] - EPS)
    h2 = max(H2, H0 + (Z-R)*box[0] - EPS,
             h1 + (Z-U)*box[2] - EPS)
    return h1, h2


def profile_data(box, heights, profile, t):
    a = box[2*profile+1] if t < ANCHORS[profile] else box[2*profile]
    h = (H0, *heights)[profile]
    d = a*ANCHORS[profile] - h
    return a, d


def poly_for(family, box, heights, label, t):
    if label.startswith("B"):
        need(label in ("B1", "B61/100", "B11/25", "B17/50"),
             ("unknown baseline label", label))
        rate = F(label[1:])
        need(rate >= baseline(t), ("baseline invalid", label, t))
        return scale(T6, 7*rate)
    match = re.fullmatch(r"([012]):(N[35](?:lo|hi)|C|E1|Qbal|Mspan|Mref|H3|H3lip|Dlo|Dhi)", label)
    need(match is not None, ("unknown polynomial", label))
    profile, name = int(match[1]), match[2]
    kind = "N" if profile == 0 else family if profile == 1 else "Q"
    a, d = profile_data(box, heights, profile, t)
    c = max(F(0), a-d/t)

    def pj(j):
        return scale(linear_power(j*d, 1-j*a, 6), 7) if 1-j*c > 0 else ZERO

    if name.startswith("N"):
        q = int(name[1])
        need((kind == "N" and q == 5) or (kind == "M" and q == 3), ("wrong normal codimension", kind, name))
        if name.endswith("lo"):
            need(0 < c <= F(1, 7), ("Nlo domain", c))
            p = linear_power(-d, a, q)
            p = (F(0),)*(6-q) + p[:q+1]
            return add(scale(T6, 7), scale(p, -7*F(7*q, q+1)**q))
        need(F(1, 7) <= c < F(1, 2), ("Nhi domain", c))
        p = linear_power(-7*d, 6*q-1+7*a, q)
        p = (F(0),)*(6-q) + p[:q+1]
        return add(scale(T6, 7), scale(p, -F(7, (6*(q+1))**q)))
    if name == "C":
        need(kind in ("N", "M") and F(1, 7) <= c < F(1, 2), ("C domain", kind, c))
        return scale(pj(1), F(7, 6)**6)
    if name == "E1":
        cutoff = QHI if kind == "Q" else F(1, 2)
        need(kind in ("N", "M", "Q") and c >= cutoff, ("E1 domain", kind, c))
        return pj(1)
    if name == "Qbal":
        need(kind == "Q" and 0 < c <= QLO, ("Qbal domain", kind, c))
        return scale(mul(linear_power(2*d, 1-2*a, 5), linear_power(-10*d, 1+10*a, 1)), 7)
    if name == "Mspan":
        need(kind == "M" and c > 0, ("Mspan domain", kind, c))
        return add(pj(1), scale(pj(2), 3), scale(pj(3), -5), scale(pj(4), 2))
    if name == "Mref":
        need(kind == "M" and c >= F(2, 11), ("Mref domain", kind, c))
        return add(pj(1), pj(2), scale(pj(3), -1))
    if name == "H3":
        need(kind == "M" and c >= F(1, 6), ("H3 domain", kind, c))
        return scale(T6, 7*F(256, 729))
    if name == "H3lip":
        need(kind == "M" and 0 < c < F(1, 6), ("H3lip domain", kind, c))
        return (F(0),)*5 + (126*d, 7*(F(256, 729)+3-18*a))
    if name == "Dlo":
        need(kind == "D" and c > 0, ("Dlo domain", kind, c))
        return add(scale(pj(2), 6), scale(pj(3), -8), scale(pj(4), 3))
    need(name == "Dhi" and kind == "D" and c >= F(1, 6), ("Dhi domain", kind, c))
    return scale(T6, 7*C4)


def split_cell(box, heights, label, left, right):
    # Exact affine roots make a midpoint domain check a whole-open-interval check.
    # No sampling approximation is used: c(t)=a-d/t has constant derivative sign.
    cuts = {left, right}
    cuts.update(t for t in (R, F(99, 100), U, Z) if left < t < right)
    if not label.startswith("B"):
        profile = int(label[0])
        broad = sorted(cuts)
        for lo, hi in zip(broad, broad[1:]):
            a, d = profile_data(box, heights, profile, (lo+hi)/2)
            for value in CUT_C:
                if a != value:
                    root = d/(a-value)
                    if lo < root < hi:
                        cuts.add(root)
    return sorted(cuts)


def audit_cell(family, tag, box, heights, cell):
    line, left, right, label, upper = cell
    need(F(0) <= left < right <= T, ("bad time interval", line))
    cuts = split_cell(box, heights, label, left, right)
    value = F(0)
    polys = set()
    for lo, hi in zip(cuts, cuts[1:]):
        polynomial = poly_for(family, box, heights, label, (lo+hi)/2)
        value += integrate(polynomial, lo, hi)
        polys.add(polynomial)
        COUNTS["exact_domain_subintervals"] += 1
    need(value >= 0, ("negative selected integral", line, value))
    rounded = (SCALE*value).__floor__() + 1
    need(upper == rounded, ("wrong strict numerator", family, tag, line, upper, rounded, value))
    need(value < F(upper, SCALE), ("not strict", line))
    # A merged cell must really have one polynomial (including the zero polynomial).
    need(len(polys) == 1, ("unprinted polynomial change inside merged cell", family, tag, line))
    COUNTS["printed_cells"] += 1
    COUNTS["label:"+label] += 1
    CELL_RECORDS.append({"line": line, "family": family, "box": tag,
                         "left": str(left), "right": str(right), "label": label,
                         "integral": str(value), "strict_upper_numerator": upper,
                         "verified": True})
    return value


def parse_certificate(text):
    """Consume the complete canonical stream with a strict ordered grammar."""
    lines = [(i, line.strip()) for i, line in enumerate(text.splitlines(), 1)
             if line.strip()]
    cursor = 0

    def peek():
        return lines[cursor][1] if cursor < len(lines) else ""

    def take(pattern):
        nonlocal cursor
        need(cursor < len(lines), ("unexpected end of certificate", pattern))
        number, line = lines[cursor]
        match = re.fullmatch(pattern, line)
        need(match is not None, ("unexpected record", number, line, pattern))
        cursor += 1
        return number, match

    def literal(expected):
        take(re.escape(expected))

    def fields(count):
        nonlocal cursor
        need(cursor < len(lines), "unexpected end in record")
        number, line = lines[cursor]
        values = [s.strip() for s in line.split("|")]
        need(len(values) == count, ("wrong field count", number, count))
        cursor += 1
        return number, values

    def tag(value):
        need(re.fullmatch(r"[1-6]\.[1-5]\.[1-5][01]*", value) is not None,
             ("invalid box tag", value))
        return value

    def cell_list():
        records = []
        while re.match(r"(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)? \|", peek()):
            number, values = fields(4)
            records.append((number, rational(values[0]), rational(values[1]),
                            values[2], natural(values[3])))
        need(records, "empty integration list")
        return records

    families = {}
    for family in ("Q", "M", "D"):
        literal(f"FAMILY {family}.hmin={HMIN[family]}.")
        data = {"splits": [], "boxes": {}, "cells": {}, "tails": {}}
        _, match = take(r"Grid0=(.*);grid1=(.*);grid2=(.*)\.")
        data["grids"] = tuple(tuple(rational(v) for v in match[i].split(","))
                              for i in (1, 2, 3))
        literal("Splits: parent tag | coordinate | midpoint")
        parents = set()
        while not peek().startswith("Rectangles:"):
            number, values = fields(3)
            parent = tag(values[0])
            need(parent not in parents, ("duplicate split", family, parent))
            parents.add(parent)
            data["splits"].append((number, parent, natural(values[1]), rational(values[2])))
        literal("Rectangles: tag | l0 | u0 | l1 | u1 | l2 | u2 | h1 | h2 | total upper numerator")
        while not peek().startswith("Rectangle "):
            number, values = fields(10)
            key = tag(values[0])
            need(key not in data["boxes"], ("duplicate box", family, key))
            data["boxes"][key] = (number, tuple(map(rational, values[1:7])),
                                  tuple(map(rational, values[7:9])), natural(values[9]))
        while peek().startswith("Rectangle "):
            _, match = take(r"Rectangle ([\d.]+)\. l \| r \| polynomial \| upper numerator")
            key = tag(match[1])
            need(key not in data["cells"], ("duplicate cell list", family, key))
            data["cells"][key] = cell_list()
        literal("Unbounded a0>=16: l | r | polynomial | upper numerator")
        data["cells"]["unbounded0"] = cell_list()
        _, match = take(r"a0>=16 total upper=([0-9]+)\.")
        data["tails"][0] = natural(match[1])
        for index, anchor in ((1, 26), (2, 27)):
            _, match = take(rf"a{index}>=8 post{anchor} loss<=([0-9/]+); volume left>([0-9/]+)\.")
            data["tails"][index] = (rational(match[1]), rational(match[2]))
        nsplits, nboxes, ncells, _ = EXPECTED[family]
        need(len(data["splits"]) == nsplits, ("split count", family))
        need(len(data["boxes"]) == nboxes, ("box count", family))
        need(sum(len(v) for k, v in data["cells"].items() if k != "unbounded0") == ncells,
             ("bounded integration record count", family))
        need(len(data["cells"]["unbounded0"]) == 5, ("unbounded record count", family))
        families[family] = data
    need(cursor == len(lines), ("unconsumed or extra records", lines[cursor:cursor+1]))
    return families, [number for number, _ in lines]


def audit_family(family, data):
    grids = data["grids"]
    expected_grids = ((H0/R, F(1, 2), F(1), F(2), F(4), F(8), F(16)),
                      (HMIN[family]/U, F(1, 2), F(1), F(2), F(4), F(8)),
                      (H2/Z, F(1, 2), F(1), F(2), F(4), F(8)))
    need(grids == expected_grids, ("grid coverage", family))
    leaves = {}
    for indices in product(*(range(len(g)-1) for g in grids)):
        key = ".".join(str(i+1) for i in indices)
        leaves[key] = tuple(v for g, i in zip(grids, indices) for v in (g[i], g[i+1]))
    for line, parent, coordinate, midpoint in data["splits"]:
        need(parent in leaves and coordinate in (0, 1, 2), ("invalid split parent", line))
        parent_box = leaves.pop(parent)
        need(parent_box[2*coordinate] < midpoint < parent_box[2*coordinate+1], ("invalid split point", line))
        for child in (0, 1):
            box = list(parent_box)
            box[2*coordinate+1-child] = midpoint
            key = parent + str(child)
            need(key not in leaves, ("split collision", line))
            leaves[key] = tuple(box)
    need(set(leaves) == set(data["boxes"]), ("leaf omissions or overlaps", family))
    need(set(data["cells"]) == set(leaves) | {"unbounded0"}, ("missing or extra cell list", family))
    exact_max, total_max, worst = F(-1), -1, None
    start_count = COUNTS["printed_cells"]
    for tag, (_, box, heights, printed_total) in data["boxes"].items():
        need(box == leaves[tag], ("wrong leaf bounds", family, tag))
        need(heights == height_pair(family, box), ("wrong linked heights", family, tag))
        cells = data["cells"][tag]
        cursor, exact, total = F(0), F(0), 0
        for cell in cells:
            need(cell[1] == cursor, ("gap or overlap in time", family, tag, cell[0]))
            exact += audit_cell(family, tag, box, heights, cell)
            cursor = cell[2]
            total += cell[4]
        need(cursor == T and total == printed_total < TARGET, ("box total or endpoint", family, tag, total))
        need(exact < F(total, SCALE), ("box rounding", family, tag))
        if total > total_max:
            total_max, worst = total, tag
        exact_max = max(exact_max, exact)
    bounded_cells = COUNTS["printed_cells"] - start_count
    # All other profiles are discarded in this branch; only old slope 16 is used.
    box = (F(16), F(16), F(8), F(8), F(8), F(8))
    heights = height_pair(family, box)
    cursor, unbounded0, total0 = F(0), F(0), 0
    for cell in data["cells"]["unbounded0"]:
        need(cell[1] == cursor, ("tail0 gap", family, cell[0]))
        need(cell[3] == "B1" or cell[3].startswith("0:"), ("tail0 borrowed profile", cell[0]))
        unbounded0 += audit_cell(family, "unbounded0", box, heights, cell)
        cursor, total0 = cell[2], total0+cell[4]
    need(cursor == T and total0 == data["tails"][0] == 5020326735 < TARGET, ("tail0 total", family))
    h = HMIN[family]
    if family == "Q":
        need(h/U > QHI, "Q high slope domain")
        loss1 = (U-h)**7/7
    elif family == "M":
        need(h/U > F(2, 11), "M high slope domain")
        loss1 = (U-h)**7/7 + (U-2*h)**7/15 - (U-3*h)**7/23
    else:
        need(U-2*h > 0, "D high slope numerator")
        loss1 = 6*(U-2*h)**7/15
    need(data["tails"][1] == (loss1, F(3, 40)-loss1), ("tail1 printed fractions", family))
    need(F(3, 40)-loss1 > F(1, 1000), ("tail1 margin", family))
    need(H2/Z > QHI, "final high slope domain")
    loss2 = (Z-H2)**7/7
    need(data["tails"][2] == (loss2, F(1, 32)-loss2), ("tail2 printed fractions", family))
    need(F(1, 32)-loss2 > F(1, 1000), ("tail2 margin", family))
    expected = {"Q": (172, 9986634132), "M": (188, 9989635195), "D": (150, 9796511857)}
    need((len(leaves), total_max) == expected[family], ("claimed family maxima", family))
    return {"initial_boxes": 150, "splits": len(data["splits"]), "final_boxes": len(leaves),
            "bounded_cells": bounded_cells, "maximum_upper_numerator": total_max,
            "worst_box": worst, "gap_below_999_over_1000": str(F(TARGET-total_max, SCALE)),
            "exact_maximum_selected_integral": str(exact_max),
            "unbounded0_cells": len(data["cells"]["unbounded0"]),
            "unbounded0_upper_numerator": total0,
            "unbounded1_loss": str(loss1), "unbounded1_volume_left": str(F(3, 40)-loss1),
            "unbounded2_loss": str(loss2), "unbounded2_volume_left": str(F(1, 32)-loss2)}


def audit_budget():
    a, p, u = F(623, 500), F(39, 40), Z
    payment = (F(13, 2)-F(7, 2))/u + F(5, 4)*F(7, 2)+F(92, 125)
    lower_b = ((8-p)*u-F(13, 2))/(a*u-1)
    inv_s = u/(F(13, 2)-lower_b)
    upper_sigma = (a*u-1)*F(10, 3)+7-8*u+u*p
    need(payment == F(70999, 9000), "double payment arithmetic")
    need(F(79, 10)-payment == F(101, 9000) > F(1, 1000), "double error margin")
    need(lower_b == F(27175, 8642), "cubic b arithmetic")
    need(inv_s == F(4321, 13425) and inv_s > H2+EPS, "cubic height gap")
    need(upper_sigma == F(8479, 15000) < F(2, 3), "triple-plane exclusion arithmetic")
    gaps = {"initial": F(5, 38)-H0, "Q": F(13, 70)-HMIN["Q"],
            "M": F(104, 505)-HMIN["M"], "D": F(52, 275)-HMIN["D"],
            "final": inv_s-H2}
    need(all(g > EPS for g in gaps.values()), "height gap too small")
    # The accepted unique crossing plus these exact signs justifies the omitted gap.
    def branch_difference(c):
        return (1-2*c)**5*(1+10*c)-(1-c)**6
    need(F(0) < QLO < QHI < F(1, 2), "crossing bracket order")
    need(branch_difference(QLO) > 0 and branch_difference(QHI) < 0, "crossing bracket signs")
    reserve = 1-F(39,100)*R**7-F(17,100)*F(99,100)**7-F(11,25)*U**7
    need(reserve > F(3,40), "full earlier original-source reserve")
    return {"double_base_plus_tail": str(payment), "double_margin_below_79_over_10": str(F(79, 10)-payment),
            "cubic_b_lower": str(lower_b), "cubic_inverse_s_lower": str(inv_s),
            "cubic_sigma_upper": str(upper_sigma), "preplace_height_gaps": {k: str(v) for k, v in gaps.items()},
            "Q_crossing_bracket_width": str(QHI-QLO),
            "full_earlier_reserve": str(reserve), "reserve_gap_above_3_over_40": str(reserve-F(3,40))}


def audit_tail_rows(path):
    rows = {}
    raw_lines = path.read_text(encoding="utf-8").splitlines()
    need(len(raw_lines) == 73, "tail table line count")
    need(raw_lines[0] == "| id | d | E | m | mu | B | r | upper lines (slope,intercept) |",
         "tail table header")
    need(raw_lines[1] == "|---|---|---|---|---|---|---|---|", "tail table separator")
    for number, raw in enumerate(raw_lines[2:], 1):
        need(raw.startswith("| ") and raw.endswith(" |"), ("tail row delimiters", number))
        fields = [x.strip() for x in raw[1:-1].split("|")]
        need(len(fields) == 8, ("tail fields", number))
        index, d, e, m = map(natural, fields[:4])
        need(index == number, ("missing/duplicate/out-of-order tail row", number, index))
        mu, cap, radius = map(rational, fields[4:7])
        need(re.fullmatch(r"\([0-9/]+,[0-9/]+\)(?:; \([0-9/]+,[0-9/]+\))*", fields[7]),
             ("invalid affine list", index))
        lines = tuple((rational(a), rational(p)) for a, p in
                      re.findall(r"\(([0-9/]+),([0-9/]+)\)", fields[7]))
        need(len(set(lines)) == len(lines), ("duplicate affine line", index))
        need(len(lines) == TAIL_AFFINE_COUNTS[index-1], ("missing/extra affine witness", index))
        rows[index] = (d, e, m, mu, cap, radius, lines)
    need((F(5, 4), F(92, 125)) in rows[51][-1], "missing double tail line")
    need((F(623, 500), F(39, 40)) in rows[55][-1], "missing cubic tail line")
    # Regenerate all descendant Hilbert states by filling the lowest degrees,
    # independently of both the table's multiplicity ranges and its mu values.
    expected_states = {(1, 1, 1): F(0)}
    for d in range(2, 7):
        for e in range(d, 8):
            c = e-d
            if c == 0:
                expected_states[(d, e, 1)] = F(0)
                continue
            m = c+1
            while True:
                remaining, degree, moment = m, 0, 0
                while remaining:
                    take = min(remaining, comb(c+degree-1, degree))
                    moment += degree*take
                    remaining -= take
                    degree += 1
                if F(d)-F(2*moment, m) <= 0:
                    break
                expected_states[(d, e, m)] = F(moment)
                m += 1
    actual_states = {(d, e, m): mu for d, e, m, mu, cap, radius, lines in rows.values()}
    need(len(actual_states) == 71 and actual_states == expected_states,
         "full state list incomplete or wrong first moment")
    # The article uses the full table also at its earlier observations.
    # Extend the historical two-root audit to every state and affine line.
    audited = list(rows)
    checks = 0
    affine_checks = 0
    nonpositive_caps = 0
    for index in audited:
        d, e, m, mu, cap, radius, lines = rows[index]
        need(cap == d-2*mu/m and cap > 0 and radius**d >= m, ("tail row data", index))
        need((radius*1000).denominator == 1 and (radius-F(1, 1000))**d < m,
             ("radius is not the minimal upward thousandth", index))
        for slope, intercept in lines:
            rhs = max(F(0), radius-slope)*cap
            surplus = F(0)
            for child, (dd, ee, mm, muu, cc, rr, ll) in rows.items():
                if dd >= d or ee > e:
                    continue
                need(child in audited, ("missed dependency", index, child))
                cutoff = min(cap, cc, cap-F(mm, m)) if dd == d-1 else min(cap, cc)
                if cutoff <= 0:
                    nonpositive_caps += 1
                    continue
                points = {F(0), cutoff}
                for a0, p0 in ll:
                    for a1, p1 in ll:
                        if a0 != a1:
                            y = (p1-p0)/(a0-a1)
                            if 0 < y < cutoff:
                                points.add(y)
                for y in points:
                    candidate = min(aa*y+pp for aa, pp in ll)-max(slope, radius)*y
                    surplus = max(surplus, candidate)
                    checks += 1
            need(intercept >= rhs+surplus, ("tail certificate inequality", index, slope, intercept, rhs+surplus))
            affine_checks += 1
    need({19, 20, 24, 25}.issubset(audited), "missing multiplicity 8/9 cases")
    return {"data_sha256": sha256(path.read_bytes()).hexdigest(),
            "roots": [51, 55], "rows_checked": audited, "exact_vertex_checks": checks,
            "affine_lines_checked": affine_checks,
            "nonpositive_successor_caps": nonpositive_caps,
            "independently_regenerated_states": len(expected_states),
            "multiplicity_8_9_threefold_rows": [19, 20, 24, 25]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    COUNTS.clear()
    CELL_RECORDS.clear()
    result = {"schema": "fujita-section14-check-v1", "status": "FAIL", "families": {}}
    try:
        expected_files = {"three-profile.txt", "tail-table.txt"}
        need(args.data.is_dir(), "missing data directory")
        need({p.name for p in args.data.iterdir()} == expected_files, "missing or extra data files")
        for name in expected_files:
            need((args.data/name).is_file(), ("nonfile payload", name))
        result["data_sha256"] = {name: sha256((args.data/name).read_bytes()).hexdigest()
                                 for name in sorted(expected_files)}
        data, consumed = parse_certificate((args.data/"three-profile.txt").read_text(encoding="utf-8"))
        result["certificate_nonblank_lines_consumed"] = len(consumed)
        result["budget"] = audit_budget()
        result["tail_table"] = audit_tail_rows(args.data/"tail-table.txt")
        for family in ("Q", "M", "D"):
            result["families"][family] = audit_family(family, data[family])
            need(result["families"][family]["bounded_cells"] == EXPECTED[family][2],
                 ("validated integration count", family))
            print(f"Checked {family}: {EXPECTED[family][1]} boxes, "
                  f"{EXPECTED[family][2]} bounded integration records", flush=True)
        need(COUNTS["printed_cells"] == 4729, "total integration record count")
        result["counts"] = dict(COUNTS)
        result["coverage"] = {"parameter_boxes": 510, "binary_splits": 60,
                              "bounded_integration_records": 4714,
                              "unbounded_first_slope_records": 15,
                              "total_integration_records": 4729,
                              "closed_root_boxes_per_family": 150,
                              "unbounded_cases_per_family": 3}
        result["proof_boundary"] = ("Exact certificate inequalities and complete finite coverage; "
                                     "the article proves geometric applicability, profile bounds, "
                                     "cutting errors, and counts in every late original degree.")
        result["status"] = "PASS"
    except (AssertionError, ValueError, OSError, ZeroDivisionError, OverflowError) as exc:
        result["failure"] = str(exc)
        print("FAIL:", exc, file=sys.stderr)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
