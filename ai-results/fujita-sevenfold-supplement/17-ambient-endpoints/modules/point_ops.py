"""Independent audit of the certificate embedded in fact 9814791a3b835060.

Run only through `danus compute computation/sub12_981_independent_audit.py`.
Reads the accepted text, not the author's generator or saved computed totals.
All mathematical decisions use exact Fraction arithmetic. Float values appear
only as human-readable companions to exact output. This is audit evidence.
"""

from .common import Q as F
from itertools import product
from math import comb
from pathlib import Path
from collections import Counter
import hashlib
import json
import re


R, U, Z, T, END = map(F, ["9/10", "26/25", "27/25", "29/25", "13/5"])
EPS = F(1, 100000)
QLO, QHI = F(16020581569, 90194313216), F(8010290789, 45097156608)
FAILURES = []
CONTEXT = ""


def check(condition, message):
    if not condition:
        raise ValueError(CONTEXT + ": " + message)


def add(*polys):
    result = [F(0)] * max(map(len, polys))
    for poly in polys:
        for j, coefficient in enumerate(poly):
            result[j] += coefficient
    return result


def scale(poly, coefficient):
    return [coefficient * x for x in poly]


def multiply(first, second):
    result = [F(0)] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i + j] += a * b
    return result


def power_line(intercept, slope, degree):
    return [F(comb(degree, j)) * intercept ** (degree - j) * slope ** j
            for j in range(degree + 1)]


def positive_power(intercept, slope, degree, left, right):
    endpoint_values = (intercept + slope * left, intercept + slope * right)
    if min(endpoint_values) >= 0:
        return power_line(intercept, slope, degree)
    if max(endpoint_values) <= 0:
        return [F(0)]
    check(False, "positive-part root lies strictly inside the written cell")
    raise ValueError(CONTEXT)


def integral(poly, left, right):
    return sum((c * (right ** (j + 1) - left ** (j + 1)) / (j + 1)
                for j, c in enumerate(poly)), F(0))


def polynom(label, left, right, box, heights, anchors, kinds):
    if label == "Bzero":
        check(left >= F(125, 52), "Bzero before its valid cutoff")
        return [F(0)]
    if label == "Bhigh":
        check(left >= T and right <= F(125, 52), "Bhigh outside its branch")
        return scale(power_line(F(125, 73), F(-52, 73), 6), 7)
    if label.startswith("B"):
        coefficient = F(label[1:])
        valid_intervals = {F(1): (F(0), R), F(61, 100): (R, F(99, 100)),
                           F(11, 25): (F(99, 100), U), F(17, 50): (U, T)}
        lo, hi = valid_intervals[coefficient]
        check(left >= lo and right <= hi, "constant baseline outside its stated piece")
        return [F(0)] * 6 + [7 * coefficient]

    index_text, name = label.split(":")
    index = int(index_text)
    anchor, height = anchors[index], heights[index]
    lo, hi = box[2 * index:2 * index + 2]
    check(right <= anchor or left >= anchor or lo == hi,
          "selected profile crosses its anchor without a cell split")
    slope = hi if (left + right) / 2 < anchor else lo
    d = slope * anchor - height
    zleft, zright = slope * left - d, slope * right - d
    check(min(zleft, zright) >= 0, "profile polynomial used while height is negative")
    check(slope * (left + right) / 2 - d > 0, "zero profile requires baseline")

    def at_least(threshold):
        return min(zleft - threshold * left, zright - threshold * right) >= 0

    def at_most(threshold):
        return max(zleft - threshold * left, zright - threshold * right) <= 0

    def pj(j):
        return scale(positive_power(j * d, 1 - j * slope, 6, left, right), 7)

    if name == "N5lo":
        check(kinds[index] == "G" and at_most(F(1, 7)), "N5lo branch")
        return add([F(0)] * 6 + [F(7)],
                   scale([F(0)] + power_line(-d, slope, 5), -7 * F(35, 6) ** 5))
    if name == "N5hi":
        check(kinds[index] == "G" and at_least(F(1, 7)) and at_most(F(1, 2)),
              "N5hi branch")
        return add([F(0)] * 6 + [F(7)],
                   scale([F(0)] + power_line(-7 * d, 29 + 7 * slope, 5), -F(7, 36 ** 5)))
    if name == "C":
        check(kinds[index] == "G" and at_least(F(1, 7)) and at_most(F(1, 2)),
              "C branch")
        return scale(pj(1), F(7, 6) ** 6)
    if name == "Qbal":
        check(kinds[index] == "Q" and at_most(QLO), "Qbal exceeds lower crossing isolator")
        return scale(multiply(power_line(2 * d, 1 - 2 * slope, 5),
                              [-10 * d, 1 + 10 * slope]), 7)
    if name == "E1":
        check((kinds[index] == "Q" and at_least(QHI)) or
              (kinds[index] == "G" and at_least(F(1, 2))), "E1 branch")
        return pj(1)
    if name == "E2":
        check(kinds[index] == "T" and at_least(F(4, 25)), "E2 branch")
        return pj(2)
    if name == "Dlo":
        check(kinds[index] == "D", "Dlo profile kind")
        return add(scale(pj(2), 6), scale(pj(3), -8), scale(pj(4), 3))
    if name == "Dhi":
        check(kinds[index] == "D" and at_least(F(1, 6)), "Dhi branch")
        return [F(0)] * 6 + [7 * F(138259, 900000)]
    if name == "Tlow":
        check(kinds[index] == "T" and at_most(F(4, 25)), "Tlow branch")
        return add(scale(pj(2), 3), scale(pj(3), -2))
    if name == "Ttwo":
        check(kinds[index] == "T" and at_most(F(4, 25)), "Ttwo branch")
        return scale(pj(2), 2)
    raise ValueError("Unknown polynomial: " + label)


