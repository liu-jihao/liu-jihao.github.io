#!/usr/bin/env python3
"""Exact independent checker for the triple-fivefold portable certificate.

Python >= 3.9, standard library only. No producer, graph, TeX engine, network,
or historical proof file is imported. All predicates survive python -O.
Adapted arithmetic: sub13_55b9_check_only.py, sub8_cf557_check_only.py,
sub13_4eb_independent_check.py; provenance is documented in manifest.json.
"""
import argparse
import hashlib
import json
import re
import sys
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb
from pathlib import Path

G = F(3997, 5000)
THETA = F(803, 1000)
TARGET = F(19999, 20000)
IDS = list(range(1, 7)) + list(range(8, 21)) + list(range(28, 39))
SMALL_IDS = list(range(1, 7)) + list(range(8, 21))
FOUR_IDS = list(range(28, 38))
LINE_COUNTS = [1,1,1,1,1,1,3,4,1,4,3,1,1,4,2,1,1,1,1,4,4,3,3,7,4,5,5,3,3,3]
HEIGHTS = {
    3: (F(3,2), F(81,250), F(377,1000)),
    4: (F(283,200), F(379,1000), F(2,5)),
    5: (F(187,125), F(297,1000), F(2,5)),
    6: (F(783,500), F(269,1000), F(3,8)),
}
CEILINGS = {3:F(49884763,50000000), 4:F(24992633,25000000),
            5:F(99899271,100000000), 6:F(19998699,20000000)}
BIN_COUNTS = {3:4, 4:11, 5:11, 6:13}
PATTERNS = {3:[[(1,F(0)),(1,F(1)),(1,F(2))],
               [(1,F(0)),(2,F(3,2))],[(1,F(2)),(2,F(1,2))]],
            **{k:[[(2,F(0)),(1,F(k))]] for k in (4,5,6)}}
PREDICATES = 0


class Invalid(ValueError):
    pass


def require(condition, label):
    global PREDICATES
    PREDICATES += 1
    if not condition:
        raise Invalid(label)


def exact_keys(obj, keys, label):
    require(type(obj) is dict and set(obj) == set(keys.split()), label + ": exact fields")


def integer(x, label):
    require(type(x) is int, label + ": integer required")
    return x


def rational(x):
    require(type(x) is str and re.fullmatch(r"-?(?:0|[1-9]\d*)(?:/[1-9]\d*)?", x),
            "malformed rational: " + repr(x))
    return F(x)


def rationals(xs, length=None):
    require(type(xs) is list and (length is None or len(xs) == length), "rational tuple shape")
    return tuple(rational(x) for x in xs)


def same_state(record, s, label):
    v = record["state"]
    require(type(v) is list and len(v) == 3 and all(type(x) is int for x in v)
            and tuple(v) == (s["d"], s["e"], s["mu"]), label + ": state")


def indexed(records, key, expected, label):
    require(type(records) is list, label + ": list")
    out = {}
    for r in records:
        require(type(r) is dict and key in r, label + ": record")
        k = integer(r[key], label)
        require(k not in out, label + ": duplicate index " + str(k))
        out[k] = r
    require(set(out) == set(expected), label + ": missing or extra indices")
    require(list(out) == list(expected), label + ": noncanonical order")
    return out


def unique_rational_tuples(records, length, label):
    require(type(records) is list, label + ": list")
    out = [rationals(x, length) for x in records]
    require(len(out) == len(set(out)), label + ": duplicate tuple")
    return set(out)


def no_duplicate_keys(pairs):
    obj = {}
    for key, value in pairs:
        if key in obj:
            raise Invalid("duplicate JSON key: " + key)
        obj[key] = value
    return obj


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=no_duplicate_keys,
                      parse_float=lambda x: (_ for _ in ()).throw(Invalid("JSON float: " + x)),
                      parse_constant=lambda x: (_ for _ in ()).throw(Invalid("JSON constant: " + x)))


def moment(c, mu):
    if c == 0:
        return 0
    degree = 0
    while comb(c + degree, degree) < mu:
        degree += 1
    return degree * mu - comb(c + degree, degree - 1) if degree else 0


def tail(s, y):
    return min(a*y+c for a,c in s["lines"])


def crossings(s):
    return {(d-c)/(a-b) for (a,c),(b,d) in combinations(s["lines"],2) if a != b}


def endpoints(s, cap):
    require(cap >= 0, "negative interval passed to endpoints")
    return {F(0), cap} | {y for y in crossings(s) if 0 <= y <= cap}


def parent_cap(s):
    cap = min(F(3), s["cap"])
    return min(cap, F(9-s["mu"],3)) if s["d"] == 4 else cap


def feasible(p, hs):
    return all(a*p[0]+b*p[1] <= c for a,b,c in hs)


def meet(l1, l2):
    a,b,c = l1
    d,e,f = l2
    det = a*e-b*d
    if det == 0:
        return None
    return ((c*e-b*f)/det, (a*f-c*d)/det)


def vertices(hs):
    return {p for a,b in combinations(hs,2) if (p := meet(a,b)) is not None
            and feasible(p,hs)}


def constraints(s, gamma, cap, order=None):
    hs = [(-F(1),F(0),F(0)), (F(0),-F(1),F(0)), (F(1),F(0),F(3)),
          (-F(1),F(1),F(0)), (F(0),F(1),cap)]
    if order is not None:
        hs.append((3*gamma-order, order, gamma*(9-s["mu"])))
    return hs


def polygon_candidates(s, hs):
    pts = vertices(hs)
    for z in crossings(s):
        for line in hs:
            p = meet(line,(F(0),F(1),z))
            if p is not None and feasible(p,hs):
                pts.add(p)
    return pts


def affine_dimension(pts):
    if not pts:
        return -1
    if len(pts) == 1:
        return 0
    p, q, *rest = sorted(pts)
    return 2 if any((q[0]-p[0])*(r[1]-p[1]) != (q[1]-p[1])*(r[0]-p[0])
                    for r in rest) else 1


def price(s, gamma, b, z):
    return F(25,27)*(7-b)+(b-z)/gamma+tail(s,z)


def polygon_max(s, gamma, cap, order=None):
    hs = constraints(s,gamma,cap,order)
    pts = polygon_candidates(s,hs)
    require(bool(pts), "unexpected empty closed price domain")
    return max(price(s,gamma,*p) for p in pts), hs, pts


def check_tails(payload, view):
    exact_keys(payload, "schema rows", "tails")
    require(payload["schema"] == "triple-fivefold-tails-v1", "tail schema")
    raw = indexed(payload["rows"], "row", IDS, "tail rows")
    states = {}
    for n, count in zip(IDS,LINE_COUNTS):
        r = raw[n]
        exact_keys(r,"row state moment cap radius lines","tail")
        require(type(r["state"]) is list and len(r["state"]) == 3
                and all(type(x) is int for x in r["state"]), "tail state types")
        d,e,mu = r["state"]
        m = integer(r["moment"],"moment")
        cap, radius = rational(r["cap"]), rational(r["radius"])
        lines = unique_rational_tuples(r["lines"],2,"tail lines")
        require(len(lines) == count, "tail line coverage " + str(n))
        s = dict(d=d,e=e,mu=mu,moment=m,cap=cap,radius=radius,lines=sorted(lines))
        require(m == moment(e-d,mu), "Artinian moment " + str(n))
        require(cap == d-F(2*m,mu) > 0, "intrinsic cap " + str(n))
        require(radius > 0 and radius**d >= mu, "radius " + str(n))
        require(all(a >= 0 and c >= 0 for a,c in lines), "nonnegative tail coefficients")
        states[n] = s
    expected = {(1,1,1)}
    for d in range(2,5):
        for e in range(d,7):
            c = e-d
            if not c:
                expected.add((d,e,1))
                continue
            ell = d//2+1
            upper = F(2*ell*comb(c+ell-1,ell-1), 2*ell-d)
            for mu in range(c+1, int(upper)+1):
                if 2*moment(c,mu) < d*mu:
                    expected.add((d,e,mu))
    require({(s["d"],s["e"],s["mu"]) for s in states.values()} == expected,
            "complete independent state census")
    require([(s["d"],s["e"],s["mu"]) for s in states.values()] == sorted(expected),
            "historical row/state mapping")
    recurrences = 0
    for n,s in states.items():
        for a,c in s["lines"]:
            surplus = F(0)
            for w in states.values():
                if w["d"] >= s["d"] or w["e"] > s["e"]:
                    continue
                cap = min(s["cap"],w["cap"])
                if w["d"] == s["d"]-1:
                    cap = min(cap,s["cap"]-F(w["mu"],s["mu"]))
                if cap > 0:
                    surplus = max(surplus,*(tail(w,y)-max(a,s["radius"])*y
                                            for y in endpoints(w,cap)))
            require(c >= max(F(0),s["radius"]-a)*s["cap"]+surplus,
                    "full affine tail recurrence " + str((n,a,c)))
            recurrences += 1
    # The small printed table is a derived view of exactly these canonical rows.
    rendered = []
    for line in view.splitlines():
        if not line.startswith(r"\("):
            continue
        plain = re.sub(r"\\tfrac\{([^{}]+)\}\{([^{}]+)\}",r"\1/\2",line)
        plain = plain.replace(r"\(", "").replace(r"\)", "")
        plain = re.sub(r"\\\\\[3pt\]$", "", plain).strip()
        c = [x.strip() for x in plain.split("&")]
        require(len(c) == 7, "tail TeX row shape")
        pairs = re.findall(r"\(([^,()]+),([^,()]+)\)", c[6])
        require(re.sub(r"\([^()]+\)|[; ]","",c[6]) == "", "tail TeX pair consumption")
        rendered.append((tuple(map(int,c[:3])),int(c[3]),F(c[4]),F(c[5]),
                         sorted((F(a),F(b)) for a,b in pairs)))
    expect_view = [((s["d"],s["e"],s["mu"]),s["moment"],s["cap"],s["radius"],s["lines"])
                   for s in states.values()]
    require(rendered == expect_view, "rendered tail table equals canonical data")
    surplus = max([F(0)] + [tail(s,y)-F(217,200)*y for s in states.values()
                             if parent_cap(s)>0 for y in endpoints(s,parent_cap(s))])
    require(surplus == F(129,125), "degree-at-least-two full tail surplus")
    require(F(217,200)**5-F(3,2) == F(1170140857,320000000000), "parent radius identity")
    require(F(100,27)+3*F(217,200)+surplus == 8-F(251,27000), "degree-two price")
    return states, dict(rows=30, affine_recurrences=recurrences,
                        rendered_rows=30, degree_two_margin="251/27000")


def check_first(payload, states):
    exact_keys(payload, "schema gamma theta constrained_order initial_mass reached_mass summaries "
               "restricted_summaries roots polygons empty_domains quartic_low", "first source")
    require(payload["schema"] == "triple-fivefold-first-source-v1", "first schema")
    require(rational(payload["gamma"]) == G and rational(payload["theta"]) == THETA
            and rational(payload["constrained_order"]) == F(27,25), "source orders")
    initial, mass = 1-3*G**5, 1-G**5-2*THETA**5
    require(rational(payload["initial_mass"]) == initial, "initial source mass")
    require(rational(payload["reached_mass"]) == mass > F(1,200) > F(1,4)**5,
            "reached mass and full-degree reserve")
    expected = [("unrestricted",n) for n in IDS if n != 38] + [("restricted",n) for n in FOUR_IDS]
    require(type(payload["polygons"]) is list, "polygon list")
    polys, dimensions, witness_count = {}, {-1:0,0:0,1:0,2:0}, 0
    for p in payload["polygons"]:
        exact_keys(p,"mode row state cap vertices witnesses upper","polygon")
        n = integer(p["row"],"polygon row")
        require(type(p["mode"]) is str, "polygon mode type")
        key = (p["mode"],n)
        require(key in expected and key not in polys, "polygon index coverage")
        s = states[n]
        same_state(p,s,"polygon")
        cap = parent_cap(s)
        require(rational(p["cap"]) == cap, "polygon cap")
        order = F(27,25) if key[0] == "restricted" else None
        maximum, hs, candidates = polygon_max(s,G,cap,order)
        vs = unique_rational_tuples(p["vertices"],2,"polygon vertices")
        require(vs == vertices(hs), "all original polygon vertices " + str(key))
        ws = unique_rational_tuples(p["witnesses"],3,"polygon witnesses")
        require(len({x[:2] for x in ws}) == len(ws), "conflicting polygon witness")
        require({x[:2] for x in ws} == candidates, "complete subdivided vertices " + str(key))
        for b,z,value in ws:
            require(value == price(s,G,b,z), "altered polygon witness " + str((key,b,z)))
        require(rational(p["upper"]) == maximum, "polygon maximum " + str(key))
        polys[key] = (maximum, hs)
        dimensions[affine_dimension(vs)] += 1
        witness_count += len(ws)
    require(list(polys) == expected, "full polygon index domain")
    require(witness_count == 466, "all 466 source polygon witnesses")
    require(type(payload["empty_domains"]) is list and len(payload["empty_domains"]) == 2,
            "two zero-cap source rows")
    for e, mode in zip(payload["empty_domains"],["unrestricted","restricted"]):
        exact_keys(e,"mode row cap","empty source")
        require(e["mode"] == mode and integer(e["row"],"empty row") == 38
                and rational(e["cap"]) == parent_cap(states[38]) == 0, "zero-cap source")
        hs = constraints(states[38],G,F(0),F(27,25) if mode == "restricted" else None)
        vs = vertices(hs)
        require(all(y == 0 for _,y in vs), "empty positive-deficit closure")
        dimensions[affine_dimension(vs)] += 1
    for family, indices, mode in [("summaries",IDS[:-1],"unrestricted"),
                                 ("restricted_summaries",FOUR_IDS,"restricted")]:
        for n,r in indexed(payload[family],"row",indices,family).items():
            exact_keys(r,"row state cap upper point" if mode == "unrestricted"
                       else "row state upper point","summary")
            same_state(r,states[n],"summary")
            if mode == "unrestricted":
                require(rational(r["cap"]) == parent_cap(states[n]), "summary cap")
            p = rationals(r["point"],2)
            upper, hs = polys[(mode,n)]
            require(feasible(p,hs) and price(states[n],G,*p) == rational(r["upper"]) == upper,
                    "summary attainment")
    unpaid = {n for (mode,n),(upper,_) in polys.items()
              if mode == "unrestricted" and upper >= 8-F(1,1000)}
    require(unpaid == {15,32,33,34,35}, "all and only five scalar exceptions")
    require(max(value[0] for key,value in polys.items() if key[0] == "restricted")
            == F(430631,54000), "restricted fourfold maximum")
    low = unique_rational_tuples(payload["quartic_low"],2,"quartic interval")
    require(len({p[0] for p in low}) == len(low) and {p[0] for p in low} == endpoints(states[15],F(1)),
            "complete quartic low-deficit endpoints")
    require(all(value == price(states[15],G,F(3),y) for y,value in low), "quartic low values")
    require(max(v for _,v in low) == F(431653343,53959500) < 8-F(1,3000),
            "quartic low margin")
    require(F(100,27)+3/G == 7+F(49267,107919), "quartic entire-boundary interface")

    roots = indexed(payload["roots"],"mu",[3,4,5,6],"normal roots")
    root_results = []
    for mu,r in roots.items():
        exact_keys(r,"mu E radius cap descendants surplus upper gap","normal root")
        E = integer(r["E"],"normal moment")
        radius, cap = rational(r["radius"]), rational(r["cap"])
        require(E == {3:3,4:6,5:8,6:10}[mu], "normal cutoff identity")
        require(radius == states[mu+29]["radius"], "normal cutoff radius")
        require(cap == min(F(9-mu,3),4-F(2*E,mu)), "normal root cap")
        maxima, empty_count = [F(0)], 0
        descendants = indexed(r["descendants"],"row",SMALL_IDS,"normal descendants")
        for n,d in descendants.items():
            exact_keys(d,"row state cap empty points upper","normal descendant")
            s = states[n]
            same_state(d,s,"normal descendant")
            A = min(cap,s["cap"])
            if s["d"] == 3:
                A = min(A,F(4*mu-2*E-s["mu"],mu))
            require(rational(d["cap"]) == A and type(d["empty"]) is bool, "normal descendant cap")
            pts = unique_rational_tuples(d["points"],2,"normal endpoints")
            if A <= 0:
                require(d["empty"] and not pts and rational(d["upper"]) == 0, "nonpositive normal cap")
                empty_count += 1
                continue
            require(not d["empty"] and len({x[0] for x in pts}) == len(pts)
                    and {x[0] for x in pts} == endpoints(s,A), "all normal endpoints")
            require(all(value == tail(s,y)-radius*y for y,value in pts), "normal endpoint values")
            maximum = max(value for _,value in pts)
            require(rational(d["upper"]) == maximum, "normal descendant maximum")
            maxima.append(maximum)
        surplus = max(maxima)
        cost = F(100,27)+(3-cap)/G+radius*cap+surplus
        require(rational(r["surplus"]) == surplus and rational(r["upper"]) == cost
                and rational(r["gap"]) == 8-cost > F(1,100), "normal root cost and gap")
        # The second source has its own threshold. Increasing gamma only improves this bound.
        require(F(100,27)+(3-cap)/THETA+radius*cap+surplus < 8-F(1,100), "own theta normal cost")
        root_results.append(dict(mu=mu, descendants=19, nonpositive_caps=empty_count, gap=str(8-cost)))
    return dict(polygons=39, written_witnesses=witness_count, summaries=39, empty_positive_domains=2,
                closed_domain_dimensions={str(k):v for k,v in dimensions.items()}, normal_roots=root_results,
                initial_mass=str(initial), reached_mass=str(mass), mass_above_one_over_200=str(mass-F(1,200)),
                restricted_price_margin="1369/54000", quartic_low_margin="22657/53959500")


def certify_price(r, s, gamma, cap, order=None):
    same_state(r,s,"price row")
    require(rational(r["cap"]) == cap, "price cap")
    upper, hs, pts = polygon_max(s,gamma,cap,order)
    p = rationals(r["point"],2)
    require(feasible(p,hs) and price(s,gamma,*p) == rational(r["upper"]) == upper,
            "independent entire polygon maximum and printed maximizer")
    require(upper < F(79999,10000), "strict integral price margin")
    return upper, len(pts)


def check_integral(payload, states):
    exact_keys(payload,"schema first intervals comparisons losses terminal_level terminal_upper jump_whole_cap","integral")
    require(payload["schema"] == "triple-fivefold-integral-v1", "integral schema")
    for n,r in indexed(payload["first"],"row",SMALL_IDS,"integral first source").items():
        exact_keys(r,"row state cap upper point","integral first row")
        certify_price(r,states[n],G,F(1) if n == 15 else states[n]["cap"])
    intervals = indexed(payload["intervals"],"interval",range(1,18),"original intervals")
    losses = indexed(payload["losses"],"interval",range(1,18),"interval losses")
    require(type(payload["comparisons"]) is list, "comparison list")
    comparisons = {}
    for r in payload["comparisons"]:
        exact_keys(r,"interval row state cap upper point","integral comparison")
        key = (integer(r["interval"],"comparison interval"),integer(r["row"],"comparison row"))
        require(key not in comparisons, "duplicate integral comparison")
        comparisons[key] = r
    require(list(comparisons) == [(i,n) for i in range(1,18) for n in FOUR_IDS],
            "all 170 constrained comparisons in order")
    remaining, last, total_candidates, prices = 1-3*G**5, G, 0, []
    for i,r in intervals.items():
        exact_keys(r,"interval left right order upper point","interval")
        left,right,order = map(rational,[r["left"],r["right"],r["order"]])
        require(last == left < right <= F(23,25) and 1 < order < 3*left, "full ordered interval coverage")
        require((5000*left).denominator == 1 and (5000*right).denominator == 1, "same original degree integrality")
        values = []
        for n in FOUR_IDS:
            value,count = certify_price(comparisons[(i,n)],states[n],left,parent_cap(states[n]),order)
            values.append(value)
            prices.append(value)
            total_candidates += count
        require(rational(r["upper"]) == max(values), "cross-fourfold interval maximum")
        pt = rationals(r["point"],2)
        require(any(feasible(pt,constraints(states[n],left,parent_cap(states[n]),order))
                    and price(states[n],left,*pt) == max(values) for n in FOUR_IDS),
                "interval summary maximizing point")
        # Row 38 is excluded by z>0, but its closed face is still enumerated.
        vs = vertices(constraints(states[38],left,F(0),order))
        require(vs and all(z == 0 for _,z in vs), "multiplicity-nine zero-cap closure")
        loss = 3*(max(F(0),right-order/3)**5-max(F(0),left-order/3)**5)
        exact_keys(losses[i],"interval loss remaining","loss")
        require(loss >= 0 and rational(losses[i]["loss"]) == loss, "exact symbolic interval loss")
        remaining -= loss
        require(rational(losses[i]["remaining"]) == remaining > 0, "positive same-degree prefix reserve")
        last = right
    require(last == rational(payload["terminal_level"]) == F(23,25), "terminal level")
    require(remaining == F(51447479205198179,45000000000000000000), "final positive original mass")
    terminal = {n:polygon_max(s,last,parent_cap(s))[0] for n,s in states.items()}
    require(max(terminal.values()) == rational(payload["terminal_upper"]) == F(2481721,310500)
            < F(7997,1000), "all 30 terminal price domains")
    require(rational(payload["jump_whole_cap"]) == F(100,27)+2/G+1
            == F(777619,107919) < F(721,100), "retained historical jump cap")
    return dict(first_rows=19, intervals=17, constrained_comparisons=170,
                recomputed_candidates=total_candidates, zero_cap_closures=17, terminal_rows=30,
                final_mass=str(remaining), price_margin=str(8-max(prices)),
                terminal_max=str(max(terminal.values())), terminal_margin=str(8-max(terminal.values())))


def box_slice_vertices(caps, total):
    vs = set()
    for free in range(len(caps)):
        rest = [i for i in range(len(caps)) if i != free]
        for sides in product((0,1),repeat=len(rest)):
            v = [F(0)]*len(caps)
            for i,side in zip(rest,sides):
                v[i] = F(caps[i])*side
            v[free] = total-sum(v)
            if 0 <= v[free] <= caps[free]:
                vs.add(tuple(v))
    return vs


def positive_integral4(a,c,left,right):
    require(left <= right, "ordered integration endpoints")
    cuts = {left,right}
    if a and left < -c/a < right:
        cuts.add(-c/a)
    cuts = sorted(cuts)
    value = F(0)
    for x,y in zip(cuts,cuts[1:]):
        if a*(x+y)/2+c > 0:
            piece = 5*sum(F(comb(4,j),j+1)*a**j*c**(4-j)*(y**(j+1)-x**(j+1))
                          for j in range(5))
            require(piece >= 0, "nonnegative sign-split polynomial integral")
            value += piece
    return value


def component_integral(shift,h,left,right,test):
    onset = G-h/right
    return (positive_integral4(F(1),F(0),F(0),onset)
            + positive_integral4(1-shift*right,shift*(right*G-h),onset,G)
            + positive_integral4(1-shift*left,shift*(left*G-h),G,test))


def check_slopes(payload, states):
    exact_keys(payload,"schema heights families","slopes")
    require(payload["schema"] == "triple-fivefold-slopes-v1", "slope schema")
    hs = indexed(payload["heights"],"mu",[3,4,5,6],"height rows")
    for mu,r in hs.items():
        exact_keys(r,"mu D derivative lower buffer","height")
        a,c,h = HEIGHTS[mu]
        A,C,D = a-F(25,27), F(175,27)+c-8, a*G-1
        derivative, lower = 3*C-A*(mu-9), D/(3*A+C)
        require(A>0 and D>0 and 3*A+C>0 and derivative<0, "height derivative and denominators")
        require(3*F(21587,8642)+mu-9>0, "old-prime numerator on entire parent interval")
        require((rational(r["D"]),rational(r["derivative"]),rational(r["lower"]),rational(r["buffer"]))
                == (D,derivative,lower,lower-h), "exact height identities")
        require(lower-h>F(1,10000) and lower>F(3,8)>F(1,10) and h>=F(3,8),
                "finite order error, fixedness, and common prime mass")
        require((a,c) in states[mu+29]["lines"], "height line in canonical scalar data")
        if mu<6:
            require((F(3,4)-h)/(THETA-G)<=104, "entire first-three slope domains")
    theta_time = (8-F(50,9)-F(367,200))/(1-F(25,27)*THETA)
    require(theta_time == F(3291,1385) and 1/theta_time < F(421,1000), "own theta time")
    require((F(421,1000)-F(3,8))/(THETA-G)==F(115,9), "entire last slope domain")
    require(box_slice_vertices((2,2,2),3)==set(permutations((F(0),F(1),F(2)))),
            "all three-plane allocation vertices")
    require(box_slice_vertices((2,3),3)=={(F(0),F(3)),(F(2),F(1))},
            "all plane-quadric allocation vertices")
    # Audit both reached levels, retaining all exceptions and moment windows.
    own_rows = 0
    for gamma in (G,THETA):
        for n,s in states.items():
            cap = parent_cap(s)
            if cap <= 0:
                continue
            value = polygon_max(s,gamma,cap)[0]
            if n not in {15,32,33,34,35}:
                require(value < 8-F(1,1000), "own-source alternative at " + str((gamma,n)))
            if n == 15:
                require(polygon_max(s,gamma,F(1))[0] <= F(431653343,53959500), "own quartic low case")
            own_rows += 1
    family_results, integral_count, bin_count = [], 0, 0
    for mu,fam in indexed(payload["families"],"mu",[3,4,5,6],"slope families").items():
        exact_keys(fam,"mu height domain count ceiling bins","slope family")
        h, ceiling = rational(fam["height"]), rational(fam["ceiling"])
        lo,hi = rationals(fam["domain"],2)
        count = integer(fam["count"],"slope count")
        require(h == HEIGHTS[mu][2] and lo == h/G and hi == (104 if mu<6 else F(115,9)),
                "source-specific full slope domain")
        require(count == BIN_COUNTS[mu] and ceiling == CEILINGS[mu] < TARGET, "family count and strict ceiling")
        previous, values = lo, []
        for n,r in indexed(fam["bins"],"bin",range(1,count+1),"slope bins").items():
            exact_keys(r,"bin left right test local global_value upper","slope bin")
            left,right,test = map(rational,[r["left"],r["right"],r["test"]])
            require(previous == left < right <= hi and G <= test <= 1, "closed slope coverage/test domain")
            require(0 <= G-h/right <= G, "profile onset")
            written = rationals(r["local"],len(PATTERNS[mu]))
            calculated = []
            for pattern in PATTERNS[mu]:
                value = F(0)
                for weight,shift in pattern:
                    v = component_integral(shift,h,left,right,test)
                    if shift == 0:
                        require(v == test**5, "zero-shift component includes entire interval")
                    value += weight*v
                calculated.append(value)
            require(tuple(calculated) == written, "altered pattern integral " + str((mu,n)))
            reserve = max(F(0),1-h-left*(test-G))**5
            upper = max(calculated)+reserve
            require(reserve == rational(r["global_value"]) and upper == rational(r["upper"]),
                    "exact global reserve and all-pattern maximum")
            require(upper <= ceiling < TARGET < 1, "strict entire-bin bound")
            values.append(upper)
            integral_count += len(calculated)
            bin_count += 1
            previous = right
        require(previous == hi, "last endpoint included")
        family_results.append(dict(mu=mu,bins=count,maximum=str(max(values)),ceiling=str(ceiling),
                                   strict_target_gap=str(TARGET-ceiling),gap_below_one=str(1-ceiling)))
    require(bin_count == 39 and integral_count == 47, "full bin/pattern counts")
    return dict(bins=bin_count,pattern_integrals=integral_count,own_source_rows=own_rows,
                height_rows=4,families=family_results)


def validate(data):
    expected_files = {"tails.json","tail-view.tex","first-source.json","integral.json","slopes.json"}
    require(data.is_dir() and {x.name for x in data.iterdir()} == expected_files,
            "exact data file set")
    require(all((data/x).is_file() and not (data/x).is_symlink() for x in expected_files),
            "ordinary self-contained data files")
    tails, tail_report = check_tails(read_json(data/"tails.json"),(data/"tail-view.tex").read_text())
    first_report = check_first(read_json(data/"first-source.json"),tails)
    integral_report = check_integral(read_json(data/"integral.json"),tails)
    slope_report = check_slopes(read_json(data/"slopes.json"),tails)
    return dict(status="PASS",schema="triple-fivefold-check-report-v1",predicates=PREDICATES,
                families={"13-triple-fivefold/tails":tail_report,
                          "13-triple-fivefold/first-source":first_report,
                          "13-triple-fivefold/integral":integral_report,
                          "13-triple-fivefold/slopes":slope_report},
                input_sha256={x:hashlib.sha256((data/x).read_bytes()).hexdigest() for x in sorted(expected_files)},
                scope="Exact finite arithmetic and exhaustive data coverage, conditional on the article's geometric and asymptotic arguments.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data",type=Path,required=True)
    parser.add_argument("--report",type=Path,required=True)
    args = parser.parse_args()
    try:
        result = validate(args.data)
    except (Invalid, ValueError, TypeError, KeyError, IndexError, OSError, ZeroDivisionError) as exc:
        result = dict(status="FAIL",schema="triple-fivefold-check-report-v1",
                      predicates=PREDICATES,error=str(exc))
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
