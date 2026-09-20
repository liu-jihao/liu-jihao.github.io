#!/usr/bin/env python3
"""Validate the four delivered quartic-threefold certificate families.

Python >=3.9, standard library only. No generator or project-state imports.
The exact predicates adapt the earlier independent sub8/sub13 and sub11
checks; see manifest.json for the full provenance and article interfaces.
"""
import argparse
from ast import literal_eval
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import re
import sys

C = Q(804700, 107919)
A = Q(5000, 3997)
T = Q(799999, 100000)
K = T - C
END = Q(18, 25)
TARGET = Q(999, 1000)
SCALE = 10**8
CAP = {3: Q(2, 3), 4: Q(1, 2), 5: Q(1, 4)}
RAD = {3: Q(1733, 1000), 4: Q(2), 5: Q(2237, 1000)}
SPECS = {
    "two-quadrics": (25, 693, Q(49900103, 50000000)),
    "cubic-plane": (14, 436, Q(19979391, 20000000)),
    "simple-chains": (46, 1618, Q(24966091, 25000000)),
    "remaining": (40, 2481, Q(19975087, 20000000)),
}
CONTEXT = {}
COUNTS = Counter()


class InvalidCertificate(ValueError):
    pass


def require(test, message, **details):
    if not test:
        raise InvalidCertificate(json.dumps(dict(CONTEXT, error=message,
                                                 **details), default=str))


def integer(text):
    require(re.fullmatch(r"-?(?:0|[1-9][0-9]*)", text) is not None,
            "invalid integer", token=text)
    return int(text)


def rational(text):
    require(re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?", text)
            is not None, "invalid rational", token=text)
    return Q(text)


def ints(text):
    return tuple(integer(x) for x in text.split(","))


def rats(text):
    return tuple(rational(x) for x in text.split(","))


def fields(text, count, sep=";"):
    result = text.split(sep)
    require(len(result) == count, "record field count", record=text)
    return result


def cross(p, q):
    return p[0]*q[1] - p[1]*q[0]


def dot(p, q):
    return p[0]*q[0] + p[1]*q[1]


def minus(p, q):
    return p[0]-q[0], p[1]-q[1]


def angle_cmp(p, q):
    def half(v):
        return 0 if v[1] > 0 or (v[1] == 0 and v[0] >= 0) else 1
    if half(p) != half(q):
        return -1 if half(p) < half(q) else 1
    z = cross(p, q)
    return -1 if z > 0 else 1 if z < 0 else 0


def upward(x):
    y = x*SCALE
    return Q(-(-y.numerator // y.denominator), SCALE)


def at_piece(piece, g):
    l, u, vl, vu = piece
    return vl + (vu-vl)*(g-l)/(u-l)


def merged(intervals):
    result = []
    for l, u in sorted(intervals):
        require(l < u, "nonpositive interval")
        if result and l <= result[-1][1]:
            result[-1] = result[-1][0], max(result[-1][1], u)
        else:
            result.append((l, u))
    return result


def base_rows(b, k):
    require(K+A*b["left"]-RAD[k]*CAP[k] > 0, "domain boundedness")
    return [(Q(-1), Q(0), Q(0), Q(0)),
            (Q(0), Q(-1), Q(0), Q(0)),
            (-CAP[k], Q(1), Q(0), Q(0)),
            (-b["right"], Q(1), Q(0), Q(-1)),
            (K+A*b["left"], -RAD[k], Q(1), Q(0))]


def check_initial(b):
    lo, hi, g0 = b["left"], b["right"], b["g0"]
    require(1 <= lo < hi <= Q(3, 2), "deficit interval range")
    gate = hi/(K+A*hi)
    require(Q(1, 2) < g0 < Q(31, 50) and gate < g0, "initial source gate")
    scaled = 10000*gate
    require(g0 == Q(scaled.numerator//scaled.denominator+1, 10000),
            "smallest strict initial 1/10000 multiple")
    require(C+(1/g0-A)*hi < T and 1/g0 > A and Q(283, 200) <= 1/g0,
            "complete small-successor payment")
    require(g0/hi > Q(1, 3), "fixed-prime threshold lower bound")


def pair_rows(rows, pair):
    require(len(pair) == 2 and pair[0] != pair[1] and
            all(0 <= i < len(rows) for i in pair), "invalid row pair", pair=pair)
    return rows[pair[0]], rows[pair[1]]


def intersection(rows, pair):
    (a, b, r0, r1), (d, e, s0, s1) = pair_rows(rows, pair)
    det = a*e-b*d
    require(det != 0, "singular claimed row intersection", pair=pair)
    return ((r0*e-b*s0)/det, (r1*e-b*s1)/det,
            (a*s0-r0*d)/det, (a*s1-r1*d)/det)


def master_point(data, pair, g0):
    h0, h1, v0, v1 = intersection(data["rows"], pair)
    return h0/g0+h1, v0/g0+v1


MODES = ("PIECES", "EMPTY", "POLYGONS", "SEGMENTS", "POINTS",
         "FUNCTIONS", "CANDIDATES", "UPPER COVERS")


def parse_conductor(path, family):
    """Parse each record once, with no overwriting of keys or ignored payload."""
    bins, current, mode = [], None, None
    header = re.compile(r"BIN (\d+): ([A-Z0-9]+); \(([^,]+),([^]]+)\]; "
                        r"gamma0=([^;]+); mass ceiling=([^ ]+)\.")
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        CONTEXT.update(line=lineno)
        m = header.fullmatch(raw)
        if m:
            require(current is None or "mass_num" in current, "unterminated bin")
            number, shape, lo, hi, g0, ceiling = m.groups()
            current = dict(number=integer(number), shape=shape,
                           left=rational(lo), right=rational(hi), g0=rational(g0),
                           ceiling=rational(ceiling), sections={}, mode_order=[])
            bins.append(current)
            mode = None
            continue
        require(current is not None and "mass_num" not in current,
                "record outside a bin")
        if raw.startswith("GAMMA "):
            require("gammas" not in current and not current["sections"],
                    "missing, repeated or misplaced gamma dictionary")
            g = {}
            for entry in raw[6:].split(";"):
                i, val = fields(entry, 2, "=")
                i = integer(i)
                require(i not in g, "duplicate gamma index")
                g[i] = rational(val)
            current["gammas"] = g
        elif raw.startswith("MASS SLACK NUMERATOR "):
            require(mode == "UPPER COVERS", "premature mass record")
            current["mass_num"] = integer(raw[len("MASS SLACK NUMERATOR "):])
            mode = None
        else:
            m = re.match(r"^([A-Z ]+) \(", raw)
            if m:
                mode = m[1]
                require(mode in MODES and mode not in current["sections"],
                        "unknown or repeated section header")
                require("gammas" in current, "section precedes gamma dictionary")
                current["mode_order"].append(mode)
                current["sections"][mode] = []
            else:
                require(mode is not None, "unparsed record", record=raw)
                current["sections"][mode].append((lineno, raw))
    require(len(bins) == SPECS[family][0] and
            [b["number"] for b in bins] == list(range(1, len(bins)+1)),
            "missing, repeated or unordered bin")
    for b in bins:
        require("mass_num" in b, "missing final mass")
        wanted = [s for s in MODES if s != "POINTS" or s in b["sections"]]
        require(b["mode_order"] == wanted, "missing or unordered record sections")
        b["sections"].setdefault("POINTS", [])
    return bins


def setup_family(b):
    lo, hi, shape = b["left"], b["right"], b["shape"]
    simple = shape in ("PPC", "PQP")
    b["simple"] = simple
    if simple:
        b["m"] = (1, 1, 2) if shape == "PPC" else (1, 2, 1)
        b["kap"] = (2, 1, 3) if shape == "PPC" else (2, 2, 2)
        if shape == "PPC":
            ranges = [(Q(1), Q(4, 3), (2, Q(18,13), Q(3,2)), (6, -4)),
                      (Q(4,3), Q(26,19), (Q(22,13), Q(18,13), Q(3,2)), (9, -Q(13,2))),
                      (Q(26,19), Q(18,13), (Q(22,13), Q(18,13), Q(19,13)), (9, -Q(13,2)))]
        else:
            ranges = [(Q(1), Q(6,5), (2, Q(4,3), 2), (6, -4)),
                      (Q(6,5), Q(22,17), (2, Q(4,3), 2), (8, -6)),
                      (Q(22,17), Q(4,3), (Q(5,3), Q(4,3), Q(5,3)), (8, -6))]
        matches = [r for r in ranges if r[0] <= lo < hi <= r[1]]
        require(len(matches) == 1, "unearned strict chain boundary")
        _, _, b["uhat"], (b0, b1) = matches[0]
        b["budget"] = b0+b1*lo
        b["inner"] = (0, 1, 2)
        b["weights"] = ((0, 1, 0), (1, 0, 1), (0, 1, 0))
    else:
        ranges = {"LOW": (Q(1), Q(8,7)), "COMMON": (Q(8,7), Q(6,5)),
                  "STAR": (Q(8,7), Q(6,5)), "PATH0": (Q(8,7), Q(6,5)),
                  "PATH1": (Q(6,5), Q(9,7)), "PATH2": (Q(9,7), Q(4,3))}
        require(shape in ranges, "unknown remaining family")
        require(ranges[shape][0] <= lo < hi <= ranges[shape][1],
                "unearned remaining-family boundary")
        b["m"] = (1, 1, 2) if shape == "COMMON" else (1, 1, 1, 1)
        n = len(b["m"])
        w = [[Q(0) for _ in range(n)] for _ in range(n)]
        b["inner"] = (1, 2) if shape.startswith("PATH") else ()
        if shape == "COMMON":
            b["kap"] = (Q(5,3), Q(5,3), Q(8,3))
            b["uhat"] = (Q(8,5), Q(8,5), Q(7,5))
            b["budget"] = 6-5*lo
            for i in range(n):
                for j in range(n):
                    if i != j:
                        w[i][j] = Q(1,3)
        elif shape == "STAR":
            b["kap"], b["uhat"] = (0,2,2,2), (Q(6,5), Q(8,5), Q(8,5), Q(8,5))
            b["budget"] = 6-5*lo
            for i in range(1, n):
                w[0][i] = w[i][0] = Q(1)
        elif shape.startswith("PATH"):
            b["kap"] = (2,1,1,2)
            b["uhat"] = ((Q(5,3),Q(4,3),Q(4,3),Q(5,3)) if shape == "PATH2"
                         else (2,Q(4,3),Q(4,3),2))
            b["budget"] = 6-4*lo if shape == "PATH0" else 8-6*lo
            for i in range(n-1):
                w[i][i+1] = w[i+1][i] = Q(1)
        else:
            b["kap"] = b["uhat"] = (0,)*n
            b["budget"] = 6-4*lo
        b["weights"] = w
    specs = []
    for k in (3, 4, 5):
        if shape == "LOW":
            specs.append((k, None, None))
        else:
            for free in range(len(b["m"])):
                others = [i for i in range(len(b["m"])) if i != free]
                for ends in product(*[(0,1,2) if i in b["inner"] else (0,2) for i in others]):
                    specs.append((k, free, ends))
    counts = {"PPC":81, "PQP":81, "LOW":3, "COMMON":36,
              "STAR":96, "PATH0":180, "PATH1":180, "PATH2":180}
    require(len(specs) == counts[shape], "allocation description enumeration")
    b["specs"] = specs


def data_for(b, d):
    require(0 <= d < len(b["specs"]), "description out of range")
    k, free, ends = b["specs"][d]
    m, kap, u = b["m"], b["kap"], b["uhat"]
    n = len(m)
    r = tuple(sum(row) for row in b["weights"])
    lower = [(Q(k-sum(kap)+kap[i]-r[i]*u[i]), Q(4-m[i]+r[i]), Q(r[i])) for i in range(n)]
    upper = [(Q(kap[i]+sum(b["weights"][i][j]*u[j] for j in range(n))),
              Q(-m[i]-r[i]), Q(-r[i])) for i in range(n)]
    if b["shape"] == "LOW":
        degrees = {3:(2,1,0,0), 4:(2,2,0,0), 5:(2,2,1,0)}[k]
        xs = [(Q(d), Q(0), Q(0)) for d in degrees]
    else:
        xs = [None]*n
        for i, label in zip([i for i in range(n) if i != free], ends):
            xs[i] = ((Q(0),)*3, lower[i], upper[i])[label]
        xs[free] = tuple((Q(k) if j == 0 else Q(0))-
                         sum(xs[i][j] for i in range(n) if i != free) for j in range(3))
    require(tuple(sum(x[j] for x in xs) for j in range(3)) == (k,0,0), "allocation sum")
    rows = base_rows(b, k)
    if b["shape"] != "LOW":
        for i, (a, bb, dd) in enumerate(xs):
            e, j, nn = lower[i]
            v, w, z = upper[i]
            zero, low, up = (-a,-bb,Q(0),dd), (e-a,j-bb,Q(0),dd-nn), (a-v,bb-w,Q(0),z-dd)
            rows.extend((zero, low, up) if b["simple"] else (zero, up))
            if not b["simple"] and i in b["inner"]:
                rows.append(low)
    beta = 1-b["budget"]/k
    require(beta >= Q(1,3), "fixed-prime beta bound")
    return dict(k=k, xs=xs, rows=rows, beta=beta, D=K+A*b["left"])


def triples(b, data, kind, pair):
    if kind == "V":
        h1, v1 = master_point(data, pair, b["g0"])
        h0 = v0 = Q(0)
    else:
        require(kind == "E" and 4 in pair, "moving candidate row")
        h0, h1, v0, v1 = intersection(data["rows"], pair)
    return tuple(sorted((b["m"][i], -data["beta"]*(a*h0+bb*v0)/b["m"][i],
                         1-data["beta"]*(a*h1+bb*v1+dd)/b["m"][i])
                        for i, (a, bb, dd) in enumerate(data["xs"])))


def farkas(b, data, support, weights, numerator):
    require(len(support) == len(weights) and len(set(support)) == len(support)
            and all(0 <= i < len(data["rows"]) for i in support)
            and all(w > 0 for w in weights), "invalid Farkas support or weights")
    normal = tuple(sum(w*data["rows"][i][j] for i,w in zip(support,weights)) for j in (0,1))
    rhs = sum(w*(data["rows"][i][2]/b["g0"]+data["rows"][i][3]) for i,w in zip(support,weights))
    require(normal == (0,0) and rhs.numerator == numerator, "Farkas identity")
    COUNTS["farkas_identities"] += 1
    return rhs


def check_geometry(b, datasets):
    master, used = {}, set()
    for mode in ("EMPTY", "POLYGONS", "SEGMENTS", "POINTS"):
        for lineno, raw in b["sections"][mode]:
            CONTEXT.update(line=lineno)
            values = fields(raw, {"EMPTY":4, "POLYGONS":3, "SEGMENTS":5, "POINTS":4}[mode])
            d = integer(values[0])
            require(d in datasets and d not in used, "missing, repeated or out-of-range description")
            used.add(d)
            data = datasets[d]
            COUNTS[mode.lower()] += 1
            if mode == "EMPTY":
                require(farkas(b,data,ints(values[1]),rats(values[2]),integer(values[3])) < 0,
                        "nonnegative empty-domain contradiction")
                continue
            if mode == "POLYGONS":
                rs, nums = ints(values[1]), ints(values[2])
                require(len(rs) == len(nums) >= 3 and len(set(rs)) == len(rs), "polygon edge list")
                pairs = [(rs[j-1],rs[j]) for j in range(len(rs))]
                ps = [master_point(data,p,b["g0"]) for p in pairs]
                normals = [data["rows"][i][:2] for i in rs]
                require(all(dot(n,n) > 0 for n in normals), "zero polygon normal")
                require(all(cross(normals[j-1],normals[j]) > 0 for j in range(len(rs))), "polygon turns")
                require(sum(angle_cmp(normals[j-1],normals[j]) >= 0 for j in range(len(rs))) == 1,
                        "polygon winding")
                for j, row in enumerate(rs):
                    n = normals[j]
                    tangent = -n[1], n[0]
                    diff = minus(ps[(j+1)%len(rs)],ps[j])
                    lam = dot(diff,tangent)/dot(tangent,tangent)
                    require(lam > 0 and lam.numerator == nums[j] and diff == tuple(lam*x for x in tangent),
                            "polygon edge-length witness")
                    rhs = data["rows"][row][2]/b["g0"]+data["rows"][row][3]
                    require(all(dot(n,p) <= rhs for p in ps), "polygon selected halfplane")
                master[d] = pairs, [(pairs[j],pairs[(j+1)%len(rs)],rs[j]) for j in range(len(rs))]
            else:
                support = ints(values[1])
                require(farkas(b,data,support,rats(values[2]),0) == 0, "nonzero equality certificate")
                if mode == "POINTS":
                    pair = ints(values[3])
                    require(all(i in support for i in pair), "point normals outside equality support")
                    master_point(data,pair,b["g0"])
                    master[d] = [pair], []
                else:
                    row = integer(values[3])
                    pairs = [ints(s) for s in fields(values[4],2,":")]
                    require(row in support and all(row in p for p in pairs), "segment equality line")
                    ps = [master_point(data,p,b["g0"]) for p in pairs]
                    n = data["rows"][row][:2]
                    require(dot(n,n) > 0, "zero segment normal")
                    # Simple-chain records orient by q-p; remaining records
                    # also allow coincident endpoints, using the line tangent.
                    direction = minus(ps[1],ps[0]) if ps[0] != ps[1] else (-n[1],n[0])
                    require(ps[0] != ps[1] or not b["simple"], "simple segment collapsed to a point")
                    ends = [p[1] if p[0] == row else p[0] for p in pairs]
                    require(dot(data["rows"][ends[0]][:2],direction) < 0 <
                            dot(data["rows"][ends[1]][:2],direction), "segment endpoint orientation")
                    if ps[0] == ps[1]:
                        COUNTS["coincident_endpoint_points"] += 1
                        master[d] = [pairs[0]], []
                    else:
                        master[d] = pairs, [(pairs[0],pairs[1],row)]
    require(used == set(datasets), "omitted description")
    return master


def expected_candidates(b, datasets, master):
    expected = {}
    for d, (pairs, edges) in master.items():
        data = datasets[d]
        def inverse_j(pair):
            h, v = master_point(data,pair,b["g0"])
            j = data["D"]*h-RAD[data["k"]]*v
            require(j > 0, "nonpositive master clipping value")
            return 1/j
        candidates = []
        def include(kind, pair, left, right):
            if left < right:
                candidates.append((kind,tuple(sorted(pair)),left,right))
            elif left == right:
                COUNTS["isolated_clipping_intervals"] += 1
            else:
                COUNTS["empty_clipping_intervals"] += 1
        for pair in pairs:
            include("V",pair,b["g0"],min(END,inverse_j(pair)))
        for p, q, row in edges:
            jp, jq = inverse_j(p), inverse_j(q)
            if row != 4 and jp != jq:
                include("E",(4,row),max(b["g0"],min(jp,jq)),min(END,max(jp,jq)))
            else:
                COUNTS["parallel_or_existing_clipping_edges"] += 1
        expected[d] = Counter(candidates)
    return expected


def check_conductor_bin(b):
    check_initial(b)
    setup_family(b)
    g = b["gammas"]
    require(list(g) == list(range(len(g))) and g[0] == b["g0"] and g[len(g)-1] == END
            and all(g[i] < g[i+1] for i in range(len(g)-1)), "gamma dictionary coverage")
    pieces = []
    for lineno, raw in b["sections"]["PIECES"]:
        CONTEXT.update(line=lineno)
        for entry in raw.split(";"):
            i,j,vl,vu = ints(entry)
            require(i in g and j in g and j == i+1 and vl >= 0 and vu >= 0, "piece index or value")
            pieces.append((g[i],g[j],Q(vl,SCALE),Q(vu,SCALE)))
    require([(p[0],p[1]) for p in pieces] == [(g[i],g[i+1]) for i in range(len(g)-1)],
            "missing, repeated or unordered piece")
    datasets = {d:data_for(b,d) for d in range(len(b["specs"]))}
    master = check_geometry(b,datasets)
    expected = expected_candidates(b,datasets,master)
    functions, reps = {}, {}
    for lineno, raw in b["sections"]["FUNCTIONS"]:
        CONTEXT.update(line=lineno)
        fid,d,kind,pair = fields(raw,4)
        fid,d,pair = integer(fid),integer(d),ints(pair)
        require(fid not in functions and fid >= 0 and d in master, "function index or representative")
        functions[fid] = triples(b,datasets[d],kind,pair)
        reps[fid] = d,kind,tuple(sorted(pair))
    require(list(functions) == list(range(len(functions))), "function index coverage")
    require(len(set(functions.values())) == len(functions), "duplicate unmerged function")
    intervals, seen, representations = defaultdict(list), set(), defaultdict(set)
    for lineno, raw in b["sections"]["CANDIDATES"]:
        CONTEXT.update(line=lineno)
        ds, entries = fields(raw,2,":")
        d = integer(ds)
        require(d in master and d not in seen, "candidate description missing, repeated or out of range")
        seen.add(d)
        actual = []
        for entry in entries.split(";"):
            kind,pair,fid,left,right = fields(entry,5,"/")
            pair,fid,left,right = ints(pair),integer(fid),integer(left),integer(right)
            require(fid in functions and left in g and right in g and left < right, "candidate index")
            require(triples(b,datasets[d],kind,pair) == functions[fid], "candidate function mismatch")
            actual.append((kind,tuple(sorted(pair)),g[left],g[right]))
            intervals[fid].append((g[left],g[right]))
            representations[fid].add((d,kind,tuple(sorted(pair))))
        require(Counter(actual) == expected[d], "clipped candidate coverage mismatch",
                description=d, missing=str(expected[d]-Counter(actual)), extra=str(Counter(actual)-expected[d]))
        COUNTS["candidates"] += len(actual)
    require(seen == {d for d,items in expected.items() if items}, "omitted candidate description")
    require(all(reps[fid] in representations[fid] for fid in functions), "unlisted function representative")
    covers, cover_keys = defaultdict(list), set()
    for lineno, raw in b["sections"]["UPPER COVERS"]:
        CONTEXT.update(line=lineno)
        fid,kind,ends,bounds,lwit,rwit = fields(raw,6)
        fid = integer(fid)
        i,j = ints(ends)
        require(fid in functions and i in g and j in g and i < j, "cover index or domain")
        left,right = g[i],g[j]
        key = fid,left,right
        require(key not in cover_keys, "duplicate or conflicting upper cover")
        cover_keys.add(key)
        vals = rats(bounds)
        if kind == "C":
            require(len(vals) == 1, "constant cover syntax")
            vals = vals[0],vals[0]
            overlap = [p for p in pieces if max(left,p[0]) < min(right,p[1])]
            require(overlap and b["g0"] <= left < right <= END, "cover outside partition")
            for p in overlap:
                require(vals[0] <= min(at_piece(p,max(left,p[0])),at_piece(p,min(right,p[1]))),
                        "constant exceeds interpolant")
        else:
            require(kind == "L" and len(vals) == 2, "linear cover syntax")
            containers = [p for p in pieces if p[0] <= left < right <= p[1]]
            require(len(containers) == 1, "linear cover crosses piece boundary")
            require(vals == (at_piece(containers[0],left),at_piece(containers[0],right)),
                    "linear chord mismatch")
        for endpoint,bound,wraw in zip((left,right),vals,(lwit,rwit)):
            witness = ints(wraw)
            require(len(witness) == len(b["m"])+3, "witness field count")
            e,v = witness[:2]
            us,slack = witness[2:-1],witness[-1]
            require(0 <= e <= 30, "witness exponent")
            scale = 10**e
            aff = functions[fid]
            terms = [max(Q(0),aa+bb*endpoint) for m,aa,bb in aff]
            require(tuple(m for m,aa,bb in aff) == tuple(sorted(b["m"])), "sorted witness weights")
            require(v <= scale*bound < v+1, "floor witness mismatch")
            require(all(ui-1 < scale*x <= ui for ui,x in zip(us,terms)), "ceiling witness mismatch")
            squares = sum(m*ui*ui for (m,aa,bb),ui in zip(aff,us))
            require(2*scale*v-squares == slack and slack >= 0, "rounding slack identity")
            exact = sum(Q(m,2)*x*x for (m,aa,bb),x in zip(aff,terms))
            require(exact <= Q(squares,2*scale*scale) <= Q(v,scale) <= bound, "endpoint bound")
            COUNTS["rounding_witnesses"] += 1
        covers[fid].append((left,right))
    require(set(functions) == set(intervals) == set(covers), "function/candidate/cover mismatch")
    for fid in functions:
        # Covers may share endpoints, but must not duplicate open intervals.
        ordered = sorted(covers[fid])
        require(all(x[1] <= y[0] for x,y in zip(ordered,ordered[1:])), "overlapping covers")
        require(merged(intervals[fid]) == merged(covers[fid]), "uncovered or extraneous function interval")
    mass = 4*b["g0"]**3+3*sum((u-l)*(vl+vu) for l,u,vl,vu in pieces)
    diff = b["ceiling"]-mass
    require(diff > 0 and diff.numerator == b["mass_num"], "mass slack numerator")
    require(b["ceiling"] < TARGET, "mass ceiling fails strict target")
    return dict(bin=b["number"],shape=b["shape"],deficit=[b["left"],b["right"]],
                initial_order=b["g0"],pieces=len(pieces),descriptions=len(datasets),
                empty=len(b["sections"]["EMPTY"]),polygons=len(b["sections"]["POLYGONS"]),
                segments=len(b["sections"]["SEGMENTS"]),points=len(b["sections"]["POINTS"]),
                functions=len(functions),candidates=sum(map(len,intervals.values())),
                covers=len(cover_keys),rounding_witnesses=2*len(cover_keys),mass=mass,
                ceiling=b["ceiling"],mass_slack=diff,gap_to_target=TARGET-b["ceiling"])


def check_ranges(bins, wanted):
    require(set(b["shape"] for b in bins) == set(wanted), "missing or unexpected configuration")
    for shape,(left,right) in wanted.items():
        bs = [b for b in bins if b["shape"] == shape]
        require(bs and bs[0]["left"] == left and bs[-1]["right"] == right
                and all(x["right"] == y["left"] for x,y in zip(bs,bs[1:])),
                "missing or overlapping deficit coverage", shape=shape)


def check_conductor(path, family):
    bins = parse_conductor(path,family)
    wanted = ({"PPC":(Q(1),Q(18,13)),"PQP":(Q(1),Q(4,3))} if family == "simple-chains" else
              {"LOW":(Q(1),Q(8,7)),"COMMON":(Q(8,7),Q(6,5)),"STAR":(Q(8,7),Q(6,5)),
               "PATH0":(Q(8,7),Q(6,5)),"PATH1":(Q(6,5),Q(9,7)),"PATH2":(Q(9,7),Q(4,3))})
    check_ranges(bins,wanted)
    receipts = []
    for b in bins:
        CONTEXT.update(bin=b["number"],shape=b["shape"])
        receipts.append(check_conductor_bin(b))
    require(sum(r["pieces"] for r in receipts) == SPECS[family][1], "total piece coverage")
    require(sum(r["rounding_witnesses"] for r in receipts) ==
            (14504 if family == "simple-chains" else 15062), "total rounding witness coverage")
    require(max(r["ceiling"] for r in receipts) == SPECS[family][2], "published maximum ceiling")
    return receipts


def parse_label(text):
    obj = literal_eval(text)
    if obj is None:
        return None
    require(isinstance(obj,(tuple,list)) and len(obj) == 3 and
            type(obj[0]) is int and type(obj[1]) is int and
            isinstance(obj[2],(tuple,list)) and len(obj[2]) == 2 and
            all(type(x) is int for x in obj[2]), "malformed maximizing label")
    return obj[0],obj[1],tuple(obj[2])


def parse_two(path, family):
    lines = path.read_text(encoding="utf-8").splitlines()
    bins, pieces = [], defaultdict(list)
    phase = 0
    title = ("Exact whole-source convex certificate. Every number is rational; maxima include all feasible vertices."
             if family == "two-quadrics" else
             "Exact cubic-plus-plane whole-source convex certificate. Every number is rational; all feasible vertices are included.")
    require(lines and lines[0] == title, "two-component data title")
    for lineno, raw in enumerate(lines[1:],2):
        if not raw.strip():
            continue
        CONTEXT.update(line=lineno)
        if raw.startswith("Full convex-piece endpoint bounds."):
            require(phase == 1, "duplicate or misplaced piece preface")
            phase = 2
            continue
        require(raw.startswith("|") and raw.endswith("|"), "unknown two-component record")
        cells = [c.strip() for c in raw[1:-1].split("|")]
        if cells[0] == "bin":
            expected = ["bin", "b lower", "b upper", "initial gamma"]
            if phase == 0:
                if family == "cubic-plane":
                    expected.append("A_C")
                expected += ["normalized capacity upper","pieces"]
                phase = 1
            else:
                require(phase == 2, "duplicate table header")
                expected = ["bin","left gamma","right gamma","left bound","right bound",
                            "left maximizing vertex","right maximizing vertex","active vertices"]
                phase = 3
            require(cells == expected, "table column definitions")
            continue
        if all(c == "---" for c in cells):
            require(phase in (1,3) and len(cells) == (8 if phase == 3 else 6+(family == "cubic-plane")),
                    "separator columns")
            continue
        number = integer(cells[0])
        if phase == 1:
            require(len(cells) == 6+(family == "cubic-plane"), "aggregate column count")
            b = dict(number=number,shape=family,left=rational(cells[1]),right=rational(cells[2]),
                     g0=rational(cells[3]),ceiling=rational(cells[-2]),count=integer(cells[-1]))
            if family == "cubic-plane":
                b["AC"] = rational(cells[4])
                require(b["AC"] == 4+min(Q(2),6-3*b["left"]), "cubic component cap")
            bins.append(b)
        else:
            require(phase == 3 and len(cells) == 8, "piece column count")
            pieces[number].append((lineno,*map(rational,cells[1:5]),
                                   parse_label(cells[5]),parse_label(cells[6]),integer(cells[7])))
    require(len(bins) == SPECS[family][0] and [b["number"] for b in bins] == list(range(1,len(bins)+1)),
            "aggregate bin coverage")
    require(set(pieces) == set(range(1,len(bins)+1)), "piece bin coverage")
    return bins,pieces


def two_rows(b, k, region, family):
    rows = base_rows(b,k)
    if family == "two-quadrics":
        rows.append((Q(k,2)-Q(9,2),Q(3),Q(0),Q(-1)))
        rows.append((Q(k)-Q(9,2),Q(3),Q(0),Q(-1)) if region == 0 else
                    (Q(9,2)-k,Q(-3),Q(0),Q(1)))
    else:
        ac,ap = b["AC"],Q(24,7)
        rows.extend([(-ac,Q(4),Q(0),Q(-1)),(-ap,Q(2),Q(0),Q(-1)),
                     (k-ac-ap,Q(6),Q(0),Q(-2))])
        rows.append([(k-ac,Q(4),Q(0),Q(-1)),(k-ap,Q(2),Q(0),Q(-1)),
                     (ac-k,Q(-4),Q(0),Q(1)),(ap-k,Q(-2),Q(0),Q(1))][region])
    return rows


def enumerate_two(b, family):
    candidates = []
    for k in (3,4,5):
        beta = 1-(6-4*b["left"])/k
        require(beta >= Q(1,3), "two-component beta")
        for region in range(2 if family == "two-quadrics" else 4):
            rows = two_rows(b,k,region,family)
            COUNTS["two_component_domains"] += 1
            for pair in combinations(range(len(rows)),2):
                COUNTS["all_boundary_row_pairs"] += 1
                r,s = pair_rows(rows,pair)
                if cross(r,s) == 0:
                    COUNTS["parallel_row_pairs"] += 1
                    continue
                h0,h1,v0,v1 = intersection(rows,pair)
                left,right = b["g0"],END
                excluded = False
                for a,bb,r0,r1 in rows:
                    slope,rhs = a*h1+bb*v1-r1,r0-a*h0-bb*v0
                    if slope > 0:
                        right = min(right,rhs/slope)
                    elif slope < 0:
                        left = max(left,rhs/slope)
                    elif rhs < 0:
                        excluded = True
                if excluded or left > right:
                    COUNTS["infeasible_row_pairs"] += 1
                    continue
                COUNTS["isolated_row_pairs" if left == right else "feasible_row_pair_intervals"] += 1
                candidates.append(dict(label=(k,region,pair),left=left,right=right,
                                       coords=(h0,h1,v0,v1),rows=rows,beta=beta))
    return candidates


def two_value(b, candidate, gamma, family):
    h0,h1,v0,v1 = candidate["coords"]
    h,v = h0+h1*gamma,v0+v1*gamma
    require(h > 0 and all(a*h+bb*v <= r0+r1*gamma for a,bb,r0,r1 in candidate["rows"]),
            "two-component endpoint feasibility")
    k,region,_ = candidate["label"]
    beta = candidate["beta"]
    if family == "two-quadrics":
        if region == 0:
            return gamma**2+max(Q(0),gamma-k*beta*h/2)**2
        s1 = Q(9,2)*h-gamma-3*v
        s2 = (k-Q(9,2))*h+gamma+3*v
        require(s1 >= 0 and s2 >= 0, "two-quadric allocation signs")
        return max(Q(0),gamma-beta*s1/2)**2+max(Q(0),gamma-beta*s2/2)**2
    ac,ap = b["AC"],Q(24,7)
    sc,sp = [(k*h,Q(0)),(Q(0),k*h),(ac*h-gamma-4*v,(k-ac)*h+gamma+4*v),
             ((k-ap)*h+gamma+2*v,ap*h-gamma-2*v)][region]
    require(sc >= 0 and sp >= 0 and sc+sp == k*h and sc <= ac*h-gamma-4*v and sp <= ap*h-gamma-2*v,
            "cubic-plane full allocation")
    return Q(3,2)*max(Q(0),gamma-beta*sc/3)**2+Q(1,2)*max(Q(0),gamma-beta*sp)**2


def check_two(path, family):
    bins,pieces = parse_two(path,family)
    check_ranges(bins,{family:(Q(1),Q(3,2) if family == "two-quadrics" else Q(10,7))})
    receipts = []
    for b in bins:
        CONTEXT.update(bin=b["number"],shape=family)
        check_initial(b)
        candidates = enumerate_two(b,family)
        cuts = {b["g0"],END}
        for item in candidates:
            cuts.update((item["left"],item["right"]))
        j = 0
        while b["g0"]+Q(j,200) < END:
            cuts.add(b["g0"]+Q(j,200))
            j += 1
        cuts = sorted(cuts)
        table = pieces[b["number"]]
        require(len(table) == b["count"] and [(p[1],p[2]) for p in table] == list(zip(cuts,cuts[1:])),
                "missing, repeated or unordered source piece")
        mass = 4*b["g0"]**3
        for lineno,l,u,vl,vu,ll,lu,count in table:
            CONTEXT.update(line=lineno)
            active = [c for c in candidates if c["left"] < (l+u)/2 < c["right"]]
            require(len(active) == count and all(c["left"] <= l < u <= c["right"] for c in active),
                    "complete fixed active list")
            for endpoint,bound,label in ((l,vl,ll),(u,vu,lu)):
                values = [(two_value(b,c,endpoint,family),c["label"]) for c in active]
                maximum = max((v for v,_ in values),default=Q(0))
                require(bound == upward(maximum), "upward endpoint bound")
                require((label is None and not active) or
                        (label in [lab for v,lab in values if v == maximum]), "maximizing label")
                COUNTS["competing_endpoint_evaluations"] += len(values)
                COUNTS["two_component_endpoint_bounds"] += 1
            COUNTS["empty_source_pieces"] += not active
            mass += 3*(u-l)*(vl+vu)
        require(b["ceiling"] == upward(mass) and b["ceiling"] < TARGET, "aggregate rounding and margin")
        receipts.append(dict(bin=b["number"],shape=family,deficit=[b["left"],b["right"]],
                             initial_order=b["g0"],pieces=len(table),candidates=len(candidates),
                             isolated_pairs=sum(c["left"] == c["right"] for c in candidates),
                             mass=mass,ceiling=b["ceiling"],gap_to_target=TARGET-b["ceiling"]))
    require(sum(r["pieces"] for r in receipts) == SPECS[family][1], "total two-component pieces")
    require(max(r["ceiling"] for r in receipts) == SPECS[family][2], "maximum two-component ceiling")
    return receipts


def validate(data_dir):
    wanted = {name+".txt" for name in SPECS}
    require(data_dir.is_dir() and {p.name for p in data_dir.iterdir()} == wanted,
            "missing or unexpected data file")
    output = []
    for family in SPECS:
        CONTEXT.clear()
        CONTEXT.update(family=family)
        COUNTS.clear()
        checker = check_two if family in ("two-quadrics","cubic-plane") else check_conductor
        receipts = checker(data_dir/(family+".txt"),family)
        output.append(dict(family="11-quartic-threefold/"+family,status="PASS",bins=receipts,
                           counts=dict(COUNTS),bin_count=len(receipts),
                           pieces=sum(r["pieces"] for r in receipts),
                           maximum_mass_ceiling=max(r["ceiling"] for r in receipts),
                           minimum_gap_to_target=min(r["gap_to_target"] for r in receipts)))
        print(json.dumps({k:output[-1][k] for k in ("family","status","bin_count","pieces",
                                                   "maximum_mass_ceiling","minimum_gap_to_target")},default=str),flush=True)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data",type=Path,required=True)
    parser.add_argument("--report",type=Path,required=True)
    args = parser.parse_args()
    report = dict(status="FAIL",checker_format="quartic-exact-v1",arithmetic="fractions.Fraction",
                  python_version=sys.version.split()[0],
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  interpretation="finite certificate validation; geometric applicability is proved in the article")
    try:
        report["families"] = validate(args.data)
        report["data_files"] = {name+".txt": hashlib.sha256((args.data/(name+".txt")).read_bytes()).hexdigest()
                                for name in SPECS}
        report["status"] = "PASS"
    except Exception as exc:
        report["error"] = str(exc)
        report["context"] = dict(CONTEXT)
        print("FAIL: "+str(exc),file=sys.stderr,flush=True)
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2,default=str)+"\n",encoding="utf-8")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
