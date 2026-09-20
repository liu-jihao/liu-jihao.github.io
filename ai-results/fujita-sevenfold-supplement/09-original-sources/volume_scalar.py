"""Independent polynomial and complete scalar-tail checks, adapted from the volume-two audit."""
from fractions import Fraction as F
from math import comb

def trim(p):
    p = list(map(F, p))
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(*ps):
    out = [F(0)] * max(map(len, ps))
    for p in ps:
        for i, x in enumerate(p):
            out[i] += x
    return trim(out)


def scale(p, x):
    return trim([F(x) * c for c in p])


def mul(*ps):
    out = [F(1)]
    for p in ps:
        nxt = [F(0)] * (len(out) + len(p) - 1)
        for i, a in enumerate(out):
            for j, b in enumerate(p):
                nxt[i + j] += a * b
        out = trim(nxt)
    return out


def power(p, n):
    return mul(*([p] * n)) if n else [F(1)]


def mobius(p, lo, hi, degree):
    # (1+y)^degree P((lo+hi*y)/(1+y)).
    return add(*(scale(mul(power([lo, hi], i),
                          power([1, 1], degree - i)), c)
                 for i, c in enumerate(p)))


def shift(p, lo):
    return add(*(scale(power([lo, 1], i), c) for i, c in enumerate(p)))


def audit(volume, rows):
    checks, details = [], {}
    def check(name, condition, value=None):
        checks.append({'name':name, 'pass':bool(condition), 'value':str(value)})
        if not condition:
            raise ValueError(f'{name}: {value}')
    def eq(name, actual, expected):
        check(name, actual == expected, actual)
    r, Q = F(29, 25), F(13, 5)
    kap, beta, q0, C0, D0 = F(37, 40), F(52, 73), F(1617, 1825), F(9, 50), F(181, 100)
    pD, pH = F(43, 50), F(1267, 1500)
    Z, W, BK = [0, 1], [F(7, 60), F(7, 6)], [beta, -1]
    ZK, WK = [kap, 1], add(W, [kap])

    pd_num = add(
        scale(mul(Z, power(BK, 6), power(ZK, 6)), D0 * beta),
        scale(mul(power([pD, r], 7), power(BK, 6)), -beta * kap**6),
        mul(power([beta*pD, -q0], 7), power(ZK, 6)))
    ph_num = add(
        scale(mul(Z, W, power(BK, 6), power(WK, 6)), D0 * beta),
        scale(mul(Z, power(add([pH], scale(W, r)), 7), power(BK, 6)), -beta * kap**6),
        scale(mul(add(Z, scale(W, -1)), power(BK, 6), power(WK, 6)), beta * pH**7),
        mul(W, power([beta*pH, -q0], 7), power(WK, 6)))
    pl = add(
        scale(mul(Z, W, power(WK, 6)), D0),
        scale(mul(Z, power(add([pH], scale(W, r)), 7)), -kap**6),
        scale(mul(add(Z, scale(W, -1)), power(WK, 6)), pH**7))
    eq("PD numerator divisible by z", pd_num[0], F(0))
    eq("PH numerator divisible by z", ph_num[0], F(0))
    pd, ph = pd_num[1:], ph_num[1:]
    eq("degrees PD PH PL", (len(pd)-1, len(ph)-1, len(pl)-1), (12, 13, 8))

    coefficient_count = 0
    for item in volume['bounded']:
        kind, (lo, hi), n, multiplier = item['kind'], item['interval'], item['degree'], item['scale']
        calculated = scale(mobius(pd if kind == 'D' else ph, lo, hi, n), multiplier)
        eq(f'P{kind} on [{lo},{hi}] complete coefficient identity', calculated, item['coefficients'])
        check(f'P{kind} all coefficients positive', all(c > 0 for c in calculated))
        coefficient_count += len(calculated)
    item = volume['unbounded']
    calculated = scale(shift(pl, item['shift']), item['scale'])
    eq('PL unbounded complete coefficient identity', calculated, item['coefficients'])
    check('PL all coefficients positive', all(c > 0 for c in calculated))
    coefficient_count += len(calculated)
    details['polynomial_coefficients_independently_matched'] = coefficient_count


    def U(p, z):
        return kap**6 * (p+r*z)**7 / (z*(z+kap)**6)


    endpoint_diffs = {
        "divisor_a_le_1": F(199, 100) - (C0+7*r*pD**6-6*pD**7/kap+q0**7/beta),
        "double_a_le_1": F(199, 100) - (C0+U(pH,F(7,60))-pH**7/F(7,60)+q0**7/beta),
        "divisor_large_at_7_10": F(199, 100)-(C0+U(pD,F(7,10))),
        "divisor_large_limit": F(199, 100)-(C0+kap**6*r**7),
    }
    for name, actual in endpoint_diffs.items():
        stated = volume['endpoint_differences'][name]
        eq(name + " printed exact difference", actual, stated)
        check(name + " positive", actual > 0)
    details["endpoint_positive_differences"] = {k: str(v) for k, v in endpoint_diffs.items()}
    tau0 = F(7, 6)
    A0, B0 = (5*tau0+tau0**(-5))/6, F(5,3)*(tau0-tau0**(-5))
    eq("kappa sixth rounding", kap**6-F(61,100), F(67166409,4096000000))
    eq("prefix correction rounding", C0-(1-kap**6)*F(9,10)**7, F(53548852728321,40960000000000000))
    eq("AMGM intercept rounding", F(21,20)-A0, F(1009,1512630))
    eq("AMGM slope rounding", B0-F(7,6), F(1009,151263))
    check("all rounding directions strict", kap**6 > F(61,100) > F(14,25) and C0 > (1-kap**6)*F(9,10)**7 and A0 < F(21,20) and B0 > F(7,6))
    for label, p in [("D", pD), ("H", pH)]:
        check(label + " crossing location assumptions", F(9,10)*kap < p < kap*r and p < q0)
    eq("divisor derivative coefficient", 6*r*kap-7*pD, F(209,500))
    eq("divisor derivative constant", pD*kap, F(1591,2000))
    check("bounded denominators positive", F(7,10) < beta and F(17,25) < beta)
    eq("reduced anchor", F(9,8)*(r-2*F(257,1000)), F(2907,4000))
    check("reduced domination", F(2907,4000) < pD and F(9,8)**6 > 2)

    entry = F(9,10)**7+F(61,100)*(F(26,25)**7-F(9,10)**7)+F(14,25)*(r**7-F(26,25)**7)
    eq("entry cumulative loss", entry, F(143360580567427,78125000000000))
    eq("volume two entry floor", 2-entry, F(12889419432573,78125000000000))
    check("entry positive", 0 < 2-entry and entry < F(46,25))
    details["entry_floor"] = str(2-entry)
    details["old_independent_endpoint_loss"] = str(entry+q0**7/beta)
    details["old_early_replacement_endpoint_loss"] = str(entry+q0**7/beta-F(11,50)*(r**7-F(26,25)**7))


    byid = {x["id"]:x for x in rows}


    def moment(c, m):
        if c == 0:
            return 0
        ans, j = 0, 1
        while comb(c+j-1,c) < m:
            ans += m-comb(c+j-1,c)
            j += 1
        return ans


    expected_states = {(1,1,1)}
    for d in range(2,7):
        for e in range(d,8):
            c = e-d
            if c == 0:
                expected_states.add((d,e,1))
                continue
            m = c+1
            while 2*moment(c,m) < d*m:
                expected_states.add((d,e,m))
                m += 1
    check("table exhausts all positive-deficit states", {(x["d"],x["e"],x["m"]) for x in rows} == expected_states)
    for row in rows:
        eq(f"row {row['id']} moment", moment(row["e"]-row["d"],row["m"]), row["moment"])
        eq(f"row {row['id']} deficit cap", F(row["d"])-2*F(row["moment"],row["m"]), row["cap"])
        rr = row["radius"]
        check(f"row {row['id']} rounded radius", rr**row["d"] >= row["m"] and (rr-F(1,1000))**row["d"] < row["m"])


    def envelope(row, x):
        return min(a*x+b for a,b in row["lines"])


    def endpoints(row, cap):
        points = {F(0),cap}
        for a,p in row["lines"]:
            for aa,pp in row["lines"]:
                if a != aa:
                    q = (pp-p)/(a-aa)
                    if 0 < q < cap:
                        points.add(q)
        return points


    tail_count = 0
    for row in rows:
        for a,p in row["lines"]:
            extra = F(0)
            for nxt in rows:
                if nxt["d"] >= row["d"] or nxt["e"] > row["e"]:
                    continue
                cap = min(row["cap"],nxt["cap"])
                if nxt["d"] == row["d"]-1:
                    cap = min(cap, row["cap"]-F(nxt["m"],row["m"]))
                if cap <= 0:
                    continue
                extra = max(extra, max(envelope(nxt,y)-max(a,row["radius"])*y for y in endpoints(nxt,cap)))
            rhs = max(row["radius"]-a,0)*row["cap"]+extra
            check(f"tail recurrence row {row['id']} line {a},{p}", p >= rhs, p-rhs)
            tail_count += 1
    details["complete_tail_affine_lines_checked"] = tail_count


    def failure_max(selected, order, truncate=None):
        return max(((7-b)/order+envelope(row,b),row["id"])
                   for row in selected for b in endpoints(row,min(row["cap"],truncate) if truncate is not None else row["cap"]))


    u0 = F(26,25)
    eq("early dimensions below four exclusion", failure_max([x for x in rows if x["d"]<=3],u0)[0], F(20773,2600))
    eq("early fourfold deficit gate", failure_max([x for x in rows if x["d"]==4],u0,F(3,2))[0], F(41553,5200))
    eq("early fivefold divisor deficit gate", failure_max([x for x in rows if x["d"]>=5],u0,F(7,5))[0], F(38951,4875))
    eligible = [x for x in rows if x["d"]==4 and x["cap"]>F(3,2)]
    eq("nine early fourfold candidates", [x["id"] for x in eligible], [28,29,30,32,33,34,39,40,41])
    early_lines = {
        28:(F(1),F(757,1000)),29:(F(119,100),F(287,500)),
        30:(F(1317,1000),F(171,500)),32:(F(1317,1000),F(71,125)),
        33:(F(283,200),F(379,1000)),34:(F(187,125),F(297,1000)),
        39:(F(283,200),F(29,50)),40:(F(187,125),F(17,40)),
        41:(F(783,500),F(331,1000)),
    }
    for rid,(a,p) in early_lines.items():
        check(f"early row {rid} line exists", (a,p) in byid[rid]["lines"])
        check(f"early row {rid} source bound increasing in deficit", a*u0 > 1)
    for rid in [28,29,30,34]:
        row,a_p = byid[rid],early_lines[rid]
        total = (7-row["cap"])/u0+a_p[0]*row["cap"]+a_p[1]
        check(f"early row {rid} paid with strict gap", total<8,8-total)
    slacks = {32:F(481,1875),33:F(217,12500),39:F(2311,5000),40:F(2921,12500),41:F(2699,37500)}
    for rid,stated in slacks.items():
        row,(a,p) = byid[rid],early_lines[rid]
        cap_total=a*row["cap"]+p
        eq(f"early row {rid} exact slack", 7-row["cap"]+u0*(cap_total-8), stated)
        check(f"early row {rid} slack hypotheses", cap_total<8 and 0<stated<F(1,2))
    for rid,stated in [(32,F(9116,4621)),(33,F(11573,5895))]:
        a,p=early_lines[rid]
        b=((8-p)*u0-7)/(a*u0-1)
        eq(f"early linear-normal row {rid} lower deficit", b,stated)
        check(f"early linear-normal row {rid} strict gate and monotonicity", b>F(19,10) and 7*a+p-8>0)

    # Retained E_1=6 root: test its proposed whole-tail line directly against
    # every actual successor permitted by the root cap, including jumps.
    a, intercept = F(187,125), F(1697,5000)
    for nxt in rows:
        if nxt["d"]>=4:
            continue
        cap=min(F(8,5),nxt["cap"])
        if nxt["d"]==3:
            cap=min(cap,F(8-nxt["m"],5))
        if cap<=0:
            continue
        excess=max(envelope(nxt,y)-a*y for y in endpoints(nxt,cap))
        check(f"E1 six full root successor row {nxt['id']}", excess<=intercept,excess)
    total=F(135,26)+F(1496,625)+intercept
    eq("E1 six exact payment",total,F(103029,13000))
    eq("E1 six strict gap",8-total,F(971,13000))
    for e,m,h,first in [(6,3,[1,2],2),(6,4,[1,2,1],4),(7,4,[1,3],3),(7,5,[1,3,1],5),(7,6,[1,3,2],7)]:
        check(f"early row {e},{m} full Hilbert numerator arithmetic",sum(h)==m and h[1]==e-4 and sum(i*x for i,x in enumerate(h))==first and len(h)-1<=3)

    for rid,stated in [(51,F(28348,8321)),(55,F(28725,11134)),(56,F(3467,1328)),(57,F(37483,15020))]:
        a,p=byid[rid]["lines"][0]
        deficit=((8-p)*r-7)/(a*r-1)
        eq(f"profile source row {rid} exact deficit",deficit,stated)
        h=F(8,25) if rid==51 else F(257,1000)
        check(f"profile source row {rid} uniform height margin",r/(7-deficit)-h>F(1,100000),r/(7-deficit)-h)
    a,p=F(5,4),F(187,250)
    smooth_b=((8-p)*r-7)/(a*r-1)
    eq("smooth divisor improved lower deficit",smooth_b,F(17654,5625))
    eq("smooth divisor strict height gap",r/(7-smooth_b)-F(3,10),F(87,217210))
    check("smooth divisor uniform margin",r/(7-smooth_b)-F(3,10)>F(1,100000))
    check("all singular divisor margins",all(e*r/7-e*F(1,100000)>F(3,10) for e in range(2,7)))
    check("source r only four nondivisorial rows remain",failure_max([x for x in rows if x["d"]<=5 and x["id"] not in [51,55,56,57]],r)[0]<8)
    check("all four fivefold rows excluded at 34/25",failure_max([byid[i] for i in [51,55,56,57]],F(34,25))[0]<8)

    full_rate=6*(1-2*F(9,50))**6-8*(1-3*F(9,50))**6+3*(1-4*F(9,50))**6
    mixed_rate=(1-F(19,100))**6+(1-2*F(19,100))**6-(1-3*F(19,100))**6
    pair_rates=[(1-F(7,40))**6,(1-2*F(7,40))**5*(1+10*F(7,40))]
    eq("full fourfold exact rate",full_rate,F(5280752440,15625000000))
    eq("linear-normal exact rate",mixed_rate,F(332908409016,1000000000000))
    eq("pair endpoint exact rate",pair_rates[0],F(1291467969,4096000000))
    eq("pair balanced exact rate",pair_rates[1],F(4084223,12800000))
    check("all early rates below 17/50",max([full_rate,mixed_rate]+pair_rates)<F(17,50))
    check("full fourfold strict threshold",F(208,1155)>F(9,50))
    check("linear-normal threshold domain",F(1040,5355)>F(19,100)>F(2,11) and F(2,11)==F(4,3*7+1))
    check("two-normal threshold with both branches",F(26,147)>F(7,40) and F(7,40)<F(4,21))
    eq("endpoint full-tail largest cost",failure_max(rows,Q)[0],F(12949,1625))
    check("endpoint error and 1/100 margin",F(12949,1625)+F(1,1000)<8-F(1,100))
    details["early_rates"]={"full_embedding":str(full_rate),"linear_normal":str(mixed_rate),"pair_endpoint":str(pair_rates[0]),"pair_balanced":str(pair_rates[1])}

    eq('four delivered rounding witnesses', volume['rounding_differences'],
       [kap**6-F(61,100), C0-(1-kap**6)*F(9,10)**7, F(21,20)-A0, B0-F(7,6)])
    # The current article uses a sharper row-32 line than the old auxiliary slack.
    a, p = F(5,4), F(373,500)
    check('current row32 line exists', (a,p) in byid[32]['lines'])
    eq('current row32 sharp slack', 7-byid[32]['cap']+u0*(a*byid[32]['cap']+p-8), F(1599,6250))
    check('both e6 deficit gates exceed 39/20', all(((8-early_lines[rid][1])*u0-7)/(early_lines[rid][0]*u0-1)>F(39,20) for rid in (32,33)))
    details['states'] = len(rows)
    details['current_row32_slack'] = '1599/6250'
    return {'checks': checks, 'details': details}
