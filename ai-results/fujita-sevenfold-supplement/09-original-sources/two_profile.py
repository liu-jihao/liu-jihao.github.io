"""Exact written-cell checker adapted from sub11_all27_written_audit.py.
Only Q/M/D data are accepted by check.py; no producer or private inputs.
"""
from fractions import Fraction as F
from collections import Counter
import re

R, U, Z, END = F(9, 10), F(26, 25), F(27, 25), F(13, 5)
H0 = F(131, 1000)
QMIN, QMAX = F(16020581569, 90194313216), F(8010290789, 45097156608)
DEN = 10**9
CHECKS = Counter()
CONTEXT = ''


def check(ok, name):
    CHECKS[name] += 1
    if not ok:
        raise AssertionError(f'{CONTEXT}: {name}')


def integrator():
    nodes = [F(j, 6) for j in range(7)]
    weights = []
    for j, x in enumerate(nodes):
        poly = [F(1)]
        for m, y in enumerate(nodes):
            if m == j:
                continue
            nxt = [F(0)] * (len(poly) + 1)
            for i, coefficient in enumerate(poly):
                nxt[i] -= y * coefficient / (x - y)
                nxt[i + 1] += coefficient / (x - y)
            poly = nxt
        weights.append(sum(v / (i + 1) for i, v in enumerate(poly)))
    for power in range(7):
        check(sum(w*x**power for w, x in zip(weights, nodes)) == F(1, power+1),
              'exact_interpolation_moment')
    return nodes, weights


NODES, WEIGHTS = integrator()


def integrate(fn, left, right):
    return (right-left)*sum(w*fn(left+x*(right-left)) for x,w in zip(NODES, WEIGHTS))


def coverage(family, h1, anchor, max1, expected):
    global CONTEXT
    grid0, grid1 = family['grid0'], family['grid1']
    check(grid0[0] == H0/R and grid0[-1] == 16, 'entire_old_slope_domain')
    check(grid1[0] == h1/anchor and grid1[-1] == max1, 'entire_new_slope_domain')
    check(all(a < b for grid in (grid0, grid1) for a,b in zip(grid, grid[1:])), 'increasing_grids')
    active = {f'{i+1}.{j+1}': [a,b,c,d]
              for i,(a,b) in enumerate(zip(grid0,grid0[1:]))
              for j,(c,d) in enumerate(zip(grid1,grid1[1:]))}
    for tag, coordinate, midpoint, line_no in family['splits']:
        CONTEXT = f'split {tag} line {line_no}'
        check(tag in active and coordinate in (0,1), 'split_active_leaf')
        bounds = active.pop(tag)
        i = 2*coordinate
        check(midpoint == (bounds[i]+bounds[i+1])/2, 'printed_exact_bisection')
        low, high = bounds.copy(), bounds.copy()
        low[i+1] = midpoint
        high[i] = midpoint
        check(tag+'0' not in active and tag+'1' not in active, 'new_children')
        active[tag+'0'], active[tag+'1'] = low, high
    check(set(active) == set(family['rectangles']) == set(family['cells']), 'exact_terminal_tag_set')
    for tag, bounds in active.items():
        CONTEXT = f'terminal {tag}'
        check(bounds == family['rectangles'][tag][0], 'terminal_bounds_match_split_tree')
    check(len(active) == expected, 'expected_terminal_count')


def baseline(t):
    return F(1) if t < R else F(61,100) if t < F(99,100) else F(11,25) if t < U else F(17,50)


def cell_function(name, left, right, bounds, kind, h1, anchor, tail=None):
    check(left < right, 'positive_cell_length')
    if name.startswith('B'):
        if name == 'Bhigh':
            check(kind in ('P2','T') and F(29,25) <= left < right <= F(125,52), 'high_baseline_domain')
            return lambda t: 7*((125-52*t)/73)**6
        if name == 'Bzero':
            check(kind in ('P2','T') and left >= F(125,52), 'zero_baseline_domain')
            return lambda t: F(0)
        rate = F(name[1:])
        splits = sorted({left,right} | {v for v in (R,F(99,100),U) if left < v < right})
        check(all(rate >= baseline((a+b)/2) for a,b in zip(splits,splits[1:])), 'baseline_rate_dominates')
        return lambda t: 7*rate*t**6
    profile, label = name.split(':')
    index = int(profile)
    check(index in (0,1), 'profile_index')
    v, h = (R,H0) if index == 0 else (anchor,h1)
    if tail is not None:
        check(profile == tail[-1] and left >= v, 'tail_uses_only_observed_profile_after_anchor')
        a = F(16)
    else:
        check(not left < v < right, 'no_changed_profile_slope_inside_cell')
        a = bounds[2*index+1] if right <= v else bounds[2*index]
    d = a*v-h
    check(d >= 0, 'nonnegative_profile_offset')
    check(a*left-d >= 0 and a*right-d > 0, 'positive_unclipped_profile_on_open_cell')
    c_left = a-d/left if left else a
    c_right = a-d/right
    cmin, cmax = min(c_left,c_right), max(c_left,c_right)
    def domain(low=F(0), high=None):
        check(cmin >= low and (high is None or cmax <= high), 'legal_capacity_branch')
    needed = []
    if label.startswith('N'):
        m = re.fullmatch(r'N([35])(lo|hi)', label)
        check(m is not None, 'generic_capacity_name')
        q, branch = int(m[1]), m[2]
        check((index == 0 and q == 5) or (index == 1 and kind == 'M' and q == 3), 'generic_codimension_matches')
        domain(F(0), F(1,7)) if branch == 'lo' else domain(F(1,7), F(1,2))
        if branch == 'lo':
            return lambda t: 7*t**6-7*F(7*q,q+1)**q*t**(6-q)*(a*t-d)**q
        return lambda t: 7*t**6-7/F(6*(q+1))**q*t**(6-q)*((6*q-1+7*a)*t-7*d)**q
    if label == 'C':
        check(index == 0 or kind == 'M', 'convex_order_capacity_matches')
        domain(F(1,7),F(1,2))
        needed = [1]
    elif label == 'E1':
        if index == 0 or kind == 'M':
            domain(F(1,2))
        else:
            check(kind == 'Q', 'endpoint_capacity_family')
            domain(QMAX)
        needed = [1]
    elif label == 'E2':
        check(index == 1 and kind in ('P2','T'), 'double_capacity_family')
        domain(F(4,25) if kind == 'T' else F(0))
        needed = [2]
    elif label == 'Tlow':
        check(index == 1 and kind == 'T', 'triple_spanning_family')
        domain(F(0),F(4,25))
        needed = [2,3]
    elif label == 'Qbal':
        check(index == 1 and kind == 'Q', 'balanced_capacity_family')
        domain(F(0), QMIN)
        needed = [2]
    elif label in ('Mspan','Mref'):
        check(index == 1 and kind == 'M', 'common_linear_capacity_family')
        domain(F(2,11) if label == 'Mref' else F(0))
        needed = [1,2,3,4] if label == 'Mspan' else [1,2,3]
    elif label in ('H3','H3lip'):
        check(index == 1 and kind == 'M', 'height_three_family')
        domain(F(1,6)) if label == 'H3' else domain(F(0),F(1,6))
        return (lambda t: 7*F(256,729)*t**6) if label == 'H3' else (
            lambda t: 7*(F(256,729)+3-18*a)*t**6+126*d*t**5)
    elif label in ('Dlo','Dhi'):
        check(index == 1 and kind == 'D', 'full_embedding_seven_height_three_family')
        domain(F(1,6) if label == 'Dhi' else F(0))
        if label == 'Dhi':
            return lambda t: 7*F(138259,900000)*t**6
        needed = [2,3,4]
    else:
        raise AssertionError(f'{CONTEXT}: unknown polynomial {name}')
    signs = {}
    for j in needed:
        pleft, pright = (1-j*a)*left+j*d, (1-j*a)*right+j*d
        check(not min(pleft,pright) < 0 < max(pleft,pright), 'positive_part_has_no_interior_crossing')
        signs[j] = pleft+pright > 0
    def p(j,t):
        return 7*((1-j*a)*t+j*d)**6 if signs[j] else F(0)
    if label.startswith('E'):
        return lambda t: p(int(label[1:]),t)
    if label == 'C':
        return lambda t: F(7,6)**6*p(1,t)
    if label == 'Tlow':
        return lambda t: 3*p(2,t)-2*p(3,t)
    if label == 'Qbal':
        check(signs[2], 'balanced_positive_part_is_positive')
        return lambda t: 7*((1-2*a)*t+2*d)**5*((1+10*a)*t-10*d)
    if label == 'Mspan':
        return lambda t: p(1,t)+3*p(2,t)-5*p(3,t)+2*p(4,t)
    if label == 'Mref':
        return lambda t: p(1,t)+p(2,t)-p(3,t)
    return lambda t: 6*p(2,t)-8*p(3,t)+3*p(4,t)


def audit_cells(cells,bounds,kind,h1,anchor,endpoint,printed_total,tail=None):
    global CONTEXT
    check(cells[0][0] == 0 and cells[-1][1] == endpoint, 'full_order_domain')
    exact_total = F(0)
    details = []
    for i,(left,right,name,numerator,line_no) in enumerate(cells):
        CONTEXT = f'{kind} {tail or bounds} proof line {line_no} {name} [{left},{right}]'
        if i:
            check(cells[i-1][1] == left, 'adjacent_cells_no_gap_or_overlap')
        fn = cell_function(name,left,right,bounds,kind,h1,anchor,tail)
        exact = integrate(fn,left,right)
        check(exact >= 0, 'nonnegative_exact_integral')
        check(exact < F(numerator,DEN), 'strict_printed_cell_integral')
        check(exact.numerator*DEN//exact.denominator+1 == numerator, 'printed_integer_is_floor_plus_one')
        exact_total += exact
        details.append({'line':line_no,'interval':[str(left),str(right)],'polynomial':name,
                        'exact_integral':str(exact),'printed_upper':numerator})
    check(sum(c[3] for c in cells) == printed_total, 'sum_of_cell_bounds_matches_total')
    check(exact_total < F(printed_total,DEN), 'exact_total_below_printed_total')
    return {'exact_total':str(exact_total),'upper_numerator':printed_total,'cells':details}

