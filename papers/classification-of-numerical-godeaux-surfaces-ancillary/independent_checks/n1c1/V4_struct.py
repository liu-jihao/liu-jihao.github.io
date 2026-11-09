#!/usr/bin/env python3
"""Structural checks on the equations produced by V1_<variant>.sing (light Python, no CAS):
 (a) w-order structure: in every equation of w-order m, every term has total unknown-w-order <= m
     (this implies: unknowns of order m occur linearly, with coefficients in order-0 unknowns only);
 (b) torus: every equation is weighted-homogeneous for the grading direction (all coefficient weights 0)
     and for the one-parameter group u (w -> t^-1 w, y3 -> t^-2 y3, P = diag(t^-2,t^-1,t,t,t^-1,t^2), t^2 = u);
     weights derived from A' = P^t A(phi) P and lambda'_{ij,k} = t^(P0+Pk-Pi-Pj) lambda_{ij,k}(phi);
 (c) weights of the free A-coefficients found by V1 (effective torus dimension e).
Usage: python3 V4_struct.py slice|general  > out_V4_<variant>.txt"""
import re, sys, os
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
variant = sys.argv[1] if len(sys.argv) > 1 else 'slice'
cd = [0, 3, 3, 3, 3, 4]

# unknowns table written by gen_v.py
info = {}
for line in open(os.path.join(HERE, 'V0_unknowns_%s.txt' % variant)):
    if line.startswith('#'):
        continue
    fl = line.split()
    name, typ, mono, word = fl[0], fl[1], fl[-2], fl[-1]
    info[name] = (typ, mono, int(word))


def parse_wy(mono):
    """'w^2*y1*y3' -> (a,b,e,f)"""
    e = {'w': 0, 'y1': 0, 'y2': 0, 'y3': 0}
    if mono != '1':
        for fac in mono.split('*'):
            if '^' in fac:
                v, p = fac.split('^')
                e[v] += int(p)
            else:
                e[fac] += 1
    return (e['w'], e['y1'], e['y2'], e['y3'])


def parse_Y(name):
    """'y0y1y3' -> exponents (a0,a1,a2,a3)"""
    ex = [0, 0, 0, 0]
    for k in re.findall(r'y(\d)', name):
        ex[int(k)] += 1
    return tuple(ex)


def torus_weights(W, Y1, Y2, Y3, P):
    def wt_wy(a, b, e, f):
        return W * a + Y1 * b + Y2 * e + Y3 * f
    wts = {}
    for name, (typ, mono, word) in info.items():
        if typ.startswith('lambda'):
            i, j, k = int(typ[8]), int(typ[9]), int(typ[11])
            wts[name] = P[0] + P[k] - P[i] - P[j] + wt_wy(*parse_wy(mono))
        elif typ == 'G':
            a0, a1, a2, a3 = parse_Y(mono)
            wts[name] = 2 * P[0] + 2 * W + wt_wy(2 * a0, a1, a2, a3) + 2 * W * 0
        elif typ.startswith('q'):
            kq = int(typ[1])
            a0, a1, a2, a3 = parse_Y(mono)
            wts[name] = P[0] + P[kq] + W + wt_wy(2 * a0, a1, a2, a3)
        elif typ.startswith('h('):
            a0, a1, a2, a3 = parse_Y(mono)
            wts[name] = 2 * P[0] + wt_wy(2 * a0, a1, a2, a3)
        elif typ.startswith('L('):
            a0, a1, a2, a3 = parse_Y(mono)
            wts[name] = P[0] + P[5] + 2 * W + wt_wy(2 * a0, a1, a2, a3)
        elif typ == 'D':
            wts[name] = 2 * P[1] + Y3
    return wts


# fixed entries of the normal form must have weight 0 (checked here as well)
def fixed_entry_check(W, Y1, Y2, Y3, P):
    wt = lambda a, b, e, f: W * a + Y1 * b + Y2 * e + Y3 * f
    checks = {
        '(1,2) y2': P[1] + P[2] + wt(0, 0, 1, 0), '(1,3) y1': P[1] + P[3] + wt(0, 1, 0, 0),
        '(2,2) y3': 2 * P[2] + wt(0, 0, 0, 1), '(2,3) y0': P[2] + P[3] + wt(2, 0, 0, 0),
        '(2,4) y1': P[2] + P[4] + wt(0, 1, 0, 0), '(3,3) y3': 2 * P[3] + wt(0, 0, 0, 1),
        '(3,4) y2': P[3] + P[4] + wt(0, 0, 1, 0), '(1,5) w': P[1] + P[5] + wt(1, 0, 0, 0),
        '(0,5) y2^2': P[0] + P[5] + wt(0, 0, 2, 0), '(0,5) y1^2': P[0] + P[5] + wt(0, 2, 0, 0),
    }
    dw = {'D at (1,1)': 2 * P[1] + wt(0, 0, 0, 1), 'D at (4,4)': 2 * P[4] + wt(0, 0, 0, 1),
          'D at (0,5)': P[0] + P[5] + wt(0, 0, 0, 2)}
    return checks, dw


tok = re.compile(r'([+-]?)([^+-]+)')


def parse_poly(s):
    """Singular polynomial string -> list of (coef Fraction, {var: exp})"""
    s = s.strip()
    terms = []
    depth = 0
    cur = ''
    pieces = []
    for ch in s:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        if ch in '+-' and depth == 0 and cur.strip() != '':
            pieces.append(cur)
            cur = ch
        else:
            cur += ch
    if cur.strip():
        pieces.append(cur)
    for p in pieces:
        p = p.strip()
        sign = 1
        if p[0] in '+-':
            sign = -1 if p[0] == '-' else 1
            p = p[1:]
        coef = Fraction(sign)
        mon = {}
        for fac in p.split('*'):
            if re.fullmatch(r'\d+(/\d+)?', fac):
                coef *= Fraction(fac)
            else:
                if '^' in fac:
                    v, e = fac.rsplit('^', 1)
                    e = int(e)
                else:
                    v, e = fac, 1
                mon[v] = mon.get(v, 0) + e
        terms.append((coef, mon))
    return terms


def split_ideal(s):
    out, depth, cur = [], 0, ''
    for ch in s:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        if ch == ',' and depth == 0:
            out.append(cur)
            cur = ''
        else:
            cur += ch
    if cur.strip():
        out.append(cur)
    return out


eqs = split_ideal(open(os.path.join(HERE, 'V1_%s_eqs_ideal.txt' % variant)).read())
ords = [int(x) for x in open(os.path.join(HERE, 'V1_%s_eqs_orders.txt' % variant)).read().strip().split(',')]
assert len(eqs) == len(ords), (len(eqs), len(ords))
print('variant', variant, ': equations parsed:', len(eqs))

# (a) order structure
viol = 0
eq_cnt_terms = 0
for s, m in zip(eqs, ords):
    for coef, mon in parse_poly(s):
        eq_cnt_terms += 1
        tot = sum(info[v][2] * e for v, e in mon.items())
        if tot > m:
            viol += 1
print('(a) terms checked:', eq_cnt_terms, '; terms whose total unknown-w-order exceeds the w-order of their equation:', viol)

# (b) torus homogeneity
tori = {
    'grading (W,Y1,Y2,Y3;P)=(1,2,2,2;-4,-1,-1,-1,-1,0)': (1, 2, 2, 2, [-4, -1, -1, -1, -1, 0]),
    'u (t-units) (W,Y1,Y2,Y3;P)=(-1,0,0,-2;-2,-1,1,1,-1,2)': (-1, 0, 0, -2, [-2, -1, 1, 1, -1, 2]),
}
for label, (W, Y1, Y2, Y3, P) in tori.items():
    checks, dw = fixed_entry_check(W, Y1, Y2, Y3, P)
    wts = torus_weights(W, Y1, Y2, Y3, P)
    nonhom = 0
    for s in eqs:
        ws = set(sum(wts[v] * e for v, e in mon.items()) for coef, mon in parse_poly(s))
        if len(ws) > 1:
            nonhom += 1
    print('(b) torus', label)
    print('    fixed entries weight (must all be 0):', checks)
    print('    weight of D from its three entries (must agree):', dw)
    print('    equations not weighted-homogeneous:', nonhom)
    if label.startswith('grading'):
        print('    coefficients with nonzero grading weight:', [v for v in wts if wts[v] != 0])
    else:
        Aw = {v: Fraction(wts[v], 2) for v in wts if not info[v][0].startswith('lambda')}
        print('    u-weights (u-units = t-weight/2) of A-coefficients:', {v: str(Aw[v]) for v in sorted(Aw, key=lambda z: (z[0], int(re.sub(r'\D', '', z) or 0)))})
        # free A-coefficients from V1
        fr = open(os.path.join(HERE, 'V1_%s_free.txt' % variant)).read().split('\n')
        freeA = fr[0].split(':', 1)[1].split()
        freeX = fr[1].split(':', 1)[1].split()
        print('    free A-coefficients and u-weights:', {v: str(Fraction(wts[v], 2)) for v in freeA})
        print('    free A-coefficients with zero u-weight:', [v for v in freeA if wts[v] == 0])
        print('    free lambda with nonzero u-weight:', sum(1 for v in freeX if wts[v] != 0), 'of', len(freeX))
