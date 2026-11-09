#!/usr/bin/env python3
"""Comparison with the case's outputs (run only after V1/V2/V4 were final).
 (a) the case's 782 equations (cases/n1c1/eqs_n1c1.txt), translated to my unknowns through the monomial
     dictionary (lambda_{ij,k} monomial, G monomial, q_i monomial, D), versus my 782 equations
     (V1_slice_eqs_ideal.txt): equal as sets up to nonzero scalars, with the same w-orders?
 (b) writes V3_graph.sing: the two graphs Gamma (case: graph_expr_n1c1.txt; mine: V1_slice_param_ideal.txt)
     are compared exactly in both directions by substituting one parametrization into the other's relations.
Usage: python3 V3_compare.py > out_V3_compare.txt"""
import re, os
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(HERE, '..', '..', 'cases', 'n1c1')


def wy_exp(mono):
    e = {'w': 0, 'y1': 0, 'y2': 0, 'y3': 0, 'y0': 0}
    if mono.strip() != '1':
        for fac in mono.split('*'):
            fac = fac.strip()
            if '^' in fac:
                v, p = fac.split('^')
                e[v] += int(p)
            else:
                e[fac] += 1
    return e


# my unknowns
mine = {}
for line in open(os.path.join(HERE, 'V0_unknowns_slice.txt')):
    if line.startswith('#'):
        continue
    fl = line.split()
    name, typ, mono = fl[0], fl[1], fl[-2]
    if typ.startswith('lambda'):
        e = wy_exp(mono)
        key = ('l', typ[8], typ[9], typ[11], e['w'], e['y1'], e['y2'], e['y3'])
    elif typ == 'G':
        ex = [0, 0, 0, 0]
        for k in re.findall(r'y(\d)', mono):
            ex[int(k)] += 1
        key = ('g',) + tuple(ex)
    elif typ.startswith('q'):
        ex = [0, 0, 0, 0]
        for k in re.findall(r'y(\d)', mono):
            ex[int(k)] += 1
        key = ('q', typ[1]) + tuple(ex)
    else:
        key = ('D',)
    mine[key] = name

# case unknowns
dic = {}
for line in open(os.path.join(CASE, 'unknowns_n1c1.txt')):
    if line.startswith('#') or '|' not in line:
        continue
    fl = [x.strip() for x in line.split('|')]
    name, typ, mono = fl[0], fl[1], fl[2]
    if typ.startswith('lambda'):
        m = re.match(r'lambda_\{(\d)(\d),(\d)\}', typ)
        e = wy_exp(mono)
        key = ('l', m.group(1), m.group(2), m.group(3), e['w'], e['y1'], e['y2'], e['y3'])
    elif typ.startswith('G'):
        e = wy_exp(mono)
        key = ('g', e['y0'], e['y1'], e['y2'], e['y3'])
    elif typ.startswith('q'):
        e = wy_exp(mono)
        key = ('q', typ[1], e['y0'], e['y1'], e['y2'], e['y3'])
    elif name == 'D':
        key = ('D',)
    else:
        continue
    dic[name] = mine[key]
print('dictionary entries (case -> mine):', len(dic), '; injective:', len(set(dic.values())) == len(dic))
print('A-coefficients:', {k: v for k, v in dic.items() if not k.startswith('l(')})


def parse_poly(s):
    s = s.strip()
    depth, cur, pieces = 0, '', []
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
    terms = {}
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
                if '^' in fac and not fac.endswith(')'):
                    v, e = fac.rsplit('^', 1)
                    e = int(e)
                else:
                    v, e = fac, 1
                mon[v] = mon.get(v, 0) + e
        key = tuple(sorted(mon.items()))
        terms[key] = terms.get(key, 0) + coef
    return {k: v for k, v in terms.items() if v != 0}


def normalize(t):
    k0 = min(t.keys())
    c = t[k0]
    return frozenset((k, v / c) for k, v in t.items())


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


my_eqs = split_ideal(open(os.path.join(HERE, 'V1_slice_eqs_ideal.txt')).read())
my_ord = [int(x) for x in open(os.path.join(HERE, 'V1_slice_eqs_orders.txt')).read().strip().split(',')]
myset = {}
for s, o in zip(my_eqs, my_ord):
    myset[normalize(parse_poly(s))] = o

case_set = {}
for line in open(os.path.join(CASE, 'eqs_n1c1.txt')):
    m = re.match(r'EQ\[\d+\]=(.*); EO\[\d+\]=(\d+);', line.strip())
    if not m:
        continue
    t = parse_poly(m.group(1))
    tt = {}
    for key, c in t.items():
        nk = tuple(sorted((dic[v], e) for v, e in key))
        tt[nk] = tt.get(nk, 0) + c
    case_set[normalize(tt)] = int(m.group(2))
print('(a) my distinct equations:', len(myset), '; case distinct equations (translated):', len(case_set))
common = set(myset) & set(case_set)
print('    common up to scalar:', len(common), '; only mine:', len(set(myset) - set(case_set)), '; only case:', len(set(case_set) - set(myset)))
print('    w-order mismatches among common equations:', sum(1 for k in common if myset[k] != case_set[k]))

# (b) Singular script for the graph comparison
rels = []
for line in open(os.path.join(CASE, 'graph_expr_n1c1.txt')):
    m = re.match(r'([gb]\(\d+\)) = (.*);', line.strip())
    if m:
        rels.append((m.group(1), m.group(2)))
case_free = ['D', 'g(5)', 'b(2)', 'b(5)', 'b(6)', 'b(8)', 'b(9)', 'b(11)', 'b(12)']
caseA = ['g(%d)' % i for i in range(1, 11)] + ['b(%d)' % i for i in range(1, 13)] + ['D']
elimc = dict(rels)
assert set(elimc) | set(case_free) == set(caseA) and not (set(elimc) & set(case_free))
myvars = ['x(%d)' % i for i in range(1, 372)] + ['qq(%d)' % i for i in range(1, 13)] + ['gg(%d)' % i for i in range(1, 11)] + ['D', 'y3', 'y2', 'y1', 'w', 'tt']
idx = {v: n + 1 for n, v in enumerate(myvars)}
S = []
S.append('// V3_graph.sing -- generated by V3_compare.py: exact comparison of the case graph and my graph (both directions).')
S.append('ring R0 = 0, (x(1..371), qq(1..12), gg(1..10), D, y3, y2, y1, w, tt), dp;')
S.append('execute("ideal MM = " + read("V1_slice_param_ideal.txt") + ";");')
S.append('ring RC = 0, (g(1..10), b(1..12), D), dp;')
S.append('ideal REL = ' + ', '.join('%s - (%s)' % (l, r) for l, r in rels) + ';')
S.append('ideal PAR = ' + ', '.join((elimc[v] if v in elimc else v) for v in caseA) + ';   // case parametrization of g(1..10), b(1..12), D')
S.append('setring R0;')
imgs = []
for v in caseA:
    imgs.append('MM[%d]' % idx[dic[v]])
S.append('map phi = RC, ' + ', '.join(imgs) + ';')
S.append('ideal CHK1 = phi(REL);')
S.append('int i, bad; bad = 0; for (i=1;i<=ncols(CHK1);i++) { if (CHK1[i] != 0) { bad++; } }')
S.append('"(b1) case relations evaluated on my parametrization of Gamma: nonzero =", bad, " of", ncols(CHK1);')
S.append('ideal MYREL; int nr = 0; for (i=372;i<=394;i++) { if (MM[i] != var(i)) { nr++; MYREL[nr] = var(i) - MM[i]; } }')
S.append('"my relations (eliminated A-coefficients):", nr;')
S.append('setring RC;')
inv = {dic[v]: v for v in caseA}
imgs2 = []
for v in myvars:
    if v in inv:
        imgs2.append('PAR[%d]' % (caseA.index(inv[v]) + 1))
    else:
        imgs2.append('0')
S.append('map psi = R0, ' + ', '.join(imgs2) + ';')
S.append('ideal CHK2 = psi(MYREL);')
S.append('bad = 0; for (i=1;i<=ncols(CHK2);i++) { if (CHK2[i] != 0) { bad++; } }')
S.append('"(b2) my relations evaluated on the case parametrization of Gamma: nonzero =", bad, " of", ncols(CHK2);')
S.append('quit;')
open(os.path.join(HERE, 'V3_graph.sing'), 'w').write('\n'.join(S) + '\n')
print('(b) wrote V3_graph.sing')
