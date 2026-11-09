#!/usr/bin/env python3
"""Light bookkeeping (no CAS): per-stage summary of out_n3_elim.txt.
For every eliminated unknown: its w-order, its 'visibility stage' (order + delay; delay(k) for lambda_{ij,k} from
out_n3_aux.txt: k=5: 0, k=2,3: 1, k=1,4: 3, k=0: 4; A-coefficients: delay 0) and the stage at which it was eliminated."""
import re, os
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
tab = {}
for line in open(os.path.join(HERE, 'n3_unknowns.txt')):
    if line.startswith('#'): continue
    p = line.split()
    idx, name, typ = int(p[0]), p[1], p[2]
    order = int(p[-1])
    k = None
    m = re.match(r'lambda_\{(\d)(\d),(\d)\}', typ)
    if m: k = int(m.group(3))
    tab[name] = (idx, typ, order, k)
delay = {5: 0, 2: 1, 3: 1, 1: 3, 4: 3, 0: 4}
out = open(os.path.join(HERE, 'out_n3_elim.txt')).read()
stages = re.findall(r'eliminated at stage (\d+) \(name\[w-order\]\): ([^\n]*)', out)
late = []
summary = []
for st, lst in stages:
    st = int(st)
    names = re.findall(r'((?:lam|g|b)\(\d+\))\[(\d+)\]', lst)
    c = Counter()
    for nm, o in names:
        idx, typ, order, k = tab[nm]
        assert order == int(o)
        if k is not None:
            vis = order + delay[k]
            c[('lambda k=%d' % k, order)] += 1
            if vis != st: late.append((nm, typ, order, vis, st))
        else:
            c[(nm, order)] += 1
            if order != st: late.append((nm, typ, order, order, st))
    summary.append((st, len(names), sorted(c.items())))
for st, n, items in summary:
    print('stage %2d: %3d eliminated:' % (st, n), '; '.join('%s ord %d: %d' % (a, o, cnt) if a.startswith('lambda') else '%s (ord %d)' % (a, o) for (a, o), cnt in items))
print('\nunknowns eliminated at a stage different from their visibility stage (order + delay):')
for t in late: print('  %-9s %-18s order %d, visible at %d, eliminated at stage %d' % t)
