#!/usr/bin/env python3
"""Independent generator for the verification of case n1c0 (N1(D), c = 0).
Written from plan/case_common.md, plan/case_n1c1.md, plan/case_n1c0.md and CCO94 7.3-7.5 (page image p. 123),
without using the case's scripts.  Own conventions:
  ring variables: w, y1, y2, y3, D, b(1..12), g(1..10), L(1..371), t   (t only used for localizing at b8)
  lambda coefficients L(n) ordered by k = 5,4,3,2,1,0, then pair (ij) = 11,12,...,55, then monomial with
  w-exponent descending (then exponents of y1, y2, y3 descending).
Writes: vring.sing (ring only), vbuild.sing (ring + system + equations), v_unknowns.txt (table), v_meta.txt.
"""
import itertools, os
HERE = os.path.dirname(os.path.abspath(__file__))
c = [0, 3, 3, 3, 3, 4]
dl = [1, -1, -1, 1, 1, -1]
pairs = [(i, j) for i in range(1, 6) for j in range(i, 6)]

def monos(weight, sign):
    """monomials w^a y1^b y2^e y3^f of weight a + 2(b+e+f) and sigma-sign (-1)^(a+e+f) = sign"""
    res = []
    if weight < 0:
        return res
    for a in range(weight, -1, -1):
        r = weight - a
        if r % 2:
            continue
        n = r // 2
        for b in range(n, -1, -1):
            for e in range(n - b, -1, -1):
                f = n - b - e
                if (-1) ** (a + e + f) == sign:
                    res.append((a, b, e, f))
    return res

def mstr(m):
    a, b, e, f = m
    parts = []
    for v, x in (('w', a), ('y1', b), ('y2', e), ('y3', f)):
        if x == 1:
            parts.append(v)
        elif x > 1:
            parts.append('%s^%d' % (v, x))
    return '*'.join(parts) if parts else '1'

# lambda ansatz
lam = []          # list of (n, i, j, k, mono)
bypk = {}         # (pair index, k) -> list of n
n = 0
for k in [5, 4, 3, 2, 1, 0]:
    for p, (i, j) in enumerate(pairs):
        for m in monos(c[i] + c[j] - c[k], dl[i] * dl[j] * dl[k]):
            n += 1
            lam.append((n, i, j, k, m))
            bypk.setdefault((p, k), []).append(n)
NL = n
cnt_k = {k: sum(1 for x in lam if x[3] == k) for k in range(6)}
cnt_a = {}
for x in lam:
    cnt_a[x[4][0]] = cnt_a.get(x[4][0], 0) + 1

# equation slots
slots = {}
for (i, j) in pairs:
    for m in monos(8 + c[i] + c[j], -dl[i] * dl[j]):
        slots[m[0]] = slots.get(m[0], 0) + 1

# sigma-fixed point data
F2 = []   # pure w-power coefficient of lambda_{ij,0}
F3 = {}   # (i,j) -> list of (n, (e, f)) pure (y2,y3) coefficients of lambda_{ij,0}
for (n_, i, j, k, m) in lam:
    if k != 0:
        continue
    a, b, e, f = m
    if b == 0 and e == 0 and f == 0:
        F2.append((n_, i, j, a))
    if a == 0 and b == 0:
        F3.setdefault((i, j), []).append((n_, (e, f)))

# F3 remainders modulo y2^2 - D*y3^2 (sufficient condition for vanishing at (0:0:+-sqrt(D):1|0))
def rem_poly(coeffs, deg):
    """coeffs: dict e -> string coefficient of y2^e y3^(deg-e). Reduce y2^2 -> D y3^2.
    Returns (r1, r0) as strings: remainder r1*y2*y3^(deg-1) + r0*y3^deg."""
    r1, r0 = [], []
    for e, s in coeffs.items():
        if e % 2 == 1:
            r1.append('D^%d*(%s)' % ((e - 1) // 2, s))
        else:
            r0.append('D^%d*(%s)' % (e // 2, s))
    return ('+'.join(r1) if r1 else '0', '+'.join(r0) if r0 else '0')

F3rem = []
for (i, j), lst in sorted(F3.items()):
    deg = lst[0][1][0] + lst[0][1][1]
    co = {e: 'L(%d)' % n_ for (n_, (e, f)) in lst}
    r1, r0 = rem_poly(co, deg)
    F3rem.append(((i, j), r1, r0))

# ---------------- Singular code ----------------
ring = 'ring R = 0, (w, y1, y2, y3, D, b(1..12), g(1..10), L(1..%d), t), dp;\n' % NL
with open(os.path.join(HERE, 'vring.sing'), 'w') as fh:
    fh.write('// ring of the n1c0 verification (vgen.py)\n' + ring + 'option(redSB);\n')

S = []
S.append('// vbuild.sing: written by vgen.py. Builds A (N1(D), c = 0), the cofactors B, the rank conditions E_ij and their coefficient equations.')
S.append(ring)
S.append('option(redSB);')
S.append('int i; int j; int k; int p; int m; int ne;')
S.append('poly y0 = w^2;')
S.append('poly G = g(1)*y3^3 + g(2)*y1^2*y3 + g(3)*y2*y3^2 + g(4)*y1^2*y2 + g(5)*y2^2*y3 + g(6)*y2^3 + g(7)*y0*y1*y3 + g(8)*y0*y1*y2 + g(9)*y0^2*y3 + g(10)*y0^2*y2;')
S.append('poly q1 = b(1)*y1*y2 + b(2)*y1*y3;')
S.append('poly q2 = b(3)*y0*y2 + b(4)*y0*y3 + b(5)*y1*y2 + b(6)*y1*y3;')
S.append('poly q3 = b(7)*y0^2 + b(8)*y0*y1 + b(9)*y1^2;')
S.append('poly q4 = b(10)*y0^2 + b(11)*y2*y3 + b(12)*y3^2;')
S.append('poly Q = y2^2 - D*y3^2 - y1^2;')
S.append('matrix A[6][6];')
S.append('A[1,1] = y0*G; A[1,2] = w*q1; A[1,3] = w*q2; A[1,4] = w*q3; A[1,5] = w*q4; A[1,6] = Q;')
S.append('A[2,2] = D*y3; A[2,3] = y2; A[2,4] = y1; A[2,5] = 0;     A[2,6] = w;')
S.append('A[3,3] = y3;   A[3,4] = 0;  A[3,5] = y1;  A[3,6] = 0;')
S.append('A[4,4] = -y3;  A[4,5] = y2; A[4,6] = 0;')
S.append('A[5,5] = -D*y3; A[5,6] = 0;')
S.append('A[6,6] = 0;')
S.append('for (i = 1; i <= 6; i++) { for (j = i + 1; j <= 6; j++) { A[j,i] = A[i,j]; } }')
S.append('intvec cc = 0, 3, 3, 3, 3, 4; intvec dd = 1, -1, -1, 1, 1, -1;')
S.append('''proc mysub(matrix M, int r, int s)
{
  int nn = nrows(M); matrix SS[nn-1][nn-1]; int a; int bb; int ii; int jj;
  ii = 0;
  for (a = 1; a <= nn; a++) { if (a != r) { ii++; jj = 0; for (bb = 1; bb <= nn; bb++) { if (bb != s) { jj++; SS[ii,jj] = M[a,bb]; } } } }
  return(SS);
}''')
S.append('// B0[1,k] = B_{0,k-1} = adj(A)_{0,k-1}; A symmetric, so the minor deleting row 1, column k is used')
S.append('matrix B0[1][6];')
S.append('for (k = 1; k <= 6; k++) { B0[1,k] = (-1)^(1+k) * det(mysub(A, 1, k)); }')
S.append('matrix BP[1][15];')
for p, (i, j) in enumerate(pairs):
    S.append('BP[1,%d] = (-1)^(%d) * det(mysub(A, %d, %d));' % (p + 1, (i + 1) + (j + 1), i + 1, j + 1))
S.append('matrix LAM[15][6];')
for p, (i, j) in enumerate(pairs):
    for k in range(6):
        terms = ['L(%d)*%s' % (n_, mstr(m)) for (n_, ii, jj, kk, m) in lam if ii == i and jj == j and kk == k]
        if terms:
            S.append('LAM[%d,%d] = %s;' % (p + 1, k + 1, ' + '.join(terms)))
S.append('ideal EE;')
S.append('for (p = 1; p <= 15; p++) { EE[p] = BP[1,p]; for (k = 1; k <= 6; k++) { EE[p] = EE[p] - LAM[p,k]*B0[1,k]; } }')
S.append('intvec PI = %s; intvec PJ = %s;' % (','.join(str(i) for i, j in pairs), ','.join(str(j) for i, j in pairs)))
S.append('''// coefficient equations, with w-order and pair; weight/sign check of every monomial
ideal EQ; intvec ORD; intvec PR; ne = 0; int badws = 0; intvec ex; int wt; int sg;
for (p = 1; p <= 15; p++)
{
  def CC = coef(EE[p], w*y1*y2*y3);
  for (k = 1; k <= ncols(CC); k++)
  {
    if (CC[2,k] != 0)
    {
      ne++; EQ[ne] = CC[2,k]; ex = leadexp(CC[1,k]); ORD[ne] = ex[1]; PR[ne] = p;
      wt = ex[1] + 2*(ex[2] + ex[3] + ex[4]);
      sg = (-1)^(ex[1] + ex[3] + ex[4]);
      if (wt != 8 + cc[PI[p]+1] + cc[PJ[p]+1]) { badws++; }
      if (sg != -dd[PI[p]+1]*dd[PJ[p]+1]) { badws++; }
    }
  }
  kill CC;
}''')
with open(os.path.join(HERE, 'vbuild.sing'), 'w') as fh:
    fh.write('\n'.join(S) + '\n')

# meta and table
with open(os.path.join(HERE, 'v_meta.txt'), 'w') as fh:
    fh.write('lambda count %d; by k 0..5: %s; by w-order: %s\n' % (NL, [cnt_k[k] for k in range(6)], [cnt_a.get(a, 0) for a in range(max(cnt_a) + 1)]))
    fh.write('equation slots %d; by w-order 0..15: %s\n' % (sum(slots.values()), [slots.get(a, 0) for a in range(16)]))
    fh.write('F2 (pure w-power coefficients of lambda_{ij,0}): %s\n' % ['L(%d)=lam_%d%d_0[w^%d]' % (n_, i, j, a) for (n_, i, j, a) in F2])
    fh.write('F3 (pure (y2,y3) coefficients of lambda_{ij,0}): %s\n' % {('%d%d' % ij): ['L(%d):y2^%d*y3^%d' % (n_, e, f) for (n_, (e, f)) in lst] for ij, lst in F3.items()})
    fh.write('F3 remainders mod y2^2 - D*y3^2:\n')
    for (ij, r1, r0) in F3rem:
        fh.write('  %s: r1 = %s ; r0 = %s\n' % (ij, r1, r0))
with open(os.path.join(HERE, 'v_unknowns.txt'), 'w') as fh:
    fh.write('# name | type | entry/monomial | w-order\n')
    fh.write('D | D | A_11 = D*y3, A_44 = -D*y3, A_05: -D*y3^2 | 0\n')
    bmon = ['q1:y1y2', 'q1:y1y3', 'q2:y0y2', 'q2:y0y3', 'q2:y1y2', 'q2:y1y3', 'q3:y0^2', 'q3:y0y1', 'q3:y1^2', 'q4:y0^2', 'q4:y2y3', 'q4:y3^2']
    for k_, s in enumerate(bmon):
        y0p = s.count('y0') + (1 if 'y0^2' in s else 0)
        fh.write('b(%d) | b | %s | %d\n' % (k_ + 1, s, 1 + 2 * y0p))
    gmon = ['y3^3', 'y1^2y3', 'y2y3^2', 'y1^2y2', 'y2^2y3', 'y2^3', 'y0y1y3', 'y0y1y2', 'y0^2y3', 'y0^2y2']
    for k_, s in enumerate(gmon):
        y0p = s.count('y0') + (1 if 'y0^2' in s else 0)
        fh.write('g(%d) | g | G:%s | %d\n' % (k_ + 1, s, 2 + 2 * y0p))
    for (n_, i, j, k, m) in lam:
        fh.write('L(%d) | lambda | lambda_%d%d,%d : %s | %d\n' % (n_, i, j, k, mstr(m), m[0]))
# python-side data for later scripts
import json
with open(os.path.join(HERE, 'v_data.json'), 'w') as fh:
    json.dump({'NL': NL, 'lam': lam, 'F2': F2, 'F3': {('%d%d' % ij): lst for ij, lst in F3.items()},
               'F3rem': F3rem, 'pairs': pairs}, fh)
print('lambda', NL, [cnt_k[k] for k in range(6)], [cnt_a.get(a, 0) for a in range(max(cnt_a) + 1)])
print('slots', sum(slots.values()), [slots.get(a, 0) for a in range(16)])
print('F2', len(F2), 'F3', sum(len(v) for v in F3.values()))
