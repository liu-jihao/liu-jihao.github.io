#!/usr/bin/env python3
"""Independent generator for the adversarial verification of case n1c1.
Written from the spec (plan/case_n1c1.md, plan/case_common.md) and the sources
(CCO94 7.3-7.5 p.122-123, CD89 p.35, 37, 52, 54 page images), not from the case's scripts.

Usage: python3 gen_v.py slice     -> V1_slice.sing   (the normal form of the case: 12 kept q's, A00 = y0*G, A05 = Q)
       python3 gen_v.py general   -> V1_general.sing (all 20 q's, A00 = y0*G + h(y1,y2,y3), A05 = Q + y0*L)
Also writes V0_counts.txt (monomial counts) and V0_unknowns_<variant>.txt (table of unknowns).

Conventions (own):
 base ring Q[w,y1,y2,y3], weights 1,2,2,2, y0 = w^2; sigma: w->-w, y1->y1, y2->-y2, y3->-y3.
 generators v0..v5, degrees cd = (0,3,3,3,3,4), sigma-signs dl = (+,-,-,+,+,-).
 entry A_jk: weight 8-cd_j-cd_k, sign -dl_j*dl_k.  lambda_{ij,k}: weight cd_i+cd_j-cd_k, sign dl_i*dl_j*dl_k.
 Ring variable order (different from any template): x(1..NX) [lambda], A-coefficients, D, y3, y2, y1, w, tt.
 Lambda numbering: k descending (5..0), then pairs (i,j) lexicographic, then monomials by decreasing w-exponent.
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
variant = sys.argv[1] if len(sys.argv) > 1 else 'slice'
assert variant in ('slice', 'general')

cd = [0, 3, 3, 3, 3, 4]
dl = [1, -1, -1, 1, 1, -1]


def wy_monos(weight, sign):
    """all monomials w^a y1^b y2^e y3^f with a + 2(b+e+f) = weight and (-1)^(a+e+f) = sign"""
    res = []
    if weight < 0:
        return res
    for a in range(weight, -1, -1):
        if (weight - a) % 2:
            continue
        n = (weight - a) // 2
        for b in range(n, -1, -1):
            for e in range(n - b, -1, -1):
                f = n - b - e
                if (-1) ** (a + e + f) == sign:
                    res.append((a, b, e, f))
    return res


def ystr(m):
    a, b, e, f = m
    parts = []
    for v, p in (('w', a), ('y1', b), ('y2', e), ('y3', f)):
        if p == 1:
            parts.append(v)
        elif p > 1:
            parts.append('%s^%d' % (v, p))
    return '*'.join(parts) if parts else '1'


def Y_monos(deg, sign, no_y0=False):
    """monomials y0^a0 y1^a1 y2^a2 y3^a3 of degree deg with (-1)^(a2+a3) = sign"""
    res = []
    for a0 in range(deg, -1, -1):
        if no_y0 and a0 > 0:
            continue
        for a1 in range(deg - a0, -1, -1):
            for a2 in range(deg - a0 - a1, -1, -1):
                a3 = deg - a0 - a1 - a2
                if (-1) ** (a2 + a3) == sign:
                    res.append((a0, a1, a2, a3))
    return res


def Ystr(m):
    a0, a1, a2, a3 = m
    return ystr((2 * a0, a1, a2, a3))


def Yname(m):
    a0, a1, a2, a3 = m
    s = ''
    for v, p in (('y0', a0), ('y1', a1), ('y2', a2), ('y3', a3)):
        s += v * p
    return s


pairs = [(i, j) for i in range(1, 6) for j in range(i, 6)]

# ---------------- lambda ansatz ----------------
lam = []
for k in range(5, -1, -1):
    for (i, j) in pairs:
        for m in wy_monos(cd[i] + cd[j] - cd[k], dl[i] * dl[j] * dl[k]):
            lam.append((i, j, k, m))
NX = len(lam)

# ---------------- A-coefficient ansatz ----------------
Gm = Y_monos(3, -1)                    # 10 anti-invariant cubics
qm_all = {1: Y_monos(2, -1), 2: Y_monos(2, -1), 3: Y_monos(2, 1), 4: Y_monos(2, 1)}
kept = {1: [(0, 1, 1, 0), (0, 1, 0, 1)],
        2: [(1, 0, 1, 0), (1, 0, 0, 1), (0, 1, 1, 0), (0, 1, 0, 1)],
        3: [(2, 0, 0, 0), (1, 1, 0, 0), (0, 2, 0, 0)],
        4: [(2, 0, 0, 0), (0, 0, 1, 1), (0, 0, 0, 2)]}
if variant == 'slice':
    qlist = [(i, m) for i in (1, 2, 3, 4) for m in kept[i]]
else:
    qlist = [(i, m) for i in (1, 2, 3, 4) for m in qm_all[i]]
qlist = list(reversed(qlist))          # own numbering: reversed
hm = Y_monos(4, -1, no_y0=True) if variant == 'general' else []   # 6 anti-invariant quartics in y1,y2,y3
lm = [(1, 0, 0, 0), (0, 1, 0, 0)] if variant == 'general' else []  # L = l1*y0 + l2*y1
NQ, NG, NH, NL = len(qlist), len(Gm), len(hm), len(lm)
NA = NQ + NG + NH + NL

# ---------------- counts file (both variants identical for lambda and slots) ----------------
with open(os.path.join(HERE, 'V0_counts.txt'), 'w') as fo:
    byk = [sum(1 for t in lam if t[2] == k) for k in range(6)]
    maxo = max(t[3][0] for t in lam)
    byo = [sum(1 for t in lam if t[3][0] == o) for o in range(maxo + 1)]
    fo.write('lambda coefficients: %d\n' % NX)
    fo.write('by k=0..5: %s\n' % byk)
    fo.write('by w-order 0..%d: %s\n' % (maxo, byo))
    slots = {}
    for (i, j) in pairs:
        for m in wy_monos(8 + cd[i] + cd[j], -dl[i] * dl[j]):
            slots[m[0]] = slots.get(m[0], 0) + 1
    fo.write('slots (monomials of weight 8+c_i+c_j, sign -d_i d_j, over the 15 pairs): %d\n' % sum(slots.values()))
    fo.write('slots by w-order: %s\n' % [slots.get(o, 0) for o in range(max(slots) + 1)])
    # trivial syzygy count
    tri = 0
    trip = {}
    for (i, j) in pairs:
        n = 0
        for mm in range(1, 6):
            n += len(wy_monos(cd[i] + cd[j] - (8 - cd[mm]), -dl[i] * dl[j] * dl[mm]))
        trip[(i, j)] = n
        tri += n
    fo.write('trivial-syzygy parameters: %d; per pair: %s\n' % (tri, trip))
    fo.write('G monomials (%d): %s\n' % (len(Gm), [Yname(m) for m in Gm]))
    for i in (1, 2, 3, 4):
        fo.write('q%d all monomials (%d): %s\n' % (i, len(qm_all[i]), [Yname(m) for m in qm_all[i]]))

# ---------------- table of unknowns ----------------
def worder_lam(t):
    return t[3][0]


with open(os.path.join(HERE, 'V0_unknowns_%s.txt' % variant), 'w') as fo:
    fo.write('# name type monomial w-order\n')
    for n, (i, j, k, m) in enumerate(lam, 1):
        fo.write('x(%d) lambda_{%d%d,%d} %s %d\n' % (n, i, j, k, ystr(m), m[0]))
    for n, (i, m) in enumerate(qlist, 1):
        fo.write('qq(%d) q%d %s %d\n' % (n, i, Yname(m), 1 + 2 * m[0]))
    for n, m in enumerate(Gm, 1):
        fo.write('gg(%d) G %s %d\n' % (n, Yname(m), 2 + 2 * m[0]))
    for n, m in enumerate(hm, 1):
        fo.write('hh(%d) h(A00 mod y0) %s %d\n' % (n, Yname(m), 0))
    for n, m in enumerate(lm, 1):
        fo.write('ll(%d) L(A05=Q+y0*L) %s %d\n' % (n, Yname(m), 2 + 2 * m[0]))
    fo.write('D D 1 0\n')

# ---------------- Singular script ----------------
Avars = ['qq(%d)' % n for n in range(1, NQ + 1)] + ['gg(%d)' % n for n in range(1, NG + 1)] + \
        ['hh(%d)' % n for n in range(1, NH + 1)] + ['ll(%d)' % n for n in range(1, NL + 1)]
small_vars = ', '.join(Avars + ['D', 'y3', 'y2', 'y1', 'w', 'tt'])
big_vars = 'x(1..%d), ' % NX + small_vars
NV = NX + NA + 6
iD = NX + NA + 1
iw = NX + NA + 5


def poly_from(coefnames, monos, fmt):
    terms = ['%s*%s' % (c, fmt(m)) for c, m in zip(coefnames, monos)]
    return ' + '.join(terms) if terms else '0'


Gpoly = poly_from(['gg(%d)' % n for n in range(1, NG + 1)], Gm, Ystr)
qpolys = {}
for i in (1, 2, 3, 4):
    names = ['qq(%d)' % (n + 1) for n, (ii, m) in enumerate(qlist) if ii == i]
    monos = [m for (ii, m) in qlist if ii == i]
    qpolys[i] = poly_from(names, monos, Ystr)
hpoly = poly_from(['hh(%d)' % n for n in range(1, NH + 1)], hm, Ystr)
Lpoly = poly_from(['ll(%d)' % n for n in range(1, NL + 1)], lm, Ystr)

A00 = 'y0*G' + (' + h' if variant == 'general' else '')
A05 = 'Q' + (' + y0*L' if variant == 'general' else '')
Arows = [
    [A00, 'w*q1', 'w*q2', 'w*q3', 'w*q4', A05],
    ['w*q1', 'D*y3', 'y2', 'y1', '0', 'w'],
    ['w*q2', 'y2', 'y3', 'y0', 'y1', '0'],
    ['w*q3', 'y1', 'y0', '-y3', 'y2', '0'],
    ['w*q4', '0', 'y1', 'y2', '-D*y3', '0'],
    [A05, 'w', '0', '0', '0', '0'],
]

sig_imgs_small = Avars + ['D', '-y3', '-y2', 'y1', '-w', 'tt']
hom_imgs_small = Avars + ['D', 'tt^2*y3', 'tt^2*y2', 'tt^2*y1', 'tt*w', 'tt']
xs = ['x(%d)' % n for n in range(1, NX + 1)]
sig_imgs_big = xs + sig_imgs_small
hom_imgs_big = xs + hom_imgs_small
zero_x_imgs = ['0'] * NX + Avars + ['D', 'y3', 'y2', 'y1', 'w', 'tt']

L = []
P = L.append
P('// %s -- generated by gen_v.py (variant %s). Independent verification of case n1c1.' % ('V1_%s.sing' % variant, variant))
P('system("--ticks-per-sec", 1000);')
P('int T0 = rtimer;')
P('int i, j, k, r, c, rr, cc, p, t, s, e, u, bad, bad2, ne, nd, dup, crossdup, nL, nsub, best, nzero, incons, nleft, nlam, nAc, nfA, nfX;')
P('intvec cd = 0,3,3,3,3,4; intvec dl = 1,-1,-1,1,1,-1;')
P('int NX = %d; int NA = %d; int iD = %d; int iw = %d;' % (NX, NA, iD, iw))
P('ring S = 0, (%s), dp;' % small_vars)
P('poly y0 = w^2;')
P('poly Q = y2^2 - D*y3^2 - y1^2;')
P('poly G = %s;' % Gpoly)
for i in (1, 2, 3, 4):
    P('poly q%d = %s;' % (i, qpolys[i]))
if variant == 'general':
    P('poly h = %s;' % hpoly)
    P('poly L = %s;' % Lpoly)
P('matrix A[6][6] = ' + ', '.join(', '.join(r) for r in Arows) + ';')
P('map sgS = S, %s;' % ', '.join(sig_imgs_small))
P('map hmS = S, %s;' % ', '.join(hom_imgs_small))
P('poly a;')
P('bad = 0;')
P('for (j=1; j<=6; j++) { for (k=1; k<=6; k++) { a = A[j,k];')
P('  if (A[j,k] != A[k,j]) { bad++; "not symmetric", j-1, k-1; }')
P('  if (sgS(a) != -dl[j]*dl[k]*a) { bad++; "sigma-sign fails at", j-1, k-1; }')
P('  if (hmS(a) != tt^(8-cd[j]-cd[k])*a) { bad++; "weight fails at", j-1, k-1; } } }')
P('"CHECK entries of A: symmetric, sigma-sign -d_j d_k, weight 8-c_j-c_k; failures =", bad;')
P('// special net alpha (rows/cols v1..v4 at w = 0) and Q')
P('matrix al[4][4]; for (j=1;j<=4;j++) { for (k=1;k<=4;k++) { al[j,k] = subst(A[j+1,k+1], w, 0); } }')
P('"CHECK det(alpha) - Q^2 == 0:", det(al) - Q^2 == 0;')
P('ideal IQ = std(ideal(Q)); bad = 0;')
P('matrix aal[4][4]; matrix mm3[3][3];')
P('for (j=1;j<=4;j++) { for (k=1;k<=4;k++) { rr = 0; for (r=1;r<=4;r++) { if (r != k) { rr++; cc = 0; for (c=1;c<=4;c++) { if (c != j) { cc++; mm3[rr,cc] = al[r,c]; } } } } aal[j,k] = (-1)^(j+k)*det(mm3); if (reduce(aal[j,k], IQ) != 0) { bad++; } } }')
P('"CHECK Q divides every entry of adj(alpha): failures =", bad;')
P('matrix I4[4][4]; for (j=1;j<=4;j++) { I4[j,j] = 1; }')
P('"CHECK alpha*adj(alpha) == det(alpha)*I:", al*aal == det(al)*I4;')
P('// adjugate of A by cofactors: B[i,j] = (-1)^(i+j) det(A without row j and column i)')
P('matrix B[6][6]; matrix m5[5][5];')
P('for (i=1;i<=6;i++) { for (j=1;j<=6;j++) { rr = 0; for (r=1;r<=6;r++) { if (r != j) { rr++; cc = 0; for (c=1;c<=6;c++) { if (c != i) { cc++; m5[rr,cc] = A[r,c]; } } } } B[i,j] = (-1)^(i+j)*det(m5); } }')
P('poly dA = det(A);')
P('matrix I6[6][6]; for (j=1;j<=6;j++) { I6[j,j] = 1; }')
P('"CHECK A*B == det(A)*I:", A*B == dA*I6;')
P('"CHECK B symmetric:", B == transpose(B);')
P('"det A sigma-invariant:", sgS(dA) == dA, " weight 16:", hmS(dA) == tt^16*dA;')
P('"leading terms: B05 + Q^3 == 0 mod w:", subst(B[1,6] + Q^3, w, 0) == 0;')
P('for (k=2;k<=5;k++) { "B0k mod w == 0 for k =", k-1, ":", subst(B[1,k], w, 0) == 0; }')
P('"time (ms) after adjugate:", rtimer - T0;')
P('// ---------------- big ring ----------------')
P('ring R = 0, (%s), dp;' % big_vars)
P('matrix B = imap(S, B);')
P('map sgR = R, %s;' % ', '.join(sig_imgs_big))
P('map hmR = R, %s;' % ', '.join(hom_imgs_big))
P('map zx = R, %s;' % ', '.join(zero_x_imgs))
P('intvec PPI = %s;' % ','.join(str(i) for (i, j) in pairs))
P('intvec PPJ = %s;' % ','.join(str(j) for (i, j) in pairs))
P('matrix LAM[15][6]; poly ff, gg0;')
# lambda polys
lamterms = {}
for n, (i, j, k, m) in enumerate(lam, 1):
    lamterms.setdefault((pairs.index((i, j)) + 1, k + 1), []).append('x(%d)*%s' % (n, ystr(m)))
for (pp, kk), terms in sorted(lamterms.items()):
    P('LAM[%d,%d] = %s;' % (pp, kk, ' + '.join(terms)))
P('bad = 0;')
P('for (p=1;p<=15;p++) { for (k=0;k<=5;k++) { i = PPI[p]; j = PPJ[p];')
P('  ff = LAM[p,k+1]; if (ff != 0) { if (sgR(ff) != dl[i+1]*dl[j+1]*dl[k+1]*ff) { bad++; } if (hmR(ff) != tt^(cd[i+1]+cd[j+1]-cd[k+1])*ff) { bad++; } } } }')
P('"CHECK lambda polys: sign d_i d_j d_k and weight c_i+c_j-c_k; failures =", bad;')
P('// rank-condition polynomials E_ij = B_ij - sum_k lambda_{ij,k} B_0k and their coefficients')
P('list EF; ideal EQ; intvec EO; intvec EP; matrix C; poly E; ne = 0; bad = 0;')
P('for (p=1;p<=15;p++) { i = PPI[p]; j = PPJ[p]; E = B[i+1,j+1];')
P('  for (k=0;k<=5;k++) { E = E - LAM[p,k+1]*B[1,k+1]; }')
P('  EF[p] = E;')
P('  if (hmR(E) != tt^(8+cd[i+1]+cd[j+1])*E) { bad++; "E weight fails", p; }')
P('  if (sgR(E) != -dl[i+1]*dl[j+1]*E) { bad++; "E sign fails", p; }')
P('  if (E != 0) { C = coef(E, w*y1*y2*y3);')
P('    for (t=1;t<=ncols(C);t++) { if (C[2,t] != 0) { ne++; EQ[ne] = C[2,t]; EO[ne] = leadexp(C[1,t])[iw]; EP[ne] = p; } } } }')
P('"CHECK E_ij homogeneous of weight 8+c_i+c_j with sign -d_i d_j; failures =", bad;')
P('"nonzero coefficient equations:", ne;')
P('intvec cnt = 0:20; for (t=1;t<=ne;t++) { cnt[EO[t]+1] = cnt[EO[t]+1] + 1; }')
P('"nonzero equations by w-order 0..19:", cnt;')
P('// dedup up to scalars (global)')
P('ideal EN; for (t=1;t<=ne;t++) { EN[t] = EQ[t]/leadcoef(EQ[t]); }')
P('ideal ED; intvec EDO; intvec EDP; nd = 0; crossdup = 0;')
P('for (t=1;t<=ne;t++) { dup = 0; for (u=1;u<=nd;u++) { if (size(ED[u]) == size(EN[t])) { if (ED[u] == EN[t]) { dup = u; break; } } }')
P('  if (dup == 0) { nd++; ED[nd] = EN[t]; EDO[nd] = EO[t]; EDP[nd] = EP[t]; } else { if (EDO[dup] != EO[t]) { crossdup++; } } }')
P('"distinct equations up to scalar (global):", nd, " duplicates across different w-orders:", crossdup;')
P('cnt = 0:20; for (t=1;t<=nd;t++) { cnt[EDO[t]+1] = cnt[EDO[t]+1] + 1; }')
P('"distinct equations by w-order 0..19:", cnt;')
P('// order structure: max total w-order of unknowns in each term <= w-order of the equation is checked in Python')
P('write(":w V1_%s_eqs_ideal.txt", string(ED));' % variant)
P('write(":w V1_%s_eqs_orders.txt", string(EDO));' % variant)
P('"time (ms) after equations:", rtimer - T0;')
# elimination procedure
P('''proc cand(poly f, int nx0, int na0)
{
  poly lin = jet(f,1) - jet(f,0);
  int s, idx, cls; poly m, dd;
  int bcls = 0; poly bm = 0; number bc = 0;
  for (s=1; s<=size(lin); s++)
  {
    m = leadmonom(lin[s]); idx = rvar(m); cls = 0;
    if ((idx >= 1) && (idx <= nx0)) { cls = 1; }
    if ((idx > nx0) && (idx <= nx0+na0)) { cls = 2; }
    if (cls > 0)
    {
      if ((bcls == 0) || (cls < bcls))
      {
        dd = diff(f, m);
        if (deg(dd) == 0) { bcls = cls; bm = m; bc = leadcoef(dd); if (cls == 1) { break; } }
      }
    }
  }
  return(list(bcls, bm, bc));
}
proc applymap(poly f, ideal MM)
{
  def R0 = basering; map phi = R0, MM; return(phi(f));
}''')
P('// global constant-pivot elimination; priority: lambda before A-coefficient (D never), then fewest terms, then first found')
P('ideal Lc = ED; nL = nd;')
P('intvec ccls = 0:nL; list cvar; list ccoef; intvec csz = 0:nL; list cres;')
P('for (e=1;e<=nL;e++) { cres = cand(Lc[e], NX, NA); ccls[e] = cres[1]; cvar[e] = cres[2]; ccoef[e] = cres[3]; csz[e] = size(Lc[e]); }')
P('list SUBV; list SUBE; list SUBC; intvec SUBF; intvec SUBK; intvec SUBO; intvec SUBS; nsub = 0; poly uu, ex; number cf; nzero = 0; incons = 0;')
P('while (1) { best = 0;')
P('  for (e=1;e<=nL;e++) { if (ccls[e] > 0) { if (best == 0) { best = e; } else { if ((ccls[e] < ccls[best]) || ((ccls[e] == ccls[best]) && (csz[e] < csz[best]))) { best = e; } } } }')
P('  if (best == 0) { break; }')
P('  uu = cvar[best]; cf = ccoef[best];')
P('  if (deg(diff(Lc[best], uu)) != 0) { "ERROR: non-constant pivot"; quit; }')
P('  ex = uu - Lc[best]/cf;')
P('  if (diff(ex, uu) != 0) { "ERROR: expression contains pivot variable"; quit; }')
P('  nsub++; SUBV[nsub] = uu; SUBE[nsub] = ex; SUBC[nsub] = cf; SUBF[nsub] = best; SUBK[nsub] = ccls[best]; SUBO[nsub] = EDO[best]; SUBS[nsub] = csz[best];')
P('  Lc[best] = 0; ccls[best] = 0;')
P('  for (e=1;e<=nL;e++) { if (Lc[e] != 0) { if (diff(Lc[e], uu) != 0) { Lc[e] = subst(Lc[e], uu, ex);')
P('        if (Lc[e] == 0) { ccls[e] = 0; nzero++; } else { if (deg(Lc[e]) == 0) { incons++; "INCONSISTENT equation", e; }')
P('          cres = cand(Lc[e], NX, NA); ccls[e] = cres[1]; cvar[e] = cres[2]; ccoef[e] = cres[3]; csz[e] = size(Lc[e]); } } } } }')
P('"elimination done: substitutions =", nsub, " equations reduced to zero =", nzero, " inconsistent =", incons;')
P('nleft = 0; for (e=1;e<=nL;e++) { if (Lc[e] != 0) { nleft++; } }')
P('"LEFTOVER equations:", nleft;')
P('if (nleft > 0) { for (e=1;e<=nL;e++) { if (Lc[e] != 0) { "left eq", e, "w-order", EDO[e], ":", Lc[e]; "  factors:", factorize(Lc[e]); } } }')
P('"time (ms) after elimination:", rtimer - T0;')
P('// pivot statistics')
P('nlam = 0; nAc = 0; for (k=1;k<=nsub;k++) { if (SUBK[k] == 1) { nlam++; } else { nAc++; } }')
P('"eliminated lambda:", nlam, " eliminated A-coefficients:", nAc;')
P('string pivs = ""; for (k=1;k<=nsub;k++) { pivs = pivs + string(SUBC[k]) + " "; }')
P('"pivot constants:", pivs;')
P('string elimA = ""; for (k=1;k<=nsub;k++) { if (SUBK[k] == 2) { elimA = elimA + string(SUBV[k]) + "[step " + string(k) + ", eq w-order " + string(SUBO[k]) + ", pivot " + string(SUBC[k]) + "] "; } }')
P('"eliminated A-coefficients (with step):", elimA;')
P('link lk3 = ":w V1_%s_subs.txt"; for (k=1;k<=nsub;k++) { write(lk3, string(k) + " | " + string(SUBV[k]) + " | class " + string(SUBK[k]) + " | pivot " + string(SUBC[k]) + " | eq " + string(SUBF[k]) + " (w-order " + string(SUBO[k]) + ", terms " + string(SUBS[k]) + ") | " + string(SUBE[k])); } close(lk3);' % variant)
P('// free unknowns')
P('intvec iselim = 0:nvars(R); for (k=1;k<=nsub;k++) { iselim[rvar(SUBV[k])] = 1; }')
P('string freeA = ""; string freeX = ""; nfA = 0; nfX = 0;')
P('for (k=1;k<=NX;k++) { if (iselim[k] == 0) { nfX++; freeX = freeX + string(var(k)) + " "; } }')
P('for (k=NX+1;k<=NX+NA+1;k++) { if (iselim[k] == 0) { nfA++; freeA = freeA + string(var(k)) + " "; } }')
P('"FREE A-coefficients (f, D included):", nfA, ":", freeA;')
P('"FREE lambda-coefficients:", nfX;')
P('link lk1 = ":w V1_%s_free.txt"; write(lk1, "freeA: " + freeA); write(lk1, "freeX: " + freeX); close(lk1);' % variant)
P('// back-substitution')
P('ideal MM = maxideal(1);')
P('for (k=nsub;k>=1;k--) { MM[rvar(SUBV[k])] = applymap(SUBE[k], MM); }')
P('"time (ms) after back-substitution:", rtimer - T0;')
P('// sanity: final expressions contain only free variables')
P('ideal ZE = maxideal(1); for (k=1;k<=nsub;k++) { ZE[rvar(SUBV[k])] = 0; } map zel = R, ZE;')
P('bad = 0; for (k=1;k<=nsub;k++) { ff = MM[rvar(SUBV[k])]; if (zel(ff) != ff) { bad++; } }')
P('"CHECK final expressions involve only free unknowns: failures =", bad;')
P('// graph check: eliminated A-coefficients free of lambda')
P('bad = 0; link lk2 = ":w V1_%s_graph.txt";' % variant)
P('for (k=NX+1;k<=NX+NA;k++) { if (iselim[k] == 1) { ff = MM[k]; if (zx(ff) != ff) { bad++; } write(lk2, string(var(k)) + " = " + string(ff) + ";"); } }')
P('close(lk2);')
P('"GRAPH CHECK (eliminated A-coefficients free of lambda): failures =", bad;')
P('write(":w V1_%s_param_ideal.txt", string(MM));' % variant)
P('// verification (ii): the parametrization kills all coefficient equations and all 15 full E_ij')
P('map PHI = R, MM;')
P('bad = 0; for (t=1;t<=ne;t++) { ff = EQ[t]; if (PHI(ff) != 0) { bad++; } }')
P('"VERIFY (ii-a) all nonzero coefficient equations vanish on the parametrization: failures =", bad;')
P('"time (ms) after (ii-a):", rtimer - T0;')
P('// verification (iii): each pivot equation is the image of its original equation under the earlier substitutions')
P('bad = 0; bad2 = 0;')
P('for (k=1;k<=nsub;k++) { ff = ED[SUBF[k]]; for (u=1;u<k;u++) { ff = subst(ff, SUBV[u], SUBE[u]); }')
P('  if (ff != SUBC[k]*(SUBV[k] - SUBE[k])) { bad++; }')
P('  for (u=1;u<=k;u++) { if (diff(SUBE[k], SUBV[u]) != 0) { bad2++; } } }')
P('"VERIFY (iii) pivot k = image of original equation under substitutions 1..k-1, equal to c_k*(u_k - e_k): failures =", bad;')
P('"VERIFY (iii) e_k free of u_1..u_k: failures =", bad2;')
P('"time (ms) after (iii):", rtimer - T0;')
P('bad = 0; for (p=1;p<=15;p++) { ff = EF[p]; if (PHI(ff) != 0) { bad++; "E fails", p; } }')
P('"VERIFY (ii-b) the 15 polynomials E_ij vanish identically on the parametrization: failures =", bad;')
P('"time (ms) total:", rtimer - T0;')
P('"memory (bytes) memory(1), memory(2):", memory(1), memory(2);')
P('quit;')

with open(os.path.join(HERE, 'V1_%s.sing' % variant), 'w') as fo:
    fo.write('\n'.join(L) + '\n')
print('variant', variant, 'NX', NX, 'NQ', NQ, 'NG', NG, 'NH', NH, 'NL', NL, 'NA', NA, 'nvars', NV)
