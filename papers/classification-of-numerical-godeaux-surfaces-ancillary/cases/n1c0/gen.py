#!/usr/bin/env python3
"""Generator for the Singular scripts of case n1c0 (lifting of N1(D), c = 0).
Notation of plan/case_common.md: coordinates w (wt 1, sigma -), y1 (wt 2, +), y2, y3 (wt 2, -); y0 = w^2.
Generators v0..v5 of degrees c = (0,3,3,3,3,4), sigma-signs delta = (+,-,-,+,+,-).
Unknowns (ring variables 1..394): l(1..371) lambda coefficients, g(1..10), b(1..12), D.
Usage: python3 gen.py   (writes unknowns_table.txt and the .sing scripts listed in SCRIPTS)
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
c = [0, 3, 3, 3, 3, 4]
dl = [1, -1, -1, 1, 1, -1]
r = [8 - x for x in c]

def mons(n, s):
    """monomials w^a y1^b y2^e y3^f of weight n = a + 2(b+e+f), sigma-sign s = (-1)^(a+e+f);
    ordered by a ascending, then (b,e,f) descending lexicographically"""
    out = []
    for a in range(n + 1):
        if (n - a) % 2:
            continue
        m = (n - a) // 2
        lst = []
        for bb in range(m + 1):
            for e in range(m - bb + 1):
                f = m - bb - e
                if (-1) ** (a + e + f) == s:
                    lst.append((a, bb, e, f))
        lst.sort(key=lambda x: (-x[1], -x[2], -x[3]))
        out.extend(lst)
    return out

def mstr(ex):
    a, bb, e, f = ex
    parts = []
    for v, p in (('w', a), ('y1', bb), ('y2', e), ('y3', f)):
        if p == 1:
            parts.append(v)
        elif p > 1:
            parts.append('%s^%d' % (v, p))
    return '*'.join(parts) if parts else '1'

# ---------------- unknowns ----------------
LAM = {}      # (i,j,k) -> list of (index, monomial exponent)
UNK = []      # list of dicts: idx, name, type, mon, order
idx = 0
for i in range(1, 6):
    for j in range(i, 6):
        for k in range(6):
            ms = mons(c[i] + c[j] - c[k], dl[i] * dl[j] * dl[k])
            LAM[(i, j, k)] = []
            for ex in ms:
                idx += 1
                LAM[(i, j, k)].append((idx, ex))
                UNK.append(dict(idx=idx, name='l(%d)' % idx, type='lambda', mon='lambda_{%d%d,%d}: %s' % (i, j, k, mstr(ex)), order=ex[0]))
assert idx == 371, idx
Gmons = [('y3^3', (0, 0, 0, 3)), ('y1^2*y3', (0, 2, 0, 1)), ('y2*y3^2', (0, 0, 1, 2)), ('y1^2*y2', (0, 2, 1, 0)),
         ('y2^2*y3', (0, 0, 2, 1)), ('y2^3', (0, 0, 3, 0)), ('y0*y1*y3', (1, 1, 0, 1)), ('y0*y1*y2', (1, 1, 1, 0)),
         ('y0^2*y3', (2, 0, 0, 1)), ('y0^2*y2', (2, 0, 1, 0))]
for n, (s, ex) in enumerate(Gmons, 1):
    idx += 1
    UNK.append(dict(idx=idx, name='g(%d)' % n, type='g', mon='G: %s' % s, order=2 + 2 * ex[0]))
Bmons = [('q1', 'y1*y2', 0), ('q1', 'y1*y3', 0), ('q2', 'y0*y2', 1), ('q2', 'y0*y3', 1), ('q2', 'y1*y2', 0), ('q2', 'y1*y3', 0),
         ('q3', 'y0^2', 2), ('q3', 'y0*y1', 1), ('q3', 'y1^2', 0), ('q4', 'y0^2', 2), ('q4', 'y2*y3', 0), ('q4', 'y3^2', 0)]
for n, (q, s, k0) in enumerate(Bmons, 1):
    idx += 1
    UNK.append(dict(idx=idx, name='b(%d)' % n, type='b', mon='%s: %s' % (q, s), order=1 + 2 * k0))
idx += 1
UNK.append(dict(idx=idx, name='D', type='D', mon='parameter of N1', order=0))
assert idx == 394
TYPCODE = {'lambda': 1, 'g': 2, 'b': 3, 'D': 4}

# special lambda coefficients for the sigma-fixed-point tests (case_common section 5)
PUREW = []   # pure w-power coefficient of lambda_{ij,0}
PUREY = []   # coefficients of monomials y2^e y3^f (no w, no y1) of lambda_{ij,0}
for i in range(1, 6):
    for j in range(i, 6):
        for (ix, ex) in LAM[(i, j, 0)]:
            a, bb, e, f = ex
            if bb == 0 and e == 0 and f == 0:
                PUREW.append((ix, i, j, ex))
            if a == 0 and bb == 0:
                PUREY.append((ix, i, j, ex))

def table():
    L = ['# Unknowns of case n1c0 (ring variables 1..394, in this order). order = w-order (case_common section 2).',
         '# idx  name      type    order  monomial']
    for u in UNK:
        L.append('%4d  %-8s  %-6s  %d      %s' % (u['idx'], u['name'], u['type'], u['order'], u['mon']))
    from collections import Counter
    byk = Counter(); byo = Counter()
    for (i, j, k), lst in LAM.items():
        byk[k] += len(lst)
        for _, ex in lst:
            byo[ex[0]] += 1
    L.append('# lambda count %d; by k: %s; by w-order: %s' % (sum(byk.values()), dict(sorted(byk.items())), dict(sorted(byo.items()))))
    ao = Counter((u['type'], u['order']) for u in UNK if u['type'] != 'lambda')
    L.append('# A-unknowns by (type, order): %s' % dict(sorted(ao.items())))
    L.append('# pure w-power coefficients of lambda_{ij,0} (F2 point test): ' + ', '.join('l(%d)=lambda_{%d%d,0}[%s]' % (ix, i, j, mstr(ex)) for ix, i, j, ex in PUREW))
    L.append('# pure (y2,y3) coefficients of lambda_{ij,0} (F3 point test): ' + ', '.join('l(%d)=lambda_{%d%d,0}[%s]' % (ix, i, j, mstr(ex)) for ix, i, j, ex in PUREY))
    # equation slots
    eqo = Counter(); tot = 0
    for i in range(1, 6):
        for j in range(i, 6):
            for ex in mons(8 + c[i] + c[j], -dl[i] * dl[j]):
                eqo[ex[0]] += 1; tot += 1
    L.append('# equation slots (monomials of weight 8+c_i+c_j, sign -d_i d_j, i<=j): %d; by w-order: %s' % (tot, dict(sorted(eqo.items()))))
    syz = 0
    for i in range(1, 6):
        for j in range(i, 6):
            for m in range(1, 6):
                n = c[i] + c[j] - r[m]
                if n >= 0:
                    syz += len(mons(n, -dl[i] * dl[j] * dl[m]))
    L.append('# trivial-syzygy directions: %d' % syz)
    return '\n'.join(L) + '\n', eqo

# ---------------- Singular code blocks ----------------
def lam_poly(i, j, k):
    return ' + '.join('l(%d)*%s' % (ix, mstr(ex)) for ix, ex in LAM[(i, j, k)])

def hdr(title):
    return '''// %s
// Case n1c0 (N1(D), c = 0). Generated by gen.py. Exact computation over Q unless stated.
system("--ticks-per-sec", 1000); int T0 = rtimer;
option(redSB);
ring R = 0, (l(1..371), g(1..10), b(1..12), D, w, y1, y2, y3, t), dp;
''' % title

def build(delays=False, verbose=True):
    S = []
    S.append('poly y0 = w^2;')
    S.append('poly G = ' + ' + '.join('g(%d)*%s' % (n, s) for n, (s, ex) in enumerate(Gmons, 1)) + ';')
    S.append('poly q1 = b(1)*y1*y2 + b(2)*y1*y3;')
    S.append('poly q2 = b(3)*y0*y2 + b(4)*y0*y3 + b(5)*y1*y2 + b(6)*y1*y3;')
    S.append('poly q3 = b(7)*y0^2 + b(8)*y0*y1 + b(9)*y1^2;')
    S.append('poly q4 = b(10)*y0^2 + b(11)*y2*y3 + b(12)*y3^2;')
    S.append('poly Q = y2^2 - D*y3^2 - y1^2;')
    S.append('''matrix A[6][6] =
  y0*G, w*q1, w*q2, w*q3, w*q4, Q,
  w*q1, D*y3, y2,   y1,   0,    w,
  w*q2, y2,   y3,   0,    y1,   0,
  w*q3, y1,   0,    -y3,  y2,   0,
  w*q4, 0,    y1,   y2,   -D*y3, 0,
  Q,    w,    0,    0,    0,    0;''')
    S.append('"A symmetric:", A == transpose(A);')
    S.append('proc cofac(matrix M, int rr, int cc) { matrix S5[5][5]; int a, bb, p, q; a = 0; for (p=1; p<=6; p++) { if (p != rr) { a++; bb = 0; for (q=1; q<=6; q++) { if (q != cc) { bb++; S5[a,bb] = M[p,q]; } } } } return((-1)^(rr+cc) * det(S5)); }')
    S.append('matrix B0[1][6]; int k; for (k=1;k<=6;k++){ B0[1,k] = cofac(A,k,1); }   // B0[1,k+1] = B_{0k}')
    # orders and types
    S.append('intvec ORD = ' + ','.join(str(u['order']) for u in UNK) + ';')
    S.append('intvec TYP = ' + ','.join(str(TYPCODE[u['type']]) for u in UNK) + ';')
    if delays:
        S.append('''// ---- delays (case_common section 2): leading terms of the row-0 cofactors
matrix al[4][4]; int ii, jj; for (ii=1;ii<=4;ii++){ for (jj=1;jj<=4;jj++){ al[ii,jj] = subst(A[ii+1,jj+1], w, 0); } }
"alpha = A[v1..v4] mod w:"; print(al);
"det(alpha) - Q^2 =", det(al) - Q^2;
"B_05 mod w minus (-Q^3):", subst(B0[1,6], w, 0) + Q^3, "   minus (+Q^3):", subst(B0[1,6], w, 0) - Q^3;
proc cof4(matrix M, int rr, int cc) { matrix S3[3][3]; int a, bb, p, q; a = 0; for (p=1; p<=4; p++) { if (p != rr) { a++; bb = 0; for (q=1; q<=4; q++) { if (q != cc) { bb++; S3[a,bb] = M[p,q]; } } } } return((-1)^(rr+cc) * det(S3)); }
matrix CW; poly adjk1;
for (k=1;k<=4;k++) {
  CW = coeffs(B0[1,k+1], w); adjk1 = cof4(al, 1, k);   // adj(alpha)_{k1} = cofactor (1,k)
  "B_0" + string(k) + ": w^0-coefficient =", CW[1,1], "; w^1-coefficient - Q*adj(alpha)_{k1} =", CW[2,1] - Q*adjk1, "; w^1-coefficient + Q*adj(alpha)_{k1} =", CW[2,1] + Q*adjk1;
}
matrix a4[4][4]; for (ii=1;ii<=4;ii++){ for (jj=1;jj<=4;jj++){ a4[ii,jj] = A[ii+1,jj+1]; } }
"B_00 + y0*adj(a)_{11} =", B0[1,1] + y0*cof4(a4,1,1);
"adj(alpha)_{k1}, k=1..4:", cof4(al,1,1), cof4(al,1,2), cof4(al,1,3), cof4(al,1,4);
int lw;
for (k=1;k<=5;k++) { CW = coeffs(B0[1,k+1], w); lw = -1;
  for (ii=1; ii<=nrows(CW); ii++) { if (lw < 0) { if (CW[ii,1] != 0) { lw = ii-1; } } }
  "lowest power of w in B_0" + string(k) + ": w^" + string(lw) + "; its coefficient (factorized):"; factorize(CW[lw+1,1]);
}
"ms:", rtimer-T0;''')
    # equations
    S.append('ideal EQ; intvec EQO; list EQT; int i, j, u, nq; poly pol; matrix C; nq = 0;')
    for i in range(1, 6):
        for j in range(i, 6):
            terms = ' + '.join('(%s)*B0[1,%d]' % (lam_poly(i, j, k), k + 1) for k in range(6))
            S.append('pol = cofac(A,%d,%d) - (%s);' % (j + 1, i + 1, terms))
            S.append('if (pol != 0) { C = coef(pol, w*y1*y2*y3); for (u=1;u<=ncols(C);u++){ if (C[2,u] != 0) { nq++; EQ[nq] = C[2,u]; EQO[nq] = leadexp(C[1,u])[395]; EQT[nq] = "E_%d%d[" + string(C[1,u]) + "]"; } } }' % (i, j))
    S.append('"nonzero coefficient equations (before removing duplicates):", nq, " ms:", rtimer-T0;')
    # per order dedup
    S.append('''// remove duplicates up to a scalar, order by order (an order-m equation is a sum of products of
// coefficients whose w-orders add up to at most m, so equations of different orders are never proportional
// unless zero; we nevertheless dedup within each order only and report)
intvec nraw = 0:16; intvec nded = 0:16; ideal Im; ideal Imd; int m; int nd;
ideal EQD; intvec EQDO; list EQDT; nd = 0;
for (m=0; m<=15; m++) {
  Im = 0; int cnt = 0; list tg = list();
  for (u=1;u<=nq;u++){ if (EQO[u]==m) { cnt++; Im[cnt] = EQ[u]; tg[cnt] = EQT[u]; } }
  nraw[m+1] = cnt;
  if (cnt > 0) {
    // keep the first of each class of proportional equations
    intvec keep = 1:cnt;
    for (u=1;u<=cnt;u++) { if (keep[u]) { for (j=u+1;j<=cnt;j++) { if (keep[j]) { if (leadmonom(Im[j]) == leadmonom(Im[u])) { if (Im[j]*leadcoef(Im[u]) - Im[u]*leadcoef(Im[j]) == 0) { keep[j] = 0; } } } } } }
    for (u=1;u<=cnt;u++) { if (keep[u]) { nd++; EQD[nd] = Im[u]; EQDO[nd] = m; EQDT[nd] = tg[u]; nded[m+1] = nded[m+1] + 1; } }
    kill keep;
  }
  kill cnt; kill tg;
}
"nonzero coefficient equations per w-order 0..15:", nraw;
"distinct up to scalar per w-order 0..15:     ", nded;
"total distinct:", nd, " ms:", rtimer-T0;
// filtration check: every order-m equation has weighted degree <= m for the order weights of the unknowns
intvec OW = ORD, 0,0,0,0,0;
int bad = 0; for (u=1;u<=nd;u++) { if (deg(EQD[u], OW) > EQDO[u]) { bad++; } }
"equations violating the order filtration (should be 0):", bad;''')
    if delays:
        S.append('''// first occurrence of each lambda coefficient (delay = first equation order - own order)
intvec FO = -1:394; ideal vv;
for (u=1;u<=nd;u++) { vv = variables(EQD[u]); for (j=1;j<=ncols(vv);j++) { k = rvar(vv[j]); if (k>=1) { if (k<=394) { if (FO[k] < 0 || EQDO[u] < FO[k]) { FO[k] = EQDO[u]; } } } } }
intvec KOF = ''' + ','.join(str(next(kk for (ii, jj, kk), lst in LAM.items() if any(ix == u['idx'] for ix, _ in lst))) for u in UNK if u['type'] == 'lambda') + ''';
int d0 = 0; int d1 = 0; int d2 = 0; int dother = 0; int never = 0;
for (k=1;k<=371;k++) {
  if (FO[k] < 0) { never++; "lambda never occurring:", k; }
  else {
    if (KOF[k]==5) { if (FO[k]-ORD[k]==0) { d0++; } else { dother++; "exception k=5:", k, FO[k], ORD[k]; } }
    if (KOF[k]>=1 && KOF[k]<=4) { if (FO[k]-ORD[k]==1) { d1++; } else { dother++; "exception k=1..4:", k, FO[k], ORD[k]; } }
    if (KOF[k]==0) { if (FO[k]-ORD[k]==2) { d2++; } else { dother++; "exception k=0:", k, FO[k], ORD[k]; } }
  }
}
"lambda delays: k=5 with delay 0:", d0, "; k=1..4 with delay 1:", d1, "; k=0 with delay 2:", d2, "; exceptions:", dother, "; never occurring:", never;
intvec FOA = FO[372..394]; "first occurrence order of g1..g10, b1..b12, D:", FOA;''')
    S.append('''proc cand(poly f, intvec SC, intvec ORD, int mmax)
{
  // best-scored unknown u (index <= 394, order <= mmax) with d f/d u a nonzero rational constant
  poly lin = jet(f, 1); int kk; int vi; int n = 0; intvec cv; poly tt;
  for (kk = 1; kk <= size(lin); kk++) {
    tt = lin[kk];
    if (deg(tt) == 1) { vi = rvar(leadmonom(tt)); if (vi >= 1) { if (vi <= 394) { if (ORD[vi] <= mmax) { n++; cv[n] = vi; } } } }
  }
  intvec res = 0, 0;
  if (n == 0) { return(res); }
  intvec tested = 0:n; int bi; int bs; int s;
  for (s = 1; s <= n; s++) {
    bi = 0; bs = 0;
    for (kk = 1; kk <= n; kk++) { if (tested[kk] == 0) { if (bi == 0 || SC[cv[kk]] < bs) { bi = kk; bs = SC[cv[kk]]; } } }
    tested[bi] = 1;
    if (deg(diff(f, var(cv[bi]))) == 0) { res = cv[bi], bs; return(res); }
  }
  return(res);
}
list TN = "lambda", "g", "b", "D";''')
    return '\n'.join(S)

def trk_block():
    """tracked unknowns: all A-unknowns, the pure-w and pure-(y2,y3) coefficients of lambda_{ij,0}"""
    names = ['g(%d)' % n for n in range(1, 11)] + ['b(%d)' % n for n in range(1, 13)] + ['D']
    names += ['l(%d)' % ix for ix, i, j, ex in PUREW]
    names += ['l(%d)' % ix for ix, i, j, ex in PUREY]
    labels = ['g%d' % n for n in range(1, 11)] + ['b%d' % n for n in range(1, 13)] + ['D']
    labels += ['lambda_%d%d,0[%s]' % (i, j, mstr(ex)) for ix, i, j, ex in PUREW]
    labels += ['lambda_%d%d,0[%s]' % (i, j, mstr(ex)) for ix, i, j, ex in PUREY]
    S = ['ideal TRK = ' + ', '.join(names) + ';',
         'list TRKN = ' + ', '.join('"%s"' % s for s in labels) + ';']
    return '\n'.join(S), len(names)

def elim_orderwise():
    return '''// ---- order-by-order constant-pivot elimination (case_common section 3(b))
// score: lambda first, then highest w-order, then smallest index
intvec SC; for (k=1;k<=394;k++) { SC[k] = (TYP[k]!=1)*100000 + (20-ORD[k])*1000 + k; }
ideal E = EQD; intvec EO = EQDO; int nE = nd;
intvec CV = 0:nE; intvec CS = 0:nE; intvec cvs;
intvec SV; ideal SE; intvec SST; int nel = 0;
int vb; int eb; int bs; poly piv; poly cc; poly ex; int nadd; int nst; int nrem; int nobs; int nobsnew; int mx;
ideal vv2; int typ;
intvec ordcount;
for (m = 0; m <= 15; m++) {
  nadd = 0; for (u=1;u<=nE;u++) { if (EO[u]==m) { if (E[u]!=0) { nadd++; } } }
  for (u=1;u<=nE;u++) { if (EO[u] <= m) { if (E[u] != 0) { cvs = cand(E[u], SC, ORD, m); CV[u] = cvs[1]; CS[u] = cvs[2]; } else { CV[u] = 0; } } }
  nst = 0;
  while (1) {
    vb = 0; eb = 0; bs = 0;
    for (u=1;u<=nE;u++) { if (EO[u] <= m) { if (CV[u] > 0) { if (vb == 0 || CS[u] < bs) { vb = CV[u]; eb = u; bs = CS[u]; } } } }
    if (vb == 0) { break; }
    piv = E[eb]; cc = diff(piv, var(vb));
    ex = var(vb) - piv/cc;
    E[eb] = 0; CV[eb] = 0;
    nel++; SV[nel] = vb; SE[nel] = ex; SST[nel] = m; nst++;
    for (u=1;u<=nE;u++) {
      if (E[u] != 0) { if (diff(E[u], var(vb)) != 0) {
        E[u] = subst(E[u], var(vb), ex);
        if (EO[u] <= m) { if (E[u] != 0) { cvs = cand(E[u], SC, ORD, m); CV[u] = cvs[1]; CS[u] = cvs[2]; } else { CV[u] = 0; } }
      } }
    }
    TRK = subst(TRK, var(vb), ex);
  }
  // stage log
  nrem = 0; nobs = 0; nobsnew = 0;
  for (u=1;u<=nE;u++) { if (EO[u] <= m) { if (E[u] != 0) {
    nrem++; vv2 = variables(E[u]); mx = -1;
    for (j=1;j<=ncols(vv2);j++) { k = rvar(vv2[j]); if (k>=1) { if (k<=394) { if (ORD[k] > mx) { mx = ORD[k]; } } } }
    if (mx < m) { nobs++; if (EO[u]==m) { nobsnew++; } }
  } } }
  // unknowns of order m by type, and eliminated at this stage by (type, order)
  intvec nu = 0:4; for (k=1;k<=394;k++) { if (ORD[k]==m) { nu[TYP[k]] = nu[TYP[k]] + 1; } }
  string elog = ""; intvec cnt2 = 0:40;   // index (type-1)*10 + order + 1
  for (k=1;k<=nel;k++) { if (SST[k]==m) { typ = TYP[SV[k]]; cnt2[(typ-1)*10 + ORD[SV[k]] + 1] = cnt2[(typ-1)*10 + ORD[SV[k]] + 1] + 1; } }
  for (typ=1;typ<=4;typ++) { for (k=0;k<=9;k++) { if (cnt2[(typ-1)*10+k+1] > 0) { elog = elog + TN[typ] + "/o" + string(k) + ":" + string(cnt2[(typ-1)*10+k+1]) + " "; } } }
  "STAGE", m, ": eqs of order m (nonzero after earlier substitutions):", nadd, "; unknowns of order m (lambda,g,b,D):", nu, "; eliminated now:", nst, " [", elog, "]; remaining eqs of order <= m:", nrem, "; involving only unknowns of order < m:", nobs, "(new at this order:", nobsnew, "); ms:", rtimer-T0;
  kill nu; kill elog; kill cnt2;
}
"order-by-order elimination finished: eliminated", nel, "; remaining equations", size(E), "; ms:", rtimer-T0;'''

def elim_global(mode='lambda_index', divide=None, maxpass=12, label='global', ordname='ORD', sstoff=100):
    """global constant-pivot elimination on ideal E (zero entries allowed); mode 'lambda_index': score = index
    (lambda first, by index; then g, b, D by index); 'orderwise_score': lambda first, highest order, index.
    ordname: the intvec of orders passed to cand (a variable with order 999 is never used as a pivot).
    divide: name of a variable assumed invertible (e.g. 'b(8)'): after each pass divide every equation by the
    largest power of it that divides it, and repeat while this changes something."""
    sc = 'SC[k] = k;' if mode == 'lambda_index' else 'SC[k] = (TYP[k]!=1)*100000 + (20-ORD[k])*1000 + k;'
    S = ['''// ---- %s constant-pivot elimination, priority %s%s
intvec SC; for (k=1;k<=394;k++) { %s }
int nE = ncols(E); intvec CV = 0:nE; intvec CS = 0:nE; intvec cvs;
int vb; int eb; int bs; poly piv; poly cc; poly ex; int pass = 0; int progress = 1; int ndiv;
for (u=1;u<=nE;u++) { if (E[u] != 0) { cvs = cand(E[u], SC, %s, 99); CV[u] = cvs[1]; CS[u] = cvs[2]; } }
while (progress) {
  pass++; progress = 0;
  while (1) {
    vb = 0; eb = 0; bs = 0;
    for (u=1;u<=nE;u++) { if (CV[u] > 0) { if (vb == 0 || CS[u] < bs) { vb = CV[u]; eb = u; bs = CS[u]; } } }
    if (vb == 0) { break; }
    piv = E[eb]; cc = diff(piv, var(vb));
    ex = var(vb) - piv/cc;
    E[eb] = 0; CV[eb] = 0;
    nel++; SV[nel] = vb; SE[nel] = ex; SST[nel] = %d + pass;
    for (u=1;u<=nE;u++) {
      if (E[u] != 0) { if (diff(E[u], var(vb)) != 0) {
        E[u] = subst(E[u], var(vb), ex);
        if (E[u] != 0) { cvs = cand(E[u], SC, %s, 99); CV[u] = cvs[1]; CS[u] = cvs[2]; } else { CV[u] = 0; }
      } }
    }
    TRK = subst(TRK, var(vb), ex);
    if (nel mod 50 == 0) { "  eliminated", nel, "; equations left", size(E), "; ms:", rtimer-T0; }
  }''' % (label, mode, (', division by ' + divide + ' (assumed invertible)') if divide else '', sc, ordname, sstoff, ordname)]
    if divide:
        S.append('''  ndiv = 0;
  for (u=1;u<=nE;u++) { if (E[u] != 0) { while (subst(E[u], %s, 0) == 0) { E[u] = E[u] / %s; ndiv++; progress = 1; } } }
  for (u=1;u<=nE;u++) { if (E[u] != 0) { cvs = cand(E[u], SC, %s, 99); CV[u] = cvs[1]; CS[u] = cvs[2]; } }
  "  pass", pass, ": eliminated so far", nel, "; equations left", size(E), "; divisions by %s in this pass:", ndiv;
  if (pass >= %d) { progress = 0; }''' % (divide, divide, ordname, divide, maxpass))
    S.append('}')
    S.append('"%s elimination finished: eliminated", nel, "; remaining equations", size(E), "; ms:", rtimer-T0;' % label)
    return '\n'.join(S)

def report(maxfact=80, tags=True):
    tagexpr = '", tag", EQDT[u], ' if tags else ''
    return '''// ---- result
intvec ISEL = 0:394; for (k=1;k<=nel;k++) { ISEL[SV[k]] = 1; }
intvec nelt = 0:4; intvec nfree = 0:4; string sfreeA = ""; string selA = "";
for (k=1;k<=394;k++) { if (ISEL[k]) { nelt[TYP[k]] = nelt[TYP[k]] + 1; if (TYP[k] > 1) { selA = selA + string(var(k)) + " "; } } else { nfree[TYP[k]] = nfree[TYP[k]] + 1; if (TYP[k] > 1) { sfreeA = sfreeA + string(var(k)) + " "; } } }
"eliminated unknowns by type (lambda, g, b, D):", nelt;
"eliminated A-unknowns:", selA;
"free unknowns by type (lambda, g, b, D):", nfree;
"free A-unknowns:", sfreeA;
string sfreeL = ""; for (k=1;k<=371;k++) { if (ISEL[k]==0) { sfreeL = sfreeL + string(var(k)) + " "; } } "free lambda-coefficients:", sfreeL;
int nleft = size(E); "remaining equations:", nleft;
int nz0 = 0; int nlin = 0; poly e0; poly tmpE; ideal lz = maxideal(1); for (k=1;k<=371;k++) { lz[k] = 0; }
map zl = R, lz;
for (u=1;u<=ncols(E);u++) { if (E[u] != 0) {
  if (deg(E[u], intvec(1:371)) > 1) { nlin++; }
  tmpE = E[u]; e0 = zl(tmpE);
  if (e0 != 0) { nz0++; }
} }
"remaining equations of degree > 1 in the lambda's:", nlin, "; with nonzero lambda-free part:", nz0;
int pr = 0;
for (u=1;u<=ncols(E);u++) { if (E[u] != 0) { pr++;
  "LEFT", pr, " (orig. order", EO[u], %s", size", size(E[u]), "):";
  if (size(E[u]) <= %d) { factorize(E[u]); } else { "  too long to factorize; leading terms:", lead(E[u]); }
} }''' % (tagexpr, maxfact)

def save_state(fname):
    return '''// ---- save state: remaining equations (with original order), substitutions, tracked values
write(":w %s", "// state written by the n1c0 scripts; execute(read(...)) in ring R");
write(":a %s", "ideal LEFT; intvec LEFTO; ideal SE0; intvec SV0; intvec SST0; ideal TRK0;");
int ns = 0;
for (u=1;u<=ncols(E);u++) { if (E[u] != 0) { ns++; write(":a %s", "LEFT[" + string(ns) + "]=" + string(E[u]) + "; LEFTO[" + string(ns) + "]=" + string(EO[u]) + ";"); } }
for (k=1;k<=nel;k++) { write(":a %s", "SV0[" + string(k) + "]=" + string(SV[k]) + "; SST0[" + string(k) + "]=" + string(SST[k]) + "; SE0[" + string(k) + "]=" + string(SE[k]) + ";"); }
for (k=1;k<=ncols(TRK);k++) { write(":a %s", "TRK0[" + string(k) + "]=" + string(TRK[k]) + ";"); }
"state saved to %s: remaining", ns, "; substitutions", nel;''' % ((fname,) * 6)

def load_state(fname, extra=()):
    """load a saved state; extra: additional equations (strings) appended with order marker 99"""
    S = ['execute(read("%s"));' % fname,
         '"state %s loaded: remaining", ncols(LEFT), "; substitutions", size(SV0);' % fname,
         'ideal E = LEFT; intvec EO = LEFTO;']
    for n, e in enumerate(extra, 1):
        S.append('E[ncols(LEFT)+%d] = %s; EO[ncols(LEFT)+%d] = 99;' % (n, e, n))
    S.append('intvec SV = SV0; ideal SE = SE0; intvec SST = SST0; int nel = size(SV0); TRK = TRK0;')
    S.append('"equations now:", size(E);')
    return '\n'.join(S)

def print_trk():
    return '''// ---- tracked values (expressions in the free unknowns, modulo the remaining equations)
"tracked list consistent:", ncols(TRK) == size(TRKN);
for (k=1;k<=ncols(TRK);k++) { if (size(TRK[k]) <= 60) { TRKN[k], "=", TRK[k]; } else { TRKN[k], "= (", size(TRK[k]), "terms ) lead", lead(TRK[k]); } }'''

def tail():
    return '"total ms:", rtimer-T0, "; memory (KB) in use / allocated:", memory(1) div 1024, memory(2) div 1024;\nquit;\n'

# ---------------- scripts ----------------
def s1():
    tb, ntrk = trk_block()
    return '\n'.join([hdr('s1_orderwise: build the rank-condition system, counts, delays, order-by-order elimination'),
                      build(delays=True), tb, elim_orderwise(), report(), print_trk(), save_state('state_s1.txt'), tail()])

def s2():
    tb, ntrk = trk_block()
    return '\n'.join([hdr('s2_global: global constant-pivot elimination from scratch, lambda first by index'),
                      build(delays=False), tb,
                      'ideal E = EQD; intvec EO = EQDO; intvec SV; ideal SE; intvec SST; int nel = 0;',
                      elim_global('lambda_index', None, label='global'), report(), print_trk(), save_state('state_s2.txt'), tail()])


def purey_forms():
    """Singular code defining the binary forms lambda_{ij,0}(0,0,y2,y3) from the tracked values TRK[34..54]"""
    from collections import OrderedDict
    base = 33   # TRK[1..23] = A-unknowns, TRK[24..32] = PUREW, TRK[33..53] = PUREY  (1-based: 24+9=33)
    forms = OrderedDict()
    for n, (ix, i, j, ex) in enumerate(PUREY):
        forms.setdefault((i, j), []).append('TRK[%d]*%s' % (base + n, mstr(ex)))
    S = ['ideal FY; list FYN; int nfy = 0;']
    for (i, j), terms in forms.items():
        S.append('nfy++; FY[nfy] = %s; FYN[nfy] = "lambda_%d%d,0(0,0,y2,y3)";' % (' + '.join(terms), i, j))
    return '\n'.join(S)

def s3():
    return '\n'.join([hdr('s3_rabinowitsch: is (b8-b3)(b8+b3) in the radical of the rank-condition ideal?'),
        'execute(read("state_s1.txt"));',
        '"loaded state_s1.txt: remaining equations", ncols(LEFT);',
        'ideal VL = variables(LEFT); string sva = ""; int u; for (u=1;u<=ncols(VL);u++) { if (rvar(VL[u]) > 371) { sva = sva + string(VL[u]) + " "; } } "A-unknowns occurring in the remaining equations:", sva;',
        'ideal JR = LEFT, 1 - t*(b(8)-b(3))*(b(8)+b(3));',
        'int T1 = rtimer; ideal SJR = std(JR);',
        '"1 in (remaining equations, 1 - t*(b8-b3)*(b8+b3)):", reduce(poly(1), SJR) == 0, "; std ms:", rtimer-T1;',
        tail()])

def branch_common(title, extra, exclude_b8, divide, sstoff):
    tb, ntrk = trk_block()
    S = [hdr(title), build(delays=False), tb, load_state('state_s1.txt', extra=extra), 'intvec ORDX = ORD;']
    if exclude_b8:
        S.append('ORDX[389] = 999;   // b8 is the localizing element: never used as a pivot')
    S.append(elim_global('lambda_index', divide, label=title.split(':')[0], ordname='ORDX', sstoff=sstoff))
    S.append(report(tags=False))
    S.append(print_trk())
    S.append('ideal J; int nj = 0; for (u=1;u<=ncols(E);u++) { if (E[u] != 0) { nj++; J[nj] = E[u]; } }')
    return S

def s4():
    S = branch_common('s4_branch_i: branch b3 = b8 = 0, test of the sigma-fixed point (1:0:...:0) in F2',
                      ['b(3)', 'b(8)'], False, None, 200)
    S.append(r"""// ---- F2 point (1:0:...:0): conditions b7 = b10 = 0 and pure w-power coefficient of every lambda_{ij,0} = 0
ideal TV = TRK[17], TRK[20], TRK[24..32]; list TVN = TRKN[17], TRKN[20], TRKN[24..32];
for (k=1;k<=ncols(TV);k++) { TVN[k], "=", TV[k]; }
"all F2 conditions identically 0 after substitution:", size(TV) == 0;
if (size(TV) > 0) {
  if (nj > 0) { int T1 = rtimer; ideal SJ = std(J); "std of remaining equations: ms", rtimer-T1, "; 1 in it:", reduce(poly(1), SJ) == 0;
    ideal NTV = reduce(TV, SJ); "normal forms of the F2 conditions modulo the branch ideal:"; NTV;
    "all F2 conditions in the branch ideal:", size(NTV) == 0; }
}""")
    S.append(tail())
    return '\n'.join(S)

def s5():
    S = branch_common('s5_branch_ii: branch b3 = b8, b8 invertible; test of the sigma-fixed points (0:0:+-sqrt(D):1|0) in F3',
                      ['b(3)-b(8)'], True, 'b(8)', 300)
    S.append(purey_forms())
    S.append(r"""// ---- F3 points p = (0:0:y2:y3|0) with y2^2 = D*y3^2 (Q(p) = 0): conditions lambda_{ij,0}(p) = 0
for (k=1;k<=nfy;k++) { FYN[k], "=", FY[k]; }
ideal K = J, 1 - t*b(8), y2^2 - D*y3^2;
int T1 = rtimer; ideal SK = std(K); "std of (remaining equations, 1 - t*b8, y2^2 - D*y3^2): ms", rtimer-T1, "; 1 in it:", reduce(poly(1), SK) == 0;
ideal NFY = reduce(FY, SK);
for (k=1;k<=nfy;k++) { "normal form of", FYN[k], "modulo (branch ideal, y2^2-D*y3^2):", NFY[k]; }
"all F3 conditions hold on the branch:", size(NFY) == 0;
ideal KJ = J, 1 - t*b(8); ideal SKJ = std(KJ);
"b11 modulo (remaining equations, 1 - t*b8):", reduce(b(11), SKJ), "; b12:", reduce(b(12), SKJ);
"1 in (remaining equations, 1 - t*b8):", reduce(poly(1), SKJ) == 0;""")
    S.append(tail())
    return '\n'.join(S)

def lamp_matrix():
    rows = []
    for i in range(1, 6):
        for j in range(i, 6):
            rows.append(', '.join(lam_poly(i, j, k) for k in range(6)))
    return 'matrix LAMP[15][6] = ' + ',\n  '.join(rows) + ';'

def syz_block():
    """exact rank of the trivial-syzygy map at the point (matrix Apt)"""
    S = ['proc cf(poly P, poly mon) { if (P == 0) { return(poly(0)); } matrix CC = coef(P, w*y1*y2*y3); int q; for (q=1;q<=ncols(CC);q++) { if (CC[1,q] == mon) { return(CC[2,q]); } } return(poly(0)); }',
         'matrix SYZM[94][371]; poly Pk; poly chk; int nbad = 0;']
    row = 0
    for i in range(1, 6):
        for j in range(i, 6):
            for m in range(1, 6):
                n = c[i] + c[j] - r[m]
                if n < 0:
                    continue
                for mu in mons(n, -dl[i] * dl[j] * dl[m]):
                    row += 1
                    for k in range(6):
                        S.append('Pk = %s*Apt[%d,%d]; chk = Pk;' % (mstr(mu), m + 1, k + 1))
                        for ix, ex in LAM[(i, j, k)]:
                            S.append('SYZM[%d,%d] = cf(Pk, %s); chk = chk - SYZM[%d,%d]*%s;' % (row, ix, mstr(ex), row, ix, mstr(ex)))
                        S.append('if (chk != 0) { nbad++; }')
    assert row == 94
    S.append('"trivial-syzygy directions leaving the lambda-ansatz (should be 0):", nbad;')
    S.append('"exact rank over Q of the trivial-syzygy map (94 directions -> 371 lambda-coordinates):", rank(SYZM);')
    return '\n'.join(S)

def s6(seed1=17, seed2=4711):
    S = branch_common('s6_branch_iii: branch b3 = -b8, b8 invertible (expected to survive)',
                      ['b(3)+b(8)'], True, 'b(8)', 400)
    S.append(r"""// ---- graph check: the eliminated A-unknowns are polynomials in the free A-unknowns (no lambda occurs)
int ngl = 0; for (k=1;k<=23;k++) { if (deg(TRK[k], intvec(1:371)) > 0) { ngl++; "A-unknown depending on a free lambda:", TRKN[k]; } }
"graph check: A-unknowns whose expression contains a free lambda (should be 0):", ngl;
"remaining equations (should be 0):", nj;
write(":w expr_iii.txt", "// branch (iii) b3 = -b8, b8 != 0: eliminated A-unknowns as polynomials in the free A-unknowns (s6_branch_iii)");
for (k=1;k<=23;k++) { write(":a expr_iii.txt", TRKN[k] + " = " + string(TRK[k])); }
"expressions written to expr_iii.txt";""")
    S.append(purey_forms())
    S.append(r"""// ---- information: sigma-fixed-point conditions on branch (iii) (the branch is excluded only if they hold identically)
"F2 point (1:0:...:0): b7 =", TRK[17], "; b10 =", TRK[20];
for (k=24;k<=32;k++) { TRKN[k], "=", TRK[k]; }
ideal SKQ = std(ideal(y2^2 - D*y3^2)); ideal NFYQ = reduce(FY, SKQ);
for (k=1;k<=nfy;k++) { "F3 points: normal form of", FYN[k], "modulo y2^2 - D*y3^2:", NFYQ[k]; }""")
    S.append(lamp_matrix())
    S.append(r"""// ---- random rational points of the branch and checks at them
LIB "matrix.lib";   // for submat
proc randpt(intvec ISEL, ideal SE, intvec SV, int nel, int seed)
{
  // free unknowns: random nonzero integers in [-9,9]; eliminated ones by back-substitution in reverse order
  system("random", seed);
  ideal VAL = maxideal(1); int k; poly sv;
  for (k=1;k<=394;k++) { if (ISEL[k]==0) { VAL[k] = random(1,9) * (2*random(0,1) - 1); } }
  for (k=nel; k>=1; k--) { map ev = basering, VAL; sv = SE[k]; VAL[SV[k]] = ev(sv); kill ev; }
  return(VAL);
}
proc checkpt(ideal VAL, matrix A, matrix LAMP)
{
  // the 15 rank conditions E_ij = B_ij - sum_k lambda_ij,k B_0k as polynomials in (w,y1,y2,y3) at the point
  map ev = basering, VAL; matrix Apt = ev(A); matrix Lpt = ev(LAMP); int i, j, k, p, nbad; poly pol;
  matrix B0p[1][6]; for (k=1;k<=6;k++) { B0p[1,k] = cofac(Apt,k,1); }
  p = 0; nbad = 0;
  for (i=1;i<=5;i++) { for (j=i;j<=5;j++) { p++;
    pol = cofac(Apt, j+1, i+1); for (k=1;k<=6;k++) { pol = pol - Lpt[p,k]*B0p[1,k]; }
    if (pol != 0) { nbad++; } } }
  return(nbad);
}
proc jacrank(ideal EQD, ideal VAL, int nd)
{
  // ranks mod 32003 of the Jacobian of all coefficient equations w.r.t. the 394 unknowns at the point
  map evp = basering, VAL;
  matrix JQ[nd+1][394]; ideal dcol; ideal dcv; int k, u; int bad = 0;
  for (k=1;k<=394;k++) { dcol = diff(EQD, var(k)); dcv = evp(dcol); if (ncols(dcv) != nd) { bad++; } for (u=1;u<=nd;u++) { JQ[u,k] = dcv[u]; } }
  JQ[nd+1,384] = 1; JQ[nd+1,389] = 1;   // last row: the branch equation b3 + b8
  def Rsave = basering;
  ring Rp = 32003, (x), dp;
  matrix Jp = imap(Rsave, JQ);
  matrix J0 = submat(Jp, 1..nd, 1..394); matrix J0l = submat(Jp, 1..nd, 1..371);
  matrix J1l = submat(Jp, 1..nd+1, 1..371);
  intvec res = rank(J0), rank(J0l), rank(Jp), rank(J1l), bad;
  setring Rsave;
  return(res);
}
ideal VAL1 = randpt(ISEL, SE, SV, nel, """ + str(seed1) + r""");
ideal VAL2 = randpt(ISEL, SE, SV, nel, """ + str(seed2) + r""");
int nconst1 = 0; int nconst2 = 0; for (k=1;k<=394;k++) { if (deg(VAL1[k]) > 0) { nconst1++; } if (deg(VAL2[k]) > 0) { nconst2++; } }
"unknowns not evaluated to constants at points 1, 2 (should be 0):", nconst1, nconst2;
string spt = ""; for (k=372;k<=394;k++) { spt = spt + string(var(k)) + "=" + string(VAL1[k]) + " "; } "point 1, A-coefficients:", spt;
spt = ""; for (k=372;k<=394;k++) { spt = spt + string(var(k)) + "=" + string(VAL2[k]) + " "; } "point 2, A-coefficients:", spt;
"point 1: rank conditions E_ij violated (as polynomials in w,y; should be 0):", checkpt(VAL1, A, LAMP);
"point 2: rank conditions E_ij violated (as polynomials in w,y; should be 0):", checkpt(VAL2, A, LAMP);
map evp1 = R, VAL1; matrix Apt = evp1(A);
"point 1: det A at the point is nonzero:", det(Apt) != 0;
""" + syz_block() + r"""
intvec jr1 = jacrank(EQD, VAL1, nd);
"point 1 mod 32003: Jacobian rank", jr1[1], "; nullity", 394 - jr1[1], "; projection of the kernel to the 23 A-coordinates: rank", 23 - jr1[1] + jr1[2], "; size mismatches", jr1[5];
"point 1 mod 32003, with b3 + b8 added: rank", jr1[3], "; nullity", 394 - jr1[3], "; projection rank", 23 - jr1[3] + jr1[4];
intvec jr2 = jacrank(EQD, VAL2, nd);
"point 2 mod 32003: Jacobian rank", jr2[1], "; nullity", 394 - jr2[1], "; projection of the kernel to the 23 A-coordinates: rank", 23 - jr2[1] + jr2[2], "; size mismatches", jr2[5];
"point 2 mod 32003, with b3 + b8 added: rank", jr2[3], "; nullity", 394 - jr2[3], "; projection rank", 23 - jr2[3] + jr2[4];""")
    S.append(tail())
    return '\n'.join(S)


# ---------------- torus (case_common section 4) ----------------
# log-parameters: (ww, w1, w2, w3, g0, g1, g2, g3, g4, g5): w -> e^ww w, y_i -> e^wi y_i, A -> P^t A P with P = diag(e^g0..e^g5)
PNAMES = ['ww', 'w1', 'w2', 'w3', 'g0', 'g1', 'g2', 'g3', 'g4', 'g5']
def lf(**kw):
    v = [0] * 10
    for k, x in kw.items():
        v[PNAMES.index(k)] += x
    return v
def addv(*vs):
    return [sum(t) for t in zip(*vs)]
def gam(k):
    return lf(**{'g%d' % k: 1})
def mon_w(a, b1, e2, f3):
    return lf(ww=a, w1=b1, w2=e2, w3=f3)
def unknown_weights():
    W = []
    for u in UNK:
        if u['type'] == 'lambda':
            ix = u['idx']
            (i, j, k), ex = next(((key, ex) for key, lst in LAM.items() for (jx, ex) in lst if jx == ix))
            v = addv(gam(0), gam(k), [-x for x in gam(i)], [-x for x in gam(j)], mon_w(*ex))
        elif u['type'] == 'g':
            n = int(u['name'][2:-1]); s_, ex = Gmons[n - 1]
            v = addv(gam(0), gam(0), mon_w(2 + 2 * ex[0], ex[1], ex[2], ex[3]))
        elif u['type'] == 'b':
            n = int(u['name'][2:-1]); q, sm, k0 = Bmons[n - 1]
            qi = int(q[1]); ex = [0, 0, 0, 0]
            for tt in sm.split('*'):
                if '^' in tt:
                    v_, pw = tt.split('^'); ex[int(v_[1])] += int(pw)
                else:
                    ex[int(tt[1])] += 1
            v = addv(gam(0), gam(qi), mon_w(1 + 2 * ex[0], ex[1], ex[2], ex[3]))
        else:   # D: coefficient of y3 in A_11
            v = addv(gam(1), gam(1), mon_w(0, 0, 0, 1))
        W.append(v)
    return W
def constraints():
    C = []
    C.append(('A_12 = y2', addv(gam(1), gam(2), mon_w(0, 0, 1, 0))))
    C.append(('A_13 = y1', addv(gam(1), gam(3), mon_w(0, 1, 0, 0))))
    C.append(('A_22 = y3', addv(gam(2), gam(2), mon_w(0, 0, 0, 1))))
    C.append(('A_24 = y1', addv(gam(2), gam(4), mon_w(0, 1, 0, 0))))
    C.append(('A_33 = -y3', addv(gam(3), gam(3), mon_w(0, 0, 0, 1))))
    C.append(('A_34 = y2', addv(gam(3), gam(4), mon_w(0, 0, 1, 0))))
    C.append(('A_15 = w', addv(gam(1), gam(5), mon_w(1, 0, 0, 0))))
    C.append(('A_05: coefficient 1 of y2^2', addv(gam(0), gam(5), mon_w(0, 0, 2, 0))))
    C.append(('A_05: coefficient -1 of y1^2', addv(gam(0), gam(5), mon_w(0, 2, 0, 0))))
    dA11 = addv(gam(1), gam(1), mon_w(0, 0, 0, 1)); dA44 = addv(gam(4), gam(4), mon_w(0, 0, 0, 1)); dA05 = addv(gam(0), gam(5), mon_w(0, 0, 0, 2))
    C.append(('D in A_44 moves like D in A_11', [a - b for a, b in zip(dA44, dA11)]))
    C.append(('D in A_05 (y3^2) moves like D in A_11', [a - b for a, b in zip(dA05, dA11)]))
    return C

def s7(free_iii=('g(5)', 'b(2)', 'b(5)', 'b(6)', 'b(8)', 'b(9)', 'b(11)', 'b(12)', 'D')):
    W = unknown_weights(); C = constraints()
    grading = lf(ww=1, w1=2, w2=2, w3=2, g0=-4, g1=-1, g2=-1, g3=-1, g4=-1, g5=0)
    sdir = lf(ww=1, g0=1, g5=-1)
    udir = lf(ww=-1, w3=-2, g0=-2, g1=-1, g2=1, g3=1, g4=-1, g5=2)   # twice the one-parameter group u of case_n1c1 item 7
    def dot(a, b): return sum(x * y for x, y in zip(a, b))
    Ws = [dot(v, sdir) for v in W]; Wu = [dot(v, udir) for v in W]; Wg = [dot(v, grading) for v in W]
    S = [hdr('s7_torus: diagonal symmetries of the normal form N1(D), c = 0; weights; semi-invariance of all equations; effective dimension')]
    S.append('LIB "matrix.lib";   // concat, submat')
    S.append('proc iszero(matrix M) { int i, j; for (i=1;i<=nrows(M);i++) { for (j=1;j<=ncols(M);j++) { if (M[i,j] != 0) { return(0); } } } return(1); }')
    S.append('ring RT = 0, (x), dp;')
    S.append('matrix CN[%d][10] = %s;' % (len(C), ', '.join(str(x) for _, v in C for x in v)))
    S.append('"constraints (rows) on the log-parameters (ww,w1,w2,w3,g0,...,g5):"; ' + ' '.join('"  %s";' % n for n, _ in C))
    S.append('module KER = syz(CN); "dimension of the torus (kernel of the constraint matrix):", ncols(KER);')
    S.append('matrix GR[10][1] = %s; matrix SD[10][1] = %s; matrix UD[10][1] = %s;' % (', '.join(map(str, grading)), ', '.join(map(str, sdir)), ', '.join(map(str, udir))))
    S.append('"grading, s, u satisfy the constraints:", iszero(CN*GR), iszero(CN*SD), iszero(CN*UD);')
    S.append('matrix GSU = concat(GR, SD, UD); "rank of (grading, s, u):", rank(GSU);')
    S.append('matrix WL[394][10] = %s;' % ', '.join(str(x) for v in W for x in v))
    S.append('matrix KM = KER; matrix WK = WL*KM; "rank of the weight matrix on all 394 unknowns (kernel coordinates):", rank(WK);')
    S.append('"grading direction acts trivially on every unknown:", iszero(WL*GR);')
    S.append('matrix KG = concat(KM, GR); "grading direction lies in the kernel (rank of kernel + grading = 3):", rank(KG);')
    # effective dimension on free coordinates of branch (iii)
    rows = [next(u['idx'] for u in UNK if u['name'] == nm) for nm in free_iii]
    S.append('intvec FR = %s;' % ','.join(map(str, rows)))
    S.append('matrix WF = submat(WK, FR, 1..ncols(WK)); "effective dimension e = rank of the weights of the free A-coordinates of branch (iii) (%s):", rank(WF);' % ' '.join(free_iii))
    Arows = list(range(372, 395))
    S.append('matrix WA = submat(WK, %s, 1..ncols(WK)); "rank of the weights of all 23 A-coordinates:", rank(WA);' % ('intvec(' + ','.join(map(str, Arows)) + ')'))
    # weights table under s and u
    S.append('setring R;')
    S.append('"weights of the A-coordinates under (s, u) [u = twice the group of case_n1c1 item 7]:";')
    for u in UNK:
        if u['type'] != 'lambda':
            S.append('"  %s (%s): s-weight %d, u-weight %d";' % (u['name'], u['mon'], Ws[u['idx'] - 1], Wu[u['idx'] - 1]))
    S.append('"weights of the lambda-coefficients under (s, u): see the table in the output below";')
    for u in UNK:
        if u['type'] == 'lambda':
            S.append('"  %s %s: %d %d";' % (u['name'], u['mon'], Ws[u['idx'] - 1], Wu[u['idx'] - 1]))
    S.append('"grading weights nonzero on some unknown (should be 0 unknowns):", %d;' % sum(1 for x in Wg if x != 0))
    # semi-invariance of every coefficient equation
    S.append(build(delays=False))
    S.append('intvec WS = %s, 0,0,0,0,0;' % ','.join(map(str, Ws)))
    S.append('intvec WU = %s, 0,0,0,0,0;' % ','.join(map(str, Wu)))
    S.append('intvec tw1 = -1,3,0; intvec tw2 = 1,-3,0; "self-test of deg with negative weights (expect 1 and 0):", deg(l(1)^2*l(2) + l(3), tw1), deg(l(1)^2*l(2) + l(3), tw2);')
    S.append('intvec mWS = -WS; intvec mWU = -WU; int nhs = 0; int nhu = 0;')
    S.append('for (u=1;u<=nd;u++) { if (deg(EQD[u], WS) != -deg(EQD[u], mWS)) { nhs++; } if (deg(EQD[u], WU) != -deg(EQD[u], mWU)) { nhu++; } }')
    S.append('"coefficient equations that are not homogeneous for the s-weights / u-weights (should be 0, 0):", nhs, nhu;')
    S.append(tail())
    return '\n'.join(S)


def s9():
    """check that every nonzero coefficient of every E_ij sits at a monomial of weight 8+c_i+c_j and sign -d_i d_j"""
    S = [hdr('s9_slots: every monomial of E_ij has weight 8+c_i+c_j and sigma-sign -delta_i delta_j (sign rule check)'), build(delays=False)]
    S.append('int nbadw = 0; int nbads = 0; int ntot = 0; intvec ex; int a; int wt; int sg;')
    for i in range(1, 6):
        for j in range(i, 6):
            terms = ' + '.join('(%s)*B0[1,%d]' % (lam_poly(i, j, k), k + 1) for k in range(6))
            S.append('pol = cofac(A,%d,%d) - (%s);' % (j + 1, i + 1, terms))
            S.append('if (pol != 0) { C = coef(pol, w*y1*y2*y3); for (u=1;u<=ncols(C);u++){ if (C[2,u] != 0) { ntot++; ex = leadexp(C[1,u]); wt = ex[395] + 2*(ex[396]+ex[397]+ex[398]); sg = (-1)^(ex[395]+ex[397]+ex[398]); if (wt != %d) { nbadw++; } if (sg != %d) { nbads++; } } } }' % (8 + c[i] + c[j], -dl[i] * dl[j]))
    S.append('"nonzero coefficients checked:", ntot, "; at a monomial of wrong weight:", nbadw, "; of wrong sigma-sign:", nbads;')
    # also the entries of A: weight and sign of every monomial of every entry
    S.append('intvec cdeg = 0,3,3,3,3,4; intvec dsg = 1,-1,-1,1,1,-1; int jj, kk; poly ent; int nbe = 0; int nbs = 0;')
    S.append('for (jj=1;jj<=6;jj++) { for (kk=1;kk<=6;kk++) { ent = A[jj,kk]; if (ent != 0) { C = coef(ent, w*y1*y2*y3); for (u=1;u<=ncols(C);u++) { ex = leadexp(C[1,u]); wt = ex[395] + 2*(ex[396]+ex[397]+ex[398]); sg = (-1)^(ex[395]+ex[397]+ex[398]); if (wt != 8 - cdeg[jj] - cdeg[kk]) { nbe++; } if (sg != -dsg[jj]*dsg[kk]) { nbs++; } } } } }')
    S.append('"entries of A: monomials of wrong weight:", nbe, "; of wrong sigma-sign (rule sign(A_jk) = -delta_j delta_k):", nbs;')
    S.append(tail())
    return '\n'.join(S)

SCRIPTS = ['s1', 's2', 's3', 's4', 's5', 's6', 's7', 's9']


if __name__ == '__main__':
    txt, eqo = table()
    open(os.path.join(HERE, 'unknowns_table.txt'), 'w').write(txt)
    which = sys.argv[1:] or SCRIPTS
    for w in which:
        open(os.path.join(HERE, '%s.sing' % w), 'w').write(globals()[w]())
    print('written:', which)
