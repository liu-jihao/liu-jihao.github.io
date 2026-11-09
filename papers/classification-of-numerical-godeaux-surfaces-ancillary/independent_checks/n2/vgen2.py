#!/usr/bin/env python3
"""Independent verifier for case n2, version 2 (written from plan/case_n2.md, plan/case_common.md and
CCO94 7.3-7.5 read from the page image p. 123, without using the case's scripts).
Generates vn2b_c1.sing and vn2b_c0.sing:
  - ring RS: the N2 normal form A, its adjugate B (cofactors), checks (A*B = det A * I, sigma signs, weights,
    order-0 data);
  - ring RB: rank-condition equations E_ij = B_ij - sum_k lambda_{ij,k} B_0k (1<=i<=j<=5), one per monomial;
  - my own global constant-pivot elimination: equations with fewest terms first; pivot variable = lambda before A,
    then the LARGEST ring index (my ring order x(1..371), g10..g1, b12..b1, y3, y2, y1, w; my lambda numbering);
  - c = 1: graph check, substitution check, random rational point, exact Jacobian ranks over Q, membership test
    B_ij in (B_00..B_05) in Q[w,y] (no lambda), perturbation tests, tau2 checks;
  - c = 0: images of the five deciding coefficients at p = (0:0:1:0|0).
Usage: python3 vgen2.py
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
cdeg = [0, 3, 3, 3, 3, 4]
delta = [1, -1, -1, 1, 1, -1]


def monos(weight, sign):
    """all monomials w^a y1^b y2^e y3^f of weight a+2(b+e+f) and sigma-sign (-1)^(a+e+f) = sign"""
    res = []
    if weight < 0:
        return res
    for a in range(weight + 1):
        rest = weight - a
        if rest % 2:
            continue
        d = rest // 2
        for b in range(d + 1):
            for e in range(d - b + 1):
                f = d - b - e
                if (-1) ** (a + e + f) == sign:
                    res.append((a, b, e, f))
    res.sort(key=lambda m: (-m[0], -m[3], -m[2], -m[1]))
    return res


def mstr(m):
    a, b, e, f = m
    parts = []
    for v, ex in (('w', a), ('y1', b), ('y2', e), ('y3', f)):
        if ex == 1:
            parts.append(v)
        elif ex > 1:
            parts.append('%s^%d' % (v, ex))
    return '*'.join(parts) if parts else '1'


pairs = [(i, j) for i in range(1, 6) for j in range(i, 6)]
pairs_rev = list(reversed(pairs))
lam = []
for (i, j) in pairs_rev:
    for k in range(5, -1, -1):
        for m in monos(cdeg[i] + cdeg[j] - cdeg[k], delta[i] * delta[j] * delta[k]):
            lam.append((len(lam) + 1, i, j, k, m))
NL = len(lam)
decide = [(idx, i, j, k, m) for (idx, i, j, k, m) in lam if k == 0 and m[0] == 0 and m[1] == 0 and m[3] == 0]


# lambda table and counts (my numbering; identical to the numbering used in the generated scripts)
nsyz = 0
for (i, j) in pairs:
    for mm in range(1, 6):
        nsyz += len(monos(cdeg[i] + cdeg[j] - (8 - cdeg[mm]), -delta[i] * delta[j] * delta[mm]))
slots = {}
for (i, j) in pairs:
    for m in monos(8 + cdeg[i] + cdeg[j], -delta[i] * delta[j]):
        slots[m[0]] = slots.get(m[0], 0) + 1
bylam_k, bylam_w = {}, {}
for (_, i, j, k, m) in lam:
    bylam_k[k] = bylam_k.get(k, 0) + 1
    bylam_w[m[0]] = bylam_w.get(m[0], 0) + 1
with open(os.path.join(HERE, 'vn2b_lambda_table.txt'), 'w') as fh:
    fh.write('# my lambda numbering x(idx): idx i j k monomial w-order\n')
    fh.write('# total %d ; by k %s ; by w-order %s\n' % (NL, dict(sorted(bylam_k.items())), dict(sorted(bylam_w.items()))))
    fh.write('# trivial-syzygy coefficient count: %d\n' % nsyz)
    fh.write('# RC monomial slots: %d ; by w-order %s\n' % (sum(slots.values()), dict(sorted(slots.items()))))
    fh.write('# deciding coefficients at (0:0:1:0|0): %s\n' % ', '.join('x(%d)=lambda_{%d%d,0}[%s]' % (d[0], d[1], d[2], mstr(d[4])) for d in decide))
    for (idx, i, j, k, m) in lam:
        fh.write('%d %d %d %d %s %d\n' % (idx, i, j, k, mstr(m), m[0]))

Gvars = ['g%d' % t for t in range(10, 0, -1)]
Bvars = ['b%d' % t for t in range(12, 0, -1)]
Avars = Gvars + Bvars
small_vars = ','.join(Avars + ['y3', 'y2', 'y1', 'w'])
big_vars = 'x(1..%d),' % NL + small_vars


def lam_expr(i, j, k):
    terms = ['x(%d)*%s' % (idx, mstr(m)) if mstr(m) != '1' else 'x(%d)' % idx
             for (idx, ii, jj, kk, m) in lam if (ii, jj, kk) == (i, j, k)]
    return '(' + ' + '.join(terms) + ')' if terms else '0'


def tau2_weights():
    """tau2: w -> w, y1 -> t y1, y2 -> t^2 y2, y3 -> y3, A -> P A(phi) P, P = diag(t^-3,t^-1,1,1,t^-1,t).
    weight of a coefficient = P_j + P_k + (tau2-degree of its monomial), y1:1, y2:2, y3:0, y0:0."""
    P = [-3, -1, 0, 0, -1, 1]
    Gm = {1: (0, 0, 0, 3), 2: (0, 2, 0, 1), 3: (0, 0, 1, 2), 4: (0, 2, 1, 0), 5: (0, 0, 2, 1), 6: (0, 0, 3, 0),
          7: (1, 1, 0, 1), 8: (1, 1, 1, 0), 9: (2, 0, 0, 1), 10: (2, 0, 1, 0)}
    Bm = {1: (2, (1, 0, 1, 0)), 2: (2, (1, 0, 0, 1)), 3: (2, (0, 1, 1, 0)), 4: (2, (0, 1, 0, 1)),
          5: (3, (2, 0, 0, 0)), 6: (3, (1, 1, 0, 0)), 7: (3, (0, 2, 0, 0)), 8: (3, (0, 0, 2, 0)),
          9: (3, (0, 0, 1, 1)), 10: (3, (0, 0, 0, 2)), 11: (4, (2, 0, 0, 0)), 12: (4, (0, 0, 0, 2))}
    wt = {}
    for m, (e0, e1, e2, e3) in Gm.items():
        wt['g%d' % m] = 2 * P[0] + e1 + 2 * e2
    for m, (k, (e0, e1, e2, e3)) in Bm.items():
        wt['b%d' % m] = P[0] + P[k] + e1 + 2 * e2
    return wt


PROCS = r"""
proc cofac(matrix M, int i, int j)
{
  int n = nrows(M);
  matrix S[n-1][n-1];
  int r, s, rr, ss;
  rr = 0;
  for (r = 1; r <= n; r++) {
    if (r != i) {
      rr++; ss = 0;
      for (s = 1; s <= n; s++) {
        if (s != j) { ss++; S[rr, ss] = M[r, s]; }
      }
    }
  }
  int sg = 1;
  if (((i + j) mod 2) == 1) { sg = -1; }
  return(sg * det(S));
}
proc iszero(matrix M)
{
  return(size(ideal(M)) == 0);
}
"""

MATRIX = r"""
poly y0 = w^2;
poly Q = y2*y3 - y1^2;
poly G = g1*y3^3 + g2*y1^2*y3 + g3*y2*y3^2 + g4*y1^2*y2 + g5*y2^2*y3 + g6*y2^3 + g7*y0*y1*y3 + g8*y0*y1*y2 + g9*y0^2*y3 + g10*y0^2*y2;
poly q1 = 0;
poly q2 = b1*y0*y2 + b2*y0*y3 + b3*y1*y2 + b4*y1*y3;
poly q3 = b5*y0^2 + b6*y0*y1 + b7*y1^2 + b8*y2^2 + b9*y2*y3 + b10*y3^2;
poly q4 = b11*y0^2 + b12*y3^2;
int cc = @C@;
matrix A[6][6] =
  y0*G, w*q1, w*q2, w*q3, w*q4, Q,
  w*q1, y2,   0,    y1,   0,    w,
  w*q2, 0,    y3,   cc*y0, y1,  0,
  w*q3, y1,   cc*y0, y3,  0,    0,
  w*q4, 0,    y1,   0,    y2,   0,
  Q,    w,    0,    0,    0,    0;
"""

CHECKS = r"""
int i, j, k, t, u, v;
matrix B[6][6];
for (i = 1; i <= 6; i++) { for (j = 1; j <= 6; j++) { B[i, j] = cofac(A, i, j); } }
poly dA = det(A);
matrix I6[6][6]; for (i = 1; i <= 6; i++) { I6[i, i] = 1; }
"CHECK A symmetric: " + string(iszero(A - transpose(A)));
"CHECK cofactor matrix symmetric: " + string(iszero(B - transpose(B)));
"CHECK A*B = det(A)*I: " + string(iszero(A*B - dA*I6));
ideal sv = maxideal(1);
sv[nvars(RS)-3] = -y3; sv[nvars(RS)-2] = -y2; sv[nvars(RS)] = -w;
map sig = RS, sv;
intvec dl = 1, -1, -1, 1, 1, -1;
int badA = 0; int badB = 0; poly tA; poly tB;
for (i = 1; i <= 6; i++) { for (j = 1; j <= 6; j++) {
  tA = A[i, j]; tB = B[i, j];
  if (sig(tA) != -dl[i]*dl[j]*tA) { badA++; }
  if (sig(tB) != -dl[i]*dl[j]*tB) { badB++; }
} }
"CHECK sigma(A_jk) = -d_j d_k A_jk, failures: " + string(badA);
"CHECK sigma(B_ij) = -d_i d_j B_ij, failures: " + string(badB);
"CHECK det A sigma-invariant: " + string(sig(dA) == dA);
intvec cw = 0, 3, 3, 3, 3, 4;
intvec wv = 0:(nvars(RS)-4), 2, 2, 2, 1;
int badW = 0;
for (i = 1; i <= 6; i++) { for (j = 1; j <= 6; j++) {
  if (A[i, j] != 0) {
    if (jet(A[i, j], 8 - cw[i] - cw[j], wv) - jet(A[i, j], 7 - cw[i] - cw[j], wv) != A[i, j]) { badW++; }
  }
} }
"CHECK weights of A entries, failures: " + string(badW);
matrix al[4][4];
for (i = 1; i <= 4; i++) { for (j = 1; j <= 4; j++) { al[i, j] = subst(A[i+1, j+1], w, 0); } }
"CHECK det(alpha) - Q^2 = " + string(det(al) - Q^2);
matrix adjal[4][4];
for (i = 1; i <= 4; i++) { for (j = 1; j <= 4; j++) { adjal[i, j] = cofac(al, i, j); } }
int qdiv = 1;
for (i = 1; i <= 4; i++) { for (j = 1; j <= 4; j++) { if (reduce(adjal[i, j], std(ideal(Q))) != 0) { qdiv = 0; } } }
"CHECK Q divides adj(alpha): " + string(qdiv);
"CHECK B_05 + Q^3 = 0 mod w: " + string(subst(B[1, 6], w, 0) + Q^3 == 0);
int okb = 1;
for (k = 2; k <= 5; k++) {
  if (subst(B[1, k], w, 0) != 0) { okb = 0; }
  if (subst(diff(B[1, k], w), w, 0) != Q*adjal[k-1, 1]) { okb = 0; }
}
"CHECK B_0k = 0 mod w and [w^1]B_0k = Q*adj(alpha)_{k1} (k=1..4): " + string(okb);
matrix aa[4][4];
for (i = 1; i <= 4; i++) { for (j = 1; j <= 4; j++) { aa[i, j] = A[i+1, j+1]; } }
"CHECK B_00 = -y0*adj(a)_11: " + string(B[1, 1] == -y0*cofac(aa, 1, 1));
"B_00 = " + string(B[1, 1]);
"""

BIGRING_HEAD = r"""
matrix B = imap(RS, B);
poly y0 = w^2;
int nv = nvars(RB);
int NU = nv - 4;
int NL = @NL@;
ideal EQ; intvec ORD = 0:1000;
int neq = 0; int badslot = 0;
poly E; matrix C; intvec ev;
"""

EQ_BLOCK = r"""if (E != 0) {
  C = coef(E, w*y1*y2*y3);
  for (t = 1; t <= ncols(C); t++) {
    neq++; EQ[neq] = C[2, t]; ev = leadexp(C[1, t]); ORD[neq] = ev[nv];
    if (ev[nv] + 2*(ev[nv-1] + ev[nv-2] + ev[nv-3]) != @WT@) { badslot++; }
    if (((ev[nv] + ev[nv-2] + ev[nv-3]) mod 2) != @SG@) { badslot++; }
  }
}"""

COUNTS = r"""
"number of nonzero coefficient equations: " + string(neq);
"monomials with wrong weight or sign: " + string(badslot);
intvec NZ = 0:16;
for (t = 1; t <= neq; t++) { NZ[ORD[t]+1] = NZ[ORD[t]+1] + 1; }
"nonzero equations by w-order 0..15: " + string(NZ);
int nconst = 0;
for (t = 1; t <= neq; t++) { if (deg(EQ[t]) == 0) { nconst++; } }
"constant equations: " + string(nconst);
ideal KEPT; intvec KORD = 0:1000; int nk = 0; int found; poly pn;
for (u = 0; u <= 15; u++) {
  for (t = 1; t <= neq; t++) {
    if (ORD[t] == u) {
      pn = EQ[t] / leadcoef(EQ[t]);
      found = 0;
      for (v = 1; v <= nk; v++) { if (leadmonom(KEPT[v]) == leadmonom(pn)) { if (KEPT[v] == pn) { found = 1; break; } } }
      if (!found) { nk++; KEPT[nk] = pn; KORD[nk] = u; }
    }
  }
}
intvec ND = 0:16;
for (t = 1; t <= nk; t++) { ND[KORD[t]+1] = ND[KORD[t]+1] + 1; }
"distinct up to scalar: " + string(nk);
"distinct by w-order 0..15: " + string(ND);
ideal EQ0 = KEPT;
"""

ELIM = r"""
// ---------------- my constant-pivot elimination, stage by stage in the w-order of the equations ----------------
// Stage m adds the (deduplicated) equations of w-order m; within a stage, sweeps look only at equations with at most
// thr terms (thr reset to 2 after a sweep with a pivot, else thr grows), so short equations are used first.
// Pivot variable of an equation: it occurs only in one degree-1 term, with constant coefficient;
// priority: lambda (index <= NL) before A-unknowns; within a type, the LARGEST ring index (my numbering).
ideal EQW;
ideal X = maxideal(1);
intvec ELI = 0:nv;  int npiv = 0;
int best; int ok; int s; poly f; poly lin; poly hi; ideal vl; ideal vh; poly cf; poly ex;
int incons = 0; int thr; int maxsz; int progress; int sweeps = 0; int stg; int nadd; int nrem; int pst;
ideal NEWE;
int tstart = timer;
for (stg = 0; stg <= 15; stg++) {
  // add the order-stg equations, with all earlier substitutions applied
  NEWE = 0; nadd = 0;
  for (t = 1; t <= nk; t++) { if (KORD[t] == stg) { nadd++; NEWE[nadd] = EQ0[t]; } }
  map PHX = RB, X;
  NEWE = PHX(NEWE);
  kill PHX;
  EQW = simplify(EQW + NEWE, 2);
  pst = npiv;
  thr = 2;
  while (1) {
    progress = 0; maxsz = 0; sweeps++;
    for (t = 1; t <= ncols(EQW); t++) {
      f = EQW[t];
      if (f != 0) {
        if (size(f) > maxsz) { maxsz = size(f); }
        if (size(f) <= thr) {
          if (deg(f) == 0) { incons++; }
          else {
            lin = jet(f, 1) - jet(f, 0);
            if (lin != 0) {
              hi = f - jet(f, 1);
              vl = variables(lin); vh = variables(hi);
              best = 0;
              for (u = 1; u <= ncols(vl); u++) {
                if (vl[u] != 0) {
                  ok = 1;
                  for (s = 1; s <= ncols(vh); s++) { if (vh[s] == vl[u]) { ok = 0; } }
                  if (ok) {
                    v = rvar(vl[u]);
                    if (best == 0) { best = v; }
                    else {
                      if ((v <= NL) && (best > NL)) { best = v; }
                      else { if (((v <= NL) == (best <= NL)) && (v > best)) { best = v; } }
                    }
                  }
                }
              }
              if (best > 0) {
                cf = diff(f, var(best));
                if (deg(cf) != 0) { "ERROR nonconstant pivot"; }
                ex = var(best) - f / cf;
                EQW = subst(EQW, var(best), ex);
                X = subst(X, var(best), ex);
                npiv++; ELI[best] = npiv;
                progress = 1;
              }
            }
          }
        }
      }
    }
    if (progress) { thr = 2; }
    else { if (thr >= maxsz) { break; } thr = thr + 1; }
  }
  EQW = simplify(EQW, 2 + 8);
  "stage " + string(stg) + ": equations added " + string(nadd) + ", pivots " + string(npiv - pst) + ", equations remaining " + string(size(EQW)) + ", max terms " + string(maxsz) + ", CPU s " + string(timer - tstart);
}
"elimination CPU seconds: " + string(timer - tstart) + "   sweeps: " + string(sweeps);
EQW = simplify(EQW, 2 + 8);
"pivots: " + string(npiv) + "   nonzero constant equations met: " + string(incons);
int nlel = 0; int nael = 0;
for (v = 1; v <= NU; v++) { if (ELI[v] > 0) { if (v <= NL) { nlel++; } else { nael++; } } }
"eliminated lambda: " + string(nlel) + "   eliminated A: " + string(nael);
"free lambda: " + string(NL - nlel) + "   free A: " + string(NU - NL - nael);
string fa = "";
for (v = NL + 1; v <= NU; v++) { if (ELI[v] == 0) { fa = fa + " " + varstr(v); } }
"free A-unknowns:" + fa;
string ea = "";
for (v = NL + 1; v <= NU; v++) { if (ELI[v] > 0) { ea = ea + " " + varstr(v) + "(step " + string(ELI[v]) + ")"; } }
"eliminated A-unknowns (with step number):" + ea;
"leftover equations: " + string(size(EQW));
map PHI = RB, X;
ideal CHK = PHI(EQ);
"original equations (all, with scalar multiples) not identically 0 after the final substitution: " + string(size(simplify(CHK, 2)));
ideal CHKn = simplify(CHK, 1+2+4);
ideal EQWn = simplify(EQW, 1+2+4);
"distinct normalized substituted equations: " + string(size(CHKn)) + " ; leftover (normalized, distinct): " + string(size(EQWn)) + " ; union: " + string(size(simplify(EQWn + CHKn, 1+2+4)));
int graphok = 1; ideal va;
for (v = NL + 1; v <= NU; v++) {
  if (ELI[v] > 0) {
    va = variables(X[v]);
    for (u = 1; u <= ncols(va); u++) { if (va[u] != 0) { if (rvar(va[u]) <= NL) { graphok = 0; } } }
  }
}
"graph check (no lambda in the images of eliminated A-unknowns): " + string(graphok);
for (v = NL + 1; v <= NU; v++) { if (ELI[v] > 0) { varstr(v) + " = " + string(X[v]); } }
"""

POST_C1_HEAD = r"""
// ---------------- c = 1: random rational point, exact Jacobian ranks, membership test, torus ----------------
LIB "matrix.lib";
write("ssi:w vn2b_c1_state.ssi", X);
system("random", 20261109);
ideal RV = maxideal(1);
for (v = 1; v <= NU; v++) { if (ELI[v] == 0) { RV[v] = random(1, 7) * (2*random(0, 1) - 1); } }
map RND = RB, RV;
ideal PT = RND(X);
ideal PTA = PT[NL+1..NU];
"random point: values of the A-unknowns (g10..g1, b12..b1): " + string(PTA);
map ATPT = RB, PT;
ideal EVQ = ATPT(EQ);
"equations (all 720) not vanishing at the point: " + string(size(simplify(EVQ, 2)));
int tj = timer;
matrix JJ = jacob(EQ0);
matrix JP = ATPT(JJ);
matrix JPU = submat(JP, 1..nrows(JP), 1..NU);
matrix JPL = submat(JP, 1..nrows(JP), 1..NL);
int rkJ = rank(JPU);
int rkL = rank(JPL);
"exact rank over Q of the Jacobian (630 distinct eqs x 393 unknowns): " + string(rkJ) + "  -> tangent dim " + string(NU - rkJ);
"exact rank over Q of the lambda-columns: " + string(rkL) + "  -> rank of the projection of the tangent space to the 22 A-coordinates: " + string(22 - rkJ + rkL);
"Jacobian CPU seconds: " + string(timer - tj);
ideal T2 = maxideal(1);
"""

POST_C1_TAIL = r"""
map PHT2 = RB, T2;
int t2bad = 0; poly tt; poly tt2;
for (v = NL + 1; v <= NU; v++) { if (ELI[v] > 0) { tt = X[v]; tt2 = T2[v]; if (PHT2(tt) != (tt2/var(v))*tt) { t2bad++; "  not tau2-homogeneous: " + varstr(v); } } }
"images of eliminated A-unknowns that are not tau2-homogeneous of the right weight: " + string(t2bad);
// membership test without lambda: B_ij in J = (B_00..B_05) in Q[w,y] (inline blocks generated by vgen2.py)
ring RW = 0, (y3, y2, y1, w), wp(2, 2, 2, 1);
ideal PTA = imap(RB, PTA);
int nb; int mi; int mj; ideal tgt; ideal J; ideal SJ; matrix BN[6][6];
@MEM:PTA:at the random graph point@
ideal P2 = PTA; P2[5] = P2[5] + 1;
@MEM:P2:perturbation g6 -> g6 + 1@
ideal P3 = PTA; P3[20] = P3[20] + 1;
@MEM:P3:perturbation b3 -> b3 + 1@
ideal P4 = PTA; P4[1] = P4[1] + 1;
@MEM:P4:perturbation g10 -> g10 + 1@
ideal P5 = PTA; P5[15] = P5[15] + 1;
@MEM:P5:perturbation b8 -> b8 + 1@
setring RB;
system("random", 777);
ideal RV2 = maxideal(1);
for (v = NL + 1; v <= NU; v++) { if (ELI[v] == 0) { RV2[v] = random(1, 9) * (2*random(0, 1) - 1); } }
map RND2 = RB, RV2;
ideal PTB = RND2(X);
ideal PTB2 = PTB[NL+1..NU];
"second graph point (A-values g10..g1, b12..b1): " + string(PTB2);
setring RW;
ideal PTB2 = imap(RB, PTB2);
@MEM:PTB2:at the second graph point@
// torus tau2 acts on the normal form (exact, parameter tq; the name t is taken by an int loop variable)
ring RT = (0, tq), (@SMALL@), dp;
matrix A = imap(RS, A);
ideal ph = maxideal(1); ph[25] = tq*y1; ph[24] = tq^2*y2;
map PHI2 = RT, ph;
matrix Aph = PHI2(A);
matrix PP[6][6]; PP[1,1] = 1/tq^3; PP[2,2] = 1/tq; PP[3,3] = 1; PP[4,4] = 1; PP[5,5] = 1/tq; PP[6,6] = tq;
ideal sc = maxideal(1);
"""


def post_c1():
    wt = tau2_weights()
    out = [POST_C1_HEAD]
    for n, name in enumerate(Avars):
        w_ = wt[name]
        fac = ('%d' % (2 ** w_)) if w_ >= 0 else ('1/%d' % (2 ** (-w_)))
        out.append('T2[NL+%d] = (%s)*%s;' % (n + 1, fac, name))
    out.append('"tau2 weights (g10..g1, b12..b1): %s";' % ' '.join('%s:%d' % (nm, wt[nm]) for nm in Avars))
    out.append(POST_C1_TAIL.replace('@SMALL@', small_vars))
    for n, name in enumerate(Avars):
        w_ = wt[name]
        out.append('sc[%d] = tq^%d*%s;' % (n + 1, w_, name) if w_ >= 0 else 'sc[%d] = %s/tq^%d;' % (n + 1, name, -w_))
    out.append('map SC = RT, sc;')
    out.append('matrix Asc = SC(A);')
    out.append('"CHECK tau2: P*A(phi)*P = A with A-unknowns scaled by tq^weight: " + string(iszero(PP*Aph*PP - Asc));')
    txt = '\n'.join(out)
    import re
    def mem(mo):
        nm, label = mo.group(1), mo.group(2)
        return ('tgt = %s, y3, y2, y1, w;\nmap EVA = RS, tgt;\nBN = EVA(B);\n'
                'J = BN[1, 1], BN[1, 2], BN[1, 3], BN[1, 4], BN[1, 5], BN[1, 6];\nSJ = std(J);\nnb = 0;\n'
                'for (mi = 2; mi <= 6; mi++) { for (mj = mi; mj <= 6; mj++) { if (reduce(BN[mi, mj], SJ) != 0) { nb++; } } }\n'
                '"membership (%s): number of B_ij (1<=i<=j<=5) NOT in (B_00..B_05): " + string(nb);\nkill EVA;') % (nm, label)
    return re.sub(r'@MEM:([A-Za-z0-9]+):([^@]+)@', mem, txt)


def post_c0():
    out = ['// ---------------- c = 0: deciding coefficients at p = (0:0:1:0|0) ----------------']
    for (idx, i, j, k, m) in decide:
        out.append('"x(%d) = lambda_{%d%d,0}[%s] : eliminated " + string(ELI[%d] > 0) + ", image = " + string(X[%d]);'
                   % (idx, i, j, mstr(m), idx, idx))
    out.append('write("ssi:w vn2b_c0_state.ssi", X);')
    out.append('write("ssi:w vn2b_c0_leftover.ssi", EQW);')
    out.append('int dmax = 0; int dmin = 100; for (t = 1; t <= size(EQW); t++) { if (deg(EQW[t]) > dmax) { dmax = deg(EQW[t]); } if (ord(EQW[t]) < dmin) { dmin = ord(EQW[t]); } }')
    out.append('"leftover: max degree " + string(dmax) + ", min degree of a term " + string(dmin);')
    return '\n'.join(out)


for cval in (1, 0):
    body = ['// vn2b_c%d.sing -- independent verification of case n2, c = %d (generated by vgen2.py)' % (cval, cval)]
    body.append('ring RS = 0, (%s), dp;' % small_vars)
    body.append(PROCS)
    body.append(MATRIX.replace('@C@', str(cval)))
    body.append(CHECKS)
    body.append('ring RB = 0, (%s), dp;' % big_vars)
    body.append(BIGRING_HEAD.replace('@NL@', str(NL)))
    for (i, j) in pairs_rev:
        wt_ = 8 + cdeg[i] + cdeg[j]
        sg = 0 if -delta[i] * delta[j] == 1 else 1
        expr = 'B[%d, %d]' % (i + 1, j + 1)
        for k in range(6):
            le = lam_expr(i, j, k)
            if le != '0':
                expr += ' - %s*B[1, %d]' % (le, k + 1)
        body.append('E = %s;' % expr)
        body.append(EQ_BLOCK.replace('@WT@', str(wt_)).replace('@SG@', str(sg)))
    body.append(COUNTS)
    body.append(ELIM)
    body.append(post_c1() if cval == 1 else post_c0())
    body.append('"memory MB: " + string(memory(2) div 1048576);')
    body.append('quit;')
    with open(os.path.join(HERE, 'vn2b_c%d.sing' % cval), 'w') as fh:
        fh.write('\n'.join(body) + '\n')
print('written vn2b_c1.sing, vn2b_c0.sing; NL =', NL, '; deciding:', [(d[0], d[1], d[2], mstr(d[4])) for d in decide])
