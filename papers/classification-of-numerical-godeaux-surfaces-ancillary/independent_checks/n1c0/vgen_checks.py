#!/usr/bin/env python3
"""Writes v8_check_iii.sing (substitution check of the branch (ii) and (iii) parametrizations into all original equations,
exact lambda-ranks blockwise at rational points of (iii), det A, exact syzygy count), v9_jac_modp.sing (Jacobian ranks
mod 32003 at random points of (iii) and (ii)), v10_torus.sing (semi-invariance of all equations under the two torus
directions s, u derived by hand in VERIFY.md).  Torus weights, the trivial-syzygy count and the random points are
computed here."""
import json, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(HERE, 'v_data.json')))
lam = data['lam']; pairs = [tuple(p) for p in data['pairs']]
c = [0, 3, 3, 3, 3, 4]; dl = [1, -1, -1, 1, 1, -1]

def monos(weight, sign):
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

# trivial syzygy parameters: c^{(ij)}_m homogeneous of weight c_i + c_j - (8 - c_m), sign -d_i d_j d_m
ntriv = 0; triv_by_pair = {}
for (i, j) in pairs:
    t = 0
    for m in range(1, 6):
        t += len(monos(c[i] + c[j] - (8 - c[m]), -dl[i] * dl[j] * dl[m]))
    triv_by_pair[(i, j)] = t; ntriv += t

# torus: hand derivation (VERIFY.md section 6). Log-weights (ww, w1, w2, w3; g0..g5).
S_DIR = (1, 0, 0, 0, [1, 0, 0, 0, 0, -1])
U_DIR = (-1, 0, 0, -2, [-2, -1, 1, 1, -1, 2])
GR_DIR = (1, 2, 2, 2, [-4, -1, -1, -1, -1, 0])

def weights(dirn):
    ww, w1, w2, w3, g = dirn
    def mono_w(a, b, e, f):
        return a * ww + b * w1 + e * w2 + f * w3
    W = {}
    # fixed-entry constraints (must be 0)
    cons = {
        'A12=y2': g[1] + g[2] + w2, 'A13=y1': g[1] + g[3] + w1, 'A22=y3': 2 * g[2] + w3, 'A24=y1': g[2] + g[4] + w1,
        'A33=-y3': 2 * g[3] + w3, 'A34=y2': g[3] + g[4] + w2, 'A15=w': g[1] + g[5] + ww,
        'A05:y2^2': g[0] + g[5] + 2 * w2, 'A05:y1^2': g[0] + g[5] + 2 * w1,
        'D consistent (A11 vs A44)': (2 * g[1] + w3) - (2 * g[4] + w3),
        'D consistent (A11 vs A05)': (2 * g[1] + w3) - (2 * w3 - 2 * w2)}
    W['D'] = 2 * g[1] + w3
    # b's: q_k sits in entry (0,k) as w*q_k
    bdef = [(1, (0, 1, 1, 0)), (1, (0, 1, 0, 1)), (2, (2, 0, 1, 0)), (2, (2, 0, 0, 1)), (2, (0, 1, 1, 0)), (2, (0, 1, 0, 1)),
            (3, (4, 0, 0, 0)), (3, (2, 1, 0, 0)), (3, (0, 2, 0, 0)), (4, (4, 0, 0, 0)), (4, (0, 0, 1, 1)), (4, (0, 0, 0, 2))]
    for idx, (k, (a, b, e, f)) in enumerate(bdef):
        W['b(%d)' % (idx + 1)] = g[0] + g[k] + mono_w(a + 1, b, e, f)   # extra w from w*q_k; y0 = w^2 counted in a
    gdef = [(0, 0, 0, 3), (0, 2, 0, 1), (0, 0, 1, 2), (0, 2, 1, 0), (0, 0, 2, 1), (0, 0, 3, 0), (2, 1, 0, 1), (2, 1, 1, 0), (4, 0, 0, 1), (4, 0, 1, 0)]
    for idx, (a, b, e, f) in enumerate(gdef):
        W['g(%d)' % (idx + 1)] = 2 * g[0] + mono_w(a + 2, b, e, f)      # A00 = y0*G
    for (n, i, j, k, m) in lam:
        W['L(%d)' % n] = g[0] + g[k] - g[i] - g[j] + mono_w(*m)
    return W, cons

Ws, cs = weights(S_DIR); Wu, cu = weights(U_DIR); Wg, cg = weights(GR_DIR)
names = ['D'] + ['b(%d)' % k for k in range(1, 13)] + ['g(%d)' % k for k in range(1, 11)] + ['L(%d)' % n for n in range(1, 372)]
assert len(names) == 394
gr_nonzero = [nm for nm in names if Wg[nm] != 0]

# random points (fixed seed), values of the free unknowns of (iii) and (ii) are assigned in Singular from these lists
random.seed(20261109)
def rnd():
    v = 0
    while v == 0:
        v = random.randint(-9, 9)
    return v
P1 = [rnd() for _ in range(394)]
P2 = [rnd() for _ in range(394)]
P3 = [rnd() for _ in range(394)]

def ivec(lst):
    return ','.join(str(x) for x in lst)

with open(os.path.join(HERE, 'v_checks_meta.txt'), 'w') as fh:
    fh.write('trivial syzygy parameters: %d; by pair: %s\n' % (ntriv, {('%d%d' % k): v for k, v in triv_by_pair.items()}))
    fh.write('torus constraints (must all be 0):\n  s: %s\n  u: %s\n  grading: %s\n' % (cs, cu, cg))
    fh.write('grading direction: unknowns with nonzero weight: %s\n' % gr_nonzero)
    fh.write('weights (s, u) of the A-unknowns: %s\n' % {nm: (Ws[nm], Wu[nm]) for nm in names[:23]})

# ---------- v8 ----------
V8 = r'''// v8_check_iii.sing: (1) the parametrizations of branches (ii), (iii) (my own elimination, v6/v7) satisfy all 870
// original coefficient equations identically modulo t*b8 - 1; (2) at three rational points of (iii): all E_ij vanish,
// det A != 0, exact rank of the lambda-linear system (blockwise over the 15 pairs), and the count of trivial syzygies.
LIB "matrix.lib";
system("--ticks-per-sec", 1000);
int T0 = rtimer;
< "vbuild.sing";
poly detA = 0; for (k = 1; k <= 6; k++) { detA = detA + A[1,k]*B0[1,k]; }
execute(read("v_param_iii.txt"));
execute(read("v_param_ii.txt"));
ideal Gt = std(t*b(8) - 1);
int ii; int kk;
map phi3 = R, w, y1, y2, y3, TRiii, t;
ideal Z = reduce(phi3(EQ), Gt);
"branch (iii): original equations not satisfied identically by the parametrization: " + string(size(Z));
map phi2 = R, w, y1, y2, y3, TRii, t;
Z = reduce(phi2(EQ), Gt);
"branch (ii): original equations not satisfied identically by the parametrization: " + string(size(Z));
// free unknowns of (iii)
string sf = ""; int nfA = 0; int nfL = 0;
for (ii = 5; ii <= 398; ii++) { if (ELiii[ii] == 0) { if (ii <= 27) { nfA++; sf = sf + string(var(ii)) + " "; } else { nfL++; } } }
"branch (iii) free unknowns: A " + string(nfA) + " (" + sf + "), lambda " + string(nfL);
proc evalpt(intvec pv, intvec el, ideal TRX)
{
  // assign pv[i] to every free unknown var(4+i) (el[4+i] == 0), others 0; t = 1/b8; return the images of all unknowns
  int i2; ideal vals;
  for (i2 = 1; i2 <= 394; i2++) { if (el[4+i2] == 0) { vals[i2] = pv[i2]; } else { vals[i2] = 0; } }
  number b8v = leadcoef(vals[9]);
  map ev = basering, w, y1, y2, y3, vals, 1/b8v;
  ideal res = ev(TRX);
  return(res);
}
list PTS = list(intvec(''' + ivec(P1) + r'''), intvec(''' + ivec(P2) + r'''), intvec(''' + ivec(P3) + r'''));
int pt; ideal VAL; int j2; matrix Mb; int rk; int rka; int rtot; int ratot; poly dA;
for (pt = 1; pt <= 3; pt++)
{
  VAL = evalpt(PTS[pt], ELiii, TRiii);
  "point " + string(pt) + ": D = " + string(VAL[1]) + ", b3 = " + string(VAL[4]) + ", b8 = " + string(VAL[9]);
  // all E_ij at the point (values of all unknowns)
  map evall = R, w, y1, y2, y3, VAL, 0;
  ideal EEp = evall(EE);
  "  E_ij at the point, nonzero: " + string(size(EEp));
  // det A at the point
  ideal VALA = VAL; for (j2 = 24; j2 <= 394; j2++) { VALA[j2] = var(4+j2); }
  map evA = R, w, y1, y2, y3, VALA, t;
  dA = evA(detA);
  "  det A at the point: nonzero " + string(dA != 0) + ", number of terms " + string(size(dA));
  // lambda-linear system at the A-values, blockwise by pair
  ideal EQa = evA(EQ);
  rtot = 0; ratot = 0;
  for (kk = 1; kk <= 15; kk++)
  {
    ideal blk; int nb = 0;
    for (j2 = 1; j2 <= ncols(EQa); j2++) { if (PR[j2] == kk) { if (EQa[j2] != 0) { nb++; blk[nb] = EQa[j2]; } } }
    if (nb > 0)
    {
      matrix J1 = jacob(blk);
      ideal z0 = blk; for (j2 = 28; j2 <= 398; j2++) { z0 = subst(z0, var(j2), 0); }
      matrix Jaug = concat(J1, transpose(matrix(z0)));
      rtot = rtot + rank(J1); ratot = ratot + rank(Jaug);
      kill J1; kill z0; kill Jaug;
    }
    kill blk; kill nb;
  }
  "  exact rank of the lambda-system (sum over the 15 pairs): " + string(rtot) + "; augmented: " + string(ratot) + "; kernel dimension " + string(371 - rtot);
  kill evall; kill EEp; kill VALA; kill evA; kill EQa;
}
"trivial syzygy parameters (count from vgen_checks.py): ''' + str(ntriv) + r'''";
"total time (ms): " + string(rtimer - T0) + ", memory(2) = " + string(memory(2));
quit;
'''
open(os.path.join(HERE, 'v8_check_iii.sing'), 'w').write(V8)

# ---------- v9: Jacobian ranks mod 32003 ----------
V9 = r'''// v9_jac_modp.sing: Jacobian of all 752 distinct equations w.r.t. the 394 unknowns, at two points of (iii) and one
// point of (ii) (values from the exact parametrizations, reduced mod 32003); rank mod p and the rank of the projection
// of the kernel to the 23 A-coordinates. Evidence only: rank mod p <= rank over Q.
system("--ticks-per-sec", 1000);
int T0 = rtimer;
< "vring.sing";
execute(read("v_eqs.txt"));
execute(read("v_param_iii.txt"));
execute(read("v_param_ii.txt"));
ideal EQd = simplify(EQ, 1); EQd = simplify(EQd, 4); EQd = simplify(EQd, 2);
"distinct equations: " + string(ncols(EQd));
matrix JJ = jacob(EQd);   // columns = ring variables (399)
proc evalpt(intvec pv, intvec el, ideal TRX)
{
  int i2; ideal vals;
  for (i2 = 1; i2 <= 394; i2++) { if (el[4+i2] == 0) { vals[i2] = pv[i2]; } else { vals[i2] = 0; } }
  number b8v = leadcoef(vals[9]);
  map ev = basering, w, y1, y2, y3, vals, 1/b8v;
  ideal res = ev(TRX);
  return(res);
}
list PTS = list(intvec(''' + ivec(P1) + r'''), intvec(''' + ivec(P2) + r'''), intvec(''' + ivec(P3) + r'''));
list VALS = list(evalpt(PTS[1], ELiii, TRiii), evalpt(PTS[2], ELiii, TRiii), evalpt(PTS[3], ELii, TRii));
list LAB = list("(iii) point 1", "(iii) point 2", "(ii) point 3");
// keep only the columns of the 394 unknowns
matrix JU[ncols(EQd)][394];
int a; int bcol;
for (a = 1; a <= ncols(EQd); a++) { for (bcol = 1; bcol <= 394; bcol++) { JU[a,bcol] = JJ[a, bcol+4]; } }
ring Rp = 32003, (w, y1, y2, y3, D, b(1..12), g(1..10), L(1..371), t), dp;
matrix JU = imap(R, JU);
list VALS = imap(R, VALS);
int pt; int i3; int rk; int rkA;
for (pt = 1; pt <= 3; pt++)
{
  ideal vv = VALS[pt];
  map ev = Rp, 0, 0, 0, 0, vv, 0;
  matrix Jp = ev(JU);
  rk = rank(Jp);
  // rank of the projection of the kernel to the 23 A-coordinates = 23 - rk + rank(lambda-columns)
  matrix JL[nrows(Jp)][371];
  for (a = 1; a <= nrows(Jp); a++) { for (i3 = 1; i3 <= 371; i3++) { JL[a, i3] = Jp[a, i3 + 23]; } }
  rkA = rank(JL);
  LAB[pt] + ": rank mod 32003 = " + string(rk) + ", nullity = " + string(394 - rk) + ", rank of the lambda-columns = " + string(rkA) + ", rank of the kernel projection to the A-coordinates = " + string(23 - rk + rkA);
  kill vv; kill ev; kill Jp; kill JL;
}
"total time (ms): " + string(rtimer - T0) + ", memory(2) = " + string(memory(2));
quit;
'''
open(os.path.join(HERE, 'v9_jac_modp.sing'), 'w').write(V9)

# ---------- v10: torus semi-invariance ----------
M = 40
def tu_map(W):
    imgs = []
    for nm in names:
        wt = W[nm] + M
        assert wt >= 0
        imgs.append('%s*T^%d*U' % (nm, wt))
    return ', '.join(imgs)

V10 = r'''// v10_torus.sing: semi-invariance of all 870 coefficient equations under the torus directions s and u
// (weights from vgen_checks.py; hand derivation in VERIFY.md). Trick: x -> x*T^(wt(x)+40)*U; a term of degree d and
// weight W becomes T^(W+40d)*U^d, so W = deg_T - 40*deg_U must be constant on each equation.
system("--ticks-per-sec", 1000);
int T0 = rtimer;
ring R = 0, (w, y1, y2, y3, D, b(1..12), g(1..10), L(1..371), t, T, U), dp;
execute(read("v_eqs.txt"));
proc checkhom(ideal I)
{
  int bad = 0; int k; int j; int w0; int wj; intvec e;
  for (k = 1; k <= ncols(I); k++)
  {
    if (I[k] != 0)
    {
      matrix C = coef(I[k], T*U);
      e = leadexp(C[1,1]); w0 = e[nvars(basering)-1] - 40*e[nvars(basering)];
      for (j = 2; j <= ncols(C); j++) { e = leadexp(C[1,j]); wj = e[nvars(basering)-1] - 40*e[nvars(basering)]; if (wj != w0) { bad++; break; } }
      kill C;
    }
  }
  return(bad);
}
map ps = R, w, y1, y2, y3, ''' + tu_map(Ws) + r''', t, T, U;
"equations not homogeneous for the s-weights: " + string(checkhom(ps(EQ)));
map pu = R, w, y1, y2, y3, ''' + tu_map(Wu) + r''', t, T, U;
"equations not homogeneous for the u-weights: " + string(checkhom(pu(EQ)));
"total time (ms): " + string(rtimer - T0) + ", memory(2) = " + string(memory(2));
quit;
'''
open(os.path.join(HERE, 'v10_torus.sing'), 'w').write(V10)
print('trivial syzygies', ntriv)
print('constraints s', cs)
print('constraints u', cu)
print('constraints grading', cg)
print('grading nonzero weights', gr_nonzero[:10], len(gr_nonzero))
for nm in names[:23]:
    print(nm, (Ws[nm], Wu[nm]))
