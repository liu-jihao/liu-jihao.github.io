# Generator for n1c1 script D: the residual diagonal torus (plan/case_common.md section 4), c = 1 fixed.
# Logarithmic weights: w -> t^W w, y_i -> t^Yi y_i, A -> P^t A(s.) P with P = diag(t^P0, ..., t^P5); D -> t^dD D.
# A coefficient a of a monomial mu in A_jk gets weight wt(mu) + Pj + Pk; a lambda coefficient of mu in lambda_{ij,k}
# gets weight wt(mu) + P0 + Pk - Pi - Pj (from adj(P^t A P) = det(P)^2 P^-1 adj(A) P^-1).
import os
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, 'gen_n1c1.py')).read().split('# ---------------- common include')[0])
P = ['P0', 'P1', 'P2', 'P3', 'P4', 'P5']
def wtmon(a, b, e, f):
    t = []
    if a: t.append('%d*W' % a)
    if b: t.append('%d*Y1' % b)
    if e: t.append('%d*Y2' % e)
    if f: t.append('%d*Y3' % f)
    return t
def parse_y(m):  # monomial string in y0..y3 -> exponents (a(w), b, e, f) with y0 = w^2
    ex = {'y0': 0, 'y1': 0, 'y2': 0, 'y3': 0}
    for tkn in m.split('*'):
        if '^' in tkn:
            v, p = tkn.split('^'); ex[v] += int(p)
        else:
            ex[tkn] += 1
    return (2 * ex['y0'], ex['y1'], ex['y2'], ex['y3'])
forms = []   # linear form (string) of the weight of each unknown, in the ring order l(1..371), g(1..10), b(1..12), D
for (i, j, k, ex) in lam:
    t = wtmon(*ex) + ['P0', P[k], '-' + P[i], '-' + P[j]]
    forms.append(' + '.join(t).replace('+ -', '- '))
for m in Gmons:              # entry A00 = y0*G
    a, b_, e, f = parse_y(m)
    t = wtmon(a + 2, b_, e, f) + ['2*P0']
    forms.append(' + '.join(t))
qcol = {'q1': 1, 'q2': 2, 'q3': 3, 'q4': 4}
for (q, m) in bmons:         # entry A0k = w*q_k
    a, b_, e, f = parse_y(m)
    t = wtmon(a + 1, b_, e, f) + ['P0', P[qcol[q]]]
    forms.append(' + '.join(t))
forms.append('dD')
assert len(forms) == 394
# fixed entries of the normal form N1(D), c = 1 (0-indexed rows/cols): conditions on the weights
conds = [
    ('A05: y2^2', 'P0 + P5 + 2*Y2'), ('A05: y1^2', 'P0 + P5 + 2*Y1'), ('A05: D*y3^2', 'P0 + P5 + 2*Y3 - dD'),
    ('A15: w', 'P1 + P5 + W'), ('A11: D*y3', '2*P1 + Y3 - dD'), ('A12: y2', 'P1 + P2 + Y2'), ('A13: y1', 'P1 + P3 + Y1'),
    ('A22: y3', '2*P2 + Y3'), ('A23: c*y0', 'P2 + P3 + 2*W'), ('A24: y1', 'P2 + P4 + Y1'),
    ('A33: -y3', '2*P3 + Y3'), ('A34: y2', 'P3 + P4 + Y2'), ('A44: -D*y3', '2*P4 + Y3 - dD'),
]
L = []
L.append('// n1c1, script D: residual torus of the normal form N1(D), c = 1')
L.append('system("--ticks-per-sec", 1000); int T0 = rtimer;')
L.append('ring RW = 0, (W, Y1, Y2, Y3, P0, P1, P2, P3, P4, P5, dD), dp;')
L.append('ideal CND = ' + ', '.join(f for (_, f) in conds) + ';')
L.append('"conditions (fixed entries):", %d;' % len(conds))
L.append('ideal SCND = std(CND); "dimension of the solution space (torus dimension, D allowed to move):", dim(SCND);')
L.append('"reduced linear conditions:"; SCND;')
L.append('// candidate basis: grading direction gr and the one-parameter group u (t^2 = u)')
L.append('ideal GR = 1, 2, 2, 2, -4, -1, -1, -1, -1, 0, 0;   // W, Y1, Y2, Y3, P0..P5 = c_j - 4, dD')
L.append('ideal UU = -1, 0, 0, -2, -2, -1, 1, 1, -1, 2, -4;  // w -> w/t, y3 -> y3/t^2, P = diag(t^-2, t^-1, t, t, t^-1, t^2), D -> D/t^4')
L.append('map mgr = RW, GR; map muu = RW, UU;')
L.append('"grading direction satisfies all conditions:", size(mgr(CND)) == 0, " ; u satisfies all conditions:", size(muu(CND)) == 0;')
L.append('matrix BM[2][11]; int k; for (k=1;k<=11;k++){ BM[1,k] = GR[k]; BM[2,k] = UU[k]; } "rank of (gr, u):", rank(BM);')
L.append('ideal WF = ' + ', '.join(forms) + ';')
L.append('ideal WGR = mgr(WF); ideal WUU = muu(WF);')
L.append('int ngr = 0; for (k=1;k<=394;k++){ if (WGR[k] != 0) { ngr++; } } "coefficients with nonzero weight under the grading direction:", ngr;')
L.append('intvec WTU; for (k=1;k<=394;k++){ WTU[k] = int(leadcoef(WUU[k])); if (WUU[k] == 0) { WTU[k] = 0; } }')
L.append('ring R1 = 0, (l(1..371), g(1..10), b(1..12), D), dp;')
L.append('short = 0;')
L.append('int k, j, n;')
L.append('intvec FREEV; execute(read("free_n1c1.txt"));')
L.append('string s = ""; for (k=372;k<=394;k++){ s = s + string(var(k)) + ":" + string(WTU[k]) + " "; } "u-weights (t-normalization, t^2 = u) of D, g, b:", s;')
L.append('s = ""; for (k=1;k<=size(FREEV);k++){ if (FREEV[k] > 371) { s = s + string(var(FREEV[k])) + ":" + string(WTU[FREEV[k]]) + " "; } } "u-weights of the free A-coefficients:", s;')
L.append('int nzl = 0; for (k=1;k<=size(FREEV);k++){ if (FREEV[k] <= 371 && WTU[FREEV[k]] != 0) { nzl++; } } "free lambda with nonzero u-weight:", nzl, "of 94";')
L.append('s = ""; for (k=1;k<=371;k++){ s = s + string(WTU[k]) + " "; } "u-weights of l(1..371):", s;')
L.append('// effective dimension e: rank of the weight matrix (rows gr, u) on the free A-coefficients, and on all free coordinates')
L.append('int nfa = 0; for (k=1;k<=size(FREEV);k++){ if (FREEV[k] > 371) { nfa++; } }')
L.append('matrix WMA[2][nfa]; n = 0; for (k=1;k<=size(FREEV);k++){ if (FREEV[k] > 371) { n++; WMA[1,n] = 0; WMA[2,n] = WTU[FREEV[k]]; } }')
L.append('matrix WMF[2][size(FREEV)]; for (k=1;k<=size(FREEV);k++){ WMF[1,k] = 0; WMF[2,k] = WTU[FREEV[k]]; }')
L.append('"e = rank of the weight matrix on the free A-coefficients:", rank(WMA), " ; on all 103 free coordinates:", rank(WMF);')
L.append('// verification 1: the u-element maps the normal form to itself, with the predicted weights of D, g, b')
L.append('ring RT = (0,t), (l(1..371), g(1..10), b(1..12), D, w, y1, y2, y3), dp;')
L.append('short = 0; int k, j;')
L.append('< "common_n1c1_A.sing";')
L.append('matrix PP[6][6]; PP[1,1] = 1/t^2; PP[2,2] = 1/t; PP[3,3] = t; PP[4,4] = t; PP[5,5] = 1/t; PP[6,6] = t^2;')
L.append('ideal SG = maxideal(1); SG[395] = w/t; SG[398] = y3/t^2; map sgm = RT, SG;')
L.append('matrix AS = sgm(A); matrix AT = transpose(PP)*AS*PP;')
L.append('ideal SC = maxideal(1); for (k=1;k<=394;k++){ if (WTU[k] >= 0) { SC[k] = t^WTU[k]*var(k); } else { SC[k] = var(k)/t^(-WTU[k]); } } map scl = RT, SC;')
L.append('matrix AW = scl(A);')
L.append('int okm = 1; for (k=1;k<=6;k++){ for (j=1;j<=6;j++){ if (AT[k,j] != AW[k,j]) { okm = 0; } } }')
L.append('"P^t A(w/t, y1, y2, y3/t^2) P equals A with D, g, b scaled by t^weight:", okm;')
L.append('"for comparison, the fixed entries after the change (row v1, row v4, A05):", AT[2,2], AT[5,5], AT[1,6];')
L.append('// verification 2: every coefficient equation is homogeneous for these weights (so the scaling maps solutions to solutions)')
L.append('ideal EQ; intvec EO; execute(read("eqs_n1c1.txt"));')
L.append('int nh = 0; poly fs, f0;')
L.append('for (k=1;k<=ncols(EQ);k++){ f0 = EQ[k]; fs = scl(f0); if (fs*leadcoef(f0) - leadcoef(fs)*f0 != 0) { nh++; } }')
L.append('"coefficient equations (of", ncols(EQ), ") that are not weighted-homogeneous for the u-weights:", nh;')
L.append('"total ms:", rtimer-T0, " memory KB:", memory(2) div 1024;')
L.append('quit;')
open(os.path.join(HERE, 'n1c1_D_torus.sing'), 'w').write('\n'.join(L) + '\n')
# the matrix A in ring RT (same text as common, without ring declaration and lambda polys)
src = open(os.path.join(HERE, 'common_n1c1.sing')).read().split('\n')
keep = [ln for ln in src if not ln.startswith('ring ') and not ln.startswith('poly lam_') and not ln.startswith('short') and not ln.startswith('proc cofac') and not ln.startswith('< ')]
open(os.path.join(HERE, 'common_n1c1_A.sing'), 'w').write('\n'.join(keep) + '\n')
print('written n1c1_D_torus.sing, common_n1c1_A.sing')
