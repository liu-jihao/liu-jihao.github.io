# Generator for n1c1 script C: back-substitution, graph check, verification of the parametrization,
# and the exact syzygy-rank check at a random rational point (plan/case_common.md 3(e), 3(f)).
import os, random
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, 'gen_n1c1.py')).read().split('# ---------------- common include')[0])
random.seed(20261109)
pairs = [(i, j) for i in range(1, 6) for j in range(i, 6)]
# syzygy basis and lambda basis per pair
syz = {}
lbas = {}
tot = 0
for (i, j) in pairs:
    S = []
    for m in range(1, 6):
        n = c[i] + c[j] - r[m]
        if n < 0:
            continue
        for ex in mons(n, -dl[i] * dl[j] * dl[m]):
            S.append((m, mstr(ex)))
    syz[(i, j)] = S
    tot += len(S)
    K = []
    for k in range(6):
        for ex in mons(c[i] + c[j] - c[k], dl[i] * dl[j] * dl[k]):
            K.append('%s*e(%d)' % (mstr(ex), k + 1))
    lbas[(i, j)] = K
assert tot == 94, tot
freeA = ['g(5)', 'b(2)', 'b(5)', 'b(6)', 'b(8)', 'b(9)', 'b(11)', 'b(12)', 'D']
vals = {}
for v in freeA:
    x = 0
    while x == 0:
        x = random.randint(-9, 9)
    vals[v] = x
L = []
L.append('// n1c1, script C: back-substitution, graph check, verification of the parametrized family, syzygy rank')
L.append('system("--ticks-per-sec", 1000); int T0 = rtimer;')
L.append('< "common_n1c1.sing";')
L.append('list SV; list SE; intvec SSTG; intvec SEO; list SPIV; int NS;')
L.append('execute(read("subs_n1c1.txt"));')
L.append('intvec FREEV; execute(read("free_n1c1.txt"));')
L.append('"substitutions read:", NS, " free unknowns read:", size(FREEV);')
L.append('int i, j, k, u, ok, n;')
L.append('intvec isfree = 0:394; for (k=1;k<=size(FREEV);k++){ isfree[FREEV[k]] = 1; }')
L.append('// back-substitution in reverse elimination order: IMG[u] = expression of u in the free unknowns')
L.append('ideal IMG = maxideal(1);')
L.append('for (k = NS; k >= 1; k--) { map phi = R0, IMG; poly e = SE[k]; IMG[SV[k]] = phi(e); kill phi, e; }')
L.append('"back-substitution done, ms:", rtimer-T0;')
L.append('// sanity: every image involves only free unknowns (and no w, y)')
L.append('int bad = 0; int maxsz = 0; int totsz = 0;')
L.append('for (k = 1; k <= NS; k++) { u = SV[k]; poly e = IMG[u]; if (deg(e) > 0) { ideal vv = variables(e); for (j=1;j<=ncols(vv);j++){ n = rvar(vv[j]); if (n > 394) { bad++; } else { if (isfree[n] == 0) { bad++; } } } kill vv; } if (size(e) > maxsz) { maxsz = size(e); } totsz = totsz + size(e); kill e; }')
L.append('"images involving a non-free unknown or a coordinate:", bad, " ; image sizes: total", totsz, " max", maxsz;')
L.append('// graph check: eliminated g, b as polynomials in the free unknowns; no free lambda may occur')
L.append('int nlam = 0;')
L.append('write(":w graph_expr_n1c1.txt", "// eliminated g, b of case n1c1 as polynomials in the free A-coefficients D, g(5), b(2), b(5), b(6), b(8), b(9), b(11), b(12) (written by n1c1_C_graph.sing)");')
L.append('for (u = 372; u <= 393; u++) { if (isfree[u] == 0) { poly e = IMG[u]; int hasl = 0; if (deg(e) > 0) { ideal vv = variables(e); for (j=1;j<=ncols(vv);j++){ if (rvar(vv[j]) <= 371) { hasl = 1; } } kill vv; } nlam = nlam + hasl; string(var(u)), "=", e, "   [terms:", size(e), ", contains a lambda:", hasl, "]"; write(":a graph_expr_n1c1.txt", string(var(u)) + " = " + string(e) + ";"); kill e, hasl; } }')
L.append('"graph check: eliminated g, b whose expression contains a lambda-coefficient:", nlam;')
L.append('// verification: substitute the parametrization into A and the lambda_{ij,k} and recompute E_ij exactly')
L.append('map PHI = R0, IMG;')
L.append('matrix AP = PHI(A);')
L.append('matrix BP[1][6]; for (k=1;k<=6;k++){ BP[1,k] = cofac(AP,k,1); }')
L.append('poly pol; int nz = 0;')
for (i, j) in pairs:
    terms = ' + '.join('PHI(lam_%d_%d_%d)*BP[1,%d]' % (i, j, k, k + 1) for k in range(6))
    L.append('pol = cofac(AP,%d,%d) - (%s); if (pol != 0) { nz++; "E_%d%d nonzero, terms:", size(pol); }' % (j + 1, i + 1, terms, i, j))
L.append('"verification of the rank conditions E_ij = 0 (15 pairs) for the parametrized family: nonzero E_ij:", nz, " ms:", rtimer-T0;')
L.append('// exact syzygy-rank check at a random rational point of V (free A-coefficients below; lambda irrelevant)')
sub = ', '.join('%s, %d' % (v, vals[v]) for v in freeA)
L.append('"random point (free A-coefficients):", "%s";' % ', '.join('%s=%d' % (v, vals[v]) for v in freeA))
L.append('matrix AQ[6][6]; for (i=1;i<=6;i++){ for (j=1;j<=6;j++){ AQ[i,j] = subst(AP[i,j], %s); } }' % sub)
L.append('"A at the point:"; print(AQ);')
L.append('ring RS = 0, (w, y1, y2, y3, e(1..6)), dp;')
L.append('matrix AQ = imap(R0, AQ); int i, j, k, rk, rtot, okc; rtot = 0; okc = 1;')
L.append('poly dQ = det(AQ); "det A at the point is nonzero:", dQ != 0, " (terms:", size(dQ), ")";')
L.append('ideal P; ideal K; matrix M;')
for p, (i, j) in enumerate(pairs):
    S = syz[(i, j)]
    if not S:
        L.append('"pair (%d,%d): no syzygy directions";' % (i, j))
        continue
    gens = []
    for (m, mo) in S:
        gens.append(' + '.join('(%s)*AQ[%d,%d]*e(%d)' % (mo, m + 1, k + 1, k + 1) for k in range(6)))
    L.append('P = ' + ', '.join(gens) + ';')
    L.append('K = ' + ', '.join(lbas[(i, j)]) + ';')
    L.append('M = coeffs(P, K); if (matrix(K)*M != matrix(P)) { okc = 0; } rk = rank(M); rtot = rtot + rk; "pair (%d,%d): syzygy directions %d, rank", rk;' % (i, j, len(S)))
L.append('"all syzygy vectors lie in the lambda-space (coefficient extraction exact):", okc;')
L.append('"total rank of the trivial-syzygy map at the point:", rtot, " (number of directions: 94)";')
L.append('"total ms:", rtimer-T0, " memory KB:", memory(2) div 1024;')
L.append('quit;')
open(os.path.join(HERE, 'n1c1_C_graph.sing'), 'w').write('\n'.join(L) + '\n')
print('written n1c1_C_graph.sing; syzygy directions', tot, '; point', vals)
