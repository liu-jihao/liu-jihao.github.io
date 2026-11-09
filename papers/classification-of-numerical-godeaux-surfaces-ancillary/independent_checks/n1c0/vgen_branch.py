#!/usr/bin/env python3
"""Writes the three branch scripts v5_branch_i.sing, v6_branch_ii.sing, v7_branch_iii.sing (verification of n1c0).
Each starts from v_state2.txt (my own elimination: v2 lambda-only, v4 with A-unknowns, D and b8 never pivots),
imposes the branch condition, and continues the constant-pivot elimination; for (ii), (iii) b8 is inverted through
a variable t with t*b8 = 1: a pivot coefficient may be c*b8^a*t^b (a unit), all polynomials are kept in normal form
modulo t*b8 - 1, and an equation divisible by b8 or t is divided (units).  Every unknown is tracked."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(HERE, 'v_data.json')))
F2 = data['F2']
F3rem = data['F3rem']
IB8, IT = 13, 399

def tr(name):
    """index in TR of an unknown: TR[i] = var(4+i); D = var 5, b(k) = 5+k, g(k) = 17+k, L(n) = 27+n"""
    if name == 'D':
        return 1
    if name.startswith('b('):
        return 1 + int(name[2:-1])
    if name.startswith('g('):
        return 13 + int(name[2:-1])
    if name.startswith('L('):
        return 23 + int(name[2:-1])

def subL(s):
    """replace L(n) in a remainder string by the tracked image TR[23+n]"""
    import re
    return re.sub(r'L\((\d+)\)', lambda m: 'TR[%d]' % (23 + int(m.group(1))), s)

COMMON_HEAD = r'''system("--ticks-per-sec", 1000);
int T0 = rtimer;
< "vring.sing";
execute(read("v_state2.txt"));
int i; int k; int m; int ix; int best; int npiv = 0; int pass = 0; int found; int np; int ndiv; int nunit = 0;
ideal E = E2; ideal TR = TR2; intvec elim = ELIM2;
ideal Gt = std(t*b(8) - 1);
proc isunit(poly d, int loc)
{
  if (d == 0) { return(0); }
  if (loc == 0) { if (deg(d) == 0) { return(1); } return(0); }
  if (size(d) != 1) { return(0); }
  intvec e = leadexp(d); int ii;
  for (ii = 1; ii <= size(e); ii++) { if (e[ii] != 0) { if (ii != 13) { if (ii != 399) { return(0); } } } }
  return(1);
}
'''

LOOP = r'''intvec allowed = 0:nvars(basering);
for (i = 6; i <= 398; i++) { if (elim[i] == 0) { allowed[i] = 1; } }
allowed[13] = 0;
poly f; poly dfd; poly bestd; ideal V; poly s; intvec ee; poly dinv;
found = 1;
while (found)
{
  found = 0; pass++; np = 0; ndiv = 0;
  for (k = ncols(E); k >= 1; k--)
  {
    f = E[k];
    if (f != 0)
    {
      V = variables(f); best = 0;
      for (m = 1; m <= ncols(V); m++)
      {
        ix = univariate(V[m]);
        if (ix > 0) { if (allowed[ix] == 1)
        {
          dfd = diff(f, V[m]);
          if (isunit(dfd, LOC)) { if (ix > best) { best = ix; bestd = dfd; } }
        } }
      }
      if (best > 0)
      {
        ee = leadexp(bestd);
        dinv = (1/leadcoef(bestd)) * t^(ee[13]) * b(8)^(ee[399]);
        if (LOC == 0) { dinv = 1/leadcoef(bestd); }
        s = -(f - bestd*var(best))*dinv;
        if (LOC) { s = reduce(s, Gt); }
        E = subst(E, var(best), s); TR = subst(TR, var(best), s);
        if (LOC) { E = reduce(E, Gt); TR = reduce(TR, Gt); }
        allowed[best] = 0; elim[best] = 1; npiv++; np++; found = 1;
        if (deg(bestd) > 0) { nunit++; "  unit pivot " + string(var(best)) + ", coefficient " + string(bestd); }
      }
    }
  }
  if (LOC)
  {
    for (k = 1; k <= ncols(E); k++)
    {
      if (E[k] != 0) { while (subst(E[k], b(8), 0) == 0) { E[k] = E[k] / b(8); ndiv++; } }
      if (E[k] != 0) { while (subst(E[k], t, 0) == 0) { E[k] = E[k] / t; ndiv++; } }
    }
  }
  E = simplify(E, 2);
  "pass " + string(pass) + ": pivots " + string(np) + ", divisions by b8 or t " + string(ndiv) + ", nonzero equations left " + string(size(E));
  if (ndiv > 0) { found = 1; }
}
int nel = 0; int nelA = 0; string sfree = ""; int nfl = 0;
for (i = 5; i <= 398; i++) { if (elim[i] == 1) { nel++; if (i <= 27) { nelA++; } } else { if (i <= 27) { sfree = sfree + string(var(i)) + " "; } else { nfl++; } } }
"pivots with a non-constant (unit) coefficient: " + string(nunit);
"eliminated or fixed by the branch condition: " + string(nel) + " (A-unknowns " + string(nelA) + ", lambda " + string(nel - nelA) + ")";
"free A-unknowns: " + sfree + " (" + string(23 - nelA) + "); free lambda: " + string(nfl);
ideal EB = simplify(E, 2);
"leftover equations: " + string(size(EB));
for (k = 1; k <= ncols(EB); k++) { if (EB[k] != 0) { list FF = factorize(EB[k]); "  LEFT " + string(k) + ": " + string(FF[1]) + "  mult " + string(FF[2]); kill FF; } }
"images of the 23 A-unknowns:";
for (i = 1; i <= 23; i++) { "  " + string(var(i+4)) + " -> " + string(TR[i]); }
int lamA = 0; ideal VA;
for (i = 1; i <= 23; i++) { if (TR[i] != 0) { VA = variables(TR[i]); for (k = 1; k <= ncols(VA); k++) { if (univariate(VA[k]) >= 28) { lamA++; } } } }
"graph check: lambda variables occurring in the images of the A-unknowns: " + string(lamA);
int occ = 0; ideal VE;
for (i = 1; i <= ncols(TR); i++) { if (TR[i] != 0) { VE = variables(TR[i]); for (k = 1; k <= ncols(VE); k++) { ix = univariate(VE[k]); if (ix > 0) { if (ix <= 398) { if (elim[ix] == 1) { occ++; } } } } } }
"eliminated variables occurring in tracked images: " + string(occ);
'''

def f2_block():
    lines = ['"--- F2 point (1:0:...:0): b7, b10 and the pure w-power coefficients of the lambda_{ij,0} ---";',
             'int f2z = 1;']
    items = [('b7', tr('b(7)')), ('b10', tr('b(10)'))] + [('lam_%d%d_0[w^%d] = L(%d)' % (i, j, a, n), 23 + n) for (n, i, j, a) in F2]
    for (nm, idx) in items:
        lines.append('"  %s -> " + string(TR[%d]); if (TR[%d] != 0) { f2z = 0; }' % (nm, idx, idx))
    lines.append('"all F2 conditions identically 0 after substitution: " + string(f2z);')
    return '\n'.join(lines) + '\n'

def f3_block():
    lines = ['"--- F3 points (0:0:+-sqrt(D):1|0): remainders of lambda_{ij,0}(0,0,y2,y3) modulo y2^2 - D*y3^2 ---";',
             'int f3z = 1; poly rr;']
    for (ij, r1, r0) in F3rem:
        for nm, r in (('r1', r1), ('r0', r0)):
            lines.append('rr = reduce(%s, Gt); "  lambda_%d%d,0 %s -> " + string(rr); if (rr != 0) { f3z = 0; }' % (subL(r), ij[0], ij[1], nm))
    lines.append('"all F3 remainders identically 0 after substitution: " + string(f3z);')
    return '\n'.join(lines) + '\n'

def write(fname, label, subst_lines, loc, extra):
    S = ['// %s: branch %s of case n1c0 (verification). Generated by vgen_branch.py.' % (fname, label), COMMON_HEAD,
         'int LOC = %d;' % loc]
    S += subst_lines
    S.append(LOOP)
    S.append(f2_block())
    S.append(f3_block())
    S.append(extra)
    S.append('"total time (ms): " + string(rtimer - T0) + ", memory(2) = " + string(memory(2));')
    S.append('quit;')
    open(os.path.join(HERE, fname), 'w').write('\n'.join(S) + '\n')

write('v5_branch_i.sing', '(i) b3 = b8 = 0',
      ['E = subst(E, b(3), 0); E = subst(E, b(8), 0); TR = subst(TR, b(3), 0); TR = subst(TR, b(8), 0);',
       'elim[8] = 1; elim[13] = 1;', 'E = simplify(E, 2);'], 0, '')
write('v6_branch_ii.sing', '(ii) b3 = b8, b8 invertible',
      ['E = subst(E, b(3), b(8)); TR = subst(TR, b(3), b(8));', 'elim[8] = 1;', 'E = simplify(E, 2);'], 1,
      'write(":w v_param_ii.txt", "ideal TRii = " + string(TR) + ";");\nwrite(":a v_param_ii.txt", "intvec ELii = " + string(elim) + ";");\n')
write('v7_branch_iii.sing', '(iii) b3 = -b8, b8 invertible',
      ['E = subst(E, b(3), -b(8)); TR = subst(TR, b(3), -b(8));', 'elim[8] = 1;', 'E = simplify(E, 2);'], 1,
      'write(":w v_param_iii.txt", "ideal TRiii = " + string(TR) + ";");\nwrite(":a v_param_iii.txt", "intvec ELiii = " + string(elim) + ";");\nwrite(":a v_param_iii.txt", "ideal EBiii = " + string(EB) + ";");\n')
print('written')
