# v7_torus.py -- diagonal torus of the order-<=2 normal forms (case nets T6 (iv)), exact rational linear algebra.
# Unknowns: G0..G5 (log-weights of P = diag), OM (w), S1, S2, S3 (y1, y2, y3), dD, dc.
from fractions import Fraction as F
names = ['G0','G1','G2','G3','G4','G5','OM','S1','S2','S3','dD','dc']
idx = {n:i for i,n in enumerate(names)}
def row(d):
    r = [F(0)]*len(names)
    for k,v in d.items(): r[idx[k]] += F(v)
    return r
def term(j,k,a,beta,rhs=None):
    # term w^a y1^b1 y2^b2 y3^b3 in entry (j,k): G_j + G_k + a*OM + beta.S = rhs (0, dD or dc)
    d = {}
    d['G%d'%j] = d.get('G%d'%j,0)+1; d['G%d'%k] = d.get('G%d'%k,0)+1
    d['OM'] = d.get('OM',0)+a
    for i,b in enumerate(beta):
        d['S%d'%(i+1)] = d.get('S%d'%(i+1),0)+b
    if rhs: d[rhs] = d.get(rhs,0)-1
    return row(d)
def rank(M):
    M = [r[:] for r in M]; rk = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        p = next((i for i in range(rk,len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[rk],M[p] = M[p],M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c]/M[rk][c]
                M[i] = [x - f*y for x,y in zip(M[i],M[rk])]
        rk += 1
    return rk
y1,y2,y3 = (1,0,0),(0,1,0),(0,0,1)
def sq(a,b): return tuple(x+y for x,y in zip(a,b))
N1 = [term(1,1,0,y3,'dD'), term(1,2,0,y2), term(2,2,0,y3), term(1,3,0,y1), term(2,4,0,y1), term(3,3,0,y3),
      term(3,4,0,y2), term(4,4,0,y3,'dD'), term(2,3,2,(0,0,0),'dc'),
      term(0,5,0,sq(y2,y2)), term(0,5,0,sq(y3,y3),'dD'), term(0,5,0,sq(y1,y1)), term(1,5,1,(0,0,0))]
N2 = [term(1,1,0,y2), term(2,2,0,y3), term(1,3,0,y1), term(2,4,0,y1), term(3,3,0,y3), term(4,4,0,y2),
      term(2,3,2,(0,0,0),'dc'), term(0,5,0,sq(y2,y3)), term(0,5,0,sq(y1,y1)), term(1,5,1,(0,0,0)), row({'dD':1})]
N3 = [term(1,1,0,y3), term(1,2,0,y2), term(1,3,0,y1), term(2,4,0,y1), term(3,4,0,y2), term(4,4,0,y3),
      term(2,3,2,(0,0,0),'dc'), term(0,5,0,sq(y2,y2)), term(0,5,0,sq(y1,y1)), term(1,5,1,(0,0,0)), row({'dD':1})]
grading = {'G0':-4,'G1':-1,'G2':-1,'G3':-1,'G4':-1,'G5':0,'OM':1,'S1':2,'S2':2,'S3':2,'dD':0,'dc':0}
scal = {'G0':1,'G5':-1,'OM':1,'dc':2}
def satisfies(M, v):
    vec = [F(v.get(n,0)) for n in names]
    return all(sum(a*b for a,b in zip(r,vec)) == 0 for r in M)
for nm,M in (('N1(D)',N1),('N2',N2),('N3',N3)):
    r = rank(M)
    rc = rank(M + [row({'dc':1})])
    print(f'{nm}: conditions {len(M)}, rank {r}, torus dimension {len(names)-r}; with c fixed (dc = 0): {len(names)-rc}; '
          f'grading solves: {satisfies(M, grading)}; scaling s (dc = 2) solves: {satisfies(M, scal)}')
