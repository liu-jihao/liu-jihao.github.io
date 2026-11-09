# Light bookkeeping for PLAN.md (no CAS). Notation of PLAN.md:
# coordinates w (deg 1, sigma -), y1 (deg 2, +), y2, y3 (deg 2, -); y0 = w^2.
# generators v0=1 (deg 0,+), v1,v2 (deg 3,-), v3,v4 (deg 3,+), v5 (deg 4,-); relation signs = -generator signs.
from itertools import product
from collections import Counter, defaultdict
c = [0,3,3,3,3,4]; dl = [1,-1,-1,1,1,-1]; r = [8-x for x in c]
def mons(n, s):
    """monomials w^a y1^b y2^e y3^f of weight n=a+2(b+e+f) and sigma-sign s=(-1)^(a+e+f)"""
    out=[]
    for a in range(n+1):
        if (n-a)%2: continue
        m=(n-a)//2
        for b in range(m+1):
            for e in range(m-b+1):
                f=m-b-e
                if (-1)**(a+e+f)==s: out.append((a,b,e,f))
    return out
# lambda ansatz
tot=0; byk=Counter(); byorder=Counter(); rows=[]
for i in range(1,6):
    for j in range(i,6):
        for k in range(6):
            n=c[i]+c[j]-c[k]; s=dl[i]*dl[j]*dl[k]
            ms=mons(n,s); tot+=len(ms); byk[k]+=len(ms)
            for m in ms: byorder[m[0]]+=1
            rows.append((i,j,k,n,s,len(ms)))
print("lambda coefficients total:", tot, " by k:", dict(byk))
print("lambda coefficients by w-order:", dict(sorted(byorder.items())))
# RC equations: monomials of B_ij - sum lambda B_0k : weight 8+c_i+c_j, sign -d_i d_j
eqo=Counter(); eqt=0
for i in range(1,6):
    for j in range(i,6):
        ms=mons(8+c[i]+c[j], -dl[i]*dl[j]); eqt+=len(ms)
        for m in ms: eqo[m[0]]+=1
print("RC coefficient equations (monomial slots, i<=j):", eqt, " by w-order:", dict(sorted(eqo.items())))
# trivial syzygies: lambda_{ij,.} -> lambda_{ij,.} + sum_{m=1..5} c_m A_{m,.}; deg c_m = c_i+c_j-r_m, sign d_i d_j * (sign of row m) ... 
# row m sign eps_m=-d_m; A_{mk} has sign eps_m d_k; need sign(c_m)*eps_m*d_k = d_i d_j d_k -> sign(c_m)= d_i d_j / eps_m = -d_i d_j d_m
syz=0
for i in range(1,6):
    for j in range(i,6):
        for m in range(1,6):
            n=c[i]+c[j]-r[m]
            if n<0: continue
            syz+=len(mons(n,-dl[i]*dl[j]*dl[m]))
print("trivial-syzygy directions:", syz)
# A-unknowns per w-order for each normal form (kept q-monomials as in PLAN.md)
def wo(mon):  # mon as string of y0,y1,y2,y3 exponents (k0,k1,k2,k3); order of w*q-term = 1+2*k0 ; of y0*G term = 2+2*k0
    return mon
G=[(k0,k1,k2,k3) for k0 in range(4) for k1 in range(4) for k2 in range(4) for k3 in range(4) if k0+k1+k2+k3==3 and (k2+k3)%2==1]
print("G monomials (anti-invariant cubics in y0..y3):", len(G))
gord=Counter(2+2*m[0] for m in G); print("  G coefficient w-orders (A_00 = y0*G):", dict(sorted(gord.items())))
names={'y0':0,'y1':1,'y2':2,'y3':3}
def parse(s):
    e=[0,0,0,0]
    for t in s.split('*'):
        if '^' in t: v,p=t.split('^'); e[names[v]]+=int(p)
        else: e[names[t]]+=1
    return tuple(e)
kept={
 'N1': {'q1':['y1*y2','y1*y3'],'q2':['y0*y2','y0*y3','y1*y2','y1*y3'],'q3':['y0^2','y0*y1','y1^2'],'q4':['y0^2','y2*y3','y3^2']},
 'N2': {'q1':[],'q2':['y0*y2','y0*y3','y1*y2','y1*y3'],'q3':['y0^2','y0*y1','y1^2','y2^2','y2*y3','y3^2'],'q4':['y0^2','y3^2']},
 'N3': {'q1':[],'q2':['y0*y2','y0*y3'],'q3':['y0^2','y0*y1','y1^2','y2^2','y2*y3','y3^2'],'q4':['y0^2','y2^2','y2*y3','y3^2']},
}
for nf,d in kept.items():
    o=Counter()
    for q,ms in d.items():
        for s in ms: o[1+2*parse(s)[0]]+=1
    print(nf, "kept q-coefficients:", sum(o.values()), " by w-order:", dict(sorted(o.items())))
# lambda table
print("\n(i,j,k): weight, sign, #monomials")
for t in rows: print(t)
