"""Second finite computation of the paper. Lists, from the conditions in the proofs, the parameter
tuples of the three tables of the appendix, compares them with the rows in tables.json (the tables as
printed), and recomputes the rows from the definitions with exact rational arithmetic.
Usage: python3 check_tables.py tables.json"""
import json,sys
from fractions import Fraction as F
from math import gcd,ceil,floor
D=json.load(open(sys.argv[1]))
W=(139754,122195,98574,75044,51760,31572,15936,7272)
def hj(r,q):
    out=[]
    while q:
        e=(r+q-1)//q; out.append(e); r,q=q,e*q-r
    return out
def chain_data(ws):
    n=len(ws)
    P=[1]; a,b=0,1
    for w in ws: a,b=b,w*b-a; P.append(b)
    Q=[1]; a,b=0,1
    for w in reversed(ws): a,b=b,w*b-a; Q.append(b)
    Q=Q[::-1]; R=P[n]
    A=[F(P[i]+Q[i+1],R) for i in range(n)]
    for i,w in enumerate(ws):
        assert w*A[i]==(A[i-1] if i else 1)+(A[i+1] if i<n-1 else 1)
    b=[1-x for x in A]; k=[w-2 for w in ws]
    gamma=n-sum(p*q for p,q in zip(k,b))
    def c(m):
        f=[m*x-floor(m*x) for x in b]
        q=sum(w*z*z for w,z in zip(ws,f))-2*sum(f[i]*f[i+1] for i in range(n-1))
        return (sum(p*z for p,z in zip(k,f))-q)/2
    return R,n,57820*(gamma+1-F(1,R))+sum(w*c(m) for w,m in zip(W,range(2,10)))
def data(e,r,q,n,t):
    return chain_data((hj(r,q) if r>1 else [])[::-1]+[e]+(hj(n,n-t) if n>1 else []))
def T1(e,r,q,n,t):
    R,l,C=data(e,r,q,n,t); v=2*R*(C-91*(l+1)); assert v.denominator==1; return R,int(v)
def T2(e,r,q,n,t):
    R,l,C=data(e,r,q,n,t); v=2*R*(C-500*(R-1)); assert v.denominator==1; return R,int(v)
A=set();B=set();small=set()
for e in range(3,19):
    for r in range(1,19):
        if 5*(e-2)*r>92: continue
        for q in ([0] if r==1 else [q for q in range(1,r) if gcd(r,q)==1]):
            al=(e-1)*r-q
            if al<=9:
                for t in range(1,al+1):
                    for v in range(t):
                        if gcd(t,v)==1 and (v>0 or t==1):
                            A.add((e,r,q,t,v)); h=max(0,ceil(F(r-v,t))-1)
                            while al*((h+1)*t+v)+r*t<20:
                                small.add((e,r,q,t,(h+1)*t+v)); h+=1
            else:
                for t in range(1,min(al,9)+1):
                    for n in range(max(r,t),(r*(46-5*t))//(5*al-46)+1):
                        if gcd(n,t)==1 and (n!=t or n==1):
                            if al*n+r*t<20: small.add((e,r,q,t,n))
                            elif (e,r,q,t,n)!=(14,1,0,1,2): B.add((e,r,q,t,n))
small.add((14,1,0,1,2))
ok=True
def need(c,msg):
    global ok
    if not c: ok=False; print('FAIL',msg)
# Table A1 of the paper: one row for each tuple of order at most 19, and the tuple (14,1,0,1,2)
rows={tuple(x[:5]):x[5:] for x in D['A1']}
need(set(rows)==small,'table A1 tuples')
for k,(R,T) in rows.items():
    r,v=T1(k[0],k[1],k[2],k[4],k[3]); need(r==R and v==T and v>=0,('A1',k))
# Table A2: 450 headers, consecutive intervals from the first h with order at least 20, the last one unbounded
H={tuple(h['header']):h['rows'] for h in D['A2']}
need(set(H)==A,'table A2 headers')
pts=0
for (e,r,q,t,v),rs in H.items():
    al=(e-1)*r-q; start=max(0,ceil(F(r-v,t))-1)
    while al*((start+1)*t+v)+r*t<20: start+=1
    pos=start
    for lo,hi,c0,c1,c2 in rs:
        need(lo==pos,('A2 gap',e,r,q,t,v))
        pos=None if hi is None else hi+1
        need(min(c0,c1 or 0,c2)>=0,('A2 sign',e,r,q,t,v))
        top=8 if hi is None else min(hi-lo,8)
        for z in range(top+1):
            n=(lo+z+1)*t+v; pts+=1
            need(T2(e,r,q,n,t)[1]==c0+(c1 or 0)*z+c2*z*(z-1),('A2',e,r,q,t,v,lo+z))
    need(pos is None,('A2 end',e,r,q,t,v))
# Table B2: every tuple with alpha at least 10 and order at least 20, except (14,1,0,1,2)
ent={tuple(x[:5]):x[5] for x in D['B2']}
need(set(ent)==B,'table B2 tuples')
for k,T in ent.items():
    v=T2(k[0],k[1],k[2],k[4],k[3])[1]; need(v==T and v>=0,('B2',k))
print(json.dumps({'A1_rows':len(rows),'A2_headers':len(H),'A2_rows':sum(len(v) for v in H.values()),'A2_points':pts,'B2_entries':len(ent),'all_match':ok}))
sys.exit(0 if ok else 1)
