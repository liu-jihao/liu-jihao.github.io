# v9_monomials.py -- all monomials w^a y1^b y2^e y3^f of the weight 8 - c_j - c_k and sigma-sign -delta_j delta_k
# of each entry A_jk (weights w:1, y:2; sigma: w,y2,y3 -> -w,-y2,-y3), compared with the ansatz of case_common.md 1.
c = [0,3,3,3,3,4]; dl = [1,-1,-1,1,1,-1]
def mons(wt, sign):
    out = []
    for a in range(wt+1):
        if (wt - a) % 2: continue
        d = (wt - a)//2
        for b in range(d+1):
            for e in range(d-b+1):
                f = d-b-e
                if (-1)**(a+e+f) == sign: out.append((a,b,e,f))
    return out
def name(m):
    a,b,e,f = m; s = []
    if a: s.append('w^%d'%a if a>1 else 'w')
    for v,k in (('y1',b),('y2',e),('y3',f)):
        if k: s.append(v+('^%d'%k if k>1 else ''))
    return '*'.join(s) or '1'
for j in range(6):
    for k in range(j,6):
        wt = 8 - c[j] - c[k]; sg = -dl[j]*dl[k]
        ms = mons(wt, sg) if wt >= 0 else []
        print(f'A_{j}{k}: weight {wt}, sign {"+" if sg>0 else "-"}: {len(ms)} monomials: ' + ', '.join(name(m) for m in ms))
