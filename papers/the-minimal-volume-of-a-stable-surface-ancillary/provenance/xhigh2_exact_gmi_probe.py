"""Bounded discovery of exact-valid GMI/rounding cuts for the complete RR model.

Run only with [job runner]. Every saved cut has a rational derivation; a solver
infeasibility status is not a certificate. No previous unverified cuts are used.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd, lcm, floor
import sys, json, time
import numpy as np
from scipy.sparse import csr_matrix
R = Path(__file__).resolve().parent
sys.path.insert(0, str(R / 'xhigh2_highspy_deps'))
import highspy as H

P = int(sys.argv[1]) if len(sys.argv) > 1 else 45
assert 18<=P<=221
LIMIT = 240
start = time.monotonic()
rows = json.loads((R / 'xhigh2_complete_304.json').read_text())
n = len(rows); nc = n + P - 17
cc = [list(map(F, r['c'])) for r in rows]
costs = [F(r['certificate_cost']) for r in rows]
upper = [min(833, (151733*c.denominator-1)//c.numerator,
             1667//(r['ell']+1)) for r,c in zip(rows,costs)] + [1 if m<=156 else 2 for m in range(18,P+1)]
# Keep even RR-zero columns. Hence no zero-packet assumptions are needed here.
AA=[]; bb=[]; equality=[]; names=[]
def add(a,b,eq,name):
    a={j:F(v) for j,v in a.items() if v}; b=F(b)
    scale=max([abs(b),F(1)]+[abs(v) for v in a.values()])
    AA.append({j:v/scale for j,v in a.items()});bb.append(b/scale)
    equality.append(eq);names.append(name)
ss=[[c[m-2]-m*(m-1)//2*c[0] for m in range(2,P+1)] for c in cc]
assert all(v.denominator==1 for s in ss for v in s)
add({j:s[1] for j,s in enumerate(ss)},2,True,'s3')
for m in range(4,P+1):
    a={j:s[m-2]-2*s[m-3]+s[m-4] for j,s in enumerate(ss)}
    for k,v in ((m,-1),(m-1,2),(m-2,-1)):
        if k>=18:a[n+k-18]=v
    add(a,1,True,('d2',m))
add({j:c[0] for j,c in enumerate(cc)},-1-F(1,84235998),False,'volume_lower')
add({j:-c[0] for j,c in enumerate(cc)},1+F(1,48983),False,'volume_upper')
add({j:F(r['gamma'])-c[0] for j,(r,c) in enumerate(zip(rows,cc))},10,False,'rank')
add({j:F(r['gamma'])+1-1/F(r['group'])-F(4,3)*c[0]
     for j,(r,c) in enumerate(zip(rows,cc))},F(40,3),False,'BMY')
add({j:r['ell']+1 for j,r in enumerate(rows)},1667,False,'length')
add(dict(enumerate(costs)),151733,False,'cost')
for a in range(18,P+1):
    for k in range(2,P//a+1):
        add({n+a-18:1,n+k*a-18:-1},0,False,('multiple',a,k*a))
    for b in range(a+1,P+1):
        if gcd(a,b)<=17 and lcm(a,b)<=156:
            add({n+a-18:1,n+b-18:1},1,False,('gcd',a,b))
base_count=len(AA)
if len(sys.argv)>2 and sys.argv[2]=='export':
    # Extra actual multiplication inequalities used only by this new export.
    for a in range(18,P+1):
        for b in range(a,P-a+1):
            ub=1 if b<=156 else 2
            ent={n+a-18:ub,n+b-18:1,n+a+b-18:-1}
            if a==b:ent[n+a-18]=ub+1
            add(ent,ub,False,('hopf',a,b))
    out=['NAME RR45','ROWS',' N OBJ']
    out.extend((' E ' if eq else ' L ')+f'R{i}' for i,eq in enumerate(equality))
    out+=['COLUMNS'," MARK0 'MARKER' 'INTORG'"]
    for j in range(nc):
        if j<n:out.append(f' X{j} OBJ {costs[j]/151733}')
        for i,a in enumerate(AA):
            if j in a:out.append(f' X{j} R{i} {a[j]}')
    out+=[" MARK1 'MARKER' 'INTEND'",'RHS']
    out.extend(f' RHS1 R{i} {b}' for i,b in enumerate(bb) if b)
    out+=['BOUNDS']
    for j,u in enumerate(upper):
        out.extend([f' LO BND X{j} 0',f' UP BND X{j} {u}'])
    out+=['ENDATA']
    path=R/f'xhigh2_exact_prefix_{P}.mps';path.write_text('\n'.join(out)+'\n')
    meta={'prefix':P,'local_rows':n,'variables':nc,'constraints':len(AA),
          'names':names,'upper':upper,'all_input_coefficients_exact_rational':True}
    (R/f'xhigh2_exact_prefix_{P}_meta.json').write_text(json.dumps(meta)+'\n')
    print('exported',path,'variables',nc,'constraints',len(AA),flush=True);raise SystemExit
if len(sys.argv)>2 and sys.argv[2]=='audit':
    # Replay the archive using exact arithmetic only. This branch invokes no
    # LP solver and does not trust the discovery basis, signs, or rounding.
    data=json.loads((R/f'xhigh2_gmi_{P}.json').read_text())
    checked=0
    for cut in data['cuts']:
        parents=cut['parents'];assert base_count<=parents<=len(AA)
        yy={i:F(y) for i,y in cut['y']};flip=set(cut['flip'])
        coeff=[F(0)]*nc;rhs=F(0)
        for i,y in yy.items():
            assert i<parents
            rhs+=y*bb[i]
            for j,a in AA[i].items():coeff[j]+=y*a
        rhs-=sum(coeff[j]*upper[j] for j in flip)
        fraction=rhs-floor(rhs);assert fraction==F(cut['fraction']) and 0<fraction<1
        gi=[]
        for j,a in enumerate(coeff):
            z=-a if j in flip else a;f=z-floor(z)
            gi.append(f/fraction if f<=fraction else (1-f)/(1-fraction))
        lhs=[g if j in flip else -g for j,g in enumerate(gi)]
        bound=-1+sum(gi[j]*upper[j] for j in flip)
        for i,y in yy.items():
            if equality[i]:continue
            h=y/fraction if y>=0 else -y/(1-fraction)
            bound+=h*bb[i]
            for j,a in AA[i].items():lhs[j]+=h*a
        scale=F(cut['round_scale']);assert scale>0
        recorded=dict(cut['a']);ints=[recorded.get(j,0) for j in range(nc)]
        assert all(z==floor(a*scale+F(1,2)) for z,a in zip(ints,lhs))
        rounded=scale*bound+sum(max(F(0),F(z)-a*scale)*u for z,a,u in zip(ints,lhs,upper))
        assert rounded==F(cut['round_bound']) and floor(rounded)==cut['b']
        add(recorded,cut['b'],False,('gmi',checked));checked+=1
        if checked%25==0:print('rationally audited',checked,flush=True)
        if time.monotonic()-start>480:break
    result={'checked':checked,'total':len(data['cuts']),'all_valid':checked==len(data['cuts']),
            'infeasibility_proved':False,'seconds':time.monotonic()-start}
    (R/f'xhigh2_gmi_audit_{P}.json').write_text(json.dumps(result)+'\n')
    print(result,flush=True);raise SystemExit
def sparse(a):
    return csr_matrix(np.array([[float(row.get(j,0)) for j in range(nc)] for row in a]))
mat=sparse(AA).tocsc()
lp=H.HighsLp();lp.num_col_=nc;lp.num_row_=len(AA)
lp.col_cost_=[float(c)/151733 for c in costs]+[0.]*(P-17)
lp.col_lower_=[0.]*nc;lp.col_upper_=upper
lp.row_lower_=[float(b) if e else -H.kHighsInf for b,e in zip(bb,equality)]
lp.row_upper_=[float(b) for b in bb]
lp.a_matrix_.format_=H.MatrixFormat.kColwise
lp.a_matrix_.start_=mat.indptr;lp.a_matrix_.index_=mat.indices;lp.a_matrix_.value_=mat.data
h=H.Highs();h.setOptionValue('output_flag',False);h.setOptionValue('threads',1)
h.setOptionValue('solver','simplex');h.setOptionValue('presolve','off')
h.setOptionValue('primal_feasibility_tolerance',1e-8)
h.passModel(lp)
archive=[]; seen=set(); end='budget'; witnesses=[]
def cut_from(y,flip,x):
    # Sum y_i(A_i*x+s_i)=sum y_i*b_i. Equality slacks are zero.
    a=[F(0)]*nc; b=sum(v*t for v,t in zip(y,bb))
    for v,row in zip(y,AA):
        if v:
            for j,t in row.items():a[j]+=v*t
    b-=sum(a[j]*upper[j] for j in flip)
    f=b-floor(b)
    if f<F(1,100) or f>F(99,100):return None
    g=[]
    for j,v in enumerate(a):
        if j in flip:v=-v
        fj=v-floor(v)
        g.append(fj/f if fj<=f else (1-fj)/(1-f))
    hh=[F(0) if eq else v/f if v>=0 else -v/(1-f)
        for v,eq in zip(y,equality)]
    out=[g[j] if j in flip else -g[j] for j in range(nc)]
    rhs=sum(v*t for v,t in zip(hh,bb))+sum(g[j]*upper[j] for j in flip)-1
    for v,row in zip(hh,AA):
        if v:
            for j,t in row.items():out[j]+=v*t
    # Deliberately rounded integer coefficients, with exact bound compensation.
    mag=max([F(1)]+[abs(z) for z in out]); scale=F(10**11)/mag
    ints=[floor(z*scale+F(1,2)) for z in out]
    bound=rhs*scale+sum(max(F(0),F(z)-v*scale)*u for z,v,u in zip(ints,out,upper))
    ib=floor(bound)
    violation=sum(z*t for z,t in zip(ints,x))-ib
    norm=max(1,max(map(abs,ints)))
    if violation/norm<1e-6:return None
    return {'y':[[i,str(v)] for i,v in enumerate(y) if v], 'flip':sorted(flip),
            'fraction':str(f),'round_scale':str(scale),'round_bound':str(bound),
            'a':[[j,v] for j,v in enumerate(ints) if v], 'b':ib,
            'violation':violation/norm,'parents':len(AA)}

if len(sys.argv)>2 and sys.argv[2]=='refine':
    prev=json.loads((R/f'xhigh2_gmi_{P}.json').read_text())
    for c in prev['cuts']:
        archive.append(c);add(dict(c['a']),c['b'],False,('gmi',len(archive)-1))
        idx=list(AA[-1]);vals=[float(AA[-1][j]) for j in idx]
        h.addRow(-H.kHighsInf,float(bb[-1]),len(idx),np.array(idx,dtype=np.int32),np.array(vals))
        seen.add((tuple(tuple(t) for t in c['a']),c['b']))

for sweep in range(20):
    if time.monotonic()-start>LIMIT:break
    h.setOptionValue('time_limit',max(1,min(30,LIMIT-(time.monotonic()-start))))
    h.run();st=h.getModelStatus()
    print('round',sweep,'status',str(st),'cuts',len(archive),'seconds',round(time.monotonic()-start,2),flush=True)
    if st==H.HighsModelStatus.kInfeasible:
        end='floating_infeasible';break
    if st!=H.HighsModelStatus.kOptimal:end=str(st);break
    sol=h.getSolution();x=np.array(sol.col_value)
    bas=h.getBasis()
    basic=[j for j,s in enumerate(bas.col_status) if s==H.HighsBasisStatus.kBasic]
    slack=[i for i,s in enumerate(bas.row_status) if s==H.HighsBasisStatus.kBasic]
    afloat=np.array([[float(row.get(j,0)) for j in range(nc)] for row in AA])
    basis_matrix=np.column_stack([afloat[:,basic],np.eye(len(AA))[:,slack]])
    assert basis_matrix.shape==(len(AA),len(AA))
    inv=np.linalg.inv(basis_matrix)
    flip={j for j,s in enumerate(bas.col_status) if s==H.HighsBasisStatus.kUpper}
    choices=sorted([(abs(x[int(j)]-round(x[int(j)])),i,int(j))
                    for i,j in enumerate(basic) if j>=0],reverse=True)
    pending=[]
    for frac,i,j in choices[:20]:
        if frac<1e-6:continue
        yy=inv[i,:]
        y=[F(float(z)).limit_denominator(10**12) for z in yy]
        candidate=cut_from(y,flip,x)
        if candidate:
            key=(tuple(tuple(t) for t in candidate['a']),candidate['b'])
            if key not in seen:pending.append(candidate);seen.add(key)
        if time.monotonic()-start>LIMIT:break
    if not pending:
        print('fractional_basic',choices[:5],flush=True)
        end='no_violated_cut';witnesses=[(int(j),float(t)) for j,t in enumerate(x) if t>1e-8];break
    for c in pending:
        archive.append(c);add(dict(c['a']),c['b'],False,('gmi',len(archive)-1))
        idx=list(AA[-1]);vals=[float(AA[-1][j]) for j in idx]
        h.addRow(-H.kHighsInf,float(bb[-1]),len(idx),np.array(idx,dtype=np.int32),np.array(vals))
    print('added',len(pending),'max_violation',max(c['violation'] for c in pending),flush=True)
    (R/f'xhigh2_gmi_{P}.json').write_text(json.dumps({'P':P,'base_count':base_count,'cuts':archive,'status':'ongoing'})+'\n')

# If the new LP looks infeasible, seek and rationally validate a Farkas ray.
certificate=None
if end=='floating_infeasible':
    status,exists,ray=h.getDualRay()
    if exists:
        for sign in (1,-1):
            y=[F(float(sign*z)).limit_denominator(10**8) for z in ray]
            if any(v>0 and not e for v,e in zip(y,equality)):continue
            rhs=sum(v*b for v,b in zip(y,bb));res=[F(0)]*nc
            for v,row in zip(y,AA):
                if v:
                    for j,a in row.items():res[j]-=v*a
            lower=rhs+sum(v*u for v,u in zip(res,upper) if v<0)
            if lower>0:
                certificate={'y':[[i,str(v)] for i,v in enumerate(y) if v], 'margin':str(lower)}
                end='exact_infeasible';break
out={'P':P,'base_count':base_count,'cuts':archive,'status':end,'farkas':certificate,
     'witness':witnesses,'seconds':time.monotonic()-start}
(R/f'xhigh2_gmi_{P}.json').write_text(json.dumps(out)+'\n')
print('RESULT',end,'cuts',len(archive),'seconds',round(time.monotonic()-start,2),flush=True)
