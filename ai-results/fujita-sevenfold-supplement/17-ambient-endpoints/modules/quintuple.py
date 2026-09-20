"""Independent cold audit of the COMPLETE native e072 printed certificate.

Run only as: danus compute computation/sub8_e072_full_cold_audit/audit_certificate.py
Uses only standard-library exact rational arithmetic. Does not import, execute,
or consume the producer's scripts, search JSON, or arithmetic summary.
"""
from .common import Q, record_lines, require
from pathlib import Path
from math import comb, ceil
from collections import Counter
from hashlib import sha256
import json, re, time

R, U, V, T = map(Q, ['9/10', '26/25', '27/25', '29/25'])
H0, HQ, EPS = map(Q, ['1323/10000', '2017/10000', '1/1000000'])
CL, CR = map(Q, ['16020581569/90194313216', '8010290789/45097156608'])
DEN = 10**9
ZERO, ONE, X = (Q(0),), (Q(1),), (Q(0), Q(1))

def encode(a):
    if isinstance(a, Q): return str(a)
    if isinstance(a, dict): return {str(k): encode(v) for k, v in a.items()}
    if isinstance(a, (list, tuple)): return [encode(v) for v in a]
    return a

def canonical(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return tuple(p)

def plus(*ps):
    return canonical(sum((p[i] if i < len(p) else Q(0)) for p in ps)
                     for i in range(max(map(len, ps))))

def scale(p, c): return canonical(c * a for a in p)
def minus(a, b): return plus(a, scale(b, -1))

def times(a, b):
    c = [Q(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b): c[i+j] += ai * bj
    return canonical(c)

def power(p, n):
    z = ONE
    for _ in range(n): z = times(z, p)
    return z

def value(p, t):
    z = Q(0)
    for a in reversed(p): z = z * t + a
    return z

def integral7(p, a, b):
    return 7 * sum(c * (b**(i+1)-a**(i+1)) / (i+1)
                   for i, c in enumerate(p))

def positive_power(p, n, t):
    return power(p, n) if value(p, t) > 0 else ZERO

def cuts(lines, a, b):
    endpoints = {a, b}
    for p in lines:
        assert len(p) <= 2, ('nonaffine switch', p)
        if len(p) == 2 and p[1]:
            z = -p[0] / p[1]
            if a < z < b: endpoints.add(z)
    return sorted(endpoints)

def prepared(box):
    al, ah, bl, bh, ll, lh, tl, th, zl, zh, ol, oh = box
    h = max(HQ, H0+al*(U-R)-EPS)
    ah = min(ah, (lh+EPS-H0)/(V-R))
    bh = min(bh, (lh+EPS-h)/(V-U))
    bl = max(bl, h/U)
    d = 1+32*(V-ll)**7
    zl = max(zl, ll-(V-U)*d)
    failed = {'a': al > ah, 'b': bl > bh,
              'chord': h*V/U > lh+EPS, 'z': zl > zh}
    if any(failed.values()): return {'empty': failed}
    assert 0 <= tl < th <= 1 and 0 <= ol < oh <= 1
    assert 0 <= zl <= zh <= Q(1,5) < lh < 1
    ell = max(2, ceil(tl/(1-tl)))
    k2 = ell*((1-tl)+tl/ell)**2
    mu = max(Q(0), h-tl*zh-(1-tl)*lh-EPS)
    A = 1-zl
    mass = max(Q(0), A*A-4*mu*mu/k2)**6/A**5+EPS
    return {'box': (al,ah,bl,bh,ll,lh,tl,th,zl,zh,ol,oh),
            'h': h, 'ell': ell, 'mass': mass, 'mu': mu, 'K2': k2}

def joint_switches(D, F, H, beta, old):
    result = [D,F,H,minus(D,F),minus(D,H),minus(F,H)]
    for j in ((2,3,4,5) if old else (2,)):
        A = minus(D, scale(H, Q(j)/(1-beta)))
        B = j*beta/(1-beta)-1
        result.extend(plus(A, scale(e, B)) for e in (F,D,H))
    return result

def joint(D, F, H, beta, old, mid):
    # Direct symbolic integration of 6*c_j*(D-e-j*h_e)_+^5 de.
    dm, fm, hm = [value(p, mid) for p in (D,F,H)]
    if dm <= 0 or fm >= dm: return ZERO
    if hm <= 0: return positive_power(minus(D,F), 6, mid)
    lower_tail = F if fm >= hm else H
    out = positive_power(minus(D,lower_tail), 6, mid)
    upper = D if dm <= hm else H
    if value(upper,mid) <= fm: return out
    for j, c in (((2,10),(3,-20),(4,15),(5,-4)) if old else ((2,1),)):
        A = minus(D, scale(H, Q(j)/(1-beta)))
        B = j*beta/(1-beta)-1
        if B == 0:
            term = scale(times(positive_power(A,5,mid), minus(upper,F)), 6)
        else:
            top = positive_power(plus(A,scale(upper,B)), 6, mid)
            low = positive_power(plus(A,scale(F,B)), 6, mid)
            term = scale(minus(top,low), 1/B)
        out = plus(out, scale(term,c))
    return out

def selected_polynomial(label, frame, mid):
    F, D, Y0, Y1, R0, R1, ol, tl = frame
    if label.startswith('B'):
        rate = Q(1) if mid<R else Q(61,100) if mid<Q(99,100) else Q(11,25) if mid<U else Q(17,50)
        assert Q(label[1:]) == rate, ('baseline', label, mid)
        return scale(power(X,6),rate)
    if label == 'P': return positive_power(D,6,mid)
    if label in ('J0','JD0','J1','JD1'):
        args = {'J0':(X,F,Y0,ol,True), 'JD0':(D,ZERO,R0,ol,True),
                'J1':(X,F,Y1,tl,False), 'JD1':(D,ZERO,R1,tl,False)}
        return joint(*args[label], mid)
    prefix, name = label.split(':')
    Z, Y, old = {'0':(X,Y0,True), 'D0':(D,R0,True),
                 '1':(X,Y1,False), 'D1':(D,R1,False)}[prefix]
    zm, ym = value(Z,mid), value(Y,mid)
    if zm <= 0:
        assert name in ('0','span'), (label,zm,ym)
        return ZERO
    if ym <= 0:
        assert name in ('B1','span'), (label,zm,ym)
        return power(Z,6)
    c = ym/zm
    if name == 'span':
        if old:
            coeff = [(j,(-1)**(j+1)*(comb(6,j)-5*comb(5,j-1))) for j in range(1,7)]
        else: coeff = [(1,1),(2,1),(3,-1)]
        return plus(*(scale(positive_power(minus(Z,scale(Y,j)),6,mid),a) for j,a in coeff))
    if c >= 1:
        assert name == '0', (label,c)
        return ZERO
    if name == 'B1':
        assert not old and CL < c < CR, (label,c)
        return power(Z,6)
    if name == 'Nlo':
        assert old and c <= Q(1,7), (label,c)
        return minus(power(Z,6), scale(times(Z,power(Y,5)), Q(35,6)**5))
    if name == 'Nhi':
        assert old and c >= Q(1,7), (label,c)
        return minus(power(Z,6),scale(times(Z,power(plus(scale(Z,29),scale(Y,7)),5)),Q(1,36)**5))
    if name == 'C':
        assert old and c >= Q(1,7), (label,c)
        return scale(power(minus(Z,Y),6), Q(7,6)**6)
    if name == 'A':
        assert c >= (Q(1,2) if old else CR), (label,c)
        return power(minus(Z,Y),6)
    if name == 'Q':
        assert not old and c <= CL, (label,c)
        return times(power(minus(Z,scale(Y,2)),5),plus(Z,scale(Y,10)))
    raise AssertionError(('unknown label',label))

def domains(p):
    al,ah,bl,bh,ll,lh,tl,th,zl,zh,ol,oh = p['box']
    h = p['h']; d = 1+32*(V-ll)**7
    anchors = [Q(0),R,Q(99,100),U,V,T]
    ratios = [Q(1,j) for j in (7,6,5,4,3,2,1)]+[Q(2,11),CL,CR]
    for start,end in zip(anchors,anchors[1:]):
        mid=(start+end)/2; a=ah if mid<R else al; b=bh if mid<U else bl
        oldline=(H0-a*R,a); newline=(h-b*U,b)
        f26=(zl-25*(lh-zl)*U,25*(lh-zl)) if mid<U else (Q(0),zl/U)
        f27=(ll-d*V,d) if mid<V else (ll-25*(ll-zh)*V,25*(ll-zh))
        first=cuts([oldline,newline,f26,f27,minus(f26,f27)],start,end)
        for l,r in zip(first,first[1:]):
            m=(l+r)/2
            Y0=oldline if value(oldline,m)>0 else ZERO
            Y1=newline if value(newline,m)>0 else ZERO
            F=max((ZERO,f26,f27),key=lambda a:value(a,m))
            D=minus(X,F); R0=minus(Y0,scale(F,oh)); R1=minus(Y1,scale(F,th))
            switches=[]
            for Z,Y in ((X,Y0),(X,Y1),(D,R0),(D,R1)):
                switches.extend([Z,Y]+[minus(Y,scale(Z,c)) for c in ratios])
            for args in ((X,F,Y0,ol,True),(D,ZERO,R0,ol,True),
                         (X,F,Y1,tl,False),(D,ZERO,R1,tl,False)):
                switches.extend(joint_switches(*args))
            second=cuts(switches,l,r)
            frame=(F,D,Y0,Y1,R0,R1,ol,tl)
            for a,b in zip(second,second[1:]): yield a,b,frame

def validate(key,text):
    started=time.monotonic()
    rows=record_lines(text)
    root_line=next(i for i,x in enumerate(rows) if x.startswith('Root:'))
    root=tuple(map(Q,[x.strip() for x in rows[root_line+1].split('|')]))
    expected=(H0/R,(Q(2,5)+EPS-H0)/(V-R),HQ/U,(Q(2,5)+EPS-HQ)/(V-U),
              Q(7,30),Q(2,5),Q(0),Q(1),Q(0),Q(1,5),Q(0),Q(1))
    assert root==expected
    si=root_line+3
    ei=next(i for i,x in enumerate(rows) if x.startswith('Empty boxes:'))
    li=next(i for i,x in enumerate(rows) if x.startswith('Leaves:'))
    ci=next(i for i,x in enumerate(rows) if x.startswith('Leaf '))
    boxes={'R':root}; split_count=0
    for line in rows[si:ei]:
        tag,axis,mid=[x.strip() for x in line.split('|')]; axis=int(axis); mid=Q(mid)
        assert 0<=axis<6 and tag in boxes
        box=boxes.pop(tag); assert box[2*axis]<mid<box[2*axis+1]; assert mid==(box[2*axis]+box[2*axis+1])/2
        left=list(box); right=list(box); left[2*axis+1]=mid; right[2*axis]=mid
        boxes[tag+'0']=tuple(left); boxes[tag+'1']=tuple(right); split_count+=1
    terminals={}; empty=[]
    for line in rows[ei+1:li]:
        fields=[x.strip() for x in line.split('|')]; assert len(fields)==14
        tag=fields[0]; box=tuple(map(Q,fields[1:13])); reason=fields[13]
        assert tag not in terminals; terminals[tag]=box; p=prepared(box)
        if reason=='interval': assert 'empty' in p
        elif reason=='full26mass': assert 'empty' not in p and p['mass']<Q(3,40)
        else: raise AssertionError(reason)
        empty.append({'tag':tag,'reason':reason,'prepared':p})
    leaves={}
    for line in rows[li+1:ci]:
        f=[x.strip() for x in line.split('|')]; assert len(f)==18
        tag=f[0]; box=tuple(map(Q,f[1:13])); assert tag not in terminals
        terminals[tag]=box; p=prepared(box); assert 'empty' not in p
        assert p['h']==Q(f[13]) and p['ell']==int(f[14]) and p['mass']==Q(f[15])
        assert p['mass']>=Q(3,40)
        leaves[tag]={'prepared':p,'method':f[16],'printed_upper':Q(f[17]),'cells':[]}
    assert boxes==terminals and split_count+1==len(terminals)
    assert sum(Q(1,2**(len(tag)-1)) for tag in terminals)==1
    index=ci; observed=set(); totalcells=0
    while index<len(rows):
        match=re.fullmatch(r'Leaf(R[01]*):left\|right\|polynomial\|uppernumeratorwithdenominator10\^9', re.sub(r'\s+', '', rows[index]))
        assert match, (index+1,rows[index]); tag=match.group(1)
        assert tag in leaves and tag not in observed; observed.add(tag); index+=1
        while index<len(rows) and not rows[index].startswith('Cell sum '):
            a,b,label,num=[x.strip() for x in rows[index].split('|')]
            leaves[tag]['cells'].append((Q(a),Q(b),label,int(num))); index+=1; totalcells+=1
        match=re.fullmatch(r'Cell sum (\d+)\.',rows[index]); assert match
        assert int(match.group(1))==sum(c[3] for c in leaves[tag]['cells']); index+=1
    assert observed==set(leaves)
    # Independent constants and rational bracket checks.
    rho=Q(2237,1000); initial=(R*rho-1)/(7*rho-8)
    assert rho*rho>5 and initial==Q(10133,76590) and initial-H0>EPS
    heights=[]
    for a,p,expected_h in [(Q(5,4),Q(92,125),Q(150,743)),
                           (Q(623,500),Q(719,1000),Q(7396,36025)),
                           (Q(33,25),Q(76,125),Q(233,1155))]:
        assert 7*a+p-8>0
        h=(a*U-1)/(7*a+p-8); assert h==expected_h; heights.append(h)
    assert min(heights)-HQ>EPS
    difference=minus(times(power(minus(ONE,scale(X,2)),5),plus(ONE,scale(X,10))),power(minus(ONE,X),6))
    assert 0<CL<CR<Q(1,2) and value(difference,CL)>0>value(difference,CR)
    retained=[(-1)**j*(j-1)*comb(5,j) for j in range(2,6)]
    assert retained==[10,-20,15,-4]
    methods=Counter(); labels=Counter(); maxima={}; receipts=[]; branches=0; zero_B=Counter()
    for number,(tag,leaf) in enumerate(leaves.items(),1):
        p=leaf['prepared']; method=leaf['method']; assert method in ('29','27mass','26mass')
        endpoint={'29':T,'27mass':V,'26mass':U}[method]
        parts=list(domains(p)); cursor=0; end=Q(0); subtotal=0; leafbranches=0
        exactsum=Q(0); maxdegree=0
        for a,b,label,num in leaf['cells']:
            assert a==end and a<b and num>0; end=b; subtotal+=num; labels[label]+=1
            while parts[cursor][1]<=a: cursor+=1
            j=cursor; covered=a; polynomial=None; cellbranches=0
            while j<len(parts) and parts[j][0]<b:
                l,r,frame=parts[j]; left=max(a,l); right=min(b,r)
                assert left==covered and left<right
                current=selected_polynomial(label,frame,(left+right)/2)
                if polynomial is None: polynomial=current
                assert current==polynomial, (tag,label,a,b,l,r)
                covered=right; cellbranches+=1; j+=1
            assert covered==b
            exact=integral7(polynomial,a,b)
            assert 0<=exact<Q(num,DEN), (tag,label,a,b,exact,Q(num,DEN))
            exactsum+=exact; maxdegree=max(maxdegree,len(polynomial)-1)
            branches+=cellbranches; leafbranches+=cellbranches
        assert end==endpoint
        addition=p['mass'] if method=='26mass' else (1-p['box'][4])**7 if method=='27mass' else Q(0)
        upper=Q(subtotal,DEN)+addition
        assert upper==leaf['printed_upper'] and upper<(Q(999,1000) if method=='29' else Q(1))
        methods[method]+=1
        if method not in maxima or upper>maxima[method]['value']: maxima[method]={'tag':tag,'value':upper}
        receipts.append({'tag':tag,'method':method,'cells':len(leaf['cells']),
          'branches':leafbranches,'cell_numerator_sum':subtotal,'criterion':upper,
          'exact_integral_sum':exactsum,'max_polynomial_degree':maxdegree,'ell':p['ell']})
        if number%25==0:
            print(json.dumps({'leaves_checked':number,'branches':branches,'seconds':round(time.monotonic()-started,2)}),flush=True)
    summary={'status':'PASS','failures':0,'native_proof_lines':len(rows),
      'splits':split_count,'nonempty_boxes':len(leaves),'empty_boxes':len(empty),
      'cells':totalcells,'polynomial_branch_checks':branches,'methods':methods,
      'maximum_criteria':maxima,'labels':labels,'ell_values':sorted(set(x['ell'] for x in receipts)),
      'initial_source_height':initial,'fivefold_source_heights':heights,
      'Q_crossing_signs':[value(difference,CL),value(difference,CR)],
      'certificate_sections':{'prose':[1,149],'appendix_header':[150,156],
        'root':[157,160],'splits':[si+1,ei],'empty':[ei+1,li],
        'leaves':[li+1,ci],'cells':[ci+1,len(rows)]},
      'seconds':round(time.monotonic()-started,2)}
    assert (len(leaves),len(empty),totalcells)==(812,105,10922)
    assert dict(methods)=={'29':573,'27mass':236,'26mass':3}
    summary.update(key=key,status='pass',margins={method:str((Q(999,1000) if method=='29' else 1)-row['value']) for method,row in maxima.items()})
    return encode(summary)
