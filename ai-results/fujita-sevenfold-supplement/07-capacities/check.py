#!/usr/bin/env python3
"""Independent exact checker for the finite Section 07 arithmetic interface.

Python >=3.10, standard library only. No generator, optimizer, network, private
state, article path, or assertion statement is used. See README for scope.
"""
import argparse
import json
import re
import sys
from collections import Counter
from fractions import Fraction as Q
from math import comb
from pathlib import Path


class Invalid(ValueError):
    pass


def require(ok, message):
    if not ok:
        raise Invalid(message)


def keys(obj, fields):
    require(type(obj) is dict and set(obj)==set(fields), 'unexpected or missing fields: '+repr(obj))


def rational(x):
    require(type(x) is str and re.fullmatch(r'-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?',x) is not None,
            'malformed rational: '+repr(x))
    return Q(x)


def integer(x):
    require(type(x) is int, 'integer required: '+repr(x))
    return x


def vector(x, count, convert):
    require(type(x) is list and len(x)==count, 'wrong vector length')
    return [convert(y) for y in x]


def unique_object(pairs):
    out={}
    for key,value in pairs:
        require(key not in out,'duplicate JSON key: '+key)
        out[key]=value
    return out


def read_json(path):
    def reject(value):
        raise Invalid('non-exact JSON number: '+value)
    return json.loads(path.read_text(encoding='utf-8'),object_pairs_hook=unique_object,
                      parse_float=reject,parse_constant=reject)


# Ascending polynomial coefficients. Multiplication, differentiation and
# integration here are independent of the binomial witness producer.
def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0: a.pop()
    return a


def plus(*polys):
    result=[Q(0)]*max(map(len,polys))
    for p in polys:
        for j,x in enumerate(p): result[j]+=x
    return trim(result)


def scale(p,c): return trim([c*x for x in p])


def times(a,b):
    result=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): result[i+j]+=x*y
    return trim(result)


def power(p,n):
    result=[Q(1)]
    for _ in range(n): result=times(result,p)
    return result


def derivative(p): return trim([j*x for j,x in enumerate(p)][1:] or [Q(0)])


def evaluate(p,x):
    result=Q(0)
    for c in reversed(p): result=result*x+c
    return result


def compose(p,q):
    result=[Q(0)]
    for c in reversed(p): result=plus(times(result,q),[c])
    return result


def positive(x): return max(Q(0),x)


def span_coefficients(t,s):
    return [(-1)**(j+1)*(comb(t+s,j)-s*comb(t+s-1,j-1)) for j in range(1,t+s+1)]


def span(t,s,c):
    return sum((a*positive(1-j*c)**6 for j,a in enumerate(span_coefficients(t,s),1)),Q(0))


def pair(n,c):
    return [positive(1-c)**(n-1),positive(1-2*c)**(n-2)*(1+2*(n-2)*c)]


def F6(s,p): return 1-Q(10)/p+(5*s+5)/p**2-s*s/p**3
def F7(s,p): return 1-Q(20)/p+(15*s+37)/p**2-(6*s*s+12*s)/p**3+s**3/p**4


SPANS={'linear-fivefold':(1,2),'linear-fourfold':(1,3),'full-fourfold':(0,4),'full-fivefold':(0,3)}
PAIRS={'six-fifth':(6,'1/5'),'seven-sixth':(7,'1/6'),'early-pair':(7,'3/20'),'later-pair':(7,'7/40')}
CONSTANTS=('six-interior','six-boundary','seven-interior','seven-boundary','four-row-constant',
           'strip-derivative','disjoint-correction','disjoint-square','three-cut-anchor',
           'envelope-rounding','linear-endpoint','linear-later','full-later')
EARLY={
    (3,7,5):('7/5','171/100','211/1000'),
    (4,5,2):('3','119/100','287/500'),
    (4,6,3):('8/3','1317/1000','71/125'),
    (4,6,4):('2','283/200','379/1000'),
    (4,7,4):('5/2','283/200','29/50'),
    (4,7,5):('2','187/125','17/40'),
    (4,7,6):('5/3','783/500','331/1000'),
    (4,7,7):('10/7','1627/1000','11/40')}
LATER={k:v for k,v in EARLY.items() if k in ((4,6,3),(4,6,4),(4,7,4),(4,7,5),(4,7,6))}
OWN={(5,6,2):('4','1149/1000','247/250'),(5,7,3):('11/3','623/500','39/40'),
     (5,7,4):('3','33/25','77/100'),(5,7,5):('13/5','69/50','673/1000')}
DIVISORS={1:('1','198/125'),2:('1123/1000','703/500'),3:('1201/1000','223/200'),
          4:('63/50','121/125'),5:('327/250','801/1000'),6:('5/3','42/125')}
EASY={29:('3/2','1/4'),30:('3/2','1/4'),31:('3/2','1/4'),32:('5/3','229/1000'),
      33:('3/2','59/200'),34:('3/2','59/200'),35:('783/500','269/1000'),
      36:('1627/1000','61/250'),39:('5/3','133/500'),40:('5/3','131/500'),
      41:('5/3','131/500'),42:('1627/1000','11/40'),45:('5/3','42/125')}
HARD={28:('4','7/4','0'),37:('1/2','841/500','111/500'),38:('2/9','1733/1000','14/125'),
      43:('5/4','841/500','32/125'),44:('10/9','1733/1000','239/1000'),
      46:('8/11','7/4','131/500'),47:('1/2','931/500','49/250'),
      48:('4/13','1899/1000','17/100'),49:('1/7','7/4','1/10')}


def expected_ids():
    result={key:'spanning' for key in SPANS}
    result.update({key:'pairs' for key in PAIRS})
    result.update({key:'constants' for key in CONSTANTS})
    result.update({key:'polynomial' for key in ('quartic','quartic-derivative','quartic-endpoints')})
    for prefix,states in (('early',EARLY),('later',LATER),('own',OWN)):
        for d,e,m in states: result[f'{prefix}-{d}{e}{m}']='source-tables'
    result.update({f'fourfold-{j}':'fourfold-majorant' for j in range(28,50)})
    result.update({f'divisor-{j}':'baseline' for j in range(1,7)})
    result.update({f'fivefold-{j}':'baseline' for j in range(3,6)})
    result.update({f'primitive-{j}':'integrals' for j in range(5)})
    result.update({key:'crossings' for key in ('crossing-left','crossing-right')})
    result.update({f'join-{j}':'envelopes' for j in range(1,7)})
    result.update({f'five-normal-{s}':'envelopes' for s in ('sixth','third','half')})
    result.update({f'K5-{s}':'envelopes' for s in ('zero','switch','half')})
    return result


def validate(data):
    require(data.is_dir(),'missing data directory')
    require({p.name for p in data.iterdir()}=={'certificate.json','bindings.json'},'unexpected/missing data files')
    doc=read_json(data/'certificate.json')
    keys(doc,('schema','records'))
    require(doc['schema']=='fujita-capacities-v1','unknown schema')
    require(type(doc['records']) is list,'records must be a list')
    expected=expected_ids()
    records={}
    for record in doc['records']:
        require(type(record) is dict and type(record.get('id')) is str,'invalid record')
        key=record['id']
        require(key not in records,'duplicate record: '+key)
        require(key in expected and record.get('family')==expected[key],'unknown record/family: '+key)
        records[key]=record
    require(set(records)==set(expected),'missing coverage: '+','.join(sorted(set(expected)-set(records))))
    bd=read_json(data/'bindings.json')
    keys(bd,('schema','article_sha256','bindings'))
    require(bd['schema']=='fujita-capacities-bindings-v1','unknown binding schema')
    article_names=('section.tex','capacity-numerics.tex','capacity-refinements.tex','source-layers.tex','finite-certificate.tex')
    keys(bd['article_sha256'],article_names)
    for value in bd['article_sha256'].values():
        require(type(value) is str and re.fullmatch('[0-9a-f]{64}',value) is not None,'invalid article identity')
    require(type(bd['bindings']) is list,'bindings must be a list')
    seen=set()
    for b in bd['bindings']:
        keys(b,('id','file','formula','snippet'))
        require(all(type(v) is str and v for v in b.values()),'invalid binding field')
        require(b['id'] in expected and b['id'] not in seen,'unknown/duplicate binding')
        require(b['file'] in article_names,'unrecognized article binding')
        seen.add(b['id'])
    require(seen==set(expected),'missing article bindings')

    receipts=[]
    def emit(key, **values):
        receipts.append(dict(id=key,family=expected[key],values={k:str(v) for k,v in values.items()}))
    def equal(actual,witness,key):
        require(actual==witness,f'{key}: altered/incorrect witness; actual {actual}, witness {witness}')
    def fields(record,*extra): keys(record,('id','family',*extra))

    for key,(t,s) in SPANS.items():
        r=records[key]; fields(r,'n','t','s','coefficients')
        equal((integer(r['n']),integer(r['t']),integer(r['s'])),(7,t,s),key+' domain')
        cs=vector(r['coefficients'],t+s,integer)
        equal(cs,span_coefficients(t,s),key)
        emit(key,coefficient_count=len(cs),coefficients=cs)
    for key,(n,c0) in PAIRS.items():
        r=records[key]; fields(r,'n','c','values')
        c=rational(r['c']); equal((integer(r['n']),c),(n,Q(c0)),key+' domain')
        values=vector(r['values'],2,rational); equal(values,pair(n,c),key)
        bound=Q(11,25) if key=='early-pair' else Q(17,50) if key=='later-pair' else max(values)
        if key in ('early-pair','later-pair'): require(max(values)<bound,key+' strict bound')
        emit(key,coordinate=values[0],balanced=values[1],target=bound,gap=bound-max(values))

    C4=Q(13,96)+Q(2,3)**5*Q(432,3125)
    evaluated={
        'six-interior':F6(Q(25,3),Q(125,27)), 'six-boundary':F6(Q(8),Q(4)),
        'seven-interior':F7(Q(12),Q(8)), 'seven-boundary':F7(Q(45,4),Q(25,4)),
        'four-row-constant':C4,'strip-derivative':4*Q(2,5)**2*Q(3,5)**3,
        'disjoint-correction':Q(1,16)/Q(1,4)*Q(1,3)**5,
        'disjoint-square':Q(31,243)+Q(1,972), 'three-cut-anchor':Q(5,6)**5*max(pair(6,Q(1,5))),
        'envelope-rounding':Q(31,200)+Q(169,200)*Q(67,200),
        'linear-endpoint':span(1,2,Q(2,11)), 'linear-later':span(1,2,Q(19,100)),
        'full-later':span(0,4,Q(9,50))}
    for key in CONSTANTS:
        r=records[key]; fields(r,'value'); equal(evaluated[key],rational(r['value']),key)
        emit(key,value=evaluated[key])
    require(C4<Q(31,200) and Q(5,6)**6<Q(67,200),'upward rounding direction')
    require(C4+(1-C4)*Q(5,6)**6<evaluated['envelope-rounding']<Q(11,25),'rounded full envelope')
    require(evaluated['linear-endpoint']<Q(2,5) and evaluated['linear-later']<Q(17,50)
            and evaluated['full-later']<Q(17,50),'full competing spanning rates')
    require(Q(5,16)<Q(992,3125)<Q(1053,3125)<Q(11,32)<Q(256,729),'three-weight endpoint comparisons')
    require(Q(125,972)<Q(13,96),'disjoint square strict margin')
    require(Q(256,729)+Q(9,500)<Q(3,8) and Q(256,729)+Q(3,50)<Q(21,50),'level strip margins')

    r=records['quartic']; fields(r,'coefficients'); R=vector(r['coefficients'],5,integer)
    s2=[0,12,-3]; s3=[0,0,6,-2]
    numerator=plus(power(s3,4),scale(power(s3,3),-20),
        times(plus(scale(s2,15),[37]),power(s3,2)),
        scale(times(plus(scale(power(s2,2),6),scale(s2,12)),s3),-1),power(s2,3))
    differentiated=plus(times(derivative(numerator),s3),scale(times(numerator,derivative(s3)),-4))
    factored=scale(times(times(times(times([0,0,0,0,1],[2,-1]),[-1,1]),[5,-2]),compose(R,[3,-1])),12)
    equal(differentiated,factored,'quartic derivative factorization')
    equal(derivative(derivative(R)),scale(times([-3,4],[-7,4]),30),'quartic second derivative roots')
    r=records['quartic-derivative']; fields(r,'values')
    points=[Q(1,2),Q(3,4),Q(7,4),Q(2)]
    dvalues=[evaluate(derivative(R),x) for x in points]
    equal(dvalues,vector(r['values'],4,rational),'quartic-derivative')
    require(dvalues[0]>0 and dvalues[1]>0 and dvalues[2]<0 and dvalues[3]<0,'quartic monotonicity signs')
    r=records['quartic-endpoints']; fields(r,'values')
    values=[evaluate(R,Q(1,2)),evaluate(R,Q(2))]
    equal(values,vector(r['values'],2,rational),'quartic-endpoints'); require(min(values)>0,'quartic positivity')
    emit('quartic',identity='coefficient equality',degree=len(R)-1)
    emit('quartic-derivative',values=dvalues)
    emit('quartic-endpoints',values=values,minimum=min(values))

    for prefix,states,u0 in (('early',EARLY,'51/50'),('later',LATER,'26/25'),('own',OWN,'29/25')):
        for (d,e,m),inputs in states.items():
            key=f'{prefix}-{d}{e}{m}'; r=records[key]
            fields(r,'d','e','m','u','cap','a','p','beta','sigma','tail','hilbert','e1')
            equal(tuple(integer(r[x]) for x in ('d','e','m')),(d,e,m),key+' state')
            u=rational(r['u']); equal(u,Q(u0),key+' order')
            cap,a,p=[rational(r[x]) for x in ('cap','a','p')]
            equal((cap,a,p),tuple(map(Q,inputs)),key+' source parameters')
            require(a*u>1 and 7*a+p>8 and a*cap+p<8,key+' denominator/monotonicity')
            beta=((8-p)*u-7)/(a*u-1)
            sigma=7-cap+u*(a*cap+p-8)
            equal(sigma,rational(r['sigma']),key+' slack')
            if prefix!='later': equal(beta,rational(r['beta']),key+' deficit')
            else: require(r['beta'] is None,key+' unexpected beta witness')
            require(beta>1 and sigma>=0,key+' tangent applicability numbers')
            gate=Q(1) if key=='early-474' else Q(1,2)
            require(sigma<gate,key+' strict tangent gate')
            if prefix=='own':
                equal(a*cap+p,rational(r['tail']),key+' complete cap tail')
                require(type(r['hilbert']) is str,key+' Hilbert polynomial')
                terms=r['hilbert'].split('+'); coeffs={}
                for term in terms:
                    match=re.fullmatch(r'(?:(\d*)z(?:\^(\d+))?|(\d+))',term)
                    require(match is not None,key+' malformed Hilbert polynomial')
                    degree=int(match[2] or '1') if match[3] is None else 0
                    coefficient=int(match[1] or '1') if match[3] is None else int(match[3])
                    require(degree not in coeffs and degree<=2 and coefficient>0,key+' duplicate/out-of-range Hilbert term')
                    coeffs[degree]=coefficient
                require(coeffs.get(0)==1 and coeffs.get(1)==e-d and sum(coeffs.values())==m,key+' entire Hilbert numerator')
                equal(sum(j*v for j,v in coeffs.items()),rational(r['e1']),key+' first moment')
            else:
                require(r['tail'] is None and r['hilbert'] is None and r['e1'] is None,key+' unexpected moment data')
            if key=='early-474':
                split=7-Q(9,4)+u*(a*Q(9,4)+p-8)
                equal(split,Q(17161,40000),key+' reducedness split')
                require(split<Q(1,2),key+' split margin')
            emit(key,beta=beta,sigma=sigma,gate=gate,gap=gate-sigma)

    for j in range(28,50):
        key=f'fourfold-{j}'; r=records[key]; fields(r,'row','a','p','cap','value')
        equal(integer(r['row']),j,key+' index')
        a,p=rational(r['a']),rational(r['p'])
        if j in EASY:
            equal((a,p),tuple(map(Q,EASY[j])),key+' parameters')
            require(r['cap'] is None and r['value'] is None,key+' unexpected cap/value')
            require(a<=Q(5,3) and p<=Q(42,125),key+' coefficient domination')
            bound=p
        else:
            cap=rational(r['cap'])
            equal((cap,a,p),tuple(map(Q,HARD[j])),key+' parameters')
            require(cap>0,key+' interval domain')
            bound=positive(a-Q(5,3))*cap+p
            equal(bound,rational(r['value']),key+' endpoint value')
            require(bound<Q(42,125),key+' strict intercept bound')
        emit(key,intercept_bound=bound,gap=Q(42,125)-bound)

    for m,(a0,p0) in DIVISORS.items():
        key=f'divisor-{m}'; r=records[key]; fields(r,'m','a','p')
        equal(integer(r['m']),m,key+' index')
        a,p=rational(r['a']),rational(r['p']); equal((a,p),(Q(a0),Q(p0)),key+' full-tail parameters')
        den=7*a+p-8; require(den>0,key+' denominator')
        slope=m*a/den-Q(125,73); at_one=m*(a-1)/den
        require(slope>=0 and at_one>=0,key+' unbounded affine comparison')
        emit(key,denominator=den,slope_difference=slope,value_at_one=at_one)
    for m in range(3,6):
        key=f'fivefold-{m}'; r=records[key]; fields(r,'m','affine','comparison')
        equal(integer(r['m']),m,key+' index')
        A,B,C=vector(r['affine'],3,integer)
        # Derive 1-2/(7-beta(u)) from that SAME row's complete price.
        _,a,p=map(Q,OWN[(5,7,m)])
        require(C>0 and 7*a+p-8>0,key+' denominator')
        ratio=Q(C)/(7*a+p-8)
        equal(Q(A),2*ratio,key+' constant from beta')
        equal(Q(B),(8-5*a-p)*ratio,key+' coefficient from beta')
        coefficients=[8*C*125-9*73*A,8*C*52-9*73*B]
        equal(coefficients,vector(r['comparison'],2,integer),key+' cleared affine numerator')
        margin=Q(coefficients[0])-coefficients[1]*Q(34,25)
        require(A-B*Q(34,25)>0 and coefficients[1]>0 and margin>0,key+' whole closed interval')
        emit(key,endpoint_numerator=margin,interval='[29/25,34/25]')
    require(Q(9,8)**6>2,'factor-two comparison')
    # The separate double row includes its equality endpoint.
    require(Q(55875)==44998*Q(375,302),'double-fivefold endpoint identity')
    require((1149*Q(29,25)-1000)/(1031*Q(29,25))>=Q(4,21),'double-fivefold pair endpoint domain')

    primitives={}
    for j in range(5):
        key=f'primitive-{j}'; r=records[key]; fields(r,'j','coefficients')
        equal(integer(r['j']),j,key+' knot')
        coefficients=vector(r['coefficients'],8,integer)
        density=power([-j,1],6)
        equal(derivative(coefficients),scale(density,7),key+' derivative')
        equal(evaluate(coefficients,Q(j)),Q(0),key+' matching knot')
        primitives[j]=coefficients
        emit(key,coefficient_count=8,zero_at_knot=True,active_interval=f'[{j},infinity)')
    # Shared COMPLETE cell partition [0,1],...,[3,4],[4,infinity).
    # Check coefficient identities on each cell, including the unbounded cell;
    # do not replace a polynomial identity by an endpoint/sample calculation.
    # The article's union-volume proof supplies nonnegativity and monotonicity.
    for t,s in SPANS.values():
        for lower in range(5):
            density=[Q(0)]; primitive=[Q(0)]
            for j,a in enumerate(span_coefficients(t,s),1):
                if j<=lower:
                    density=plus(density,scale(power([-j,1],6),a))
                    primitive=plus(primitive,scale(primitives[j],a))
            equal(derivative(primitive),scale(density,7),'full signed polynomial on knot cell')

    for key,c0 in (('crossing-left','16020581569/90194313216'),('crossing-right','8010290789/45097156608')):
        r=records[key]; fields(r,'c','difference'); c=rational(r['c']); equal(c,Q(c0),key+' bracket')
        values=pair(7,c); delta=values[1]-values[0]
        equal(delta,rational(r['difference']),key+' sign witness')
        require(delta>0 if key=='crossing-left' else delta<0,key+' crossing sign')
        emit(key,c=c,balanced_minus_coordinate=delta)
    require(0<Q('16020581569/90194313216')<Q('8010290789/45097156608')<Q(1,2),'crossing bracket order')
    for j in range(1,7):
        key=f'join-{j}'; r=records[key]; fields(r,'r','value'); equal(integer(r['r']),j,key+' height')
        c=Q(1,7); low=1-(Q(7*j,j+1)*c)**j
        high=1-((6*j+7*c-1)/(6*(j+1)))**j
        equal(low,high,key+' both branches'); equal(low,rational(r['value']),key+' value')
        emit(key,join=low)
    for tag,c0 in (('sixth','1/6'),('third','1/3'),('half','1/2')):
        key='five-normal-'+tag; r=records[key]; fields(r,'c','values')
        c=rational(r['c']); equal(c,Q(c0),key+' argument')
        require(36-42*c>0 and 1-c>0,key+' all denominators')
        coord=(1-c)**6
        chord=1-(Q(5,6)+(c-Q(1,7))/(6*(1-c-Q(1,7))))**5
        v=(8*c-1)/6
        # At c=1/2, use the limit of the call bound; do not divide by 0.
        call=Q(0) if c==Q(1,2) else (1-c-v)**7/(7*(1-2*c)*(1-c)**5*(c-v))
        b=1-((30+7*c-1)/36)**5
        moment=(Q(7,6)*(1-c))**6
        e=coord if c==Q(1,2) else min(b,moment)
        values=[coord,chord,call,b,moment,e,max(coord,chord),max(coord,call)]
        equal(values,vector(r['values'],8,rational),key+' all competing expressions')
        emit(key,values=values,full_minimum=min(values[5:]))
    require(7**6<2**17,'five-normal coordinate crossing')
    for tag,c0 in (('zero','0'),('switch','4/25'),('half','1/2')):
        key='K5-'+tag; r=records[key]; fields(r,'c','values','single_enabled')
        c=rational(r['c']); equal(c,Q(c0),key+' argument')
        one=positive(1-2*c)**6
        values=[span(0,3,c),2*one,one]
        equal(values,vector(r['values'],3,rational),key+' competitors')
        require(type(r['single_enabled']) is bool and r['single_enabled']==(c>=Q(4,25)),key+' single-pair range')
        emit(key,values=values,enabled='all' if r['single_enabled'] else 'signed and factor-two')

    require(len(receipts)==len(expected),'internal coverage error')
    margins={'early_rounded_envelope':Q(11,25)-evaluated['envelope-rounding'],
             'later_full_fourfold':Q(17,50)-evaluated['full-later'],
             'later_linear_fourfold':Q(17,50)-evaluated['linear-later'],
             'retained_three_eighths':Q(3,8)-Q(256,729)-Q(9,500),
             'early_three_normal':Q(21,50)-Q(256,729)-Q(3,50),
             'five_normal_endpoint_integer_gap':2**17-7**6}
    return dict(status='pass',schema='fujita-capacities-report-v1',record_count=len(receipts),
                family_counts=dict(sorted(Counter(expected.values()).items())),
                integral_partition=['[0,1]','[1,2]','[2,3]','[3,4]','[4,infinity)'],
                signed_integral_polynomial_cells=20,
                arithmetic='integer and Fraction only; no numerical error',
                strict_margins={k:str(v) for k,v in margins.items()},records=receipts,
                scope='Finite article identities, values, coefficient witnesses and affine endpoint reductions. '
                      'Continuum maxima, geometric applicability, full scalar census and asymptotic errors remain article proofs.')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    args=p.parse_args(argv)
    try:
        result=validate(args.data)
        code=0
    except (ValueError,TypeError,KeyError,IndexError,OSError,ZeroDivisionError) as exc:
        result=dict(status='fail',error=str(exc),error_type=type(exc).__name__)
        code=1
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
    return code


if __name__=='__main__':
    sys.exit(main())
