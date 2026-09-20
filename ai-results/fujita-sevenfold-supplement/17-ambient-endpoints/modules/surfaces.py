"""Independent checks of the two initial surface histories and distinct primes."""
from .common import Q,require,record_lines,unique_insert
from . import surface_ops as c, distinct_ops as d
import re

def surface(key,text):
    constants=c.constants();families=c.parse(text);require(set(families)=={'Q','M'},'surface family coverage');results=[]
    for family,f in families.items():
        c.verify_tree(f,family);maxima={};cells=pieces=0;minimum=None;methods={}
        for tag,leaf in f['leaves'].items():
            method=leaf['method'];methods[method]=methods.get(method,0)+1;cursor=Q(0);num=0
            for left,right,label,N,line in f['cells'][tag]:
                require(left==cursor and left<right,'surface time coverage')
                exact,parts=c.integrate_cell(family,leaf['box'],label,left,right)
                require(len({p[2]for p in parts})==1,'nonpolynomial merged surface cell')
                slack=Q(N,10**9)-exact;require(exact>=0 and slack>0,('surface cell',family,tag,label))
                minimum=slack if minimum is None else min(minimum,slack)
                cursor=right;num+=N;cells+=1;pieces+=len(parts)
            require(cursor==(c.T if method=='29'else c.v),'surface horizon')
            require(num==leaf['numerator']==f['sums'][tag],'surface sum')
            full=Q(num,10**9)+(0 if method=='29'else(1-leaf['box'][4])**7)
            target=Q(999,1000)if method=='29'else Q(1)
            require(full==leaf['total']<target,'surface full criterion')
            maxima[method]=max(maxima.get(method,Q(0)),full)
        require((len(f['leaves']),len(f['empty']),cells)=={'Q':(1075,278,12408),'M':(611,112,6722)}[family],'surface inventory')
        require(max(maxima.values())=={'Q':Q(8589931765637723389,8589934592000000000),'M':Q(29323928302434726949,29353417113600000000)}[family],'surface maximum')
        results.append(dict(group=family,leaves=len(f['leaves']),empty=len(f['empty']),cells=cells,branch_pieces=pieces,methods=methods,maximum_by_method={m:str(x)for m,x in maxima.items()},margins={m:str((Q(999,1000)if m=='29'else 1)-x)for m,x in maxima.items()},minimum_cell_slack=str(minimum)))
    return dict(key=key,status='pass',groups=results,constants=constants)

def distinct(key,text):
    families={};family=None;tag=None
    for ln in record_lines(text):
        if ln.startswith('FAMILY '):
            family=ln.split(';')[0][7:];require(family in ['NQ','NE','IQ','IM'],'distinct family')
            unique_insert(families,family,{'cells':[]});continue
        require(family is not None,'distinct record before family');f=families[family]
        if ln=='Splits: node | midpoint.':continue
        if ln.startswith('Leaf '):
            match=re.fullmatch(r'Leaf R theta=0,([01]); method=(26|27);\s*upper=([\d/]+)<999/1000\.',ln)
            require(match is not None and 'upper'not in f,'distinct leaf header')
            theta,method,upper=match.groups();require((int(theta),int(method))==((0,27)if family[0]=='N'else(1,26)),'distinct root and horizon')
            f['upper']=Q(upper);continue
        if ln=='left | right | polynomial | numerator':continue
        if ln.startswith('Sum numerator='):
            require('sum'not in f,'duplicate distinct sum');f['sum']=int(re.fullmatch(r'Sum numerator=(\d+)\.',ln)[1]);continue
        fields=[v.strip()for v in ln.split('|')];require(len(fields)==4,'unexpected distinct record')
        f['cells'].append((Q(fields[0]),Q(fields[1]),fields[2],int(fields[3])))
    require(set(families)=={'NQ','NE','IQ','IM'},'distinct missing family');reports=[]
    cross=d.minus(d.times(d.power(d.minus(d.ONE,d.scale(d.T,2)),5),d.plus(d.ONE,d.scale(d.T,10))),d.power(d.minus(d.ONE,d.T),6))
    seq=d.sturm(d.trim(cross[1:]));require(d.variations(seq,Q(0))-d.variations(seq,Q(1,2))==1 and d.value(cross,d.cm)>0>d.value(cross,d.cp),'unique capacity crossing')
    for family,f in families.items():
        cursor=Q(0);num=0;slacks=[];branches=0
        for left,right,label,N in f['cells']:
            require(left==cursor and left<right,'distinct cell cover');value=Q(0)
            for a,b,F,H0,H1,D in d.branches(family,left,right):
                p=d.polynomial(family,label,a,b,F,H0,H1,D);d.nonnegative(p,a,b)
                value+=7*d.integral(p,a,b);branches+=1
            slack=Q(N,10**9)-value;require(0<=value and slack>0,'distinct integral')
            slacks.append(slack);cursor=right;num+=N
        require(cursor==(d.v if family[0]=='N'else d.u)and num==f['sum'],'distinct horizon and sum')
        mass=Q(23,30)**7 if family[0]=='N'else Q(15263,18750)**7
        criterion=Q(num,10**9)+mass;require(criterion==f['upper']<Q(999,1000),'distinct full mass criterion')
        reports.append(dict(group=family,cells=len(f['cells']),branches=branches,mass=str(mass),maximum_upper=str(criterion),strict_margin=str(Q(999,1000)-criterion),minimum_cell_slack=str(min(slacks))))
    require(sum(len(f['cells'])for f in families.values())==31,'distinct cell inventory')
    return dict(key=key,status='pass',groups=reports)
