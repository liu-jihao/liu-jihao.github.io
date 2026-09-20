"""Scalar and rectangular smooth-source certificates; exact independent readers."""
from .common import Q, require, record_lines, unique_insert, split_box
from . import smooth_ops as s, bc74_ops as b, threefold_ops as c
import re
KEYS=('source-noncontained26','source-regular26','source-noncontained-initial','source-initial-fourfold','source-initial-threefold')
CONFIG={
 'source-noncontained26':('55cf',(Q(7,30),Q(2,5)),Q(29,25),Q(99,100)),
 'source-regular26':('e193',(Q(36,125),Q(2,5)),Q(29,25),Q(19,20)),
 'source-initial-fourfold':('d936',(Q(263,1800),Q(56167,60000),Q(7,30),Q(3,10)),Q(27,25),Q(999,1000)),
 'source-noncontained-initial':('bc74',b.ROOT,Q(29,25),Q(999,1000)),
}

def read_direct(text,root):
    families={};family=None;mode=None;tag=None
    for ln in record_lines(text):
        if ln.startswith('FAMILY '):
            family=ln.split()[1];unique_insert(families,family,dict(boxes={'R':root},leaves={},traces={},splits=0));mode=None;continue
        require(family is not None,'record before family')
        f=families[family]
        if ln.startswith('Root '):
            # The independently prescribed root is checked against the literal header.
            vectors=re.findall(r'\[([^\]]+)\]',ln)
            require(tuple(Q(x)for v in vectors for x in v.split(','))==root,'root header differs from formula')
            mode='splits'if 'Splits:'in ln else None;continue
        if ln.startswith('Splits:'):mode='splits';continue
        if ln.startswith('Empty boxes:'):mode='empty';continue
        if ln.startswith('Leaves:'):mode='leaves';continue
        if ln.startswith('Leaf '):
            tag=ln.split(':',1)[0][5:];require(tag in f['leaves'],'unknown trace leaf')
            unique_insert(f['traces'],tag,dict(cells=[],sum=None));mode='cells';continue
        if re.match(r'R[01]*\s*\|',ln):
            fields=[x.strip()for x in ln.split('|')];name=fields[0]
            if mode=='splits':
                axis=0 if len(root)==2 or fields[1]=='a'else 1
                if len(root)>2:require(fields[1]in ['a','lambda'],'split coordinate')
                require(name in f['boxes'],'split parent missing')
                lo,hi=split_box(f['boxes'].pop(name),axis,Q(fields[-1]))
                unique_insert(f['boxes'],name+'0',lo);unique_insert(f['boxes'],name+'1',hi);f['splits']+=1
            elif mode=='leaves':
                require(name in f['boxes'],'leaf outside cover')
                box=tuple(map(Q,fields[1:1+len(root)]));require(box==f['boxes'][name],'leaf parameters differ')
                unique_insert(f['leaves'],name,dict(box=box,numbers=list(map(int,fields[1+len(root):]))))
            else:raise ValueError('unexpected box record '+ln)
            continue
        if mode=='cells'and (ln.startswith('Sum ')or ln.startswith('Incoming sum ')):
            require(f['traces'][tag]['sum']is None,'duplicate sum')
            f['traces'][tag]['sum']=int(re.search(r'\d+',ln)[0]);continue
        if mode=='cells'and'|'in ln:
            fields=[x.strip()for x in ln.split('|')];require(len(fields)==4,'malformed rational cell')
            f['traces'][tag]['cells'].append((Q(fields[0]),Q(fields[1]),fields[2],int(fields[3])));continue
        raise ValueError('unconsumed direct record '+ln)
    return families

def direct(key,text):
    kind,root,end,target=CONFIG[key];families=read_direct(text,root)
    require(set(families)==({'Q'}if kind=='e193'else {'Q','M'}),'missing or extra direct family')
    results=[]
    for family,f in families.items():
        require(set(f['boxes'])==set(f['leaves'])==set(f['traces']),'incomplete direct coverage')
        cells=subcells=0;maximum=Q(0);min_slack=None
        for tag,leaf in f['leaves'].items():
            tr=f['traces'][tag];box=leaf['box'];cursor=Q(0);total=0
            for left,right,label,num in tr['cells']:
                require(left==cursor and left<right<=end,'cell coverage')
                knots=b.breaks(family,box,left,right)if kind=='bc74'else s.all_breaks(kind,family,box,left,right)
                polys=[b.polynomial(family,box,label,a,z)if kind=='bc74'else s.branch_poly(kind,family,box,label,a,z)for a,z in zip(knots,knots[1:])]
                require(all(poly==polys[0]for poly in polys),'merged cell changes polynomial')
                exact=s.integral(polys[0],left,right);slack=Q(num,10**9)-exact
                require(exact>=0 and slack>0,'failed rational cell '+key+'/'+family+'/'+tag+'/'+label)
                min_slack=min(min_slack,slack)if min_slack is not None else slack
                cursor=right;total+=num;cells+=1;subcells+=len(polys)
            require(cursor==end and total==tr['sum']==leaf['numbers'][0],'cell sum or final endpoint')
            value=Q(total,10**9)+((1-box[2])**7 if kind=='d936'else 0)
            require(value<target,'failed source criterion')
            if kind=='d936':require(value<Q(leaf['numbers'][-1],10**9)<target,'whole criterion upper')
            maximum=max(maximum,value)
        results.append(dict(group=family,leaves=len(f['leaves']),splits=f['splits'],cells=cells,branch_pieces=subcells,maximum_upper=str(maximum),target=str(target),strict_margin=str(target-maximum),minimum_cell_slack=str(min_slack)))
    return dict(family=key,groups=results,endpoint=str(end),coverage='complete closed prescribed roots, every split child and cell interval')

def threefold(text):
    families=c.parse(text);require(set(families)=={'Q','M'},'threefold family coverage');results=[]
    for family,f in families.items():
        c.verify_tree(f,family);maximum=Q(0);cells=pieces=0;min_slack=None;methods={}
        for tag,leaf in f['leaves'].items():
            require(leaf['method']in ['29','27mass'],'unknown endpoint method')
            methods[leaf['method']]=methods.get(leaf['method'],0)+1
            cursor=Q(0);num=0
            for left,right,label,N,line in f['cells'][tag]:
                require(left==cursor and left<right,'threefold time coverage')
                exact,parts=c.integrate_cell(family,leaf['box'],label,left,right)
                delta=Q(N,10**9)-exact;require(exact>=0 and delta>0,'threefold cell witness '+family+'/'+tag+'/'+label)
                min_slack=min(min_slack,delta)if min_slack is not None else delta
                cursor=right;num+=N;cells+=1;pieces+=len(parts)
            require(cursor==(c.T if leaf['method']=='29'else c.v),'threefold horizon')
            require(num==leaf['numerator']==f['sums'][tag],'threefold sum')
            full=Q(num,10**9)+(0 if leaf['method']=='29'else(1-leaf['box'][4])**7)
            require(full==leaf['total']<Q(999,1000),'threefold full criterion');maximum=max(maximum,full)
        require((len(f['leaves']),len(f['empty']),cells)=={'Q':(281,43,3065),'M':(213,36,2467)}[family],'threefold inventory')
        results.append(dict(group=family,leaves=len(f['leaves']),empty=len(f['empty']),cells=cells,branch_pieces=pieces,methods=methods,maximum_upper=str(maximum),strict_margin=str(Q(999,1000)-maximum),minimum_cell_slack=str(min_slack)))
    return dict(family='source-initial-threefold',groups=results,coverage='complete four-coordinate root, strict empty tests and all closed faces')

def validate(key,text):return threefold(text)if key=='source-initial-threefold'else direct(key,text)
