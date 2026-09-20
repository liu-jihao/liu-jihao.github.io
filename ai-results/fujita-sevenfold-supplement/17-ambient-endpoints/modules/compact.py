"""Validate all delivered linked-profile rectangles and literal cell traces.
Polynomial arithmetic adapted from the independent f08 and 981 audits.
"""
from itertools import product
from .common import Q, require, record_lines, unique_insert, split_box
from . import f08_ops as f08, point_ops as point
import re

SPECS={
 'reserve-d26':[('G','131/1000','9/10','16'),('D','189/1000','26/25','8')],
 'reserve-q26-d27':[('G','131/1000','9/10','16'),('Q','37/200','26/25','8'),('D','17/80','27/25','12')],
 'reserve-m26-d27':[('G','131/1000','9/10','16'),('M','41/200','26/25','8'),('D','17/80','27/25','12')],
 'reserve-q26-c27':[('G','131/1000','9/10','16'),('Q','37/200','26/25','8'),('Q','1609/5000','27/25','12')],
 'reserve-m26-c27':[('G','131/1000','9/10','16'),('M','41/200','26/25','8'),('Q','1609/5000','27/25','12')],
}
for key,kind,height,anchor in [('d26','D','189/1000','26/25'),('d27','D','17/80','27/25'),('c27','Q','1609/5000','27/25')]:
    SPECS['point-'+key]=[('G','131/1000','9/10','2'),(kind,height,anchor,'8'),('T','257/1000','29/25','9/4')]
EXPECTED={'reserve-d26':(133,1310,9842943211),'reserve-q26-d27':(787,8590,9843582362),'reserve-m26-d27':(237,2694,9843027297),'reserve-q26-c27':(277,2489,9843540583),'reserve-m26-c27':(371,3821,9843229660),'point-d26':(76,826,9948845674),'point-d27':(279,3538,9989986795),'point-c27':(60,606,8870441956)}

def validate(key,text):
    qdict,pdict,grids,specs={},{},{},[]
    splits,rows=[],[];mode='preamble'
    for ln in record_lines(text):
        # Whitespace in dictionary lists and field separators is documentary.
        m=re.fullmatch(r'(q\d+)\s*=\s*(\d+(?:/\d+)?)',ln)
        if m: unique_insert(qdict,m[1],Q(m[2]));continue
        m=re.fullmatch(r'(p\d+)\s*=\s*(\S+)',ln)
        if m: unique_insert(pdict,m[1],m[2]);continue
        m=re.fullmatch(r'(\d+)\s*\|\s*([GQMDT])\s*\|\s*(\S+)\s*\|\s*(\S+)\s*\|\s*(\S+)',ln)
        if m:
            require(int(m[1])==len(specs),'profile specification index')
            specs.append((m[2],Q(m[3]),Q(m[4]),Q(m[5])));continue
        m=re.fullmatch(r'Grid(\d+)=(.*)\.',ln)
        if m: unique_insert(grids,int(m[1]),[qdict[s.strip()]for s in m[2].split(',')]);continue
        if ln.startswith('Splits:'):mode='splits';continue
        if ln.startswith('Final rows:'):mode='rows';continue
        if mode=='splits':
            m=re.fullmatch(r'([0-9.]+)/([0-9]+)/(q\d+)',ln)
            require(m is not None,'malformed split '+ln)
            splits.append((m[1],int(m[2]),qdict[m[3]]));continue
        if mode=='rows':
            fields=[s.strip()for s in ln.split('|')];require(len(fields)==4,'malformed final row')
            tag,ps,n,trace=fields
            values=[qdict[s.strip()]for s in ps.split(',')]
            cells=[(qdict[a.strip()],pdict[b.strip()])for a,b in (s.split('/')for s in trace.split(','))]
            rows.append((tag,values,int(n),cells));continue
        require(ln.startswith(('COMPACT FULL CERTIFICATE','Each q-symbol','Every row','A trace','Integrate each','Rational dictionary:','Polynomial dictionary:','FAMILY ','Endpoint=','Specifications:','Heights:','Before an anchor')),'unknown certificate record '+ln)
    expected=[(k,Q(h),Q(a),Q(c))for k,h,a,c in SPECS[key]]
    require(specs==expected,'source profile specifications differ from article')
    n=len(specs);require(set(grids)==set(range(n)),'missing grid')
    for i,g in grids.items():
        require(len(g)>=2 and all(a<b for a,b in zip(g,g[1:])),'grid ordering')
        require(g[0]==specs[i][1]/specs[i][2] and g[-1]==specs[i][3],'incomplete slope root')
    for dictionary,prefix in [(qdict,'q'),(pdict,'p')]:
        require(set(dictionary)=={prefix+str(i)for i in range(1,len(dictionary)+1)},'dictionary index coverage')
    leaves={}
    for inds in product(*(range(1,len(grids[i]))for i in range(n))):
        leaves['.'.join(map(str,inds))]=tuple(x for i,j in enumerate(inds)for x in (grids[i][j-1],grids[i][j]))
    initial=len(leaves)
    for tag,axis,mid in splits:
        require(tag in leaves,'missing or repeated split parent '+tag)
        lo,hi=split_box(leaves.pop(tag),axis,mid)
        unique_insert(leaves,tag+'0',lo);unique_insert(leaves,tag+'1',hi)
    require(len(rows)==len(leaves),'incomplete final leaf count')
    endpoint=Q('29/25')if key.startswith('reserve')else Q('13/5')
    target=Q(63,64)if key.startswith('reserve')else Q(999,1000)
    seen=set();cells_count=0;maximum=0;min_rounding=None
    for tag,values,total,cells in rows:
        require(tag not in seen and tag in leaves,'duplicate or unknown leaf '+tag);seen.add(tag)
        require(len(values)==3*n,'missing row parameters')
        slopes,heights=values[:2*n],values[2*n:]
        require(tuple(slopes)==leaves[tag],'leaf box mismatch '+tag)
        require(heights==f08.derive_heights(specs,slopes),'linked height mismatch '+tag)
        left=Q(0);exact=Q(0)
        for right,label in cells:
            require(left<right<=endpoint,'cell ordering or range '+tag)
            if key.startswith('reserve'):
                area=f08.integrate_cell(left,right,label,specs,slopes,heights)
            else:
                poly=point.polynom(label,left,right,slopes,heights,[s[2]for s in specs],[s[0]for s in specs])
                area=point.integral(poly,left,right)
            require(area>=0,'negative integral '+tag);exact+=area;left=right;cells_count+=1
        require(left==endpoint,'missing terminal interval '+tag)
        upper=Q(total,10**10);require(exact<upper<target,'failed strict integral or endpoint margin '+tag)
        min_rounding=min(min_rounding,upper-exact)if min_rounding is not None else upper-exact
        maximum=max(maximum,total)
    require(seen==set(leaves),'leaf coverage')
    require((len(rows),cells_count,maximum)==EXPECTED[key],'claimed inventory or maximum mismatch')
    return dict(family=key,rectangles=len(rows),cells=cells_count,initial_boxes=initial,splits=len(splits),target=str(target),maximum_upper=str(Q(maximum,10**10)),strict_margin=str(target-Q(maximum,10**10)),minimum_rounding_slack=str(min_rounding),endpoint=str(endpoint),coverage='complete closed slope grid and all bisection children',arithmetic='exact polynomial integration with full branch checks')
