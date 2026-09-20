"""Read the delivered trees and witnesses, independently reconstruct every box.
Arithmetic follows the independent all-cells audit, with no cached assertions.
"""
from .common import Q, require, record_lines, unique_insert, split_box
from . import persistent_ops as p
import re,json

GROUPS={'persistent-fourfold':['NE','M'],'persistent-small-surface':['Q4'],'persistent-fivefold':['NQ','QT','Q5','QN']}

def parse(text):
    lines=record_lines(text);trees={};witnesses={};i=0
    while i<len(lines):
        ln=lines[i];i+=1
        m=re.fullmatch(r'Family (\w+); domain version (\d+)\.',ln)
        if m:
            group,version=m[1],int(m[2]);head=lines[i];mx=lines[i+1];tree=lines[i+2];i+=3
            counts=re.fullmatch(r'Leaves (\d+), empty leaves (\d+), internal nodes (\d+), audited rational cells (\d+)\.',head)
            maximum=re.fullmatch(r'Largest upper integer (\d+) with denominator10000000\.',mx)
            require(counts is not None and maximum is not None,'tree headers')
            require(re.fullmatch(r'[\[\],0-9E]+',tree)is not None,'malformed tree')
            node=json.loads(tree.replace('E','"E"'))
            unique_insert(trees,group,dict(version=version,counts=tuple(map(int,counts.groups())),maximum=int(maximum[1]),tree=node));continue
        m=re.fullmatch(r'FAMILY (\w+); (\d+) nonempty leaves; (\d+) empty leaves; (\d+) cells\.',ln)
        require(m is not None,'unknown persistent record '+ln)
        group=m[1];qdict={};leaves={};empties={};current=None
        header=tuple(map(int,m.groups()[1:]));require(lines[i]=='Rational dictionary:','dictionary marker');i+=1
        while i<len(lines)and not lines[i].startswith('FAMILY '):
            ln=lines[i];i+=1
            m=re.fullmatch(r'(z\d+)=(-?\d+(?:/\d+)?)',ln)
            if m:unique_insert(qdict,m[1],Q(m[2]));continue
            if ln=='Written witnesses:':continue
            def vec(s):return tuple(qdict[x]for x in s.split(','))
            # Formula records are whitespace-free for structural parsing.
            compact=re.sub(r'\s+','',ln)
            m=re.fullmatch(r'Leaf(R[01]*)B=\[(.*)\]h=(z\d+)horizon=(26|27|29)mass=(z\d+)\.',compact)
            if m:
                current=m[1];unique_insert(leaves,current,dict(effective=vec(m[2]),height=qdict[m[3]],method=m[4],mass=qdict[m[5]]));continue
            if compact.startswith('Cells'):
                require(current is not None and 'cells'not in leaves[current],'duplicate or orphan cell trace')
                cells=[];left=Q(0)
                for item in compact[5:].rstrip('.').split(';'):
                    fields=item.split(':');require(len(fields)>=3,'malformed cell')
                    right=qdict[fields[0]];name=':'.join(fields[1:-1]);num=int(fields[-1])
                    require(left<right,'nonincreasing cell trace');cells.append((left,right,name,num));left=right
                leaves[current]['cells']=cells;continue
            m=re.fullmatch(r'Sum=(\d+);test(\d+)/1000000000\+\(1-(z\d+)\)\^7<=(\d+)/10000000\.',compact)
            if m:
                require(current is not None and 'upper_n'not in leaves[current],'duplicate or orphan total')
                leaf=leaves[current];summ=sum(c[3]for c in leaf['cells'])
                require(int(m[1])==int(m[2])==summ and qdict[m[3]]==leaf['mass'],'sum or mass mismatch')
                leaf['upper_n']=int(m[4]);continue
            m=re.fullmatch(r'Empty(R[01]*)raw=\[(.*?)\]B=\[(.*?)\]h=(z\d+);(.*)\.',compact)
            if m:
                claims=[]
                for claim in m[5].split(','):
                    name,relation=claim.split(':');a,b=relation.split('>');claims.append((name,qdict[a],qdict[b]))
                unique_insert(empties,m[1],dict(raw=vec(m[2]),effective=vec(m[3]),height=qdict[m[4]],claims=claims));continue
            raise ValueError('unconsumed persistent witness '+ln)
        require(set(qdict)=={'z'+str(j)for j in range(len(qdict))},'dictionary index coverage')
        unique_insert(witnesses,group,dict(header=header,leaves=leaves,empty=empties))
    return trees,witnesses

def validate(key,text):
    trees,groups=parse(text);require(set(trees)==set(groups)==set(GROUPS[key]),'missing or extra persistent group')
    results=[]
    for group in GROUPS[key]:
        meta=trees[group];p.version=meta['version'];expected_version=1 if group in ['NQ','NE']else 5
        require(p.version==expected_version,'wrong physical/normalized domain version')
        h0,hk,old,kind=p.data[group];lmax=(Q(157,250)+Q(2,5)+p.eps)/3
        root=(h0/p.r,(lmax+p.eps-h0)/(p.v-p.r),hk/p.u,25*(lmax+p.eps-hk),Q(7,30),lmax,Q(3003,10000),Q(157,250),Q(0),Q(0)if kind.startswith('N')else Q(1),Q(0),Q(1)if old=='S5'else Q(0))
        if p.version==5:root=(Q(0),Q(1))*4+root[8:]
        source=groups[group];leaves=source['leaves'];empties=source['empty'];seen=set();empty_seen=set();cells=nodes=0;maximum=0;min_cell=None;pending=[('R',root,meta['tree'])]
        while pending:
            tag,box,node=pending.pop()
            if isinstance(node,list):
                require(len(node)==3 and isinstance(node[0],int),'malformed internal tree node')
                axis=node[0]
                if p.version==1:
                    widths=[(box[2*j+1]-box[2*j])/(root[2*j+1]-root[2*j])if root[2*j+1]>root[2*j]else Q(0)for j in range(6)]
                    require(axis==max(range(6),key=lambda j:widths[j]),'split tie convention')
                require(0<=axis<6,'split axis out of range')
                a,b=split_box(box,axis,(box[2*axis]+box[2*axis+1])/2);nodes+=1
                pending.extend([(tag+'1',b,node[2]),(tag+'0',a,node[1])]);continue
            physical=p.enclose_normalized(group,box)if p.version==5 else box
            invalid,effective,h1=p.tighten(group,physical)
            if node=='E':
                require(tag in empties and tag not in empty_seen and invalid,'missing or false empty leaf '+tag);empty_seen.add(tag)
                record=empties[tag];require(record['raw']==box and record['effective']==effective and record['height']==h1,'empty witness parameters')
                comparisons={'a_low>a_high':(effective[0],effective[1]),'b_low>b_high':(effective[2],effective[3]),'lambda_low>lambda_high':(effective[4],effective[5]),'Lambda_low>Lambda_high':(effective[6],effective[7]),'h_1*v/u>lambda_high+epsilon':(h1*p.v/p.u,effective[5]+p.eps)}
                require(record['claims'],'missing empty witness')
                for name,a,b in record['claims']:require(name in comparisons and (a,b)==comparisons[name]and a>b,'false strict empty test')
                continue
            require(isinstance(node,int)and 0<=node<=9999000,'invalid terminal numerator')
            require(not invalid and tag in leaves and tag not in seen,'missing, repeated or empty nonempty leaf '+tag);seen.add(tag)
            leaf=leaves[tag];require(leaf['effective']==effective and leaf['height']==h1 and leaf.get('upper_n')==node,'effective leaf mismatch '+tag)
            end=Q(0);summ=0
            for left,right,name,num in leaf['cells']:
                require(left==end and left<right,'cell coverage '+tag)
                val=p.cell(group,effective,h1,name,left,right);slack=Q(num,10**9)-val
                require(val>=0 and slack>0,'invalid rational witness '+group+'/'+tag+'/'+name)
                min_cell=min(min_cell,slack)if min_cell is not None else slack;end=right;summ+=num;cells+=1
            ll,lh,ml,mh=effective[4:8]
            fu=max(Q(0),ll-(mh+p.eps-ll)*(p.v-p.u)/(p.T-p.v),ml-(1+1000*(p.T-ml)**7)*(p.T-p.u))
            mass={'26':fu,'27':ll,'29':ml}[leaf['method']]
            require(end=={'26':p.u,'27':p.v,'29':p.T}[leaf['method']]and mass==leaf['mass'],'terminal time or mass')
            upper=Q(summ,10**9)+(1-mass)**7
            require(upper<=Q(node,10**7)<Q(9999,10000),'whole leaf criterion')
            maximum=max(maximum,node)
        require(seen==set(leaves)and empty_seen==set(empties),'missing or extra leaf witnesses')
        require((len(seen),len(empty_seen),nodes,cells)==meta['counts'],'tree inventory')
        require((len(seen),len(empty_seen),cells)==source['header']and maximum==meta['maximum'],'witness/header inventory')
        results.append(dict(group=group,domain_version=p.version,leaves=len(seen),empty=len(empty_seen),internal_nodes=nodes,cells=cells,maximum_upper=str(Q(maximum,10**7)),strict_margin=str(Q(9999,10000)-Q(maximum,10**7)),minimum_cell_slack=str(min_cell)))
    return dict(family=key,groups=results,coverage='entire six-coordinate root, including fixed and collapsed coordinates; all empty leaves strictly certified',arithmetic='exact integration after every affine branch switch')
