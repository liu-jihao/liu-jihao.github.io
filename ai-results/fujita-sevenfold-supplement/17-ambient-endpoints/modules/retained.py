"""Read the complete early 1099-cell Bernstein/shell certificate."""
from .common import Q,require,record_lines,unique_insert
from . import retained_ops as c
import re,json

def validate(key,text):
    rows=[];seen=set();active=False
    for lineno,line in enumerate(record_lines(text),1):
        if line.startswith('Case '):
            m=re.fullmatch(r'Case (prime3|no_prime),\s*slope row (\d+): a in \[([^,]+),\s*([^\]]+)\]\.',line)
            require(m is not None,'retained row header');case,num,left,right=m.groups();ident=(case,int(num));require(ident not in seen,'duplicate retained row');seen.add(ident)
            rows.append(dict(case=case,number=int(num),a0=Q(left),a1=Q(right),cells=[]));active=True;continue
        if not active:continue
        row=rows[-1]
        if line.startswith('Each line has:'):continue
        if line.startswith('Exact integer sum of the last column:'):
            require('printed_sum'not in row,'duplicate retained row sum')
            m=re.match(r'Exact integer sum of the last column: (\d+);',line);require(m is not None,'retained sum record');row['printed_sum']=int(m[1]);continue
        parts=[v.strip()for v in line.split('|')];require(len(parts)==4,'unknown retained cell record '+line)
        left,right,formula,D=parts;tokens=formula.split();require(tokens[0]in ['B','S'],'retained cell kind')
        cell=dict(l=Q(left),v=Q(right),kind=tokens[0],D=int(D),line=lineno)
        if cell['kind']=='B':
            cell.update(choices=tokens[1],M=[int(x)for x in ''.join(tokens[2:]).split(',')]);require(len(cell['choices'])==16 and len(cell['M'])==7,'Bernstein witness size')
        else:require(len(tokens)==2,'shell witness size');cell['rate']=Q(tokens[1])
        row['cells'].append(cell)
    require(len(rows)==14,'retained row count')
    for kind,count in [('prime3',10),('no_prime',4)]:
        selected=[row for row in rows if row['case']==kind]
        require([x['number']for x in selected]==list(range(1,count+1)),'retained row indices')
        require(selected[0]['a0']==c.HEIGHT/c.R and selected[-1]['a1']==3,'retained root slope endpoints')
        require(all(a['a1']==b['a0']for a,b in zip(selected,selected[1:])),'retained root slope cover')
    reports=[];cells=coefs=Bcells=Scells=0;minimum=None
    for row in rows:
        rebuilt,branch=c.partition(row);printed=[row['cells'][0]['l']]+[z['v']for z in row['cells']]
        require(rebuilt==printed and printed[0]==0 and printed[-1]==c.U,'retained exact time partition')
        for index,cell in enumerate(row['cells']):
            left,right=cell['l'],cell['v'];require(left==printed[index] and 0<right-left<=Q(1,40),'retained cell coverage')
            h,f=c.profiles(row,left,right);rate=Q(1)if right<=c.R else Q(61,100)if right<=Q(99,100)else Q(11,25)
            if cell['kind']=='S':require(cell['rate']==rate,'retained shell rate');bound=rate*(right**7-left**7);Scells+=1
            else:
                codes='';Bcells+=1
                for lowbeta,B in c.BETAS:
                    poly,code,_=c.interval_polynomial(h,f,B,1-lowbeta,left,right);codes+=code;bern=c.bernstein(poly,left,right)
                    require(c.at(poly,left)>=0 and c.at(poly,right)>=0,'Bernstein endpoint sign')
                    require(all(Q(M,10**6)>=coef for M,coef in zip(cell['M'],bern)),'Bernstein coefficient enclosure');coefs+=7
                require(codes==cell['choices'],'retained slice choices')
                bound=(right-left)*Q(sum(cell['M']),10**6)
            delta=Q(cell['D'],10**8)-bound;require(delta>=0,'retained contribution enclosure')
            minimum=delta if minimum is None else min(minimum,delta);cells+=1
        total=sum(z['D']for z in row['cells']);target=Q(789,1000)if row['case']=='prime3'else Q(169,200)
        require(total==row['printed_sum']and Q(total,10**8)<target,'retained total or target')
        report=dict(case=row['case'],row=row['number'],cells=len(row['cells']),upper=str(Q(total,10**8)),target=str(target),strict_margin=str(target-Q(total,10**8)));reports.append(report)
        print(json.dumps({'family':key,**report}),flush=True)
    require(cells==1099,'retained cell census')
    gaps={'height_preparation':Q(5,38)-c.HEIGHT-Q(1,40000),'height_floor':c.HEIGHT-c.R/7,'unbounded_slope':Q(5,7)-c.R**7*Q(10,7),'prime_contradiction':Q(211,1000)-Q(4,5)**7,'no_prime_contradiction':1-Q(169,200)-Q(21,25)**7/2}
    require(all(v>0 for v in gaps.values()),'retained source or analytic tail gap')
    return dict(key=key,status='pass',cells=cells,bernstein_cells=Bcells,shell_cells=Scells,coefficient_inequalities=coefs,beta_intervals=16,rows=reports,minimum_cell_rounding_slack=str(minimum),analytic_gaps={name:str(v)for name,v in gaps.items()})
