"""Written deterministic trees: exact independent leaf enclosures and coverage."""
import ast,re,json
from .common import Q,require,unique_insert
from . import reserve_ops as r, nonsmooth_ops as n, incidence_ops as i, payment_ops as p
KEYS=('initial-fivefold-reserve','initial-fivefold-nonsmooth','initial-fivefold-incidence','initial-fivefold-fivefold')

def take_tree(text):
    text=text.lstrip();require(text.startswith('('),'missing binary tree');depth=0
    for end,ch in enumerate(text):
        if ch=='(':depth+=1
        if ch==')':depth-=1
        if depth==0:break
    require(depth==0,'unterminated tree')
    token=re.sub(r'\s+','',text[:end+1]);require(re.fullmatch(r'[0-9(),]+',token)is not None,'malformed tree token')
    tree=ast.literal_eval(token)
    def valid(t):
        if type(t)is int:require(t>=0,'negative leaf');return
        require(type(t)is tuple and len(t)==2,'tree is not binary');valid(t[0]);valid(t[1])
    valid(tree);return tree,text[end+1:].lstrip()

def split_leaves(tree,box,path=''):
    if type(tree)is int:yield path,tree,box;return
    widths=[(b-a)/(1+(a+b)/2)for a,b in box];axis=max(range(len(box)),key=widths.__getitem__)
    a,b=box[axis];mid=(a+b)/2;require(a<mid<b,'noninterior tree split')
    left,right=list(box),list(box);left[axis]=(a,mid);right[axis]=(mid,b)
    yield from split_leaves(tree[0],tuple(left),path+'0');yield from split_leaves(tree[1],tuple(right),path+'1')

def named(text,pattern):
    rows=[]
    while text.strip():
        text=text.lstrip();m=re.match(pattern,text);require(m is not None,'malformed tree header '+text[:100])
        tree,text=take_tree(text[m.end():]);rows.append((m.groups(),tree))
    return rows

def reserve(key,text):
    base={}
    for _ in range(3):
        text=text.lstrip();m=re.match(r'K_([CAB]),\s*H=([\d/]+),\s*(\d+) leaves:\s*',text);require(m is not None,'reserve header')
        kind,height,count=m.groups();require(Q(height)==r.heights[kind],'reserve height')
        tree,text=take_tree(text[m.end():]);unique_insert(base,kind,(int(count),tree))
    require(set(base)=={'A','B','C'},'reserve family coverage')
    text=text.lstrip();header,text=text.split('\n',1);require(header.split()==['l0','v0','l1','v1','old','replacement'],'refinement column header')
    refinements={}
    while text.strip().strip('.'):
        text=text.lstrip();m=re.match(r'((?:\d+/\d+\s+){4})(\d+)\s+',text);require(m is not None,'malformed refinement row')
        box=tuple(map(Q,m[1].split()));old=int(m[2]);tree,text=take_tree(text[m.end():]);unique_insert(refinements,box,(old,tree))
    require(len(refinements)==15,'refinement coverage');reports=[]
    for kind in ['A','B','C']:
        count,tree=base[kind];root=(r.h0/r.r0,Q(2),r.heights[kind]/r.u,Q(4));leaves=list(r.leaves(tree,root))
        require(len(leaves)==count=={'A':47,'B':31,'C':147}[kind],'base leaf count')
        for box,N in leaves:require(0<=r.box_integral(kind,box)<=Q(N,10**6),'base reserve leaf ceiling')
        if kind=='C':
            require({b:N for b,N in leaves if N>=998000}=={b:row[0]for b,row in refinements.items()},'replacement domain exactly equals old high leaves')
            leaves=[entry for box,N in leaves for entry in (r.leaves(refinements[box][1],box)if box in refinements else [(box,N)])]
        values=[r.box_integral(kind,box)for box,N in leaves];Ns=[N for box,N in leaves]
        require(all(0<=v<=Q(N,10**6)<Q(499,500)for v,N in zip(values,Ns)),'final reserve leaf')
        require((len(leaves),max(Ns))=={'A':(47,997781),'B':(31,997970),'C':(163,997987)}[kind],'final reserve census')
        reports.append(dict(group=kind,base_leaves=count,leaves=len(leaves),maximum_upper=str(Q(max(Ns),10**6)),strict_margin=str(Q(499,500)-Q(max(Ns),10**6)),minimum_rounding_slack=str(min(Q(N,10**6)-v for v,N in zip(values,Ns)))))
    old=r.r0**7+r.cap_integral('Q',r.r0,r.h0,Q(2),r.r0,r.r);require(old<Q(631,1000),'reserve old slope tail')
    tails={'old':str(old)}
    for kind,bound in [('A','857/1000'),('B','802/1000'),('C','875/1000')]:
        v=Q(19,25)+r.cap_integral(kind,r.u,r.heights[kind],Q(4),r.u,r.r);require(v<Q(bound),'reserve new slope tail');tails[kind]=str(v)
    return dict(key=key,status='pass',groups=reports,replacement_subtrees=15,analytic_tails=tails)

def validate(key,text):
    if key=='initial-fivefold-reserve':return reserve(key,text)
    if key=='initial-fivefold-nonsmooth':
        n.scalar_checks();rows=named(text,r'(generic|four6|four7):\s*(\d+) leaves; maximum (\d+)/1000000\s*');den=10**6
        expected={'generic':(380,999894),'four6':(195,999856),'four7':(172,999856)}
    elif key=='initial-fivefold-incidence':
        i.scalar_checks();rows=named(text,r'(old|new) / (generic|four6|four7):\s*(\d+) leaves; maximum (\d+)/1000000\s*');den=10**6
        expected={'old/generic':(556,999897),'old/four6':(292,999897),'old/four7':(195,999836),'new/generic':(77,999873),'new/four6':(40,999673),'new/four7':(15,999640)}
    else:
        analytic=p.analytic_checks();rows=named(text,r'(masscap|partial):\s*(\d+) leaves; maximum (\d+)/10000000\s*');den=10**7
        expected={'masscap':(1773,9998961),'partial':(2315,9998999)}
    reports={}
    for fields,tree in rows:
        if key=='initial-fivefold-nonsmooth':
            name,count,maximum=fields;root=((n.h/n.r,Q(2)),(n.HEIGHT[name]/n.u,Q(4)),(n.H/n.T,Q(30)))
            evaluate=lambda box:n.enclosure(name,box)[0]
        elif key=='initial-fivefold-incidence':
            mode,kind,count,maximum=fields;name=mode+'/'+kind;root=((i.H/i.R,Q(2)),(i.HEIGHT[kind]/i.U,Q(4)),(i.AP_MIN,Q(60)))
            def evaluate(box):
                a,b,_=i.box_integrals(mode,kind,box);return min(a,Q(499,500))+b
        else:
            name,count,maximum=fields;hp=Q(31,100)if name=='masscap'else Q(3003,10000);theta=Q(1)if name=='masscap'else Q(43,50)
            root=((p.h0/p.r,Q(2)),(p.h1/p.u,Q(4)),((hp-Q(13,70))/(p.T-p.u),Q(200)))
            evaluate=lambda box:p.box_bound(box,hp,theta)[0]
        require(name in expected and (int(count),int(maximum))==expected[name],'tree inventory header')
        leaves=list(split_leaves(tree,root));require(len(leaves)==int(count)and max(x[1]for x in leaves)==int(maximum),'tree count and maximum')
        worst=Q(0);slack=None
        for index,(path,N,box)in enumerate(leaves,1):
            value=evaluate(box);require(0<=value<=Q(N,den)<Q(9999,10000),('leaf ceiling',key,name,path,str(value),N))
            delta=Q(N,den)-value;slack=delta if slack is None else min(slack,delta);worst=max(worst,value)
            if index%500==0:print(json.dumps({'family':key,'group':name,'checked':index,'total':len(leaves)}),flush=True)
        unique_insert(reports,name,dict(leaves=len(leaves),root=[[str(a),str(b)]for a,b in root],maximum_upper=str(Q(int(maximum),den)),maximum_exact=str(worst),strict_margin=str(Q(9999,10000)-Q(int(maximum),den)),minimum_rounding_slack=str(slack)))
    require(set(reports)==set(expected),'missing tree family')
    return dict(key=key,status='pass',groups=reports,analytic=analytic if key=='initial-fivefold-fivefold'else 'all source-specific scalar margins and unbounded slope inequalities checked')
