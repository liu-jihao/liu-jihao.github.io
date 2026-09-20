"""Independent exact replay of 98b's complete indexed written certificate.

Finite question: Does the complete printed tree cover the root, do all empty
tests follow from the printed constraints, and does each named bound apply on
its entire cell and integrate below the printed numerator? No producer code
or producer result is imported. Run only with danus compute DEFAULT.
The polynomial engine is the previous independent sub10 c2bf cold checker;
the target parser, old codimension, uniform heights, and method criteria here
are derived from 98b itself. This imports no producer code or result.
"""
from .common import Q as F, record_lines, require
from pathlib import Path
from collections import Counter
import hashlib
import json
import re
import sys

r, u, v, T = F(9,10), F(26,25), F(27,25), F(29,25)
eps = F(1,10**6)
old_heights = {'Q': F(13333,100000), 'M': F(1323,10000)}
hk = {'Q': F(2017,10000), 'M': F(129,625)}
cminus, cplus = F(16020581569,90194313216), F(8010290789,45097156608)

class Poly:
    def __init__(self, x=0):
        a = list(x.a) if isinstance(x, Poly) else list(x) if isinstance(x, (tuple,list)) else [F(x)]
        while len(a)>1 and not a[-1]: a.pop()
        self.a = tuple(F(y) for y in a)
    def __add__(self, other):
        b = Poly(other).a
        return Poly([(self.a[i] if i<len(self.a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(self.a),len(b)))])
    __radd__ = __add__
    def __neg__(self): return Poly([-x for x in self.a])
    def __sub__(self, other): return self + -Poly(other)
    def __rsub__(self, other): return Poly(other) + -self
    def __mul__(self, other):
        b = Poly(other).a
        c = [F(0)]*(len(self.a)+len(b)-1)
        for i,x in enumerate(self.a):
            for j,y in enumerate(b): c[i+j] += x*y
        return Poly(c)
    __rmul__ = __mul__
    def __truediv__(self, other): return Poly([x/F(other) for x in self.a])
    def __pow__(self, n):
        z = Poly(1)
        for _ in range(n): z = z*self
        return z
    def at(self, t):
        z = F(0)
        for x in reversed(self.a): z = z*t+x
        return z
    def integral(self, a, b):
        return sum((x*(b**(j+1)-a**(j+1))/F(j+1) for j,x in enumerate(self.a)),F(0))

t = Poly([0,1])
class Split(Exception): pass

class Interval:
    """Every branch is proved by affine sign, splitting at its exact root."""
    def __init__(self, l, z): self.l,self.z,self.mid=l,z,(l+z)/2
    def sign(self, x):
        x=Poly(x)
        assert len(x.a)<=2, ('nonaffine branch',x.a)
        if len(x.a)==2 and x.a[1]:
            root=-x.a[0]/x.a[1]
            if self.l<root<self.z: raise Split(root)
        q=x.at(self.mid)
        return (q>0)-(q<0)
    def pos(self, x): return Poly(x) if self.sign(x)>0 else Poly(0)
    def minimum(self,a,b): return Poly(a) if self.sign(Poly(a)-b)<=0 else Poly(b)
    def maximum(self,a,b): return Poly(a) if self.sign(Poly(a)-b)>=0 else Poly(b)
    def require(self, expression, label):
        assert self.sign(expression)>=0, ('invalid capacity regime',label,str(self.l),str(self.z))

def parse(text):
    families={}; rationals={}; labels={}; used_r=Counter(); used_p=Counter()
    family=mode=tag=None; parsed_lines=[]
    def rational(index):
        i=int(index); used_r[i]+=1
        return rationals[i]
    for line_no,line in enumerate(['RATIONAL DICTIONARY: index:value'] + record_lines(text),1):
        if line == 'RATIONAL DICTIONARY: index:value':
            mode='rational'; parsed_lines.append(line_no); continue
        if line == 'POLYNOMIAL DICTIONARY: index:label':
            assert len(rationals)==5891 and set(rationals)==set(range(5891))
            mode='polynomial'; parsed_lines.append(line_no); continue
        if line.startswith('FAMILY '):
            family=line.split()[1]
            assert family in ('Q','M') and family not in families
            assert len(labels)==23 and set(labels)==set(range(23))
            families[family]={'splits':{},'empty':{},'leaves':{},'cells':{},'sums':{}}
            mode='family';parsed_lines.append(line_no)
            continue
        if mode is None: continue
        parsed_lines.append(line_no)
        if not line.strip(): continue
        if mode=='rational':
            a,b=line.split(':'); i=int(a)
            assert i==len(rationals) and str(F(b))==b
            rationals[i]=F(b)
            continue
        if mode=='polynomial':
            a,b=line.split(':',1); i=int(a)
            assert i==len(labels)
            labels[i]=b
            continue
        f=families[family]; a=line.split()
        if a[0]=='ROOT':
            assert len(a)==2 and 'root' not in f
            f['root']=tuple(rational(i) for i in a[1].split(','))
            assert len(f['root'])==8
        elif a[0]=='S':
            assert len(a)==4 and a[1] not in f['splits']
            assert 0<=int(a[2])<=3
            f['splits'][a[1]]=(int(a[2]),rational(a[3]),line_no)
        elif a[0]=='E':
            assert len(a)==2 and a[1] not in f['empty']
            f['empty'][a[1]]=None
        elif a[0]=='L':
            assert len(a)==6 and a[1] not in f['leaves'] and a[2] in ('29','27mass')
            tag=a[1]
            f['leaves'][tag]={'height':rational(a[3]),'method':a[2],'numerator':int(a[4]),'total':rational(a[5]),'line':line_no}
            f['cells'][tag]=[];f['sums'][tag]=0
        elif a[0]=='C':
            assert len(a)==4
            left=f['cells'][tag][-1][1] if f['cells'][tag] else F(0)
            polyindex=int(a[2]); used_p[polyindex]+=1
            f['cells'][tag].append((left,rational(a[1]),labels[polyindex],int(a[3]),line_no))
            f['sums'][tag]+=int(a[3])
        else: raise AssertionError(('unparsed appendix line',line_no,line))
    assert rationals[0]==0
    assert set(used_r)==set(rationals) and set(used_p)==set(labels)
    return families

def clipped(box,family):
    h0=old_heights[family]
    al,ah,bl,bh,ll,lh,tl,th=box
    height=max(hk[family],h0+al*(u-r)-eps)
    ahi=min(ah,(lh+eps-h0)/(v-r))
    bhi=min(bh,(lh+eps-height)/(v-u))
    blo=max(bl,height/u)
    return height,ahi,blo,bhi

def empty(box,family):
    height,ahi,blo,bhi=clipped(box,family)
    return box[0]>ahi or blo>bhi or height*v/u>box[5]+eps

def verify_tree(f,family):
    h0=old_heights[family]
    expected_root=(h0/r,(F(2,5)+eps-h0)/(v-r),hk[family]/u,(F(2,5)+eps-hk[family])/(v-u),F(7,30),F(2,5),F(0),F(1))
    assert f['root']==expected_root,('root',family)
    pending={'R':f['root']}
    for tag,(axis,mid,line) in f['splits'].items():
        assert tag in pending,('split absent',family,line)
        box=pending.pop(tag)
        assert box[2*axis] < mid < box[2*axis+1]
        assert mid==(box[2*axis]+box[2*axis+1])/2
        left,right=list(box),list(box)
        left[2*axis+1]=mid;right[2*axis]=mid
        pending[tag+'0']=tuple(left);pending[tag+'1']=tuple(right)
    assert set(f['leaves']).isdisjoint(f['empty'])
    assert set(pending)==set(f['leaves'])|set(f['empty'])
    assert set(f['cells'])==set(f['leaves'])==set(f['sums'])
    for tag in f['empty']:
        box=pending[tag];f['empty'][tag]=box
        assert empty(box,family),('false empty',family,tag)
    for tag,leaf in f['leaves'].items():
        box=pending[tag];leaf['box']=box
        assert not empty(box,family)
        assert leaf['height']==clipped(box,family)[0]
        assert box[6]<1 and box[6]<box[7]<=1

def joint(iv,family,D,Fp,H,theta):
    if iv.sign(D)<=0: return Poly(0)
    if iv.sign(H)<=0: return iv.pos(D-Fp)**6
    upper=iv.minimum(D,H)
    result=iv.pos(D-iv.maximum(Fp,H))**6
    if iv.sign(upper-Fp)<=0: return result
    pairs=[(2,1)] if family=='Q' else [(2,3),(3,-2)]
    for j,c in pairs:
        A=D-j*H/(1-theta)
        B=j*theta/(1-theta)-1
        if B:
            result+=c*(iv.pos(A+B*upper)**6-iv.pos(A+B*Fp)**6)/B
        else:
            result+=6*c*(iv.pos(A)**5)*(upper-Fp)
    return result

def named_poly(iv,family,box,name):
    h0=old_heights[family]
    # Fix all geometric anchor sides independently of the printed cell cuts.
    for anchor in (r,F(99,100),u,v):
        if iv.l<anchor<iv.z: raise Split(anchor)
    al,ah,bl,bh,ll,lh,tl,th=box
    height,ahi,blo,bhi=clipped(box,family)
    Y0=iv.pos(h0+(ahi if iv.mid<r else al)*(t-r))
    Y1=iv.pos(height+(bhi if iv.mid<u else blo)*(t-u))
    forced=iv.pos(ll+12*(t-v)) if iv.mid<v else ll+25*(ll-F(1,5))*(t-v)
    remain=iv.pos(t-forced)
    debit=Y1-th*forced
    if name.startswith('B'):
        constant=F(name[1:])
        baseline=F(1) if iv.mid<r else F(61,100) if iv.mid<F(99,100) else F(11,25) if iv.mid<u else F(17,50)
        assert constant>=baseline,('baseline',name)
        return constant*t**6
    if name=='P': return remain**6
    if name=='J': return joint(iv,family,t,forced,Y1,tl)
    if name=='JD': return joint(iv,family,remain,Poly(0),debit,tl)
    prefix,kind=name.split(':')
    if prefix=='0': Z,Y,q=t,Y0,5
    elif prefix=='1': Z,Y,q=t,Y1,3
    elif prefix=='D': Z,Y,q=remain,iv.pos(debit),3
    else: raise AssertionError(('prefix',prefix))
    if iv.sign(Z)<=0: return Poly(0)
    if iv.sign(Y)<=0: return Z**6
    if kind=='B1': return Z**6
    if kind=='0':
        iv.require(Y-Z,name)
        return Poly(0)
    if kind=='Nlo':
        iv.require(Z/7-Y,name)
        return Z**6-F(7*q,q+1)**q*Z**(6-q)*Y**q
    if kind in ('Nhi','C'):
        iv.require(Y-Z/7,name);iv.require(Z-Y,name)
        if kind=='C': return F(7,6)**6*(Z-Y)**6
        return Z**6-Z**(6-q)*((6*q-1)*Z+7*Y)**q/F(6*(q+1))**q
    if kind=='A':
        iv.require(Y-(cplus if family=='Q' and prefix!='0' else F(1,2))*Z,name)
        iv.require(Z-Y,name)
        return (Z-Y)**6
    if kind=='Q':
        assert family=='Q' and prefix!='0'
        iv.require(cminus*Z-Y,name)
        return (Z-2*Y)**5*(Z+10*Y)
    if kind in ('span','R'):
        if family=='M': iv.require(Y-F(2,11)*Z,name)
        return iv.pos(Z-Y)**6+iv.pos(Z-2*Y)**6-iv.pos(Z-3*Y)**6
    if kind=='S':
        assert family=='M'
        return iv.pos(Z-Y)**6+3*iv.pos(Z-2*Y)**6-5*iv.pos(Z-3*Y)**6+2*iv.pos(Z-4*Y)**6
    if kind=='H':
        assert family=='M';iv.require(Y-Z/6,name)
        return F(256,729)*Z**6
    if kind=='Hlip':
        assert family=='M';iv.require(Z/6-Y,name)
        return (F(256,729)+3)*Z**6-18*Z**5*Y
    raise AssertionError(('unknown polynomial',name))

def integrate_cell(family,box,name,left,right):
    pending=[(left,right)]
    value=F(0);pieces=[]
    while pending:
        l,z=pending.pop()
        try: poly=named_poly(Interval(l,z),family,box,name)
        except Split as e:
            mid=e.args[0]
            pending.extend([(l,mid),(mid,z)])
            continue
        value+=7*poly.integral(l,z)
        pieces.append((l,z,poly.a))
    return value,pieces

def constants():
    rho=F(2237,1000)
    assert rho**2==F(5004169,10**6)>5
    assert rho**3==F(11194326053,10**9)>9
    assert (r*rho-1)/(7*rho-8)==F(10133,76590)
    assert F(10133,76590)-old_heights['M']==F(143,76590000)>eps
    assert (2*r-1)/(14-8)==F(2,15)
    assert F(2,15)-old_heights['Q']==F(1,300000)>eps
    table=[(F(5,4),F(92,125),F(150,743)),(F(623,500),F(719,1000),F(7396,36025)),(F(33,25),F(76,125),F(233,1155)),(F(1317,1000),F(71,125),F(9242,44675)),(F(283,200),F(379,1000),F(1179,5710))]
    for a,p,bound in table:
        assert 7*a+p-8>0 and (a*u-1)/(7*a+p-8)==bound
    assert min(row[2] for row in table[:3])-hk['Q']==F(73,2310000)>eps
    assert min(row[2] for row in table[3:])-hk['M']==F(57,713750)>eps
    assert 32*127**7==17052027525296096<18794531250000000==11*150**7
    # Balanced minus endpoint is positive at 0<c<c7, negative above c7.
    cross=lambda c:(1-2*c)**5*(1+10*c)-(1-c)**6
    assert cross(cminus)>0 and cross(cplus)<0
    return {'rho_square':str(rho**2),'old_Q_height_gap':str(F(2,15)-old_heights['Q']),'old_M_height_gap':str(F(10133,76590)-old_heights['M']),'Q_height_gap':str(F(73,2310000)),'M_height_gap':str(F(57,713750)),'crossover_bracket_signs':['positive','negative']}

