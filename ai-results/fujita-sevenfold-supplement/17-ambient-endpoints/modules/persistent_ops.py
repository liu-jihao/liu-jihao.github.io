"""DEFAULT independent all-cells and full-domain rational audit.
Imports no generator code. This validates evidence, not a fact.
"""
from .common import Q as R
from pathlib import Path
from math import comb
from functools import lru_cache
import json,sys,time

def add(a,b):return tuple((a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(max(len(a),len(b))))
def scale(a,c):return tuple(x*c for x in a)
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):
    out=[R(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return tuple(out)
@lru_cache(maxsize=10000)
def power(a,n):
    out=(R(1),)
    for _ in range(n):out=mul(out,a)
    return out
def ev(a,t):return sum(x*t**i for i,x in enumerate(a))
@lru_cache(maxsize=25000)
def monomial_integral(n,l,z):return (z**n-l**n)/R(n)
def integral(a,l,z):return 7*sum(x*monomial_integral(i+1,l,z)for i,x in enumerate(a)if x)
def roots(a,l,z):return[-a[0]/a[1]]if len(a)>1 and a[1]and l<-a[0]/a[1]<z else[]
zero=(R(0),);tt=(R(0),R(1));r,u,v,T=R(9,10),R(26,25),R(27,25),R(29,25);eps=R(1,1000000)
cl,cr=R(16020581569,90194313216),R(8010290789,45097156608)
version=1
data={'NQ':(R(263,2000),R(1857,10000),'GEN','NQ'),'NE':(R(263,2000),R(189,1000),'GEN','NE'),
 'M':(R(263,2000),R(129,625),'GEN','IM'),'QT':(R(263,2000),R(2017,10000),'T','IQ'),
 'Q4':(R(13333,100000),R(2017,10000),'GEN','IQ'),'Q5':(R(1323,10000),R(2017,10000),'S5','IQ'),
 'QN':(R(263,2000),R(2017,10000),'NON','IQ')}
jcoeff={'IQ':[(2,1)],'IM':[(2,3),(3,-2)],'S5':[(2,10),(3,-20),(4,15),(5,-4)]}

def tighten(group,box):
    h0,hk,old,kind=data[group];al,ah,bl,bh,ll,lh,ml,mh,ql,qh,ol,oh=box
    if version==3:
        hpre=max(hk,h0+al*(u-r)-eps);bl=max(bl,hpre/u)
        ll=max(ll,h0+al*(v-r)-eps,hpre+bl*(v-u)-eps)
        ml=max(ml,h0+al*(T-r)-eps,hpre+bl*(T-u)-eps)
    ml=max(ml,3*ll-R(2,5)-eps);lh=min(lh,(mh+R(2,5)+eps)/3)
    h1=max(hk,h0+al*(u-r)-eps)
    ah=min(ah,(lh+eps-h0)/(v-r),(mh+eps-h0)/(T-r))
    bh=min(bh,(lh+eps-h1)/(v-u),(mh+eps-h1)/(T-u));bl=max(bl,h1/u)
    invalid=al>ah or bl>bh or ll>lh or ml>mh or h1*v/u>lh+eps
    return invalid,(al,ah,bl,bh,ll,lh,ml,mh,ql,qh,ol,oh),h1

def enclose_normalized(group,box):
    h0,hk,old,kind=data[group];ar,az,br,bz,lr,lz,mr,mz,qr,qz,orr,oz=box
    lmax=(R(157,250)+R(2,5)+eps)/3
    ll,lh=[R(7,30)+t*(lmax-R(7,30))for t in[lr,lz]]
    ml,mh=[(1-t)*max(R(3003,10000),3*l-R(2,5)-eps)+t*R(157,250)for l,t in[(ll,mr),(lh,mz)]]
    def cap_a(l,m):
        caps=[(l+eps-h0)/(v-r),(m+eps-h0)/(T-r)]
        caps.extend(((l+eps)*u/v-h0+eps)/(u-r)for _ in[0])
        caps.append(((m+eps)*u/T-h0+eps)/(u-r));return min(caps)
    al,ah=[(1-t)*h0/r+t*cap_a(l,m)for t,l,m in[(ar,ll,ml),(az,lh,mh)]]
    hl,hh=[max(hk,h0+(u-r)*a-eps)for a in[al,ah]]
    def value_b(h,l,m,t):return h/u+t*(min((l+eps-h)/(v-u),(m+eps-h)/(T-u))-h/u)
    lower=max(hl/u,min(value_b(h,ll,ml,t)for h in[hl,hh]for t in[br,bz]))
    kink=(3*lh-mh)/2+eps
    heights=[hl,hh]+([kink]if hl<kink<hh else[])
    upper=max(value_b(h,lh,mh,t)for h in heights for t in[br,bz])
    return(al,ah,lower,upper,ll,lh,ml,mh,qr,qz,orr,oz)

def raw(group,box,h1,t):
    h0,hk,old,kind=data[group];al,ah,bl,bh,ll,lh,ml,mh,ql,qh,ol,oh=box
    a=ah if t<r else al;b=bh if t<u else bl
    d=(mh+eps-ll)/(T-v)if t<v else 25*(ll-R(1,5));e=1+1000*(T-ml)**7
    return(h0-a*r,a),(h1-b*u,b),(ll-d*v,d),(ml-e*T,e)

def joint(D,F,H,q,kind,mid):
    positive=lambda p:power(p,6)if ev(p,mid)>0 else zero
    if ev(D,mid)<=0 or ev(F,mid)>=ev(D,mid):return zero
    M=F if ev(F,mid)>=ev(H,mid)else H;U=D if ev(D,mid)<=ev(H,mid)else H
    out=positive(sub(D,M))
    if ev(U,mid)<=ev(F,mid):return out
    for j,c in jcoeff[kind]:
        a=sub(D,scale(H,R(j)/(1-q)));b=R(j)*q/(1-q)-1
        if b:p=scale(sub(positive(add(a,scale(U,b))),positive(add(a,scale(F,b)))),1/b)
        else:p=scale(mul(power(a,5)if ev(a,mid)>0 else zero,sub(U,F)),6)
        out=add(out,scale(p,c))
    return out

def polynomial(group,box,name,H0,H1,F,l,z):
    h0,hk,old,kind=data[group];ql,qh,ol,oh=box[8:];D=sub(tt,F);R1=sub(H1,scale(F,qh));R0=sub(H0,scale(F,oh));mid=(l+z)/2
    if name.startswith('B'):
        rate=R(name[1:]);allowed=R(1)if mid<r else R(61,100)if mid<R(99,100)else R(11,25)if mid<u else R(17,50)
        assert rate==allowed
        return scale(power(tt,6),rate)
    if name=='P':return power(D,6)if ev(D,mid)>0 else zero
    if name in['J1','JD1','J0','JD0']:
        if name in['J0','JD0']:assert old=='S5'
        else:assert kind in['IQ','IM']
        return joint(tt,F,H0 if name=='J0'else H1,ol if name=='J0'else ql,'S5'if name=='J0'else kind,mid)if name in['J0','J1']else joint(D,zero,R0 if name=='JD0'else R1,ol if name=='JD0'else ql,'S5'if name=='JD0'else kind,mid)
    prefix,tag=name.split(':')
    if prefix=='0':Z,Y=tt,H0
    elif prefix=='1':Z,Y=tt,H1
    elif prefix=='D1':Z,Y=D,R1
    else:
        assert prefix=='D0'and old in['NON','S5'];Z,Y=D,R0
    q=4 if old=='T'and prefix in['0','D0']else 5 if prefix in['0','D0']else 3
    family='Q'if prefix in['1','D1']and kind in['IQ','NQ']else'M'if prefix in['1','D1']and kind=='IM'else'N'
    def ge(p):assert ev(p,l)>=0 and ev(p,z)>=0,(group,name,l,z,p)
    zm,ym=ev(Z,mid),ev(Y,mid);pp=lambda p,n:power(p,n)if ev(p,mid)>0 else zero
    if tag=='0':assert zm<=0 or ym>=zm;return zero
    if tag=='B1':
        assert zm>0 and(ym<=0 or family=='Q'and cl<=ym/zm<=cr)
        return power(Z,6)
    assert zm>0 and ym>0;ge(Z);ge(Y)
    if tag=='Nlo':
        assert family!='Q';ge(sub(Z,scale(Y,7)))
        return sub(power(Z,6),scale(mul(power(Z,6-q),power(Y,q)),R(7*q,q+1)**q))
    if tag in['Nhi','C']:
        assert family!='Q';ge(sub(scale(Y,7),Z));ge(sub(Z,Y))
        return scale(power(sub(Z,Y),6),R(7,6)**6)if tag=='C'else sub(power(Z,6),scale(mul(power(Z,6-q),power(add(scale(Z,6*q-1),scale(Y,7)),q)),R(1,6*(q+1))**q))
    if tag=='Q':
        assert family=='Q';ge(sub(scale(Z,cl),Y))
        return mul(power(sub(Z,scale(Y,2)),5),add(Z,scale(Y,10)))
    if tag=='A':
        ge(sub(Y,scale(Z,cr if family=='Q'else R(1,2))))
        return pp(sub(Z,Y),6)
    if tag in['H','Hlip']:
        assert family=='M'
        ge(sub(scale(Y,6),Z)if tag=='H'else sub(Z,scale(Y,6)))
        return scale(power(Z,6),R(256,729))if tag=='H'else sub(scale(power(Z,6),R(256,729)+3),scale(mul(power(Z,5),Y),18))
    if tag in['S','R']:
        assert family=='M'
        if tag=='R':ge(sub(Y,scale(Z,R(2,11))))
        coeff=[(1,1),(2,1),(3,-1)]if tag=='R'else[(1,1),(2,3),(3,-5),(4,2)]
    elif tag=='span':
        assert(prefix in['0','D0']and old=='S5')or family=='Q'and kind=='IQ'
        a,s=(1,5)if prefix in['0','D0']else(1,2)
        coeff=[(j,(-1)**(j+1)*(comb(a+s,j)-s*comb(a+s-1,j-1)))for j in range(1,a+s+1)]
    else:raise AssertionError(name)
    out=zero
    for j,c in coeff:out=add(out,scale(pp(sub(Z,scale(Y,j)),6),c))
    return out

def cell(group,box,h1,name,l,z):
    total=R(0);basic=sorted(set([l,z]+[x for x in[r,R(99,100),u,v]if l<x<z]))
    for aa,bb in zip(basic,basic[1:]):
        H0,H1,F1,F2=raw(group,box,h1,(aa+bb)/2)
        cuts=sorted(set([aa,bb]+sum((roots(p,aa,bb)for p in[H0,H1,F1,F2,sub(F1,F2)]),[])))
        for cc,dd in zip(cuts,cuts[1:]):
            mid=(cc+dd)/2;Y0=H0 if ev(H0,mid)>0 else zero;Y1=H1 if ev(H1,mid)>0 else zero;F=max([zero,F1,F2],key=lambda p:ev(p,mid));D=sub(tt,F)
            ql,qh,ol,oh=box[8:];R1=sub(Y1,scale(F,qh));R0=sub(Y0,scale(F,oh));expr=[]
            for Z,Y in[(tt,Y0),(tt,Y1),(D,R1),(D,R0)]:
                expr.extend([Z,Y]);expr.extend(sub(Y,scale(Z,c))for c in[R(1,7),R(1,6),R(1,5),R(1,4),R(1,3),R(1,2),R(1),R(2,11),cl,cr])
            old,kind=data[group][2:]
            joints=[]
            if kind in['IQ','IM']:joints.extend([(tt,F,Y1,ql,kind),(D,zero,R1,ql,kind)])
            if old=='S5':joints.extend([(tt,F,Y0,ol,'S5'),(D,zero,R0,ol,'S5')])
            for Z,P,H,q,k in joints:
                expr.extend([Z,P,H,sub(Z,P),sub(Z,H),sub(P,H)])
                for j,c in jcoeff[k]:
                    a=sub(Z,scale(H,R(j)/(1-q)));b=R(j)*q/(1-q)-1
                    expr.extend(add(a,scale(E,b))for E in[P,Z,H])
            fine=sorted(set([cc,dd]+sum((roots(p,cc,dd)for p in expr),[])))
            for x,y in zip(fine,fine[1:]):total+=integral(polynomial(group,box,name,Y0,Y1,F,x,y),x,y)
    return total
