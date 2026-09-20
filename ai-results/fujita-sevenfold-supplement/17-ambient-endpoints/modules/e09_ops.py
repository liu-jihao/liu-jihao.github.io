#!/usr/bin/env python3
"""Independent exact audit of the WRITTEN e09 certificate, not producer data.

The portable reader supplies the delivered text. Native project executions
use danus compute DEFAULT. No producer code or producer JSON is consumed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import time
from functools import lru_cache

from .common import Q as F, ENGINE, require, unwrap

ZERO, ONE = F(0), F(1)
R, U, T, Q = F(9, 10), F(27, 25), F(29, 25), F(13, 5)
H0, H1, HP = F(131, 1000), F(59, 250), F(3003, 10000)
ROOT = ((H0/R, F(2)), (H1/U, F(4)), (F(8021, 5600), F(200)))
INC = [(a, b, c, d) for a, b in [(ZERO, F(3,4)), (F(3,4), ONE)]
       for c, d in [(ZERO,F(1,2)), (F(1,2),F(5,8)),
                    (F(5,8),F(3,4)), (F(3,4),ONE)]]
GRID = sorted({ZERO, R, U, T, Q} | {F(j,50) for j in range(1,130)})
CSTAR, ETA = F(71,400), F(1,5000)
LINE_T, LINE_ZERO = (ONE, ZERO), (ZERO, ZERO)
COUNTERS = {}

def hit(name):
    COUNTERS[name] = COUNTERS.get(name, 0) + 1

def value(line, t):
    return line[0]*t + line[1]

def add(a,b):
    return (a[0]+b[0], a[1]+b[1])

def scale(a,s):
    return (a[0]*s, a[1]*s)

def sub(a,b):
    return (a[0]-b[0], a[1]-b[1])

def positive(x):
    return max(ZERO,x)

def pieces(lo,hi,*lines):
    points = {lo,hi}
    for a,b in lines:
        if a:
            z = -b/a
            if lo < z < hi:
                points.add(z)
                hit("interior_affine_roots")
    s = sorted(points)
    return zip(s,s[1:])

@lru_cache(maxsize=100000)
def r6(line,lo,hi):
    """Exact 7 integral (a t+b)_+^6, including negative/zero a."""
    a,b = line
    if not a:
        hit("r6_zero_slope")
        return 7*positive(b)**6*(hi-lo)
    return (positive(value(line,hi))**7-positive(value(line,lo))**7)/a

@lru_cache(maxsize=100000)
def f5(line,linear,lo,hi):
    """Independent binomial integration of 7(line)_+^5 * linear."""
    total = ZERO
    for l,v in pieces(lo,hi,line):
        if value(line,(l+v)/2)<=0:
            hit("f5_nonpositive_piece")
            continue
        length = v-l
        z, dz = value(line,l), line[0]*length
        w, dw = value(linear,l), linear[0]*length
        total += 7*length*sum((F(math.comb(5,j))*z**(5-j)*dz**j
                    *(w/F(j+1)+dw/F(j+2)) for j in range(6)), ZERO)
    return total

def capacity(i,degree,height,lo,hi):
    """Printed C_i, with exact threshold splitting; no sampling bound."""
    c = CSTAR if i==0 else F(2,11)
    result = ZERO
    for l,v in pieces(lo,hi,degree,height,sub(height,scale(degree,c))):
        mid = (l+v)/2
        d,h = value(degree,mid),value(height,mid)
        if d<=0:
            hit("capacity_zero_degree")
            continue
        if h<=0:
            hit("capacity_nonpositive_height")
            result += r6(degree,l,v)
            continue
        if i==0:
            if h<c*d:
                hit("q_balanced_branch")
                result += f5(sub(degree,scale(height,2)),add(degree,scale(height,10)),l,v)
            else:
                hit("q_axis_eta_branch")
                result += r6(sub(degree,height),l,v)+ETA*r6(degree,l,v)
        else:
            if h<c*d:
                hit("a_below_jump")
                terms = ((1,1),(2,3),(3,-5),(4,2))
            else:
                hit("a_above_jump")
                terms = ((1,1),(2,1),(3,-1))
            result += sum((coef*r6(sub(degree,scale(height,j)),l,v)
                           for j,coef in terms),ZERO)
    assert result>=0, ("negative integral",i,degree,height,lo,hi,result)
    return result

@lru_cache(maxsize=180000)
def independent(i,slope,lo,hi):
    anchor,height = (R,H0) if i==0 else (U,H1)
    return capacity(i,LINE_T,(slope,height-slope*anchor),lo,hi)

@lru_cache(maxsize=60000)
def prime(slope,lo,hi):
    p = (slope,HP-slope*T)
    result = ZERO
    for l,v in pieces(lo,hi,p):
        result += r6(LINE_T if value(p,(l+v)/2)<=0 else sub(LINE_T,p),l,v)
    return result

@lru_cache(maxsize=180000)
def partial(i,slope,pslope,upper,capped,lo,hi):
    anchor,height = (R,H0) if i==0 else (U,H1)
    h,p = (slope,height-slope*anchor),(pslope,HP-pslope*T)
    result = ZERO
    for l,v in pieces(lo,hi,h,p,sub(h,scale(p,upper))):
        mid = (l+v)/2
        hh = h if value(h,mid)>0 else LINE_ZERO
        pp = p if value(p,mid)>0 else LINE_ZERO
        q = pp
        if capped and upper*value(pp,mid)>value(hh,mid):
            hit("partial_capped")
            q = scale(hh,1/upper)
        result += capacity(i,sub(LINE_T,q),sub(hh,scale(q,upper)),l,v)
    assert result>=0
    return result

@lru_cache(maxsize=180000)
def joint(i,slope,pslope,alpha,lo,hi):
    anchor,height = (R,H0) if i==0 else (U,H1)
    h,p = (slope,height-slope*anchor),(pslope,HP-pslope*T)
    tail = ZERO
    js = (1,) if i==0 else (2,3)
    sums = {j:ZERO for j in js}
    for l,v in pieces(lo,hi,h,p,sub(h,LINE_T),sub(p,LINE_T),sub(p,h)):
        mid = (l+v)/2
        hh = h if value(h,mid)>0 else LINE_ZERO
        pp = p if value(p,mid)>0 else LINE_ZERO
        maximum = hh if value(hh,mid)>=value(pp,mid) else pp
        end = hh if value(hh,mid)<=mid else LINE_T
        tail += r6(sub(LINE_T,maximum),l,v)
        if value(pp,mid)>=value(end,mid):
            hit("joint_empty_slice")
            continue
        for j in js:
            c = j*alpha/(1-alpha)-1
            g = sub(LINE_T,scale(hh,F(j)/(1-alpha)))
            if c:
                term = (r6(add(scale(end,c),g),l,v)-r6(add(scale(pp,c),g),l,v))/c
            else:
                hit("joint_zero_z_slope")
                term = 6*f5(g,sub(end,pp),l,v)
            assert term>=0, ("negative B integral",i,j,alpha,l,v,term)
            sums[j] += term
    result = tail+sums[1] if i==0 else tail+min(2*sums[2],3*sums[2]-2*sums[3])
    assert result>=0, ("negative joint integral",i,result)
    return result

@lru_cache(maxsize=256)
def baseline(lo,hi):
    return hi**7-lo**7

@lru_cache(maxsize=80000)
def cell_minimum(row,a0,a1,ap,l,v):
    """Cache a pure exact cell calculation; every leaf still checks its sum."""
    b0,B0,b1,B1 = INC[row-1]
    bounds = [baseline(l,v),independent(0,a0,l,v),independent(1,a1,l,v),prime(ap,l,v),
              joint(0,a0,ap,b0,l,v),joint(1,a1,ap,b1,l,v),
              partial(0,a0,ap,B0,False,l,v),partial(0,a0,ap,B0,True,l,v),
              partial(1,a1,ap,B1,False,l,v),partial(1,a1,ap,B1,True,l,v)]
    if l>=T:
        bounds.append(r6((F(-52,73),F(125,73)),l,v))
    winner = min(range(len(bounds)),key=bounds.__getitem__)
    return winner,bounds[winner]

def evaluate(box,row):
    incoming,outgoing = ZERO,ZERO
    winners = [0]*11
    for l,v in zip(GRID,GRID[1:]):
        a0 = box[0][1 if v<=R else 0]
        a1 = box[1][1 if v<=U else 0]
        ap = box[2][1 if v<=T else 0]
        winner,bound = cell_minimum(row,a0,a1,ap,l,v)
        winners[winner] += 1
        if l<T:
            assert v<=T
            incoming += bound
        else:
            outgoing += bound
    a = incoming+(1-HP)**7
    b = min(incoming,F(499,500))+outgoing
    return incoming,outgoing,min(a,b),"prime" if a<=b else "point",winners

def parse_tree(source):
    """Parse every token, without eval or a producer assertion."""
    tokens = re.findall(r"\d+|[(),]",source)
    assert ''.join(tokens)==re.sub(r"\s+","",source)
    pos = 0
    def parse():
        nonlocal pos
        if tokens[pos].isdigit():
            n = int(tokens[pos]); pos += 1
            return n
        assert tokens[pos]=='('; pos += 1
        left = parse()
        assert tokens[pos]==','; pos += 1
        right = parse()
        assert tokens[pos]==')'; pos += 1
        return (left,right)
    tree = parse()
    assert pos==len(tokens)
    return tree

def read_leaves(text):
    raw = text.encode()
    compact = re.sub(r"\s+", "", text)
    pattern = r"Row(\d+):\(b_0,B_0,b_1,B_1\)=\(([^)]+)\)\.(\d+)leaves;largestupperinteger(\d+)\.(.*?)(?=Row|\Z)"
    matches = re.findall(pattern, compact)
    assert "".join(m.group(0) for m in re.finditer(pattern,compact)) == compact
    assert len(matches)==8
    all_leaves = {}
    coverage = []
    for row,inc,count,maximum,tree_text in matches:
        row=int(row)
        assert row not in all_leaves and 1 <= row <= 8
        assert tuple(F(x) for x in inc.split(','))==INC[row-1]
        tree = parse_tree(tree_text)
        leaves = []
        nodes,ties,depth = 0,0,0
        def descend(t,box,path):
            nonlocal nodes,ties,depth
            depth=max(depth,len(path))
            if isinstance(t,int):
                assert 0<=t<=9999000
                leaves.append((path,box,t))
                return
            nodes+=1
            widths=[(b-a)/(1+(a+b)/2) for a,b in box]
            k=max(range(3),key=widths.__getitem__)
            ties+=sum(w==widths[k] for w in widths)>1
            a,b=box[k]; m=(a+b)/2
            assert a<m<b
            left=list(box); right=list(box)
            left[k]=(a,m); right[k]=(m,b)
            descend(t[0],tuple(left),path+'L')
            descend(t[1],tuple(right),path+'R')
        descend(tree,ROOT,'')
        assert len(leaves)==int(count)
        assert max(x[2] for x in leaves)==int(maximum)
        assert nodes+1==len(leaves)
        assert len({x[0] for x in leaves})==len(leaves)
        all_leaves[row]=leaves
        coverage.append(dict(row=row,leaves=len(leaves),internal_nodes=nodes,
                             maximum_written=int(maximum),max_depth=depth,split_ties=ties,
                             tree_sha256=hashlib.sha256(tree_text.encode()).hexdigest()))
    assert sum(len(v) for v in all_leaves.values())==10756
    assert max(x['maximum_written'] for x in coverage)==9998671
    return raw,all_leaves,coverage

def qstring(q):
    return str(q)

def self_checks():
    checks={}
    def check(name,condition):
        assert condition,name
        checks[name]=True
    check('old_margin',F(5,38)-H0==F(11,19000))
    check('new_margin',F(10559,44675)-H1==F(157,446750))
    check('rounding_margin',F(1,200000)<min(F(11,19000),F(157,446750)))
    check('prime_margin',F(225,749)-HP==F(753,7490000)>F(1,10000))
    check('prime_slope',(HP-F(13,70))/(T-U)==ROOT[2][0])
    check('big_neighborhood',F(1,500)-F(7,10000)*F(7,6)**6>0)
    check('height_two_ratio',2*F(10,11)**5<F(4,3))
    excess=(1-2*CSTAR)**5*(1+10*CSTAR)-(1-CSTAR)**6
    check('q_switch_majorant',0<excess<ETA)
    old=R**7+F(4,3)*(R-H0)**7
    new=F(19,25)+(U-H1)**7/3+(U-2*H1)**7/7
    p=F(499,500)+(T-HP)**7/199
    check('unbounded_old',old<F(691,1000))
    check('unbounded_new',new<F(867,1000))
    check('unbounded_prime',p==F(1989490791196134952369386003613,
                               1990000000000000000000000000000)<F(9999,10000))
    check('new_tail_above_jump',H1/U>F(2,11))
    costs=[F(43,13),F(51,13),F(12879,2600),F(76841,13000),F(11234,1625),F(12949,1625)]
    check('full_point_cost',max(costs)==F(12949,1625))
    check('both_error_budget',max(costs)+F(1,1000)<8-F(1,100))
    check('global_written_margin',F(9999,10000)-F(9998671,10000000)==F(329,10000000))
    # Compare the binomial implementation with the printed primitive, on
    # all sign combinations and zero slopes. This tests clipping, too.
    def printed_f5(line,other,l,v):
        a,b=line; c,d=other
        total=ZERO
        for lo,hi in pieces(l,v,line):
            if value(line,(lo+hi)/2)<=0: continue
            if a:
                def primitive(t):
                    z=value(line,t)
                    return c*z**7/a**2+F(7,6)/a*(d-c*b/a)*z**6
            else:
                def primitive(t):
                    return 7*b**5*(c*t*t/2+d*t)
            total+=primitive(hi)-primitive(lo)
        return total
    for a in [-3,-1,0,1,3]:
        for b in [-2,0,2]:
            for c in [-2,0,3]:
                ln=(F(a),F(b)); ot=(F(c),F(1,3))
                check(f'f5_{a}_{b}_{c}',f5(ln,ot,F(-1,2),F(5,3))==printed_f5(ln,ot,F(-1,2),F(5,3)))
    # Exercise the zero-slope B_j formula on actual legal slice inputs.
    for i,alpha in [(0,F(1,2)),(1,F(1,3)),(1,F(1,4))]:
        result=joint(i,F(1,2),F(2),alpha,F(1),T)
        check(f'joint_singular_formula_{i}_{alpha}',result>=0)
    return dict(checks=checks,exact_values=dict(old_tail=str(old),new_tail=str(new),
                prime_tail=str(p),q_switch_excess=str(excess),endpoint=str(max(costs)+F(1,1000))))
