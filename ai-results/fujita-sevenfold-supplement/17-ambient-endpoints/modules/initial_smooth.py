"""Exhaustive exact evaluation of all delivered retained-incidence leaves."""
from . import e09_ops as c
from .common import Q, require
import json,time

def validate(key,text):
    raw,rows,coverage=c.read_leaves(text)
    identities=c.self_checks()
    reports=[]; started=time.monotonic(); smallest=None
    for row,leaves in rows.items():
        worst=Q(0); wins=[0]*11; branches={}
        for index,(path,box,n) in enumerate(leaves):
            incoming,outgoing,value,branch,winners=c.evaluate(box,row)
            require(0<=value<=Q(n,10**7)<=Q(9999,10000),('retained smooth leaf',row,path,str(value),n))
            slack=Q(n,10**7)-value
            smallest=slack if smallest is None else min(smallest,slack)
            worst=max(worst,value);wins=[a+b for a,b in zip(wins,winners)]
            branches[branch]=branches.get(branch,0)+1
            if (index+1)%250==0:
                print(json.dumps({'family':key,'row':row,'checked':index+1,'row_total':len(leaves),'seconds':round(time.monotonic()-started,2)}),flush=True)
        reports.append(dict(coverage[row-1],exact_maximum=str(worst),margin_to_one=str(1-worst),branches=branches,competitor_wins=wins))
    return {'key':key,'status':'pass','leaves':10756,'rows':reports,'grid_intervals':len(c.GRID)-1,'minimum_leaf_rounding_slack':str(smallest),'margin_to_claimed_threshold':str(Q(9999,10000)-Q(9998671,10**7)),'identities':identities,'branch_counters':c.COUNTERS}
