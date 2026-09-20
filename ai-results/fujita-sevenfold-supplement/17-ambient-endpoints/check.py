#!/usr/bin/env python3
"""Check every ambient-endpoint data object using exact standard-library arithmetic.
Usage: python3 check.py --data data --report report.json
No candidate generator, private state, network service or cached verdict is used.
"""
from pathlib import Path
import argparse,gc,hashlib,json,platform,sys,time,traceback
from modules.common import ENGINE,require
from modules import compact,persistent,sources,initial_smooth,surfaces,quintuple,trees,retained

KEYS=tuple(sorted(compact.SPECS))+('retained-parameter','initial-fivefold-reserve','initial-fivefold-nonsmooth','initial-fivefold-incidence','initial-fivefold-fivefold','initial-fivefold-smooth','source-noncontained26','source-regular26','source-noncontained-initial','source-initial-fourfold','source-initial-threefold','source-initial-surface','source-quintuple-surface','distinct-divisors','persistent-fourfold','persistent-small-surface','persistent-fivefold')

def dispatch(key,text):
    if key in compact.SPECS:return compact.validate(key,text)
    if key in persistent.GROUPS:return persistent.validate(key,text)
    if key in sources.KEYS:return sources.validate(key,text)
    if key=='initial-fivefold-smooth':return initial_smooth.validate(key,text)
    if key=='source-initial-surface':return surfaces.surface(key,text)
    if key=='source-quintuple-surface':return quintuple.validate(key,text)
    if key=='distinct-divisors':return surfaces.distinct(key,text)
    if key in trees.KEYS:return trees.validate(key,text)
    if key=='retained-parameter':return retained.validate(key,text)
    raise ValueError('Unknown data family: '+key)

def inventory(data):
    require(len(KEYS)==25 and len(set(KEYS))==25,'checker family inventory')
    found={p.name for p in data.iterdir()if p.is_file()}
    expected={k+'.txt'for k in KEYS}
    require(found==expected,'missing or unexpected data files: '+str(sorted(found^expected)))

def clear_caches():
    # Independent families need no shared arithmetic cache; keep peak memory bounded.
    for name,module in tuple(sys.modules.items()):
        if name.startswith('modules.') and module is not None:
            for obj in tuple(vars(module).values()):
                if callable(obj)and hasattr(obj,'cache_clear'):obj.cache_clear()
    gc.collect()

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data',type=Path,required=True);ap.add_argument('--report',type=Path,required=True)
    args=ap.parse_args(argv);started=time.monotonic()
    report={'schema':'ambient-endpoints-check/1','section':'17-ambient-endpoints','status':'running','scope':'complete','engine':ENGINE,'python':platform.python_version(),'families':[]}
    try:
        require(__debug__,'Run without -O: exact checker assertions must be enabled')
        inventory(args.data)
        for key in KEYS:
            payload=(args.data/(key+'.txt')).read_bytes()
            result=dispatch(key,payload.decode('utf-8'))
            result.update(id='17-ambient-endpoints/'+key,key=key,status='pass',data_sha256=hashlib.sha256(payload).hexdigest())
            report['families'].append(result);clear_caches()
            report['seconds']=round(time.monotonic()-started,3)
            args.report.parent.mkdir(parents=True,exist_ok=True)
            args.report.write_text(json.dumps(report,indent=2)+'\n')
            print(json.dumps({'family':key,'status':'pass','seconds':round(time.monotonic()-started,3)}),flush=True)
        require({r['key']for r in report['families']}==set(KEYS),'incomplete result coverage')
        report['status']='pass';code=0
    except Exception as exc:
        report['status']='fail';report['error']=str(exc);report['traceback']=traceback.format_exc();code=1
        print(report['traceback'],file=sys.stderr)
    report['seconds']=round(time.monotonic()-started,3)
    args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(report,indent=2)+'\n')
    return code
if __name__=='__main__':sys.exit(main())
