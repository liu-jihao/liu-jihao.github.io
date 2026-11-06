"""Compress the immutable definitive proof artifact; preserve the original.
Run through [job runner] to cap CPU and I/O.
"""
from pathlib import Path
import gzip,hashlib,time,json
R=Path(__file__).resolve().parent
source=R/'xhigh2_scip_exact_45.complete.vipr.opt';dest=Path(str(source)+'.gz')
expected='42a2669f4670ad31f1ae1869db820462ae847312b2ec7a530261bd614e743df5'
start=time.monotonic();h=hashlib.sha256()
with source.open('rb') as inp,dest.open('wb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,compresslevel=6,mtime=0) as out:
    for block in iter(lambda:inp.read(1<<20),b''):h.update(block);out.write(block)
assert h.hexdigest()==expected
check=hashlib.sha256()
with gzip.open(dest,'rb') as inp:
    for block in iter(lambda:inp.read(1<<20),b''):check.update(block)
assert check.hexdigest()==expected
print({'compressed_bytes':dest.stat().st_size,'uncompressed_sha256':expected,
       'seconds':time.monotonic()-start},flush=True)
