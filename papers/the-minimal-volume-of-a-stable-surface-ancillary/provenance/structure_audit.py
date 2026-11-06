"""One pass over the decompressed proof file large_mld_proof.vipr. Checks that every (asm) line has
integral coefficients and an integral right-hand side, and that every (lin), (rnd), and (uns) line
refers only to earlier lines, each at most once. Prints the SHA256 hash and the counts of the file,
and ends with the line "all_asm_integral_and_references_backward: true" if both properties hold.
Usage: python3 structure_audit.py large_mld_proof.vipr"""
import hashlib, sys, time
from collections import Counter
from fractions import Fraction

start = time.monotonic()
sha = hashlib.sha256()
size = 0
with open(sys.argv[1], 'rb') as f:
    def read_line():
        global size
        line = f.readline()
        sha.update(line)
        size += len(line)
        return line
    nv = nc = None
    while True:
        line = read_line()
        assert line, 'no DER section'
        if line.startswith(b'VAR '):
            nv = int(line.split()[1])
        elif line.startswith(b'CON '):
            nc = int(line.split()[1])
            for _ in range(nc):
                read_line()
        elif line.startswith(b'DER '):
            nd = int(line.split()[1])
            break
    kinds = Counter()
    bad_asm = bad_ref = 0
    for i in range(nd):
        line = read_line()
        lhs, rest = line.split(b'{', 1)
        reason = rest.split(b'}', 1)[0].split()
        kind = reason[0].decode()
        kinds[kind] += 1
        here = nc + i
        if kind in ('lin', 'rnd'):
            n = int(reason[1])
            assert len(reason) == 2 + 2 * n
            refs = [int(x) for x in reason[2::2]]
            if len(set(refs)) != len(refs) or any(not 0 <= j < here for j in refs):
                bad_ref += 1
        elif kind == 'uns':
            assert len(reason) == 5
            if any(not 0 <= int(j) < here for j in reason[1:]):
                bad_ref += 1
        elif kind == 'asm':
            t = lhs.split()
            assert t[1] in (b'L', b'G')
            n = int(t[3])
            assert len(t) == 4 + 2 * n
            if Fraction(t[2].decode()).denominator != 1:
                bad_asm += 1
            for j in range(n):
                if not 0 <= int(t[4 + 2 * j]) < nv or Fraction(t[5 + 2 * j].decode()).denominator != 1:
                    bad_asm += 1
        else:
            raise ValueError(kind)
    while True:
        line = read_line()
        if not line:
            break
        assert not line.strip()
print('bytes', size)
print('sha256', sha.hexdigest())
print('variables', nv, 'constraints', nc, 'derivations', nd)
print('kinds', dict(sorted(kinds.items())))
print('asm lines with a non-integral number', bad_asm, 'lines with a bad reference', bad_ref)
print('seconds', round(time.monotonic() - start, 1))
print('all_asm_integral_and_references_backward:', 'true' if bad_asm == bad_ref == 0 else 'false')
sys.exit(0 if bad_asm == bad_ref == 0 else 1)
