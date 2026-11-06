"""Removes constraints that no derivation of a VIPR 1.0 file refers to, and renumbers the
references (the reasons and the trailing index of each derivation). A reference to a removed
constraint stops the script. Usage: drop_unused_constraints.py in.vipr out.vipr i1 i2 ...
The output has to be checked again with viprchk; this script is not part of any argument."""
import sys
src,dst=sys.argv[1:3]; drop=sorted(int(x) for x in sys.argv[3:])
with open(src,'rb') as f,open(dst,'wb') as g:
    while True:
        line=f.readline(); assert line
        if line.startswith(b'CON '):
            h=line.split(); nc=int(h[1]); g.write(b'CON %d %s\n'%(nc-len(drop),h[2]))
            new=[]; k=0
            for i in range(nc):
                c=f.readline()
                if k<len(drop) and drop[k]==i: new.append(None); k+=1
                else: new.append(i-k); g.write(c)
            assert k==len(drop)
            continue
        if line.startswith(b'DER '):
            nd=int(line.split()[1]); g.write(line); break
        g.write(line)
    def m(j):
        if j>=nc: return j-len(drop)
        assert new[j] is not None,j
        return new[j]
    for i in range(nd):
        line=f.readline(); lhs,rest=line.split(b'{',1); body,after=rest.split(b'}',1)
        t=body.split(); kind=t[0]
        if kind in(b'lin',b'rnd'):
            n=int(t[1]); assert len(t)==2+2*n
            for j in range(n): t[2+2*j]=b'%d'%m(int(t[2+2*j]))
        elif kind==b'uns':
            assert len(t)==5
            for j in range(1,5): t[j]=b'%d'%m(int(t[j]))
        elif kind!=b'asm': raise ValueError(kind)
        a=int(after.strip()); a=a if a<0 else m(a)
        g.write(lhs+b'{ '+b' '.join(t)+b' } %d\n'%a)
    rest=f.read(); assert not rest.strip(),rest[:80]
