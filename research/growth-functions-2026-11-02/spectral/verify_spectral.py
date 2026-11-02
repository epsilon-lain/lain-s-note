"""Exact finite Fourier growth-family kernel checks; no trigonometric floats."""
from pathlib import Path
from random import Random
from math import ceil
from itertools import zip_longest
import json

def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return p

def multiply(p,q):
    out=[0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):out[i+j]+=a*b
    return trim(out)

def subtract(p,q):
    out=[0]*max(len(p),len(q))
    for i,a in enumerate(p):out[i]+=a
    for i,a in enumerate(q):out[i]-=a
    return trim(out)

def divide_unit(p,r):
    """Monic division by w**r-1."""
    work=list(p);quo=[0]*max(1,len(p)-r)
    for j in range(len(p)-1,r-1,-1):
        c=work[j];quo[j-r]=c;work[j]=0;work[j-r]+=c
    return trim(quo),trim(work[:r] or [0])

def residue(p,r):
    out=[0]*r
    for j,c in enumerate(p):out[j%r]+=c
    return trim(out)

def run():
    rng=Random(20261102);count=0;generated=0
    for r in range(1,25):
        divisor=[-1]+[0]*(r-1)+[1]
        for N in range(0,55):
            for _ in range(3):
                p=[rng.randint(-20,20) for j in range(N+1)]
                q,rem=divide_unit(p,r)
                assert rem==residue(p,r)
                assert trim(p)==trim([a+b for a,b in zip_longest(multiply(q,divisor),rem,fillvalue=0)])
                count+=1
            if N>=r:
                q=[rng.randint(-20,20) for j in range(N-r+1)]
                p=multiply(q,divisor)
                recovered,rem=divide_unit(p,r)
                assert recovered==trim(q) and rem==[0]
                assert max(map(abs,recovered))<=ceil(len(p)/r)*max(map(abs,p))
                # The coordinates in the monomial-times-divisor basis are recovered q.
                rebuilt=[0]
                for j,c in enumerate(recovered):
                    term=[0]*j+[c*x for x in divisor]
                    rebuilt=[a+b for a,b in zip_longest(rebuilt,term,fillvalue=0)]
                assert trim(rebuilt)==trim(p)
                generated+=1
    p3=[-1,0,0,1];p2=[-1,0,1]
    assert subtract(p3,[0]+p2)==[-1,1]
    lcm=[-1,-1,0,1,1]
    assert multiply(lcm,[-1,1])==multiply(p3,p2)
    assert divide_unit(lcm,2)[1]==[0] and divide_unit(lcm,3)[1]==[0]
    assert divide_unit(lcm,6)[1]!=[0]
    out={'status':'PASS','random_division_and_residue_checks':count,
         'integer_basis_reconstruction_and_height_checks':generated,
         'two_skeleton_identity_checks':4,'uses_pi_digits':False,
         'scope':'Exact finite Fourier kernel arithmetic only; no pi approximation claim.'}
    Path(__file__).with_name('checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':run()
