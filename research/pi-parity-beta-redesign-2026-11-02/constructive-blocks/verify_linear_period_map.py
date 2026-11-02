from math import comb,gcd
from fractions import Fraction as F
from pathlib import Path
import json

def odd(n):
    n=abs(n)
    while n and n%2==0:n//=2
    return n

def run():
    checks=targets=0
    for m in range(1,13):
        A,C=comb(2*m,m),comb(4*m,2*m)
        S,R,T,D=2*m-1,6*m+1,2**(6*m+2),A*C
        al0=F(A,64**m);al1=F(R,S)*al0;b1=F(-1,S*C)
        for h0 in range(-15,16):
            for h1 in range(-15,16):
                if gcd(h0,h1)!=1:continue
                v=S*h0+R*h1
                if not v:continue
                alpha=h0*al0+h1*al1;beta=h1*b1
                assert -4*beta/alpha==F(T*h1,D*v)
                g=gcd(T*h1,D*v)
                assert (T*D*S)%g==0
                assert odd(D*gcd(h1,S))%odd(g)==0
                assert odd(gcd(h1,D)*gcd(h1,S))%odd(g)==0
                checks+=1
        for p in range(-5,6):
            for q in range(1,8):
                h1=D*S*p;h0=T*q-R*D*p;v=S*h0+R*h1
                assert F(T*h1,D*v)==F(p,q)
                targets+=1
    result=dict(status='passed',primitive_linear_blocks_checked=checks,
                arbitrary_target_reconstructions=targets,
                scope='Exact finite audit of the proved linear period and gcd identities; no approximation limit inferred.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':run()
