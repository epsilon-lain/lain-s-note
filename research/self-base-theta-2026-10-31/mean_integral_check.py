from fractions import Fraction as F
from math import comb, lcm, gcd
import json
import mpmath as mp
mp.mp.dps=250
rows=[]
for n in [2,4,8,16,32,64]:
    a=F(0)
    for j in range(n+1):
        k=n+j
        a+=(-1)**j*comb(n,j)*sum((F((-1)**r,2*(k-1-r)+1) for r in range(k)),F(0))
    D=lcm(*range(1,4*n,2))
    assert (D*a).denominator==1
    A=int(D*a); B=(-1)**n*D*2**(n-2); g=gcd(A,B)
    J=mp.quad(lambda x:(x*x*(1-x*x))**n/(1+x*x),[0,1/mp.sqrt(2),1])
    direct=mp.mpf(A)/D+(-1)**n*mp.mpf(2)**(n-2)*mp.pi
    assert abs(direct/J-1)<mp.mpf("1e-150")
    rows.append(dict(n=n,gcd=str(g),q=str(abs(B)//g),
                     normalized_error_numeric=mp.nstr(abs(B)*mp.mpf(D)*J/g**2,30)))
for beta in [mp.sqrt(2),mp.sqrt(3),mp.pi]:
    for m in map(mp.mpf,[".001",".01",".1",".4"]):
        c=mp.mpf(".5")
        H=mp.acosh(c/m)/(2*mp.log(beta))
        value=mp.quad(lambda t:1/(m*mp.cosh(2*t*mp.log(beta))),[-H,0,H])
        exact=mp.acos(m/c)/(m*mp.log(beta))
        assert abs(value/exact-1)<mp.mpf("1e-150")
print(json.dumps(dict(status="Numerical checks, not interval certificates; P2848 unresolved",
                     excursion_integral_checks=12,pi_integral_samples=rows),indent=2))
