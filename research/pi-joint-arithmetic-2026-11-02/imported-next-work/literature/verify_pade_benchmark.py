import math
from fractions import Fraction
import mpmath as mp
mp.mp.dps=180
rows=[]
for m in (3,5,7,9,13,21,31):
    coeff=[math.factorial(2*m-k)//(math.factorial(k)*math.factorial(m-k)) for k in range(1,m+1,2)]
    assert coeff[-1]==1
    assert max(coeff)==math.factorial(2*m-1)//math.factorial(m-1)
    assert all(x>y for x,y in zip(coeff,coeff[1:]))
    for a,b in ((13,2),(22,7),(710,113)):
        if math.gcd(a,b)>1: continue
        x=Fraction(-a*a,b*b)
        val=sum(c*x**j for j,c in enumerate(coeff))
        assert val.denominator==b**(m-1)
    T=2*mp.pi; omega=1j*T
    direct=sum(mp.mpf(c)*omega**(2*j) for j,c in enumerate(coeff))
    integ=mp.quad(lambda t:t**m*(1-t)**m*mp.exp(omega*t),[0,mp.mpf('.5'),1])
    exact_remainder=(-1)**(m+1)*omega**(2*m)/(2*math.factorial(m))*integ
    rel=abs((direct-exact_remainder)/direct)
    assert rel<mp.mpf('1e-120')
    chi=-integ/(mp.factorial(m)**2/mp.factorial(2*m+1))
    assert 1-mp.pi**2/(2*(2*m+3))<mp.re(chi)<1
    p=lambda Z:sum(mp.mpf(math.factorial(2*m-k)*math.factorial(m))/(math.factorial(2*m)*math.factorial(k)*math.factorial(m-k))*Z**k for k in range(m+1))
    g=lambda Z:p(Z)-p(-Z)
    root=mp.findroot(lambda y:mp.im(g(1j*y)),(T,T+mp.mpf('.1')),solver='secant')
    d=root-T
    rows.append({'m':m,'degree':m-1,'height':max(coeff),'root_above_2pi':mp.nstr(d,15),'log_root_distance':mp.nstr(mp.log(d),15),'remainder_relative_error':mp.nstr(rel,4)})
print(rows)
