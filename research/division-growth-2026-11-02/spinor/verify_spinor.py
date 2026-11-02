"""Exact algebra and arithmetic checks. Requires SymPy; no pi digits as inputs."""
from fractions import Fraction
from pathlib import Path
import json
import sympy as s

t = s.symbols('t', real=True)
p, q, dp, dq = s.symbols('p q dp dq', real=True)
D = p*p+q*q
X, Y = (p*p-q*q)/D, 2*p*q/D
dx = s.diff(X,p)*dp+s.diff(X,q)*dq
dy = s.diff(Y,p)*dp+s.diff(Y,q)*dq
assert s.cancel(X*X+Y*Y-1) == 0
assert s.cancel(X*dy-Y*dx-2*(p*dq-q*dp)/D) == 0
c, dc = s.symbols('c dc')
assert s.expand(c*p*(dc*q+c*dq)-c*q*(dc*p+c*dp)-c*c*(p*dq-q*dp)) == 0

def height(poly):
    return max(map(abs, s.Poly(poly,t).all_coeffs()))

count = 0
examples = []
for A in [s.Integer(0), s.Integer(1), 1+t, 1-t]:
    for B in [s.Integer(-1), s.Integer(0), s.Integer(1), 1+t]:
        P = s.expand(1+t*t*A*A)
        Q = s.expand(t*P.subs(t,1)+t*(1-t)*B)
        W = s.expand(P*s.diff(Q,t)-Q*s.diff(P,t))
        den = s.expand(P*P+Q*Q)
        assert Q.subs(t,0) == 0 and Q.subs(t,1) == P.subs(t,1)
        d = max(s.degree(P,t),s.degree(Q,t),1)
        H = max(height(P),height(Q),1)
        assert W == 0 or s.degree(W,t) <= 2*d-2
        assert height(W) <= 2*d*(d+1)*H*H
        assert height(den) <= 2*(d+1)*H*H
        C = 1+t+t*t
        assert s.expand((C*P)*s.diff(C*Q,t)-(C*Q)*s.diff(C*P,t)-C*C*W) == 0
        for z in [Fraction(0),Fraction(1,7),Fraction(1,3),Fraction(1,2),Fraction(1)]:
            n,m = z.numerator,z.denominator
            value = s.Rational(n,m)
            Ph,Qh = m**d*P.subs(t,value),m**d*Q.subs(t,value)
            Wh = m**(2*d-2)*W.subs(t,value)
            assert all(x.is_Integer for x in [Ph,Qh,Wh])
            raw = s.Rational(4*m*m*Wh,Ph*Ph+Qh*Qh)
            assert raw == s.cancel(4*W/den).subs(t,value)
            assert s.denom(raw) <= 2*(d+1)**2*H*H*max(abs(n),m)**(2*d)
            count += 1
        examples.append({'P':str(P),'Q':str(Q),'degree':int(d),'height':int(H)})
P,Q=1+t*t,2*t
assert s.cancel(4*(P*s.diff(Q,t)-Q*s.diff(P,t))/(P*P+Q*Q)-8*(1-t*t)/(1+6*t*t+t**4)) == 0
report={'status':'PASS','symbolic_identities':4,'integer_constructions':len(examples),
        'exact_sample_cost_checks':count,'examples':examples,
        'scope':'Exact interface verification; no new pi approximation or asymptotic quality claim.'}
Path(__file__).with_name('checks.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print('PASS: circle identities, common-factor cancellation, 16 integer constructions, 80 exact sample costs')
