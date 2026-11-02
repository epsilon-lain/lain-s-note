"""Recover the normalized cusp equation from the source's explicit curve."""
import sympy as s
from pathlib import Path
import json
x,t,z=s.symbols('x t z')
g2=3*t**3*(t-1);g3=t**4*(t-1)**2
P=4*x**3-g2*x-g3
assert s.cancel(1728*g2**3/(g2**3-27*g3**2))==1728*t
aa,bb,A,B,C=s.symbols('a b A B C')
R=A*x*x+B*x+C
connection=[]
for k in (0,1):
 rhs=x**k*(s.diff(g2,t)*x+s.diff(g3,t))/2
 residual=s.Poly(s.expand((aa+bb*x)*P+s.diff(R,x)*P-R*s.diff(P,x)/2-rhs),x)
 sol=s.solve(residual.all_coeffs(),(aa,bb,A,B,C))
 assert sol
 assert s.cancel(residual.as_expr().subs(sol))==0
 connection.append([s.factor(sol[aa]),s.factor(sol[bb])])
a,b=connection[0];c,d=connection[1]
r1=s.factor(a+d+s.diff(b,t)/b)
r0=s.factor(s.diff(a,t)+b*c-a*d-a*s.diff(b,t)/b)
ell=3/(4*t)+1/(4*(t-1))
# W(t)=F(t)/(t^(3/4)*(t-1)^(1/4)).
F1=s.factor(-(2*ell+r1))
F0=s.factor(ell**2-s.diff(ell,t)+r1*ell-r0)
# Change t=1/z; F'(t)=-z² f'(z), F''=z⁴ f''+2z³ f'.
coef2=z**4
coef1=s.factor(2*z**3-F1.subs(t,1/z)*z**2)
coef0=s.factor(F0.subs(t,1/z))
normal=s.factor(z*(1-z)/coef2)
assert s.cancel(coef1*normal-(1-s.Rational(3,2)*z))==0
assert s.cancel(coef0*normal+s.Rational(5,144))==0
out={'status':'PASS','curve':'y^2=4*x^3-3*t^3*(t-1)*x-t^4*(t-1)^2',
 'j':'1728*t','gauss_manin_matrix':[[str(q) for q in row] for row in connection],
 'gauge':'F(t)=t^(3/4)*(t-1)^(1/4)*period(dx/y), up to a nonzero constant',
 'cusp_variable':'z=1/t',
 'normalized_equation':"z*(1-z)*f''+(1-3*z/2)*f'-5*f/144=0",
 'checks':'Two exact polynomial differential reductions, j invariant, algebraic gauge and variable change.'}
Path(__file__).with_name('gauss_manin_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
