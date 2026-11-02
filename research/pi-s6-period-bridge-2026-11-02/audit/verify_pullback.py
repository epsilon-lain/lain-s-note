"""Exact rational-operator audit of the elliptic-period bridge."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb
import json
import sympy as s
u=s.symbols('u');h=s.Function('h')(u)
z=27*u**2/(4-u)**3
A=s.factor(z/s.diff(z,u))
ell=-1/(2*(4-u))
seq=[h]
for _ in range(3):
 seq.append(s.expand(A*(s.diff(seq[-1],u)+ell*seq[-1])))
pull=(1-z)*seq[3]-3*z*seq[2]/2-23*z*seq[1]/36-5*z*seq[0]/72
q=[h]
for _ in range(3):q.append(s.expand(u*s.diff(q[-1],u)))
target=(1-u)*q[3]-3*u*q[2]/2-3*u*q[1]/4-u*q[0]/8
factor=s.factor(s.expand(pull).coeff(s.diff(h,u,3))/s.expand(target).coeff(s.diff(h,u,3)))
checks=[]
for j in range(4):
 term=s.diff(h,u,j) if j else h
 assert s.cancel(s.expand(pull-factor*target).coeff(term))==0
 checks.append(j)
# Verify the ordinary-point normalizations independently.
H=[s.Rational(comb(2*n,n)**3,64**n) for n in range(20)]
G=[s.Rational(factorial(6*n),factorial(3*n)*factorial(n)**3*1728**n) for n in range(10)]
left=s.series(sum(G[n]*z**n for n in range(10)),u,0,18).removeO()
right=s.series(s.sqrt(1-u/4)*sum(H[n]*u**n for n in range(20)),u,0,18).removeO()
assert s.expand(left-right)==0
assert z.subs(u,s.Rational(1,64))==s.Rational(1728,255**3)
out={'status':'PASS','coefficientwise_differential_operator_checks':checks,
     'pullback':'z=27*u^2/(4-u)^3','gauge':'G(z(u))=sqrt(1-u/4)*H(u)',
     'rational_operator_factor':str(factor),'series_agreement_degrees':[0,17],
     'evaluation':'u=1/64 gives j=255^3',
     'scope':'An exact elliptic-period connection; not a rational-approximation quality theorem.'}
Path(__file__).with_name('pullback_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
