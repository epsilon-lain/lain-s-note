"""Exact finite demonstrations; these do not validate any irrationality exponent.
Run with Python 3 and SymPy: python check_classical_models.py
"""
import json
import sympy as s

x, z, epsilon = s.symbols('x z epsilon')
checks = {}

# First-type Hermite-Pade: one exponential-polynomial combination has a high zero.
a,b,c,d = s.symbols('a b c d')
F=(z-2)*s.exp(z)+z+2
coeffs=[s.expand(s.series(F,z,0,8).removeO()).coeff(z,k) for k in range(8)]
assert coeffs[:3]==[0,0,0] and coeffs[3]==s.Rational(1,6)
checks['type_I_zero_order']=3
checks['type_I_coefficients_z3_to_z7']=[str(c) for c in coeffs[3:]]

# A residue construction in the m=n=1 case, computed independently.
Qz=z**2*(z+1)**2
R=[]
for h in [0,1]:
    integrand=z**h*x*s.exp(z*s.log(x))/Qz
    R.append(s.simplify(sum(s.residue(integrand,z,p) for p in [0,-1])))
assert s.simplify(R[0]-((x+1)*s.log(x)+2*(1-x)))==0
assert s.simplify(R[1]-(x-1-s.log(x)))==0
A=s.Matrix([[2*(1-x),x+1],[x-1,-1]])
assert s.factor(A.det())==-(x-1)**2
assert s.simplify(A.det().subs(x,s.I))==2*s.I
checks['residue_remainders']=[str(r) for r in R]
checks['residue_coefficient_determinant']=str(s.factor(A.det()))
checks['determinant_at_i']=str(s.simplify(A.det().subs(x,s.I)))

# Second-type simultaneous approximation: one common Q for two logarithms.
q1,q2=s.symbols('q1 q2')
common=1+q1*x+q2*x**2
conditions=[s.expand(s.series(common*s.log(1+sign*x),x,0,4).removeO()).coeff(x,3)
            for sign in [1,-1]]
solution=s.solve(conditions,[q1,q2])
assert solution=={q1:0,q2:-s.Rational(1,3)}
Q=6-2*x**2
remainders=[]
for sign in [1,-1]:
    P=sign*6*x-3*x**2
    r=s.expand(s.series(Q*s.log(1+sign*x)-P,x,0,8).removeO())
    assert all(r.coeff(x,k)==0 for k in range(4))
    assert r.coeff(x,4)==-s.Rational(1,2)
    assert r.coeff(x,5)==sign*s.Rational(8,15)
    remainders.append(str(r))
assert Q.subs(x,s.I)==8
assert (6*x-3*x**2-(-6*x-3*x**2)).subs(x,s.I)==12*s.I
checks['type_II_normalized_denominator_solution']={str(k):str(v) for k,v in solution.items()}
checks['type_II_remainders_through_x7']=remainders

# Determinant cancellation: three distinct Taylor degrees are required.
M=s.Matrix([[1,j*epsilon,s.exp(j*epsilon)] for j in range(3)])
D=s.factor(M.det())
assert s.simplify(D-epsilon*(s.exp(epsilon)-1)**2)==0
assert s.expand(s.series(D,epsilon,0,4).removeO())==epsilon**3
checks['collision_determinant']=str(D)
checks['collision_zero_order']=3

# The arithmetic penalty increases the exponent when raw rates are fixed.
alpha,beta,delta=s.symbols('alpha beta delta',positive=True)
ratio=1+(beta+delta)/(alpha-delta)
assert s.simplify(s.diff(ratio,delta)-(alpha+beta)/(alpha-delta)**2)==0
checks['denominator_penalty_derivative']=str(s.simplify(s.diff(ratio,delta)))
checks['scope']='Exact finite identities only; no proof or computational validation of mu(pi)=2.'
print(json.dumps(checks,indent=2))
