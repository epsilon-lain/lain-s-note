"""Exact integral-local-system extraction from the S6 period matrices.

No numerical periods or value of pi is used by this verifier.
"""
from pathlib import Path
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form

I=s.eye(4)
T1=s.Matrix([[1,0,-6,2],[0,-1,1,1],[0,-1,0,1],[0,0,0,1]])
T2=s.Matrix([[1,6,0,-3],[0,0,-1,1],[0,1,0,0],[0,0,0,1]])
assert T1**3==I and T2**4==I
assert T1.det()==T2.det()==1
T0=(T1*T2).inv();N=T0-I
assert N*N==s.zeros(4)
assert N.rank()==2
SN=smith_normal_form(N,domain=s.ZZ)
assert [abs(SN[j,j]) for j in range(4)]==[1,1,0,0]
# Primitive invariant flag: Z gamma < span(gamma,u,w) < Z^4.
for T in (T1,T2):
 assert T[:,0]==s.Matrix([1,0,0,0])
 assert all(T[3,j]==0 for j in range(3))
 assert T[3,3]==1
A1=T1[1:3,1:3];A2=T2[1:3,1:3];A0=T0[1:3,1:3]
assert A1==s.Matrix([[-1,1],[-1,0]])
assert A2==s.Matrix([[0,-1],[1,0]])
assert A0==s.Matrix([[1,-1],[0,1]])
assert A1**3==s.eye(2) and A2**2==-s.eye(2)
assert A2**4==s.eye(2)
for m in range(-20,21):
 assert T0**m==I+m*N
# The rank-two quotient does not split off as an invariant rational subspace.
z=s.symbols('z0:8');X=s.Matrix(4,2,z)
relations=list(T1*X-X*A1)+list(T2*X-X*A2)
assert s.linsolve(relations,z)==s.FiniteSet(tuple([0]*8))
# No nonconstant entire period column is independent of both extension variables.
mu,beta,tau=s.symbols('mu beta tau')
P=s.Matrix([[6*mu,tau,1,0],[beta,mu,0,1]])
v=s.symbols('v0:4');period=P*s.Matrix(v)
coeff_relations=[]
for entry in period:
 coeff_relations.extend([s.diff(entry,mu),s.diff(entry,beta)])
assert s.linsolve(coeff_relations,v)==s.FiniteSet((0,0,v[2],v[3]))
result={
 'status':'PASS',
 'source':'https://alpo.ge/s6.pdf, introductory setup on PDF page 2 (numbered page 1)',
 'T1':T1.tolist(),'T2':T2.tolist(),'T0':T0.tolist(),'N':N.tolist(),
 'N_squared_zero':True,'N_rank':2,'N_smith_invariants':[1,1,0,0],
 'primitive_invariant_flag':'Z gamma < Z gamma + Z u + Z w < Z^4',
 'rank_two_quotient_matrices':{'A1':A1.tolist(),'A2':A2.tolist(),'A0':A0.tolist()},
 'rank_two_quotient_maps':['tau -> (tau-1)/tau','tau -> -1/tau','tau -> tau-1'],
 'cusp_powers_checked':[-20,20],
 'rational_equivariant_embedding_dimension':0,
 'whole_period_column_mu_beta_free':'v_gamma=v_u=0; only constant columns remain',
 'scope':'An exact integral rank-two subquotient and cusp-difference mechanism, not a rational approximation to pi or a claimed elliptic quotient of the general two-torus.'
}
Path(__file__).with_name('monodromy_checks.json').write_text(json.dumps(result,indent=2,default=int)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('T1','T2','T0','N')},indent=2,default=int))
