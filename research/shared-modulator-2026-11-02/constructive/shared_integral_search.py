"""Shared modulation constraints inserted into the exact pi endpoint kernel.
All candidate generation and Gram ranking use rational/integer arithmetic.
"""
from fractions import Fraction as F
from math import comb,lcm,gcd,isqrt
from pathlib import Path
import json,sys,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'support'))
from integer_kernel import data
from integral_certificate import certificate
T=s.symbols('t');S=T*T-10*T+5

def sf(q):
 q=s.Rational(q);return F(int(q.p),int(q.q))

def Mmatrix(c):
 M=s.zeros(c+1,c+1)
 for j in range(c):M[j,j]=2*j+1;M[j+1,j]=2*j-2*c+3
 for j in range(c):M[j,c]=comb(c-1,j)
 return M

def hmatrix(M,B,a,b,c):
 rows=[]
 for row in B.tolist():
  N=M*s.Matrix(row);p=sum(v*T**j for j,v in enumerate(N))
  h,r=s.div(p,T**a*(1-T)**b);assert r==0
  rows.append([int(s.expand(h).coeff(T,j)) for j in range(c-a-b+1)])
 return s.Matrix(rows)

def saturated_core(M,B,c):
 # Additional conditions N mod S5=0 are equivalent to H mod S5=0,
 # since the endpoint factor is coprime with S5.
 rows=[]
 for d in range(2):
  rows.append([s.rem(sum(M[j,k]*T**j for j in range(c+1)),S,T).coeff(T,d) for k in range(c+1)])
 C=s.Matrix(rows)*B.T
 normal,left,right=smith_normal_decomp(DomainMatrix.from_Matrix(C).convert_to(s.ZZ))
 rank=C.rank();assert rank==2
 U=right.to_Matrix();assert abs(U.det())==1
 K=U[:,rank:].T*B
 assert C*U[:,rank:]==s.zeros(2,K.rows)
 return K.lll(delta=F(3,4))

def mapped_basis(M,B,c,lam):
 cn=c+6;Mn=Mmatrix(cn);inv=Mn.inv()
 multiplier=s.Poly((1+T)**6-4*s.Rational(lam.numerator,lam.denominator)*T*S*S,T)
 rows=[]
 for row in B.tolist():
  N=M*s.Matrix(row);p=sum(v*T**j for j,v in enumerate(N))
  out=s.Poly(p*multiplier.as_expr(),T)
  V=inv*s.Matrix([out.nth(j) for j in range(cn+1)])
  D=lcm(*(int(v.q) for v in V));ints=[int(v*D) for v in V]
  content=gcd(*ints);rows.append([v//content for v in ints])
 Bn=s.Matrix(rows).lll(delta=F(3,4))
 return Mn,Bn,cn

def root_upper(v,bits=80):
 v=sf(v);n=isqrt((v.numerator<<(2*bits))//v.denominator)
 return F(n+1,1<<bits)

def select(M,B,a,b,c,m):
 HB=hmatrix(M,B,a,b,c);d=HB.cols-1
 G=s.Matrix([[s.Rational(1,2*i+2*j+1) for j in range(d+1)] for i in range(d+1)])
 L,D=G.LDLdecomposition(hermitian=False)
 roots=[root_upper(D[j,j]) for j in range(d+1)]
 metric=L*s.diag(*[s.Rational(v.numerator,v.denominator) for v in roots])
 FE=HB*metric*s.Rational(1,16**m)
 seen=set();candidates=[]
 for budget in [0,8*m,12*m,16*m,20*m]:
  cost=s.Rational(1,2**budget)
  feature=FE.row_join(B[:,-1:]*cost).row_join(s.Matrix([sum(row[:-1]) for row in B.tolist()])*cost)
  D0=lcm(*(int(v.q) for v in feature));extended=(feature*D0).row_join(B)
  reduced=extended.lll(delta=F(3,4))
  for row in reduced.tolist():
   R=list(map(int,row[-(c+1):]));alpha=R[-1]
   if not alpha:continue
   rat=F(-sum(R[:-1]),2**(c-3)*alpha)
   if rat<=0 or rat in seen:continue
   seen.add(rat);H=hmatrix(M,s.Matrix([R]),a,b,c)[0,:].T
   energy=sf((H.T*G*H)[0]);factor=F(4*rat.denominator**2,abs(alpha)*16**m)
   candidates.append({'p':str(rat.numerator),'q':str(rat.denominator),'alpha':str(alpha),'R':R[:-1],
        'H':[int(v) for v in H],'energy':str(energy),'gram_Q_bound_squared':str(factor*factor*energy),
        'gram_Q_bound_display':float(factor)*float(energy)**.5,'selection_budget':budget,
        'max_construct_coefficient_bits':max(abs(v).bit_length() for v in R)})
 candidates.sort(key=lambda r:F(r['gram_Q_bound_squared']))
 if candidates:
  for r in candidates[:2]:
   cert=certificate(c,r['R'],int(r['alpha']),K=180)
   assert F(cert['quality_upper'])**2<F(r['gram_Q_bound_squared'])
   r['quality_interval']=[cert['quality_lower'],cert['quality_upper']]
   r['quality_display']=cert['quality_upper_display']
 period_coordinates=[[int(row[-1]),int(sum(row[:-1]))] for row in B.tolist()]
 period_map=s.Matrix(period_coordinates).T
 return {'basis_rank':B.rows,'period_map_rank':period_map.rank(),'basis_period_coordinates':period_coordinates,'c':c,'candidate_count':len(candidates),'best':candidates[:2]}

def run():
 start=time.monotonic();rows=[]
 for m in range(1,5):
  a,b,c=m,2*m,4*m+1;M,B=data(c,a,b)
  for name,MM,BB,cc,lam in [('unconstrained',M,B,c,None),('core_S5',M,saturated_core(M,B,c),c,None)]:
   item={'m':m,'method':name,**select(MM,BB,a,b,cc,m)}
   rows.append(item);print(m,name,[(r['p'],r['q'],r['quality_display']) for r in item['best']],flush=True)
  # The literal shared-g product requires six extra denominator powers.
  for lam in [F(1,4),F(1,2),F(1),F(2)]:
   Mn,Bn,cn=mapped_basis(M,B,c,lam)
   item={'m':m,'method':'literal_shared_g','lambda':str(lam),**select(Mn,Bn,a,b,cn,m)}
   rows.append(item);print(m,'literal',str(lam),[(r['p'],r['q'],r['quality_display']) for r in item['best']],flush=True)
 out={'status':'PASS','pi_or_CF_input_in_generation':False,
      'scope':'Finite fixed-Gram-LLL comparison; constrained spaces are subspaces, not a global improvement theorem.',
      'core_factor':'t^2-10*t+5','literal_multiplier':'(1+t)^6-4*lambda*t*(t^2-10*t+5)^2',
      'shape_budget':'16^(-m); the six extra literal-g denominator powers are retained and bounded by 1 for the conservative Gram metric.',
      'rows':rows,'elapsed_seconds_display':time.monotonic()-start}
 (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':'PASS','cases':len(rows),'seconds':out['elapsed_seconds_display']},indent=2))

if __name__=='__main__':run()
