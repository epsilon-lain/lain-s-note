"""Independent checks against the exact primitive and Machin intervals."""
from fractions import Fraction as F
from random import Random
from pathlib import Path
import json
from integral_certificate import numerator,integral_interval

def atan(a,K=150):
 s=sum((F((-1)**k,(2*k+1)*a**(2*k+1)) for k in range(K)),F(0))
 return s,s+F(1,(2*K+1)*a**(2*K+1))
a,b=atan(5);u,v=atan(239);pl,ph=16*a-4*v,16*b-4*u
rng=Random(20261102);count=0
for c in range(1,21):
 for sample in range(4):
  r=[rng.randrange(-10000,10001) for _ in range(c)]
  alpha=0 if sample==0 else rng.randrange(-1000,1001)
  lo,hi=integral_interval(numerator(c,r,alpha),c,100)
  base=F(sum(r),2**(c-1))
  exactlo=base+F(alpha,4)*(pl if alpha>=0 else ph)
  exacthi=base+F(alpha,4)*(ph if alpha>=0 else pl)
  assert lo<=exactlo<=exacthi<=hi
  count+=1
out={'status':'PASS','independent_primitive_checks':count,'c_range':[1,20],'zero_period_cases':20,'seed':20261102,'main_certificate_uses_pi':False,'independent_check_uses_machin_interval':True}
Path(__file__).with_name('independent_certificate_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
