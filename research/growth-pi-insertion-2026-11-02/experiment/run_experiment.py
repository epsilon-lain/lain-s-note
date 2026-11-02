"""Fixed rational feedback grid; float proposals, exact positive certificates and pi scoring.
This bounded discovery routine can miss cycles. Uncertified parameters remain unresolved.
"""
from fractions import Fraction as Q
from math import floor
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def point(v):return Q(v),Q(v)
def atan_reciprocal(q,terms):
 s=sum((Q((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(terms)),Q(0))
 r=Q(1,(2*terms+1)*q**(2*terms+1))
 return (s,s+r) if terms%2==0 else (s-r,s)
A=atan_reciprocal(5,50);B=atan_reciprocal(239,16)
PI=(16*A[0]-4*B[1],16*A[1]-4*B[0])
assert Q(314,100)<PI[0]<PI[1]<Q(22,7)

def quality(r):
 d=r.denominator
 if r<PI[0]:return ((PI[0]-r)*d*d,(PI[1]-r)*d*d)
 if r>PI[1]:return ((r-PI[1])*d*d,(r-PI[0])*d*d)
 return (Q(0),max(r-PI[0],PI[1]-r)*d*d)

def FF(x,om,kap):
 u=x-floor(x)
 psi=u if u<=.25 else (.5-u if u<=.75 else u-1)
 return x+om-kap*psi

def branch(x,om,kap):
 j=x.numerator//x.denominator;u=x-j
 if u<=Q(1,4):return 1-kap,om+kap*j,(Q(j)-Q(1,4),Q(j)+Q(1,4))
 if u<=Q(3,4):return 1+kap,om-kap*(j+Q(1,2)),(Q(j)+Q(1,4),Q(j)+Q(3,4))
 return 1-kap,om+kap*(j+1),(Q(j)+Q(3,4),Q(j)+Q(5,4))

def bitheight(v):return max(abs(v.numerator).bit_length(),v.denominator.bit_length())

def certify(om,kap,xfloat,p,q):
 # Float only chooses a finite branch word, which is then validated exactly.
 A=Q(1);B=Q(0);xf=xfloat;oms=float(om);ks=float(kap)
 for _ in range(q):
  j=floor(xf);u=xf-j
  if u<=.25:a,b=1-kap,om+kap*j
  elif u<=.75:a,b=1+kap,om-kap*(j+Q(1,2))
  else:a,b=1-kap,om+kap*(j+1)
  A,B=a*A,a*B+b;xf=FF(xf,oms,ks)
 if not 0<A<1:return None
 x=(Q(p)-B)/(A-1)
 if not Q(-1,100)<x<Q(101,100):return None
 y=x;r=Q(1,100);Ai=Q(1);maxbits=bitheight(y)
 for _ in range(q):
  a,b,(lo,hi)=branch(y,om,kap)
  if not lo<y<hi:return None
  r=min(r,(y-lo)/(2*Ai),(hi-y)/(2*Ai))
  y=a*y+b;Ai*=a;maxbits=max(maxbits,bitheight(y))
 if y!=x+p or Ai!=A:return None
 # Every interval stays within its branch; F^q(I)-p=x+A(I-x) strictly inside I.
 return {'omega':str(om),'kappa':str(kap),'lift_p':p,'iterate_q':q,
  'rotation':str(Q(p,q)),'periodic_point':str(x),'trapping_radius':str(r),
  'q_step_multiplier':str(A),'strict_margin':str((1-A)*r),
  'max_orbit_rational_bits':maxbits,
  'cycle_equation_residual':'0','positive_certificate':'PASS'}

def discover(om,kap,burn=6000,maxq=256):
 x=.123456789;oms=float(om);ks=float(kap)
 for _ in range(burn):x=FF(x,oms,ks)
 x0=x-floor(x);y=x0
 for q in range(1,maxq+1):
  y=FF(y,oms,ks);p=round(y-x0)
  if abs(y-x0-p)<.0001:
   c=certify(om,kap,x0,p,q)
   if c:return c
 return None

def run():
 om=Q(1,7);input_r=3+om;input_lo,input_hi=quality(input_r)
 summaries=[];certificates=[];unknown=[]
 for j in range(1,601):
  kap=Q(j,5000);c=discover(om,kap)
  if not c:unknown.append(j);continue
  out=3+Q(c['rotation']);lo,hi=quality(out)
  c.update({'output':str(out),'actual_reduced_denominator':out.denominator,
            'quality_interval':[str(lo),str(hi)]})
  certificates.append(c)
  summaries.append({'grid_index':j,'kappa':str(kap),'output':str(out),
                    'denominator':out.denominator,'quality_display':float((lo+hi)/2),
                    'raw_cycle_length':c['iterate_q']})
 best=min(certificates,key=lambda c:(Q(c['quality_interval'][1]),c['actual_reduced_denominator']))
 assert all(Q(best['quality_interval'][1]) < Q(c['quality_interval'][0]) for c in certificates if c['output'] != best['output'])
 special=next(c for c in certificates if c['kappa']=='493/5000')
 assert special['output']=='355/113'
 qlo,qhi=map(Q,special['quality_interval'])
 assert qhi<input_lo/18 and special['actual_reduced_denominator']>input_r.denominator
 out={'status':'PASS','scope':'All 600 declared rational parameters attempted; only positive exact certificates accepted. Cycle discovery is bounded and incomplete. Best means among certified outputs, not a global optimum.',
      'parameter_selection':'This full fixed-grid reproduction was specified after target-guided preliminary exploration near kappa=.0986. Its execution uses no pi digits or continued fractions for branch proposals or cycle certificates; pi is evaluated only after exact certificates for scoring. This is not an independent blind parameter discovery.',
      'grid':{'omega':'1/7','kappa':'j/5000, j=1,...,600','burn_in':6000,'max_tested_period':256,'near_return_tolerance':'1e-4 (proposal only)'},
      'attempted_parameters':600,'positive_certificates':len(certificates),'unresolved_parameters':len(unknown),
      'unresolved_grid_indices':unknown,'certified_rows':summaries,
      'input':{'rational':str(input_r),'quality_interval':[str(input_lo),str(input_hi)],'quality_display':float((input_lo+input_hi)/2)},
      'best_certified':best,'explicit_joint_improvement':special,
      'joint_improvement':{'denominator_ratio':str(Q(113,7)),
                          'quality_gain_certified_lower':str(input_lo/qhi),
                          'quality_gain_lower_display':float(input_lo/qhi),
                          'quality_gain_greater_than_18':True},
      'pi_interval':[str(x) for x in PI],
      'source_known_benchmark':'355/113 is a classical known approximation, not a new pi record or infinite quality-zero construction.'}
 (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
 (HERE/'certificates.json').write_text(json.dumps(certificates,indent=2)+'\n')
 print(json.dumps({k:out[k] for k in ['status','attempted_parameters','positive_certificates','unresolved_parameters','joint_improvement']},indent=2))
 print(json.dumps({'best_output':best['output'],'best_kappa':best['kappa'],'special_max_orbit_bits':special['max_orbit_rational_bits']},indent=2))

if __name__=='__main__':run()
