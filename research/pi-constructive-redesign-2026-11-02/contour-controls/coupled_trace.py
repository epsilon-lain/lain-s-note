"""Genuinely coupled conjugate kernels: independent small/large endpoint powers."""
from quadratic_trace_contour import K,one,zero,v,mul,power,from_numerator
from endpoint_design import pibounds
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as sy

def coupled(a,b,d,c):
 big=K.from_sympy(-sy.sqrt(2)-1)
 Q1=[one+v*v,-2*one,one];Q2=[one+big*big,-2*one,one]
 N=mul(mul(power(Q1,a),power(Q2,b)),power([-one,one],2*d))
 return from_numerator(N,c)

def main():
 pl,ph=pibounds(300)
 rows=[];records=[];best=None
 for cost in range(1,19):
  for a in range(5):
   for b in range(5):
    if a+b==0:continue
    for d in range(3):
     c=cost-a-b-d
     if not 1<=c<=8:continue
     value=coupled(a,b,d,c)
     if value is None:continue
     r,r0,B=value
     if r<=0:continue
     p,q=r.numerator,r.denominator
     if r<pl:lo,hi=(pl-r)*q*q,(ph-r)*q*q
     elif r>ph:lo,hi=(r-ph)*q*q,(r-pl)*q*q
     else:raise AssertionError('pi interval too wide')
     scale=10**18
     row={'a':a,'b':b,'d':d,'c':c,'p':str(p),'q':str(q),'rho':str(K.to_sympy(r0)),'B':str(K.to_sympy(B)),
          'quality_display':float((lo+hi)/2),'quality_strict_outward_bounds':[str(F(lo.numerator*scale//lo.denominator,scale)),str(F(-((-hi.numerator*scale)//hi.denominator),scale))]}
     rows.append(row)
     if best is None or hi<best:best=hi;records.append(row)
 result={'status':'PASS','candidate_count':len(rows),'scope':'Finite coupled algebraic-contour experiment; no target rational or continued fraction input, no new pi record or infinite quality-zero theorem.', 'bounds':{'a_b':[0,4],'a_plus_b_positive':True,'d':[0,2],'c':[1,8]}, 'record_improvements':records,
         'best_10':sorted(rows,key=lambda r:r['quality_display'])[:10],
         'best_with_both_factors':min((r for r in rows if r['a']>0 and r['b']>0),key=lambda r:r['quality_display']),
         'best_by_pole_exponent':{str(c):min((r for r in rows if r['c']==c),key=lambda r:r['quality_display']) for c in range(1,9)}}
 Path(__file__).with_name('coupled_trace_checks.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'status':'PASS','count':len(rows),'records':records,'best_by_c':{c:{k:v for k,v in row.items() if k in ['a','b','d','c','p','q','quality_display']} for c,row in result['best_by_pole_exponent'].items()}},indent=2))
if __name__=='__main__':main()
