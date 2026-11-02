"""Validate the predeclared growing family directly from integer integrands."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json
from integral_certificate import certificate,numerator

def outward(value,up=False):
 scale=10**24
 n=-((-value.numerator*scale)//value.denominator) if up else value.numerator*scale//value.denominator
 return str(F(n,scale))

def run(source):
 data=json.loads(source.read_text());rows=[]
 for case in data['per_m']:
  m=case['m'];row=case['best_observed'];c=case['c']
  assert (case['a'],case['b'],c)==(m,2*m,4*m+1)
  alpha=int(row['alpha']);r=list(map(int,row['R_coefficients']))
  assert numerator(c,r,alpha)==list(map(int,row['N_coefficients']))
  result=certificate(c,r,alpha,max(100,5*c+40))
  assert (result['p'],result['q'])==(row['p'],row['q'])
  lower,upper=F(result['quality_lower']),F(result['quality_upper'])
  checklo,checkhi=map(F,row['quality_rational_bounds'])
  assert max(lower,checklo)<=min(upper,checkhi)
  rows.append({'m':m,'c':c,'p':result['p'],'q':result['q'],'quality_strict_outward_bounds':[outward(lower),outward(upper,True)],'quality_upper_display':float(upper),'certified_below_one':upper<1,'K':result['K']})
 out={'status':'PASS','family':'a=m,b=2m,c=4m+1','completed_m':data['completed_m'],'pi_input_in_integral_certificates':False,'selection_scope':'Best among five predeclared cost weights per m; not a globally optimal integer vector.','certificates':rows,'scope':'Finite constructive family validation, known rationals; no quality-zero subsequence proof.'}
 Path(__file__).with_name('new_family_integral_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,default=Path(__file__).parent.parent/'integer-kernel/new_family_results.json');run(parser.parse_args().input)
