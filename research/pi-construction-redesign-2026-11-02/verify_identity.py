from fractions import Fraction as F
from pathlib import Path
import json
import sympy as s
x=s.symbols('x')
k=x**4*(1-x*x)**3*(10*x*x-3)/(1+x*x)**4
assert s.simplify(s.integrate(k,(x,0,1))-(113*s.pi-355))==0
row=json.loads(Path(__file__).with_name('combination_results.json').read_text())['best']
assert (row['p'],row['q'])==(355,113)
l,h=map(F,row['quality_interval'])
assert F('0.00340631193013')<l<h<F('0.00340631193015')
result={'status':'PASS','independent_symbolic_integral':'113*pi-355','quality_strict_bounds':['0.00340631193013','0.00340631193015'],'scope':'Known rational recovered by bounded constructive search; no infinite-subsequence theorem.'}
Path(__file__).with_name('identity_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
