from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import sympy as s
r=json.loads(Path(__file__).with_name('antiderivative_search_results.json').read_text());t=s.symbols('t');out=[]
for row in r['records']:
 c,a,b=row['c'],row['a'],row['b'];v=list(map(int,row['R_coefficients']));alpha=int(row['alpha']);M=s.zeros(c+1,c+1)
 for j in range(c):M[j,j]=2*j+1;M[j+1,j]=2*j-2*c+3
 for j in range(c):M[j,c]=comb(c-1,j)
 N=M*s.Matrix(v+[alpha]);P=sum(int(z)*t**j for j,z in enumerate(N));H,rem=s.div(P,t**a*(1-t)**b);assert rem==0
 h=sum(abs(int(x)) for x in s.Poly(H,t).all_coeffs());q=int(row['q']);base=F(a**a*b**b,2**a*(a+b)**(a+b));upper=F(4*q*q*h,abs(alpha))*base
 lo,hi=map(F,row['quality_rational_bounds']);assert hi<upper
 out.append({'c':c,'a':a,'b':b,'p':row['p'],'q':row['q'],'H_coefficients_low_to_high':[str(s.expand(H).coeff(t,j)) for j in range(c-a-b+1)],'H_l1':str(h),'Q_actual_display':row['quality_display'],'Q_upper_certificate':str(upper),'Q_upper_display':float(upper)})
Path(__file__).with_name('constructive_upper_certificates.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS:',len(out),'exact constructive upper certificates')
