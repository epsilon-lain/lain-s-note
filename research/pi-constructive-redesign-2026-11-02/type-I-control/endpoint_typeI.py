from fractions import Fraction as F
from math import gcd,lcm
from functools import reduce
from pathlib import Path
import sympy as s
import json,time

def coefficient(a,k):return s.Rational((-1)**k,(2*k+1)*a**(2*k+1)) if k>=0 else s.Rational(0)
def atan_interval(a,K=150):
 partial=sum((F((-1)**k,(2*k+1)*a**(2*k+1)) for k in range(K)),F(0))
 return partial,partial+F(1,(2*K+1)*a**(2*K+1))
a,b=atan_interval(5);c,d=atan_interval(239);pilo,pihi=16*a-4*d,16*b-4*c

def solve(m0,m2,m3):
 degrees=(m0,m2,m3);offsets=(0,m0+1,m0+m2+2);M=sum(degrees)+3;L=M-2
 rows=[]
 for k in range(L):
  row=[s.Rational(0)]*M
  if k<=m0:row[k]=1
  for a,m,o in ((2,m2,offsets[1]),(3,m3,offsets[2])):
   for j in range(m+1):row[o+j]=coefficient(a,k-j)
  rows.append(row)
 rows.append([0]*(m0+1)+[1]*(m2+1)+[-1]*(m3+1))
 ns=s.polys.matrices.DomainMatrix.from_Matrix(s.Matrix(rows)).nullspace().to_Matrix()
 if ns.rows!=1:return {'degrees':degrees,'nullity':ns.rows}
 vec=list(ns[0,:]);den=lcm(*(int(v.q) for v in vec));vec=[int(v*den) for v in vec];g=reduce(gcd,map(abs,vec));vec=[v//g for v in vec]
 polys=[vec[o:o+m+1] for o,m in zip(offsets,degrees)]
 A=sum(polys[0]);Q=sum(polys[1]);assert Q==sum(polys[2])
 assert s.Matrix(rows)*s.Matrix(vec)==s.zeros(len(rows),1)
 if not Q:return {'degrees':degrees,'nullity':1,'Q':0}
 r=F(-4*A,Q);p,q=r.numerator,r.denominator;g=gcd(abs(Q),abs(4*A))
 if r<pilo:lo,hi=(pilo-r)*q*q,(pihi-r)*q*q
 elif r>pihi:lo,hi=(r-pihi)*q*q,(r-pilo)*q*q
 else:raise AssertionError('pi bound too wide')
 bound=sum((F(abs(poly[j]),(2*(L-j)+1)*a**(2*(L-j)+1)) for a,poly in ((2,polys[1]),(3,polys[2])) for j in range(len(poly))),F(0))
 upper=F(4*abs(Q),g*g)*bound
 assert lo<=upper
 return {'degrees':degrees,'order':L,'nullity':1,'coefficients':polys,'A':A,'Q':Q,'reduction_gcd':g,'p':p,'q':q,'quality_interval':[str(lo),str(hi)],'quality_display':float((lo+hi)/2),'simple_integral_upper':str(upper),'simple_integral_upper_display':float(upper)}

def main():
 started=time.monotonic();results=[]
 for total in range(0,16):
  for m0 in range(total+1):
   for m2 in range(total-m0+1):
    m3=total-m0-m2
    if max(m0,m2,m3)>7:continue
    results.append(solve(m0,m2,m3))
 valid=[r for r in results if r.get('Q')]
 valid.sort(key=lambda r:r['quality_display'])
 out={'status':'PASS','candidate_count':len(results),'valid_count':len(valid),'scope':'Finite endpoint-coupled type-I experiment, no asymptotic quality theorem','elapsed_seconds':time.monotonic()-started,'best_10':valid[:10],'balanced':[r for r in results if len(set(r['degrees']))==1]}
 Path(__file__).with_name('endpoint_typeI_results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({**out,'best_10':[{k:v for k,v in row.items() if k not in ['quality_interval','coefficients','simple_integral_upper']} for row in valid[:10]],'balanced':[{k:v for k,v in row.items() if k not in ['quality_interval','coefficients','simple_integral_upper']} for row in out['balanced']]},indent=2))
if __name__=='__main__':main()
