"""Rational Schur-complement projection at fixed reduced rational output."""
from fractions import Fraction as F
from pathlib import Path
from math import lcm,gcd
import json,sys,time
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).resolve().parent))
from gram_wave import data,hcoeff,sf

def main():
    results=[]; r=F(355,113)
    for m in range(1,5):
        a,b,c=m,2*m,4*m+1; M,B=data(c,a,b); d=m+1
        HB=s.Matrix([hcoeff(M,row,a,b,c) for row in B.tolist()])
        G=s.Matrix([[s.Rational(1,2*i+2*j+1) for j in range(d+1)] for i in range(d+1)])
        K=HB*G*HB.T
        P=s.Matrix([list(B[:,-1]),[sum(row[:-1]) for row in B.tolist()]])
        assert P.rank()==2 and K.det()>0
        S=(P*K.inv()*P.T).inv()
        target=s.Matrix([1,-2**(c-3)*s.Rational(r.numerator,r.denominator)])
        y=K.inv()*P.T*S*target
        assert P*y==target
        R=B.T*y
        den=lcm(*(int(z.q) for z in R))
        R=[int(z*den) for z in R]
        alpha=R[-1]; rat=F(-sum(R[:-1]),2**(c-3)*alpha)
        assert rat==r
        en=sf((target.T*S*target)[0])
        assert en>0
        factor=F(4*r.denominator**2)*F(2,27)**m
        results.append({'m':m,'p':str(r.numerator),'q':str(r.denominator),
          'alpha':str(alpha),'R':list(map(str,R[:-1])),
          'max_coefficient_bits':max(abs(z).bit_length() for z in R),
          'effective_metric_2x2':[[str(z) for z in list(S.row(i))] for i in range(2)],
          'optimal_normalized_energy':str(en),
          'gram_Q_bound_squared':str(factor*factor*en),
          'gram_Q_bound_display':float(factor)*float(en)**.5})
    out={'status':'PASS','scope':'Fixed known fraction benchmark; projection changes certificates, not the approximation or its true Q.',
      'pi_or_cf_input':False,'fixed_rational_input':'355/113','results':results}
    Path(__file__).with_name('fiber_projection_short_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'PASS','results':[{k:z[k] for k in ('m','q','max_coefficient_bits','gram_Q_bound_display')} for z in results]},indent=2))
if __name__=='__main__':main()
