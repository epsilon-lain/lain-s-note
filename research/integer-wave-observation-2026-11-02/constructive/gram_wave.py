"""Exact rational continuum-wave metric; no pi input in selection/certification."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,lcm,isqrt
import json,sys,time
import sympy as s
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'support'))
from integer_kernel import data
from integral_certificate import certificate

def sf(v): return F(int(v.p),int(v.q))

def root_upper(v,bits=80):
    v=sf(v)
    n=isqrt((v.numerator << (2*bits))//v.denominator)
    return F(n+1,1<<bits)

def hcoeff(M,row,a,b,c):
    t=s.Symbol('t'); N=M*s.Matrix(row)
    p=sum(v*t**j for j,v in enumerate(N))
    h,rem=s.div(p,t**a*(1-t)**b)
    assert rem==0
    return [int(s.expand(h).coeff(t,j)) for j in range(c-a-b+1)]

def quadratic(h,G):
    return sf((s.Matrix([h])*G*s.Matrix(h))[0])

def main():
    start=time.monotonic(); out=[]; existing=[]
    old=json.loads((ROOT/'support/existing_family_results.json').read_text())
    for entry in old['per_m']:
        m=entry['m']; r=entry['best_observed']; h=list(map(int,r['H_coefficients']))
        G=s.Matrix([[s.Rational(1,2*i+2*j+1) for j in range(len(h))] for i in range(len(h))])
        v=quadratic(h,G); l1=sum(abs(x) for x in h)
        # Compare squared bounds, so all comparisons remain rational.
        improve=F(l1*l1)/v
        factor=F(4*int(r['q'])**2,abs(int(r['alpha'])))*F(2,27)**m
        existing.append({'m':m,'p':r['p'],'q':r['q'],'H':h,'G_energy':str(v),
          'old_shape_upper':r['shape_quality_upper'],
          'gram_shape_upper_display':float(factor)*float(v)**.5,
          'guaranteed_improvement_squared':str(improve),
          'improvement_display':float(improve)**.5})
    for m in range(1,7):
        a,b,c=m,2*m,4*m+1; M,B=data(c,a,b); d=c-a-b
        HB=s.Matrix([hcoeff(M,row,a,b,c) for row in B.tolist()])
        G=s.Matrix([[s.Rational(1,2*i+2*j+1) for j in range(d+1)] for i in range(d+1)])
        L,D=G.LDLdecomposition(hermitian=False)
        roots=[root_upper(D[i,i]) for i in range(d+1)]
        metric=L*s.diag(*[s.Rational(v.numerator,v.denominator) for v in roots])
        shape=F(2,27)**m
        FE=HB*metric*s.Rational(shape.numerator,shape.denominator)
        seen=set(); candidates=[]
        for budget in [0,8*m,12*m,16*m,20*m]:
            cost=s.Rational(1,2**budget)
            feature=FE.row_join(B[:,-1:]*cost).row_join(s.Matrix([sum(row[:-1]) for row in B.tolist()])*cost)
            den=lcm(*(int(v.q) for v in feature))
            extended=(feature*den).row_join(B)
            reduced=extended.lll(delta=F(3,4))
            for row in reduced.tolist():
                R=list(map(int,row[-(c+1):])); alpha=R[-1]
                if not alpha: continue
                rat=F(-sum(R[:-1]),2**(c-3)*alpha)
                if rat<=0 or rat in seen: continue
                seen.add(rat); h=hcoeff(M,R,a,b,c); en=quadratic(h,G)
                factor=F(4*rat.denominator**2,abs(alpha))*shape
                candidates.append({'p':str(rat.numerator),'q':str(rat.denominator),
                  'alpha':str(alpha),'R':R[:-1],'H':h,'energy':str(en),
                  'gram_Q_bound_squared':str(factor*factor*en),
                  'gram_Q_bound_display':float(factor)*float(en)**.5,
                  'cost_budget':budget})
        candidates.sort(key=lambda r:F(r['gram_Q_bound_squared']))
        for r in candidates[:3]:
            cert=certificate(c,r['R'],int(r['alpha']),K=160)
            r['direct_Q_interval']=[cert['quality_lower'],cert['quality_upper']]
            r['direct_Q_upper_display']=cert['quality_upper_display']
            assert F(cert['quality_upper'])**2 < F(r['gram_Q_bound_squared'])
        out.append({'m':m,'a':a,'b':b,'c':c,'integer_kernel_rank':B.rows,
          'distinct_candidates':len(candidates),'top_by_gram_certificate':candidates[:3]})
    result={'status':'PASS','pi_or_cf_input':False,
      'scope':'Finite metric and generator audit. No new approximation record or quality-zero theorem.',
      'metric':'Gram G_ij=1/(2i+2j+1) for H(x^2); exact LDL with rational upper square roots.',
      'existing_candidate_bound_comparison':existing,'search':out,'elapsed_seconds':time.monotonic()-start}
    (Path(__file__).with_name('gram_wave_checks.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'elapsed_seconds':result['elapsed_seconds'],'comparisons':existing},indent=2))
    for x in out:
        print('m',x['m'],'distinct',x['distinct_candidates'],'top',[(r['p'],r['q'],r['gram_Q_bound_display'],r['direct_Q_upper_display']) for r in x['top_by_gram_certificate']])

if __name__=='__main__': main()
