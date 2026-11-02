"""Exact global minimization of the gcd-charged Gram score on an integer chain.

The candidate generator uses only rational moments and integer periods. It does
not read the saved inputs' direct pi-quality intervals or any continued fraction.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd,lcm,comb
import json
import sympy as s
from integral_certificate import certificate

HERE=Path(__file__).resolve().parent
k,t=s.symbols('k t')

def sf(x):
    x=s.Rational(x)
    return F(int(x.p),int(x.q))

def linear_congruence(v,w,d):
    # v*k+w == 0 mod d.
    h=gcd(v,d)
    if w%h:return None
    mod=d//h
    return (0,1) if mod==1 else ((-w//h)*pow(v//h,-1,mod)%mod,mod)

def combine(x,y):
    if x is None or y is None:return None
    a,m=x;b,n=y;h=gcd(m,n)
    if (b-a)%h:return None
    modulus=lcm(m,n);nn=n//h
    shift=0 if nn==1 else ((b-a)//h)*pow(m//h,-1,nn)%nn
    return (a+m*shift)%modulus,modulus

def critical_integer_candidates(poly,variable):
    """Cover every integer global minimizer of a coercive polynomial,
    after deletion of the possible integer alpha=0 point.
    Root intervals are exact rational isolating intervals, not floats.
    """
    der=s.Poly(s.diff(poly,variable),variable)
    intervals=der.intervals(eps=s.Rational(1,16))
    points=set()
    for (left,right),multiplicity in intervals:
        # The extra neighbor covers deleting a forbidden integral alpha-root.
        for n in range(int(s.floor(left))-1,int(s.ceiling(right))+2):
            points.add(n)
    return sorted(points),[{'lo':str(a),'hi':str(b),'multiplicity':mult}
                           for (a,b),mult in intervals]

def run():
    data=json.loads((HERE/'inputs.json').read_text())
    m,c=data['m'],data['c'];a,b=data['rows']
    alpha0,alpha1=int(a['alpha']),int(b['alpha'])
    R0,R1=list(map(int,a['R'])),list(map(int,b['R']))
    H0,H1=list(map(int,a['H'])),list(map(int,b['H']))
    x=s.symbols('x')
    for rr,al,hh in [(R0,alpha0,H0),(R1,alpha1,H1)]:
        R=x*sum(v*x**(2*j) for j,v in enumerate(rr))
        N=(1+x*x)*s.diff(R,x)-2*(c-1)*x*R+al*(1+x*x)**(c-1)
        expected=x**(2*m)*(1-x*x)**(2*m)*sum(v*x**(2*j) for j,v in enumerate(hh))
        assert s.expand(N-expected)==0
    dyad=2**(c-3)
    raw=[-sum(R0),dyad*alpha0,-sum(R1),dyad*alpha1]
    content=gcd(*raw)
    P0,A0,P1,A1=[x//content for x in raw]
    determinant=P0*A1-P1*A0
    assert determinant!=0
    alpha=alpha0+k*alpha1;AA=A0+k*A1
    HH=[x+k*y for x,y in zip(H0,H1)]
    energy=s.expand(sum(HH[i]*HH[j]/s.Integer(2*i+2*j+1)
                        for i in range(len(HH)) for j in range(len(HH))))
    # q=abs(A)/g and alpha=(content/dyad)*A, so no division by a polynomial.
    prefactor=s.Rational(16*dyad**2,content**2)*s.Rational(2,27)**(2*m)
    raw_score=s.Poly(s.expand(prefactor*AA**2*energy),k).as_expr()
    assert s.degree(raw_score,k) in (2,4)
    assert s.Poly(raw_score,k).LC()>0
    classes=[];candidate_set=set()
    for divisor in s.divisors(abs(determinant)):
        d=int(divisor)
        ap=combine(linear_congruence(A1,A0,d),linear_congruence(P1,P0,d))
        if ap is None:continue
        residue,modulus=ap
        objective=s.expand(raw_score.subs(k,residue+modulus*t)/d**4)
        points,roots=critical_integer_candidates(objective,t)
        valid=[n for n in points if alpha0+(residue+modulus*n)*alpha1!=0]
        assert valid
        best=min(valid,key=lambda n:(sf(objective.subs(t,n)),n))
        bestk=residue+modulus*best
        score=sf(objective.subs(t,best))
        for n in valid:candidate_set.add(residue+modulus*n)
        classes.append({'d':d,'residue':residue,'modulus':modulus,
                        'derivative_root_intervals':roots,'candidate_count':len(valid),
                        'best_k':bestk,'best_guaranteed_bound_squared':str(score)})
    chosen=min(classes,key=lambda row:(F(row['best_guaranteed_bound_squared']),row['best_k']))
    chosen_k=chosen['best_k']
    def actual_bound(kk):
        al=alpha0+alpha1*kk
        if not al:return None
        den=A0+A1*kk;num=P0+P1*kk;g=gcd(den,num)
        rat=F(num,den)
        en=sf(energy.subs(k,kk))
        score=sf(raw_score.subs(k,kk))/g**4
        factor=F(4*rat.denominator**2,abs(al))*F(2,27)**m
        assert score==factor**2*en
        return score,rat,g
    best_score,best_rat,best_g=actual_bound(chosen_k)
    # Every divisor score majorizes the actual score, and at d=g(k) equality
    # holds. Thus the global minimum over these classes is the true minimum.
    assert chosen['d']==best_g
    assert best_score==F(chosen['best_guaranteed_bound_squared'])
    audit_count=0
    # Independent wide-window arithmetic regression, not the global proof.
    for kk in range(-2000,2001):
        if alpha0+alpha1*kk==0:continue
        result=actual_bound(kk)
        assert result[0]>=best_score
        assert abs(determinant)%result[2]==0
        audit_count+=1
    R=[x+chosen_k*y for x,y in zip(R0,R1)]
    al=alpha0+chosen_k*alpha1
    check=certificate(c,R,al,K=180)
    assert (check['p'],check['q'])==(str(best_rat.numerator),str(best_rat.denominator))
    assert F(check['quality_upper'])**2<best_score
    initial=[]
    for row in (a,b):
        exact=certificate(c,list(map(int,row['R'])),int(row['alpha']),K=180)
        initial.append({key:exact[key] for key in ['p','q','quality_lower','quality_upper','quality_upper_display']})
    actual_improvement_ratio=F(check['quality_upper'])/F(initial[0]['quality_lower'])
    assert actual_improvement_ratio<F(1,100)
    output={'status':'PASS','scope':data['scope'],'m':m,'c':c,
            'period_pair':[P0,A0,P1,A1],'common_raw_content':content,
            'determinant':determinant,'divisor_count':len(s.divisors(abs(determinant))),
            'nonempty_progressions':len(classes),'distinct_candidates':len(candidate_set),
            'energy_polynomial':str(energy),'raw_score_polynomial':str(raw_score),
            'classes':classes,'global_best':{'k':chosen_k,'p':str(best_rat.numerator),
             'q':str(best_rat.denominator),'gcd':best_g,'alpha':str(al),'R':R,
             'H':[x+chosen_k*y for x,y in zip(H0,H1)],
             'Q_bound_squared':str(best_score),'Q_bound_display':float(best_score)**.5,
             'direct_integral_certificate':check},
            'input_integral_certificates':initial,
            'selected_step':{'old_reduced_q':initial[0]['q'],'new_reduced_q':check['q'],
              'gcd_gain':str(F(best_g,gcd(P0,A0))),
              'raw_denominator_ratio':str(F(abs(A0+A1*chosen_k),abs(A0))),
              'true_Q_ratio_certified_upper':str(actual_improvement_ratio),
              'true_Q_ratio_upper_display':float(actual_improvement_ratio),
              'is_growing_denominator_progress':False},
            'wide_window_regression_cases':audit_count,
            'selection_uses_pi_or_CF':False}
    (HERE/'checks.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':'PASS','period_pair':output['period_pair'],
                      'divisor_count':output['divisor_count'],
                      'nonempty_progressions':len(classes),'candidates':len(candidate_set),
                      'best_k':chosen_k,'best_rational':str(best_rat),
                      'Gram_Q_upper':float(best_score)**.5,
                      'integral_Q_upper':check['quality_upper_display'],
                      'window_checks':audit_count},indent=2))

if __name__=='__main__':run()
