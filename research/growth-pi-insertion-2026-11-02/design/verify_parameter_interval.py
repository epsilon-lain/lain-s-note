"""Independent exact replay and a simple continuum of parameters for 16/113.
Uses only the standard library. No numerical trig or approximate comparison.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import gcd
import json
from tongue_design import branch,pi_interval,quality

HERE=Path(__file__).resolve().parent

def run():
    cert=json.loads((HERE/'target_certificate.json').read_text())
    omega=Q(cert['input_omega']);kappa=Q(cert['input_kappa'])
    p,q=cert['p'],cert['q']
    assert (omega,kappa,p,q)==(Q(1,7),Q(493,5000),16,113)
    assert gcd(p,q)==1 and cert['integer_lift']==p and cert['tested_period']==q
    x0=Q(cert['period_point']);lo,hi=map(Q,cert['phase_interval'])
    r=min(x0-lo,hi-x0)/2
    assert r>0
    I0,I1=x0-r,x0+r
    a,b=I0,I1;v=x0;A=Q(1)
    for _ in range(q):
        left,right,s,c=branch((a+b)/2)
        assert left<=a<=b<=right
        assert left<=v<=right
        m=1-kappa*s
        a,b=m*a+omega-kappa*c,m*b+omega-kappa*c
        v=m*v+omega-kappa*c
        A*=m
    assert v==x0+p and A==Q(cert['multiplier']) and 0<A<1
    a-=p;b-=p
    gap=min(a-I0,I1-b)
    assert gap==(1-A)*r>0
    # Every real kappa in this closed interval yields the same rotation number.
    delta=Q(1,10**10);M=Q(11,10)
    assert 0<kappa-delta<kappa+delta<Q(1,10)<1
    S=sum((M**j for j in range(q)),Q(0))
    perturbation=delta*S/4
    assert perturbation<gap
    remaining=gap-perturbation
    # Joint denominator/quality comparison; oracle is posthoc, not construction.
    baseline=quality(3+omega);output=quality(3+Q(p,q))
    assert output['q']>baseline['q']
    assert Q(output['Q_upper'])<Q(baseline['Q_lower'])/18
    report={'status':'PASS','parameter_interval':{'omega':str(omega),'kappa_center':str(kappa),
      'kappa_radius':str(delta),'lower':str(kappa-delta),'upper':str(kappa+delta),
      'all_real_parameters_included':True},
      'output_rotation':str(Q(p,q)),'output_fraction_for_pi':str(3+Q(p,q)),
      'primitive_period':q,'fixed_point':str(x0),'trapping_radius':str(r),
      'multiplier':str(A),'q_step_gap':str(gap),'global_lipschitz_bound':str(M),
      'q_step_perturbation_bound':str(perturbation),'remaining_strict_margin':str(remaining),
      'displays_only':{'radius':float(r),'multiplier':float(A),'gap':float(gap),
                      'perturbation':float(perturbation),'remaining_margin':float(remaining)},
      'pi_quality':{'input':baseline,'output':output,'improvement_factor_certified_lower':str(Q(baseline['Q_lower'])/Q(output['Q_upper']))},
      'scope':'Exact fixed-orbit replay, trapping and continuum perturbation certificate. Known 355/113 benchmark, not a new pi record or an infinite improvement theorem.'}
    (HERE/'parameter_interval_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','omega':str(omega),'kappa_center':str(kappa),'kappa_radius':str(delta),
                     'rotation':str(Q(p,q)),'fraction':str(3+Q(p,q)),
                     'input_q':baseline['q'],'output_q':output['q'],
                     'input_Q':baseline['Q_display'],'output_Q':output['Q_display'],
                     'improvement_factor_lower_display':float(Q(baseline['Q_lower'])/Q(output['Q_upper'])),
                     'displays_only':report['displays_only']},indent=2))

if __name__=='__main__':run()
