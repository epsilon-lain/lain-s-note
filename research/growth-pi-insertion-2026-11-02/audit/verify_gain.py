"""Independent arithmetic and absolute-coordinate replay of the insertion step.
Construction engine is used only as an input certificate provider; the checks below
use a different lift formula and an independent Machin interval implementation.
"""
from fractions import Fraction as Q
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'algorithm'))
from itinerary_certificate import propose

def atan_interval(q,N):
    a=sum((Q((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N)),Q(0))
    b=a+Q((-1)**N,(2*N+1)*q**(2*N+1))
    return min(a,b),max(a,b)

def pi_interval():
    a,b=atan_interval(5,80),atan_interval(239,24)
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]

def quality(r,interval):
    p,q=r.numerator,r.denominator
    lo,hi=interval[0]-r,interval[1]-r
    assert lo*hi>0
    return min(abs(lo),abs(hi))*q*q,max(abs(lo),abs(hi))*q*q

def true_lift(x,omega,kappa):
    cell=x.numerator//x.denominator;u=x-cell
    if u<=Q(1,4):wave=u
    elif u<=Q(3,4):wave=Q(1,2)-u
    else:wave=u-1
    return x+omega-kappa*wave

def run():
    omega,kappa=Q(1,7),Q(493,5000)
    cert=propose(omega,kappa,max_period=256,burn=20000)
    assert cert and (cert['p'],cert['q'])==(16,113)
    x0=Q(cert['period_point']);x=x0
    for j in range(113):x=true_lift(x,omega,kappa)
    assert x==x0+16
    A=Q(cert['multiplier']);L,R=map(Q,cert['phase_interval'])
    assert 0<A<1 and L<x0<R
    radius=min(x0-L,R-x0)/2
    images=[]
    for start in (x0-radius,x0+radius):
        x=start
        for j in range(113):x=true_lift(x,omega,kappa)
        images.append(x-16)
    assert images==[x0-A*radius,x0+A*radius]
    assert x0-radius<images[0]<images[1]<x0+radius
    pi=pi_interval();assert pi[1]-pi[0]<Q(1,10**100)
    old=Q(3)+omega;new=Q(3)+Q(cert['p'],cert['q'])
    oldQ,newQ=quality(old,pi),quality(new,pi)
    improvement=(oldQ[0]/newQ[1],oldQ[1]/newQ[0])
    assert new.denominator>old.denominator and improvement[0]>18
    qratio=Q(new.denominator,old.denominator)
    relative_error_upper=(newQ[1]/new.denominator**2)/(oldQ[0]/old.denominator**2)
    assert relative_error_upper<1/qratio**2
    out={'status':'PASS','input':{'omega':str(omega),'kappa':str(kappa),'numerator_denominator_total_bits':sum(v.numerator.bit_length()+v.denominator.bit_length() for v in (omega,kappa))},
        'old':{'fraction':str(old),'quality_interval':list(map(str,oldQ)),'quality_display':float(oldQ[1])},
        'new':{'fraction':str(new),'quality_interval':list(map(str,newQ)),'quality_display':float(newQ[1])},
        'denominator_growth':str(qratio),'strict_quality_improvement_lower':str(improvement[0]),'improvement_display':float(improvement[0]),
        'relative_absolute_error_upper':str(relative_error_upper),
        'independent_absolute_lift_checks':113,'independent_trapping_endpoint_checks':226,
        'period_point_numerator_bits':x0.numerator.bit_length(),'period_point_denominator_bits':x0.denominator.bit_length(),
        'trapping_radius':str(radius),'trapping_margin':str((1-A)*radius),
        'pi_interval_width':str(pi[1]-pi[0]),
        'scope':'A certified finite gain from 22/7 to known 355/113. The parameter was suggested in exploration; no novel record or asymptotic theorem.'}
    (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'PASS','input':out['input'],'old_fraction':str(old),'new_fraction':str(new),'old_Q':float(oldQ[1]),'new_Q':float(newQ[1]),'q_ratio':str(qratio),'gain_factor_lower_display':out['improvement_display'],'period_point_bits':[x0.numerator.bit_length(),x0.denominator.bit_length()]},indent=2))

if __name__=='__main__':run()
