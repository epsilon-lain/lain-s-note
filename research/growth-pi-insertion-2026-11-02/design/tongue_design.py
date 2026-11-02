"""Fixed-feedback symbolic drive intervals for an exact rational circle orbit.
No target digits are used to construct the examples. Standard library only.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import gcd
import json

def floor(v):return v.numerator//v.denominator

def branch(v):
    k=floor(v);u=v-k
    if u<Q(1,4):return Q(k)-Q(1,4),Q(k)+Q(1,4),Q(1),-Q(k)
    if u<Q(3,4):return Q(k)+Q(1,4),Q(k)+Q(3,4),Q(-1),Q(k)+Q(1,2)
    return Q(k)+Q(3,4),Q(k)+Q(5,4),Q(1),-Q(k+1)

def F(v,omega,kappa):
    lo,hi,s,c=branch(v)
    return (1-kappa*s)*v+omega-kappa*c

def symbolic_drive_interval(omega,kappa,p,q,x0):
    A=Q(1);B=Q(0);C=Q(0);x=x0
    affine=[];trace=[]
    for j in range(q):
        lo,hi,s,c=branch(x)
        affine.append((A,B,C));trace.append((lo,hi,s,c))
        m=1-kappa*s
        A,B,C=m*A,m*B+1,m*C-kappa*c
        x=F(x,omega,kappa)
    assert x==x0+p
    assert A!=1
    # For a fixed itinerary x0 = slope*Omega + intercept.
    slope=B/(1-A);intercept=(C-p)/(1-A)
    lower=None;upper=None;constraints=[]
    for (Aj,Bj,Cj),(lo,hi,s,c) in zip(affine,trace):
        aj=Aj*slope+Bj;bj=Aj*intercept+Cj
        if aj==0:
            assert lo<=bj<=hi;continue
        l,h=(lo-bj)/aj,(hi-bj)/aj
        if l>h:l,h=h,l
        lower=l if lower is None else max(lower,l)
        upper=h if upper is None else min(upper,h)
        constraints.append({'x_slope':str(aj),'x_intercept':str(bj),'branch':[str(lo),str(hi)],'omega_interval':[str(l),str(h)]})
    assert lower<=omega<=upper
    # Verify all endpoint orbit identities, allowing either compatible breakpoint branch.
    def psi(v):
        u=v-floor(v)
        if u<=Q(1,4):return u
        if u<=Q(3,4):return Q(1,2)-u
        return u-1
    for om in (lower,omega,upper):
        xx=slope*om+intercept;start=xx
        for _ in range(q):xx=xx+om-kappa*psi(xx)
        assert xx==start+p
    return {'omega':str(omega),'kappa':str(kappa),'p':p,'q':q,'primitive':gcd(p,q)==1,
      'initial_point':str(x0),'orbit_multiplier':str(A),
      'omega_interval':[str(lower),str(upper)],'width':str(upper-lower),
      'initial_point_affine':[str(slope),str(intercept)],'constraints':constraints,
      'parameter_bit_lengths':{'omega_numerator':abs(omega.numerator).bit_length(),'omega_denominator':omega.denominator.bit_length(),
                               'kappa_numerator':kappa.numerator.bit_length(),'kappa_denominator':kappa.denominator.bit_length()}}

def atan_inv(n,K):
    x=sum((Q((-1)**j,(2*j+1)*n**(2*j+1)) for j in range(K)),Q(0))
    rem=Q(1,(2*K+1)*n**(2*K+1))
    return (x,x+rem) if K%2==0 else (x-rem,x)

def pi_interval():
    a,b=atan_inv(5,60),atan_inv(239,20)
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]

def quality(r):
    pilo,pihi=pi_interval();p,q=r.numerator,r.denominator
    lo,hi=pilo-r,pihi-r
    low=Q(0) if lo<=0<=hi else min(abs(lo),abs(hi))
    high=max(abs(lo),abs(hi))
    return {'p':p,'q':q,'Q_lower':str(q*q*low),'Q_upper':str(q*q*high),'Q_display':float(q*q*(low+high)/2)}

def run():
    cases=[(Q(51,100),Q(1,2),1,2,Q(1,10)),(Q(17,50),Q(1,2),1,3,Q(28,125)),(Q(7,40),Q(1,2),1,7,Q(89,940))]
    results=[symbolic_drive_interval(*case) for case in cases]
    baseline=quality(3+Q(7,40));output=quality(3+Q(1,7))
    assert Q(output['Q_upper'])<Q(baseline['Q_lower'])
    out={'status':'PASS','scope':'Closed symbolic drive intervals from fixed rational feedback and exact itineraries; no new pi record or infinite construction.',
      'examples':results,'pi_posthoc_benchmark':{'input':baseline,'output':output,'input_used_pi':False,'known_fraction_only':True,'denominator_increased':False}}
    Path(__file__).with_name('checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'intervals':[{'rho':f"{r['p']}/{r['q']}",'drive_interval':r['omega_interval'],'width':r['width']} for r in results],'benchmark':out['pi_posthoc_benchmark']},indent=2))

if __name__=='__main__':run()
