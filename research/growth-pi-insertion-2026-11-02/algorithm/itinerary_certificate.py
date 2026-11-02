"""Linear-pass exact certificates from floating-point itinerary proposals.
Floats only choose candidate periods and branches; every certificate is checked
again by exact rational inequalities. No target irrational is used by this module.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import floor,gcd
import json

def _floor(v):return v.numerator//v.denominator

def float_step(u,omega,kappa):
    if u<.25:psi=u
    elif u<.75:psi=.5-u
    else:psi=u-1
    v=u+omega-kappa*psi
    j=floor(v)
    return v-j,j

def _branch(u):
    # Closed branch boundaries are safe since the two formulae agree there.
    if u<.25:return Q(0),Q(1,4),0
    if u<.75:return Q(1,4),Q(3,4),1
    return Q(3,4),Q(1),2

def _affine(omega,kappa,branch):
    if branch==0:return 1-kappa,omega
    if branch==1:return 1+kappa,omega-kappa/2
    return 1-kappa,omega+kappa

def certify_itinerary(omega,kappa,p,q,start):
    """Return exact certified cycle or None. start may be any float phase.

    This is O(q) rational operations, not O(q) bit time. The certificate can be
    reverified without using start or any float: replay the exact initial point.
    """
    omega,kappa=Q(omega),Q(kappa)
    if not (0<kappa<1) or q<1:return None
    # The rounded orbit determines only a proposal itinerary.
    u=float(start)%1.;A=Q(1);B=Q(0);lo=Q(0);hi=Q(1);lift=0;word=[]
    for _ in range(q):
        l,r,branch=_branch(u)
        lo=max(lo,(l-B)/A);hi=min(hi,(r-B)/A)
        if lo>hi:return None
        a,b=_affine(omega,kappa,branch)
        # Both use the same proposed branch, even if rounding lands on a corner.
        v=float(a)*u+float(b);j=floor(v);u=v-j
        A,B=a*A,a*B+b-j;lift+=j
        # The integer-cell part of the proposal is also checked exactly.
        lo=max(lo,-B/A);hi=min(hi,(1-B)/A)
        if lo>hi:return None
        word.append((branch,j))
    if lift!=p:return None
    if A==1:
        if B!=0:return None
        x=(lo+hi)/2
    else:x=-B/(A-1)
    if not lo<=x<=hi:return None
    u=x;actual_p=0;maxbits=0
    for branch,j in word:
        l,r=((Q(0),Q(1,4)),(Q(1,4),Q(3,4)),(Q(3,4),Q(1)))[branch]
        if not l<=u<=r:return None
        a,b=_affine(omega,kappa,branch)
        nxt=a*u+b-j
        if not 0<=nxt<=1:return None
        actual_p+=j;u=nxt
        maxbits=max(maxbits,u.numerator.bit_length(),u.denominator.bit_length())
    if u!=x or actual_p!=p:return None
    g=gcd(p,q)
    return {'p':p//g,'q':q//g,'tested_period':q,'period_point':str(x),
      'phase_interval':[str(lo),str(hi)],'multiplier':str(A),'strict_local_attractor':A<1 and lo<x<hi,
      'integer_lift':p,'word':word,'maximum_orbit_rational_bits':maxbits,
      'input_omega':str(omega),'input_kappa':str(kappa),'status':'CERTIFIED'}

def verify_certificate(cert):
    """Independent exact replay; accepts a JSON round-tripped certificate.
    It never consults floating-point data or a target irrational.
    """
    omega=Q(cert['input_omega']);kappa=Q(cert['input_kappa'])
    assert 0<kappa<1
    x=Q(cert['period_point']);lo,hi=map(Q,cert['phase_interval'])
    assert lo<=x<=hi
    A=Q(1);B=Q(0);p=0
    for branch,j in cert['word']:
        assert branch in (0,1,2) and isinstance(j,int)
        l,r=((Q(0),Q(1,4)),(Q(1,4),Q(3,4)),(Q(3,4),Q(1)))[branch]
        assert l<=A*lo+B<=A*hi+B<=r
        a,b=_affine(omega,kappa,branch)
        A,B=a*A,a*B+b-j
        assert 0<=A*lo+B<=A*hi+B<=1
        p+=j
    assert A*x+B==x
    assert A==Q(cert['multiplier'])
    assert p==cert['integer_lift']
    assert cert['tested_period']==len(cert['word'])
    ratio=Q(p,len(cert['word']))
    assert ratio==Q(cert['p'],cert['q'])
    assert cert['strict_local_attractor']==(A<1 and lo<x<hi)
    return True

def propose(omega,kappa,max_period=1000,burn=20000,top=16):
    """Generate candidate p/q from near-returns, never from pi or a target fraction."""
    om,kap=float(omega),float(kappa);u=.123456789
    for _ in range(burn):u,_=float_step(u,om,kap)
    start=u;lift=0;candidates=[]
    for q in range(1,max_period+1):
        u,j=float_step(u,om,kap);lift+=j
        # Compare lifts modulo one; an end phase near 1 can match a start near 0.
        delta=u-start;shift=round(delta)
        dist=abs(delta-shift);p=lift+shift
        candidates.append((dist,q,p))
    for dist,q,p in sorted(candidates)[:top]:
        cert=certify_itinerary(omega,kappa,p,q,start)
        if cert is not None:return cert
    return None

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('omega');ap.add_argument('kappa')
    ap.add_argument('--max-period',type=int,default=1000);ap.add_argument('--burn',type=int,default=20000)
    args=ap.parse_args()
    cert=propose(Q(args.omega),Q(args.kappa),args.max_period,args.burn)
    if cert:cert={k:v for k,v in cert.items() if k!='word'}
    print(json.dumps(cert,indent=2))
