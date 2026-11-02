"""Independent certificates for exact rational stepping and changed period ratios.
No imports from the construction modules. Standard library only.
"""
from fractions import Fraction as F
from pathlib import Path
import json

def floor(v):return v.numerator//v.denominator

def psi(v):
    u=v-floor(v)
    if u<=F(1,4):return u
    if u<=F(3,4):return F(1,2)-u
    return u-1

def step(v,omega):return v+omega-psi(v)/2

def rational_step(v,a=F(1)):
    n=floor(v/a+F(1,2));e=v/a-n
    return a*(n+1+4*e**3)

def run():
    cases=[{'omega':F(51,100),'p':1,'q':2,'orbit':[F(1,10),F(14,25),F(11,10)],'I':(F(9,100),F(11,100)),'multiplier':F(3,4),'margin':F(1,400)},
           {'omega':F(17,50),'p':1,'q':3,'orbit':[F(28,125),F(113,250),F(96,125),F(153,125)],'I':(F(107,500),F(117,500)),'multiplier':F(3,8),'margin':F(1,160)}]
    report=[]
    for c in cases:
        om=c['omega'];orbit=c['orbit']
        assert all(step(a,om)==b for a,b in zip(orbit,orbit[1:]))
        assert orbit[-1]==orbit[0]+c['p']
        lo,hi=c['I'];derivative=F(1)
        # Check that each entire interval lies inside a single affine branch.
        for _ in range(c['q']):
            breakpoints=[F(k)+d for k in range(floor(lo)-1,floor(hi)+2) for d in (F(1,4),F(3,4))]
            assert not any(lo<t<hi for t in breakpoints)
            newlo,newhi=step(lo,om),step(hi,om)
            derivative*=(newhi-newlo)/(hi-lo)
            lo,hi=newlo,newhi
        lo-=c['p'];hi-=c['p'];I0,I1=c['I']
        assert derivative==c['multiplier'] and 0<derivative<1
        assert min(lo-I0,I1-hi)==c['margin']>0
        report.append({'omega':str(om),'rho':str(F(c['p'],c['q'])),'exact_orbit':[str(x) for x in orbit],
           'trapping_interval':[str(I0),str(I1)],'q_step_contraction':str(derivative),'margin':str(c['margin'])})
    plateau=0
    for j in range(201):
        om=F(19,40)+F(j,4000)
        x=10*om-5
        assert step(step(x,om),om)==x+1
        plateau+=1
    # Independently check the cubic orbit closed form and denominator growth.
    rational_cases=0
    for a in (F(1),F(2),F(3),F(7)):
        for u0 in (F(-2,5),F(-1,7),F(0),F(1,7),F(2,5)):
            v=a*u0
            for k in range(6):
                err=v/a-k
                assert err==(2*u0)**(3**k)/2
                if u0==F(1,7):assert err.denominator==7**(3**k)
                v=rational_step(v,a)
                rational_cases+=1
    # A rational bound for the Fejer smoothing perturbation, including three iterates.
    smooth_error=(F(8,1025)+F(1,1024))/9
    assert smooth_error<F(1,1000)
    for c in cases:
        perturbation=sum((F(3,2)**k for k in range(c['q'])),F(0))*smooth_error
        assert perturbation<c['margin']
    out={'status':'PASS','arithmetic':'standard-library Fraction only','period_certificates':report,
         'closed_locking_plateau_checks':plateau,'rational_cubic_closed_form_checks':rational_cases,
         'fejer_step_error_upper':str(smooth_error),'scope':'Independent exact orbit/trapping/denominator controls. No improved approximation to pi is asserted.'}
    Path(__file__).with_name('independent_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':run()
