#!/usr/bin/env python3
"""Check platform positions and feedback certificates; no new-fraction search."""
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp

from verify import balanced_energy, phi, screened_width

mp.mp.dps = 75
ROOT = Path(__file__).resolve().parent


def endpoints(a,b,xi):
    rho=mp.mpf(a)/b
    s=1/mp.mpf(xi)
    z=mp.exp(-s/rho)
    prefactor=2*mp.exp(s)
    count=max(1,int(mp.ceil(mp.log(mp.mpf('1e-52')*(1-z)/prefactor)/mp.log(z))))
    total=mp.mpf(0)
    for k in range(1,count+1):
        d=(k*b)//a
        vd=mp.exp(-s*d)/d
        vn=mp.exp(-s*(d+1))/(d+1)
        total+=(d+1)*vd-d*vn
    tail=prefactor*z**(count+1)/(1-z)
    assert tail < mp.mpf('1e-50')
    width,wtail=screened_width(b,xi)
    return total,total+width,tail+wtail


def finite_part(a,b):
    rho=mp.mpf(a)/b
    ds=[(j*b)//a for j in range(1,a+1)]
    return (rho-2*rho*(mp.log(b)+mp.euler)
            -mp.fsum(mp.digamma(mp.mpf(d)/b)+mp.digamma(mp.mpf(d+1)/b)
                     for d in ds)/b)


def integral_s(z,fermion=False):
    # At z>=0.15, the omitted t>8 tail is below exp(-1400).
    z=mp.mpf(z)
    a=2*mp.pi*z
    if fermion:
        f=lambda t: 1/(mp.exp(a*mp.cosh(t))+1)
    else:
        f=lambda t: 1/mp.expm1(a*mp.cosh(t))
    return mp.quad(f,[0,1,2,4,8])


def run():
    result={}
    rows=[]
    for a,b in [(1,2),(1,3),(2,5),(7,22)]:
        lower,upper,tail=endpoints(a,b,5)
        rho=mp.mpf(a)/b
        eps=mp.mpf('1e-17')
        left=(balanced_energy(rho-eps,5)-balanced_energy(rho-2*eps,5))/eps
        right=(balanced_energy(rho+2*eps,5)-balanced_energy(rho+eps,5))/eps
        assert abs(lower-left) < mp.mpf('1e-40')
        assert abs(upper-right) < mp.mpf('1e-40')
        rows.append({'density':f'{a}/{b}','mu_minus':str(lower),'mu_plus':str(upper)})
    result['position_vs_energy_slopes']={'passed':True,'rows':rows}

    convergence=[]
    for a,b in [(1,2),(1,3),(2,5),(7,22)]:
        target=finite_part(a,b)
        errors=[]
        for xi in [10,40,160]:
            lo,hi,tail=endpoints(a,b,xi)
            renormalized=lo-2*mp.mpf(a)/b*mp.log(xi)
            error=abs(renormalized-target)
            errors.append(error)
            convergence.append({'density':f'{a}/{b}','xi':xi,
                                'renormalized_lower':str(renormalized),
                                'predicted_limit':str(target),'absolute_error':str(error)})
            if a==1 and b==2:
                center=-mp.log1p(-mp.exp(-mp.mpf(1)/xi))
                assert abs((lo+hi)/2-center) < mp.mpf('1e-45')
        assert errors[2]<errors[1]<errors[0]
    assert abs(finite_part(1,2)+mp.mpf('.5')) < mp.mpf('1e-65')
    result['renormalized_positions']={'passed':True,'rows':convergence}

    # Independent finite sums with a rigorous elementary tail bound.
    series=[]
    for a,b in [(1,3),(2,5),(7,22)]:
        rho=mp.mpf(a)/b
        nmax=1000
        partial=mp.fsum(1/mp.mpf((k*b)//a)+1/mp.mpf((k*b)//a+1)-2*rho/k
                        for k in range(1,nmax+1))
        predicted=finite_part(a,b)-rho-2*rho*mp.log(rho)
        tail=rho*rho/(nmax*(1-rho/(nmax+1)))
        assert abs(partial-predicted)<tail
        series.append({'density':f'{a}/{b}','difference':str(partial-predicted),
                       'absolute_tail_bound':str(tail)})
    rational_bound=Fraction(25,28)-Fraction(21,22)+Fraction(150,2573)
    assert rational_bound < 0
    pi_bound=mp.mpf(25)/28-3/mp.pi+3/(2*mp.pi*(3*mp.pi-1))
    result['bias_coefficient']={'passed':True,'series_checks':series,
                                'pi_two_term_upper_bound':str(pi_bound),
                                'negative_rational_upper_bound':str(rational_bound)}

    # Finite screening identity at two existing historical benchmark fractions.
    alpha=1/mp.pi
    bridge=[]
    for p,q in [(22,7),(355,113)]:
        r=mp.mpf(q)/p
        c=q*(p-mp.pi*q)
        assert c>0
        for xi in [10,100,p]:
            lo,hi,tail=endpoints(q,p,xi)
            kminus=lo/(alpha-r)
            kplus=hi/(alpha-r)
            width,wtail=screened_width(p,xi)
            measured=(kplus-kminus)*c
            predicted=mp.pi*p*q*width
            assert abs(measured-predicted)<mp.mpf('1e-45')
            bridge.append({'known_fraction':f'{p}/{q}','xi':xi,
                           'pi_quality':str(c),'kappa_minus':str(kminus),
                           'kappa_plus':str(kplus),'kappa_width':str(kplus-kminus),
                           'width_times_quality':str(measured),
                           'exact_normalization':str(predicted)})
    result['pi_quality_bridge']={'passed':True,'not_a_blind_discovery':True,'rows':bridge}

    # Check the actual total feedback energy against nearby densities.
    a,b=2,7
    rho=mp.mpf(a)/b
    lo,hi,tail=endpoints(a,b,5)
    kappa=(lo+hi)/(2*(alpha-rho))
    objective=lambda x: balanced_energy(x,5)+kappa*(x-alpha)**2/2
    f0=objective(rho)
    for step in ['-0.01','-0.001','-0.00001','0.00001','0.001','0.01']:
        assert objective(rho+mp.mpf(step))>f0
    result['feedback_objective']={'passed':True,'density':'2/7','kappa':str(kappa),
                                  'nearby_competitors_checked':6}

    bessel=[]
    for z in ['0.15','0.5','1','2']:
        sz=integral_s(z)
        dz=integral_s(z,fermion=True)
        bound=1/(mp.pi**4*mp.mpf(z)**4)
        assert 0<dz<sz<bound
        bessel.append({'z':z,'S':str(sz),'D':str(dz),'upper_bound':str(bound)})
    lam=mp.mpf('.01')
    threshold=mp.sqrt(8*lam)/mp.pi**2
    quality_22=7*(22-7*mp.pi)
    quality_355=113*(355-113*mp.pi)
    # Choose h=1/113 only to check the known state, not to claim generation.
    scaled_energy=quality_355**2/2-4*lam*integral_s(1,fermion=True)
    assert quality_22>threshold>quality_355
    assert scaled_energy<0
    result['negative_energy_filter']={'passed':True,'bessel_rows':bessel,
                                     'lambda':str(lam),'uniform_quality_threshold':str(threshold),
                                     'known_355_scaled_min_energy':str(scaled_energy),
                                     'not_a_dynamics_or_discovery_test':True}
    return result


if __name__=='__main__':
    data=run()
    (ROOT/'feedback_verification.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({key: {'passed':value['passed']} for key,value in data.items()},indent=2))
