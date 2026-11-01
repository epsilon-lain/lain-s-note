"""Exact parameter checks for RADIAL_PACKET_EMBEDDING.md.

Checks finite members of the analytically proved high-dimensional family.
Does not numerically certify the uploaded interpolation theorem.
"""
from fractions import Fraction as F
from pathlib import Path
import json


def main():
    cases=[]
    for m in (3,4,8,16,32,64):
        K=3**m//2**m
        theta=F(1,2)
        w0=F(2**m)
        v0=F(2*K)
        T=theta*(1+F(1,4*m))
        theta_new=theta*(1+F(1,2*m))
        v0_new=F(9,10)*v0
        assert v0<w0 and K*theta**m<1
        assert K*w0/v0*theta**m==F(1,2)
        assert theta<T<theta_new<1 and v0_new<v0
        assert K*theta_new**m<1
        assert K*w0/v0_new*theta_new**m<1
        kappa=max(v0_new/v0,T/theta_new)
        inclusion_constant=T/(theta*theta_new)/(1-kappa)
        assert inclusion_constant==8*m+2
        u=theta/T
        linear=2*F(K-1,K)*m*(T-theta)/(m+2)
        collision=F(2,K*v0*(m+2))*(1+(K-1)*u**m*(K+m+1-m*u))
        b=m*theta/(m+2)
        assert linear<F(1,4*(m+2))
        assert collision<=F(1,m+2)
        # c=ln(2)/4<1/4, giving an exact rational cost upper bound.
        upper=linear+collision/4
        assert upper<F(1,2*(m+2))<1-b
        cases.append({'m':m,'K':K,'w0':str(w0),'v0':str(v0),
                      'theta':str(theta),'T':str(T),
                      'theta_new':str(theta_new),'v0_new':str(v0_new),
                      'v0_over_w0':str(v0/w0),
                      'new_fiber_volume_ratio':str(K*theta_new**m),
                      'new_full_packet_volume_ratio':str(K*w0/v0_new*theta_new**m),
                      'H_required_greater_than_sum_weights_times':str(inclusion_constant),
                      'normalized_linear_cost':str(linear),
                      'normalized_collision_coefficient':str(collision),
                      'strict_cost_upper_bound_using_c_less_than_one_quarter':str(upper),
                      'arithmetic_threshold_without_errors':str(1-b)})
    result={'status':'passed','arithmetic':'exact Fraction',
            'cases':cases,
            'scope':'Finite parameter verification; theorem-based full rank remains conditional on the uploaded interpolation theorem and simultaneous weight thresholds. This family has B=theta, not the manuscript critical B>theta path.'}
    Path(__file__).with_name('radial_parameter_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(cases)}))


if __name__=='__main__':
    main()
