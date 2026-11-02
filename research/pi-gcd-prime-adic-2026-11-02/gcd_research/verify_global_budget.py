"""Exact finite audits of the global odd-surplus budget; no extrapolation.

Only diagnostic log summaries use floating point. Every inequality and
valuation assertion is checked with integers. No evaluation of pi is used.
"""
from pathlib import Path
from math import gcd, log
import json
import sys


def v2(x):
    assert x > 0
    return (x & -x).bit_length() - 1


def predicted_v2(n):
    k, r = divmod(n, 4)
    return (2*k, 2*k, 2*k+2, 2*k+3+v2(k+1))[r]


def main(N=5000):
    a0, a, b0, b, fac = 1, 1, 0, 1, 1
    best_h = (1, 1, 1)  # numerator, denominator of h^2, index
    max_density = (0.0, 1, 1)
    max_E = (1, 1)
    checkpoints, excess_count = [], 0
    all_excesses = []
    first_excess = []
    for n in range(1, N+1):
        if n > 1:
            k = n-1
            a0, a = a, (2*k+1)*a+k*k*a0
            b0, b = b, (2*k+1)*b+k*k*b0
            fac *= n
        gn = gcd(a, 4*b)
        assert v2(a) == predicted_v2(n)
        vf, vg = v2(fac), v2(gn)
        numerator = gn // gcd(gn, fac)
        E = numerator >> v2(numerator)
        # Exact refined necessary bound h <= sqrt(n)*2^(v2(A)-v2(n!))*E.
        exponent = 2*(v2(a)-vf)
        lhs, rhs = gn*gn, fac*fac*E*E
        if exponent >= 0:
            rhs <<= exponent
        else:
            lhs <<= -exponent
        assert lhs <= rhs
        # Coarser all-index bound h <= 2(n+1)^2 sqrt(n) 2^(-n/2) E.
        assert (gn*gn << n) <= 4*(n+1)**4*fac*fac*E*E
        if n >= 3:
            # Existing exact Legendre error bound implies Q_n>1/2 here.
            assert 4*fac*fac > (2*n+1)*gn*gn
        if gn*gn*n*best_h[1] > best_h[0]*fac*fac:
            best_h = (gn*gn*n, fac*fac, n)
        if E > 1:
            excess_count += 1
            all_excesses.append({'n':n,'E':str(E)})
            if len(first_excess) < 10:
                first_excess.append({'n': n, 'E': str(E)})
        density = log(E)/n
        if density > max_density[0]:
            max_density = (density, n, E)
        if E > max_E[0]:
            max_E = (E, n)
        if n in (72,107,166,500,1000,2000,N):
            checkpoints.append({'n':n, 'odd_surplus':str(E),
                                'log_h':log(gn)+0.5*log(n)-log(fac),
                                'v2_g':vg,'v2_factorial':vf,
                                'required_compensation_log_slope':log(2)/2})
    common = gcd(best_h[0],best_h[1])
    result={'status':'PASS','checked_n':N,'all_index_two_adic_checks':N,
            'exact_budget_checks':2*N,
            'maximum_h_squared':{'n':best_h[2],
                'numerator':str(best_h[0]//common),
                'denominator':str(best_h[1]//common)},
            'indices_with_odd_surplus':excess_count,
            'all_odd_surpluses':all_excesses,
            'first_odd_surpluses':first_excess,
            'maximum_observed_log_E_over_n':{'value':max_density[0],
                'n':max_density[1],'E':str(max_density[2])},
            'maximum_observed_E':{'n':max_E[1],'E':str(max_E[0])},
            'quality_above_one_half_certified_for_n':[3,N],
            'checkpoints':checkpoints,
            'scope':'Finite exact verification only. Diagnostic logarithms '
                    'are floating-point summaries. No infinite odd-surplus '
                    'bound or assertion about pi follows from the scan.'}
    out=Path(__file__).with_name('global_budget_checks.json')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 5000)
