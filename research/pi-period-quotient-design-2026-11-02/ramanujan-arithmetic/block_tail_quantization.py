"""A factor-preserving tail correction of the Ramanujan period series.

No pi value or continued fraction enters generation. The only external analytic
input is Ramanujan's identity, Guillera2018 formula(2). A published integer
divisibility is checked on the bounded controls but not inferred from them.
"""
from fractions import Fraction as F
from math import comb, gcd
from pathlib import Path
import json

BBASE = 4096


def term(n):
    return F((42*n+5)*comb(2*n,n)**3, BBASE**n)


def valuation_two(x):
    assert x > 0
    return (x & -x).bit_length()-1


def generate(N):
    assert N >= 2
    B = comb(2*N,N)
    d = 42*N+5
    # A geometric bound certifies an interval of width <=1/4 on the
    # quantization variable. Selection uses only B and exact integers.
    L = 1
    while 4*B > BBASE*63*64**(L-1):
        L += 1
    a = term(N)
    G = sum((term(N+j)/a for j in range(L)), F(0))
    # Use the proved geometric envelope, not an additional source coefficient.
    U = F(1,63*64**(L-1))
    lower, upper = F(B,BBASE)*G, F(B,BBASE)*(G+U)
    assert 0 < lower < upper and upper-lower <= F(1,4)
    middle = (lower+upper)/2
    u = (middle+F(1,2)).numerator // (middle+F(1,2)).denominator
    assert u >= 0
    assert max(abs(F(u)-lower),abs(F(u)-upper)) <= F(5,8)
    prefix = sum((term(k) for k in range(N)), F(0))
    C = prefix*BBASE**(N-1)
    assert C.denominator == 1
    C = int(C)
    assert C % (2*N*B) == 0  # He2017 theorem, finite implementation audit.
    D = C+u*d*B*B
    v = B-BBASE*u
    Cnext = BBASE*C+d*B**3
    unreduced_period = u*Cnext+v*C
    unreduced_integer = 16*BBASE**(N-1)*(BBASE*u+v)
    assert BBASE*u+v == B
    assert unreduced_period == B*D
    assert unreduced_integer == B*(16*BBASE**(N-1))
    rat = F(unreduced_integer,unreduced_period)
    assert rat == F(16*BBASE**(N-1),D)
    dyadic_gcd = gcd(16*BBASE**(N-1),D)
    assert dyadic_gcd == 2**valuation_two(D)
    assert rat.denominator == D//dyadic_gcd
    S_tilde = F(D,BBASE**(N-1))
    assert 5 <= S_tilde < 6
    error_upper = F(2,5)*d*B*B/BBASE**(N-1)
    joint_upper = rat.denominator**2*error_upper
    return dict(N=N,tail_terms=L,source_coefficient_max_index=N+L-1,
                B=B,u=u,v=v,prefix_integer=C,corrected_integer=D,
                p=rat.numerator,q=rat.denominator,
                guaranteed_common_factor=B,
                final_dyadic_gcd=dyadic_gcd,
                tail_coordinate_lower=lower,tail_coordinate_upper=upper,
                error_upper=error_upper,joint_upper=joint_upper,
                direct_prefix=F(16,prefix),candidate=rat)


def pi_interval(last=140):
    S = sum((term(k) for k in range(last+1)),F(0))
    a = term(last+1)
    return F(16)/(S+F(64,63)*a),F(16)/(S+a)


def quality(r,lo,hi):
    assert not lo <= r <= hi
    er = (lo-r,hi-r) if r < lo else (r-hi,r-lo)
    return tuple(x*r.denominator**2 for x in er)


def run():
    # Polynomial ratio proof: a_(n+1)/a_n <1/64 for n>=1.
    # The cross-multiplied difference is168n^3+312n^2+132n-7>0.
    for n in range(1,80):
        assert term(n+1)/term(n) < F(1,64)
        assert 168*n**3+312*n*n+132*n-7 > 0
    plo,phi = pi_interval()
    rows = []
    for N in range(2,13):
        r = generate(N)
        qlo,qhi = quality(r['candidate'],plo,phi)
        assert qhi <= r['joint_upper']
        M = r['source_coefficient_max_index']
        direct = F(16)/sum((term(k) for k in range(M+1)),F(0))
        dlo,dhi = quality(direct,plo,phi)
        rows.append(dict(N=N,tail_terms=r['tail_terms'],
                         source_coefficient_max_index=M,
                         u=r['u'],p=str(r['p']),q=str(r['q']),
                         guaranteed_common_factor=str(r['B']),
                         final_dyadic_gcd=str(r['final_dyadic_gcd']),
                         strict_quality_lower=str(qlo),strict_quality_upper=str(qhi),
                         certified_joint_upper=str(r['joint_upper']),
                         displayed_quality=float(qhi),
                         same_input_direct_q=str(direct.denominator),
                         same_input_quality_improvement_lower=str(dlo/qhi),
                         displayed_same_input_improvement=float(dlo/qhi),
                         changed_from_prefix=(r['candidate'] != r['direct_prefix'])))
    out=dict(status='passed',controls='N=2,...,12; no parameter search',
             source='Guillera2018 formula2; He2017 published divisibility theorem',
             rows=rows,
             scope='The common factor is genuine but coefficient height is charged. The method is a binomial-grid tail quantization; no quadratic-quality limit is inferred.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ('N','tail_terms','displayed_quality','displayed_same_input_improvement','changed_from_prefix')} for r in rows],indent=2))


if __name__ == '__main__':
    run()
