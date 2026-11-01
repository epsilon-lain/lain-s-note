"""Exact finite checks for ENVELOPE_BARRIER.md; standard library only."""
from fractions import Fraction as F
from math import floor
from pathlib import Path
import json


def indices(weights, bound):
    if not weights:
        if bound > 0:
            yield ()
        return
    k = 0
    while k * weights[0] < bound:
        for rest in indices(weights[1:], bound-k*weights[0]):
            yield (k,)+rest
        k += 1


def finite_case(m, K, H, v0, w0, weights, theta, T, vp):
    rows=list(indices([v0]+[w/theta for w in weights],F(H)))
    M=K*len(rows)
    types=list(indices(weights,T*H))
    L=len(types)
    packet=K*sum(1 for _ in indices([vp]+[w/T for w in weights],F(H)))
    # Count columns with the non-strict boundary exactly, on a common denominator.
    from math import lcm
    den=lcm(*(w.denominator for w in [w0]+weights))
    cols=sum(1 for _ in indices([w0]+weights,F(H)+F(1,den)))
    occupancy={}
    for row in rows:
        beta=row[1:]
        assert sum(w*b for w,b in zip(weights,beta)) < T*H
        assert vp*row[0]+sum(w*b for w,b in zip(weights,beta))/T < H
        occupancy[beta]=occupancy.get(beta,0)+K
    ss=sum(n*n for n in occupancy.values())
    assert L*ss >= M*M
    return {'m':m,'K':K,'H':H,'original_rows':M,'columns':cols,
            'envelope_rows':packet,'possible_types':L,
            'sum_squared_occupancies':ss,'finite_cauchy_verified':True,
            'envelope_dimension_count_only':packet<=cols}


def main():
    finite=[]
    for H in [8,16,24,40]:
        finite.append(finite_case(2,2,H,F(1,2),F(1),[F(1),F(3)],F(1,4),F(3,10),F(2,5)))
    critical=[]
    theta=F(1,2); B=F(3,5); C=F(13,10)
    for m in [3,10,25,50,100]:
        K=floor(C**m); w0=B**(-m); v0=2*K*theta**m*w0
        assert K*theta**m<1
        assert K*w0*theta**m/v0==F(1,2)
        coefficient=K*K*theta**m*w0/(v0*v0*(m+1))
        assert coefficient==(B/theta)**m/(4*(m+1))
        lower=coefficient/6  # c=ln(2)/4 > 1/6, from ln(2)>2/3.
        if m==50: assert lower>7
        critical.append({'m':m,'K':K,'v0_over_w0':str(v0/w0),
                         'cost_lower_bound_using_c_gt_1_over_6':str(lower)})
    result={'status':'passed','finite_count_checks':finite,'critical_path_checks':critical,
            'scope':'Finite arithmetic checks support the proved asymptotic envelope restriction. No certification of interpolation thresholds or pi badly approximable.'}
    Path(__file__).with_name('verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','finite_cases':len(finite),'critical_cases':len(critical)}))


if __name__=='__main__':
    main()
