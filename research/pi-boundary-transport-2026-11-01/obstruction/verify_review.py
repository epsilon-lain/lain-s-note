"""Independent exact checks for the boundary coefficient budget and columns.

Uses the sibling standard-library Gaussian-rational helper module. It does
not rely on the shear verifier's matrix or boundary implementations.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import json

from verify_gram_displacement import QI, determinant


def prod(values):
    out = F(1)
    for value in values:
        out *= value
    return out


def coefficient(beta, gamma, costs, gap):
    alt = sum((-1)**sum(kappa)
              * prod(comb(g,k) for g,k in zip(gamma,kappa))
              for kappa in product(*(range(g+1) for g in gamma))
              if sum(c*k for c,k in zip(costs,kappa)) < gap)
    # lambda_i = 2^{-cost_i}; powers remain exact rational numbers.
    lam_pow = F((-1)**sum(gamma),
                2**sum(c*g for c,g in zip(costs,gamma)))
    return prod(comb(b+g,b) for b,g in zip(beta,gamma))*lam_pow*alt


def run():
    count = 0
    row_counts = 0
    for weights in ((2,), (2,3)):
        m = len(weights)
        height = 12
        costs = tuple(2*w-1 for w in weights)
        types = [a for a in product(*(range(height//w+1) for w in weights))
                 if sum(w*x for w,x in zip(weights,a)) <= height]
        rows = [(s,beta) for s in range(height)
                for beta in types
                if s+2*sum(w*b for w,b in zip(weights,beta)) < height]
        for s,beta in rows:
            gap = height-s-2*sum(w*b for w,b in zip(weights,beta))
            row_sum = F(0)
            for alpha in types:
                gamma = tuple(a-b for a,b in zip(alpha,beta))
                if min(gamma)<0 or sum(gamma)>s:
                    continue
                if sum(c*g for c,g in zip(costs,gamma))<gap:
                    continue
                coeff = coefficient(beta,gamma,costs,gap)
                # Stronger degree-specific estimate underlying 3^{H/w_*}.
                bound = F(3**sum(alpha),2**gap)
                assert abs(coeff) <= bound
                row_sum += abs(coeff)
                count += 1
            assert row_sum <= len(types)*F(3**(height//min(weights)),2**gap)
            row_counts += 1

    compatible = []
    points = [QI(), QI(F(2,3),F(-1,7)), QI(0,F(601,100))]
    for n in range(1,7):
        for r in points:
            rows = []
            for j in range(2):
                for s in range(n):
                    y = [QI(comb(h,s) if s<=h else 0) for h in range(n)]
                    xy = []
                    for h in range(n):
                        # Direct coefficient extraction of
                        # (1+t)^h (j*r + log(1+t)).
                        value = j*r*(comb(h,s) if s<=h else 0)
                        value += sum(F((-1)**(ell+1),ell)*comb(h,s-ell)
                                     for ell in range(1,s+1)
                                     if s-ell<=h)
                        xy.append(value)
                    rows.append(y+xy)
            assert determinant(rows)==r**n
            compatible.append({'n':n,'r_real':str(r.re),'r_imag':str(r.im),
                               'determinant_equals_r_power_n':True})
    result = {'status':'passed','outside_coefficient_checks':count,
              'row_sum_checks':row_counts,
              'compatible_column_cases':compatible,
              'scope':'Finite exact identities, not a high-dimensional determinant bound.'}
    out = Path(__file__).with_name('review_verification.json')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','outside_coefficient_checks':count,
                      'row_sum_checks':row_counts,
                      'compatible_column_checks':len(compatible)}))


if __name__=='__main__':
    run()
