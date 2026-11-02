"""Select saved period designs using the proved rational upper bound only.

This does not read the post-construction pi quality intervals. Candidate
generation itself uses integer/rational kernel data, not a target fraction.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import reduce
from math import gcd
import json


def main():
    root = Path(__file__).resolve().parent
    selected = []
    for name in ('hilbert_quotient_lift.json','hilbert_forced_period_lift.json'):
        data = json.loads((root/'period-main-kernel'/name).read_text())
        for block in data['per_m']:
            row = min(block['candidates'],
                      key=lambda r: F(r['joint_quality_squared_upper']))
            alpha = int(row['alpha'])
            rr = [int(x) for x in row['R_coefficients']]
            hh = [int(x) for x in row['H_coefficients']]
            z = sum(rr)
            candidate = F(-4*z,2**(block['c']-1)*alpha)
            assert candidate == F(int(row['p']),int(row['q']))
            bound = F(row['joint_quality_squared_upper'])
            factor = block['forced_period_factor']
            assert alpha % factor == z % factor == 0
            selected.append(dict(source=name,m=block['m'],
                                 p=str(candidate.numerator),q=str(candidate.denominator),
                                 forced_period_factor=factor,
                                 pi_period=str(alpha),boundary_period=str(z),
                                 primitive_coordinate_content=reduce(gcd,map(abs,rr+[alpha])),
                                 multiplier_content=reduce(gcd,map(abs,hh)),
                                 common_period_factor=gcd(abs(alpha),abs(z)),
                                 joint_quality_squared_upper=str(bound),
                                 certified_quality_below_one=bound<1))
    # Post-selection regression checks; these do not choose a candidate.
    example = next(r for r in selected
                   if r['source']=='hilbert_forced_period_lift.json' and r['m']==2)
    assert example['p']=='355' and example['q']=='113'
    assert example['forced_period_factor']==3
    assert example['common_period_factor']==3
    assert example['primitive_coordinate_content']==example['multiplier_content']==1
    assert 20*F(example['joint_quality_squared_upper'])<1
    result = dict(status='passed',selection='Minimum exact rational quality-squared upper bound; no pi intervals read.',
                  per_scale=selected,
                  exact_benchmark_certificate='Q^2 < 217375834098688/4366350858338469 < 1/20',
                  scope='A bounded certified construction benchmark; no new approximation record or infinite-sequence quality claim.')
    (root/'certificate_selected_periods.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
