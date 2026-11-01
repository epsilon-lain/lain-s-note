"""Exact finite illustration of the countermodel theorem.

The infinite classification is proved in the accompanying note, not inferred
from these finite checks. Only the standard library is used.
"""
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import sys


def convergents(coefficients):
    pp, p = 0, 1
    qq, q = 1, 0
    rows = []
    for a in coefficients:
        pp, p = p, a*p+pp
        qq, q = q, a*q+qq
        rows.append((p,q))
    return rows


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--pi-certificate',type=Path)
    args=parser.parse_args()
    sys.set_int_max_str_digits(0)
    # Taken from the already certified pi prefix in the parent project.
    prefix = [3,7,15,1,292,1,1,1,2,1,3,1,14]
    rows = convergents(prefix)
    p,q = rows[-1]
    r,s = rows[-2]
    assert abs(p*s-r*q)==1
    lower,upper = sorted((Fraction(p,q),Fraction(p+r,q+s)))
    assert upper-lower == Fraction(1,q*(q+s))
    decimals = 10
    assert lower*10**decimals//1 == upper*10**decimals//1
    Q=113
    profile=[]
    thresholds=[Fraction(1,2),Fraction(1,10),Fraction(1,100),Fraction(1,1000)]
    for denominator in range(1,Q+1):
        numerator=int(denominator*lower+Fraction(1,2))
        assert numerator==int(denominator*upper+Fraction(1,2))
        residual_lower=denominator*lower-numerator
        residual_upper=denominator*upper-numerator
        assert residual_lower*residual_upper>0
        kl,ku=sorted((denominator*abs(residual_lower),denominator*abs(residual_upper)))
        for threshold in thresholds:
            assert ku<threshold or kl>threshold
        profile.append((kl,ku,denominator,numerator))
    profile.sort()
    assert all(left[1]<right[0] for left,right in zip(profile,profile[1:]))
    tail=[]
    for k in range(1,257):
        j=k.bit_length()-1
        tail.append(j+1 if k>=2 and k==2**j else 1)
    all_rows=convergents(prefix+tail)
    spikes=[]
    N=len(prefix)-1
    for j in range(1,9):
        n=N+2**j
        P,Qn=all_rows[n-1]
        _,Qprev=all_rows[n-2]
        ratio=Fraction(Qprev,Qn)
        kl=1/(j+2+ratio)
        ku=1/(j+1+ratio)
        assert Fraction(1,j+3)<kl<ku<=Fraction(1,j+1)
        spikes.append({'j':j,'convergent_index':n-1,
                       'denominator':str(Qn),
                       'quality_strict_lower':str(kl),
                       'quality_strict_upper':str(ku)})
    result={'status':'passed','pi_prefix':prefix,
            'cylinder_lower':str(lower),'cylinder_upper':str(upper),
            'cylinder_width':str(upper-lower),
            'shared_decimal_prefix_10_places':str(int(lower*10**decimals)),
            'finite_profile_denominator_bound':Q,
            'all_nearest_numerators_fixed':True,
            'all_quality_rankings_fixed':True,
            'rational_thresholds':[str(x) for x in thresholds],
            'all_threshold_decisions_fixed':True,
            'best_visible_fraction':f'{profile[0][3]}/{profile[0][2]}',
            'sparse_tail_spikes':spikes,
            'scope':'Finite exact cylinder checks. Infinite BA/non-BA and mu=2 claims rely on the proof, not these checks.'}
    if args.pi_certificate:
        source_bytes=args.pi_certificate.read_bytes()
        previous=json.loads(source_bytes)
        previous_prefix=previous['certified_cf_coefficients']
        assert previous_prefix[:len(prefix)]==prefix
        assert previous['maximum_denominator_power_of_ten']==5000
        quality_lower=Fraction(previous['finite_global_minimum_certificate']['record_less_than_half'][0])
        assert quality_lower>Fraction(1,21000)
        large_rows=convergents(previous_prefix)
        _,last_q=large_rows[-1]
        _,prior_q=large_rows[-2]
        diameter_denominator=last_q*(last_q+prior_q)
        assert diameter_denominator>10**11994
        assert diameter_denominator>21000*10**10000
        assert 1000*(10**2498)**2<=10**5000
        result['existing_certificate_transfer']={
            'source_filename':args.pi_certificate.name,
            'source_sha256':hashlib.sha256(source_bytes).hexdigest(),
            'certified_prefix_length':len(previous_prefix),
            'maximum_prefix_partial_quotient':max(previous_prefix[1:]),
            'last_denominator_digits':len(str(last_q)),
            'cylinder_diameter_denominator_digits':len(str(diameter_denominator)),
            'strict_target_distance_bound':'10^-11994',
            'normalized_quality_error_bound_for_q_at_most_10pow5000':'10^-1994',
            'rational_side_comparisons_preserved_to_denominator':'10^5000',
            'nearest_grid_numerators_preserved_to_denominator':'10^5000/2',
            'entire_quality_order_and_four_thresholds_preserved_to_denominator':'10^2498',
            'dependency':'Uses the independently verified earlier pi certificate; does not re-prove that source certificate.'}
    out=Path(__file__).with_name('verification.json')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='sparse_tail_spikes'},indent=2))


if __name__=='__main__':
    main()
