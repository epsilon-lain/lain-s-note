"""Variable-base rational-approximation experiment; Python standard library only.

All floor decisions use exact rational enclosures, never float rounding.
The output is a finite certificate, not a proof of infinitely many hits.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from math import isqrt, gcd
from pathlib import Path
import json
import random


def atan_bounds(d, terms):
    total = sum((F((-1)**k, (2*k+1)*d**(2*k+1))
                 for k in range(terms)), F(0))
    following = F((-1)**terms, (2*terms+1)*d**(2*terms+1))
    return min(total, total+following), max(total, total+following)


def pi_bounds():
    a, b = atan_bounds(5, 160)
    c, d = atan_bounds(239, 50)
    return 16*a-4*d, 16*b-4*c


def rounded_enclosure(lo, hi, places=200):
    scale = 10**places
    lower = lo.numerator*scale // lo.denominator
    upper = -(-hi.numerator*scale // hi.denominator)
    return lower, upper, scale


def certified_digits(lower, upper, scale, base, length=3):
    """Both interval endpoints must belong to the same length-digit cell."""
    factor = base**length
    prefix = lower*factor // scale
    if prefix != upper*factor // scale:
        raise ArithmeticError('Insufficient enclosure precision')
    result = []
    for _ in range(length):
        result.append(prefix % base)
        prefix //= base
    assert prefix == 0
    return list(reversed(result))


def decimal_string(value, places=16):
    with localcontext() as context:
        context.prec = places
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def test_exact_identities():
    rng = random.Random(2848)
    for _ in range(12000):
        base = rng.randrange(2, 1000)
        denominator = rng.randrange(2, 10000)
        alpha = F(rng.randrange(1, denominator), denominator)
        q = base-1
        first = (base*alpha).__floor__()
        actual = q*q*(alpha-F(first,q))
        for length in (1,2,3,4):
            power = base**length
            prefix = (power*alpha).__floor__()
            tail = power*alpha-prefix
            repeated_prefix = first*(power-1)//q
            difference = prefix-repeated_prefix
            factor = F(q*q,power)
            reconstructed = factor*(difference+tail-F(first,q))
            assert actual == reconstructed
            lower = factor*(difference-F(first,q))
            upper = factor*(difference+1-F(first,q))
            assert lower <= actual < upper
            assert factor*max(0,abs(difference)-1) <= abs(actual)
            assert abs(actual) < factor*(abs(difference)+1)
    return 48000


def analyze_pi(limit=1000000):
    lo, hi = pi_bounds()
    lower, upper, scale = rounded_enclosure(lo,hi)
    alpha_lower, alpha_upper = lower-3*scale, upper-3*scale
    repeated2, repeated3, record_rows, repeated_rows = 0, [], [], []
    best = None
    for base in range(4,limit+1):
        q = base-1
        d1,d2,d3 = certified_digits(alpha_lower,alpha_upper,scale,base)
        if d1 == d2:
            repeated2 += 1
            factor=F(q*q,base**3)
            cert_lo=factor*(d3-d1-F(d1,q))
            cert_hi=factor*(d3-d1+1-F(d1,q))
            actual_lo=F(q*(q*alpha_lower-d1*scale),scale)
            actual_hi=F(q*(q*alpha_upper-d1*scale),scale)
            assert cert_lo <= actual_lo <= actual_hi < cert_hi
            repeated_rows.append({'base':base,'q':q,'digits':[d1,d2,d3],
                'signed_three_digit_bound_exact':[str(cert_lo),str(cert_hi)]})
            if d2 == d3:
                repeated3.append(base)
        # q times the distance of q*pi to its nearest integer.
        nearest = (2*q*lower+scale)//(2*scale)
        assert nearest == (2*q*upper+scale)//(2*scale)
        gap_lower, gap_upper = q*lower-nearest*scale,q*upper-nearest*scale
        assert gap_lower*gap_upper > 0
        candidate_lower = q*min(abs(gap_lower),abs(gap_upper))
        candidate_upper = q*max(abs(gap_lower),abs(gap_upper))
        if best is None or candidate_upper < best[0]:
            best = candidate_lower,candidate_upper
            record_rows.append({
                'base':base,'q':q,'numerator':nearest,
                'reduced_denominator':q//gcd(nearest,q),
                'digits':[d1,d2,d3],
                'normalized_error_lower':decimal_string(F(candidate_lower,scale),25),
                'normalized_error_upper':decimal_string(F(candidate_upper,scale),25),
                'normalized_error_exact_interval':{
                    'lower_numerator':str(candidate_lower),
                    'upper_numerator':str(candidate_upper),'denominator':str(scale)},
                'relative_third_gap':decimal_string(F(abs(d3-d1),base)),
            })
        elif best[0] <= candidate_upper and candidate_lower <= best[1]:
            raise ArithmeticError('Unresolved comparison of record errors')
    # Certify selected convergents beyond the exhaustive scan.
    cf_lo,cf_hi=F(lower,scale),F(upper,scale)
    p_prev2,p_prev,q_prev2,q_prev=0,1,1,0
    convergents=[]
    for index in range(40):
        a=cf_lo.__floor__()
        assert a == cf_hi.__floor__()
        p=a*p_prev+p_prev2
        q=a*q_prev+q_prev2
        if q>2:
            base=q+1
            digits=certified_digits(alpha_lower,alpha_upper,scale,base)
            error1=F(q*(q*lower-p*scale),scale)
            error2=F(q*(q*upper-p*scale),scale)
            assert error1*error2 > 0
            convergents.append({'index':index,'numerator':p,'q':q,'base':base,
                'digits':digits,'normalized_error':decimal_string(abs((error1+error2)/2),20),
                'relative_third_gap':decimal_string(F(abs(digits[2]-digits[0]),base),20)})
        cf_lo,cf_hi=1/(cf_hi-a),1/(cf_lo-a)
        p_prev2,p_prev,q_prev2,q_prev=p_prev,p,q_prev,q
    return {'base_min':4,'base_max':limit,'tested_bases':limit-3,
            'pi_enclosure_places':200,'pi_enclosure_width_units':upper-lower,
            'repeated_two_count':repeated2,'repeated_three_bases':repeated3,
            'repeated_two_rows':repeated_rows,
            'records':record_rows,'convergents':convergents}


def sqrt2_control():
    scale=10**200
    root=isqrt(2*scale*scale)
    lower,upper=root-scale,root+1-scale
    # sqrt(2)-1=[0;2,2,...].
    p_prev2,p_prev,q_prev2,q_prev=0,1,1,0
    rows=[]
    for index in range(21):
        a=0 if index==0 else 2
        p,q=a*p_prev+p_prev2,a*q_prev+q_prev2
        if q>=5:
            digits=certified_digits(lower,upper,scale,q+1)
            assert digits[0]==digits[1]
            error=abs(F(q*(q*(lower+upper)-2*p*scale),2*scale))
            assert error>F(1,4)
            rows.append({'q':q,'base':q+1,'digits':digits,
                'normalized_error':decimal_string(error,20),
                'relative_third_gap':decimal_string(F(abs(digits[2]-digits[0]),q+1),20)})
        p_prev2,p_prev,q_prev2,q_prev=p_prev,p,q_prev,q
    return rows


if __name__ == '__main__':
    result={'identity_checks':test_exact_identities(),
            'pi':analyze_pi(),'sqrt2_control':sqrt2_control()}
    destination=Path(__file__).with_name('results.json')
    destination.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('Exact identity checks:',result['identity_checks'])
    print('Pi bases tested:',result['pi']['tested_bases'])
    print('First two equal:',result['pi']['repeated_two_count'])
    print('First three equal:',result['pi']['repeated_three_bases'])
    print('First two equal bases:',[row['base'] for row in result['pi']['repeated_two_rows']])
    print('Record rows:')
    for row in result['pi']['records']:
        print({key:value for key,value in row.items() if key!='normalized_error_exact_interval'})
    print('Selected further convergents:')
    for row in result['pi']['convergents']:
        if row['digits'][0]==row['digits'][1] and float(row['normalized_error'])<0.1:
            print(row)
    print('Sqrt2 controls:',result['sqrt2_control'][:4])
