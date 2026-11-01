"""Certified pi residual search using only Python integer arithmetic.

All approximation decisions use a Machin interval, never floating-point pi.
The output is finite evidence, not a proof of liminf q*|q*pi-p| = 0.
"""
from pathlib import Path
import argparse
import json
import math
import sys

def arctan_scaled_bounds(n, scale):
    """Return integers L,U with L/scale < atan(1/n) < U/scale."""
    lower = upper = 0
    power = n
    k = 0
    while True:
        denominator = (2*k+1)*power
        term_floor = scale//denominator
        if term_floor == 0:
            # The remaining alternating series has this sign and size < 1/scale.
            if k % 2:
                lower -= 1
            else:
                upper += 1
            return lower, upper, k
        if k % 2:
            lower -= term_floor+1
            upper -= term_floor
        else:
            lower += term_floor
            upper += term_floor+1
        k += 1
        power *= n*n

def pi_interval(digits):
    scale = 10**digits
    l5,u5,n5 = arctan_scaled_bounds(5,scale)
    l239,u239,n239 = arctan_scaled_bounds(239,scale)
    return 16*l5-4*u239,16*u5-4*l239,scale,{'atan5_terms':n5,'atan239_terms':n239}

def continued_fraction_prefix(lower, upper, denominator):
    ln,ld,un,ud=lower,denominator,upper,denominator
    coefficients=[]
    while True:
        a,b=ln//ld,un//ud
        if a != b:
            return coefficients
        coefficients.append(a)
        lr,ur=ln-a*ld,un-a*ud
        if lr == 0 or ur == 0:
            return coefficients
        ln,ld,un,ud=ud,ur,ld,lr

def fixed_decimal_bounds(lo,hi,denominator,places=24):
    scale=10**places
    ilo=lo*scale//denominator
    ihi=-((-hi*scale)//denominator)
    def render(integer):
        sign='-' if integer<0 else ''
        integer=abs(integer)
        whole,frac=divmod(integer,scale)
        return f'{sign}{whole}.{frac:0{places}d}'
    return [render(ilo),render(ihi)]

def scientific_decimal_bounds(lo,hi,denominator,places=24):
    """Outward-rounded scientific notation; all rounding uses integers."""
    if not (0 < lo <= hi or lo <= hi < 0):
        raise ValueError('Scientific display requires a strictly signed interval')
    absolute_nearest_zero=lo if lo>0 else -hi
    exponent=len(str(absolute_nearest_zero))-len(str(denominator))
    if exponent>=0:
        if absolute_nearest_zero < denominator*10**exponent:
            exponent-=1
    elif absolute_nearest_zero*10**(-exponent) < denominator:
        exponent-=1
    if exponent>=0:
        mantissas=fixed_decimal_bounds(lo,hi,denominator*10**exponent,places)
    else:
        mantissas=fixed_decimal_bounds(lo*10**(-exponent),hi*10**(-exponent),denominator,places)
    return [f'{bound}e{exponent:+d}' for bound in mantissas]

def residual_certificate(p,q,L,U,M):
    lo,hi=q*L-p*M,q*U-p*M
    if hi<0:
        sign=-1; alo,ahi=-hi,-lo
    elif lo>0:
        sign=1; alo,ahi=lo,hi
    else:
        raise ValueError('Residual sign not separated by supplied pi interval')
    return {'p':str(p),'q':str(q),'q_digits':len(str(q)),'beta_sign':sign,
        'beta_interval':fixed_decimal_bounds(lo,hi,M),
        'q_abs_beta_interval':fixed_decimal_bounds(q*alo,q*ahi,M)},q*alo,q*ahi

def scan(digits=4000, uncertainty_places=24, max_denominator_digits=None):
    if digits<=0 or uncertainty_places<=0:
        raise ValueError('digits and uncertainty_places must be positive')
    if max_denominator_digits is not None and max_denominator_digits<=0:
        raise ValueError('max_denominator_digits must be positive')
    L,U,M,counts=pi_interval(digits)
    coefficients=continued_fraction_prefix(L,U,M)
    prev2p,prevp,prev2q,prevq=0,1,1,0
    minlo=minhi=None
    records=[]
    positive_records=[]
    pos_minlo=pos_minhi=None
    quality_count=0
    last_q=0
    first_above=None
    max_q=None if max_denominator_digits is None else 10**max_denominator_digits
    for index,a in enumerate(coefficients):
        p,q=a*prevp+prev2p,a*prevq+prev2q
        assert math.gcd(p,q)==1
        if index:
            assert abs(p*prevq-prevp*q)==1
        prev2p,prevp,prev2q,prevq=prevp,p,prevq,q
        if max_q is not None and q>max_q:
            first_above={'convergent_index':index,'q':str(q),'q_digits':len(str(q))}
            break
        if q*q*(U-L)*10**uncertainty_places > M:
            continue
        item,lo,hi=residual_certificate(p,q,L,U,M)
        item.update({'convergent_index':index,
                     'next_partial_quotient':coefficients[index+1] if index+1<len(coefficients) else None})
        quality_count+=1
        last_q=q
        if minlo is None or hi<minlo:
            records.append(item)
        elif lo<minhi:
            raise RuntimeError('A possible quality record is unresolved; increase digits')
        minlo=lo if minlo is None else min(minlo,lo)
        minhi=hi if minhi is None else min(minhi,hi)
        if item['beta_sign']>0:
            if pos_minlo is None or hi<pos_minlo:
                positive_records.append(item)
            elif lo<pos_minhi:
                raise RuntimeError('A possible positive-side record is unresolved')
            pos_minlo=lo if pos_minlo is None else min(pos_minlo,lo)
            pos_minhi=hi if pos_minhi is None else min(pos_minhi,hi)
    if not records or not positive_records:
        raise RuntimeError('Pi interval too coarse to certify quality records')
    if max_q is not None:
        if first_above is None:
            raise RuntimeError('Pi interval insufficient to reach first convergent above Q')
        if quality_count != first_above['convergent_index']:
            raise RuntimeError('Some convergent quality below Q lacks requested precision')
        if 2*minhi >= M or 2*pos_minhi >= M:
            raise RuntimeError('Legendre coverage requires certified records below 1/2')
    # Each target is an exact signed combination of the original positive
    # Leibniz residuals gamma_1=3*pi-8 and gamma_2=105*pi-304.
    # The raw numerator/denominator both have gcd 72 for primitive p/q.
    bestp,bestq=int(records[-1]['p']),int(records[-1]['q'])
    seedu,seedv=304*bestq-105*bestp,3*bestp-8*bestq
    rawp,rawq=8*seedu+304*seedv,3*seedu+105*seedv
    assert rawp==72*bestp and rawq==72*bestq
    assert math.gcd(rawp,rawq)==72
    best_index=records[-1]['convergent_index']
    local=[]
    if best_index>=2:
        p2,p1,q2,q1=0,1,1,0
        for index,a in enumerate(coefficients[:best_index+1]):
            p,q=a*p1+p2,a*q1+q2
            p2,p1,q2,q1=p1,p,q1,q
            if index>=best_index-2:
                item=residual_certificate(p,q,L,U,M)[0]
                item['convergent_index']=index
                item['partial_quotient']=a
                local.append(item)
        a=coefficients[best_index]
        assert int(local[2]['p'])==a*int(local[1]['p'])+int(local[0]['p'])
        assert int(local[2]['q'])==a*int(local[1]['q'])+int(local[0]['q'])
    pairs=[(22,7),(333,106),(355,113),(103993,33102)]
    examples=[residual_certificate(p,q,L,U,M)[0] for p,q in pairs]
    assert 22+333==355 and 7+106==113
    assert 292*355+333==103993 and 292*113+106==33102
    for item in records+positive_records+examples+local:
        p,q=int(item['p']),int(item['q'])
        lo,hi=q*L-p*M,q*U-p*M
        item['beta_scientific_interval']=scientific_decimal_bounds(lo,hi,M)
        item['pi_minus_p_over_q_scientific_interval']=scientific_decimal_bounds(lo,hi,M*q)
    return {'status':'passed','method':'Exact integer Machin interval and interval-stable continued fractions',
        'decimal_scale_digits':digits,'pi_scaled_lower':str(L),'pi_scaled_upper':str(U),
        'pi_scale_denominator_exponent':digits,'pi_width_units':U-L,
        'arctan_counts':counts,'certified_cf_coefficients':coefficients,
        'certified_cf_length':len(coefficients),'quality_checked_convergents':quality_count,
        'quality_absolute_uncertainty_max':f'1e-{uncertainty_places}',
        'last_quality_checked_q':str(last_q),'last_quality_checked_q_digits':len(str(last_q)),
        'strict_quality_records':records,'positive_beta_quality_records':positive_records,
        'maximum_denominator_power_of_ten':max_denominator_digits,
        'first_convergent_above_denominator_limit':first_above,
        'finite_global_minimum_certificate':None if max_q is None else {
            'principle':'Legendre: abs(x-p/q)<1/(2q^2) for reduced p/q implies a continued-fraction convergent.',
            'record_less_than_half':records[-1]['q_abs_beta_interval'],
            'coverage':'Every convergent with denominator <= Q checked; first convergent above Q certified.',
            'conclusion':'The last strict record is the unique minimum of q*abs(q*pi-p) among reduced fractions with 1<=q<=Q.'},
        'best_local_cancellation_certificate':local,
        'best_leibniz_seed_certificate':{
            'gamma1':'3*pi-8','gamma2':'105*pi-304',
            'u':str(seedu),'v':str(seedv),'common_factor':72,
            'identity':'beta=(u*gamma1+v*gamma2)/72',
            'normalization':'8*u+304*v=72*p; 3*u+105*v=72*q; gcd(rawp,rawq)=72'},
        'cancellation_examples':examples,
        'scope':'Finite certified convergent scan. It does not prove an infinite subsequence with q*abs(beta) tending to zero.'}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--digits',type=int,default=4000)
    parser.add_argument('--max-denominator-digits',type=int)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('residual_search_results.json'))
    args=parser.parse_args()
    if args.digits<=0:
        parser.error('--digits must be positive')
    if args.max_denominator_digits is not None and args.max_denominator_digits<=0:
        parser.error('--max-denominator-digits must be positive')
    if hasattr(sys,'set_int_max_str_digits'):
        sys.set_int_max_str_digits(max(10000,3*args.digits))
    result=scan(args.digits,max_denominator_digits=args.max_denominator_digits)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','decimal_scale_digits','certified_cf_length',
        'quality_checked_convergents','last_quality_checked_q_digits']},indent=2))
    for row in result['strict_quality_records']:
        print('record',row['convergent_index'],'q_digits',row['q_digits'],
              'quality',row['q_abs_beta_interval'],'next_a',row['next_partial_quotient'])
    print('positive-side records:',[(r['convergent_index'],r['q_digits'],r['q_abs_beta_interval']) for r in result['positive_beta_quality_records']])

if __name__=='__main__':
    main()
