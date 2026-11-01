"""Exact finite experiments accelerating the original Leibniz tail.

No floating-point arithmetic is used for constructions or comparisons.
This is a finite search, not a claim about all possible accelerations.
Dependency: search_residuals.py, which supplies a certified Machin interval.
"""
from fractions import Fraction
from math import comb, gcd
from pathlib import Path
import argparse
import hashlib
import json
import sys

from search_residuals import pi_interval, scientific_decimal_bounds


def leibniz_prefixes(last):
    values = [Fraction(0)]
    for n in range(last):
        values.append(values[-1] + Fraction(4 * (-1)**n, 2*n+1))
    return values


def shifted_euler(prefix, n, depth):
    """S_n + 4*(-1)^n times an Euler-transformed alternating tail.

    Its jth unsigned correction is
        2*j! / product_{r=0}^j (2*n+2*r+1).
    Successive corrections have ratio (j+1)/(2*n+2*j+3).
    """
    value = prefix
    term = Fraction(2, 2*n+1)
    for m in range(1, depth+1):
        value += (-1)**n * term
        yield m, value
        term *= Fraction(m, 2*n+2*m+1)


def shanks_rows(prefixes, order):
    """Wynn epsilon recurrence; even rows are Shanks transforms.

    epsilon[-1,n]=0; epsilon[0,n]=S_n;
    epsilon[k+1,n]=epsilon[k-1,n+1]
                   +1/(epsilon[k,n+1]-epsilon[k,n]).
    A zero divisor aborts rather than silently inventing a value.
    """
    previous = [Fraction(0)] * len(prefixes)
    current = prefixes
    for level in range(1, 2*order+1):
        following = []
        for n in range(len(current)-1):
            difference = current[n+1] - current[n]
            if not difference:
                raise ZeroDivisionError((level, n))
            following.append(previous[n+1] + 1/difference)
        previous, current = current, following
        if level % 2 == 0:
            yield level//2, current


def determinant(matrix):
    a = [row[:] for row in matrix]
    sign = 1
    result = Fraction(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            sign *= -1
        result *= a[i][i]
        for j in range(i+1, len(a)):
            ratio = a[j][i]/a[i][i]
            for k in range(i+1, len(a)):
                a[j][k] -= ratio*a[i][k]
    return sign*result


def check_construction_identities():
    # Independent finite-difference formula checks the Euler recurrence.
    for n in range(8):
        term = Fraction(2, 2*n+1)
        for j in range(10):
            delta = sum((Fraction((-1)**k*comb(j,k), 2*n+2*k+1)
                         for k in range(j+1)), Fraction(0))
            assert term == 4*delta/Fraction(2**(j+1))
            term *= Fraction(j+1, 2*n+2*j+3)
    # Shanks's determinant definition checks a different implementation.
    values = leibniz_prefixes(16)
    difference = [b-a for a,b in zip(values,values[1:])]
    for order, row in shanks_rows(values, 4):
        for n in range(5):
            tail = [[difference[n+i+j] for j in range(order+1)]
                    for i in range(order)]
            numerator = determinant([values[n:n+order+1]] + tail)
            denominator = determinant([[Fraction(1)]*(order+1)] + tail)
            assert row[n] == numerator/denominator
    assert next(shanks_rows(values,1))[1][0] == 3
    return {'euler_independent_finite_difference_checks':80,
            'shanks_independent_determinant_checks':20}


def original_tail_denominator_checks(count=1000):
    values = leibniz_prefixes(2*count)
    for N in range(1, count+1):
        q = values[2*N].denominator
        D = 4*N-1
        for prime in (3,5):
            power = 1
            while power*prime <= D:
                power *= prime
            assert q % power == 0 and (q//power) % prime != 0
        assert 15*q > D*D
        assert 3*q >= D+2
    return {'checked_N':count,
            'identities':'v_3(q_N)=floor(log_3(4N-1)); v_5(q_N)=floor(log_5(4N-1))',
            'scope':'Finite checks support, but do not replace, the separate elementary proof.'}


def run(max_shift=100, euler_depth=128, shanks_order=12, digits=1500):
    checks = check_construction_identities()
    checks['original_tail_denominators'] = original_tail_denominator_checks()
    L,U,M,counts = pi_interval(digits)
    prefixes = leibniz_prefixes(max_shift+2*shanks_order)
    digest = hashlib.sha256()
    examined = {'euler':0,'shanks':0,'baseline':0}
    low_quality = {'euler':[],'shanks':[]}
    max_q_digits = {'euler':0,'shanks':0,'baseline':0}
    summaries = {'euler_by_shift':[], 'shanks_by_order':[]}
    samples = []
    all_best = {}
    positive_prefix_best = {}

    def candidate(family, parameters, value):
        p,q = value.numerator,value.denominator
        assert gcd(p,q) == 1
        lo,hi = q*L-p*M,q*U-p*M
        if lo > 0:
            alo,ahi = lo,hi
        elif hi < 0:
            alo,ahi = -hi,-lo
        else:
            raise RuntimeError('Increase pi precision to separate residual signs')
        klo,khi = q*alo,q*ahi
        if (khi-klo)*10**24 >= M:
            raise RuntimeError('Increase precision for certified quality comparisons')
        if klo <= M <= khi:
            raise RuntimeError('K<1 is unresolved')
        if 2*klo <= M <= 2*khi:
            raise RuntimeError('K<1/2 is unresolved')
        row = {'family':family, **parameters, 'p':str(p),'q':str(q),
               'q_digits':len(str(q)), 'beta_sign':1 if lo>0 else -1,
               'absolute_error_interval':scientific_decimal_bounds(alo,ahi,M*q),
               'signed_beta_interval':scientific_decimal_bounds(lo,hi,M),
               'K_interval':scientific_decimal_bounds(klo,khi,M),
               'K_less_than_one':khi<M, 'K_less_than_half':2*khi<M}
        examined[family] += 1
        max_q_digits[family] = max(max_q_digits[family],len(str(q)))
        digest.update(json.dumps([family,parameters,str(p),str(q)],sort_keys=True).encode())
        digest.update(b'\n')
        if family != 'baseline' and khi < M:
            low_quality[family].append(row)
        return {'row':row,'q':q,'alo':alo,'ahi':ahi,'klo':klo,'khi':khi}

    def better(first, second, metric='K'):
        if first is None:
            return second
        if first['row']['p']==second['row']['p'] and first['q']==second['q']:
            return first
        if metric == 'K':
            if second['khi'] < first['klo']:
                return second
            if second['klo'] > first['khi']:
                return first
        else:
            if second['ahi']*first['q'] < first['alo']*second['q']:
                return second
            if second['alo']*first['q'] > first['ahi']*second['q']:
                return first
        raise RuntimeError('An ordering is unresolved; increase pi precision')

    for n in range(max_shift+1):
        best_K = best_error = None
        for m,value in shifted_euler(prefixes[n],n,euler_depth):
            item = candidate('euler',{'prefix_terms':n,'tail_terms':m},value)
            best_K = better(best_K,item)
            best_error = better(best_error,item,'error')
            all_best['euler'] = better(all_best.get('euler'),item)
            if n>=2 and n%2==0:
                positive_prefix_best['euler'] = better(positive_prefix_best.get('euler'),item)
            if n in (0,2,10,100) and m in (1,2,4,8,16,32,64,128):
                samples.append(item['row'])
        summaries['euler_by_shift'].append({'prefix_terms':n,
                                            'best_K':best_K['row'],
                                            'best_absolute_error':best_error['row']})

    for order,values in shanks_rows(prefixes,shanks_order):
        best_K = best_error = None
        for n,value in enumerate(values[:max_shift+1]):
            item = candidate('shanks',{'prefix_terms':n,'order':order},value)
            best_K = better(best_K,item)
            best_error = better(best_error,item,'error')
            all_best['shanks'] = better(all_best.get('shanks'),item)
            if n>=2 and n%2==0:
                positive_prefix_best['shanks'] = better(positive_prefix_best.get('shanks'),item)
            if n in (0,2,10,100) and order in (1,2,3,4,8,12):
                samples.append(item['row'])
        summaries['shanks_by_order'].append({'order':order,
                                             'best_K':best_K['row'],
                                             'best_absolute_error':best_error['row']})

    baselines = []
    for n in range(2,max_shift+1,2):
        item = candidate('baseline',{'prefix_terms':n,'original_N':n//2},prefixes[n])
        baselines.append(item['row'])

    return {'status':'passed',
            'scope':'Exact finite scans of two fixed acceleration families; no universal failure or infinite-subsequence claim.',
            'parameters':{'max_shift':max_shift,'euler_depth':euler_depth,
                          'shanks_order':shanks_order,'pi_scale_digits':digits},
            'construction':{
                'base':'S_n=4 sum_{k=0}^{n-1} (-1)^k/(2k+1); original R_N=S_{2N}',
                'euler':'E_{n,m}=S_n+(-1)^n sum_{j=0}^{m-1} 2*j!/product_{r=0}^j(2*n+2*r+1)',
                'shanks':'Wynn epsilon_{2r,n} constructed from S_n,...,S_{n+2r}',
                'quality':'beta=q*pi-p; K=q*abs(beta)=q^2*abs(pi-p/q); p/q always reduced'},
            'pi_certificate':{'scaled_lower':str(L),'scaled_upper':str(U),
                              'denominator_power_of_ten':digits,'width_units':U-L,
                              'machin_arctan_counts':counts},
            'checks':checks,'examined':examined,'max_q_digits':max_q_digits,
            'all_constructed_fraction_sha256':digest.hexdigest(),
            'best_K_all':{k:v['row'] for k,v in all_best.items()},
            'best_K_from_nonempty_positive_original_prefix':{k:v['row'] for k,v in positive_prefix_best.items()},
            'all_candidates_with_K_less_than_one':low_quality,
            **summaries,'selected_samples':samples,'original_positive_tail_baselines':baselines,
            'data_retention':'Every candidate was checked. This JSON retains group extrema, every K<1 candidate, selected samples, and baseline rows; rerun code reproduces the entire finite scan.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-shift',type=int,default=100)
    parser.add_argument('--euler-depth',type=int,default=128)
    parser.add_argument('--shanks-order',type=int,default=12)
    parser.add_argument('--digits',type=int,default=1500)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('acceleration_results.json'))
    args = parser.parse_args()
    if min(args.max_shift,args.euler_depth,args.shanks_order,args.digits)<=0:
        parser.error('all numerical parameters must be positive')
    if hasattr(sys,'set_int_max_str_digits'):
        sys.set_int_max_str_digits(max(10000,3*args.digits))
    result = run(args.max_shift,args.euler_depth,args.shanks_order,args.digits)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','parameters','examined','max_q_digits',
                    'best_K_all','best_K_from_nonempty_positive_original_prefix',
                    'all_candidates_with_K_less_than_one')},ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
