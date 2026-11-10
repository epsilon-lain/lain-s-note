"""Signed integer filters of the Beukers/Apéry zeta(2) moment remainder.

Generation uses the recurrence and a rigorous scalar maximum bound only.
Pi enters after every source and bounded conic candidate has been generated.
No fixed decimal pi value or pi continued-fraction table is an input.
"""
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path
import argparse
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent
LEGACY = ROOT.parent / 'pi-small-norm-construction-2026-11-04'
sys.path.insert(0, str(LEGACY))
from conic_attempt import source, root_convergents
from construct import display_interval
from verify import atan_taylor

sys.set_int_max_str_digits(0)


def polynomial_powers(kmax):
    h = [1]
    ans = []
    for k in range(kmax + 1):
        ans.append(h)
        nxt = [0] * (len(h) + 2)
        for j, c in enumerate(h):
            nxt[j] += c
            nxt[j+1] -= 11*c
            nxt[j+2] -= c
        h = nxt
    return ans


def max_bound(n, k, bits=96):
    """Upper bound on max t^n(1-11t-t^2)^k, 0<=t<=mu.

    The maximum is at the positive root of
    (n+2k)t^2 + 11(n+k)t - n = 0. The power with exponent zero is 1.
    """
    a, b = n + 2*k, 11*(n+k)
    D = b*b + 4*a*n
    scale = 1 << bits
    root = isqrt(D * scale * scale)
    low = F(root - b*scale, 2*a*scale)
    high = F(root + 1 - b*scale, 2*a*scale)
    assert 0 < low <= high < F(1, 11)
    Hlow = 1 - 11*low - low*low
    assert Hlow > 0
    upper = high**n * Hlow**k
    simple = F(n**n * (k**k if k else 1), 11**n * (n+k)**(n+k))
    return min(upper, simple), upper, simple


def digest(record):
    return hashlib.sha256(json.dumps(record, sort_keys=True,
                         separators=(',', ':')).encode()).hexdigest()


def run(order_max=128):
    us, vs = source(order_max)
    polys = polynomial_powers((order_max-1)//2)
    raw, hits = [], []
    for n in range(1, order_max+1):
        for k in range((order_max-n)//2 + 1):
            h = polys[k]
            w = [(-1)**j * c for j, c in enumerate(h)]
            U = sum(wj*us[n+j] for j, wj in enumerate(w))
            V = sum((wj*vs[n+j] for j, wj in enumerate(w)), F(0))
            assert U.denominator == 1 and U != 0
            U = U.numerator
            aa, bb = 6*V.numerator, V.denominator*U
            if bb < 0:
                aa, bb = -aa, -bb
            common = gcd(aa, bb)
            A, B = aa//common, bb//common
            R = F(A, B)
            M, Msharp, Msimple = max_bound(n, k)
            error_bound = 12*M/abs(U)
            error_sign = (1 if n % 2 == 0 else -1) * (1 if U > 0 else -1)
            witness = {'n':n, 'k':k, 'order':n+2*k, 'U':str(U), 'V':str(V),
                       'raw_A':str(aa), 'raw_B':str(bb), 'gcd':str(common),
                       'A':str(A), 'B':str(B), 'M':str(M),
                       'source_error_bound':str(error_bound), 'error_sign':error_sign}
            row = {'n':n, 'k':k, 'order':n+2*k,
                   'B_bits':B.bit_length(), 'U_bits':abs(U).bit_length(),
                   'gcd_bits':common.bit_length(),
                   'max_weight_bits':max(abs(x).bit_length() for x in w),
                   'source_bound_decimal':display_interval(B*B*error_bound, B*B*error_bound),
                   'exact_witness_sha256':digest(witness),
                   'root_convergents_checked':0, 'budget_hits':0}
            if 9 < R < 16:
                budgets = {'fixed25':min(25,B), 'square_root':isqrt(B)}
                limits = {name:isqrt((H*error_bound.denominator)//(B*error_bound.numerator))
                          for name,H in budgets.items()}
                for p,q in root_convergents(A,B,max(limits.values())):
                    if not 3*q <= p <= 4*q:
                        continue
                    row['root_convergents_checked'] += 1
                    norm = B*p*p-A*q*q
                    for name,H in budgets.items():
                        if q > limits[name] or abs(norm) > H:
                            continue
                        base = F(norm,B)
                        # Exact one-sided residual range; cancellation is retained.
                        signed_bound = max(abs(base),abs(base-error_sign*q*q*error_bound))/6
                        triangle_bound = (abs(base)+q*q*error_bound)/6
                        assert signed_bound <= triangle_bound
                        hits.append({'case_index':len(raw),'n':n,'k':k,'budget':name,
                                     'p':str(p),'q':str(q),'norm':str(norm),
                                     'H':str(H),'q_limit':str(limits[name]),
                                     'signed_quality_upper':str(signed_bound),
                                     'triangle_quality_upper':str(triangle_bound),
                                     'source_witness':witness})
                        row['budget_hits'] += 1
            raw.append((row, witness))
        if n % 16 == 0:
            print(json.dumps({'n_completed':n,'cases':len(raw),'budget_hits':len(hits)}),flush=True)

    bits = 4*max(int(w['B']).bit_length() for _,w in raw)+256
    l5,h5,_ = atan_taylor(5,bits)
    l239,h239,_ = atan_taylor(239,bits)
    pl,ph = 16*l5-4*h239,16*h5-4*l239
    l2,h2 = pl*pl,ph*ph
    for row,w in raw:
        R = F(int(w['A']),int(w['B']))
        if w['error_sign'] > 0:
            lo,hi = l2-R,h2-R
        else:
            lo,hi = R-h2,R-l2
        assert 0 < lo < hi < F(w['source_error_bound'])
        B = int(w['B'])
        row['actual_source_quality_decimal'] = display_interval(B*B*lo,B*B*hi)
        row['actual_source_quality_interval'] = [str(B*B*lo),str(B*B*hi)]
    for hit in hits:
        p,q = int(hit['p']),int(hit['q'])
        if F(p,q) < pl:
            ql,qh = q*(q*pl-p),q*(q*ph-p)
        else:
            assert F(p,q) > ph
            ql,qh = q*(p-q*ph),q*(p-q*pl)
        assert 0 < ql < qh < F(hit['signed_quality_upper'])
        hit['actual_quality_interval'] = [str(ql),str(qh)]
        hit['actual_quality_decimal'] = display_interval(ql,qh)

    # Keep giant rational score intervals only for selected witnesses.
    best_source_index = min(range(len(raw)),key=lambda i:F(raw[i][0]['actual_source_quality_interval'][1]))
    best_filtered_index = min((i for i,(r,w) in enumerate(raw) if r['k']>0),
                              key=lambda i:F(raw[i][0]['actual_source_quality_interval'][1]))
    byorder = {row['n']:i for i,(row,w) in enumerate(raw) if row['k']==0}
    improved_over_same_order = []
    improved_over_n = []
    for i,(row,w) in enumerate(raw):
        if not row['k']:
            continue
        upper=F(row['actual_source_quality_interval'][1])
        if upper<F(raw[byorder[row['order']]][0]['actual_source_quality_interval'][0]):
            improved_over_same_order.append(i)
        if upper<F(raw[byorder[row['n']]][0]['actual_source_quality_interval'][0]):
            improved_over_n.append(i)
    selected = {best_source_index,best_filtered_index}
    comparisons = []
    for order in [8,16,32,64,128]:
        if order > order_max:
            continue
        base_i = byorder[order]
        filt_i = min((i for i,(row,w) in enumerate(raw)
                      if row['order']==order and row['k']>0),
                     key=lambda i:F(raw[i][0]['actual_source_quality_interval'][1]))
        base_interval = list(map(F,raw[base_i][0]['actual_source_quality_interval']))
        filt_interval = list(map(F,raw[filt_i][0]['actual_source_quality_interval']))
        selected.update([base_i,filt_i])
        comparisons.append({'order':order,'baseline_case_index':base_i,
                            'filtered_case_index':filt_i,
                            'filtered_n':raw[filt_i][0]['n'],'filtered_k':raw[filt_i][0]['k'],
                            'baseline_B_bits':raw[base_i][0]['B_bits'],
                            'filtered_B_bits':raw[filt_i][0]['B_bits'],
                            'baseline_quality_decimal':display_interval(*base_interval),
                            'filtered_quality_decimal':display_interval(*filt_interval),
                            'gain_decimal':display_interval(base_interval[0]/filt_interval[1],
                                                            base_interval[1]/filt_interval[0])})
    if improved_over_same_order:
        selected.add(min(improved_over_same_order,key=lambda i:F(raw[i][0]['actual_source_quality_interval'][1])))
    selected.update(hit['case_index'] for hit in hits)
    witness_records=[]
    for i in sorted(selected):
        row,w=raw[i]
        witness_records.append({'case_index':i,'witness':w,
                                'actual_source_quality_decimal':row['actual_source_quality_decimal'],
                                'costs':[row['B_bits'],row['U_bits'],row['gcd_bits'],row['max_weight_bits']],
                                'root_convergents_checked':row['root_convergents_checked'],
                                'budget_hits':row['budget_hits']})
    rows=[]
    for row,w in raw:
        exact={'witness':w,
               'costs':[row['B_bits'],row['U_bits'],row['gcd_bits'],row['max_weight_bits']],
               'root_convergents_checked':row['root_convergents_checked'],
               'budget_hits':row['budget_hits']}
        rows.append([row['n'],row['k'],digest(exact)])
    distinct=sorted({(int(h['p']),int(h['q'])) for h in hits},key=lambda t:(t[1],t[0]))
    besthit=min(hits,key=lambda h:F(h['actual_quality_interval'][1])) if hits else None
    summary={'cases':len(rows),'filtered_cases':sum(r[1]>0 for r in rows),
             'root_convergents_checked':sum(r['root_convergents_checked'] for r,w in raw),
             'budget_hits':len(hits),'distinct_fractions':[f'{p}/{q}' for p,q in distinct],
             'best_source_case':best_source_index,'best_filtered_source_case':best_filtered_index,
             'filtered_improve_same_highest_input_order':len(improved_over_same_order),
             'filtered_improve_same_n':len(improved_over_n),
             'best_actual_pi_quality':None if besthit is None else besthit['actual_quality_decimal'],
             'best_actual_pi_fraction':None if besthit is None else besthit['p']+'/'+besthit['q'],
             'scoring_bits':bits}
    out={'schema':2,'parameters':{'order_max':order_max,'n_min':1,'k_min':0,
                                 'max_bound_dyadic_bits':96,
                                 'norm_budgets':['min(25,B)','floor(sqrt(B))']},
         'construction':'H(t)=(1-11t-t^2)^k signed integer filters of Beukers moments',
         'case_columns':['n','k','exact_record_sha256'],
         'case_hash_record':'{witness, costs:[B_bits,U_bits,gcd_bits,max_weight_bits],root_convergents_checked,budget_hits}; canonical JSON sort_keys=True separators=(comma,colon)',
         'cases':rows,'selected_source_witnesses':witness_records,
         'same_highest_order_comparisons':comparisons,'hits':hits,'summary':summary,
         'scope':'Finite signed-remainder source construction with actual gcd and final pi conic checks. No new record or infinite-success claim.'}
    (ROOT/'moment_filter_results.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(json.dumps(summary,indent=2),flush=True)
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--order-max',type=int,default=128)
    args=parser.parse_args()
    assert args.order_max >= 3
    run(args.order_max)
