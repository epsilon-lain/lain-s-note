"""Independent reconstruction of signed Apéry moment-filter experiments.

Sources use binomial sums and a telescoping companion, polynomial coefficients
use closed multinomial sums, and root CFs come from strict dyadic intervals.
No import of moment_filter.py or either source constructor.
"""
from fractions import Fraction as F
from math import comb, gcd, isqrt
from pathlib import Path
import hashlib
import importlib.util
import json
import sys

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parent
helper_path = ROOT.parent/'pi-joint-selection-2026-11-10'/'verify.py'
spec = importlib.util.spec_from_file_location('independent_joint_helpers', helper_path)
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)


def coefficients(k):
    # Coefficients of (1+11z-z^2)^k, without polynomial-power convolution.
    return [sum(comb(k, j-c)*comb(j-c, c)*11**(j-2*c)*(-1)**c
                for c in range(max(0, j-k), j//2+1))
            for j in range(2*k+1)]


def bound(n, k, bits):
    a, b = n+2*k, 11*(n+k)
    D = b*b+4*a*n
    scale = 1 << bits
    m = isqrt(D*scale*scale)
    low, high = F(m-b*scale, 2*a*scale), F(m+1-b*scale, 2*a*scale)
    # Certify that these rationals enclose the maximizer by its polynomial.
    assert 0 < low < high < F(1, 11)
    assert a*low*low+b*low-n <= 0 < a*high*high+b*high-n
    H = lambda t: 1-11*t-t*t
    assert H(low) > 0
    sharp = high**n*H(low)**k
    simple = F(n**n*(k**k if k else 1), 11**n*(n+k)**(n+k))
    return min(sharp, simple)


def digest(record):
    return hashlib.sha256(json.dumps(record, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def decimal_contains(display, lo, hi):
    dl, dh = map(F, display)
    assert dl <= lo < hi <= dh


def run():
    path = ROOT/'moment_filter_results.json'
    data = json.loads(path.read_text())
    assert data['schema'] == 2
    par = data['parameters']
    limit = par['order_max']
    assert par['norm_budgets'] == ['min(25,B)', 'floor(sqrt(B))']
    us, rs = helper.sources(limit)
    vs = [r*u/6 for r,u in zip(rs,us)]
    polys = {k: coefficients(k) for k in range((limit-1)//2+1)}
    # Independent, tighter pi scoring precision than the constructor.
    score_bits = data['summary']['scoring_bits']+512
    pl, ph, terms = helper.pi_interval(score_bits)
    pi2l, pi2h = pl*pl, ph*ph
    records, qualities, expected_hits = [], [], []
    cursor = 0
    for n in range(1, limit+1):
        for k in range((limit-n)//2+1):
            saved_n, saved_k, saved_hash = data['cases'][cursor]
            assert (saved_n, saved_k) == (n,k)
            weights = polys[k]
            U = sum(c*us[n+j] for j,c in enumerate(weights))
            V = sum((c*vs[n+j] for j,c in enumerate(weights)), F(0))
            assert U != 0 and isinstance(U,int)
            aa, bb = 6*V.numerator, V.denominator*U
            if bb < 0:
                aa, bb = -aa, -bb
            common = gcd(aa,bb)
            A, B = aa//common, bb//common
            assert gcd(A,B) == 1 and B > 0
            M = bound(n,k,par['max_bound_dyadic_bits'])
            error = 12*M/abs(U)
            sign = (-1 if n%2 else 1)*(1 if U>0 else -1)
            w = {'n':n,'k':k,'order':n+2*k,'U':str(U),'V':str(V),
                 'raw_A':str(aa),'raw_B':str(bb),'gcd':str(common),
                 'A':str(A),'B':str(B),'M':str(M),
                 'source_error_bound':str(error),'error_sign':sign}
            costs = [B.bit_length(),abs(U).bit_length(),common.bit_length(),
                     max(abs(c).bit_length() for c in weights)]
            count, hits = 0,0
            R = F(A,B)
            if 9 < R < 16:
                budgets = {'fixed25':min(25,B),'square_root':isqrt(B)}
                caps = {name:isqrt(H*error.denominator//(B*error.numerator))
                        for name,H in budgets.items()}
                for p,q in helper.root_convergents(A,B,max(caps.values())):
                    if not 3*q <= p <= 4*q:
                        continue
                    count += 1
                    norm = B*p*p-A*q*q
                    for name,H in budgets.items():
                        if q > caps[name] or abs(norm) > H:
                            continue
                        assert gcd(p,q)==1 and B*q*q*error <= H
                        base = F(norm,B)
                        signed = max(abs(base),abs(base-sign*q*q*error))/6
                        triangle = (abs(base)+q*q*error)/6
                        assert signed <= triangle
                        expected_hits.append({'case_index':cursor,'n':n,'k':k,'budget':name,
                            'p':str(p),'q':str(q),'norm':str(norm),'H':str(H),
                            'q_limit':str(caps[name]),'signed_quality_upper':str(signed),
                            'triangle_quality_upper':str(triangle),'source_witness':w})
                        hits += 1
            rec = {'witness':w,'costs':costs,'root_convergents_checked':count,'budget_hits':hits}
            assert digest(rec) == saved_hash, (n,k,'source record mismatch')
            if sign > 0:
                el, eh = pi2l-R, pi2h-R
            else:
                el, eh = R-pi2h,R-pi2l
            assert 0 < el < eh < error
            qualities.append((B*B*el,B*B*eh))
            records.append(rec)
            cursor += 1
        if n%32==0:
            print(json.dumps({'verified_n':n,'records':cursor}),flush=True)
    assert cursor == len(data['cases'])
    assert len(expected_hits)==len(data['hits'])
    fractions = set()
    for expected, saved in zip(expected_hits,data['hits']):
        for key,value in expected.items():
            assert saved[key]==value,(key,'hit mismatch')
        p,q = int(saved['p']),int(saved['q'])
        ql,qh = helper.quality_interval(p,q,pl,ph)
        sl,sh = map(F,saved['actual_quality_interval'])
        assert sl <= ql < qh <= sh < F(saved['signed_quality_upper'])
        decimal_contains(saved['actual_quality_decimal'],ql,qh)
        fractions.add((p,q))
    selected = data['selected_source_witnesses']
    for saved in selected:
        i=saved['case_index']
        for key,value in records[i].items():
            assert saved[key]==value,(i,key)
        ql,qh=qualities[i]
        decimal_contains(saved['actual_source_quality_decimal'],ql,qh)
    byorder={r['witness']['n']:i for i,r in enumerate(records) if r['witness']['k']==0}
    best=min(range(cursor),key=lambda i:qualities[i][1])
    filtered=[i for i,r in enumerate(records) if r['witness']['k']>0]
    bestf=min(filtered,key=lambda i:qualities[i][1])
    def better(i,j):
        if qualities[i][1]<qualities[j][0]:
            return True
        assert qualities[j][1]<qualities[i][0] or records[i]['witness']['A']==records[j]['witness']['A'] and records[i]['witness']['B']==records[j]['witness']['B'], (i,j,'rank unresolved')
        return False
    improve_order=[i for i in filtered if better(i,byorder[records[i]['witness']['order']])]
    improve_n=[i for i in filtered if better(i,byorder[records[i]['witness']['n']])]
    summary=data['summary']
    assert summary['cases']==cursor
    assert summary['filtered_cases']==len(filtered)
    assert summary['root_convergents_checked']==sum(r['root_convergents_checked'] for r in records)
    assert summary['budget_hits']==len(expected_hits)
    assert summary['distinct_fractions']==[f'{p}/{q}' for p,q in sorted(fractions,key=lambda t:(t[1],t[0]))]
    assert summary['best_source_case']==best and summary['best_filtered_source_case']==bestf
    assert summary['filtered_improve_same_highest_input_order']==len(improve_order)
    assert summary['filtered_improve_same_n']==len(improve_n)
    for saved in data['same_highest_order_comparisons']:
        order=saved['order']
        baseline=byorder[order]
        eligible=[i for i in filtered if records[i]['witness']['order']==order]
        chosen=min(eligible,key=lambda i:qualities[i][1])
        assert saved['baseline_case_index']==baseline and saved['filtered_case_index']==chosen
        assert saved['filtered_n']==records[chosen]['witness']['n']
        assert saved['filtered_k']==records[chosen]['witness']['k']
        assert saved['baseline_B_bits']==records[baseline]['costs'][0]
        assert saved['filtered_B_bits']==records[chosen]['costs'][0]
        bl,bh=qualities[baseline]
        fl,fh=qualities[chosen]
        decimal_contains(saved['baseline_quality_decimal'],bl,bh)
        decimal_contains(saved['filtered_quality_decimal'],fl,fh)
        decimal_contains(saved['gain_decimal'],bl/fh,bh/fl)
    best_hit=min(data['hits'],key=lambda hit:helper.quality_interval(int(hit['p']),int(hit['q']),pl,ph)[1])
    assert summary['best_actual_pi_fraction']==best_hit['p']+'/'+best_hit['q']
    ql,qh=helper.quality_interval(int(best_hit['p']),int(best_hit['q']),pl,ph)
    decimal_contains(summary['best_actual_pi_quality'],ql,qh)
    report={'checks':'PASS','file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'source_records_reconstructed':cursor,'filtered_records':len(filtered),
        'independent_root_convergents_checked':summary['root_convergents_checked'],
        'certified_hits':len(expected_hits),'distinct_fractions':summary['distinct_fractions'],
        'source_error_signs_and_bounds_checked':cursor,
        'better_same_highest_order':len(improve_order),'better_same_n':len(improve_n),
        'same_highest_order_comparisons_checked':len(data['same_highest_order_comparisons']),
        'independent_pi_bits':score_bits,'independent_taylor_terms':terms,
        'method':'Binomial and telescoping sources; closed multinomial filter coefficients; dyadic CF; direct integer gcd; outward Taylor score.',
        'scope':'Finite exact reconstruction; source moment identity reviewed separately; no infinite-success or pi-record theorem.'}
    (ROOT/'moment_filter_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    run()
