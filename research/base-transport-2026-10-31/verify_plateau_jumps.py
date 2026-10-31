"""Independent exact checks of survival-window denominator jumps.

No stored pi digits or continued-fraction coefficients are inputs.
All certification decisions use integers or Fraction. Decimal output is display only.
"""
from __future__ import annotations
from fractions import Fraction as F
from decimal import Decimal, localcontext
from pathlib import Path
from math import gcd
import json
import hashlib

HERE = Path(__file__).resolve().parent


def atan_interval(inv: int, terms: int) -> tuple[F, F]:
    partial = sum((F((-1)**j, (2*j+1)*inv**(2*j+1)) for j in range(terms)), F(0))
    nxt = F((-1)**terms, (2*terms+1)*inv**(2*terms+1))
    return min(partial, partial+nxt), max(partial, partial+nxt)


def ceil_rat(x: F) -> int:
    return -((-x.numerator)//x.denominator)


def certify_digits(count: int = 256, bits: int = 1024):
    a, b = atan_interval(5, 260)
    c, d = atan_interval(239, 80)
    exact_lo, exact_hi = 16*a-4*d, 16*b-4*c
    scale = 1 << bits
    pl = (exact_lo*scale).__floor__()
    ph = ceil_rat(exact_hi*scale)
    assert ph-pl == 1
    lo = hi = scale
    digits, remainders = [], [(F(1), F(1))]
    for _ in range(count):
        lo = pl*lo//scale
        hi = (ph*hi+scale-1)//scale
        digit = lo//scale
        assert digit == hi//scale
        lo -= digit*scale
        hi -= digit*scale
        assert 0 < lo <= hi < scale and 0 <= digit <= 3
        digits.append(digit)
        remainders.append((F(lo,scale), F(hi,scale)))
    return (F(pl,scale), F(ph,scale)), digits, remainders


def homogeneous(p: int, q: int, digits: list[int], depth: int):
    numerator = den = 1
    for digit in digits[:depth]:
        den *= q
        numerator = p*numerator-digit*den
    return numerator, den


class Windows:
    def __init__(self, digits: list[int]):
        self.digits = digits
        self.queries = 0

    def side(self, p: int, q: int, depth: int) -> int:
        self.queries += 1
        numerator = den = 1
        for digit in self.digits[:depth]:
            den *= q
            numerator = p*numerator-digit*den
            # Escape is absorbing on [3,4], so early exit preserves the side.
            if numerator <= -2*den:
                return -1
            if numerator >= 2*den:
                return 1
        return 0

    def least(self, depth: int) -> dict:
        a,b,c,d = 3,1,4,1
        while True:
            sign = self.side(a+c,b+d,depth)
            if sign == 0:
                p,q = a+c,b+d
                assert b*c-a*d == 1
                assert self.side(a,b,depth) == -1
                assert self.side(c,d,depth) == 1
                numerator, den = homogeneous(p,q,self.digits,depth)
                assert abs(numerator) < 2*den and gcd(p,q) == 1
                ln, ld = homogeneous(a,b,self.digits,depth)
                rn, rd = homogeneous(c,d,self.digits,depth)
                assert ln < -2*ld and rn > 2*rd
                return {'N':depth,'p':p,'q':q,'left':[a,b],'right':[c,d]}
            def query(k: int) -> int:
                if sign < 0:
                    return self.side(a+k*c,b+k*d,depth)
                return self.side(c+k*a,d+k*b,depth)
            low, high = 1,2
            while query(high) == sign:
                low, high = high,2*high
            while high-low > 1:
                mid = (high+low)//2
                if query(mid) == sign:
                    low = mid
                else:
                    high = mid
            if sign < 0:
                a,b = a+low*c,b+low*d
            else:
                c,d = c+low*a,d+low*b

    def endpoint(self, depth: int, value: int, bits: int = 500) -> tuple[F,F]:
        left,right = F(3),F(4)
        for _ in range(bits):
            mid = (left+right)/2
            num,den = homogeneous(mid.numerator,mid.denominator,self.digits,depth)
            if num < value*den:
                left=mid
            else:
                right=mid
        assert right-left == F(1,2**bits)
        return left,right


def absolute_bounds(lo: F, hi: F) -> tuple[F,F]:
    if lo > 0: return lo,hi
    if hi < 0: return -hi,-lo
    raise AssertionError('Insufficient precision to separate rational from pi')


def display(x: F) -> str:
    with localcontext() as ctx:
        ctx.prec = 26
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def interval(lo:F,hi:F) -> dict:
    assert lo <= hi
    return {'exact':[str(lo),str(hi)],'display':[display(lo),display(hi)]}


def main():
    (pl,ph),digits,remainders = certify_digits()
    win=Windows(digits)
    rows=[win.least(n) for n in range(4,161)]
    lookup={r['N']:r for r in rows}
    expected={4:F(19,6),5:F(22,7),8:F(245,78),9:F(311,99),10:F(355,113),14:F(355,113),15:F(53938,17169)}
    for n,v in expected.items():
        assert F(lookup[n]['p'],lookup[n]['q']) == v
    platforms=[]
    for row in rows:
        if platforms and platforms[-1]['p']==row['p'] and platforms[-1]['q']==row['q']:
            platforms[-1]['end']=row['N']
        else:
            platforms.append({'start':row['N'],'end':row['N'],'p':row['p'],'q':row['q']})
    good=[]
    jumps=[]
    for old,new in zip(platforms,platforms[1:]):
        e=new['start']; p,q=old['p'],old['q']; u,v=new['p'],new['q']
        rat=F(p,q)
        error=absolute_bounds(rat-ph,rat-pl)
        kl,kh=q*q*error[0],q*q*error[1]
        jump=F(v,q)
        assert jump > 1
        # The universal lower bound needs no small-kappa hypothesis.
        assert 29*kl*jump > 1
        rec={'exit_depth':e,'old_fraction':str(rat),'new_fraction':str(F(u,v)),
             'plateau_length':old['end']-old['start']+1,'jump':str(jump),
             'jump_display':display(jump),'kappa':interval(kl,kh),
             'cross_determinant':u*q-p*v}
        jumps.append(rec)
        if kh <= F(1,64):
            assert jump < 1/kh+1
            assert abs(u*q-p*v)==1
            above = rat > ph
            parent=lookup[e-1]['left' if above else 'right']
            a,b=parent
            assert v>=b and (v-b)%q==0
            k=(v-b)//q
            assert k>=1 and u==k*p+a
            # Isolate the endpoint on the far side of pi from the old fraction.
            low,high=win.endpoint(e,-2 if above else 2)
            dl,dh=(rat-high,rat-low) if above else (low-rat,high-rat)
            assert dl>0
            k_low=((1/(q*dh)-b)/q).__floor__()+1
            k_high=((1/(q*dl)-b)/q).__floor__()+1
            assert k_low==k_high==k
            # Exact rounding defect theta is strictly between zero and one.
            theta_low=jump-1/(q*q*dl)
            theta_high=jump-1/(q*q*dh)
            assert 0<theta_low<=theta_high<1
            # k-1 lies beyond the far endpoint; k is inside the new window.
            assert win.side((k-1)*p+a,(k-1)*q+b,e)==(-1 if above else 1)
            assert win.side(u,v,e)==0
            # Record exact geometric coefficients; finite checks do not establish an asymptotic rate.
            rlo,rhi=remainders[e]
            rec=dict(rec,parent=parent,ray_index=k,
                     far_endpoint=interval(low,high),
                     far_distance=interval(dl,dh),
                     rounding_defect=interval(theta_low,theta_high),
                     kappa_times_jump=interval(kl*jump,kh*jump),
                     exact_geometric_coefficient=interval(error[0]/dh,error[1]/dl),
                     native_remainder_at_exit=interval(rlo,rhi))
            good.append(rec)
    out={'status':'all_exact_checks_passed',
         'input':'No stored pi digits or continued fractions; Machin identity with rational alternating bounds.',
         'native_digits':len(digits),'precision_bits':1024,
         'digit_prefix':''.join(map(str,digits[:60])),
         'pi_interval':interval(pl,ph),'depth_range':[4,160],
         'layers':len(rows),'platform_count':len(platforms),
         'complete_platform_transitions':len(jumps),
         'maximum_observed_platform_length':max(p['end']-p['start']+1 for p in platforms),
         'all_jump_lower_bounds_passed':True,
         'small_kappa_threshold':'1/64','small_kappa_transition_count':len(good),
         'small_kappa_exact_successor_certificates':good,
         'max_jump':max(jumps,key=lambda x:F(x['jump'])),
         'integer_queries':win.queries,
         'scope':'Finite exact checks of proofs in PLATEAU_EXIT_JUMPS.md. No all-depth bound, no new pi approximation record, no P2848 conclusion.',
         'platforms':platforms,'transitions':jumps}
    (HERE/'plateau_jump_results.json').write_text(json.dumps(out,indent=2)+'\n')
    summary = {
        'status': out['status'],
        'scope': out['scope'],
        'source_snapshot': 'ed23900554caf17b68a0d5903cbd564839e513dd',
        'digit_count': len(digits),
        'digit_bytes_sha256': hashlib.sha256(bytes(digits)).hexdigest(),
        'arithmetic': 'Standard-library integers and Fraction; 1024-bit outward dyadic bounds',
        'pi_exact_interval': [str(pl),str(ph)],
        'depth_range': [4,160], 'layers': len(rows),
        'minimality_certificates': len(rows), 'platform_count': len(platforms),
        'complete_transitions': len(jumps),
        'universal_jump_lower_bound_checks': len(jumps),
        'determinant_one_transitions': sum(abs(r['cross_determinant'])==1 for r in jumps),
        'small_kappa_threshold': '1/64',
        'small_kappa_exact_successor_checks': len(good),
        'max_observed_plateau_length': out['maximum_observed_platform_length'],
        'integer_queries': win.queries,
        'selected_certificates': [{
            'exit_depth': r['exit_depth'],
            'old_fraction': r['old_fraction'], 'new_fraction': r['new_fraction'],
            'parent': r['parent'], 'ray_index': r['ray_index'],
            'plateau_length': r['plateau_length'], 'jump': r['jump'],
            'kappa_approx': r['kappa']['display'][0],
            'far_endpoint_exact_interval': r['far_endpoint']['exact'],
            'rounding_defect_approx': r['rounding_defect']['display'][0],
            'kappa_times_jump_approx': r['kappa_times_jump']['display'][0],
            'geometric_coefficient_approx': r['exact_geometric_coefficient']['display'][0]
        } for r in good],
        'display_warning': 'Fields ending in _approx are rounded, not interval certificates. Exact rational bounds are retained in plateau_jump_results.json and recomputed by this script.'
    }
    (HERE/'plateau_jump_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['pi_exact_interval','selected_certificates']},indent=2))
    for r in good:
        print('GOOD', r['old_fraction'], '->', r['new_fraction'], 'ray index',r['ray_index'])

if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Run without -O: assertions are the certificate checks.')
    main()
