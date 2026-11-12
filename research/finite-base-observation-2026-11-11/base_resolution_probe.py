"""Independent finite base-grid / continued-fraction consistency probe.

Pure integer arithmetic. No imported experiment code, files, pi digits, or API.
This checks finite instances of the Legendre-based theorem in FIXED_DEPTH.md.
"""
from math import gcd
import argparse
import json
from pathlib import Path

if not __debug__:
    raise RuntimeError('Assertions required')

S = 10**600

def atan_bounds(x, terms):
    lo = hi = 0
    for j in range(terms):
        floor = S // ((2*j+1)*x**(2*j+1))
        if j % 2:
            lo -= floor+1; hi -= floor
        else:
            lo += floor; hi += floor+1
    den = (2*terms+1)*x**(2*terms+1)
    tail = (S+den-1)//den
    return lo-tail, hi+tail

a,b=atan_bounds(5,450); c,d=atan_bounds(239,150)
PL,PU=16*a-4*d,16*b-4*c
assert 3*S < PL < PU < 7*S//2

def convergents_through(Q):
    ln,ld,un,ud=PL,S,PU,S
    pm2,pm1,qm2,qm1=0,1,1,0
    result=[]
    for n in range(2000):
        term=ln//ld
        assert term==un//ud, 'More pi precision needed'
        p,q=term*pm1+pm2,term*qm1+qm2
        result.append((n,p,q))
        if q>Q:
            return result
        lr,ur=ln-term*ld,un-term*ud
        assert lr>0 and ur>0
        ln,ld,un,ud=ud,ur,ld,lr
        pm2,pm1,qm2,qm1=pm1,p,qm1,q
    raise AssertionError('CF cap reached')

def quality(p,q):
    # Endpoints of S*q^2*abs(pi-p/q).
    lo,hi=q*PL-p*S,q*PU-p*S
    if lo>=0:return q*lo,q*hi
    if hi<=0:return -q*hi,-q*lo
    return 0,q*max(-lo,hi)

def certified_winner(candidates):
    values=[(quality(p,q),p,q) for p,q in candidates]
    best=min(values,key=lambda x:x[0][1])
    for val,p,q in values:
        if (p,q)!=(best[1],best[2]):assert best[0][1]<val[0]
    assert best[0][1]<S//2
    return best[1],best[2]

def run(b,k,brute=False):
    M=b**k
    allcf=convergents_through(M)
    selected=[(p,q) for n,p,q in allcf if q<=M and M%q==0]
    winner=certified_winner(selected)
    if brute:
        grid=[]
        for m in range(3*M,4*M+1):
            g=gcd(m,M)
            grid.append((m//g,M//g))
        assert certified_winner(grid)==winner
    lo,hi=quality(*winner)
    return dict(base=b,max_digits=k,denominator_upper_bound=str(M),
                first_certified_denominator_above_limit=str(allcf[-1][2]),
                cf_terms_including_first_denominator_above_limit=len(allcf),
                admissible_convergents=len(selected),minimizer_p=str(winner[0]),minimizer_q=str(winner[1]),
                admissible_convergent_list=[{'index':n,'p':str(p),'q':str(q)} for n,p,q in allcf if q<=M and M%q==0],
                quality_lower_numerator=str(lo),quality_upper_numerator=str(hi),scale_decimal_power=600,
                brute_grid_comparison=brute)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='Optional full certificate report JSON path')
    args=parser.parse_args()
    cases=[run(b,k,True) for b,k in [(2,8),(3,6),(7,3),(10,3),(11,3),(113,1)]]
    cases += [run(b,100) for b in (2,7,10)]
    result={'format':'external-base-resolution-probe/1','scope':'Finite grids only, not infinite-base-power or pi bad-approximability theorem',
            'uploaded_code_executed':False,'provider_requests':0,'pi_interval_width_scaled':PU-PL,'cases':cases}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
