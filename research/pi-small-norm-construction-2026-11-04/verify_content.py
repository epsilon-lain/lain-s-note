"""Verify prescribed-content witnesses and all ensuing bounded conic searches.

Uses closed binomial and telescoping sums for sources and dyadic-interval CF
extraction for quadratic roots. Does not import either constructor.
"""
from fractions import Fraction as F
from math import comb,gcd,isqrt
from pathlib import Path
import hashlib,json,time
from verify_conic import root_convergents_interval
from verify import atan_taylor

ROOT=Path(__file__).resolve().parent

def run():
    start=time.monotonic();path=ROOT/'content_lift_results.json'
    data=json.loads(path.read_text());nmax=data['n_max']
    u=[sum(comb(n,j)**2*comb(n+j,j) for j in range(n+1)) for n in range(nmax+2)]
    v=[F(0)];ratio=F(0)
    for n in range(1,nmax+2):
        ratio+=F(5*(-1)**(n-1),n*n*u[n]*u[n-1]);v.append(u[n]*ratio)
    expected=[];total=0;max_gap_bits=0
    for case in data['cases']:
        n=case['n'];r0,r1=6*v[n]/u[n],6*v[n+1]/u[n+1]
        A0,B0,A1,B1=r0.numerator,r0.denominator,r1.numerator,r1.denominator
        w=case['content_witness'];G=int(w['modulus']);k=int(w['kernel_slope_numerator'])
        x,y,content=int(w['weight0']),int(w['weight1']),int(w['content'])
        assert G==abs(A0*B1-A1*B0)>0 and 0<=k<G
        assert (A0+k*A1)%G==(B0+k*B1)%G==0
        assert x>=0 and y>=0 and x+y>0 and gcd(x,y)==1
        aa,bb=x*A0+y*A1,x*B0+y*B1
        assert gcd(aa,bb)==content and G%content==0
        A,B=int(case['A']),int(case['B'])
        assert (A*content,B*content)==(aa,bb) and gcd(A,B)==1
        assert min(r0,r1)<=F(A,B)<=max(r0,r1)
        assert B<=min(B0,B1)
        delta=F(30,(n+1)**2*u[n]*u[n+1]);assert delta==abs(r0-r1)
        assert F(case['source_width'])==delta
        assert (int(case['source0_B']),int(case['source1_B']))==(B0,B1)
        assert F(case['height_gain_over_source0'])==F(B0,B)
        assert F(case['source0_B_squared_width'])==B0*B0*delta
        assert F(case['new_B_squared_width'])==B*B*delta
        max_gap_bits=max(max_gap_bits,delta.denominator.bit_length()-delta.numerator.bit_length())
        budgets={'fixed25':min(25,B),'square_root':isqrt(B)}
        limits={name:isqrt((H*delta.denominator)//(B*delta.numerator)) for name,H in budgets.items()}
        checked=0;best=None
        for p,q in root_convergents_interval(A,B,max(limits.values())):
            if not 3*q<=p<=4*q:continue
            checked+=1;total+=1;norm=B*p*p-A*q*q
            best=abs(norm) if best is None else min(best,abs(norm))
            for name,H in budgets.items():
                if q<=limits[name] and abs(norm)<=H:
                    expected.append((n,name,p,q,norm,H,limits[name],(F(abs(norm),B)+q*q*delta)/6))
        assert case['root_convergents_checked']==checked
        assert case['min_absolute_norm']==(str(best) if best is not None else None)
    bits=max_gap_bits+256
    l5,h5,_=atan_taylor(5,bits);l239,h239,_=atan_taylor(239,bits)
    lo,hi=16*l5-4*h239,16*h5-4*l239
    for n in range(1,nmax+1):
        r0,r1=6*v[n]/u[n],6*v[n+1]/u[n+1]
        assert min(r0,r1)<lo*lo<hi*hi<max(r0,r1)
    assert len(expected)==len(data['hits'])
    distinct=set()
    for values,hit in zip(expected,data['hits']):
        n,name,p,q,norm,H,qlim,bound=values
        assert (n,name)==(hit['n'],hit['budget'])
        assert [p,q,norm,H,qlim]==[int(hit[z]) for z in ['p','q','norm','H','q_limit']]
        assert F(hit['quality_upper_from_conic'])==bound and gcd(p,q)==1
        ql,qh=map(F,hit['quality_interval']);dl,du=map(F,hit['quality_decimal_enclosure'])
        assert dl<=ql<qh<=du and qh<bound
        if F(p,q)<lo:tl,tu=q*(q*lo-p),q*(q*hi-p)
        else:
            assert F(p,q)>hi
            tl,tu=q*(p-q*hi),q*(p-q*lo)
        assert ql<=tl<tu<=qh
        distinct.add((p,q))
    summary=data['summary']
    assert summary['sources_preconditioned']==len(data['cases'])==nmax
    assert summary['root_convergents_checked']==total
    assert summary['hits_with_budgets']==len(expected)
    assert summary['distinct_fractions']==[f'{p}/{q}' for p,q in sorted(distinct,key=lambda z:z[1])]
    assert summary['last_source_original_B_digits']==len(data['cases'][-1]['source0_B'])
    assert summary['last_source_new_B_digits']==len(data['cases'][-1]['B'])
    assert summary['max_hit_n']==max(h['n'] for h in data['hits'])
    report={'checks':'PASS','file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'content_witnesses_checked':len(data['cases']),'quadratic_candidates_rechecked':total,
            'certified_hits':len(expected),'distinct_fractions':summary['distinct_fractions'],
            'last_source_original_B_digits':summary['last_source_original_B_digits'],
            'last_source_new_B_digits':summary['last_source_new_B_digits'],
            'height_claim_includes_weight_cost':True,'max_hit_n':summary['max_hit_n'],
            'scope':'Finite positive-cone weight selection, no global minimal-height or infinite-success theorem.',
            'elapsed_seconds_display':round(time.monotonic()-start,6)}
    (ROOT/'content_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':run()
