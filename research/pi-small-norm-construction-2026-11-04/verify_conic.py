"""Independent small-norm checks using binomial sums and dyadic root intervals.

No import of conic_attempt.py. Its three-term source recurrence is replaced by
the closed binomial sum plus a telescoping Casoratian sum. Its quadratic-surd
CF recurrence is replaced by exact Euclidean extraction from dyadic intervals.
"""
from fractions import Fraction as F
from math import comb,gcd,isqrt
from pathlib import Path
import hashlib,json,sys,time
from verify import atan_taylor

sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parent

def root_convergents_interval(A,B,qmax):
    D=A*B;s=isqrt(D)
    pm2,pm1,qm2,qm1=0,1,1,0
    if s*s==D:
        num,den=s,B
        while den:
            a,rem=divmod(num,den)
            p,q=a*pm1+pm2,a*qm1+qm2
            if q>qmax:return
            yield p,q
            pm2,pm1,qm2,qm1=pm1,p,qm1,q
            num,den=den,rem
        return
    bits=4*qmax.bit_length()+2*B.bit_length()+128
    scale=1<<bits
    lower=isqrt((A*scale*scale)//B)
    lo,hi=F(lower,scale),F(lower+1,scale)
    assert lo*lo<F(A,B)<hi*hi
    while True:
        a=lo.numerator//lo.denominator
        au=-((-hi.numerator)//hi.denominator)-1
        assert a==au,'Increase exact dyadic precision; no branch may be guessed.'
        p,q=a*pm1+pm2,a*qm1+qm2
        if q>qmax:return
        yield p,q
        pm2,pm1,qm2,qm1=pm1,p,qm1,q
        assert lo>a
        lo,hi=1/(hi-a),1/(lo-a)

def verify():
    start=time.monotonic();path=ROOT/'conic_results.json'
    data=json.loads(path.read_text());nmax=data['parameters']['n_max']
    u=[sum(comb(n,j)**2*comb(n+j,j) for j in range(n+1)) for n in range(nmax+2)]
    v=[F(0)];ratio=F(0)
    for n in range(1,nmax+2):
        ratio+=F(5*(-1)**(n-1),n*n*u[n]*u[n-1])
        v.append(u[n]*ratio)
    assert (u[:2],v[:2])==([1,3],[F(0),F(5)])
    widths=[]
    for n,record in enumerate(data['sources'],1):
        assert record['n']==n and F(record['u'])==u[n] and F(record['v'])==v[n]
        r0,r1=F(6*v[n],u[n]),F(6*v[n+1],u[n+1])
        assert record['pi_squared_endpoints']==list(map(str,sorted([r0,r1])))
        delta=F(30,(n+1)**2*u[n]*u[n+1])
        assert abs(r0-r1)==delta
        widths.append(delta)
    bits=max(d.denominator.bit_length()-d.numerator.bit_length()+256 for d in widths)
    l5,h5,terms5=atan_taylor(5,bits);l239,h239,terms239=atan_taylor(239,bits)
    pl,pu=16*l5-4*h239,16*h5-4*l239
    assert 3<pl<pu<4
    for record in data['sources']:
        lo,hi=map(F,record['pi_squared_endpoints'])
        assert lo<pl*pl<pu*pu<hi
    expected=[];checked=0
    for case in data['cases']:
        n,k=case['n'],case['mix'];delta=widths[n-1]
        R=6*(v[n]+k*v[n+1])/(u[n]+k*u[n+1]);A,B=R.numerator,R.denominator
        budgets={'fixed25':min(25,B),'square_root':isqrt(B)}
        limits={name:isqrt((H*delta.denominator)//(B*delta.numerator)) for name,H in budgets.items()}
        total=0;best=None;hit_count=0
        for p,q in root_convergents_interval(A,B,max(limits.values())):
            if not 3*q<=p<=4*q:continue
            total+=1;checked+=1
            norm=B*p*p-A*q*q
            best=abs(norm) if best is None else min(best,abs(norm))
            for name,H in budgets.items():
                if q<=limits[name] and abs(norm)<=H:
                    expected.append((n,k,name,A,B,p,q,norm,H,limits[name],delta))
                    hit_count+=1
        assert case['A_digits']==len(str(A)) and case['B_digits']==len(str(B))
        assert case['convergents_checked']==total and case['hit_count_with_budgets']==hit_count
        assert case['min_absolute_norm']==(str(best) if best is not None else None)
        assert case['q_limit_digits']=={name:len(str(qmax)) for name,qmax in limits.items()}
    assert len(expected)==len(data['hits'])
    distinct=set()
    for values,hit in zip(expected,data['hits']):
        n,k,name,A,B,p,q,norm,H,qlim,delta=values
        assert (n,k,name)==(hit['n'],hit['mix'],hit['budget'])
        assert [A,B,p,q,norm,H,qlim]==[int(hit[z]) for z in ['A','B','p','q','norm','H','q_limit']]
        assert gcd(p,q)==1 and abs(norm)<=H<=B and q*q*delta<=F(H,B)
        assert F(hit['source_width'])==delta
        bound=(F(abs(norm),B)+q*q*delta)/6
        assert bound==F(hit['quality_upper_from_conic'])<=F(H,3*B)
        lo,hi=map(F,hit['quality_interval']);dl,du=map(F,hit['quality_decimal_enclosure'])
        assert dl<=lo<hi<=du and hi<bound
        if F(p,q)<pl:tl,tu=q*(q*pl-p),q*(q*pu-p)
        else:
            assert F(p,q)>pu
            tl,tu=q*(p-q*pu),q*(p-q*pl)
        assert lo<=tl<tu<=hi
        assert hit['q_digits']==len(str(q))
        distinct.add((p,q))
    assert data['summary']['conics_tried']==len(data['cases'])==nmax*len(data['parameters']['mixes'])
    assert data['summary']['root_convergents_checked']==checked
    assert data['summary']['hits_with_budgets']==len(expected)
    assert data['summary']['max_hit_n']==max(h['n'] for h in data['hits'])
    assert {(int(d['p']),int(d['q'])) for d in data['summary']['distinct_fractions']}==distinct
    report={'checks':'PASS','file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'binomial_source_values':len(u),'telescoping_companion_values':len(v),
            'pi_squared_brackets_checked':nmax,'conics_rechecked':len(data['cases']),
            'independent_root_convergents_checked':checked,'certified_hits':len(expected),
            'distinct_fractions':[f'{p}/{q}' for p,q in sorted(distinct,key=lambda t:t[1])],
            'max_hit_n':max(h['n'] for h in data['hits']),
            'independent_pi_grid_bits':bits,'taylor_terms':[terms5,terms239],
            'no_new_fraction_better_than_355_over_113':distinct=={(3,1),(22,7),(355,113)},
            'scope':'Finite experiment only. Neither an infinite successful family nor an exclusion of larger or modified families.',
            'elapsed_seconds_display':round(time.monotonic()-start,6)}
    (ROOT/'conic_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':verify()
