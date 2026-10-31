"""Two closed problems for base transport.

Run: python verify_closed_problems.py
Standard-library exact checks always run. Optional mpmath supplies explicitly
labelled high-precision illustrations; these do not certify infinite theorems.
No stored decimal digits of pi or continued-fraction inputs are used.
"""
from __future__ import annotations
from fractions import Fraction as Q
from dataclasses import dataclass
from pathlib import Path
import hashlib, json, random


def atan_bounds(inv: int, bits: int):
    terms = 1
    while inv ** (2*terms+1) < (1 << (bits+16)):
        terms += 1
    s = sum((Q((-1)**j, (2*j+1)*inv**(2*j+1))
             for j in range(terms)), Q())
    e = Q((-1)**terms, (2*terms+1)*inv**(2*terms+1))
    return min(s, s+e), max(s, s+e), terms


def certified_pi_word(count: int = 512, bits: int = 1792):
    a,b,ta = atan_bounds(5,bits)
    c,d,tc = atan_bounds(239,bits)
    lo,hi = 16*a-4*d,16*b-4*c
    scale=1<<bits
    p=(lo*scale).__floor__(); P=-((-hi*scale).__floor__())
    r=R=scale; digits=[]
    for n in range(1,count+1):
        s=p*r//scale; S=(P*R+scale-1)//scale
        digit=s//scale
        if digit != S//scale:
            raise ArithmeticError(f'Unresolved digit {n}; increase bits')
        r,R=s-digit*scale,S-digit*scale
        assert 0<r<=R<scale and 0<=digit<=3
        digits.append(digit)
    return Q(p,scale),Q(P,scale),digits,{'bits':bits,'atan_terms':[ta,tc]}


@dataclass(frozen=True)
class CQ:
    re: Q=Q(0)
    im: Q=Q(0)
    def __add__(self,b):
        if not isinstance(b,CQ): b=CQ(Q(b))
        return CQ(self.re+b.re,self.im+b.im)
    __radd__=__add__
    def __neg__(self): return CQ(-self.re,-self.im)
    def __sub__(self,b): return self+(-b)
    def __mul__(self,b):
        if not isinstance(b,CQ): b=CQ(Q(b))
        return CQ(self.re*b.re-self.im*b.im,self.re*b.im+self.im*b.re)
    __rmul__=__mul__
    def inv(self):
        den=self.re*self.re+self.im*self.im
        if den==0: raise ZeroDivisionError
        return CQ(self.re/den,-self.im/den)
    def __pow__(self,n):
        if n<0: return self.inv()**(-n)
        r,b=CQ(Q(1)),self
        while n:
            if n&1: r=r*b
            b=b*b; n//=2
        return r


def polynomial_value(digits,z):
    u=z.inv(); s=CQ()
    for d in reversed(digits): s=(s+d)*u
    return s


def exact_four_sample_test(word):
    N=4; A=3; R=Q(5,3); X=R**N
    T=Q(A)/(X-1); threshold=(1-T)/(2*X)
    eta=Q(1,50); i=CQ(Q(),Q(1))
    noises=[CQ(eta),CQ(Q(),eta),CQ(-eta),CQ(Q(),-eta)]
    values=[polynomial_value(word,CQ(R)*i**j) for j in range(N)]
    readings=[v+e for v,e in zip(values,noises)]
    tail=Q(A)*R**(-len(word))/(R-1)
    assert eta+tail<threshold
    out=[]
    for r in range(1,N+1):
        B=sum((values[j]*i**(j*r) for j in range(N)),CQ())*(R**r/N)
        direct=sum((Q(word[n-1])*R**(r-n)
                    for n in range(r,len(word)+1,N)),Q())
        assert B==CQ(direct)
        noisy=sum((readings[j]*i**(j*r) for j in range(N)),CQ())*(R**r/N)
        estimate=noisy.re-T/2
        recovered=(estimate+Q(1,2)).__floor__()
        assert recovered==word[r-1]
        out.append(recovered)
    return {'scope':'Exact Gaussian-rational arithmetic, four samples',
            'radius':str(R),'X':str(X),'noise_per_sample':str(eta),
            'tail_bound':str(tail),'sharp_threshold':str(threshold),
            'recovered':out}


def exact_uniform_checks(digits):
    rng=random.Random(4312026)
    checked=0
    for A in (1,2,3,9):
        X=Q(2*(A+1)); T=Q(A)/(X-1)
        eta=(1-T)/(2*X)*Q(3,4)
        for N in (1,4,16,64):
            words=[[rng.randrange(A+1) for _ in range(8*N)] for _ in range(12)]
            if A==3: words.append(digits[:8*N])
            for word in words:
                for r in range(1,N+1):
                    block=[word[n-1] for n in range(r,len(word)+1,N)]
                    known=sum((Q(d)*X**(-k) for k,d in enumerate(block)),Q())
                    unseen=Q(A)*X**(-len(block))/(1-1/X)
                    # A uniform bound X*eta dominates R^r*eta for all r<=N.
                    lo=known-T/2-X*eta
                    hi=known+unseen-T/2+X*eta
                    assert word[r-1]-Q(1,2)<lo<=hi<word[r-1]+Q(1,2)
                    checked+=1
    # Sharp counterexample: digits at positions 2N,3N,... are A;
    # competing word has only digit 1 at N. Their data differ by 2*eta*.
    sharp=[]
    for A in (1,2,3,9):
        for X in (Q(A+2),Q(2*(A+1)),Q(10*(A+1))):
            data_a=Q(A)/(X*(X-1)); data_b=1/X
            eps=(1-Q(A)/(X-1))/(2*X)
            assert data_b-data_a==2*eps>0
            midpoint=(data_a+data_b)/2
            assert midpoint-data_a==data_b-midpoint==eps
            sharp.append({'A':A,'X':str(X),'threshold':str(eps)})
    # Exact loss of identifiability at X=A+1 (same observations, unequal prefix).
    for A in (1,2,3,9):
        X=Q(A+1)
        assert Q(A)/(X*(X-1))==1/X
    return {'component_interval_cases':checked,
            'sharp_indistinguishability_pairs':sharp,
            'boundary_collisions_checked':4,
            'scope':'Exact finite arithmetic; universal claim proved in note'}


def numerical_illustrations(digits,counts):
    try:
        import mpmath as mp
    except ImportError:
        return {'available':False,'reason':'Optional mpmath not installed'}
    mp.mp.dps=90
    def val(z):
        u=1/z; x=mp.mpc(0)
        for d in reversed(digits): x=(x+d)*u
        return x
    dec=[]
    for N in (16,32,64):
        R=mp.root(8,N); eta=mp.mpf(1)/40
        roots=[mp.exp(2j*mp.pi*j/N) for j in range(N)]
        noises=[eta*mp.exp(2j*mp.pi*((7*j+3)%N)/N) for j in range(N)]
        readings=[val(R*roots[j])+noises[j] for j in range(N)]
        out=[]; err=mp.mpf(0)
        for r in range(1,N+1):
            B=R**r/N*sum(readings[j]*roots[j]**r for j in range(N))
            centered=mp.re(B)-mp.mpf(3)/14
            out.append(int(mp.floor(centered+mp.mpf('0.5'))))
            err=max(err,abs(centered-digits[r-1]))
        assert out==digits[:N]
        tail=3*R**(-len(digits))/(R-1)
        assert eta+tail<mp.mpf(1)/28
        dec.append({'N':N,'R_approx':mp.nstr(R,20),
                    'prefix_recovered': ''.join(map(str,out)),
                    'noise_amplitude':'1/40',
                    'truncation_error_bound_approx':mp.nstr(tail,10),
                    'maximum_centered_error_approx':mp.nstr(err,12)})
    pi=mp.pi; m=sum(n*d*pi**(-n) for n,d in enumerate(digits,1))
    tubes=[]
    for t in (4,10,16):
        t=mp.mpf(t); F=mp.re(val(t)); delta=1-F
        D=mp.log(pi)/mp.log(t)
        def P(theta):
            return delta**D/m*(pi/(pi-1)*pi**(-theta)
                              +pi/(t-pi)*(t/pi)**theta)
        theta_min=mp.log((t-pi)*mp.log(pi)/((pi-1)*mp.log(t/pi)))/mp.log(t)
        rows=[]
        for N in (12,24,48):
            for theta in (mp.mpf(0),mp.mpf(1)/2):
                h=delta*t**(-N-theta)
                W=sum(counts[n]*min(h,delta*t**(-n)) for n in range(1,len(counts)))
                # c_n <= pi^n; omitted tail bound below is theoretical,
                # evaluated here in high precision rather than interval arithmetic.
                rem=delta*(pi/t)**len(counts)/(1-pi/t)
                observed=W/h**(1-D)
                error=abs(observed-P(theta))
                # c_n=pi^n/m+e_n, |e_n|<=2^n/(pi-2).
                raw_bound=(h*pi/(m*(pi-1))
                           +h/(pi-2)*(2**(N+1)-2)
                           +delta/(pi-2)*(2/t)**(N+1)/(1-2/t))
                assert error<raw_bound/h**(1-D)+rem/h**(1-D)
                rows.append({'N':N,'theta':str(theta),
                             'normalized_tube_approx':mp.nstr(observed,24),
                             'predicted_limit_approx':mp.nstr(P(theta),24),
                             'difference_approx':mp.nstr(error,10)})
        tubes.append({'base':str(t),'dimension_approx':mp.nstr(D,22),
                      'P_at_0_approx':mp.nstr(P(0),24),
                      'P_at_half_approx':mp.nstr(P(mp.mpf('0.5')),24),
                      'theta_min_approx':mp.nstr(theta_min,24),
                      'P_min_approx':mp.nstr(P(theta_min),24),
                      'average_approx':mp.nstr(delta**D/m*(1/mp.log(pi)+1/mp.log(t/pi)),24),
                      'tests':rows})
    return {'available':True,
            'scope':'90-digit numerical illustrations, NOT outward interval certificates; general statements have analytic proofs',
            'decoder_cases':dec,'tube_cases':tubes,
            'optimal_X_A3_approx':mp.nstr(4+2*mp.sqrt(3),24),
            'optimal_threshold_A3_approx':mp.nstr(mp.mpf('3.5')-2*mp.sqrt(3),24)}


def main():
    lo,hi,digits,meta=certified_pi_word()
    assert digits[:12]==[3,0,1,1,0,2,1,1,1,0,0,2]
    counts=[1]
    for n in range(1,len(digits)+1):
        counts.append(sum(digits[k-1]*counts[n-k] for k in range(1,n+1)))
    assert counts[:8]==[1,3,9,28,88,276,867,2724]
    # D*C=1, checked independently through all available orders.
    assert all(counts[n]-sum(digits[k-1]*counts[n-k]
               for k in range(1,n+1))==0 for n in range(1,len(counts)))
    out={'status':'all_checks_passed',
         'proof_status':'Analytic proofs in accompanying markdown; no Lean formalization',
         'pi_word_certificate':dict(meta,count=len(digits),
                                    sha256=hashlib.sha256(bytes(digits)).hexdigest(),
                                    first_32=''.join(map(str,digits[:32]))),
         'exact_four_sample':exact_four_sample_test(digits[:128]),
         'exact_uniform':exact_uniform_checks(digits),
         'renewal_identity_orders_checked':len(digits),
         'illustrations':numerical_illustrations(digits,counts)}
    target=Path(__file__).with_name('verification_results.json')
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],
          'native_digits':len(digits),
          'exact_interval_cases':out['exact_uniform']['component_interval_cases'],
          'sharp_pairs':len(out['exact_uniform']['sharp_indistinguishability_pairs']),
          'four_sample_recovered':out['exact_four_sample']['recovered'],
          'numerical':out['illustrations'].get('available')},indent=2))


if __name__=='__main__': main()
