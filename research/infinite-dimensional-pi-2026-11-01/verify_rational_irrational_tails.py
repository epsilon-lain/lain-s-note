"""Exact denominator checks for ordinary Basel partial sums."""
from fractions import Fraction
from pathlib import Path
import json

def v2(n):
    count = 0
    while n % 2 == 0:
        count += 1
        n //= 2
    return count

def main():
    total = Fraction(0)
    examples = []
    for N in range(1, 257):
        total += Fraction(1, N*N)
        B = total.denominator
        k = N.bit_length()-1
        assert v2(B) == 2*k
        assert B >= 4**k and 4*B > N*N
        # Integral tail bounds, kept rational; no floating-point pi is used.
        tail_lower, tail_upper = Fraction(1,N+1), Fraction(1,N)
        assert Fraction(B,N+1) > Fraction(N*N,4*(N+1))
        assert Fraction(B*B,N+1) > Fraction(N**4,16*(N+1))
        if N in [1,2,3,4,8,16,32,64,128,256]:
            examples.append({'N':N, 'numerator':str(total.numerator),
                'denominator':str(B), 'v2_denominator':v2(B),
                'tail_lower':str(tail_lower), 'tail_upper':str(tail_upper),
                'denominator_scaled_tail_lower':str(Fraction(B,N+1))})
    # Contrast: infinitely many rational squares can sum to a rational number.
    for N in range(1,50):
        partial = sum((Fraction(1,4**n) for n in range(1,N+1)), Fraction(0))
        assert Fraction(1,3)-partial == Fraction(1,3*4**N)
    for N in range(1,21):
        # Polynomial identity (1+t²) Σ_{k=0}^{2N-1} (-1)^k t^(2k) = 1-t^(4N).
        coefficients=[0]*(4*N+1)
        for k in range(2*N):
            coefficients[2*k]+=(-1)**k
            coefficients[2*k+2]+=(-1)**k
        assert coefficients == [1]+[0]*(4*N-1)+[-1]
    result={'status':'passed','basel_prefixes':256,'geometric_prefixes':49,
            'pi_tail_polynomial_identities':20,
            'method':'Exact Fraction arithmetic; analytic tail inequalities proved in the note.',
            'examples':examples}
    out=Path(__file__).with_name('rational_irrational_tails_results.json')
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('Passed: 256 Basel denominator checks, 49 geometric sums, and 20 pi-tail polynomial identities.')

if __name__=='__main__':
    main()
