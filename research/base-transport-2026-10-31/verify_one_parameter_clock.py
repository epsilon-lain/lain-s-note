"""One integer input enumerates all positive rational candidates."""
from fractions import Fraction as F
from itertools import product
from math import isqrt
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('transport', HERE/'verify_transport.py')
transport = importlib.util.module_from_spec(spec)
spec.loader.exec_module(transport)


def decode(t):
    assert t >= 1
    k = (isqrt(8*(t-1)+1)-1)//2
    j = t-1-k*(k+1)//2
    return j+1, k-j+1


def encode(p, q):
    assert p >= 1 and q >= 1
    return (p+q-2)*(p+q-1)//2+p


def evaluate(digits, b):
    result = F(0)
    for digit in reversed(digits):
        result = (result+digit)/b
    return result


def main():
    pair_checks = 0
    for p, q in product(range(1, 101), repeat=2):
        assert decode(encode(p, q)) == (p, q)
        pair_checks += 1
    scale_checks = 0
    for q in range(1, 101):
        for p in range(3*q, 4*q+1):
            t = encode(p, q)
            b = F(p, q)
            assert F(t, q*q) == (b+1)**2/2-(b+3)/(2*q)+F(1, q*q)
            assert 6*q*q <= t < F(25, 2)*q*q
            scale_checks += 1
    _, _, digits = transport.certify_digits(128)
    examples = []
    for p, q in [(22, 7), (333, 106), (355, 113), (103993, 33102)]:
        t = encode(p, q)
        assert decode(t) == (p, q)
        L = 1
        while 3**L < t*t:
            L += 1
        assert L <= len(digits)
        b = F(p, q)
        short = evaluate(digits[:L], b)
        full_low = evaluate(digits, b)
        full_high = full_low+3*b**(-len(digits))/(b-1)
        assert short <= full_low <= full_high
        tail_bound = F(3, 2*t)
        # This interval controls both absolute-value deviations by Lipschitz.
        assert t*(full_high-short) <= tail_bound
        value = t*abs(short-1)
        examples.append({'p': p, 'q': q, 't': t, 'digit_depth': L,
                         'V_approx': float(value),
                         'V_exact': str(value),
                         'error_upper': str(t*(full_high-short)),
                         'permitted_error': str(tail_bound)})
    result = {'status': 'exact finite checks passed; P2848 unresolved',
              'pair_inverse_checks': pair_checks, 'scale_checks': scale_checks,
              'native_digits_certified': len(digits),
              'first_ten_pairs': [list(decode(t)) for t in range(1, 11)],
              'examples': examples}
    (HERE/'one_parameter_clock_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({**{k:v for k,v in result.items() if k != 'examples'},
                      'examples': [{k:v for k,v in x.items() if k not in
                                    ('V_exact', 'error_upper')} for x in examples]},
                     ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
