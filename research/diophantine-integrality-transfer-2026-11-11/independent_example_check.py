"""Check the root examples without importing certificates.py."""
from fractions import Fraction as F
from math import gcd, lcm, prod
from pathlib import Path
import hashlib
import importlib.util
import json

ROOT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    'prior_independent_pi', ROOT.parent / 'pi-joint-selection-2026-11-10' / 'verify.py')
PI = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PI)
FORMS = [(51, 2504), (51, 2503), (51, 103), (1632, 80471)]


def v(n, prime):
    n, count = abs(n), 0
    assert n
    while n % prime == 0:
        n //= prime
        count += 1
    return count


def run():
    path = ROOT / 'examples.json'
    data = json.loads(path.read_text())
    normalizers = 0
    for row in data['normalizer_examples']:
        x = F(row['x'])
        p, q, N = x.numerator, x.denominator, row['order']
        P = prod(a*p*p + b*q*q for a, b in FORMS)
        num, den = 149**(2*N+1)*P, q**8 + 149**2*P
        g = gcd(num, den)
        assert F(row['R']) == F(num, den)
        assert F(row['beta']) == F(51*p*p + 2504*q*q, q*q)
        assert int(row['raw_gcd']) == int(row['expected_gcd']) == g
        assert int(row['input_denominator']) == q
        assert int(row['output_denominator']) == den//g
        assert row['v149_R'] == v(num, 149) - v(den, 149)
        normalizers += 1

    quartics = 0
    for name in ['fixed_false_input_residuals', 'sharp_quartic_gap_examples']:
        for row in data[name]:
            x, b, s, t = map(F, row['coordinates'])
            # Expanded integer polynomial, instead of the root SOS evaluator.
            loss = x*x - 2*x*b + b*b + s*s*b*b + s*s*t*t - 2*s*t + 1
            L = lcm(*(r.denominator for r in [x, b, s, t]))
            assert int(row['L']) == L
            assert F(row['loss']) == loss
            assert F(row['normalized_loss']) == L**4*loss
            assert loss > 0 and (L**4*loss).denominator == 1
            if name == 'sharp_quartic_gap_examples':
                assert loss == F(1, L**4) and L**4*loss == 1
            else:
                assert x == b == F(1, 2) and s*t == 1
                assert loss == F(1, 4*t*t)
            quartics += 1

    lower, upper, terms = PI.pi_interval(256)
    old_lower, old_upper = map(F, data['pi_interval'])
    assert old_lower < lower < upper < old_upper
    cert = data['finite_pi_linear_certificate']
    forms = cert['forms']
    assert forms == [[7, 22], [113, 355]]
    assert forms[0][0]*forms[1][1] - forms[1][0]*forms[0][1] == -1
    assert all(A != 0 for A, B in forms)
    bounds = [max(abs(A*lower-B), abs(A*upper-B)) for A, B in forms]
    for (A, B), stored, bound in zip(forms, cert['eta_upper'], bounds):
        assert F(stored) == max(abs(A*old_lower-B), abs(A*old_upper-B))
        assert F(stored) >= bound
    largest = max(bounds)
    assert 56*largest < F(1, 2)
    # The next denominator really fails this pair's sufficient precision
    # condition, rather than merely exceeding an imprecisely rounded cap.
    largest_lower = max(min(abs(A*lower-B), abs(A*upper-B)) for A, B in forms)
    assert 57*largest_lower > F(1, 2)
    assert cert['denominator_min'] == 1 and cert['denominator_max'] == 56
    assert cert['error_lower_bound'] == '1/(226*q)'
    fraction_checks = 0
    for q in range(1, 57):
        c = F(1, 226*q)
        for p in range(3*q, 4*q+1):
            r = F(p, q)
            assert r < lower or r > upper
            actual_lower = lower-r if r < lower else r-upper
            assert actual_lower > c
            fraction_checks += 1
    # Fractions outside [3,4] have even larger distance; this covers all p.
    assert min(lower-3, 4-upper) > F(1, 226)
    result = {'checks': 'PASS', 'file': path.name,
              'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'normalizer_examples_reconstructed': normalizers,
              'expanded_quartic_examples_reconstructed': quartics,
              'independent_pi_bits': 256, 'dyadic_taylor_terms': terms,
              'finite_pi_range_verified': [1, 56],
              'pi_fraction_instances_checked': fraction_checks,
              'all_integer_numerators_covered_by_range_and_outside_bound': True,
              'method': 'Homogeneous integer forms, expanded quartic, independent outward-dyadic Machin intervals.',
              'scope': 'Finite known-input illustration only; no new pi result or verification of the source paper.'}
    (ROOT/'independent_examples.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    run()
