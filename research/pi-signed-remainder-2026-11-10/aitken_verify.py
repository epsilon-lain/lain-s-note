"""Check Aitken diagnostics from independent binomial/telescoping sources."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import importlib.util
import json

ROOT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    'joint_independent', ROOT.parent / 'pi-joint-selection-2026-11-10' / 'verify.py')
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)


def error(value, lower, upper):
    if value < lower:
        return lower - value, upper - value
    assert value > upper
    return value - upper, value - lower


def run():
    path = ROOT / 'aitken_results.json'
    data = json.loads(path.read_text())
    _, rs = CHECK.sources(max(data['parameters']['n']) + 3)
    pl, pu, terms = CHECK.pi_interval(data['parameters']['scoring_bits'])
    lower, upper = pl * pl, pu * pu
    for row in data['cases']:
        n = row['n']
        r0, r1, r2, r3 = rs[n:n + 4]
        # Determinant formula, independently of local-ratio construction.
        transformed = (r0 * r2 - r1 * r1) / (r0 - 2 * r1 + r2)
        next_transformed = (r1 * r3 - r2 * r2) / (r1 - 2 * r2 + r3)
        assert F(row['source']) == r0
        assert F(row['transformed_source']) == transformed
        d = r1 - r0
        rho = -(r2 - r1) / d
        next_rho = -(r3 - r2) / (r2 - r1)
        variation = d * rho * (1 / (1 + rho) - 1 / (1 + next_rho))
        assert variation == next_transformed - transformed
        el, eh = error(r0, lower, upper)
        tl, th = error(transformed, lower, upper)
        gl, gh = tl / eh, th / el
        factor = F(transformed.denominator ** 2, r0.denominator ** 2)
        cl, ch = factor * gl, factor * gh
        ml, mh = error(r2, lower, upper)
        matched_factor = F(transformed.denominator ** 2, r2.denominator ** 2)
        matched_l, matched_h = matched_factor * tl / mh, matched_factor * th / ml
        assert matched_h < 1 or matched_l > 1
        assert gh < 1 < cl
        assert row['analytic_error_strictly_improved']
        assert row['arithmetic_cost_strictly_worse']
        for key, lo, hi in [('error_ratio_decimal_enclosure', gl, gh),
                            ('B_squared_error_ratio_decimal_enclosure', cl, ch)]:
            dl, du = map(F, row[key])
            assert dl <= lo < hi <= du
        assert row['source_denominator_digits'] == len(str(r0.denominator))
        assert row['transformed_denominator_digits'] == len(str(transformed.denominator))
        assert F(row['same_highest_input_source']) == r2
        assert row['same_highest_input_denominator_digits'] == len(str(r2.denominator))
        dl, du = map(F, row['same_highest_input_cost_ratio_decimal'])
        assert dl <= matched_l < matched_h <= du
        assert row['same_highest_input_cost_improved'] == (matched_h < 1)
        exact = {'n': n, 'R': str(r0), 'T': str(transformed), 'rho': str(rho),
                 'rho_next': str(next_rho), 'T_next_minus_T': str(variation),
                 'source_error': list(map(str, (el, eh))),
                 'transformed_error': list(map(str, (tl, th))),
                 'error_ratio': list(map(str, (gl, gh))),
                 'arithmetic_cost_ratio': list(map(str, (cl, ch))),
                 'same_highest_input_source': str(r2),
                 'same_highest_input_cost_ratio': list(map(str, (matched_l, matched_h)))}
        digest = hashlib.sha256(json.dumps(exact, sort_keys=True,
                                          separators=(',', ':')).encode()).hexdigest()
        assert digest == row['exact_record_sha256']
    assert [r['n'] for r in data['cases']] == data['parameters']['n']
    expected = {'candidates': 9, 'analytic_improvements': 9,
                'arithmetic_cost_improvements': 0,
                'same_highest_input_cost_improvements': 1,
                'last_source_denominator_digits': 238,
                'last_transformed_denominator_digits': 371}
    assert data['summary'] == expected
    report = {'checks': 'PASS', 'file': path.name,
              'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'source_formula': 'binomial sum and telescoping companion',
              'candidate_formula': 'rational determinant rather than local-tail weights',
              'cases_verified': len(data['cases']), 'taylor_terms': terms,
              'summary': expected,
              'scope': 'Nine finite diagnostics; no general failure theorem.'}
    (ROOT / 'aitken_verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    run()
