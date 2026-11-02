"""Exact complete-period polynomials and complete small exterior sums.

Uses integer Bareiss determinants and rational polynomial interpolation.
The computation is formal in Z; Z=1 never substitutes for the actual period.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb, factorial, prod, gcd, lcm
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


@lru_cache(None)
def stirling(n, h):
    if n == 0:
        return int(h == 0)
    if h <= 0 or h > n:
        return 0
    return stirling(n-1, h-1)+h*stirling(n-1, h)


def selection(hmax, k, v0, v1):
    rows = [(j, s, b) for j in range(k) for b in range(hmax//v1+1)
            for s in range(hmax//v0+1) if v0*s+v1*b<hmax]
    counts = Counter(s+b for j, s, b in rows if j == 0)
    packed = {n for n, count in counts.items() if 2*count-1 <= n}
    targets = []
    for j, s, b in rows:
        n = s+b
        a = b+counts[n] if n in packed else b
        targets.append((j, n-a, a))
    return rows, targets, counts, packed


def entry(row, column, center, shifted=False):
    j, s, b = row
    h, a = column
    if b > a:
        return 0
    value = 0
    for ell in range(min(s, a-b)+1):
        n = s-ell
        time = factorial(h)*stirling(n, h) if shifted else h**n
        value += (comb(a, b)*comb(a-b, ell)*(factorial(s)//factorial(n))
                  * time*(j*center)**(a-b-ell))
    return value


def det_bareiss(matrix):
    a = [row[:] for row in matrix]
    n = len(a)
    sign, previous = 1, 1
    if not n:
        return 1
    for k in range(n-1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        diagonal = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = a[i][j]*diagonal-a[i][k]*a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator//previous
            a[i][k] = 0
        previous = diagonal
    return sign*a[-1][-1]


def trim(poly):
    while len(poly)>1 and not poly[-1]:
        poly.pop()
    return poly


def add(first, second):
    result = [Fraction(0)]*max(len(first), len(second))
    for i, value in enumerate(first):
        result[i] += value
    for i, value in enumerate(second):
        result[i] += value
    return trim(result)


def scale(poly, value):
    return trim([x*value for x in poly])


def eval_poly(poly, value):
    out = 0
    for coefficient in reversed(poly):
        out = out*value+coefficient
    return out


def interpolate(values):
    """Forward differences in the falling-factorial polynomial basis."""
    basis = [Fraction(1)]
    result = [Fraction(0)]
    differences = [Fraction(x) for x in values]
    for k in range(len(values)):
        result = add(result, scale(basis, differences[0]))
        differences = [b-a for a, b in zip(differences, differences[1:])]
        next_basis = [Fraction(0)]*(len(basis)+1)
        for i, value in enumerate(basis):
            next_basis[i] -= value*k/(k+1)
            next_basis[i+1] += value/(k+1)
        basis = next_basis
    return trim(result)


def period_polynomial(rows, columns, shifted=False):
    degree_bound = sum(a for h, a in columns)-sum(b for j, s, b in rows)
    if degree_bound < 0:
        return [Fraction(0)]
    denominator = prod(factorial(s) for j, s, b in rows)
    def evaluate(z):
        matrix = [[entry(row, col, z, shifted) for col in columns] for row in rows]
        return Fraction(det_bareiss(matrix), denominator)
    values = [evaluate(z) for z in range(degree_bound+1)]
    out = interpolate(values)
    # Independent evaluation points beyond the interpolation samples.
    assert eval_poly(out, -1) == evaluate(-1)
    assert eval_poly(out, degree_bound+2) == evaluate(degree_bound+2)
    return out


def sparse(poly):
    return {str(i): str(x) for i, x in enumerate(poly) if x}


def primitive_positive_scale(poly):
    """Clear denominators and positive content without changing signs."""
    poly = trim([Fraction(x) for x in poly])
    denominator = lcm(*(x.denominator for x in poly))
    integers = [int(x*denominator) for x in poly]
    content = gcd(*integers)
    return [Fraction(x//content) for x in integers] if content else [Fraction(0)]


def remainder(dividend, divisor):
    out = list(dividend)
    while len(out) >= len(divisor) and any(out):
        shift = len(out)-len(divisor)
        factor = out[-1]/divisor[-1]
        for i, value in enumerate(divisor):
            out[i+shift] -= factor*value
        trim(out)
    return out


def sturm_certificate(poly):
    initial = primitive_positive_scale(poly)
    derivative = primitive_positive_scale([i*initial[i] for i in range(1, len(initial))])
    chain = [initial, derivative]
    while True:
        next_poly = primitive_positive_scale(scale(remainder(chain[-2], chain[-1]), -1))
        if not any(next_poly):
            break
        chain.append(next_poly)
    def infinity_signs(negative):
        return [(-1 if p[-1]<0 else 1)*((-1)**(len(p)-1) if negative else 1)
                for p in chain]
    minus, plus = infinity_signs(True), infinity_signs(False)
    variations = lambda signs: sum(a != b for a, b in zip(signs, signs[1:]))
    real_count = variations(minus)-variations(plus)
    return {'polynomial_degree': len(initial)-1,
            'sturm_chain_degrees': [len(p)-1 for p in chain],
            'signs_at_minus_infinity': minus, 'signs_at_plus_infinity': plus,
            'variations_at_minus_infinity': variations(minus),
            'variations_at_plus_infinity': variations(plus),
            'real_root_count': real_count,
            'gcd_with_derivative_degree': len(chain[-1])-1,
            'nonreal_roots_with_multiplicity': len(initial)-1-real_count}


def gaussian_evaluation(poly, p, q):
    """Evaluate at 2ip/q exactly, as a pair of rational components."""
    x = Fraction(2*p, q)
    real, imaginary = Fraction(0), Fraction(0)
    for d, value in enumerate(poly):
        if d % 2:
            imaginary += value*x**d*(-1)**((d-1)//2)
        else:
            real += value*x**d*(-1)**(d//2)
    return real, imaginary


def aggregate_coefficients(poly, degree_surplus):
    result = defaultdict(lambda: [Fraction(0)])
    for d, value in enumerate(poly):
        if not value:
            continue
        assert d <= degree_surplus
        for gamma in range(degree_surplus-d+1):
            term = [Fraction(0)]*d+[value*comb(degree_surplus-d, gamma)*(-1)**gamma]
            result[gamma] = add(result[gamma], term)
    return dict(result)


def small_complete_sum():
    rows = [(j, s, 0) for j in range(2) for s in range(3)]
    columns = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 2)]
    degree_surplus = 4
    poly = period_polynomial(rows, columns)
    expected = [Fraction(0), Fraction(0), Fraction(0), Fraction(3), Fraction(1)]
    assert poly == expected
    summed = defaultdict(lambda: [Fraction(0)])
    terms = []
    for choices in product(*(range(s+1) for j, s, b in rows)):
        targets = [(j, s-a, a) for (j, s, b), a in zip(rows, choices)]
        gamma = sum(choices)
        target_poly = period_polynomial(targets, columns)
        contribution = scale(target_poly, (-1)**gamma)
        summed[gamma] = add(summed[gamma], contribution)
        if any(target_poly):
            terms.append({'target_rows': targets, 'gamma': gamma,
                          'outside_rows': sum(a>0 for j, t, a in targets),
                          'period_minor': sparse(target_poly),
                          'signed_contribution': sparse(contribution)})
    predicted = aggregate_coefficients(poly, degree_surplus)
    for gamma in set(summed)|set(predicted):
        assert summed.get(gamma, [Fraction(0)]) == predicted.get(gamma, [Fraction(0)])
    return {'source_rows': rows, 'columns': columns, 'D': degree_surplus,
            'period_polynomial': sparse(poly),
            'aggregate_coefficients': {str(g): sparse(p) for g, p in sorted(summed.items())},
            'target_selection_count': 36, 'nonzero_target_minor_count': len(terms),
            'nonzero_terms': terms}


def genuine_case():
    case = json.loads((HERE/'genuine_32_fixture.json').read_text())
    hmax = case['H']
    rows, targets, counts, packed = selection(hmax, case['K'], case['v0'], case['v1'])
    assert len(rows) == len(set(targets)) == 32
    assert len(set(targets)-set(rows)) == case['outside_rows'] == 22
    assert all(h*case['w0']+a*case['w1']<=hmax for h, a in case['columns'])
    output = []
    for shifted in (False, True):
        columns = [tuple(col) for col in case['columns']]
        poly = period_polynomial(rows, columns, shifted)
        target_poly = period_polynomial(targets, columns, shifted)
        assert any(poly) and any(target_poly)
        source_at_r = gaussian_evaluation(poly, *case['rational_approximation'])
        assert any(source_at_r)
        surplus = sum(a for h, a in columns)-sum(b for j, s, b in rows)
        valuation = next(i for i, x in enumerate(poly) if x)
        gamma = case['gamma']
        aggregate = aggregate_coefficients(poly, surplus)
        sturm = sturm_certificate(poly[valuation:])
        assert sturm['real_root_count'] == 2
        assert sturm['gcd_with_derivative_degree'] == 0
        assert sturm['polynomial_degree'] == 16
        output.append({'basis': 'shifted (Y-1)^h X^a' if shifted else 'raw Y^h X^a',
                       'H': hmax, 'rows': len(rows),
                       'outside_rows': case['outside_rows'], 'gamma': gamma,
                       'D': surplus, 'source_polynomial': sparse(poly),
                       'source_polynomial_valuation': valuation,
                       'source_polynomial_degree': len(poly)-1,
                       'rational_source_exactly_nonzero': any(source_at_r),
                       'quotient_sturm_certificate': sturm,
                       'aggregate_degree_in_zeta': surplus-valuation,
                       'selected_target_period_polynomial': sparse(target_poly),
                       'full_aggregate_coefficient_at_gamma': sparse(aggregate.get(gamma, [Fraction(0)])),
                       'complete_coefficient_at_gamma_is_zero':
                           not any(aggregate.get(gamma, [Fraction(0)]))})
        print(json.dumps({k: v for k, v in output[-1].items()
                          if k not in {'source_polynomial', 'selected_target_period_polynomial',
                                       'full_aggregate_coefficient_at_gamma'}}), flush=True)
    return output


def controlled_replacements():
    """A complete weighted packet with two separated staircase exchanges."""
    hmax, k, v0, v1, w0, w1 = 25, 2, 2, 8, 1, 1
    rows = [(j, s, b) for j in range(k) for b in range(hmax//v1+1)
            for s in range(hmax//v0+1) if v0*s+v1*b<hmax]
    rows.sort(key=lambda row: (row[1], row[0], row[2]))
    counts = {s: (hmax-v0*s+v1-1)//v1 for s in range((hmax+v0-1)//v0)}
    original = [(s, a) for s, count in counts.items() for a in range(k*count)]
    times = [6, 10]
    changed = list(original)
    multiplier = 1
    for t in times:
        count = counts[t]
        assert counts[t-1] == count
        deleted = (t, k*count-1)
        inserted = (t-1, k*count)
        assert inserted not in changed
        changed[changed.index(deleted)] = inserted
        multiplier *= k*count
    assert all(w0*h+w1*a<=hmax for h, a in changed)
    assert all(b-a>=2 for a, b in zip(times, times[1:]))
    expected_degree = sum(count*count for count in counts.values())
    original_poly = period_polynomial(rows, original, shifted=True)
    changed_poly = period_polynomial(rows, changed, shifted=True)
    assert original_poly == [Fraction(0)]*expected_degree+[Fraction(1)]
    assert changed_poly == scale(original_poly, multiplier)
    degree_surplus = sum(a for h, a in changed)-sum(b for j, s, b in rows)
    assert degree_surplus == expected_degree+len(times)
    aggregate = aggregate_coefficients(changed_poly, degree_surplus)
    assert aggregate == {g: scale(changed_poly, (-1)**g*comb(len(times), g))
                         for g in range(len(times)+1)}
    return {'H': hmax, 'K': k, 'v0': v0, 'v1': v1, 'w0': w0, 'w1': w1,
            'rows': len(rows), 'replacement_times': times,
            'replacement_count': len(times), 'multiplier': multiplier,
            'original_degree': expected_degree, 'new_D': degree_surplus,
            'aggregate_factor': '(1-zeta)^2',
            'original_polynomial': sparse(original_poly),
            'changed_polynomial': sparse(changed_poly)}


def main():
    small = small_complete_sum()
    print('Small complete sum verified:', json.dumps({
        k: small[k] for k in ['period_polynomial', 'aggregate_coefficients',
                              'target_selection_count', 'nonzero_target_minor_count']}), flush=True)
    genuine = genuine_case()
    controlled = controlled_replacements()
    result = {'status': 'passed', 'small_complete_sum': small, 'genuine_32_row_case': genuine,
              'controlled_positive_time_surplus_example': controlled,
              'scope': 'Exact formal polynomial identities; evaluation at omega=2*pi*i '
                       'is literal. No period specialization or magnitude lower bound.'}
    Path(__file__).with_name('aggregate_period_verification.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
