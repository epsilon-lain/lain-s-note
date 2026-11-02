"""Exact checks of the joint time/transverse raw-column cost law.

The general limits are proved in SPECIAL_TIME_COST.md. This verifier checks
polynomial density integrals and rational enclosures; it does not certify a
nonzero determinant or a feasible approximation tuple.
"""

from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import json


def cdf(m, a, c):
    n = m + 1
    if c <= 0:
        return F(0)
    if c >= max(F(1), a):
        return F(1)
    if a == 1:
        return c ** n
    if a > 1:
        if c <= 1:
            return c ** n / a
        return 1 - (a - c) ** n / (a * (a - 1) ** (n - 1))
    return (c ** n - max(F(0), c - a) ** n / (1 - a) ** (n - 1)) / a


def stop_loss(m, a, c):
    """Integral_0^c F_a(u) du = E[(c-C_a)_+]."""
    n = m + 1
    mean = F(m) / (m + 2) + a / (m + 2)
    if c <= 0:
        return F(0)
    if c >= max(F(1), a):
        return c - mean
    if a == 1:
        return c ** (n + 1) / (n + 1)
    if a > 1:
        if c <= 1:
            return c ** (n + 1) / (a * (n + 1))
        return c - mean + (a - c) ** (n + 1) / (
            a * (a - 1) ** (n - 1) * (n + 1))
    return (c ** (n + 1) - max(F(0), c - a) ** (n + 1)
            / (1 - a) ** (n - 1)) / (a * (n + 1))


def cost_mean(m, a, cutoff):
    mass = cdf(m, a, cutoff)
    assert mass > 0
    return cutoff - stop_loss(m, a, cutoff) / mass


def integrate_power_polynomial(m, coefficients, lower, upper):
    """Integrate x^(m-1) times the supplied polynomial."""
    return sum(coefficient * (upper ** (m + i) - lower ** (m + i)) / (m + i)
               for i, coefficient in enumerate(coefficients))


def joint_integrals(m, a, cutoff):
    """Independent piecewise integration of the joint triangle density."""
    right = min(F(1), cutoff)
    breaks = {F(0), right}
    if a != 1:
        crossing = (a - cutoff) / (a - 1)
        if 0 < crossing < right:
            breaks.add(crossing)
    breaks = sorted(breaks)
    mass = time = spatial = F(0)
    prefactor = m * (m + 1)
    for left, right in zip(breaks, breaks[1:]):
        midpoint = (left + right) / 2
        if 1 - midpoint <= (cutoff - midpoint) / a:
            intercept, slope = F(1), F(-1)
        else:
            intercept, slope = cutoff / a, -1 / a
        cap = [intercept, slope]
        cap_squared = [intercept ** 2, 2 * intercept * slope, slope ** 2]
        mass += prefactor * integrate_power_polynomial(m, cap, left, right)
        time += F(prefactor, 2) * integrate_power_polynomial(m, cap_squared, left, right)
        spatial += prefactor * integrate_power_polynomial(m, [F(0)] + cap, left, right)
    return mass, time, spatial


def exact_integral_checks():
    checked = lower_coordinates = upper_coordinates = 0
    for m, a in product(range(1, 9), (F(1, 4), F(1, 2), F(1), F(3, 2), F(2), F(5), F(20))):
        maximum = max(F(1), a)
        cutoffs = {maximum * x for x in (F(1, 4), F(1, 2), F(3, 4), F(1))}
        cutoffs.add(min(F(1), a))
        for cutoff in sorted(cutoffs):
            mass, time, spatial = joint_integrals(m, a, cutoff)
            assert mass == cdf(m, a, cutoff)
            assert spatial + a * time == cutoff * mass - stop_loss(m, a, cutoff)
            assert cost_mean(m, a, cutoff) == (spatial + a * time) / mass
            slack = mass - time - spatial
            assert time >= 0 and spatial >= 0 and slack >= 0
            if cutoff <= min(F(1), a):
                assert time / mass == cutoff / (a * (m + 2))
                assert spatial / mass == m * cutoff / (m + 2)
                assert cost_mean(m, a, cutoff) == F(m + 1, m + 2) * cutoff
                lower_coordinates += 1
            if a > 1 and cutoff >= 1:
                tail = 1 - mass
                delta = a - cutoff
                k = delta / (a - 1)
                assert spatial == F(m, m + 2) * (1 - tail * k)
                assert slack == (1 - tail * delta / a) / (m + 2)
                expected_cost = (F(m, m + 2) + a / (m + 2)
                                 - tail * (a - F(m + 1, m + 2) * delta)) / mass
                assert cost_mean(m, a, cutoff) == expected_cost
                upper_coordinates += 1
            checked += 1
    return {"joint_density_integrals": checked,
            "lower_branch_coordinate_checks": lower_coordinates,
            "upper_branch_coordinate_checks": upper_coordinates}


def quantile_bounds(m, a, mass=F(1, 2), steps=96):
    lower, upper = F(0), max(F(1), a)
    for _ in range(steps):
        midpoint = (lower + upper) / 2
        if cdf(m, a, midpoint) < mass:
            lower = midpoint
        else:
            upper = midpoint
    assert cdf(m, a, lower) <= mass <= cdf(m, a, upper)
    return lower, upper


def time_mean_at_cutoff(m, cutoff):
    n = m + 1
    tail = (1 - cutoff) ** n
    mass = 1 - tail
    assert mass > 0
    return (F(1, n + 1) - tail * (1 - F(n, n + 1) * (1 - cutoff))) / mass


def time_mean_bounds(m, mass=F(1, 2), steps=96):
    n = m + 1
    lower, upper = F(0), F(1)
    for _ in range(steps):
        midpoint = (lower + upper) / 2
        if 1 - (1 - midpoint) ** n < mass:
            lower = midpoint
        else:
            upper = midpoint
    return time_mean_at_cutoff(m, lower), time_mean_at_cutoff(m, upper)


def uniform_bound_checks():
    mass = F(1, 2)
    checked = 0
    for m in (1, 2, 4, 8, 16, 32):
        time_lower, time_upper = time_mean_bounds(m)
        for a in sorted({F(1, 3), F(1), F(3, 2), F(2), F(4), F(m), F(m * m), F(2 ** min(m, 20))}):
            q_lower, q_upper = quantile_bounds(m, a)
            value_lower = cost_mean(m, a, q_lower)
            value_upper = cost_mean(m, a, q_upper)
            # Use the harder interval endpoints to get actual rational certificates.
            assert value_lower >= 1 + a * time_upper - F(2) / (mass * (m + 2))
            coefficient = a - F(m, m + 1)
            time_endpoint = time_lower if coefficient >= 0 else time_upper
            upper_bound = F(m, m + 1) + coefficient * time_endpoint
            assert value_upper <= upper_bound
            if a >= 1:
                assert value_lower >= 1 + (a - 1) * time_upper - F(1) / (mass * (m + 2))
            checked += 1
    return checked


def compact_enclosure(lower, upper, digits=12):
    scale = 10 ** digits
    lo = (lower * scale).numerator // (lower * scale).denominator
    hi = -((-upper * scale).numerator // (-upper * scale).denominator)
    return [str(F(lo, scale)), str(F(hi, scale))]


def log_two_bounds(terms=32):
    lower = sum(F(2, (2 * k + 1) * 3 ** (2 * k + 1)) for k in range(terms))
    error = F(9, 4 * (2 * terms + 1) * 3 ** (2 * terms + 1))
    return lower, lower + error


def scaling_examples():
    records = []
    log_lower, log_upper = log_two_bounds()
    assert F(2, 3) < log_lower < log_upper < F(7, 10)
    d_lower, d_upper = 1 - log_upper, 1 - log_lower
    for m in (8, 16, 32, 64):
        for label, a in (("linear", F(m)), ("quadratic", F(m * m)), ("exponential", F(2 ** m))):
            q_lower, q_upper = quantile_bounds(m, a, steps=128)
            lower = cost_mean(m, a, q_lower)
            upper = cost_mean(m, a, q_upper)
            if label == "linear":
                normalized = (lower, upper)
                limit = (1 + d_lower, 1 + d_upper)
                meaning = "Q_m; proved limit 2-ln(2)"
            else:
                normalized = (m * lower / a, m * upper / a)
                limit = (d_lower, d_upper)
                meaning = "m Q_m/a; proved limit 1-ln(2)"
            records.append({"m": m, "a": str(a), "scale": label,
                            "normalized_cost_enclosure": compact_enclosure(*normalized),
                            "proved_limit_enclosure": compact_enclosure(*limit),
                            "meaning": meaning})
    return records


def raw_columns(w0, weights, height):
    return [(h, alpha)
            for alpha in product(*(range(height // w + 1) for w in weights))
            for h in range(height // w0 + 1)
            if w0 * h + sum(w * x for w, x in zip(weights, alpha)) <= height]


def finite_selection_checks():
    checked = 0
    for w0, weights, height in ((1, (2,), 3), (1, (1, 1), 2)):
        columns = raw_columns(w0, weights, height)
        number = len(columns) // 2
        for a in (F(1, 2), F(1), F(2), F(5)):
            costs = [a * F(w0 * h, height)
                     + F(sum(w * x for w, x in zip(weights, alpha)), height)
                     for h, alpha in columns]
            best = sum(sorted(costs)[:number])
            observed = min(sum(costs[i] for i in selected)
                           for selected in combinations(range(len(columns)), number))
            assert best == observed
            checked += 1
    return checked


if __name__ == "__main__":
    result = {
        "status": "passed", "integral_checks": exact_integral_checks(),
        "uniform_coefficient_bound_checks": uniform_bound_checks(),
        "finite_column_selection_checks": finite_selection_checks(),
        "scaling_examples": scaling_examples(),
        "scope": "Exact joint positive-cost law and lower-quantile column optimization. "
                 "No source-degree feasibility, time matching, determinant rank, "
                 "determinant lower bound, or pi consequence is certified.",
    }
    Path(__file__).with_name("special_time_cost_verification.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: len(v) if isinstance(v, list) else v
                      for k, v in result.items()}, indent=2))
