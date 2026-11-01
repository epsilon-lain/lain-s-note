#!/usr/bin/env python3
"""Exact rational controls for quantized position, apparent velocity and q² scale.

These examples test an observation mechanism, not pi's actual approximation
tail or the validity of any claimed irrationality-exponent theorem.
"""
from fractions import Fraction as F
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def quantize(x, h):
    """Nearest h-grid point; a midpoint is assigned to the larger grid index."""
    x, h = F(x), F(h)
    assert h > 0
    z = x/h + F(1, 2)
    return h*(z.numerator//z.denominator)


def circle_point(t):
    """A smooth rational parametrization of a unit-circle arc."""
    t = F(t)
    return ((1-t*t)/(1+t*t), 2*t/(1+t*t))


def squared_length(point):
    return sum((x*x for x in point), F(0))


def serialize(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serialize(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(v) for v in value]
    return value


def main():
    position_checks = 0
    for h in (F(1, 2), F(1, 3), F(1, 7), F(1, 10)):
        for numerator in range(-50, 51):
            x = F(numerator, 13)
            assert abs(quantize(x, h)-x) <= h/2
            position_checks += 1

    increment_checks = 0
    forward_velocity_checks = 0
    for h in (F(1, 3), F(1, 7), F(1, 16)):
        for numerator in range(-10, 11):
            t = F(numerator, 7)
            for k in range(1, 10):
                dt = F(k, 11)
                # Smooth trajectory y(t)=t², with derivative 2t and L=2.
                y0, y1 = t*t, (t+dt)*(t+dt)
                q0, q1 = quantize(y0, h), quantize(y1, h)
                assert abs((q1-q0)-(y1-y0)) <= h
                assert abs((q1-q0)/dt-(y1-y0)/dt) <= h/dt
                assert abs((q1-q0)/dt-2*t) <= h/dt+dt
                increment_checks += 1
                forward_velocity_checks += 1

    fixed_resolution_spikes = []
    for k in range(2, 22):
        h = F(1)
        dt = F(1, k)
        left, right = F(1, 2)-dt/2, F(1, 2)+dt/2
        observed = (quantize(right, h)-quantize(left, h))/dt
        assert observed == k
        assert (right-left)/dt == 1
        fixed_resolution_spikes.append(dict(h=h, dt=dt,
                                            true_velocity=F(1), observed_velocity=observed))

    improving_resolution_spikes = []
    for k in range(1, 21):
        h = F(1, 2*k+1)
        dt = h*h
        left, right = F(1, 2)-dt/2, F(1, 2)+dt/2
        q0, q1 = quantize(left, h), quantize(right, h)
        assert q0 == k*h and q1 == (k+1)*h
        observed = (q1-q0)/dt
        assert observed == 1/h
        assert max(abs(q0-left), abs(q1-right)) <= h/2
        improving_resolution_spikes.append(dict(h=h, dt=dt, center=F(1, 2),
                                                true_velocity=F(1), observed_velocity=observed,
                                                position_error_bound=h/2))

    circle_checks = 0
    circle_samples = []
    for h in (F(1, 3), F(1, 7), F(1, 16)):
        for k in range(-26, 27):
            t = F(k, 13)
            point = circle_point(t)
            observed = tuple(quantize(x, h) for x in point)
            assert squared_length(point) == 1
            assert squared_length(tuple(a-b for a, b in zip(point, observed))) <= h*h/2
            circle_checks += 1
    for t in (F(1, 2), F(2, 3)):
        point = circle_point(t)
        observed = tuple(quantize(x, F(1, 7)) for x in point)
        circle_samples.append(dict(t=t, true_point=point, observed_point=observed,
                                   true_radius_squared=F(1),
                                   observed_radius_squared=squared_length(observed)))
    assert circle_samples[0]['observed_radius_squared'] == F(52, 49)
    assert circle_samples[1]['observed_radius_squared'] == F(45, 49)

    # This cusp belongs to an exact statistic and requires no rounded numbers.
    nearest = lambda x: min(abs(x), abs(1-x))
    landmark_checks = []
    for k in range(2, 22):
        eps = F(1, 2*k)
        left, center, right = F(1, 2)-eps, F(1, 2), F(1, 2)+eps
        left_slope = (nearest(center)-nearest(left))/eps
        right_slope = (nearest(right)-nearest(center))/eps
        assert left_slope == 1 and right_slope == -1
        landmark_checks.append(dict(epsilon=eps, exact_left_slope=left_slope,
                                    exact_right_slope=right_slope))

    normalized_error_checks = 0
    false_hit_examples = []
    for q in range(1, 81):
        p = 3*q+1
        a = F(p, q)
        h = F(1, q*q)
        target = a+F(1, 10*q*q)
        observed = quantize(target, h)
        true_kappa = q*q*abs(target-a)
        observed_kappa = q*q*abs(observed-a)
        assert observed == a
        assert true_kappa == F(1, 10) and observed_kappa == 0
        assert abs(true_kappa-observed_kappa) <= q*q*h/2
        normalized_error_checks += 1
        if q in (1, 7, 20, 80):
            false_hit_examples.append(dict(p=p, q=q, h=h, target=target, candidate=a,
                                           true_kappa=true_kappa, observed_kappa=observed_kappa))
    # Keep the user's familiar candidate as a separate exact control.
    candidate = F(22, 7)
    h = F(1, 49)
    target = candidate+F(1, 490)
    assert quantize(target, h) == candidate

    output = dict(
        status='passed',
        scope=__doc__,
        quantizer='Q_h(x)=h floor(x/h+1/2), midpoint ties toward larger index',
        position_checks=position_checks,
        increment_checks=increment_checks,
        forward_velocity_checks=forward_velocity_checks,
        fixed_resolution_spikes=fixed_resolution_spikes,
        improving_resolution_spikes=improving_resolution_spikes,
        exact_unit_circle_checks=circle_checks,
        circle_inward_outward_samples=circle_samples,
        exact_landmark_cusp_checks=landmark_checks,
        normalized_error_checks=normalized_error_checks,
        false_normalized_hits=false_hit_examples,
        candidate_22_over_7_control=dict(candidate=candidate, h=h, target=target,
                                        true_kappa=F(1, 10), observed_kappa=F(0)),
        limitations=[
            'The targets in the false-hit family vary with q and are rational controls, not pi.',
            'The circle is a mathematical observation control, not a planetary dynamics model.',
            'No numerical rank, irrationality measure or infinite pi tail is inferred.',
            'Exact non-smooth landmark statistics are distinct from rounded measurement error.',
        ],
    )
    path = HERE/'results_quantization.json'
    path.write_text(json.dumps(serialize(output), ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(dict(status='passed', position_checks=position_checks,
                          increment_checks=increment_checks, circle_checks=circle_checks,
                          spike_sequences=40, normalized_error_checks=normalized_error_checks)))


if __name__ == '__main__':
    main()
