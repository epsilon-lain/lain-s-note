"""Reproducible checks of the joint-bound note; no claim about pi approximants."""
from fractions import Fraction as F
from heapq import heapify, heappop, heappush
from math import exp, factorial, log, prod
from pathlib import Path
import json
import random


def lattice(weights, limit):
    """Exact lattice points in a weighted simplex."""
    if not weights:
        yield (), F(0)
        return
    w, *rest = weights
    for j in range(int(limit // w) + 1):
        for a, y in lattice(rest, limit - w*j):
            yield (j,) + a, w*j + y


def exact_minimum(weights, mass, cap):
    # Independent primal calculation: allocate each next unit where its
    # incremental quadratic cost is smallest. Separable convexity proves
    # optimality over this finite set.
    c = F(1, 6)
    points = list(lattice(weights, F(cap)))
    count = [0] * len(points)
    heap = [(2*y + c, i) for i, (_, y) in enumerate(points)]
    heapify(heap)
    energy = F(0)
    for _ in range(mass):
        increment, i = heappop(heap)
        energy += increment
        count[i] += 1
        y = points[i][1]
        heappush(heap, (2*y + c*(2*count[i]+1), i))
    # Any excluded point has y>cap, hence its first cost is >2*cap+c.
    # It could not improve the constructed allocation.
    assert increment < 2*cap+c
    H = F(10)
    t = (increment + heap[0][0]) / (4*H)
    assert t*H < cap
    square_sum = sum((t*H-y)**2 for _, y in lattice(weights, t*H))
    dual = 2*t*H*mass - square_sum/c
    assert energy >= dual
    direct = sum(c*n*n+2*y*n for n, (_, y) in zip(count, points))
    assert direct == energy
    return {"weights": list(map(str, weights)), "mass": mass,
            "exact_minimum": str(energy), "exact_dual": str(dual),
            "nonnegative_gap": str(energy-dual)}


def riemann_checks():
    out = []
    for weights in ((F(1),), (F(1), F(2)), (F(1), F(2), F(3))):
        m = len(weights)
        t = F(3, 4)
        ratios = []
        for H in (8, 16, 32, 64):
            exact = sum((t*H-y)**2 for _, y in lattice(weights, t*H))
            leading = 2*(t*H)**(m+2)/(factorial(m+2)*prod(weights))
            ratios.append({"H": H, "ratio": float(exact/leading)})
        assert abs(ratios[-1]["ratio"]-1) < abs(ratios[0]["ratio"]-1)
        out.append({"weights": list(map(str, weights)), "ratios": ratios})
    return out


def bound_checks():
    rng = random.Random(2848)
    c = log(2)/4
    worst = -float("inf")
    for _ in range(600):
        m = rng.choice((1, 2, 3, 10, 30, 100))
        K = rng.randint(1, 12)
        theta = rng.uniform(.2, .99)
        rho = rng.uniform(.1, .95)
        w0 = exp(rng.uniform(-2, 24))
        v0 = K*theta**m*w0/rho
        ss = [rng.uniform(0, 2) for _ in range(m)]
        tau = max(ss)
        logP = sum(log(1+s/2) for s in ss)
        T = exp((log(c*rho/w0)+logP)/(m+1))
        b = m/(m+2)
        # U>=b and tau/F0>=0. Replacing U by b and dropping tau/F0
        # makes the LHS larger; thus this is the harder numerical check.
        lhs = 2*(m+1)*T/(m+2)-theta*sum(1+s for s in ss)/(m+2)-b
        lhs -= 100*K/w0 + log(2)/v0
        delta = exp(log(c*rho/(400*(m+2)**2))/m)
        rhs = b*(delta*(1+tau*tau/(4*(1+tau)))-1)
        slack = lhs-rhs
        assert slack <= 1e-10*max(1, abs(lhs), abs(rhs))
        worst = max(worst, slack)
    thresholds = []
    for m in (1, 10, 100, 1000, 10000):
        rho = .5
        d = exp(log(400*(m+2)**2/(c*rho))/m)-1
        # tau^2/(4(1+tau))=d; equality is already insufficient.
        tau_min = 2*d+2*(d*d+d)**.5
        thresholds.append({"m": m, "rho": rho,
                           "necessary_tau_strictly_above": tau_min})
    return {"sample_count": 600, "largest_lhs_minus_rhs": worst,
            "thresholds": thresholds}


def monotone_constraint():
    m, theta, T = 3, F(4, 5), F(2, 5)
    R = (theta+T)/2
    x = R/theta
    source_cdf = (m+1)*x**m-m*x**(m+1)
    target_cdf = F(1)
    assert target_cdf > source_cdf
    return {"m": m, "theta": str(theta), "T": str(T), "R": str(R),
            "source_fraction": str(source_cdf),
            "free_target_fraction": str(target_cdf),
            "violated_mass_fraction": str(target_cdf-source_cdf)}


if __name__ == "__main__":
    result = {
        "status": "checks passed; not a proof or disproof of P2848",
        "finite_check_c": "1/6 (the finite inequality holds for every c>0)",
        "exact_primal_dual": [
            exact_minimum((F(1),), 100, 20),
            exact_minimum((F(1), F(2)), 250, 20),
            exact_minimum((F(1), F(2), F(3)), 1000, 20)],
        "riemann_sums": riemann_checks(),
        "paper_constant_bound": bound_checks(),
        "omitted_monotone_constraint": monotone_constraint(),
    }
    destination = Path(__file__).with_name("results.json")
    destination.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n")
    print(json.dumps(result, indent=2, ensure_ascii=False))
