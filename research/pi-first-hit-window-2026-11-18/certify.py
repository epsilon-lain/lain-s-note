"""Generate a finite first-record certificate using exact continued fractions."""
import json
from fractions import Fraction as F

DIGITS = 350
SCALE = 10 ** DIGITS


def atan_bounds(d, terms):
    partial = sum((F((-1)**k, (2*k+1)*d**(2*k+1)) for k in range(terms)), F(0))
    following = partial + F((-1)**terms, (2*terms+1)*d**(2*terms+1))
    return min(partial, following), max(partial, following)


def pi_bounds():
    a, b = atan_bounds(5, 260)
    c, d = atan_bounds(239, 80)
    lower, upper = 16*a-4*d, 16*b-4*c
    return F(lower*SCALE // 1, SCALE), F(-(-upper*SCALE // 1), SCALE)


def cost(p, q, lower, upper):
    left, right = q*(q*lower-p), q*(q*upper-p)
    return max(left, -right, 0), max(abs(left), abs(right))


def record(n, p, q, interval):
    lower, upper = interval
    grid = 10**15
    return dict(n=n, p=p, q=q, cost_floor_1e15=int(lower*grid//1),
                cost_ceil_1e15=int(-(-upper*grid//1)))


def generate():
    lower, upper = pi_bounds()
    left, right = lower, upper
    old_p, p, old_q, q = 0, 1, 1, 0
    best = None
    records = []
    for n in range(1000):
        digit = left // 1
        if digit != right // 1:
            raise ArithmeticError("Increase the precision of the pi enclosure")
        old_p, p = p, digit*p+old_p
        old_q, q = q, digit*q+old_q
        score = cost(p, q, lower, upper)
        if best is None or score[1] < best[0]:
            records.append(record(n, p, q, score))
            best = score
            if q > 113:
                return dict(pi_grid_digits=DIGITS, pi_lower_integer=int(lower*SCALE),
                            pi_upper_integer=int(upper*SCALE), records=records)
        elif score[0] <= best[1]:
            raise ArithmeticError("Ambiguous score comparison")
        if left == digit:
            raise ArithmeticError("Insufficient enclosure for reciprocal")
        left, right = 1/(right-digit), 1/(left-digit)
    raise ArithmeticError("No new record within the search budget")


if __name__ == "__main__":
    print(json.dumps(generate(), indent=2))
