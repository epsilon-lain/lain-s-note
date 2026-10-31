"""Independent exact-rational baseline and search regressions (stdlib only)."""
import json
import random
import unittest
from fractions import Fraction as F
from pathlib import Path
import search


def independent_pi():
    # More terms than production 160-bit arithmetic; independent direct sums.
    def bounds(a):
        s = sum((F((-1)**k, (2*k+1)*a**(2*k+1)) for k in range(100)), F(0))
        return s, s+F(1, 201*a**201)
    a, b = bounds(5)
    c, d = bounds(239)
    return 16*a-4*d, 16*b-4*c


def independent_error(p, q, lo, hi):
    x = F(p, q)
    assert x < lo or x > hi
    return (lo-x, hi-x) if x < lo else (x-hi, x-lo)


class SearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lo, cls.hi = independent_pi()

    def test_pi_certificate(self):
        for bits in (8, 16, 32, 80, 160, 256):
            lo, hi, scale = search.pi_bounds(bits)
            self.assertLess(F(lo, scale), self.lo)
            self.assertGreater(F(hi, scale), self.hi)

    def test_random_rectangles_against_all_pairs(self):
        rng = random.Random(2848)
        cases = [(1, 30, 1, 120), (50, 120, 150, 380),
                 (8, 13, 1, 10), (1, 5, 50, 57), (226, 226, 710, 710)]
        for _ in range(75):
            qmin, pmin = rng.randint(1, 50), rng.randint(1, 180)
            cases.append((qmin, qmin+rng.randint(0, 12), pmin, pmin+rng.randint(0, 18)))
        for qmin, qmax, pmin, pmax in cases:
            fast = search.scan(qmax, qmin, pmin, pmax)
            slow = search.scan(qmax, qmin, pmin, pmax, exhaustive=True)
            self.assertEqual(fast["pairs_evaluated"], qmax-qmin+1)
            self.assertEqual(slow["pairs_evaluated"], (qmax-qmin+1)*(pmax-pmin+1))
            for metric in ("absolute", "normalized"):
                self.assertEqual(fast[metric], slow[metric])
                p, q = fast[metric][-1][:2]
                low, high = independent_error(p, q, self.lo, self.hi)
                weight = q*q if metric == "normalized" else 1
                for den in range(qmin, qmax+1):
                    for num in range(pmin, pmax+1):
                        if (metric == "absolute" and num*q == p*den) or (num, den) == (p, q):
                            continue
                        a, _ = independent_error(num, den, self.lo, self.hi)
                        other_weight = den*den if metric == "normalized" else 1
                        self.assertLess(high*weight, a*other_weight)

    def test_precision_retry(self):
        coarse = search.scan(500, bits=8)
        fine = search.scan(500)
        self.assertGreater(coarse["precision_retries"], 0)
        for metric in ("absolute", "normalized"):
            self.assertEqual([x[:2] for x in coarse[metric]], [x[:2] for x in fine[metric]])

    def test_known_records_and_ties(self):
        result = search.scan(1000)
        self.assertEqual([x[:2] for x in result["normalized"]], [(3, 1), (22, 7), (355, 113)])
        self.assertEqual(result["absolute"][-1][:2], (355, 113))
        # A clipped rectangle can contain only an unreduced representation.
        only = search.scan(226, 226, 710, 710)
        data = search.report(only, 0)["normalized_records"][-1]
        self.assertEqual((data["reduced_p"], data["reduced_q"]), (355, 113))
        _, _, a, b = only["normalized"][-1]
        _, _, c, d = result["normalized"][-1]
        self.assertEqual((a, b), (2*c, 2*d))
        self.assertEqual((a*226, b*226), (4*c*113, 4*d*113))

    def test_native_digits_and_migration(self):
        digits = search.native_digits(64)
        self.assertEqual(digits[:6], [3, 0, 1, 1, 0, 2])
        # Independent interval propagation with exact Fraction, without rounding.
        lo = hi = F(1)
        for digit in digits:
            lo, hi = lo*self.lo, hi*self.hi
            self.assertEqual(lo.__floor__(), digit)
            self.assertEqual(hi.__floor__(), digit)
            lo, hi = lo-digit, hi-digit
        for p, q, sign in [(22, 7, -1), (355, 113, -1), (3, 1, 1)]:
            item = search.migration_interval(p, q, digits)
            if sign < 0:
                self.assertLess(int(item["upper_numerator"]), 0)
            else:
                self.assertGreater(int(item["lower_numerator"]), 0)
        self.assertIn("not_evaluated", search.migration_interval(1, 2, digits))

    def test_invalid_ranges(self):
        for kwargs in ({"qmax": 0}, {"qmax": 5, "qmin": 6},
                       {"qmax": 1, "pmin": 3, "pmax": 2},
                       {"qmax": 3, "exhaustive": True}):
            with self.assertRaises(ValueError):
                search.scan(**kwargs)


if __name__ == "__main__":
    unittest.main()
