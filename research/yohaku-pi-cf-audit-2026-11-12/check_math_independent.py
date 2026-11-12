"""Independent, self-contained, standard-library-only runA witness checker.

Run without inputs: python check_math_independent.py
Optionally inspect old fractions too: python check_math_independent.py --artifact PATH
Optional JSON output file: --output PATH

No project code is imported or executed. No network, third-party packages, or
precomputed pi digits are needed. This is an external audit, not an action
performed by the experiment's model.

Mathematical premises and enclosure argument:
* For theta=atan(1/5), tan(2 theta)=5/12 and tan(4 theta)=120/119.
  If phi=atan(1/239), tan(4 theta-phi)=1 because
  120*239-119 == 119*239+120. Since 0 < 4 theta-phi < pi/2,
  pi = 16 atan(1/5) - 4 atan(1/239) (Machin's identity).
* atan(1/x) = sum(k>=0) (-1)^k/((2k+1)*x^(2k+1)). The alternating
  remainder has magnitude smaller than the first omitted term. Every retained
  term is enclosed by integer floor/ceiling at scale 10^600. Signs are handled
  outward, and a ceiling of the remainder is added to both endpoints.
* The CF transform t -> 1/(t-a) reverses endpoints. A term is accepted only
  when both rational endpoints have the same floor and positive remainders.
  Thus all 432 extracted CF terms apply to pi, with a_0=3. The usual integer
  convergent recurrence gives p_430/q_430; its denominator has 216 digits.
* We then bound v=q_430*(q_430*pi-p_430) directly from the independent pi
  interval. We also check 0 < q_430*pi-p_430 < 1/2, so p_430 is the nearest
  integer and v=q_430*||q_430*pi||. Strict integer cross-products prove the
  claimed forty-decimal-place enclosure. There is no floating-point math.

This establishes a finite witness and hence an infimum upper bound. It proves
nothing about pi's asymptotic liminf or whether pi is badly approximable.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

if not __debug__:
    raise RuntimeError("Run without -O/-OO; this checker requires its assertions.")
sys.set_int_max_str_digits(0)
EXPECTED_HASH = "324a5a829cd5b95bd972e80cf904a69ee1b7e3bed3ea44e5653c8a663200f4d4"


def arctan_inverse_bounds(x, terms, scale):
    """Integer enclosure of scale*atan(1/x), using outward rounding.

    Each rational series term lies in [floor, floor+1]. The alternating
    remainder has absolute size strictly below the first omitted term.
    Conservatively add that remainder to both ends of the interval.
    """
    lower = upper = 0
    for k in range(terms):
        denominator = (2 * k + 1) * x ** (2 * k + 1)
        term_floor = scale // denominator
        if k % 2 == 0:
            lower += term_floor
            upper += term_floor + 1
        else:
            lower -= term_floor + 1
            upper -= term_floor
    denominator = (2 * terms + 1) * x ** (2 * terms + 1)
    remainder_ceiling = (scale + denominator - 1) // denominator
    return lower - remainder_ceiling, upper + remainder_ceiling


def decimal_grid(numerator, denominator, digits):
    quotient, remainder = divmod(numerator * 10**digits, denominator)
    upper = quotient + (remainder != 0)
    return {
        "digits": digits,
        "floor_integer": quotient,
        "ceil_integer": upper,
        "remainder_nonzero": remainder != 0,
        "floor_decimal": f"{quotient // 10**digits}.{quotient % 10**digits:0{digits}d}",
        "ceil_decimal": f"{upper // 10**digits}.{upper % 10**digits:0{digits}d}",
        "strict_crossmultiply": quotient * denominator < numerator * 10**digits < upper * denominator,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--artifact", type=Path, help="Optional historical pi-cf-certified-summary.json; verifies its fixed SHA-256 first")
    parser.add_argument("--output", type=Path, help="Optional output JSON path (stdout is always printed)")
    args = parser.parse_args()
    digest = None
    fraction_checks = {}
    if args.artifact:
        raw = args.artifact.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        assert digest == EXPECTED_HASH
        artifact = json.loads(raw)
        for name in ["runA", "runB"]:
            bound = artifact[name]["min_markov_upper_bound"]
            numerator, denominator = int(bound["upper_num"]), int(bound["upper_den"])
            assert numerator > 0 and denominator > 0
            fraction_checks[name] = {
                "numerator_digits": len(bound["upper_num"]),
                "denominator_digits": len(bound["upper_den"]),
                "at_n": bound["at_n"],
                "decimal_18": decimal_grid(numerator, denominator, 18),
                "decimal_40": decimal_grid(numerator, denominator, 40),
                "stored_decimal_40_matches_floor": decimal_grid(numerator, denominator, 40)["floor_decimal"] == bound["decimal_40"],
            }
    assert 120*239 - 119 == 119*239 + 120
    # pi = 16*atan(1/5) - 4*atan(1/239). All endpoints are integers / 10^600.
    scale = 10**600
    l5, h5 = arctan_inverse_bounds(5, 450, scale)
    l239, h239 = arctan_inverse_bounds(239, 150, scale)
    pi_lower, pi_upper = 16*l5 - 4*h239, 16*h5 - 4*l239
    assert 3*scale < pi_lower < pi_upper < 4*scale
    # Derive the first 432 CF terms while checking both rational endpoints.
    lower_num, lower_den, upper_num, upper_den = pi_lower, scale, pi_upper, scale
    p_minus2, p_minus1, q_minus2, q_minus1 = 0, 1, 1, 0
    terms, p430, q430 = [], None, None
    for n in range(432):
        a = lower_num // lower_den
        assert a == upper_num // upper_den
        terms.append(a)
        p = a*p_minus1 + p_minus2
        q = a*q_minus1 + q_minus2
        if n == 430:
            p430, q430 = p, q
        p_minus2, p_minus1 = p_minus1, p
        q_minus2, q_minus1 = q_minus1, q
        rem_lower = lower_num - a*lower_den
        rem_upper = upper_num - a*upper_den
        assert rem_lower > 0 and rem_upper > 0
        lower_num, lower_den, upper_num, upper_den = upper_den, rem_upper, lower_den, rem_lower
    assert terms[431] == 20776 and terms.index(max(terms)) == 431
    # n=430 is even: this convergent is below pi. Check this directly too.
    error_lo = q430*pi_lower - p430*scale
    error_hi = q430*pi_upper - p430*scale
    assert 0 < error_lo < error_hi and 2*error_hi < scale
    v_lo, v_hi = q430*error_lo, q430*error_hi
    floor_40 = int("481305612462894046281170790752902942")
    # Both pi endpoints are inside the claimed one-unit decimal grid.
    lower_ok = floor_40*scale < v_lo*10**40
    upper_ok = v_hi*10**40 < (floor_40+1)*scale
    assert lower_ok and upper_ok
    report = {
        "format": "external-v05-math-integer-check/1",
        "uploaded_code_executed": False,
        "artifact_sha256": digest,
        "artifact_fractions": fraction_checks,
        "independent_runA_witness": {
            "pi_enclosure_method": "Machin identity; independent outward-rounded integer arctangent series",
            "scale_decimal_digits": 600,
            "atan_inverse5_terms": 450,
            "atan_inverse239_terms": 150,
            "pi_enclosure_width_scaled_integer": pi_upper-pi_lower,
            "certified_cf_terms": len(terms),
            "a431": terms[431],
            "q430_digits": len(str(q430)),
            "q431_digits": len(str(q)),
            "p430": str(p430),
            "q430": str(q430),
            "q431": str(q),
            "p430_is_nearest_integer_to_q430_pi": 2*error_hi < scale,
            "strict_40digit_lower": "0.0000481305612462894046281170790752902942",
            "strict_40digit_upper": "0.0000481305612462894046281170790752902943",
            "strict_18digit_upper": "0.000048130561246290",
            "lower_integer_inequality_pass": lower_ok,
            "upper_integer_inequality_pass": upper_ok,
            "scope": "A finite witness for an infimum bound, not a theorem on pi's tail liminf or badly approximable status.",
        },
    }
    if args.output:
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
