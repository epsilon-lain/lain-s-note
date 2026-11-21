// Build: g++ -std=c++17 -O3 verify_million_gmp.cpp -lgmpxx -lgmp -o verify_million
// Run:   ./verify_million [quotients_file]
// With no input, generate and certify the coefficients from pi. An optional
// Zuuv or B2 text file is read only after generation to check it term by term.
#include <gmpxx.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Integer = mpz_class;
using Digits = std::vector<unsigned long>;
constexpr std::size_t public_sample_size = 1000000;

// D. and G. Chudnovsky's formula for pi, in the integer recurrence given at
// https://flintlib.org/doc/algorithms/special-functions/constants.html
// pi = 426880 sqrt(10005) / S. In S each hypergeometric coefficient is
// multiplied by intercept + slope*k, and the signs alternate.
const Integer intercept = 13591409;
const Integer slope = 545140134;
const Integer cubic_factor = Integer("10939058860032000");  // 640320^3 / 24

struct SeriesBlock {
    Integer product_numerator;
    Integer product_denominator;
    Integer sum_numerator;
};

SeriesBlock split_series(unsigned long begin, unsigned long end) {
    if (end - begin == 1) {
        if (begin == 0) {
            return {1, 1, intercept};
        }
        Integer k = begin;
        Integer numerator = (6*k - 5) * (2*k - 1) * (6*k - 1);
        Integer denominator = cubic_factor * k * k * k;
        Integer term = numerator * (intercept + slope*k);
        if (begin % 2 == 1) {
            term = -term;
        }
        return {numerator, denominator, term};
    }
    auto middle = begin + (end - begin)/2;
    auto left = split_series(begin, middle);
    auto right = split_series(middle, end);
    return {
        left.product_numerator * right.product_numerator,
        left.product_denominator * right.product_denominator,
        left.sum_numerator * right.product_denominator
            + left.product_numerator * right.sum_numerator
    };
}

struct Interval {
    Integer lower, upper, denominator;
};

Interval enclose_pi() {
    // Consecutive magnitudes in S have ratio less than 72*42/cubic_factor < 1.
    // The partial sum with an even number of terms is below S. Including the
    // next positive term gives an upper bound. All operations below are exact.
    constexpr unsigned long count = 100000;
    const Integer scale = Integer(1) << 4200000;
    auto lower_sum = split_series(0, count);
    auto next = split_series(count, count+1);
    Integer upper_sum_numerator = lower_sum.sum_numerator * next.product_denominator
                                + lower_sum.product_numerator * next.sum_numerator;
    Integer upper_sum_denominator = lower_sum.product_denominator * next.product_denominator;
    Integer radicand = 10005 * scale * scale;
    Integer root;
    mpz_sqrt(root.get_mpz_t(), radicand.get_mpz_t());

    Integer lower = 426880 * root * upper_sum_denominator / upper_sum_numerator;
    Integer upper = 426880 * (root+1) * lower_sum.product_denominator
                  / lower_sum.sum_numerator + 1;
    if (upper <= lower || upper - lower > 20) {
        throw std::runtime_error("invalid or too wide an enclosure of pi");
    }
    return {lower, upper, scale};
}

struct Matrix {
    Integer a, b, c, d;
};

Matrix multiply(const Matrix& left, const Matrix& right) {
    return {
        left.a*right.a + left.b*right.c, left.a*right.b + left.b*right.d,
        left.c*right.a + left.d*right.c, left.c*right.b + left.d*right.d
    };
}

// This matrix is the product of [[a_i,1],[1,0]]. After the prefix with
// matrix M, x=(M.a*t+M.b)/(M.c*t+M.d). Its open CF cylinder is exactly the
// image of t>1. The returned pair is the inverse image of numerator/denominator.
std::pair<Integer,Integer> positive_tail(const Matrix& matrix, const Integer& numerator,
                                         const Integer& denominator) {
    Integer top = matrix.b*denominator - matrix.d*numerator;
    Integer bottom = matrix.c*numerator - matrix.a*denominator;
    if (bottom < 0) {
        top = -top;
        bottom = -bottom;
    }
    if (bottom <= 0 || top <= bottom) {
        throw std::runtime_error("point outside the proposed open CF cylinder");
    }
    return {top, bottom};
}

Matrix propose_digits(const Integer& input_numerator, const Integer& input_denominator,
                      std::size_t count, Digits& digits) {
    // Retain leading bits for a fast proposal. Approximate arithmetic never
    // establishes that a digit is correct; the final two cylinder tests do.
    auto used_bits = mpz_sizeinbase(input_numerator.get_mpz_t(), 2);
    auto keep_bits = 6*count + 256;
    auto shift = used_bits > keep_bits ? used_bits - keep_bits : 0;
    Integer numerator = input_numerator >> shift;
    Integer denominator = input_denominator >> shift;

    if (count <= 128) {
        Matrix prefix{1, 0, 0, 1};
        for (std::size_t i = 0; i < count; ++i) {
            if (denominator <= 0) {
                throw std::runtime_error("insufficient proposal precision");
            }
            Integer quotient, remainder;
            mpz_fdiv_qr(quotient.get_mpz_t(), remainder.get_mpz_t(),
                        numerator.get_mpz_t(), denominator.get_mpz_t());
            if (quotient < 1 || !mpz_fits_ulong_p(quotient.get_mpz_t())) {
                throw std::runtime_error("invalid proposed coefficient");
            }
            digits.push_back(quotient.get_ui());
            prefix = {prefix.a*quotient + prefix.b, prefix.a,
                      prefix.c*quotient + prefix.d, prefix.c};
            numerator = std::move(denominator);
            denominator = std::move(remainder);
        }
        return prefix;
    }
    auto half = count/2;
    auto left = propose_digits(numerator, denominator, half, digits);
    auto tail = positive_tail(left, numerator, denominator);
    auto right = propose_digits(tail.first, tail.second, count-half, digits);
    return multiply(left, right);
}

void compare_input(const char* path, const Digits& digits) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open the comparison file");
    }
    input >> std::ws;
    if (input.peek() == '#') {
        std::string header;
        std::getline(input, header);
    }
    for (std::size_t i = 0; i < public_sample_size; ++i) {
        unsigned long value;
        if (!(input >> value) || value != digits[i]) {
            throw std::runtime_error("coefficient mismatch at a[" + std::to_string(i) + "]");
        }
    }
    input >> std::ws;
    if (!input.eof() || input.bad()) {
        throw std::runtime_error("extra data or read error after the final coefficient");
    }
    std::cout << "1000000 submitted quotients match\n";
}

void report_records(const Digits& digits) {
    // For pi, 1/E(p_(i-1)/q_(i-1)) lies in (a_i,a_i+2) for i>=1:
    // the forward correction is in (0,1) and the backward one in [0,1).
    // In this sample every competing coefficient is separated from the
    // current record digit by at least 2. Any ambiguous case aborts.
    unsigned long current = 0;
    for (std::size_t i = 1; i < digits.size(); ++i) {
        unsigned long digit = digits[i];
        if (i == 1 || digit >= current + 2) {
            std::cout << "record n=" << i-1 << " next=" << digit << '\n';
            current = digit;
        } else if (digit + 2 > current) {
            throw std::runtime_error("two-tail comparison required at a[" + std::to_string(i) + "]");
        }
    }
    auto reserved = std::max_element(digits.begin()+900000, digits.begin()+public_sample_size);
    auto after_record = std::max_element(digits.begin()+453294, digits.begin()+public_sample_size);
    std::cout << "reserved a[" << reserved-digits.begin() << "]=" << *reserved << '\n';
    std::cout << "postrecord a[" << after_record-digits.begin() << "]=" << *after_record << '\n';
    std::cout << "extra a[1000000]=" << digits[public_sample_size] << '\n';
}

int main(int argc, char** argv) {
    if (argc > 2) {
        std::cerr << "Usage: verify_million [quotients_file]\n";
        return 1;
    }
    try {
        auto pi = enclose_pi();
        Digits digits;
        digits.reserve(public_sample_size+1);
        // The extra coefficient decides the error record at n=999999.
        auto prefix = propose_digits(pi.lower, pi.denominator, public_sample_size+1, digits);
        if (prefix.a*prefix.d - prefix.b*prefix.c != -1) {
            throw std::runtime_error("incorrect full matrix determinant");
        }
        positive_tail(prefix, pi.lower, pi.denominator);
        positive_tail(prefix, pi.upper, pi.denominator);
        if (argc == 2) {
            compare_input(argv[1], digits);
        }
        report_records(digits);
        std::cout << "1000001 digits independently certified\n";
    } catch (const std::exception& problem) {
        std::cerr << problem.what() << '\n';
        return 1;
    }
}
