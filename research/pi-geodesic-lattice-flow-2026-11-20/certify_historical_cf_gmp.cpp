// Build: g++ -std=c++17 -O3 certify_historical_cf_gmp.cpp -lgmpxx -lgmp -o certify_historical_cf_gmp
// Run: ./certify_historical_cf_gmp 460000
#include <gmpxx.h>

#include <iostream>
#include <string>
#include <utility>

struct Interval {
    mpz_class lo;
    mpz_class hi;
};


// If a_j = (-1/base^2)^j/(2j+1), then for a block [l,h)
// P/Q=a_h/a_l and T/Q=sum_{j=l}^{h-1} a_j/a_l.
struct SeriesBlock {
    mpz_class p;
    mpz_class q;
    mpz_class t;
};

SeriesBlock split_series(unsigned long l, unsigned long h, unsigned long base) {
    if (h-l == 1) {
        mpz_class p = -mpz_class(2*l+1), q = mpz_class(base*base)*(2*l+3);
        return {p,q,q};
    }
    const unsigned long mid = l+(h-l)/2;
    const auto left = split_series(l,mid,base);
    const auto right = split_series(mid,h,base);
    return {left.p*right.p, left.q*right.q, left.t*right.q + left.p*right.t};
}

Interval arctan_reciprocal(const mpz_class& scale, unsigned long bits,
                           unsigned long base, unsigned long lower_log2_base) {
    const unsigned long terms = bits/(2*lower_log2_base)+1;
    const auto block = split_series(0,terms,base);
    const mpz_class numerator = block.t*scale;
    const mpz_class denominator = base*block.q;
    mpz_class lo, hi;
    mpz_fdiv_q(lo.get_mpz_t(),numerator.get_mpz_t(),denominator.get_mpz_t());
    mpz_cdiv_q(hi.get_mpz_t(),numerator.get_mpz_t(),denominator.get_mpz_t());
    // The next alternating term has magnitude < 1/scale and determines
    // the sign of the unknown remainder.
    if (terms%2 == 0) ++hi;
    else --lo;
    return {lo,hi};
}

void print_decimal_bounds(const mpq_class& lo, const mpq_class& hi) {
    const mpz_class decimal_scale("1000000000000000000");
    mpz_class lo_numerator = lo.get_num() * decimal_scale;
    mpz_class hi_numerator = hi.get_num() * decimal_scale;
    mpz_class lower, upper;
    mpz_fdiv_q(lower.get_mpz_t(),lo_numerator.get_mpz_t(),lo.get_den_mpz_t());
    mpz_cdiv_q(upper.get_mpz_t(),hi_numerator.get_mpz_t(),hi.get_den_mpz_t());
    std::cout << "error [" << lower << "," << upper << "]*1e-18";
}

int main(int argc, char** argv) {
    const unsigned long total = argc == 2 ? std::stoul(argv[1]) : 30000;
    const unsigned long bits = 5*total + 512;
    const mpz_class scale = mpz_class(1) << bits;
    const auto a = arctan_reciprocal(scale, bits, 5, 2);
    const auto b = arctan_reciprocal(scale, bits, 239, 7);
    mpz_class ln = 16*a.lo - 4*b.hi, ld = scale;
    mpz_class hn = 16*a.hi - 4*b.lo, hd = scale;
    mpz_class q_before = 1, q = 0;
    mpq_class best_lo = 10, best_hi = 10;
    for (unsigned long i = 0; i < total; ++i) {
        mpz_class l, h, lrem, hrem;
        mpz_fdiv_qr(l.get_mpz_t(), lrem.get_mpz_t(), ln.get_mpz_t(), ld.get_mpz_t());
        mpz_fdiv_qr(h.get_mpz_t(), hrem.get_mpz_t(), hn.get_mpz_t(), hd.get_mpz_t());
        if (l != h || lrem == 0) return 1;
        if (i >= 1 && mpq_class(l+2) >= 1/best_hi) {
            // Current complete quotient is in [ln/ld,hn/hd].
            // Invert after adding the backward complete quotient q_before/q.
            mpq_class lo(hd*q, hn*q + hd*q_before);
            mpq_class hi(ld*q, ln*q + ld*q_before);
            lo.canonicalize();
            hi.canonicalize();
            if (hi < best_lo) {
                best_lo = lo;
                best_hi = hi;
                std::cout << "index " << i-1 << ", next digit " << l << ", ";
                print_decimal_bounds(lo,hi);
                std::cout << ", denominator digits " << q.get_str().size() << "\n";
            } else if (lo <= best_hi) {
                std::cerr << "Insufficient interval precision at " << i-1 << "\n";
                return 2;
            }
        }
        q_before.swap(q);
        q += l * q_before;
        ln = std::move(hd);
        hn = std::move(ld);
        ld = std::move(hrem);
        hd = std::move(lrem);
    }
    std::cout << "Certified " << total << " partial quotients\n";
}
