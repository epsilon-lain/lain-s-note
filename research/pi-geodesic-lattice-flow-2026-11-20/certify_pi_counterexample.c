#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>

static void atan_interval(mpz_t lo, mpz_t hi, const mpz_t scale, unsigned long m)
{
    mpz_t total, power, term;
    mpz_inits(total, power, term, NULL);
    mpz_fdiv_q_ui(power, scale, m);
    unsigned long k = 0;
    for (; mpz_sgn(power); ++k) {
        mpz_fdiv_q_ui(term, power, 2*k + 1);
        if (k % 2 == 0) mpz_add(total, total, term);
        else mpz_sub(total, total, term);
        mpz_fdiv_q_ui(power, power, m*m);
    }
    mpz_sub_ui(lo, total, k/2 + 1);
    mpz_add_ui(hi, total, (k+1)/2 + 1);
    mpz_clears(total, power, term, NULL);
}

static int weighted_error(mpz_t lower, mpz_t upper, const mpz_t p,
                          const mpz_t q, const mpz_t lo, const mpz_t hi,
                          const mpz_t scale)
{
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mul(tmp, p, scale);
    mpz_mul(lower, q, lo); mpz_sub(lower, lower, tmp);
    mpz_mul(upper, q, hi); mpz_sub(upper, upper, tmp);
    if (mpz_sgn(lower) * mpz_sgn(upper) <= 0) {
        mpz_clear(tmp);
        return 0;
    }
    mpz_abs(lower, lower); mpz_abs(upper, upper);
    if (mpz_cmp(lower, upper) > 0) mpz_swap(lower, upper);
    mpz_mul(lower, lower, q); mpz_mul(upper, upper, q);
    mpz_ui_pow_ui(tmp, 10, 26);
    mpz_mul(lower, lower, tmp); mpz_fdiv_q(lower, lower, scale);
    mpz_mul(upper, upper, tmp); mpz_cdiv_q(upper, upper, scale);
    mpz_clear(tmp);
    return 1;
}

static int verify_witness(const char *filename, unsigned long places)
{
    if (places < 100 || places > 10000000) return 1;
    mpz_t scale, p, q, l18, u18, l57, u57, l239, u239, lo, hi, lower, upper, divisor;
    mpz_inits(scale, p, q, l18, u18, l57, u57, l239, u239, lo, hi, lower, upper, divisor, NULL);
    FILE *file = fopen(filename, "r");
    if (!file) return 4;
    int read = gmp_fscanf(file, "%Zd%Zd", p, q);
    fclose(file);
    if (read != 2 || mpz_sgn(q) <= 0) return 4;
    mpz_gcd(divisor, p, q);
    if (mpz_cmp_ui(divisor, 1) != 0) return 4;

    mpz_ui_pow_ui(scale, 10, places);
    atan_interval(l18, u18, scale, 18);
    atan_interval(l57, u57, scale, 57);
    atan_interval(l239, u239, scale, 239);
    /* Gauss identity: pi=48 atan(1/18)+32 atan(1/57)-20 atan(1/239). */
    mpz_mul_ui(lo, l18, 48); mpz_addmul_ui(lo, l57, 32); mpz_submul_ui(lo, u239, 20);
    mpz_mul_ui(hi, u18, 48); mpz_addmul_ui(hi, u57, 32); mpz_submul_ui(hi, l239, 20);
    if (!weighted_error(lower, upper, p, q, lo, hi, scale)) return 3;
    gmp_printf("Gauss interval for the primitive witness: [%Zd, %Zd] * 10^(-26)\n", lower, upper);
    mpz_clears(scale, p, q, l18, u18, l57, u57, l239, u239, lo, hi, lower, upper, divisor, NULL);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc == 4 && strcmp(argv[1], "verify") == 0)
        return verify_witness(argv[2], strtoul(argv[3], NULL, 10));
    unsigned long places = argc > 1 ? strtoul(argv[1], NULL, 10) : 250000;
    if (places < 100 || places > 10000000 || argc > 3) return 1;
    unsigned long *digits = calloc(2*places, sizeof *digits);
    if (!digits) return 1;

    mpz_t scale, lo5, hi5, lo239, hi239, lo, hi;
    mpz_inits(scale, lo5, hi5, lo239, hi239, lo, hi, NULL);
    mpz_ui_pow_ui(scale, 10, places);
    atan_interval(lo5, hi5, scale, 5);
    atan_interval(lo239, hi239, scale, 239);
    mpz_mul_ui(lo, lo5, 16); mpz_submul_ui(lo, hi239, 4);
    mpz_mul_ui(hi, hi5, 16); mpz_submul_ui(hi, lo239, 4);

    mpz_t a, b, c, d, q1, q2, r1, r2;
    mpz_inits(a,b,c,d,q1,q2,r1,r2,NULL);
    mpz_set(a,lo); mpz_set(b,scale);
    mpz_set(c,hi); mpz_set(d,scale);
    unsigned long count = 0, max_index = 0;
    while (mpz_sgn(b) && mpz_sgn(d) && count < 2*places) {
        mpz_fdiv_qr(q1,r1,a,b); mpz_fdiv_qr(q2,r2,c,d);
        if (mpz_cmp(q1,q2) != 0) break;
        if (!mpz_fits_ulong_p(q1)) return 2;
        digits[count] = mpz_get_ui(q1);
        if (digits[count] > digits[max_index]) {
            max_index = count;
            printf("digit record: a[%lu]=%lu\n",count,digits[count]);
        }
        ++count;
        mpz_swap(a,b); mpz_swap(b,r1);
        mpz_swap(c,d); mpz_swap(d,r2);
    }
    printf("Certified terms: %lu\n",count);

    mpz_t p,q,p_prev,q_prev,tmp,err_lo,err_hi;
    mpz_inits(p,q,p_prev,q_prev,tmp,err_lo,err_hi,NULL);
    mpz_set_ui(p,1); mpz_set_ui(q_prev,1);
    for (unsigned long i = 0; i < max_index; ++i) {
        mpz_set(tmp,p_prev); mpz_addmul_ui(tmp,p,digits[i]);
        mpz_swap(p_prev,p); mpz_swap(p,tmp);
        mpz_set(tmp,q_prev); mpz_addmul_ui(tmp,q,digits[i]);
        mpz_swap(q_prev,q); mpz_swap(q,tmp);
    }
    if (!weighted_error(err_lo, err_hi, p, q, lo, hi, scale)) return 3;
    gmp_printf("weighted error bracket: [%Zd, %Zd] * 10^(-26)\n",err_lo,err_hi);
    char *q_text = mpz_get_str(NULL, 10, q);
    printf("target p/q is convergent %lu; q has %zu decimal digits\n", max_index-1, strlen(q_text));
    free(q_text);
    if (argc == 3) {
        FILE *witness = fopen(argv[2], "w");
        if (!witness) return 4;
        int written = gmp_fprintf(witness, "%Zd\n%Zd\n",p,q);
        if (fclose(witness) || written < 0) return 4;
    }
    mpz_clears(scale, lo5, hi5, lo239, hi239, lo, hi, a, b, c, d, q1, q2, r1, r2,
               p, q, p_prev, q_prev, tmp, err_lo, err_hi, NULL);
    free(digits);
    return 0;
}
