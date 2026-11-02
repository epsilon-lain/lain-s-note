# Exact finite search for combined odd-prime surplus

2026-11-02. For a prime p, call n a surplus index when
v_p(g_n)>v_p(n!). This is only a local gain; it does not say that g_n/n!
is large after all prime contributions are included.

The modular search covers every n<=50000 and every odd prime p<=500.
Only 43 of those primes need scanning: when a_r has no zero modulo p for
0<=r<p, the Lucas property makes a_n a unit for all n, so no surplus can
ever occur at p. For each remaining p, the recurrence was evaluated modulo
p^(v_p(50000!)+3), which preserves every tested divisibility by
p^(v_p(n!)+1). All arithmetic is exact.

There were 115 surplus events among these primes, and no index with two
surplus primes in this range. In particular p=17 and p=23, or p=17 and
p=41, never occurred together through n=50000.

The complete unreduced integers A_n,B_n and the exact gcd(A_n,4B_n) were
then computed at all 115 candidate indices. None had any other odd surplus
prime, even above 500. This second check ran in about 13 seconds. Thus no
simultaneous pair having at least one prime<=500 occurs through n=50000.
The computation does not exclude pairs with both primes above 500, nor
make any statement about later indices or boundedness of g_n sqrt(n)/n!.

Examples of the exact surplus sets through 50000:

- p=17: 72,216,1224,2524,28317,42123,47859.
- p=23: 176,352,3296,3580,4048,8586,8870,11645,12159,29355.

Reproduce with `scan_combined_surplus.py`, then
`certify_known_candidates.py`. Outputs are
`combined_surplus_checks.json` and `full_candidate_surplus_checks.json`.
No candidate arose from this finite search that certifies a combined lift.
