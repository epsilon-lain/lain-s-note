# Exact experiments on the upper-half first-layer anomaly set

2026-11-02. The fixed recurrence is

\[
A_0=A_1=1,\qquad A_{r+1}=(2r+1)A_r+r^2A_{r-1}.
\]

For an integer n>=1, let

\[
\mathcal A(n)=\{p:\ p\text{ odd prime},\ n/2<p<n,
\ p\mid A_{n-p}\},
\quad P(n)=\prod_{p\in\mathcal A(n)}p,
\quad S(n)=\ln P(n).
\]

The empty product is one. This experiment records first-layer modulo-p anomalies. It does not equate them with positive excess beyond n! in the actual gcd.

## 1. Reflection gives a complete sparse event enumeration

Write ell_r=A_r/r!, D(z)=1-2z-z² and h=(p-1)/2. The established generating-function congruence gives

\[
D(z)^{-1/2}\equiv D(z)^hD(z^p)^{-1/2}\pmod p.
\]

For 0<=r<p, ell_r is therefore the coefficient q_r of D(z)^h modulo p. The reciprocal identity

\[
z^{2h}D(-1/z)^h=(-1)^hD(z)^h
\]

implies

\[
\ell_{p-1-r}\equiv(-1)^{h+r}\ell_r\pmod p.
\]

Since r! and (p-1-r)! are p-units,

\[
\boxed{p\mid A_r\iff p\mid A_{p-1-r}.}
\]

Thus a root at r produces anomaly events at the paired indices

\[
n=p+r,\qquad n'=2p-1-r=3p-1-n.
\]

The scan only evaluates A_r modulo p for
r<=min((p-1)/2,N-p). If a visible event has r above the midpoint, its reflected index is smaller than r and is inside this computed interval. This proves the enumeration covers every event with n<=N. It does not skip primes or require factorization of a large A_n.

The code also independently computes every relevant r without reflection through n=10000 and compares the complete event lists exactly. Agreement is a finite implementation audit; the congruence proves the general reflection step.

## 2. Finite certification through n=100000

The exact scan certifies

\[
\boxed{S(n)<n/6\quad(1\le n\le100000).}
\]

Its unique maximum S(n)/n remains

\[
\boxed{S(21)/21=\ln17/21=0.1349149211\ldots.}
\]

Neither statement uses floating-point logarithms. For each index other than 21, the integer certificates establish both

\[
P(n)^6<2^n<e^n,
\qquad P(n)^{21}<16^n<17^n.
\]

The first proves S(n)<n/6; the second proves its ratio is strictly below ln17/21. The comparisons against powers of two use integer bit lengths. At n=21, P(21)=17 and

\[
3^{21}17^6<8^{21},\qquad e>8/3,
\]

give the remaining strict n/6 comparison. The lower bound for e follows from its positive series. Floating logs are only labeled diagnostics for ranking the stored samples.

The scan records 7026 anomaly events at 6800 indices. The cluster-size histogram is:

| Number of first-layer anomalous primes at one index | Indices through 100000 |
|---|---:|
| 0 | 93200 |
| 1 | 6580 |
| 2 | 214 |
| 3 | 6 |

The ascending exact event stream has SHA-256
`4f51a87309af96483f97affe51bbb5734a2ab8bb4039e7a4d8ac8aca60c94f3f`,
using ASCII `n:p` followed by newline, ordered first by n, then by p.

## 3. Dangerous clusters are genuine simultaneous first-layer roots

The strongest two-prime cluster by ratio is n=68, with primes 37 and 41. Both arise by reflecting earlier single events: the partners are n=42 for p=37 and n=54 for p=41. The fact that each prime's roots pair does not prevent different primes' events from colliding at the same n.

There are exactly six three-prime clusters in the tested range:

| n | First-layer anomalous primes |
|---:|---|
| 2772 | 2383, 2671, 2689 |
| 9559 | 5741, 6841, 7459 |
| 29158 | 17359, 23327, 26497 |
| 29259 | 21139, 23027, 25343 |
| 40219 | 28393, 29723, 37813 |
| 73219 | 37813, 44959, 53231 |

For all eighteen displayed primes, the established shifted modulo-p² lift criterion gives **v_p(g_n)=1**, exactly equal to v_p(n!) in this top-half range. Thus none supplies a positive excess factor, despite contributing to S(n). The same is true for the ratio-maximizing event (n,p)=(21,17).

The shifted test computes C_0(p)=1, C_1(p)=2p+1 and

\[
C_{r+1}(p)=(2p+2r+1)C_r(p)+(p+r)^2C_{r-1}(p)
\pmod {p^2}.
\]

For an anomalous r=n-p, C_r(p) is divisible by p. Its nonzero residue modulo p² proves valuation one; residue zero would give valuation two by the existing top-half Wronskian cap. The displayed clusters all have the first outcome. This bounded computation uses no large A_n, B_n or gcd. It illustrates that S is an intentionally coarse upper budget, not an observed excess product.

## 4. Structural diagnostics and their limitations

Consider a graph whose vertices are anomalous indices n<=100000, joining n and n'=3p-1-n when both are in range. Internal edges include 21 midpoint loops. After removing loops there are no cycles in this finite graph; its largest connected components have six vertices. For example, one component contains 40219 and 73219, which share the reflected p=37813 event. These are finite graph observations. Even a global forest statement would not bound the degree of a vertex, so it would not alone prove the desired pointwise S(n) estimate.

The 21 midpoint-anomaly primes visible in this range are all congruent to 31 modulo 32. This is recorded as a diagnostic, not proved as a congruence restriction. A midpoint event has p=(2n+1)/3, so there is at most one such prime at any fixed n; controlling it cannot alone bound the remaining anomalies.

Across the entire event list, all four odd residue classes modulo 8 occur: 1744, 1683, 1714, and 1885 events for classes 1, 3, 5, and 7 respectively. A simple global exclusion of any one of these residue classes is therefore contradicted by exact examples. Complete digit-root sets can also have more than a single reflected pair: p=56807 has ten roots. Thus a strategy based on “at most two roots for every prime” is false.

The useful remaining arithmetic task is a uniform bound on collisions between different moduli at a fixed index. A bound on each individual prime's root count, reflection alone, or the observed small cluster size does not establish that bound. No eventual S(n)<=n/6 theorem is claimed here.

## Reproduction and dependencies

Run `python -B scan_global_s_clusters.py --limit 100000` from this directory. The standard-library script performs about 154 million small modular recurrence steps, independently cross-checks the direct enumeration through 10000, and produces `global_s_cluster_checks.json`. A reference run takes about 25 seconds before the small diagnostic lift queries.

The proof inputs are the previously established generating-function digit congruence and top-half shifted modulo-p² lift criterion. Their repository checkpoints are the odd-prime blocks in `pi-gcd-blast-2026-11-02` and `gcd_blocks/PRIME_BLOCKS_AND_LIFTING.md` in `pi-gcd-prime-adic-2026-11-02`. All numerical results above have finite scope; no conclusion about all subsequences of the recurrence family or all rational approximations to pi follows from them.
