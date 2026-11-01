# Exact finite optima after quotienting by the forced anchor block

Four previously tested cases now have certified **global minimum cost over every allowed type choice** in the subclass `a >= beta`, `w*a <= H`, `ell=s`, `d=k=0`, with the `j=0` rows fixed. Each optimum survives the anchor quotient with full row rank, and the same monomial columns have a nonzero original rational-center minor. The optimization does not discard any assignment on the basis of a vanishing modular determinant.

These remain finite structural results. They do not verify the manuscript's unknown separated-weight thresholds, establish an all-large-`H` statement, or optimize over terms with `d != 0` and `ell < s`.

## Exact results

All four cases have `K=2`, `w0=1`, `v0=1/2`, `nu=2`, and `c=log(2)/4`. The cost is

\[
C=\frac{2\sum_{\rm rows}w(a-\beta)+c\sum_a n_a^2}{MH}.
\]

| Dimension, weights | H | theta | Rows x columns | Global minimum C | Previous cost, approximately | Minimum, approximately |
|---|---:|---:|---:|---|---:|---:|
| 2, (1,5) | 8 | 1/4 | 48 x 55 | 13/96 + (63/32)c | 0.487612 | 0.476575 |
| 2, (1,5) | 16 | 1/3 | 208 x 262 | 217/1664 + (183/104)c | 0.456648 | 0.435327 |
| 2, (2,5) | 12 | 1/4 | 64 x 73 | 59/384 + (757/384)c | 0.499665 | 0.495255 |
| 3, (1,2,5) | 8 | 1/3 | 68 x 108 | 1/8 + (495/272)c | 0.446621 | 0.440356 |

Every displayed cost strictly improves its previous fixture and is less than `1-bbar`. The original two numerical parameter conditions hold. The fixed anchor ranks are respectively 24, 104, 32, and 34; the selected rows have quotient ranks equal to those same numbers.

## Global optimality without any rank-deficiency pruning

Let `b_a` be the fixed anchor count, and `x_a` the number of free rows assigned type `a`. Relax the problem by discarding **all** conditions on which row receives each type, including `a>=beta`, distinct labels, and rank. Keep only nonnegative integer counts with total equal to the number of free rows, and all types satisfying `wa<=H`.

Apart from a constant independent of the counts, the objective is

\[
2\sum_a wa\,x_a+c\sum_a(b_a+x_a)^2.
\]

The kth unit assigned to a type costs

\[
\Delta_a(k)=2wa+c(2b_a+2k-1).
\]

These marginal costs increase strictly with k. For each of the four certificates, **every selected marginal is strictly below every unselected marginal**. Therefore the certified counts are the unique optimal aggregate counts in the relaxed problem. A feasible full-rank assignment with those counts attains this lower bound, proving global optimality in the actual subclass.

| Case | Number of allowed types | Largest selected marginal | Smallest unselected marginal |
|---|---:|---|---|
| H=8, (1,5) | 13 | 41c | 6+7c |
| H=16, (1,5) | 38 | 4+63c | 12+17c |
| H=12, (2,5) | 13 | 8+17c | 4+41c |
| H=8, (1,2,5) | 31 | 6+9c | 4+21c |

`verify.py` enumerates **every type**, checks every last-selected/first-unselected pair, and verifies each strict comparison at both rational endpoints `14/81` and `25/144`. These endpoints enclose c: the first two positive terms of

\[
c=\frac12\sum_{k\ge0}\frac{1}{(2k+1)3^{2k+1}}
\]

give the lower bound, and replacing all omitted denominators after the first term by 3 gives the upper bound. Thus the comparison needs no floating-point assumption.

The optimizer also supplies an independently checked minimum-flow certificate. Its network has one unit per free source row, one capacity-one slot per `(j,s,a)`, a node for each type, and unit edges with the above convex marginal costs. Stored dual potentials have the form `A+B*c`. The verifier reconstructs every residual arc from the assignment and checks its reduced cost with rational bounds on c, including all reverse arcs. This is a second optimum certificate; the aggregate proof is already sufficient.

Consequently the aggregate optimum, row-feasible optimum, and full-rank optimum **coincide in all four cases**. There is no finite extra cost from rank in these examples.

## Anchor quotient and determinant certificates

Write the test matrix in anchor/free row blocks. Select anchor pivot columns J, and let L be the remaining columns. With the original row order fixed, the anchor block A has invertible `A_J`. The exact Schur projection of a candidate free row v is

\[
\bar v=v_L-v_JA_J^{-1}A_L.
\]

For selected quotient columns B,

\[
\det T_{J\cup B}=\det A_J\det\bar T_B,
\]

where the column order is J followed by B. The implementation keeps the anchor determinant, all pivot indices, and the quotient determinant. The independent verifier rebuilds the matrix, repeats this projection, and checks the factorization against a direct determinant on the same columns.

All arithmetic in this check is exact modulo the verified prime `1,000,000,009`, at formal specialization `W=2`.

| Case | Anchor determinant | Quotient determinant | Product/direct test determinant | Original determinant, same ordered columns |
|---|---:|---:|---:|---:|
| H=8, (1,5) | 56639847 | 510464783 | 947794387 | 918208015 |
| H=16, (1,5) | 372489085 | 385826179 | 91311864 | 590560636 |
| H=12, (2,5) | 462358753 | 783797555 | 872688011 | 473934798 |
| H=8, (1,2,5) | 70634668 | 936751591 | 33252555 | 494028796 |

The formal entries are

\[
\binom\alpha a[t^s](1+t)^h(jW+\log(1+t))^{|\alpha|-|a|}.
\]

Every relevant coefficient denominator is invertible modulo the prime: all logarithm indices are at most 31. A nonzero modular minor therefore certifies a nonzero polynomial over `Q[W]`. Using the **known transcendence of pi**, that same minor is nonzero at `W=2*pi*i`. This use of known transcendence is part of the finite certificate; it is not a new proof of transcendence or an irrationality measure bound.

The original matrices use the genuinely different rational centers from `(p,q)=(6,2),(22,7),(355,113)` for weights 1, 2, 5. Their entries are in `Q(i)`, with `i` mapped to `430477711` modulo the prime. The verifier checks this square root of -1, the original determinants on exactly the same monomial columns, `ceil(log q)=w`, and `|pi-p/q|<q^-2` using rational exponential and Machin-series bounds. Taking `F0=17` gives `Ti=34wi>31`, so the truncated original `Gi` and logarithm have identical required jets. All selected `d=k=0` terms have nonzero scalar factors: at the anchor `a=beta`; elsewhere the rational approximation errors cannot vanish.

## A short explicit witness

`SIMPLE_OPTIMUM.md` gives a hand-written H=8, weights (1,5) assignment with the same global minimum. Its total type counts are `(21,15,9,3)` for types `(0,0),(1,0),(2,0),(3,0)`. The first 48 columns already work. The independently rechecked determinants are `146986659` for its translated test and `918208015` for the original matrix. This distinct witness makes the smallest result directly reviewable without inspecting a network flow output.

## Reproduction

All code uses the Python standard library and is self-contained in this directory:

```sh
python3 optimize.py
python3 verify.py
```

`optimize.py` rebuilds the four optima and Schur certificates in `optimal_fixtures.json`. `verify.py` performs the independent feasibility, exhaustive marginal, residual-dual, determinant, and actual-center checks, then writes `verification.json`. No modular rank-zero observation is used as a characteristic-zero upper bound.
