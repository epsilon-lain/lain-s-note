# Active-prime capacity: a distinct conditional route

This note gives a new sufficient target for excluding bounded-quality
subsequences in this particular approximation family. It does not prove
the target's required prime-density bound. The attached scan is finite
evidence, and the arithmetic classifications in that scan are exact.

## Sequence and available literature

The sequence is

\[
U_n=[x^n](1+2x+2x^2)^n,
\qquad U_0=1,\quad U_1=2,
\]

with recurrence

\[
(n+1)U_{n+1}=2(2n+1)U_n+4nU_{n-1}.
\]

Tony D. Noe, *On the Divisibility of Generalized Central Trinomial
Coefficients*, Journal of Integer Sequences **9** (2006), Article 06.2.7,
[official primary PDF](https://cs.uwaterloo.ca/journals/JIS/VOL9/Noe/noe35.pdf),
identifies this sequence in Table 1 as the case `(a,b,c,d)=(1,2,2,-4)`,
OEIS A006139. Equation (13) proves digit multiplication modulo odd primes;
equation (14) gives the Holt reflection used by the code. Theorem 8.5
states the resulting digit divisibility criterion; Theorem 8.6 gives the
zero symmetry. Noe's discussion immediately before Theorem 8.2 explicitly
leaves the general case `p` not dividing `bd` beyond the digit criterion.
Every odd prime in the present family is in that general case.

Those results do not supply an asymptotic density, or a logarithmic mass
bound, for the set of primes defined below. In particular, a theorem
about the density of divisible *indices for one fixed prime* is different
from a theorem about the distribution of these primes.

## Exact prime-capacity reduction

Define

\[
\mathcal Z=\{p\text{ odd prime}:\exists\,0\le r<p,
\ U_r\equiv0\pmod p\},
\qquad
\vartheta_{\mathcal Z}(x)=\sum_{\substack{p\le x\\p\in\mathcal Z}}\log p.
\]

Lucas digit multiplication proves the equivalence

\[
p\mid U_n\text{ for some }n\ge0
\quad\Longleftrightarrow\quad p\in\mathcal Z.
\]

The forward direction uses a zero digit factor in the base-`p` expansion;
the reverse direction is immediate from the definition. As a consequence,
every prime contributing to the existing multi-band mass `W(n)` belongs
to `Z`, so

\[
\boxed{W(n)\le\vartheta_{\mathcal Z}(n).}
\]

Thus a proof, for some fixed `eta>0`, of the new global target

\[
\boxed{\vartheta_{\mathcal Z}(x)
\le\left(\frac{2+\log2}{4}-\eta\right)x
\quad\text{for all sufficiently large }x}
\]

would combine with the prior multi-band theorem to give

\[
\log Q_n\ge4\eta n+o(n),
\]

and hence exclude every bounded-quality subsequence in this family.
Here `(2+log 2)/4 = 0.6732867951...`. This implication is rigorous; the
boxed prime-capacity target is unproved.

This target avoids needing to synchronize the exceptional prime residues
at the same index `n`. It asks whether enough primes are even eligible.
It is therefore a different attack, but it remains a substantial
cross-prime arithmetic problem.

## Reproducible finite scan

Run from the directory containing this note:

```sh
python scan_active_primes.py --limit 20000
```

For each odd prime `p`, the code evaluates the recurrence exactly modulo
`p` for `0<=r<=(p-1)/2`. Holt reflection

\[
U_{p-1-r}\equiv(-4)^{(p-1)/2-r}U_r\pmod p
\]

then recovers all roots. The full congruence, not only its zero pattern,
is independently checked for every odd prime through 503: 11,066 exact
checks pass. Logarithms enter only the descriptive summaries, never the
prime classifications.

| `x` | Active odd primes through `x` | `theta_Z(x)/x` | `theta_Z(x)/theta_odd(x)` |
|---:|---:|---:|---:|
| 1,000 | 69 / 167 | 0.391316 | 0.409519 |
| 2,000 | 126 / 302 | 0.402987 | 0.415633 |
| 5,000 | 277 / 668 | 0.407169 | 0.414547 |
| 10,000 | 506 / 1,228 | 0.407388 | 0.411698 |
| 20,000 | 915 / 2,261 | 0.399696 | 0.403639 |

The largest observed prefix ratio `theta_Z(x)/x` through 20,000 occurs
at `x=283` and is approximately 0.54771844. Restricting to
`1,000<=x<=20,000`, the largest observed ratio is approximately
0.42261444 at `x=1327`. These are numerical summaries of finite exact
root scans; they establish no eventual inequality.

## No simple residue-class exclusion appears

All four odd classes modulo 8, and all eight odd classes modulo 16,
contain both active and inactive primes in this scan. Here are the
weighted active fractions through 20,000:

| `p mod 8` | Active / all | `theta_Z,class / theta_class` |
|---:|---:|---:|
| 1 | 223 / 556 | 0.396732 |
| 3 | 223 / 571 | 0.392123 |
| 5 | 225 / 569 | 0.398043 |
| 7 | 244 / 565 | 0.427760 |

| `p mod 16` | Active / all | `theta_Z,class / theta_class` |
|---:|---:|---:|
| 1 | 113 / 279 | 0.399995 |
| 3 | 107 / 290 | 0.375296 |
| 5 | 103 / 283 | 0.367388 |
| 7 | 105 / 282 | 0.373525 |
| 9 | 110 / 277 | 0.393440 |
| 11 | 116 / 281 | 0.409507 |
| 13 | 122 / 286 | 0.428279 |
| 15 | 139 / 283 | 0.481673 |

In particular, excluding entire residue classes modulo 8 or 16 cannot
describe the active-prime set. For instance, modulo 16 the first active
primes in the eight odd classes are respectively
`17,131,37,23,41,107,157,31`; the first inactive examples are
`193,3,5,7,89,11,13,79`.

The observed capacity is well below the sufficient threshold, making
this worth targeted theoretical study. A finite root table and an
apparent limiting fraction cannot be promoted to the required uniform
prime-mass theorem. No claim about the bad approximability of pi follows.
