# A finite-resolution countermodel: identical visible prefix, opposite bad-approximation behavior

This is a standard continued-fraction construction, not a counterexample concerning pi itself and not a new theorem about pi. Its purpose is to test what conclusions can follow from finite denominator observations plus an irrationality exponent of 2.

## Precise claim

Given any finite initial continued-fraction segment of pi, any positive tolerance, and any finite list of strict rational threshold/comparison observations, there are two irrational numbers inside that tolerance of pi, extending the prescribed initial segment, such that:

* both have irrationality exponent exactly 2;
* one is badly approximable and the other is not;
* both preserve the chosen finite observations, provided those observations have nonzero margins at pi.

This does not say that their exact distances to rationals are identical: an exact real-valued distance oracle would contain unlimited information. The preserved data are finite prefixes, finite-precision measurements, nearest-grid indices, strict rankings, and strict threshold decisions.

## Constructions and proof

Write the given prefix as [a_0; a_1,...,a_N], and let p_n/q_n be its convergents. Set p_{-1}=1, q_{-1}=0. For a complete quotient T>1 define

F_N(T)=(p_N T+p_{N-1})/(q_N T+q_{N-1}).

The determinant identity |p_N q_{N-1}-p_{N-1}q_N|=1 shows that the open continued-fraction cylinder has endpoints

p_N/q_N and (p_N+p_{N-1})/(q_N+q_{N-1}),

and diameter 1/[q_N(q_N+q_{N-1})]. Every irrational completion lies strictly inside it. Increasing N along pi's actual expansion makes this diameter tend to zero.

For the first target, append only 1s:

x_B=[a_0;...;a_N,1,1,1,...]=F_N((1+sqrt(5))/2).

Let A=max(1,a_1,...,a_N). For an irrational x with complete quotient t_{n+1}=[a_{n+1};a_{n+2},...] and convergent p_n/q_n, direct subtraction gives

K_n(x):=q_n^2 |x-p_n/q_n|=1/(t_{n+1}+q_{n-1}/q_n),

so

1/(a_{n+1}+2)<K_n(x)<1/a_{n+1}.

Thus every convergent of x_B satisfies K_n(x_B)>1/(A+2). Legendre's theorem says that a reduced nonconvergent p/q satisfies q^2|x_B-p/q|>=1/2. Consequently every reduced rational satisfies

q^2|x_B-p/q|>1/(A+2).

Unreduced denominators only increase the left side. Hence x_B is badly approximable. Infinitely many convergents have error <1/q_n^2, and the uniform lower bound prevents every exponent above 2, so mu(x_B)=2.

For the second target, append a sparse sequence of unbounded partial quotients. Number the new positions k=1,2,..., and define

b_k=j+1 if k=2^j for some integer j>=1; b_k=1 otherwise.

Set x_U=[a_0;...;a_N,b_1,b_2,...]. At the spike n_j=N+2^j, take the preceding convergent, with index n_j-1. The same identity yields

1/(j+3)<K_{n_j-1}(x_U)<1/(j+1).

Therefore inf_{p,q}q^2|x_U-p/q|=0: x_U is not badly approximable.

Nevertheless, mu(x_U)=2. Indeed q_n>=F_{n+1}, where F_n is the Fibonacci sequence (apart from irrelevant initial indices), so log q_n grows at least linearly in n. Meanwhile a_{n+1}<=max(A,n+2). For every epsilon>0, eventually a_{n+1}+2<q_n^epsilon. The convergent lower estimate then gives

|x_U-p_n/q_n|>q_n^(-2-epsilon).

Every reduced nonconvergent has error >=1/(2q^2)>q^(-2-epsilon) for sufficiently large q. This excludes infinitely many approximations of every exponent above 2. The infinitely many convergents below q_n^(-2) establish the reverse inequality, so mu(x_U)=2.

There is a quantitative interpretation of the spike subsequence. The recurrence q_n=a_n q_{n-1}+q_{n-2} implies q_n<=(a_n+1)q_{n-1}; combining this with Fibonacci growth gives

log q_{n_j-1}=Theta(2^j),

because the logarithmic extra cost of the first j spikes is O(j log j). Hence

K_{n_j-1}(x_U)=Theta(1/log log q_{n_j-1}).

The quality factor tends to zero, yet the improvement is slower than every fixed power of the denominator. This explicitly separates the endpoint exponent mu=2 from a positive uniform bad-approximation constant.

Even adding the statement "the target is transcendental" does not determine which alternative occurs. To see this without a transcendence theorem, allow every tail digit of the bounded construction to be independently 1 or 2; there are uncountably many distinct completions, all badly approximable and of exponent 2. For the unbounded construction, keep the prescribed spikes and independently choose 1 or 2 at the other positions; again there are uncountably many distinct completions, all not badly approximable and of exponent 2. Algebraic numbers form a countable set, so each family contains uncountably many transcendental targets. This is an existence statement; the explicit all-1 completion above is quadratic, and no explicit transcendence claim about the sparse completion is needed.

## Preserving a finite denominator observation profile

Fix Q>=1. For 1<=q<=Q define the nearest integer P_q(x) to qx, and define

K_q(x)=q|qx-P_q(x)|.

The nearest integer can change only when x=(k+1/2)/q, and the sign of qx-P_q(x) can change only at x=k/q. These boundaries are rational, so none equals pi. There is therefore a neighborhood of pi in which every P_q and every sign is fixed.

On that neighborhood, each K_q is affine with slope either q^2 or -q^2. For distinct q and s, the equality K_q(x)=K_s(x) has at most one solution, which is rational. Thus pi is not a tie, and the strict ranking of all Q quality values is constant on a sufficiently small neighborhood.

Likewise a comparison K_q(x)<c with a fixed rational c has a rational equality boundary. Any finite collection of these threshold decisions is locally constant at pi. For finite-precision values rather than decisions, the global Lipschitz estimate

|K_q(x)-K_q(pi)|<=q^2|x-pi|<=Q^2|x-pi|

provides an explicit error bound. Exact agreement of arbitrary real-valued measurements is not claimed.

Choose a sufficiently long actual pi prefix so that its whole cylinder lies inside this neighborhood, inside the desired decimal-prefix cell, and inside the requested tolerance. Append the two tails above. Both targets preserve all these observations, while their bad-approximation behavior differs.

For an explicit neighborhood certificate, first fix the nearest numerators and signs at pi. Take a positive number smaller than all relevant distances to nearest-integer switch boundaries and rational-grid crossings; also require it to be smaller than

|K_q(pi)-K_s(pi)|/(q^2+s^2)

for every compared pair, and smaller than |K_q(pi)-c|/q^2 for each threshold comparison. All these margins are positive. Certified rational intervals for pi can verify a rational lower bound for this minimum.

## What this rules out, and what it does not

It rules out a general deduction of bad approximability from an irrationality exponent of 2 plus any finite initial-resolution experiment. Increasing a fixed finite denominator cutoff cannot repair that logical implication. A proof about pi needs a property of pi that controls arbitrarily remote tails.

It does not prove that pi is or is not badly approximable, does not rebut a pi-specific theorem, and does not imply that finite experiments are useless for discovering such a property.

## Transferring the project's existing large finite certificate

The existing `research/irrational-numerator-blast-2026-11-01/results_q_10pow5000.json` contains 11747 certified initial partial quotients of pi and a certificate that

q^2|pi-p/q|>1/21000 for every reduced rational with 1<=q<=10^5000.

The final denominator of this prefix has 5998 digits. The denominator of its cylinder diameter has 11995 digits. Exact integer comparison verifies that the cylinder diameter is smaller than 10^(-11994). Appending either tail above therefore gives targets satisfying

|x_B-pi|<10^(-11994),  |x_U-pi|<10^(-11994).

For every q<=10^5000, including all choices of numerator, the corresponding normalized error differs by less than 10^(-1994):

|q^2|x-p/q|-q^2|pi-p/q||<=q^2|x-pi|<10^(-1994).

This is approximate agreement, not equality of real-valued error measurements.

Furthermore, every rational r with reduced denominator at most 10^5000 lies at distance greater than 1/(21000*10^10000) from pi. Since the whole cylinder diameter is smaller than this margin, neither replacement crosses any such rational r. Therefore all left/right comparisons against these rationals are exactly preserved. All nearest-grid numerators are preserved for q<=10^5000/2, because their switch points have denominator at most 2q.

For an especially strong finite observation statement, take q<=10^2498. The rational boundaries for pairwise normalized-quality comparisons have denominators at most 2q^2; the boundaries for thresholds 1/2, 1/10, 1/100, and 1/1000 have denominators at most 1000q^2. These are all below 10^5000. Consequently the entire strict quality ranking, all nearest-grid numerators, and all four threshold decision profiles up to q=10^2498 are preserved by both completions. This follows from the prior global finite certificate and the cylinder margin; it requires no enumeration up to this enormous cutoff.

The supplied verifier checks the new integer inequalities and the copied prefix metadata. It does not repeat the independent verification of the previous pi certificate; that is an explicit dependency.

## Primary references

The determinant/error identities and Legendre theorem are standard. See Yann Bugeaud, *Continued fractions*, CIMPA Bamako 2010 lecture notes: https://irma.math.unistra.fr/~bugeaud/travaux/Bamako2.pdf . The same error identities and Legendre criterion are also derived in Ben Lynn's Stanford notes: https://crypto.stanford.edu/pbc/notes/contfrac/converge.html . All tail choices and finite-observation arguments used above are stated explicitly here.
