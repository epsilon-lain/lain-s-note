# Exact constant terms and sharp logarithmic height scale

Research scratch note, 2026-11-02. Inputs are the reflection identity, exact even coefficient denominator, monic integer factors, and nearest-root scale in epsilon-lain/lain-s-note at befa7cf7a515a18c3bd00ef17a4c6f53dd9ae6df. This note does not prove irreducibility or a rational approximation theorem.

Let n=2r>=4, A_j=product_{a=0}^{j-1}a!, B_n=A_n/A_r^2, and

L_±(Z)=det([z^s](Z/2+z)^a(e^{z/2}±(-1)^a e^{-z/2}))_{0<=s,a<n}.

The established normalization is F_+=(-1)^{r(r-1)/2}B_n Z^{-r}L_+ and F_-=(-1)^{r(r+1)/2}B_n Z^{-r}L_-. These are monic even integer polynomials of degree Δ=r(r-1).

## 1. Exact constant terms

Claim: [Z^r]L_+=(-1)^{r-1}(n-1), and [Z^r]L_-=(-1)^r. Consequently |F_+(0)|=(n-1)B_n and |F_-(0)|=B_n.

For L_+, every odd row vanishes at Z=0; for L_-, every even row vanishes. Each has exactly r such rows. To compute coefficient Z^r, differentiate each vanishing row once and leave the others at Z=0. No factorial r! appears: this is coefficient extraction from a multilinear determinant, each row contributing its linear term.

At a nonvanishing row s the entries are 2*2^{a-s}/(s-a)! (zero if a>s). At a differentiated vanishing row they are a*2^{a-s-1}/(s-a+1)! (zero if a>s+1). After extracting 2^{-s}/s! from each row and 2^a a! from each column, all total factorials and powers of two cancel. The remaining nonvanishing rows are 2*binom(s,a); differentiated rows are (1/2)*binom(s,a-1). The r factors of 2 and r factors of 1/2 also cancel. Thus the desired coefficients are determinants Q_± with rows binom(s,a) or binom(s,a-1), chosen by the parity above.

Let P=(binom(s,a)) be the unit lower triangular Pascal matrix and J the upper shift, J_{k,k+1}=1. Original rows of QP^{-1} are coordinate rows e_s. Shifted rows are those of C=PJP^{-1}. For s<n-1, Pascal's identity gives row C_s=-e_s+e_{s+1}; for the last row C_{n-1,n-1}=n-1. The last assertion follows directly from P^{-1}_{a,j}=(-1)^{a-j}binom(a,j): only k=n-2 contributes to the final entry and gives binom(n-1,n-2)=n-1.

For Q_- all shifted rows have even index and therefore avoid the last row. Deleting odd coordinate pivots gives diagonal entries -1, hence det Q_-=(-1)^r. For Q_+, deleting even coordinate pivots gives a lower triangular odd-index submatrix: its first r-1 rows have only diagonal entries -1, while its final row may have earlier entries and has final diagonal entry n-1. Hence det Q_+=(-1)^{r-1}(n-1). This proves the claim.

## 2. Coefficient and Mahler heights

For a polynomial f write ∥f∥_1=sum|coefficients| and H(f)=max|coefficients|. Every reflection matrix entry satisfies ∥entry∥_1<=2(3/2)^a: sum its explicit binomial expansion and use 2^{-(s-k)}/(s-k)!<=1. Determinant expansion and submultiplicativity of ∥·∥_1 give

∥F_±∥_1 <= B_n n! 2^n (3/2)^{n(n-1)/2}.

For a monic polynomial the Mahler measure is M(f)=product_roots max(1,|α|). Hence M(f)>=|f(0)| and M(f)<=∥f∥_2<=∥f∥_1. Likewise H(f)>=|f(0)|. Combining these bounds with Section 1 yields, for either sign,

log H(F_±)=log B_n+O(n²),  log M(F_±)=log B_n+O(n²).

Stirling summation gives log A_n=(n²/2)log n-(3/4)n²+O(n log n), so

log B_n=(n²/4)log n+((log2)/4-3/8)n²+O(n log n).

In particular both log heights have the sharp common leading term (n²/4)log n.

## 3. Consequence for the geometric nearest root

The nearest nonzero root belongs to one of F_±. Using the already established log δ_n=-2n log n+O(n), the height of that full reflection factor satisfies

(-log δ_n)/(log H(F_±))=(8+o(1))/n ->0.

Thus the spectacular absolute root distance is not a positive fixed power saving relative to the full factor's coefficient height. The same holds with Mahler measure. This precisely quantifies one arithmetic cost in this family; it is not a universal obstruction, and it does not determine the height or degree of the irreducible factor containing the nearest root. An undiscovered much smaller irreducible factor is the logical escape route.

The familiar monic-factor rational denominator remains b^d for any monic even integer factor of degree d evaluated at ia/b, gcd(a,b)=1. Merely rewriting a factor as a polynomial in Z² does not improve the b exponent.
