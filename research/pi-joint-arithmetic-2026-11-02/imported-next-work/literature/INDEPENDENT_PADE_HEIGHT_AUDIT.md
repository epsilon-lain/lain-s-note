# Independent audit of the Padé benchmark and reflection height

2026-11-02. This agent review checked `PADE_LINEAR_DEGREE_BENCHMARK.md` and `../arithmetic_route/REFLECTION_HEIGHT.md` directly. It is not external peer review.

## Padé benchmark

No mathematical error was found in the stated formulas.

- With the displayed normalization the remainder is `q_m exp(Z)-p_m=(-1)^m Z^(2m+1)/(2m)!` times the beta integral. At omega=2*pi*i the integral is `-B(m+1,m+1)*chi_m`, and `(-1)^m omega^(2m+1)=i*(2*pi)^(2m+1)`. Consequently `g_m(omega)` has positive imaginary sign, and formula (5) gives positive `C_m(omega)`.
- Stirling expansion independently gives the linear terms in (6) and (10): respectively `1+2 log pi` and `2+2 log(pi/2)`.
- Since `g'_m(omega)=-1+O(1/m)` and the root is close, formula (9) has displacement `+g_m(omega)*(1+O(1/m))`, with the stated positive imaginary direction.
- The positive coefficients of `C_m` strictly decrease through their odd indices, so its exact naive height is `(2m-1)!/(m-1)!` and its logarithm has linear term `2 log 2-1`.
- The condition `b^(m-1)*|C_m(omega)|<=1/2`, divided logarithmically by m along `m/b -> c`, has limiting logarithm `log(e*pi^2/c)`. The threshold `e*pi^2` is correct for this particular baseline-value condition.
- The complete coefficient bounds in (7), large and small Rouche disks, and the rational-value denominator are consistent. No irreducibility statement for the odd part follows from irreducibility of the full Padé polynomial; the draft correctly avoids that transfer.

## Reflection-factor height

The exact constant-term argument, entry coefficient-norm bound, and leading height term `(n^2/4)*log n` check out. The height statement concerns the full reflection factors, not necessarily the minimal polynomial of the selected root.

One harmless orientation correction was communicated to the author: after deleting the even coordinate pivots in the plus case, the retained odd-index matrix has diagonal first r-1 rows and a final row that may contain earlier entries. It is lower triangular rather than upper triangular. Its determinant remains `(-1)^(r-1)*(n-1)`.

The benchmark supplies a substantially simpler classical family exhibiting the same leading root-distance scale; it does not constitute a rational approximation theorem or prove bad approximability.
