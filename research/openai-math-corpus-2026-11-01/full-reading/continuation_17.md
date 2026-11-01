# Queue 52 — continuation 17

Pinned commit: `4eaee87548250b2682affda7d05d947c05554bae`. Direct fresh reading: **4 files, 2,333 lines, 174,545 bytes**. This is scoped reading coverage, not independent theorem validation.

| Assigned source | Git blob SHA | Actual full range | Untruncated reading chunks |
|---|---|---|---|
| `build/sections/12-alternation-imaginary.tex` | `ae1c101994fabbc64496a4db7c09cbd33b88583e` | 1–411 | 1–150, 151–290, 291–411 |
| `build/sections/13-real.tex` | `0011d7030ae2e6ec2efa7c46a43f68cf07c3a310` | 1–1101 | 1–220, 221–440, 441–650, 651–850, 851–1020, 1021–1101 |
| `build/sections/14-cyclic.tex` | `d35821150a3a2db475588a41cedb24ee52f8486d` | 1–692 | 1–240, 241–460, 461–692 |
| `build/sections/15-conclusion.tex` | `639c81ea0ea05c862047608a63105bddb05406b5` | 1–129 | 1–129 |

All Git blob SHAs, exact byte counts and line counts match the cache closure and batch 29 base receipt. The JSON records absolute cache paths, SHA-256 checksums and the exact actual-read ranges.

Entire 12-alternation-imaginary, 1–411: conjugate central extension and scalar Bockstein alternation including quadratic ramification; elementary-divisor parity; imaginary-S3 CM comparator and simultaneous Selmer seed; exact central finite/singular switches and common descent multiplier; compatible residual/characteristic-zero evaluations and nonzero traces; central mixed coefficient, integral quotient and ordinary rank-one center.

Entire 13-real, 1–1101: real quadratic ring-class level; elliptic Selmer seed; real-circle Hecke detector; integral rank-two Hecke packet with multiplicities; residual transversality, integral alternating presentation and Pfaffian cofactor lift; repeated-factor divisibility; coefficient-order descent; paired/Pfaffian division; horizontal logarithms and separate integral input-clearing; Heegner identity, denominator intersection and detector primitivity; exact rank-one product anchor.

Entire 14-cyclic, 1–692: cyclic-cubic Selmer seed with F4-linear evaluation prescriptions; Eisenstein comparator; addressed conductor variables and integral curve divisibility; primitive ray-unit determinant from integral regular-to-primitive quotient sequences and regulator/class-number comparison; exact unramified logarithm lattice; Eisenstein primitive and complete order-to-ray label bijection; quotient unit, rank-one center and exact trace-factor cancellation.

Entire 15-conclusion, 1–129: exhaustive residual anchor classification; nonnegativity plus zero sum forces both twist discrepancies to vanish, and the positive twist comparison transfers equality to the original curve; direct CM branch; usual low-rank consequences and exact two-primary quotient; signed squarefree density corollary with separate manuscript dependency.

## Conditional construction ideas

**Retain repeated factors by a weighted Pfaffian cofactor identity.**

Source anchors: `build/sections/13-real.tex:374-468,473-506,508-629`.

An integral alternating presentation supplies a Pfaffian adjugate lift. For A=Delta B Delta, v=Delta b and L=H Delta, every Pfaffian matching uses each index once and Delta #A Delta=(product d_i)#B. The full product of auxiliary traces divides the projected lift even when traces share factors; regularity and coprimality are unnecessary for the polynomial identity.

Possible interface: If an arithmetic approximation matrix admits a proved alternating enlargement, expose factors on rows, the input and the output projection in one common integral basis, then use the identity to divide the exact witness by the full repeated product without losing multiplicity.

Required premises: An actual alternating enlargement with integral comparison maps; all three Delta factorizations, including double factors on entries between weighted indices; exact relation of the resulting cofactor or Pfaffian to the intended integer determinant; independently established nonzero specialization.

Missing bridge: The source provides Selmer-complex and compact-cup inputs, not a pi approximation matrix. It proves neither row rank over Q(Z) nor nonvanishing of a proposed pi-dependent integer determinant. Its factor extraction provides no archimedean smallness, column-height bound or quadratic-occupancy estimate for that same determinant.

**Combine horizontal integrality with a residually nonzero clearing series.**

Source anchors: `build/sections/13-real.tex:651-696,852-938,978-1051,1058-1076`.

Distinguished Weierstrass division makes coefficient-order quotients inject into their integer components. Separately the global cofactor class lies in N[1/2] and c times that class is integral with c nonzero modulo the dyadic uniformizer. Coordinatewise residual injectivity gives N[1/uniformizer] intersect N[1/c]=N, removing the final dyadic denominator.

Possible interface: For a parameterized determinant family, distinguish divisibility on coefficient components from descent to the original lattice. Seek two independent denominator controls and a genuine residual nonzero multiplier before claiming an integral specialized witness.

Required premises: A fixed finite free common lattice; a distinguished denominator or proved quotient-injection theorem; integrality away from the designated prime; a clearing scalar with nonzero reduction; compatibility of every specialization and basis.

Missing bridge: No common lattice or coefficient family for a pi construction is furnished. Knowing formal divisibility, or changing bases separately at each specialization, cannot establish the nonzero integer and analytic estimate needed for BA/non-BA. Rank over Q(Z), column heights and quadratic occupancy remain separate obligations.

**Prove primitivity using an independent detector or a fully normalized comparator.**

Source anchors: `build/sections/13-real.tex:134-263,290-302,704-775,1053-1092; build/sections/14-cyclic.tex:275-433,446-527,577-629,648-658`.

In the real case a weighted real-circle functional detects the Heegner trace nontrivially; an integral factorization Y=R_b(0)p rules out a nonunit quotient. In the cyclic case full unnormalized ray-unit sums, primitive character quotient lattices, exact logarithm indices and a one-global-sign label bijection determine the comparator determinant valuation; residual congruences transfer a unit only after integral curve divisibility is proved.

Possible interface: Build a primitive arithmetic comparator or an independent residue functional that certifies a proposed witness before transferring it to the target family. Preserve full sums, lattice indices and label multiplicities in the exact comparison, with integrality established before cancellation.

Required premises: An independent nonzero detector or explicitly primitive comparator; an integral factorization through a single common basis; literal congruence and determinant identities; exact multiplicities and index/normalization factors; separate target nonvanishing and analytic bounds.

Missing bridge: The component functional, ray-class units, class-number formula and CM label bijection are specialized inputs. They supply no detector for a pi approximation determinant and no rational-approximation inequality. Transferring a residual unit does not bound the absolute value of the same nonzero integer determinant; Q(Z) row rank, controlled common bases, column height and quadratic occupancy are absent.

## Scope and limits

- This is a complete first reading of the four explicitly assigned source files, not independent theorem validation. No assertion here certifies the central-extension, Pfaffian transport, cochain/limit, local logarithm, regulator, CM label, or cited classical and collection-manuscript arguments.
- This reader directly covered only sections 12–15: 4 files, 2333 lines, 174545 bytes. The parent base receipt covers README/root/references and 01–07; a different reader owns 08–11. Those receipts are not silently merged here, and this receipt alone does not establish whole-manuscript full-body status.
- The theorem scope is usual two-power Selmer corank at most one and the exact prime-two leading-term formula, with total real periods, full torsion-containing indices and the stated height convention. The density claim has an external collection-manuscript dependency.
- The BA/non-BA status of pi remains unresolved. BA would require a uniform positive lower bound for q^2|pi-p/q| for all integers p and q>0; non-BA would require a sequence of such quantities tending to zero. None of the read mechanisms supplies either.
- Outstanding transfer obligations are nonvanishing of an exact integer determinant versus analytic smallness of that same determinant, row rank over Q(Z), a common basis/lattice with controlled denominators, and quantitative column-height plus quadratic-occupancy bounds. No source analogy fills these gaps.

Fresh continuation for queue position 52, sections 12–15 only. Every line was actually read through the 13 listed untruncated chunks; exact Git blob SHAs, byte counts and line counts were recomputed against source_closure.json and batch_29/records.json. The full 721-manuscript nonselective queue remains the parent task. No source, base receipt, ledger or publication was changed.
