# 第52项初读：正文01–07与完整源码映射

本笔记由已完成初读的原始收据转写，保留该阅读分支的范围；全稿的08–15由continuation_16和continuation_17补齐，最终状态以合并台账为准。

## 原稿主张

[
  "Theorem1.1 as stated in01-introduction: for any elliptic curve E/Q with s_2(E)=corank_Z2 Sel_(2∞)(E/Q)≤1 using the usual local Kummer conditions including the real place, rank E(Q)=ord_(s=1)L(E,s)=s_2(E), the whole Tate–Shafarevich group is finite, and the positive rational leading-term quotient Q_E=L^(r)(E,1)|E(Q)_tors|²/(r! Ω_E Reg_E∏c_l(E)) has v_2(Q_E)=v_2(#Sha). There is no restriction on E[2] or reduction at two. Ω_E is the total real period and the regulator uses the specified Poincaré height diagonal H.",
  "The any-prime corollary additionally invokes the collection manuscript The Selmer converse for elliptic curves at every prime to reduce s_p≤1 to s_2≤1. The asserted density-one signed squarefree twist consequence, with ranks0 and1 each of density1/2, additionally invokes the separate Goldfeld/Selmer-distribution manuscript. Those references are dependencies, not independently checked here.",
  "Positive-comparison result read in full here: for non-CM A with analytic rank≤1 the discrepancy X(A)≥0. A positive allowed twist with analytic rank≤1 and X=0 transfers a central determinant unit; under s_2(A)≤1 this forces the corresponding analytic rank and X(A)=0. Split-pair anchors and the direct CM comparison are stated in the introduction but their full proofs lie in the still-unread08–15 modules."
]

## 证明架构

- Fix actual integral modular/Betti/Tate lattices and all period, real-component, torsion, Manin and Tamagawa factors. Kill cusp differences geometrically to obtain an integral relative-symbol lattice with primitive evaluation, rather than using a denominator-bearing rational Hecke projector. Construct smoothed classes by unnormalized integral Siegel-unit traces; use full-level reciprocity, cusp residues and archimedean unfolding to identify their rational exponential coordinates.
- Retain arithmetic maps, classes, cup products, nullhomotopies and local conditions in finite free models; pass moving conductor systems to matrix limits at all fixed Artin precisions. Establish exact dualities and determinant volumes, explicitly separating ordinary positive real cohomology over Q from modified-real Poitou–Tate and from imaginary-base duality.
- Build U over Z_2[[t,u_1,…,u_k]]. Residual auxiliary-prime/Legendre-bit tests concentrate positive cohomology at the divisor(2), controlling its vertical valuation. Characteristic-zero finite/singular line switches preserve determinant volumes and lower surplus fiber dimension, controlling horizontal primes away from fixed factors. Rectangular systems with active characters and bounded-series pole/Weierstrass tests remove remaining fixed divisors. Normality then yields integral nonzero U. Common residual twist data transfer closed-point unit status only after integrality is proved.
- Compute the center in analytic rank0 with compact/discrete Poitou–Tate and exact local log/Haar/real volumes. For analytic rank1 first establish the tame logarithm-squared identity: choose horizontal characters detecting a self-cup primitive, compare cup reciprocity with a theta–Eisenstein spectral/intersection height calculation, keep every field-degree and lattice scalar, then use exact Gross–Zagier only at the stage where rank1 is known.
- For an initially unknown s_2=1 center, use the compact Selmer line and its scalar self-cup primitive before asserting a rational spanning point. Choose a direction with eta≠0; a hypothetical first-order Selmer lift contradicts the global sum-of-local-invariants law, giving ord_v δ=1. The transferred unit and local-condition determinant identity force a simple dual-exponential zero. The spectral implication torsion Heegner trace⇒zero derivative, together with Gross–Zagier and odd sign, gives a nonzero analytic derivative. Only then apply the classical low-rank theorem and central equality.
- The introduction states the final non-CM assembly: a split-pair anchor gives X(A^h)+X(A^(hk))=0, previously proved nonnegativity makes both zero, and the positive member transfers the unit to A. Separate reducible/imaginary-S3/real-S3/cyclic and CM branches are mapped but remain for the fresh reader; this receipt does not substitute the introduction for reading their proofs.

## 条件构造线索

### 1. Separate divisorwise integrality, closed-point primitivity, and generic nonvanishing.

mechanism：Separate divisorwise integrality, closed-point primitivity, and generic nonvanishing.

source_anchor：build/sections/01-introduction.tex:189-202; build/sections/05-positive-determinant.tex:34-77,237-259,424-497; build/sections/07-positive-completion.tex:87-106

source_construction：The source controls(2) by a residual concentration argument, horizontal primes by exact switches and bounded character division, obtains U in a normal integral power-series ring, and separately uses a known central valuation0 plus common residual twist reductions to prove U is a unit. Nonnegative divisor orders alone do not make the central value a unit.

possible_interface：Audit a candidate integer-determinant construction in separate stages: exact integrality of the chosen witness, an independent nonvanishing or primitivity certificate, and analytic smallness of that literal witness.

required_premises：An explicitly identified integer determinant in one common basis, with all denominator clearing and valuations controlled; exact correspondence between every algebraic specialization and the analytic expression.

missing_bridge：This dyadic determinant line is not a determinant of a pi approximation family. Residual units or generic nonzero power series provide no row-rank theorem over Q(Z), no archimedean smallness, and no lower bound for the same pi-dependent integer. A new arithmetic family and common-basis identity are required.

status：conditional proof-design interface; no established transfer to pi

### 2. Preserve exact determinant volume while switching local conditions and lowering a finite obstruction dimension.

mechanism：Preserve exact determinant volume while switching local conditions and lowering a finite obstruction dimension.

source_anchor：build/sections/03-foundations.tex:207-319; build/sections/05-positive-determinant.tex:331-421

source_construction：Complementary finite/transverse conditions are organized by exact triangles and duality. A line switch or square switch carries a distinguished class via a unit-volume local identity, preserves the relevant determinant valuation, and reduces a surplus fiber dimension. Iteration ends at a free line where integrality is visible. The derivative comparison includes a depth margin so fixed clearing factors are harmless only at horizontal DVRs where2 is invertible.

possible_interface：Seek exact elementary transformations of an arithmetic matrix that remove obstruction directions while retaining the same determinant up to explicitly controlled integer units, then certify rank at the terminal form.

required_premises：Exact local-to-global identities; a proved dimension decrease; nonzero evaluation minors; one common rational/integer lattice with unit-volume transformations and uniform control of every clearing factor.

missing_bridge：Cohomological duality, Frobenius selection and local finite/singular coordinates are specialized arithmetic inputs. A combinatorial occupancy count does not supply their analogue, and this source gives no candidate row space over Q(Z), column-height estimate, or quadratic-occupancy bound for pi.

status：conditional rank/volume construction interface

### 3. Use bounded-series tests to turn many characterwise divisibilities into one integral division identity.

mechanism：Use bounded-series tests to turn many characterwise divisibilities into one integral division identity.

source_anchor：build/sections/03-foundations.tex:360-416; build/sections/05-positive-determinant.tex:441-497

source_construction：Bounded coefficients permit Weierstrass division and uniqueness from sufficiently high two-power characters. The pole test uses integral numerator and denominator convergence on smaller disks to preserve zero multiplicities even without a uniform bound on quotient denominators. Distinguished fixed factors and enough admissible active character tests then force the division remainder to vanish.

possible_interface：For a parameterized exact determinant identity, first place numerator and denominator in a coefficient ring with one common quantitative bound; choose a test set that determines a bounded remainder; prove the identity before specializing the target.

required_premises：A genuine bounded family, an exact algebraic divisor with the required regularity, test specializations avoiding only proved nonidentical conditions, and a uniqueness/division theorem in the relevant coefficient ring.

missing_bridge：Dyadic character uniqueness is not an archimedean estimate at pi. Without a bounded parameter family and a literal determinant identity it yields neither nonvanishing nor analytic smallness of an integer. Independent height, rank and occupancy estimates are still missing.

status：conditional interpolation/identity interface

### 4. Detect a first derivative using an obstruction pairing before identifying the class with a known geometric generator.

mechanism：Detect a first derivative using an obstruction pairing before identifying the class with a known geometric generator.

source_anchor：build/sections/06-tame.tex:148-238,1289-1298; build/sections/07-positive-completion.tex:124-193,198-403

source_construction：The compact Selmer class has a scalar self-cup primitive without a rational point being known. A direction is chosen where eta=k(g)-β(g)g_T(g_T-1)^-1x(g)≠0. A first-order lift would violate cup reciprocity, so the Bockstein has order1. The determinant unit then forces a simple singular-localization zero, and an independently proved spectral torsion implication detects the analytic derivative.

possible_interface：In an existence/nonvanishing strategy, identify an independent first-order obstruction that prevents a chosen witness from vanishing to excessive order before assuming a particular basis or geometric representative exists.

required_premises：A nonzero class established independently; an exact deformation and obstruction identity; an admissible parameter direction satisfying all constraints; a separate bridge from the obstruction derivative to the same arithmetic witness.

missing_bridge：There is no corresponding cohomological class, self-cup primitive, reciprocity law or derivative-to-integer-determinant identity for pi in this source. It does not resolve rational approximation, and using the desired nonzero determinant as the assumed class would be circular.

status：conditional noncircularity template only

## 范围限制

- This receipt reports complete first reading of the explicitly marked source ranges, not independent theorem or dependency validation. The new contraction/limit, residual concentration, derivative, intersection and Bockstein arguments and all cited classical inputs have not been independently proved or computationally certified here.
- The complete compiled source closure plus README is18 text files, each cached and exact Git-blob-SHA verified. This lead actually read10 files/4526 lines. Sections08–15 remain wholly unread in this base receipt, with precise full-file intervals retained. Downloads and source mapping are not reading coverage.
- The theorem concerns the two-primary leading-term formula in usual Selmer corank≤1. It does not assert higher-corank or other prime-part leading-term formulas. The any-prime input criterion and density corollary have separately cited collection-manuscript dependencies. The ordinary/modified real terms, total versus connected periods, factors2, full torsion-containing indices and height convention are material.
- The positive comparison is read in full here; split-pair anchors, residual-image branches, direct CM comparison and final assembly are not covered beyond their introductory statements. These unread modules must be read through a continuation before any full-body label is supported.
- No construction resolving the BA/non-BA status of pi is supplied. BA would require a uniform positive lower bound for q²|pi-p/q| for all integer p and q>0; non-BA would require a sequence with this quantity tending to0. Both remain unresolved in this work.
- Outstanding transfer obligations remain exact nonzero integer determinant versus analytic smallness of that same determinant, row rank over Q(Z), a common basis/lattice with controlled denominators, and quantitative column-height plus quadratic-occupancy bounds. The source dyadic determinant volumes, residual-unit transfers and cohomological rank tests do not establish these for pi.

Lead basepartial receipt at queue52, following the nonselective721-manuscript queue. README/root/reference list and01–07 read in full without truncation;08–15 are mapped/cached for a fresh reader. No publishing or ledger changes.
