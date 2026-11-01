# Batch 10 full-reading notes

Pinned repository: `openai/math`, commit `4eaee87548250b2682affda7d05d947c05554bae`.

This batch contains complete first source readings of queue positions 22 and 31 in ascending order. Reading is not independent validation of the manuscripts' theorems. Position 36 was reassigned by the parent and is not included as an unread placeholder.

## Coverage and byte audit

- Position 22: 15 source blobs, 2,946 lines; README, main, nine body sections, two appendices, one included figure and bibliography. All cached blobs match pinned Git blob SHA.
- Position 31: 13 source blobs, 4,446 lines; README, main, seven body sections, three included figures and bibliography. All cached blobs match pinned Git blob SHA.
- Sources cached at `source_context/round02_group4/{22,31}`. All input/include/bibliography edges resolve. Added terminal LF was stripped only when removal exactly matched the expected blob SHA.
- Position 22 main is `build/main.tex`; position 31 main is `build/paper.tex`. Neither README supplies a compile command; the recursive preflight tree plus document/include structure identify the actual mains.

## 22 — A Direct Proof of Optimal Max-Cut Hardness

The manuscript claims unconditional NP-hardness of every approximation factor strictly above the Goemans–Williamson constant, including simple unweighted graphs. The exact weighted gap compares `(1+t)/2` against `(1+B(t))/2`, with `B(t)=2 arcsin(t)/pi` and fixed error. The external inputs are the imperfect-completeness 2-to-1 game construction/Grassmann machinery and Majority Is Stablest. Unique Games is historical context in this text.

The source-game quantifier order is crucial: fix soundness, then an alphabet depending only on soundness, then arbitrary completeness error. Appendix A works through affine folding, validity conditioning, completed tables, invariant lists and tagged quotient lists, rather than treating this order as an informal numerical game promise.

The tree gate makes each Walsh monomial select one child in every pair, while a random slice selects both children of one pair. Their exact-one intersection gives precisely one correlation factor. Odd/even decomposition turns an overly good cut into excess stability. Row symmetries impose common child-pair sums; gate parity makes odd sums actual labels, not arbitrary nonzero functionals. Depth-averaged local unsplit projection costs `s/(2e m mu)` through weighted support-growth telescoping. A vector-valued affine hash then converts diffuse label masses into small maximum table-entry influence, without ambient dimension loss. Its fourth moment forces a common character in all four factors and produces a collision square-sum.

Decoding freezes most rows, plants zero blocks at a uniformly hidden product coordinate, controls the unmarked law by a binomial likelihood ratio, and compares only `2s` free shifts across a codimension-one restriction. The resulting loss is `2^(2s)`, independent of the alphabet. The final parameter order avoids a circular alphabet/soundness choice. Appendix B explicitly uses finite-field dot-product rectangles to remove rational weights while controlling cuts that split vertex clusters; loops are deleted without renormalization.

Tentative interfaces for the pi construction: exact-one transversal intersections; support-surplus telescoping; square-sum collision control for occupancy; restriction costs depending only on a bounded active block. Each requires an actual arithmetic matrix realization, Q(Z) noncancellation and a height budget. None by itself produces a nonzero integer determinant with absolute value below one.

## 31 — A boundary-only obstruction to four-dimensional disk embedding

The manuscript claims a finite family of proper discs with disjoint boundaries and framed algebraic dual spheres in a compact smooth four-manifold, satisfying exact equivariant duality/vanishing conditions, whose boundary circles nevertheless admit no disjoint locally flat filling. The conclusion permits arbitrary filling homotopy classes and framings. It derives nongoodness of F2 by subgroup/quotient closure, without identifying the constructed fundamental group with F2, and excludes the conjunction of unrestricted surgery and s-cobordism.

The finite geometry glues Hopf product balls to a one-handlebody, obtaining a graph homotopy model and a longitude–meridian interchange. All independently conjugated triple commutators are realized by capped stages. Only the base tubes are removed. Product tori compressed twice give individually embedded framed dual candidates; the unused upper cap is retained and kills the exact comparison loop in the exterior, so opposite signs carry equal Wall labels. Hypothetical actual topological disc images permit a graph homotopy supported away from smooth regular cuts, with no smoothing of those discs and no prescribed relative homotopy class.

The two-sheeted Mayer–Vietoris calculation gives complementary wall restriction tangents and a degree-two detection isomorphism. Finite transition-cochain charts retain the full nonreduced defect ideal. A finite adjoint spanning set, using diagonal spectral projectors and two unipotents, converts independently conjugated leading triple commutators into every scalar cubic product; Nakayama upgrades the first-layer vanishing to an exact cube-zero end ideal.

Polarized Chern–Simons wall functions recover that entire defect ideal through their internal derivatives. A missing generator would define a minimal small-extension obstruction invisible to all wall variations; duality and paired H2 then kill its curvature, permitting a forbidden flat lift. A separate closed-boundary gluing argument distinguishes connection changes from gauge changes and applies Stokes to curvature squared, yielding a square-ideal value identity.

All finite gradient jets are controlled in one fixed complex algebraic chart, via finite smooth/algebraic sums and explicit inverses. The algebraic replacement theorem descends a whole Taylor tuple by faithful flatness, then checks every primary component, including embedded ones, using differential powers. Ninth-order replacement preserves exact gradient generators and both cubic memberships. Algebraic witnesses and a nonzero complementary-tangent determinant are placed in one finite coefficient ring. A denominator-free recursive graph correction makes the potential vanish exactly and transports all witnesses by an invertible gradient change. Only this algebraic/formal output is reduced to positive characteristic.

The final obstruction works in a coordinate Frobenius quotient and keeps ordinary differentiation distinct from the wedge Koszul differential. Complementary zero graphs produce a nonzero socle/residue pairing. The negative cubic law and a tensor joint-kernel calculation force enough kernel factors; chain retractions with explicit homotopies make three factors liftable while preserving the pairing. Their deformation over three square-zero parameters retains all eight mixed terms. A moving normal form and the positive cubic law give an explicit ambient primitive for `s1*s2*s3` times that form, in the very same finite complex. Pairing with the actual lifted cycle is a unit, so the boundary identity forces the nonzero squarefree parameter monomial to vanish.

Tentative interfaces for the pi construction: exact ideals and membership witnesses rather than only radicals; whole-tuple descent rather than separately chosen approximations; denominator control by a fixed initial Jacobian; explicitly preserve one nonzero determinant through every correction; all triple annihilations give an at-most-two-active tensor occupancy rule; and prove nonvanishing and annihilation for the same certificate in the same ambient object. These are distant construction principles, not an arithmetic proof. Neither finite-ring existence nor preservation of a localized unit supplies a coefficient-height estimate or a small nonzero integer determinant.

## Verification boundary

Both manuscripts' complete source closures were read. No claim is made that their major new results have been independently verified. No compilation or PDF page audit was performed, and external references were not independently reread. The structured records retain explicit hypotheses, mechanisms, possible interfaces and missing target steps for every extracted idea.
