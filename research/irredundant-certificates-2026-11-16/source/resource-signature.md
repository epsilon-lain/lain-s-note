# Certified macro path signatures (native elementary audit)
External dependency: external/macro-expansion-interface@1. This extends its exact-expansion interface with explicit resource accounting; not a literature novelty claim or new AC shortest-distance claim.

Let G be a primitive directed state graph with state resource r(v)>=0; states encode all primitive-legality information. An instantiated certified macro is a fixed legal finite walk P=(v0,...,vk), with bound source/target and every primitive step verified.
Define S(P)=(k,H,Q), H=max_{0<=i<=k}r(vi), Q=sum_{i=1}^k r(vi).
The empty walk at v has S=(0,r(v),0). Endpoint types remain part of the macro interface.
For matching endpoints define (k,H,Q)*(l,J,R)=(k+l,max(H,J),Q+R).

Theorem. S(P concatenated with T)=S(P)*S(T), associatively. If every primitive edge remains available, augmenting by any such exact certified macros preserves, for every source/target, the entire attainable signature set, hence its Pareto frontier and all optimization/feasibility problems depending only on (expanded steps,peak,post-state work). Macro-only subsets yield soundness, not necessarily equality.

Proof. Concatenation counts the junction once: steps add, all visited resources combine by max, and the post-state lists concatenate without duplication. Associativity follows from that of addition and maximum. Inductively expand every macro in an augmented walk; bound matching endpoints give a primitive walk with the identical signature. Conversely use each primitive edge as itself in the augmented graph. These two maps give both attainable-set inclusions, including unreachable endpoints. Finite heuristics are not asserted to enumerate these sets.

For work including initial state, W(P)=r(v0)+Q(P), so W(PT)=W(P)+W(T)-r(junction). Thus a fixed-source AC path may use Q internally and add the initial resource once. This covers the convention reported from the external ac-08491 archived work49 path, without claiming official-core validation.

At a fixed endpoint, componentwise dominance of (k,H,Q) is stable under identical future extensions, since max and addition are monotone. This supports dropping dominated labels only for upper resource bounds and nondecreasing objectives, with future availability determined by that endpoint and the tracked budgets. No assertion for arbitrary lower/exact bounds or omitted history-dependent guards.

Endpoint-only counterexample: path 0->1->2 with resource (0,2,0) has S=(2,2,2). At peak budget0 the path is forbidden but both endpoints pass; a macro checking only endpoints would illegally bypass the middle state. Endpoint resource alone also fails to determine work.
