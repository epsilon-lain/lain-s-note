## Fixed object and action language
State is an ordered pair of freely reduced words on signed integers ±1,±2. Start s=((-2,-1),(-2,-2,-1,-1,2,2,1)); target t=((1,),(2,)).
IDs 0,1 invert r0,r1. IDs 2,3 replace r0 by red(r0 r1), red(r0 inverse(r1)); IDs 4,5 replace r1 by red(r1 r0), red(r1 inverse(r0)). IDs 6..9 conjugate r0 by g=(1,-1,2,-2); IDs 10..13 conjugate r1 by the same ordered g, using red(g r_i (-g)). Inverse IDs are (0,1,3,2,5,4,7,6,9,8,11,10,13,12). No sorting, cyclic reduction, conjugacy quotient or length/resource pruning.
L(s) is total word length. For trajectory (s0,...,sn), peak=max_i L(si), work=sum_i L(si), including both endpoints and initial state.

## Actual finite native computation
Complete forward B4 layers: (1,13,129,1186,10656), total11985, 18606 generated edges.
Complete target B4 layers: (1,10,77,574,4453), total5115, 9268 generated edges.
All 27874 edges checked against a second explicitly cased reduction implementation and their inverses checked in both implementations. These are author's internal crosschecks, not independent peer review and not official source/CLI verification.
B3(s)∩B4(t)=∅. B4(s)∩B4(t) consists of 2 states, both at distance4 from each end. All shortest-parent chains retained, including distinct action IDs with the same resulting state.
Exactly 8 labelled 8-step sequences; for each m∈{7,8} they are:
A: (m,4,1,12,12,2,5,5), lengths (9,9,7,7,5,5,4,3,2), peak9, work51.
B: (m,4,12,1,12,2,5,5), lengths (9,9,7,5,5,5,4,3,2), peak9, work49.
C: (m,4,12,12,1,2,5,5), lengths (9,9,7,5,5,5,4,3,2), peak9, work49.
D: (m,4,12,12,3,1,5,5), lengths (9,9,7,5,5,4,4,3,2), peak9, work48.
The geodesic DAG has layer sizes (1,1,1,2,2,2,1,1,1) and15 labelled edges.
Thus at raw length8 the complete peak/work distribution is (9,48):2, (9,49):4, (9,51):2. Its Pareto frontier is the single pair(9,48), attained by D. The archived official witness B with m7 has work49, so D improves this declared work by1 without changing raw length or peak.
One D witness's exact states are:
((-2,-1),(-2,-2,-1,-1,2,2,1));
((-1,-2),(-2,-2,-1,-1,2,2,1));
((-1,-2),(-2,-2,-1,-1,2));
((-1,-2),(-2,-1,-1));
((-1,-2),(-1,-1,-2));
((1,),(-1,-1,-2));
((1,),(2,1,1));
((1,),(2,1));
((1,),(2,)).

## Why enumeration is complete
Every <=7-step path would have a vertex within3 of s and4 of t (use the endpoint if shorter than3); the complete empty intersection rules this out. Any8-step path is consequently geodesic: all prefixes and suffixes are shortest, otherwise replacing one would shorten the whole path. Its step4 midpoint belongs to both complete B4 balls. Conversely, every shortest forward parent chain and reversed inverse target parent chain at a common midpoint concatenates to a legal8-step path. Enumerating all such labelled chains covers all8-step sequences. A separate DAG path count agrees. Historical distance8 is external_seed; this is a new small-ball reproduction plus new complete resource refinement, not a first discovery of that distance.

## Precisely limited macro-description result
Language: output is an ordered list of the eight path strings A7,B7,C7,D7,A8,B8,C8,D8. It may use zero or one flat macro M with a body of2..4 primitive IDs. Primitive occurrence and M reference each cost1; defining M costs its body length; name/rule header overhead is0. No nesting, parameters, output deduplication, cross-output continuation, inverse abbreviations or automatic state equality rewrites.
All55 distinct length2..4 contiguous substrings are candidates. A useful macro must be one of them; a suffix DP considers literal or matching M at every position, so tokenization is exact. No-macro is explicitly allowed.
Every single path separately has optimum8 including definition costs: no improvement.
For the whole corpus, optimum55 instead of64. A minimizer M=(2,5,5) costs3 and occurs once in each A/B/C (six uses); these six outputs cost6 each and two D outputs remain8, total3+6*6+2*8=55. Candidate enumeration and DP establish the lower bound only in this declared one-flat-macro language, not for arbitrary dictionaries/SLPs or storage bytes. The dictionary is selected post hoc from this finite corpus, not a held-out search-speed experiment.
All eight parsed outputs were expanded to their exact original ID sequences and replayed, preserving raw lengths, endpoints, peak and work.

## Resource-summary application and concrete warning
For a path fragment record (source,target,raw_length,peak,tail_work), tail_work=sum of L over states after its source. Compatible composition adds raw lengths/tail_work and takes max peak; full work is L(initial)+tail_work. The program checks every4/4 split and compressed token segmentation against direct trace.
In the minimizing shared macro's actual context u=((-1,-2),(2,1,1)), M=(2,5,5) ends at t, raw length3, peak5, tail_work9 (full local work14).
A and B have identical full endpoints, raw length8 and peak9, but work51 vs49. Therefore endpoint+raw-length+peak alone does not recover cumulative work, even on shortest paths of a real frozen AC instance. This does not refute using a correctly resource-enriched state or proved dominance.
