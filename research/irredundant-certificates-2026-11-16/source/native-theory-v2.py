"""Native extension: all inclusion-minimal positive-edge type certificates.
External baseline: external/positive-edge-general-m@1 (fixed README seed-read).
Standard library only; exhaustive partitions and an independent-in-method
candidate-graph check are NOT independent authorship/review.
"""
import hashlib
import itertools
import json
import math
import time
from collections import Counter
from pathlib import Path

EXTERNAL = {
    "id": "external/positive-edge-general-m",
    "revision": 1,
    "origin": "external_chatgpt_work",
    "fixed_readme_sha256": "e08afb8d4beb43140d53021fcf3cf790c5a2cb06571e35412e9f6e9110369fd5",
    "used": ["U_m and positive-edge cost", "forest subset-sum completion criterion", "minimum cardinality m"],
}

def partitions(n, lo=1):
    if n == 0:
        yield ()
    else:
        for a in range(lo, n + 1):
            for tail in partitions(n - a, a):
                yield (a,) + tail

def subset_bits(parts):
    bits = 1
    for a in parts:
        bits |= bits << a
    return bits

def certifies(parts, m):
    return not ((subset_bits(parts) >> m) & 1)

def irredundant_split(parts, m):
    if not certifies(parts, m):
        return False
    for i, a in enumerate(parts):
        for t in range(1, a):
            q = parts[:i] + parts[i+1:] + (t, a-t)
            if certifies(q, m):
                return False
    return True

def irredundant_interval(parts, m):
    if not certifies(parts, m):
        return False
    for a in set(parts):
        if a == 1:
            continue
        i = parts.index(a)
        other = subset_bits(parts[:i] + parts[i+1:])
        lo = m - a + 1
        hi = m - 1
        if lo < 0:
            return False
        need = ((1 << (hi-lo+1)) - 1) << lo
        if other & need != need:
            return False
    return True

def multiplicity_on_cycle(parts):
    n = sum(parts)
    den = math.prod(math.factorial(v) for v in Counter(parts).values())
    num = n * math.factorial(len(parts)-1)
    assert num % den == 0
    return num // den

def canonical_edge_subset(parts):
    cursor = 0
    edges = []
    for a in parts:
        edges.extend([cursor+i, cursor+i+1] for i in range(a-1))
        cursor += a
    assert cursor == sum(parts)
    return edges

def enumerate_partitions(m):
    seen = 0
    certs = 0
    patterns = []
    by_cost = Counter()
    labeled_by_cost = Counter()
    for p in partitions(2*m):
        seen += 1
        if not certifies(p, m):
            continue
        certs += 1
        ok = irredundant_interval(p, m)
        # Different implementation, both authored by this same native node.
        if m <= 12:
            assert ok == irredundant_split(p, m), (m,p)
        if ok:
            cost = 2*m-len(p)
            count = multiplicity_on_cycle(p)
            by_cost[cost] += 1
            labeled_by_cost[cost] += count
            patterns.append({"parts": p, "cost": cost, "cycle_edge_subset_count": count})
    assert min(by_cost) == m
    maximum = max(by_cost)
    maxima = [p["parts"] for p in patterns if p["cost"] == maximum]
    return {
        "m":m, "integer_partitions_examined":seen,
        "no_half_subset_partitions":certs,
        "inclusion_minimal_patterns":len(patterns),
        "cost_spectrum":dict(sorted(by_cost.items())),
        "labeled_cycle_subset_counts":dict(sorted(labeled_by_cost.items())),
        "minimum_cost":min(by_cost), "maximum_cost":maximum,
        "maximum_patterns":maxima,
        "one_maximum_edge_witness":canonical_edge_subset(maxima[0]),
        "all_irredundant_patterns":patterns,
    }

def cycle_edge_set(order):
    return frozenset(tuple(sorted((order[i-1],order[i]))) for i in range(len(order)))

def cycles_on(vertices):
    vertices = tuple(sorted(vertices))
    first = vertices[0]
    for tail in itertools.permutations(vertices[1:]):
        if tail[0] < tail[-1]:
            yield cycle_edge_set((first,)+tail)

def component_sizes(n, selected):
    parent = list(range(n))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for a,b in selected:
        parent[find(a)] = find(b)
    return tuple(sorted(Counter(find(a) for a in range(n)).values()))

def candidate_check(m, summary):
    n = 2*m
    vertices = tuple(range(n))
    target = [tuple(sorted((i,(i+1)%n))) for i in range(n)]
    counts = {"single":[0]*(1<<n),"double":[0]*(1<<n)}
    totals = Counter()
    def record(kind, edges):
        totals[kind] += 1
        mask = sum(1<<i for i,e in enumerate(target) if e in edges)
        counts[kind][mask] += 1
    for c in cycles_on(vertices):
        record("single",c)
    for rest in itertools.combinations(range(1,n),m-1):
        a = (0,)+rest
        aa = set(a)
        b = tuple(v for v in vertices if v not in aa)
        for ca in cycles_on(a):
            for cb in cycles_on(b):
                record("double",ca|cb)
    # Superset zeta: exact candidate counts for every observation edge subset.
    for freq in counts.values():
        for i in range(n):
            for mask in range(1<<n):
                if not (mask & (1<<i)):
                    freq[mask] += freq[mask | (1<<i)]
    full = (1<<n)-1
    certificates = [counts["double"][s] == 0 for s in range(1<<n)]
    irredundant = []
    by_cost = Counter()
    for s in range(1<<n):
        assert counts["single"][s] > 0
        if s != full:
            chosen = [e for i,e in enumerate(target) if s & (1<<i)]
            p = component_sizes(n,chosen)
            assert certificates[s] == certifies(p,m), (m,s,p)
            minimal_formula = irredundant_interval(p,m)
        else:
            minimal_formula = False
        minimal_graph = certificates[s] and all(
            not certificates[s ^ (1<<i)] for i in range(n) if s & (1<<i)
        )
        assert minimal_formula == minimal_graph, (m,s)
        if minimal_graph:
            irredundant.append(s)
            by_cost[s.bit_count()] += 1
    expected = {int(k):v for k,v in summary["labeled_cycle_subset_counts"].items()}
    assert dict(by_cost) == expected
    return {
        "m":m,"all_labeled_candidates":dict(totals),
        "all_target_edge_subsets":1<<n,
        "candidate_count_sha256":hashlib.sha256(json.dumps(counts,sort_keys=True).encode()).hexdigest(),
        "irredundant_counts_by_cost":dict(sorted(by_cost.items())),
        "irredundant_masks":irredundant,
        "component_subset_sum_and_minimality_agree_for_all_subsets":True,
    }

def family_check(t):
    m = 4*t
    p = (2,)*(2*t-1)+(2*t+1,2*t+1)
    assert sum(p) == 2*m
    assert irredundant_interval(p,m)
    assert irredundant_split(p,m)
    cost = 2*m-len(p)
    assert cost == 6*t-1
    return {"t":t,"m":m,"parts":p,"cost":cost,"optimum_external_baseline":m,"excess":2*t-1}


def interval_present(bits, lo, hi):
    assert lo >= 0 and hi >= lo
    need = ((1 << (hi-lo+1))-1) << lo
    return bits & need == need

def near_double_family(r, s):
    assert r >= 1 and 0 <= s <= r+1
    p = (1,)*(r-1)+(r+1,)
    half, width = r, r
    stages = []
    for j in range(1,s+1):
        b = (1<<j)*r+1
        assert b == 2*width+1 and b > half
        p = tuple(sorted(p+(b,b)))
        half += b
        width *= 2
        bits = subset_bits(p)
        assert not ((bits>>half)&1)
        assert interval_present(bits,half-width,half-1)
        assert interval_present(bits,half+1,half+width)
        assert irredundant_interval(p,half)
        if half <= 512:
            assert irredundant_split(p,half)
        stages.append({"j":j,"half":half,"near_half_interval_width":width,"new_part":b})
    assert half == ((1<<(s+1))-1)*r+s
    assert len(p) == r+2*s
    if s == r+1:
        assert half == r*(1<<(r+2))+1
    return {"r":r,"s":s,"m":half,"component_count":len(p),"parts":p,
            "cost":2*half-len(p),"optimal_cost_external":half,
            "ratio_as_exact_fraction":[2*half-len(p),half],"stages":stages}

def logarithmic_component_bound_check(rows):
    tested = 0
    for row in rows:
        for pattern in row["all_irredundant_patterns"]:
            p = pattern["parts"]
            k = len(p)
            assert all(a-1 <= (1<<(k-1)) for a in p)
            assert sum(p) <= k*((1<<(k-1))+1)
            tested += 1
    return {"all_irredundant_patterns_tested":tested,
            "bound":"2m <= k*(2^(k-1)+1)","assertions_passed":True}

def main():
    start = time.monotonic()
    rows = []
    for m in range(3,21):
        row = enumerate_partitions(m)
        rows.append(row)
        print(json.dumps({k:row[k] for k in ["m","integer_partitions_examined","inclusion_minimal_patterns","cost_spectrum","minimum_cost","maximum_cost","maximum_patterns"]}),flush=True)
    graph = [candidate_check(m,rows[m-3]) for m in (3,4)]
    family = [family_check(t) for t in range(1,51)]
    near_double = [near_double_family(r,s) for r in range(1,13) for s in range(1,r+2)]
    component_bound = logarithmic_component_bound_check(rows)
    result = {
        "format":"native-irredundant-cycle-certificates/2",
        "claim_scope":"integer partitions for 3<=m<=20; graph candidates only m=3,4; family checks only 1<=t<=50 and near-double 1<=r<=12, 1<=s<=r+1; general statements require written proofs",
        "external_dependencies":[EXTERNAL],
        "native_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "algorithm":"all nondecreasing integer partitions of 2m; bitset subset sums; every one-edge split criticality",
        "independent_authorship_review":False,
        "partition_spectra":rows,
        "candidate_graph_checks":graph,
        "infinite_family_finite_checks":family,
        "near_double_family_finite_checks":near_double,
        "component_bound_finite_checks":component_bound,
        "native_math_dependencies":[{"id":"M-05ac7d1d25ae1353","revision":1}],
        "elapsed_seconds":time.monotonic()-start,
        "assertions_passed":True,
    }
    Path("/output/certificate-results.json").write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"all_assertions_passed":True,"seconds":result["elapsed_seconds"],"graph_checks":graph,"near_double_final_case":near_double[-1]}),flush=True)

if __name__ == "__main__":
    main()
