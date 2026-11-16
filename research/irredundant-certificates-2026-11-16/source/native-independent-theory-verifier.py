"""Root independent literal-deletion and direct opposite-graph certificate audit."""
from collections import Counter
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
import json

def parts(total,lower=1):
    if total==0:
        yield ()
    for a in range(lower,total+1):
        for rest in parts(total-a,a):
            yield (a,)+rest

def bits(A):
    b=1
    for a in A:
        b|=b<<a
    return b

def balanced(A,m):
    return bool((bits(A)>>m)&1)

def minimal_by_all_deletions(A,m):
    if balanced(A,m):
        return False
    for i,a in enumerate(A):
        other=A[:i]+A[i+1:]
        for t in range(1,a):
            if not balanced(other+(t,a-t),m):
                return False
    return True

def labelled_count(A):
    numerator=sum(A)*factorial(len(A)-1)
    denominator=1
    for mult in Counter(A).values():
        denominator*=factorial(mult)
    assert numerator%denominator==0
    return numerator//denominator

def component_sizes(mask,n):
    parent=list(range(n))
    def find(x):
        while parent[x]!=x:
            x=parent[x]
        return x
    for i in range(n):
        if mask&(1<<i):
            a,b=find(i),find((i+1)%n)
            parent[a]=b
    return tuple(sorted(Counter(find(i) for i in range(n)).values()))

def cycle_edges(group):
    first=min(group)
    for rest in permutations(sorted(set(group)-{first})):
        if rest[0]<rest[-1]:
            order=(first,)+rest
            yield frozenset(tuple(sorted((order[i],order[(i+1)%len(order)])))
                            for i in range(len(order)))

def direct_check(m):
    n=2*m
    target=[tuple(sorted((i,(i+1)%n))) for i in range(n)]
    support=[]
    graphs=set()
    for tail in combinations(range(1,n),m-1):
        first=(0,)+tail
        second=tuple(sorted(set(range(n))-set(first)))
        for c in cycle_edges(first):
            for d in cycle_edges(second):
                graph=c|d
                assert graph not in graphs
                graphs.add(graph)
                support.append(sum(1<<i for i,e in enumerate(target) if e in graph))
    expected={3:10,4:315}[m]
    assert len(graphs)==expected
    compatible=[0]*(1<<n)
    for mask in support:
        sub=mask
        while True:
            compatible[sub]+=1
            if sub==0:break
            sub=(sub-1)&mask
    partition_counts=Counter()
    minimal_costs=Counter()
    for mask in range(1<<n):
        direct_cert=compatible[mask]==0
        direct_min=direct_cert and all(compatible[mask^(1<<i)]>0
                    for i in range(n) if mask&(1<<i))
        if mask==(1<<n)-1:
            assert direct_cert and not direct_min
            continue
        A=component_sizes(mask,n)
        partition_counts[A]+=1
        assert direct_cert==(not balanced(A,m))
        assert direct_min==minimal_by_all_deletions(A,m)
        if direct_min:
            minimal_costs[mask.bit_count()]+=1
    for A in parts(n):
        assert partition_counts[A]==labelled_count(A)
    return dict(m=m,opposite_graphs=len(graphs),target_edge_masks=1<<n,
                minimal_labelled_by_cost=sorted(minimal_costs.items()),
                proper_partition_counts_checked=len(partition_counts))

def main():
    spectra=[]
    for m in range(3,21):
        pc=0
        proper_count=0
        minimal=[]
        weighted=Counter()
        pattern=Counter()
        for A in parts(2*m):
            pc+=1
            count=labelled_count(A)
            proper_count+=count
            if minimal_by_all_deletions(A,m):
                cost=2*m-len(A)
                minimal.append(A)
                pattern[cost]+=1
                weighted[cost]+=count
        assert proper_count==(1<<(2*m))-1
        spectra.append(dict(m=m,partitions_checked=pc,spectrum=sorted(pattern),
            minimal_patterns=len(minimal),patterns_by_cost=sorted(pattern.items()),
            labelled_by_cost=sorted(weighted.items()),
            max_cost_patterns=[A for A in minimal if 2*m-len(A)==max(pattern)]))
    index={r["m"]:r for r in spectra}
    assert index[6]["spectrum"]==[6,8]
    assert index[12]["spectrum"]==[12,17,18]
    assert index[20]["spectrum"]==[20,26,29,30,31,32,33,34]
    assert (2,6,6,7,9,10) in index[20]["max_cost_patterns"]
    result=dict(format="root-independent-certificate-audit/1",status="passed",
        method="literal every-edge split, no interval shortcut; all partitions through2m=40; separate direct labelled opposite graphs for m3/4",
        reviewed_math=dict(id="M-05ac7d1d25ae1353",revision=1),
        spectra=spectra,direct_graph_checks=[direct_check(3),direct_check(4)],
        caveat="Finite reproduction is not the general proof; counts assume fixed labelled target cycle and positive-edge observations.")
    Path("/output/root-certificate-verify.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(status=result["status"],
        spectra=[dict(m=r["m"],spectrum=r["spectrum"]) for r in spectra],
        direct_graph_checks=result["direct_graph_checks"]),sort_keys=True))

if __name__=="__main__":
    main()
