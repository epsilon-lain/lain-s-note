"""Native umbrella audit: exact path summaries; no AC official-validator claim."""
import itertools
import json
from pathlib import Path

def signature(path, resource):
    return (len(path)-1, max(resource[v] for v in path),
            sum(resource[v] for v in path[1:]))

def compose(a, b):
    return (a[0]+b[0], max(a[1], b[1]), a[2]+b[2])

def enumerate_paths(edges, start, max_length):
    out = [(start,)]
    layer = [(start,)]
    for _ in range(max_length):
        layer = [p+(v,) for p in layer for v in range(3)
                 if (p[-1], v) in edges]
        out.extend(layer)
    return out

def reachable_signatures(edges, resource, start, macros):
    layers = [set() for _ in range(5)]
    layers[0].add((start, 0, resource[start], 0))
    outgoing = {v: [] for v in range(3)}
    for u, v in edges:
        outgoing[u].append((v, signature((u,v), resource)))
    for p in macros:
        outgoing[p[0]].append((p[-1], signature(p, resource)))
    for length in range(5):
        for vertex, _, peak, work in list(layers[length]):
            for dest, summary in outgoing[vertex]:
                nxt = compose((length,peak,work), summary)
                if nxt[0] <= 4:
                    layers[nxt[0]].add((dest,)+nxt)
    return layers

def main():
    possible = [(u,v) for u in range(3) for v in range(3) if u != v]
    counts = dict(graphs=0, graph_resource_pairs=0, path_checks=0,
                  cut_checks=0, layered_set_equalities=0,
                  raw_signature_records=0, macro_instances=0)
    for mask in range(1 << len(possible)):
        edges = {e for i,e in enumerate(possible) if mask & (1 << i)}
        paths_by_start = {s:enumerate_paths(edges,s,4) for s in range(3)}
        macros = [p for paths in paths_by_start.values() for p in paths if len(p)==3]
        counts["graphs"] += 1
        for resource in itertools.product(range(3), repeat=3):
            counts["graph_resource_pairs"] += 1
            counts["macro_instances"] += len(macros)
            for start in range(3):
                for path in paths_by_start[start]:
                    direct = signature(path,resource)
                    counts["path_checks"] += 1
                    for i in range(len(path)):
                        for j in range(i,len(path)):
                            a=signature(path[:i+1],resource)
                            b=signature(path[i:j+1],resource)
                            c=signature(path[j:],resource)
                            assert compose(compose(a,b),c)==direct
                            assert compose(a,compose(b,c))==direct
                            counts["cut_checks"] += 1
                raw=reachable_signatures(edges,resource,start,[])
                augmented=reachable_signatures(edges,resource,start,macros)
                for a,b in zip(raw,augmented):
                    assert a==b
                    counts["layered_set_equalities"] += 1
                    counts["raw_signature_records"] += len(a)
    resource=(0,2,0)
    counterexample=dict(edges=[[0,1],[1,2]], resource=resource,
        path=[0,1,2], exact_signature=signature((0,1,2),resource),
        endpoint_only_peak=0, peak_budget=0,
        raw_budget_feasible=False, endpoint_only_macro_feasible=True)
    assert counterexample["exact_signature"]==(2,2,2)
    result=dict(format="native-root-resource-signature-audit/1", status="passed",
        counts=counts, scope=dict(vertices=3,self_loops=False,
            resource_values=[0,1,2], max_expanded_length=4,
            macros="all legal length-two walks; all primitive edges retained"),
        convention="work sums resource of post-step vertices; peak includes initial vertex",
        counterexample=counterexample,
        caveat="Finite check supports but does not prove the general theorem; not an official AC validation.")
    Path("/output/root-audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))

if __name__ == "__main__":
    main()
