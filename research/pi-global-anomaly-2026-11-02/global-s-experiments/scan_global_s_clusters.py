"""Exact reflected-root scan for the first-layer upper-half anomaly budget.

The experiment does not establish an eventual bound. All comparisons against
n/6 and the observed maximum use integer certificates, not floating logs.
Floating logs are used only to rank/describe the saved diagnostic samples.
Python standard library only.
"""
from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from math import isqrt, log, prod
from pathlib import Path
import json
import time


def primes_upto(limit):
    sieve=bytearray(b'\x01')*(limit+1)
    sieve[:2]=b'\x00\x00'
    for p in range(2,isqrt(limit)+1):
        if sieve[p]:
            sieve[p*p:limit+1:p]=b'\x00'*((limit-p*p)//p+1)
    return [p for p in range(3,limit+1,2) if sieve[p]]


def reflected_scan(limit):
    """A_r mod p for only the lower member of each reflection pair."""
    factors=[[] for _ in range(limit+1)]
    prime_stats=[]
    steps=0
    for p in primes_upto(limit-1):
        end=min((p-1)//2,limit-p)
        previous,current=1,1
        low_roots=[]
        visible_roots=[]
        for r in range(2,end+1):
            previous,current=current,((2*r-1)*current+(r-1)**2*previous)%p
            steps+=1
            if not current:
                low_roots.append(r)
                factors[p+r].append(p)
                visible_roots.append(r)
                mirror=p-1-r
                if mirror!=r and p+mirror<=limit:
                    factors[p+mirror].append(p)
                    visible_roots.append(mirror)
        complete=end==(p-1)//2
        all_roots=sorted(set(low_roots+[p-1-r for r in low_roots])) if complete else None
        prime_stats.append({'p':p,'lower_roots':low_roots,
                            'visible_roots':sorted(visible_roots),
                            'full_digit_roots_by_reflection':all_roots,
                            'whole_digit_range_scanned':complete})
    return factors,prime_stats,steps


def direct_independent_scan(limit):
    factors=[[] for _ in range(limit+1)]
    for p in primes_upto(limit-1):
        previous,current=1,1
        for r in range(2,min(p-1,limit-p)+1):
            previous,current=current,((2*r-1)*current+(r-1)**2*previous)%p
            if not current:
                factors[p+r].append(p)
    return factors


def below_sixth(n,product):
    # product^6 < 2^n < e^n is enough, and this bit test avoids huge 2^n.
    value=product**6
    if value.bit_length()<=n:
        return 'product^6 < 2^n < e^n'
    # e>8/3 follows from 1+1+1/2+1/6 plus positive remaining terms.
    assert value*3**n<8**n
    return '3^n product^6 < 8^n < (3e)^n'


def below_reference(n,product):
    # S(n)/n < ln(17)/21 iff product^21 < 17^n.
    if n==21:
        assert product==17
        return 'equality at n=21'
    value=product**21
    if value.bit_length()<=4*n:
        return 'product^21 < 2^(4n) < 17^n'
    assert value<17**n
    return 'product^21 < 17^n directly'


@lru_cache(maxsize=None)
def top_half_lift(n,p):
    assert p<n<2*p
    r=n-p
    modulus=p*p
    previous,current=1,(2*p+1)%modulus
    for j in range(1,r):
        previous,current=current,((2*p+2*j+1)*current+(p+j)**2*previous)%modulus
    assert current%p==0
    return {'p':p,'shifted_factor_mod_p_squared':current,
            'gcd_valuation_at_p_by_established_lift_criterion':2 if not current else 1,
            'positive_excess_beyond_factorial':not current}


def sample(n,ps):
    return {'n':n,'anomalous_primes':ps,'prime_product':str(prod(ps)),
            'S_over_n_diagnostic':sum(log(p) for p in ps)/n,
            'p_over_n_diagnostic':[p/n for p in ps],
            'offsets_r':[n-p for p in ps],
            'reflected_partner_indices':[3*p-1-n for p in ps],
            'exact_top_half_lift_classes':[top_half_lift(n,p) for p in ps]}


def main():
    parser=ArgumentParser()
    parser.add_argument('--limit',type=int,default=100000)
    parser.add_argument('--direct-check-limit',type=int,default=10000)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    assert args.limit>=21
    started=time.monotonic()
    factors,stats,steps=reflected_scan(args.limit)
    check_limit=min(args.direct_check_limit,args.limit)
    direct=direct_independent_scan(check_limit)
    for n in range(check_limit+1):
        assert sorted(factors[n])==direct[n]
    classifications=Counter()
    reference_classes=Counter()
    clusters=Counter()
    residue_events=Counter()
    interval_events=Counter()
    midpoint_primes=[]
    internal_edges=[]
    nonempty=[]
    event_hash=sha256()
    for n,ps in enumerate(factors):
        if not n:
            continue
        ps.sort()
        assert len(ps)==len(set(ps))
        assert all(n//2<p<n for p in ps)
        clusters[len(ps)]+=1
        product=prod(ps)
        classifications[below_sixth(n,product)]+=1
        reference_classes[below_reference(n,product)]+=1
        if ps:
            nonempty.append((n,ps))
        for p in ps:
            event_hash.update(f'{n}:{p}\n'.encode())
            residue_events[p%8]+=1
            # Half-open bins with exact integer decisions, no float edges.
            interval=('(1/2,2/3]' if 3*p<=2*n else
                      '(2/3,3/4]' if 4*p<=3*n else
                      '(3/4,4/5]' if 5*p<=4*n else
                      '(4/5,9/10]' if 10*p<=9*n else '(9/10,1)')
            interval_events[interval]+=1
            partner=3*p-1-n
            if partner==n:
                midpoint_primes.append(p)
            if n<=partner<=args.limit:
                internal_edges.append((n,partner,p))
    ranked=sorted(nonempty,key=lambda x:sum(log(p) for p in x[1])/x[0],reverse=True)
    top_clusters=sorted(nonempty,key=lambda x:(len(x[1]),sum(log(p) for p in x[1])/x[0]),reverse=True)
    thresholds=(100,1000,10000,50000,100000)
    bands=[]
    lo=0
    for hi in thresholds:
        hi=min(hi,args.limit)
        if hi<=lo:
            continue
        rows=[x for x in nonempty if lo<x[0]<=hi]
        winner=max(rows,key=lambda x:sum(log(p) for p in x[1])/x[0]) if rows else None
        bands.append({'range':[lo+1,hi],'nonempty_indices':len(rows),
                      'events':sum(len(ps) for _,ps in rows),
                      'max_cluster_size':max((len(ps) for _,ps in rows),default=0),
                      'diagnostic_ratio_winner':sample(*winner) if winner else None})
        lo=hi
    complete=[x for x in stats if x['whole_digit_range_scanned']]
    root_heavy=sorted(complete,key=lambda x:(len(x['full_digit_roots_by_reflection']),x['p']),reverse=True)[:20]
    parent={n:n for n,_ in nonempty}
    def find(n):
        while parent[n]!=n:
            parent[n]=parent[parent[n]]
            n=parent[n]
        return n
    for n,m,p in internal_edges:
        parent[find(n)]=find(m)
    groups={}
    for n,ps in nonempty:
        groups.setdefault(find(n),[]).append(n)
    largest_groups=sorted(groups.values(),key=lambda ns:(len(ns),sum(len(factors[n]) for n in ns)),reverse=True)[:12]
    components=[{'vertices':ns,
                 'anomalies_at_vertices':[sample(n,factors[n]) for n in ns],
                 'internal_reflection_edges':[[n,m,p] for n,m,p in internal_edges if n in ns]}
                for ns in largest_groups]
    out={'status':'passed','limit':args.limit,
         'arithmetic':'Exact residues and integer inequality certificates; floats only in labeled diagnostics',
         'reflection':'A_r=0 mod p iff A_(p-1-r)=0 mod p for 0<=r<p',
         'modular_recurrence_steps_with_reflection':steps,
         'independent_direct_scan_through':check_limit,
         'nonempty_indices':len(nonempty),'anomaly_events':sum(len(ps) for _,ps in nonempty),
         'maximum_cluster_size':max(map(len,factors)),
         'all_S_n_less_than_n_over_6':True,
         'strict_maximum_ratio':{'index':21,'primes':[17],'exact_value':'ln(17)/21',
                                 'diagnostic_decimal':log(17)/21},
         'all_bound_certificate_counts':dict(classifications),
         'maximum_certificate_counts':dict(reference_classes),
         'event_stream_format':'Ascending n then ascending p, ASCII n:p followed by newline',
         'event_stream_sha256':event_hash.hexdigest(),
         'cluster_size_histogram':{str(k):v for k,v in sorted(clusters.items())},
         'event_prime_mod_8_histogram':{str(k):v for k,v in sorted(residue_events.items())},
         'event_p_over_n_intervals':dict(interval_events),
         'midpoint_anomalous_primes':midpoint_primes,
         'midpoint_prime_mod_32_histogram':dict(Counter(str(p%32) for p in midpoint_primes)),
         'reflection_graph':{'vertices_with_anomalies':len(parent),
                             'internal_edges_including_loops':len(internal_edges),
                             'components':len(groups),
                             'cycle_rank_after_removing_loops':len(internal_edges)-len(midpoint_primes)-len(parent)+len(groups),
                             'largest_component_vertex_count':max(map(len,groups.values())),
                             'largest_components':components},
         'range_summary':bands,
         'top_diagnostic_ratios':[sample(*x) for x in ranked[:30]],
         'top_cluster_samples':[sample(*x) for x in top_clusters[:30]],
         'root_heavy_complete_prime_ranges':root_heavy,
         'runtime_seconds_diagnostic':time.monotonic()-started,
         'scope':'Finite exact certification and structural diagnostics only. No eventual S(n)<=n/6 bound, subexponential anomaly product, or global pi approximation conclusion is proved.'}
    target=args.output or Path(__file__).with_name('global_s_cluster_checks.json')
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('top_diagnostic_ratios','top_cluster_samples','root_heavy_complete_prime_ranges','reflection_graph')},indent=2))


if __name__=='__main__':
    main()
