"""Full integer gcd audit of indices from the small-prime modular scan."""
import json
from math import gcd
from pathlib import Path
from time import perf_counter


def main():
    folder = Path(__file__).parent
    data = json.loads((folder/'combined_surplus_checks.json').read_text())
    known = {n:r['p'] for r in data['prime_scans'] for n in r['hits']}
    targets = set(known)
    a0=a1=1
    b0=0
    b1=1
    fac=1
    records=[]
    combined=[]
    started=perf_counter()
    for n in range(1,max(targets)+1):
        if n>=2:
            k=n-1
            a0,a1=a1,(2*k+1)*a1+k*k*a0
            b0,b1=b1,(2*k+1)*b1+k*k*b0
        fac*=n
        if n in targets:
            g=gcd(a1,4*b1)
            e=g//gcd(g,fac)
            while e%2==0:
                e//=2
            p=known[n]
            z=e
            while z%p==0:
                z//=p
            rec={'n':n,'known_surplus_prime':p,'odd_excess_factor':str(e),
                 'remaining_after_known_prime_powers':str(z)}
            records.append(rec)
            if z>1:
                combined.append(rec)
    output={'status':'PASS','candidate_count':len(records),'max_index':max(targets),
            'combined_events':combined,'records':records,
            'seconds':round(perf_counter()-started,3),
            'scope':'Complete exact gcd at all indices with a surplus prime<=500 through n50000; does not exclude pairs where both primes exceed500.'}
    (folder/'full_candidate_surplus_checks.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='records'},indent=2))


if __name__=='__main__':
    main()
