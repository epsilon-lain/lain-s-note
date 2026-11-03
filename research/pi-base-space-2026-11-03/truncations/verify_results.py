#!/usr/bin/env python3
"""Independent radix-prefix checks and exact recomputation of stored candidates."""
import json
from pathlib import Path
import math
from scan_bases import g, floor_exact, log10

HERE=Path(__file__).resolve().parent


def decode(pair):
    return g.mpq(g.mpz(pair[0],16),g.mpz(pair[1],16))


def main():
    data=json.loads((HERE/'base_scan_results.json').read_text())
    length=data['parameters']['digits']+data['parameters']['lookahead']
    bits=data['binary_precision_pi_sqrt2']+256
    targets={name:(decode(z['lower']),decode(z['upper'])) for name,z in data['target_exact_bounds_hex'].items()}
    with g.context(precision=bits,round=g.RoundDown):
        pl=g.mpq(g.const_pi())
        sl=g.mpq(g.sqrt(g.mpfr(2)))
    with g.context(precision=bits,round=g.RoundUp):
        ph=g.mpq(g.const_pi())
        sh=g.mpq(g.sqrt(g.mpfr(2)))
    fresh={'pi':(pl,ph),'sqrt2':(sl,sh)}
    prefix_checks=0
    for name,(reflo,refhi) in fresh.items():
        lo,hi=targets[name]
        assert lo<=reflo<refhi<=hi
        for b in range(2,data['parameters']['max_base']+1):
            scale=g.mpz(b)**length
            assert floor_exact(lo*scale)==floor_exact(hi*scale)==floor_exact(reflo*scale)==floor_exact(refhi*scale)
            prefix_checks+=1
    candidate_checks=0
    with g.context(precision=100):
        for name,rows in data['targets'].items():
            lo,hi=targets[name]
            for z in rows:
                b=z['base']
                for key in ['best_reduced_C','best_reduced_C_q_gt_1','best_reduced_C_m_ge_32']:
                    r=z[key]
                    p,q=g.mpz(r['p_hex'],16),g.mpz(r['q_hex'],16)
                    assert g.gcd(p,q)==1
                    assert floor_exact(lo*g.mpz(b)**r['m'])==p*(g.mpz(b)**r['m']//q)
                    cl,ch=q*q*(lo-g.mpq(p,q)),q*q*(hi-g.mpq(p,q))
                    assert cl>0
                    assert abs(log10(cl)-r['log10_C_lower'])<1e-10
                    assert abs(log10(ch)-r['log10_C_upper'])<1e-10
                    candidate_checks+=1
    partial=sum((g.mpq(1,g.mpz(10)**math.factorial(k)) for k in range(1,7)),g.mpq(0))
    compatible=[]
    for z in data['targets']['liouville_decimal']:
        r=z['best_reduced_C']
        if r['log10_C_upper']<-3500:
            assert g.mpq(g.mpz(r['p_hex'],16),g.mpz(r['q_hex'],16))==partial
            assert g.mpz(r['q_hex'],16)==g.mpz(10)**720
            compatible.append(z['base'])
    out={'independent_full_prefix_checks':prefix_checks,
         'fresh_methods':{'pi':'MPFR const_pi directed endpoints, distinct from scan Machin formula',
                          'sqrt2':'MPFR sqrt directed endpoints, distinct from scan integer-square-root bounds'},
         'saved_candidate_checks':candidate_checks,
         'liouville_S6_compatible_bases':compatible,
         'all_checks_passed':True}
    (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
