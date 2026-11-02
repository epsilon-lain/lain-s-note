"""Recompute old and improved rational norm bounds on existing candidates.

This does not run a candidate search and does not use a value of pi.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb, gcd
import hashlib
import json


def main():
    root = Path(__file__).resolve().parent
    source = root.parent/'pi-constructive-redesign-2026-11-02/integer-kernel/new_family_results.json'
    raw = source.read_bytes()
    git_blob = hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    assert git_blob == '3cd2fd50b288c470977f776c9c34e6b85c1dccde'
    data = json.loads(raw)
    rows = []
    for entry in data['per_m']:
        m = entry['m']
        best = entry['best_observed']
        alpha = int(best['alpha'])
        rr = [int(x) for x in best['R_coefficients']]
        pp, q = int(best['p']), int(best['q'])
        h = sum(abs(int(x)) for x in best['H_coefficients'])
        assert F(-sum(rr), 2**(4*m-2)*alpha) == F(pp,q)
        assert gcd(pp,q) == 1 and q > 0 and alpha != 0
        old = F(4*q*q*h,abs(alpha))*F(2,27)**m
        new = F(22*q*q*h,7*abs(alpha))*F(comb(2*m,m),64**m)
        assert new < old
        rows.append(dict(m=m,p=pp,q=q,
                         old_quality_upper=str(old),
                         sharper_quality_upper=str(new),
                         improvement_factor=str(old/new),
                         displayed_new_upper=float(new),
                         displayed_improvement=float(old/new)))
    result = dict(status='passed',new_candidate_searches=0,
                  input_commit='db1290c827479ff85beaeca8606e8258d72f0c69',
                  source_git_blob=git_blob,samples=rows,
                  scope='All upper bounds are exact rational inequalities; displays are descriptive. No asymptotic height bound or quality-zero conclusion follows.')
    (root/'existing_candidate_bound_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
