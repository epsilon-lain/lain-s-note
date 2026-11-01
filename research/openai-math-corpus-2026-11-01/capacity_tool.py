"""Exact signature capacities and finite jet-saturation certificates.

Requires SymPy. Supports the pure-period functions used in rank_bounds.md.
This gives necessary capacities for nonzero expansion terms, not a sufficient
condition and not an asymptotic Diophantine approximation theorem.
"""
from collections import defaultdict
from math import factorial
import json
import sympy as sp

Z, t, z = sp.symbols('Z t z')


def signature_capacity(columns, a):
    """Columns are (h, alpha_tuple); a is a displacement multi-index.

    r is the number of distinct (h, total-degree-minus-|a|) signatures.
    D is the dimension of their downward degree closure; D jets suffice
    to distinguish functions in that closure at a single center.
    """
    if any(len(alpha) != len(a) for _, alpha in columns):
        raise ValueError('dimension mismatch')
    if any(x < 0 for x in a) or any(h < 0 or any(x < 0 for x in alpha) for h, alpha in columns):
        raise ValueError('nonnegative integer indices required')
    signatures = sorted({(h, sum(alpha)-sum(a)) for h, alpha in columns
                         if all(x <= y for x, y in zip(a, alpha))})
    degrees = defaultdict(lambda: -1)
    for h, d in signatures:
        degrees[h] = max(degrees[h], d)
    return {'capacity': len(signatures),
            'jet_cutoff': sum(d+1 for d in degrees.values()),
            'signatures': signatures,
            'degree_closure': [(h, k) for h in sorted(degrees) for k in range(degrees[h]+1)]}


def coefficient_matrix(signatures, depth, use_log=True):
    """Exact first `depth` Taylor coefficients, one function per column."""
    matrix = []
    variable = t if use_log else z
    for ell in range(depth):
        values = []
        for h, d in signatures:
            if use_log:
                f = (1+t)**h * (Z+sp.log(1+t))**d
            else:
                f = sp.exp(h*z) * (Z+z)**d
            values.append(sp.expand(sp.diff(f, variable, ell).subs(variable, 0)/factorial(ell)))
        matrix.append(values)
    return sp.Matrix(matrix)


def main():
    cases = [
        ([(0,(0,0)),(0,(0,1)),(0,(0,2)),(0,(1,0))],(0,0)),
        ([(0,(0,0)),(0,(0,1)),(0,(0,2)),(0,(1,0))],(1,0)),
        ([(0,(0,0)),(0,(0,1)),(0,(0,2)),(0,(1,0))],(0,1)),
        ([(0,(0,0)),(0,(2,0)),(1,(0,0)),(1,(0,2))],(0,0)),
        ([(0,(1,1)),(2,(2,1)),(2,(1,2)),(0,(2,0))],(1,0)),
        ([(0,(1,0,1)),(1,(0,2,0)),(1,(1,0,1))],(0,0,0)),
    ]
    results = []
    for columns, a in cases:
        info = signature_capacity(columns,a)
        r,D = info['capacity'],info['jet_cutoff']
        coefficients = coefficient_matrix(info['signatures'], D)
        assert coefficients.rank() == r
        direct = coefficient_matrix(info['signatures'], D, use_log=False)
        assert direct.rank() == r
        # The complete downward-closed space has an invertible D-jet map.
        closure_matrix = coefficient_matrix(info['degree_closure'], D)
        det = sp.factor(closure_matrix.det())
        groups = defaultdict(int)
        for h,k in info['degree_closure']:
            groups[h] += 1
        predicted = sp.Integer(1)
        for h,size in groups.items():
            for k in range(size):
                predicted *= factorial(k)
        hs = sorted(groups)
        for i,h in enumerate(hs):
            for g in hs[i+1:]:
                predicted *= (g-h)**(groups[h]*groups[g])
        for ell in range(D):
            predicted /= factorial(ell)
        assert det != 0 and sp.simplify(det-predicted)==0
        results.append(dict(columns=columns,a=a,capacity=r,jet_cutoff=D,
                            signature_rank=coefficients.rank(),
                            closure_jet_determinant=str(det)))
    assert results[0]['capacity']==3
    no_shift_rows=4
    assert no_shift_rows>results[0]['capacity']
    print(json.dumps(dict(status='passed',sympy_version=sp.__version__,
                          exact_capacity_cases=len(cases),
                          zero_shift_example={'rows':4,'capacity':3,'excluded':True},
                          cases=results,
                          scope='Exact necessary capacities and finite jet saturation; '
                                'not sufficient independence, not an asymptotic lower displacement.'),indent=2))


if __name__=='__main__':
    main()
