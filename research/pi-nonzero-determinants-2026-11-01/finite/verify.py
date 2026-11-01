"""Exact certificates for actual translated minors and a Hall counterexample.

Uses only standard-library rational/modular arithmetic. Finite checks do not
certify the unknown weight-separation threshold or a uniform asymptotic claim.
"""
import json
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import experiment as ex

MODULUS=1_000_000_009
ROOT_I=430_477_711


def determinant_mod(matrix):
    a=[row[:] for row in matrix]
    n=len(a)
    assert all(len(row)==n for row in a)
    d=1
    for j in range(n):
        pivot=next((i for i in range(j,n) if a[i][j]),None)
        if pivot is None:return 0
        if pivot!=j:a[pivot],a[j]=a[j],a[pivot];d=-d
        value=a[j][j]
        d=d*value%MODULUS
        inv=pow(value,-1,MODULUS)
        for i in range(j+1,n):
            factor=a[i][j]*inv%MODULUS
            for k in range(j+1,n):
                a[i][k]=(a[i][k]-factor*a[j][k])%MODULUS
    return d%MODULUS


def main():
    ex.PRIME=MODULUS
    assert all(MODULUS%d for d in range(2,isqrt(MODULUS)+1))
    assert ROOT_I*ROOT_I%MODULUS==MODULUS-1
    source=Path(__file__).with_name('experiment.json')
    examples=json.loads(source.read_text())
    certificates=[]
    for example in examples:
        m,H=example['m'],example['H']
        weights=example['weights']
        cols=[(h,tuple(alpha)) for h,alpha in example['columns']]
        rows=[(j,s,tuple(beta)) for j,s,beta in example['rows']]
        assignment=[tuple(a) for a in example['assignment']]
        indices=example['minor_column_indices']
        base=[row[2] for row in rows]
        assert cols==ex.columns(m,H,weights,F(example['w0']))
        assert rows==ex.rows(m,H,weights,example['K'],F(example['theta']))
        assert len(indices)==len(rows) and len(set(indices))==len(indices)
        assert all(0<=i<len(cols) for i in indices)
        assert len(assignment)==len(rows)
        assert max(row[1] for row in rows)<MODULUS
        # The previous finite fixtures use theta=1/8 or 1/4, v0=1,
        # q=2 and p=6. Thus wi=ceil(ln2)=1 as in the manuscript.
        assert m==1 and weights==[1]
        assert example['K']*F(example['theta'])**m<1
        assert example['K']*F(example['w0'])*F(example['theta'])**m<1
        assert any(a!=r[2] for r,a in zip(rows,assignment))
        assert all(a==r[2] for r,a in zip(rows,assignment) if r[0]==0)
        assert all(all(x>=b for x,b in zip(a,r[2])) for r,a in zip(rows,assignment))
        assert all(sum(w*x for w,x in zip(weights,a))<=H for a in assignment)
        assert len({(j,s,a) for (j,s,beta),a in zip(rows,assignment)})==len(rows)
        testing=ex.vectors(cols,rows,assignment,W=2)
        translated_minor=[[row[c] for c in indices] for row in testing]
        translated_det=determinant_mod(translated_minor)
        assert translated_det
        # Original interpolation matrix at r=2ip/q=6i, before translating.
        original=ex.vectors(cols,rows,base,W=6*ROOT_I%MODULUS)
        original_minor=[[row[c] for c in indices] for row in original]
        original_det=determinant_mod(original_minor)
        assert original_det
        before=ex.cost_parts(rows,base,H,weights)
        after=ex.cost_parts(rows,assignment,H,weights)
        assert after[1]<before[1]
        # ln2/4 > 1/6. Because collision count decreases, this certifies
        # a strict cost decrease for the actual constant c=ln2/4.
        assert after[0]-before[0]+F(1,6)*(after[1]-before[1])<0
        bbar=F(sum(sum(w*b for w,b in zip(weights,row[2])) for row in rows),len(rows)*H)
        assert after[0]+F(1,4)*after[1] < 1-bbar
        certificates.append({
            'm':m,'H':H,'K':example['K'],'theta':example['theta'],
            'w0':example['w0'],'v0':1,'w1':1,'p':6,'q':2,
            'rows':len(rows),'columns':len(cols),'shifted_rows':example['shifted_rows'],
            'selected_column_indices':indices,
            'translated_minor_at_W2_mod_prime':translated_det,
            'same_original_minor_at_R6i_mod_prime':original_det,
            'base_cost_parts':[str(x) for x in before],
            'shifted_cost_parts':[str(x) for x in after],
            'finite_bbar':str(bbar),
            'cost_below_1_minus_bbar_certified':True,
            'strict_cost_improvement_certified':True,
            'distinct_test_rows':True,
        })

    # A complete column support matching exists, yet the actual test matrix
    # is singular for all W. Verify the seven polynomial monomials exactly.
    fcoeff=(11,27,-27,-11)
    gcoeff=(3,27,27,3)
    residuals=[]
    for degree in range(8):
        value=sum(f*j**degree for j,f in enumerate(fcoeff))
        if degree:
            value+=sum(g*degree*j**(degree-1) for j,g in enumerate(gcoeff))
        residuals.append(value)
    assert residuals[:7]==[0]*7 and residuals[7]==108
    cols=ex.columns(2,6,(1,1),F(1,4))
    rows=ex.rows(2,6,(1,1),K=4,theta=F(1,3))
    assert len(rows)==48 and len(cols)==252
    matching=[]
    counts={}
    for j,s,beta in rows:
        k=counts.get(beta,0);counts[beta]=k+1
        column=(k,beta)
        assert column in cols
        assert k>=s
        matching.append(cols.index(column))
    assert len(set(matching))==len(rows)
    ranks=[]
    for W in (1,2,3,7):
        matrix=ex.vectors(cols,rows,[r[2] for r in rows],W)
        rank,_=ex.rank_pivots(matrix)
        assert rank<48
        ranks.append({'W':W,'rank_mod_prime':rank})
    report={
        'status':'passed','modulus':MODULUS,'sqrt_minus_one':ROOT_I,
        'nonzero_minors':certificates,
        'hall_counterexample':{'m':2,'H':6,'K':4,'theta':'1/3','w0':'1/4',
                              'v0':1,'w1':1,'w2':1,'rows':48,'columns':252,
                              'full_column_support_matching':matching,
                              'monomial_relation_residuals_degrees_0_to_7':residuals,
                              'finite_specialization_ranks':ranks,
                              'universal_singularity_reason':'Explicit degree<=6 Hermite relation; not inferred from specialization ranks.'},
        'scope':'Finite structural certificates. Does not verify the unknown separated-weight threshold or prove a uniform H limit. Nonzero Q[W] minors transfer to W=2pi*i using the known transcendence of pi.',
    }
    Path(__file__).with_name('verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
