"""Certify an actual-period nonzero grouped eta^2 zeta_1^2 coefficient.

Computations use finite fields; the coefficient is a polynomial in the
formal period, so a nonzero residue at 1 certifies a nonzero polynomial.
"""
import json
from pathlib import Path
from verify_fiber_actual_period import P, rows, columns, entry, det_mod


def solve_matrix(left,right):
    n=len(left);a=[l[:]+r[:] for l,r in zip(left,right)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if a[i][j])
        if pivot!=j:a[j],a[pivot]=a[pivot],a[j]
        inverse=pow(a[j][j],-1,P)
        a[j]=[(x*inverse)%P for x in a[j]]
        for i in range(n):
            if i==j:continue
            factor=a[i][j]
            if factor:
                a[i]=[(x-factor*y)%P for x,y in zip(a[i],a[j])]
    return [line[n:] for line in a]


def run():
    cert=json.loads(Path(__file__).with_name('fiber_actual_period_verification.json').read_text())
    cols=[columns[i] for i in cert['selected_column_indices']]
    b=[[entry(row,col,(1,1)) for col in cols] for row in rows]
    db=det_mod(b);assert db==21281
    e=[];active=[]
    for index,(j,s,beta) in enumerate(rows):
        gap=20-s-4*beta[0]-10*beta[1]
        if s>=1 and gap<=3:
            target=(j,s-1,(beta[0]+1,beta[1]))
            e.append([(-s*(beta[0]+1)*entry(target,col,(1,1)))%P for col in cols])
            active.append(index)
        else:e.append([0]*len(cols))
    v=solve_matrix(b,e)
    # Independent residual check of the solved relation B V = E.
    assert all(sum(b[i][k]*v[k][j] for k in range(len(rows)))%P==e[i][j]
               for i in range(len(rows)) for j in range(len(rows)))
    tr=sum(v[i][i] for i in range(len(rows)))%P
    tr2=sum(v[i][j]*v[j][i] for i in range(len(rows)) for j in range(len(rows)))%P
    c2=((tr*tr-tr2)*pow(2,-1,P))%P
    coefficient=db*c2%P
    assert coefficient
    # Check rank-bound alternative: coefficient2 from low rank row support.
    # E=I_active E_active; Sylvester reduces the determinant to an active block.
    identity_columns=[[int(i==j) for j in active] for i in range(len(rows))]
    b_inverse_columns=solve_matrix(b,identity_columns)
    small=[[sum(e[i][k]*b_inverse_columns[k][l] for k in range(len(rows)))%P
            for l in range(len(active))] for i in active]
    small_coefficient=sum(small[i][i]*small[j][j]-small[i][j]*small[j][i]
                          for i in range(len(active)) for j in range(i+1,len(active)))%P
    assert small_coefficient==c2
    result={'status':'passed','prime':P,'formal_period':1,'normalization':'each source row multiplied by s!',
            'source_period_determinant_residue':db,'active_first_order_boundary_rows':len(active),
            'trace_inverse_boundary':tr,'trace_inverse_boundary_squared':tr2,
            'degree_two_relative_coefficient':c2,'eta2_zeta1_2_coefficient_residue':coefficient,
            'independent_low_rank_coefficient':small_coefficient,
            'scope':'The fully aggregated eta^2*zeta_1^2 coefficient is a nonzero polynomial in the period and hence nonzero at 2*pi*i by known transcendence. No lower bound for its value and no claim about cancellation at specified nonzero zeta_2 or in eta=1 full determinant.'}
    Path(__file__).with_name('aggregate_degree_two_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':run()
