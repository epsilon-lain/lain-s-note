"""Exact finite polynomial certificates, using only the standard library."""
from fractions import Fraction as F
from math import comb, factorial, prod, gcd, lcm
from itertools import product
from pathlib import Path
import json


def fall(a,k): return prod(range(a-k+1,a+1))


def determinant(matrix):
    a=[list(map(F,row)) for row in matrix]
    old,sign=F(1),1
    for k in range(len(a)-1):
        p=next((i for i in range(k,len(a)) if a[i][k]),None)
        if p is None:return F(0)
        if p!=k:a[p],a[k]=a[k],a[p];sign=-sign
        diag=a[k][k]
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                a[i][j]=(a[i][j]*diag-a[i][k]*a[k][j])/old
            a[i][k]=F(0)
        old=diag
    return sign*a[-1][-1]


def scaled_entry(row,col,z):
    j,s,b=row;h,a=col
    if any(x<y for x,y in zip(a,b)):return 0
    degree=sum(a)-sum(b)
    return prod(fall(x,y) for x,y in zip(a,b))*sum(
        comb(s,k)*h**(s-k)*fall(degree,k)*(j*z)**(degree-k)
        for k in range(min(s,degree)+1))


def source(H,v):
    types=[b for b in product(*(range(H//x+1) for x in v))
           if sum(x*y for x,y in zip(v,b))<H]
    # Explicit common order used by the determinant statements.
    order={(0,0):0,(1,0):1,(0,1):2}
    types.sort(key=order.__getitem__)
    return [(j,s,b) for j in (0,1) for b in types
            for s in range(H) if s+sum(x*y for x,y in zip(v,b))<H]


def evaluate(poly,z):return sum(c*z**k for k,c in poly.items())


def make_case(name,H,v,theta,w0,alphas,expected):
    rows=source(H,v)
    cols=[(h,a) for h in (0,1) for a in alphas]
    weights=[theta*x for x in v]
    assert len(rows)==len(cols)
    assert 2*theta**2<1 and 2*w0*theta**2<1
    slacks=[H-w0*h-sum(x*y for x,y in zip(weights,a)) for h,a in cols]
    assert min(slacks)>0
    D=sum(sum(a) for h,a in cols)-sum(sum(b) for j,s,b in rows)
    row_den=prod(factorial(s)*prod(factorial(x) for x in b) for j,s,b in rows)
    checks=0
    for z in list(range(D+1))+[F(-1,2),F(5,3)]:
        actual=determinant([[scaled_entry(r,c,z) for c in cols] for r in rows])
        assert actual==row_den*evaluate(expected,z)
        checks+=1
    d=lcm(*(c.denominator for c in expected.values()))
    integer={k:int(d*c) for k,c in expected.items()}
    return {'name':name,'H':H,'v0':1,'v':list(v),'theta':str(theta),
            'w0':str(w0),'w':[str(x) for x in weights],
            'rows':[[j,s,list(b)] for j,s,b in rows],
            'raw_columns':[[h,list(a)] for h,a in cols],
            'minimum_column_slack':str(min(slacks)),
            'proved_degree_bound':D,'integer_identity_sample_points':D+1,
            'additional_rational_checks':2,'identity_checks':checks,
            'row_factorial_multiplier':str(row_den),
            'polynomial':{str(k):str(c) for k,c in expected.items()},
            'least_integer_coefficient_denominator':d,
            'integer_coefficient_content':gcd(*integer.values())}


def f(x):return 13*x**3+720*x*x+6480*x-63936


def main():
    alphas=[(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(1,3),(1,4),(1,6)]
    mainpoly={27:F(13,864),25:F(720,864),23:F(6480,864),21:F(-63936,864)}
    maincase=make_case('new_residual',6,(4,5),F(1,10),F(1),alphas,mainpoly)
    assert maincase['minimum_column_slack']=='8/5'
    assert maincase['least_integer_coefficient_denominator']==864
    assert maincase['integer_coefficient_content']==1
    comp=alphas[:-1]+[(1,5)]
    # (Z²+24)(Z⁴+36 Z²−144)=Z⁶+60 Z⁴+720 Z²−3456.
    comppoly={25:F(1,864),23:F(60,864),21:F(720,864),19:F(-3456,864)}
    comparator=make_case('consecutive_exponent_comparator',6,(4,5),F(1,10),F(1),comp,comppoly)
    controlalphas=[(a,0) for a in range(5)]+[(a,1) for a in range(2)]
    control=make_case('positive_surplus_control',4,(2,3),F(1,3),F(1,6),controlalphas,{13:F(-2)})
    assert control['minimum_column_slack']=='7/6'
    assert sum(a[1][0] for a in maincase['raw_columns'])-sum(r[2][0] for r in maincase['rows'])==8
    assert sum(a[1][1] for a in control['raw_columns'])-sum(r[2][1] for r in control['rows'])==2
    residues=[int(f(x))%11 for x in range(11)]
    assert all(residues)
    # Norm(alpha)=63936/13 has odd 13-valuation, hence cannot be a square.
    assert 63936%13!=0
    lo,hi=F(6302063,10**6),F(6302064,10**6)
    assert f(-lo*lo)>0>f(-hi*hi)
    assert 39<lo*lo<hi*hi<40
    assert 13*39**2-480*39+2160>0 and 26*39-480>0
    clo,chi=F(62954,10**4),F(62956,10**4)
    g=lambda t:t**4-36*t*t-144
    assert g(clo)<0<g(chi)
    assert F(44,7)<clo<chi<lo
    denominator_checks=0
    for b in range(2,41):
        for a in range(-2*b,2*b+1):
            if gcd(a,b)!=1:continue
            value=f(-F(a,b)**2)
            assert value.denominator==b**6//gcd(b,13)
            imaginary_component=F(a,b)**21*value  # i^21=i
            assert imaginary_component.denominator==b**27//gcd(b,13)
            assert value and imaginary_component
            denominator_checks+=2
    result={'status':'passed','cases':[maincase,comparator,control],
            'total_exact_determinant_checks':sum(c['identity_checks'] for c in [maincase,comparator,control]),
            'cubic_mod_11_residues':residues,
            'new_root_bracket':[str(lo),str(hi)],
            'comparator_root_bracket':[str(clo),str(chi)],
            'new_certified_root_is_farther_than_comparator':True,
            'exact_rational_value_denominator_checks':denominator_checks,
            'scope':'Finite complete-source raw minors. Polynomial identities follow from proved degree bounds plus exact evaluations. No infinite coupled family, approximation improvement, or pi bad-approximability assertion.'}
    Path(__file__).with_name('coupled_raw_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))


if __name__=='__main__':main()
