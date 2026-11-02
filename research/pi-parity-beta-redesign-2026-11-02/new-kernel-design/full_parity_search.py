"""Bounded parity-free integer-antiderivative construction.

The declared generation features use rational samples and coefficient costs.
Only after generation are candidates checked against a rational pi interval.
"""
from fractions import Fraction as F
from math import comb,gcd,lcm
from pathlib import Path
import json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
from sympy.matrices.normalforms import hermite_normal_form
import sympy.polys.matrices.lll as lll_module
lll_module.mfloor=lambda q:q.numerator//q.denominator


def kernel(c,A,B):
    # Coordinates r_1,...,r_(2c-1),alpha; R(0)=0.
    M=s.zeros(2*c+1,2*c)
    for k in range(1,2*c):
        M[k-1,k-1]=k
        M[k+1,k-1]=k-2*c+2
    for j in range(c):M[2*j,2*c-1]=comb(c-1,j)
    rows=[list(M[j,:]) for j in range(A)]
    rows += [[sum(comb(j,k)*M[j,v] for j in range(k,2*c+1))
              for v in range(2*c)] for k in range(B)]
    C=s.Matrix(rows)
    normal,left,right=smith_normal_decomp(DomainMatrix.from_Matrix(C).convert_to(s.ZZ))
    V=right.to_Matrix()
    assert abs(V.det())==1
    rank=A+B
    K=V[:,rank:].T
    assert C*K.T==s.zeros(rank,2*c-rank)
    return M,C,K.lll(delta=F(3,4))


def features(c):
    out=[]
    count=max(9,c+2)
    for i in range(1,count+1):
        x=F(i,count+1);den=(1+x*x)**c
        row=[(k*x**(k-1)+(k-2*c+2)*x**(k+1))/den
             for k in range(1,2*c)]
        row.append(1/(1+x*x));out.append(row)
    return out


def atan(q,n=180):
    value=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(n)),F(0))
    return value,value+F(1,(2*n+1)*q**(2*n+1))


def run():
    a,b=atan(5);u,v=atan(239);plo,phi=16*a-4*v,16*b-4*u
    out=[];started=time.monotonic()
    for m in range(1,4):
        c,A,B=4*m+1,2*m,2*m
        M,C,K=kernel(c,A,B)
        period_matrix=s.Matrix([[int(r[-1]) for r in K.tolist()],
                                [sum(map(int,r[:-1])) for r in K.tolist()]])
        period_image=hermite_normal_form(period_matrix).tolist()
        ff=features(c);seen={};examined=zero=0
        for e in [0,8*m,12*m,16*m,20*m]:
            cost=F(1,2**e)
            feats=ff+[[F(0)]*(2*c-1)+[cost],
                      [cost]*(2*c-1)+[F(0)]]
            D=lcm(*(z.denominator for r in feats for z in r))
            fmat=s.Matrix([[int(z*D) for z in row] for row in feats]).T
            reduction=(K*fmat).row_join(K).lll(delta=F(3,4))
            for row in reduction.tolist():
                vec=list(map(int,row[-2*c:]));alpha=vec[-1];rr=vec[:-1]
                examined+=1
                assert C*s.Matrix(vec)==s.zeros(A+B,1)
                if not alpha:zero+=1;continue
                rat=F(-4*sum(rr),2**(c-1)*alpha)
                if rat<=0:continue
                pp,q=rat.numerator,rat.denominator
                n=list(map(int,M*s.Matrix(vec)))
                # Divide N=x^A(1-x)^B H exactly in the full x-coordinate.
                poly=s.Poly(sum(z*s.Symbol('x')**j for j,z in enumerate(n)),s.Symbol('x'))
                quotient,remainder=s.div(poly,s.Poly(s.Symbol('x')**A*(1-s.Symbol('x'))**B,s.Symbol('x')))
                assert remainder.is_zero
                hh=[int(quotient.nth(j)) for j in range(quotient.degree()+1)]
                norm=sum(map(abs,hh))
                cert=F(22,7)*F(q*q*norm,abs(alpha))*F(4,121)**m
                if rat<plo:L,U=(plo-rat)*q*q,(phi-rat)*q*q
                elif rat>phi:L,U=(rat-phi)*q*q,(rat-plo)*q*q
                else:raise AssertionError('Pi interval too wide')
                seen[(pp,q)]=dict(m=m,c=c,A=A,B=B,cost_exponent=e,
                     p=str(pp),q=str(q),alpha=str(alpha),R_coefficients=[str(z) for z in rr],
                     N_coefficients=[str(z) for z in n],H_coefficients=[str(z) for z in hh],
                     quality_lower=str(L),quality_upper=str(U),displayed_quality=float(U),
                     rational_norm_quality_upper=str(cert),H_l1=str(norm))
        rows=sorted(seen.values(),key=lambda r:F(r['quality_upper']))
        out.append(dict(m=m,integer_kernel_dimension=K.rows,rows_examined=examined,
                        rational_sample_count=max(9,c+2),
                        integer_kernel_basis=[[int(z) for z in row] for row in K.tolist()],
                        period_image_hnf=[[int(z) for z in row] for row in period_image],
                        zero_alpha_rows=zero,distinct_positive_candidates=len(rows),
                        best=rows[0],all_candidates=rows))
    result=dict(status='passed',generation='Declared rational samples and five fixed dyadic cost weights; no pi or target fraction enters candidate generation.',
       coverage='m=1,2,3 only',per_m=out,elapsed_seconds=time.monotonic()-started,
       scope='Finite construction comparison; sample selection is not an integral certificate. Quality is checked afterwards using a Machin rational interval.')
    Path(__file__).with_name('full_parity_search.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([dict(m=r['m'],dimension=r['integer_kernel_dimension'],fraction=r['best']['p']+'/'+r['best']['q'],quality=r['best']['displayed_quality'],norm_bound=r['best']['rational_norm_quality_upper']) for r in out],indent=2))


if __name__=='__main__':run()
