"""Exact rank-two period quotient and nearest-plane integrand lifts.

Pi is absent from generation. A rational Machin interval is used afterwards.
"""
from fractions import Fraction as F
from math import comb, gcd, lcm
from pathlib import Path
import importlib.util
import json
import time
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
import sympy.polys.matrices.lll as lll_module
lll_module.mfloor=lambda q:q.numerator//q.denominator


for ancestor in Path(__file__).resolve().parents:
    options=[ancestor/'new-kernel-design/full_parity_search.py',
             ancestor/'pi-parity-beta-redesign-2026-11-02/new-kernel-design/full_parity_search.py']
    dependency=next((p for p in options if p.exists()),None)
    if dependency is not None:break
else:raise RuntimeError('Previous parity-free checkpoint is required.')
spec=importlib.util.spec_from_file_location('full_primitive',dependency)
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)


def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def nearest(x):return (2*x.numerator+x.denominator)//(2*x.denominator)


def gram_schmidt(rows):
    stars=[];norms=[]
    for row in rows:
        v=list(map(F,row))
        for u,nu in zip(stars,norms):
            coefficient=dot(v,u)/nu
            v=[x-coefficient*y for x,y in zip(v,u)]
        n=dot(v,v)
        assert n>0
        stars.append(v);norms.append(n)
    return stars,norms


def orthogonal_part(v,stars,norms):
    out=list(map(F,v))
    for u,nu in zip(stars,norms):
        coefficient=dot(out,u)/nu
        out=[x-coefficient*y for x,y in zip(out,u)]
    return out


def gauss_two(S):
    def inner(a,b):return sum(F(a[i])*S[i][j]*b[j] for i in range(2) for j in range(2))
    a,b=[1,0],[0,1]
    for _ in range(10000):
        if inner(b,b)<inner(a,a):a,b=b,a
        mu=nearest(inner(a,b)/inner(a,a))
        if mu==0:return [a,b]
        b=[x-mu*y for x,y in zip(b,a)]
    raise RuntimeError('Gauss reduction did not terminate')


def babai_lift(a,pp,Z,stars,norms,ncoord):
    v=[a[0]*pp[0][j]+a[1]*pp[1][j] for j in range(len(pp[0]))]
    coefficients=[]
    for row,u,nu in reversed(list(zip(Z,stars,norms))):
        k=nearest(-dot(v,u)/nu)
        coefficients.append(k)
        v=[x+k*y for x,y in zip(v,row)]
    original=list(map(int,v[-ncoord:]))
    return v,original,list(reversed(coefficients))


def run(maximum=3):
    lo,hi=base.atan(5);u,v=base.atan(239)
    pi_lo,pi_hi=16*lo-4*v,16*hi-4*u
    started=time.monotonic();blocks=[]
    for m in range(1,maximum+1):
        c=4*m+1;ncoord=2*c
        M,C,K=base.kernel(c,2*m,2*m)
        P=sp.Matrix([[int(row[-1]) for row in K.tolist()],
                     [sum(map(int,row[:-1])) for row in K.tolist()]])
        normal,left,right=smith_normal_decomp(DomainMatrix.from_Matrix(P).convert_to(sp.ZZ))
        assert normal==left*DomainMatrix.from_Matrix(P).convert_to(sp.ZZ)*right
        U=right.to_Matrix();B=U.T*K
        assert abs(U.det())==1
        periods=P*U
        assert periods[:,2:]==sp.zeros(2,K.rows-2)
        rows=[];seen={};count=0
        sample_count=max(9,K.rows+3,c+2)
        sample=[]
        for j in range(1,sample_count+1):
            x=F(j,sample_count+1);den=(1+x*x)**c
            sample.append([(k*x**(k-1)+(k-2*c+2)*x**(k+1))/den
                           for k in range(1,2*c)]+[1/(1+x*x)])
        for e in [0,8*m,12*m,16*m,20*m]:
            weight=F(1,2**e)
            ff=sample+[[F(0)]*(ncoord-1)+[weight],
                       [weight]*(ncoord-1)+[F(0)]]
            D=lcm(*(x.denominator for row in ff for x in row))
            E=sp.Matrix([[int(x*D) for x in row] for row in ff]).T
            embedded=(B*E).row_join(B)
            pp=[[int(x) for x in embedded.row(j)] for j in range(2)]
            Zmatrix=embedded[2:,:].lll(delta=F(3,4))
            Z=[[int(x) for x in row] for row in Zmatrix.tolist()]
            for row in Z:
                assert row[-1]==0 and sum(row[-ncoord:-1])==0
            stars,norms=gram_schmidt(Z)
            residuals=[orthogonal_part(p,stars,norms) for p in pp]
            S=[[dot(x,y) for y in residuals] for x in residuals]
            assert S[0][0]>0 and S[0][0]*S[1][1]-S[0][1]**2>0
            reduced=gauss_two(S);aa,bb=reduced
            choices=[aa,bb,[aa[0]+bb[0],aa[1]+bb[1]],
                     [aa[0]-bb[0],aa[1]-bb[1]]]
            for a in choices:
                lifted,vec,co= babai_lift(a,pp,Z,stars,norms,ncoord)
                assert list(sp.Matrix([vec])*E)+vec==lifted
                assert C*sp.Matrix(vec)==sp.zeros(4*m,1)
                alpha=vec[-1];r1=sum(vec[:-1])
                assert sp.Matrix([alpha,r1])==periods[:,:2]*sp.Matrix(a)
                if not alpha:continue
                rat=F(-4*r1,2**(c-1)*alpha)
                if rat<=0:continue
                count+=1
                n=[int(x) for x in M*sp.Matrix(vec)]
                # Exact endpoint quotient in Q[x] (known integral by monicity).
                x=sp.Symbol('x')
                quotient,remainder=sp.div(sp.Poly(sum(t*x**j for j,t in enumerate(n)),x),
                                           sp.Poly(x**(2*m)*(1-x)**(2*m),x))
                assert remainder.is_zero
                h=[int(quotient.nth(j)) for j in range(quotient.degree()+1)]
                q=rat.denominator
                if rat<pi_lo:L,Uq=(pi_lo-rat)*q*q,(pi_hi-rat)*q*q
                elif rat>pi_hi:L,Uq=(rat-pi_hi)*q*q,(rat-pi_lo)*q*q
                else:raise AssertionError('Need finer interval')
                seen[(rat.numerator,q)]={'m':m,'c':c,'cost_exponent':e,
                    'quotient_coordinates':a,'zero_period_correction':co,
                    'p':str(rat.numerator),'q':str(q),'alpha':str(alpha),
                    'R_coefficients':list(map(str,vec[:-1])),'N_coefficients':list(map(str,n)),
                    'H_coefficients':list(map(str,h)),
                    'quality_lower':str(L),'quality_upper':str(Uq),'displayed_Q':float(Uq),
                    'projected_cost_squared':str(dot([F(t) for t in lifted],[F(t) for t in lifted]))}
        candidates=sorted(seen.values(),key=lambda r:F(r['quality_upper']))
        blocks.append({'m':m,'c':c,'kernel_rank':K.rows,'zero_period_rank':len(Z),
                      'sample_count':sample_count,'nonzero_positive_lifts':count,
                      'distinct_candidates':len(candidates),'best':candidates[0] if candidates else None,
                      'candidates':candidates})
    result={'status':'passed','rule':'Smith period quotient, rational Schur metric, Gauss quotient reduction, Babai zero-period lift',
            'generation':'All matrices, sample points and cost weights are rational; pi is used only for subsequent quality checking.',
            'scales':list(range(1,maximum+1)),'elapsed_seconds':time.monotonic()-started,'per_m':blocks,
            'scope':'Bounded generation comparison, not global optimality or an asymptotic quality conclusion.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([dict(m=b['m'],q=b['best']['q'] if b['best'] else None,
                           Q=b['best']['displayed_Q'] if b['best'] else None,
                           candidates=b['distinct_candidates']) for b in blocks],indent=2))


if __name__=='__main__':run()
