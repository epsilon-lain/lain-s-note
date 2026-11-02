"""Exact algebra and interval controls for shared multiplicative modulators.
Standard library only. No numerical trig or float accepts a certificate.
"""
from fractions import Fraction as Q
from math import factorial,comb,lcm
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def coeff(j):return Q(comb(2*j,j),4**j*(2*j+1))
def mulpoly(a,b,N):
    out=[Q(0)]*(N+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=N:out[i+j]+=x*y
    return out

def power(a,n,N):
    b=[Q(1)]+[Q(0)]*N
    for _ in range(n):b=mulpoly(b,a,N)
    return b

def det(a):
    a=[row[:] for row in a];out=Q(1)
    for i in range(len(a)):
        k=next(k for k in range(i,len(a)) if a[k][i])
        if k!=i:a[k],a[i]=a[i],a[k];out=-out
        p=a[i][i];out*=p
        for j in range(i+1,len(a)):
            f=a[j][i]/p
            for k in range(i,len(a)):a[j][k]-=f*a[i][k]
    return out

def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def mul(a,b):
    s=[x*y for x in a for y in b];return min(s),max(s)
def pt(x):return Q(x),Q(x)
def scale(a,b):return mul(a,pt(b))
def div(a,b):
    assert b[0]>0
    return mul(a,(1/b[1],1/b[0]))
def maxabs(a):return max(abs(a[0]),abs(a[1]))

def atan(q,N):
    a=sum((Q((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N)),Q(0))
    b=a+Q((-1)**N,(2*N+1)*q**(2*N+1))
    return min(a,b),max(a,b)
A,B=atan(5,70),atan(239,24)
EXACT_PI=add(scale(A,16),scale(B,-4))
# Outward dyadic rounding keeps interval arithmetic tractable and remains exact.
D=2**160
PI=(Q((EXACT_PI[0]*D).numerator//(EXACT_PI[0]*D).denominator,D),
    Q(-((-EXACT_PI[1]*D).numerator//(-EXACT_PI[1]*D).denominator),D))
assert PI[0]<=EXACT_PI[0]<EXACT_PI[1]<=PI[1]
TWOPI=scale(PI,2)
assert PI[1]<Q(22,7)

def sine(x,N=22):
    sq=mul(x,x);v=pt(Q((-1)**N,factorial(2*N+1)))
    for j in range(N-1,-1,-1):v=add(pt(Q((-1)**j,factorial(2*j+1))),mul(sq,v))
    v=mul(x,v);E=maxabs(x)**(2*N+3)/factorial(2*N+3)
    return add(v,(-E,E))

def error(e,T,r,eps=Q(0)):
    theta=scale(TWOPI,e*T);s=sine(theta);ss=mul(s,s)
    v=pt(coeff(r-1))
    for j in range(r-2,-1,-1):v=add(pt(coeff(j)),mul(ss,v))
    P=mul(s,v)
    # Shared gain perturbation +eps*sin(theta)**(2r+1)/(2*pi*T).
    h=add(pt(e),neg(div(P,scale(TWOPI,T))))
    if eps:
        sp=s
        for _ in range(r):sp=mul(sp,ss)
        h=add(h,scale(div(sp,scale(TWOPI,T)),eps))
    return h

def run():
    formal=[]
    for r in range(1,9):
        N=2*r+1
        S=[Q(0)]*(N+1)
        for j in range((N-1)//2+1):S[2*j+1]=Q((-1)**j,factorial(2*j+1))
        powers=[power(S,2*j+1,N) for j in range(r)]
        matrix=[[powers[j][2*i+1] for j in range(r)] for i in range(r)]
        assert det(matrix)==1
        residual=[Q(0),Q(1)]+[Q(0)]*(N-1)
        for j in range(r):
            for k in range(N+1):residual[k]-=coeff(j)*powers[j][k]
        assert all(residual[k]==0 for k in range(N))
        assert residual[N]==coeff(r)
        D=lcm(*(coeff(j).denominator for j in range(r)))
        formal.append({'r':r,'correction_order':2*r+1,'coefficients':[str(coeff(j)) for j in range(r)],
          'independent_condition_determinant':'1','first_error_coefficient':str(coeff(r)),
          'coefficient_common_denominator':str(D),'denominator_bits':D.bit_length(),
          'highest_sine_harmonic':2*r-1})
    # q=1, shared T in {1,2,3}; all errors lie in |e|<=1/(12*Tmax).
    Tmax=3;R=Q(1,12*Tmax);interval_checks=0;free_checks=0
    contraction=[]
    for r in range(1,5):
        C=Q(4,3)*coeff(r)*(2*Q(22,7)*Tmax)**(2*r)
        lam=C*R**(2*r)
        assert lam==Q(4,3)*coeff(r)*Q(11,21)**(2*r) and lam<1
        for T in (1,2,3):
            for j in range(1,13):
                e=R*Q(j,12);h=error(e,T,r)
                assert 0<h[0]<=h[1]<e
                assert maxabs(h)<=C*e**(2*r+1)
                interval_checks+=1
                for eps in (-coeff(r)/3,coeff(r)/3):
                    H=error(e,T,r,eps)
                    Cfree=(Q(4,3)*coeff(r)+abs(eps))*(2*Q(22,7)*Tmax)**(2*r)
                    assert maxabs(H)<=Cfree*e**(2*r+1)
                    assert Cfree*R**(2*r)<1
                    free_checks+=1
        contraction.append({'r':r,'radius':str(R),'uniform_contraction_bound':str(lam),
            'noise_budget':str(R*(1-lam))})
    # Switching T and integer jumps uses the same error equation; rational bounds
    # can therefore be iterated without trig. The indices are symbolic integers.
    switching=0
    for r in range(1,5):
        C=Q(4,3)*coeff(r)*(2*Q(22,7)*Tmax)**(2*r)
        a=R
        for k in range(5):
            expected=R*(C*R**(2*r))**(((2*r+1)**k-1)//(2*r))
            assert a==expected
            a=C*a**(2*r+1);switching+=1
    out={'status':'PASS','arithmetic':'standard-library exact rational formal series and interval bounds',
      'jet_designs':formal,'true_sine_interval_controls':interval_checks,
      'surviving_free_parameter_controls':free_checks,'switching_envelope_controls':switching,
      'uniform_tubes':contraction,
      'scope':'Shared modulator, finite commensurate scales and high-order all-time error bounds; no improved pi approximation theorem.'}
    (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='jet_designs'},indent=2))

if __name__=='__main__':run()
