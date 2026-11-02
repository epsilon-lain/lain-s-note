"""Exact rank and Fourier-kernel checks for one shared periodic modulator.
Standard library only: no numerical roots, pi digits, or sampled trigonometry.
"""
from fractions import Fraction as F
from math import gcd,prod
from itertools import combinations
from pathlib import Path
from random import Random
import json
from dataclasses import dataclass

@dataclass(frozen=True,eq=False)
class GaussianInteger:
    real: int = 0
    imag: int = 0
    @staticmethod
    def coerce(v):
        return v if isinstance(v,GaussianInteger) else GaussianInteger(v,0)
    def __add__(self,v):
        v=self.coerce(v);return GaussianInteger(self.real+v.real,self.imag+v.imag)
    __radd__=__add__
    def __neg__(self):return GaussianInteger(-self.real,-self.imag)
    def __sub__(self,v):return self+-self.coerce(v)
    def __rsub__(self,v):return self.coerce(v)+-self
    def __mul__(self,v):
        v=self.coerce(v);return GaussianInteger(self.real*v.real-self.imag*v.imag,self.real*v.imag+self.imag*v.real)
    __rmul__=__mul__
    def __eq__(self,v):
        v=self.coerce(v);return self.real==v.real and self.imag==v.imag
    def conjugate(self):return GaussianInteger(self.real,-self.imag)

G=GaussianInteger


def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)

def power(a,h):
    out=[1]
    for _ in range(h):out=mul(out,a)
    return out

def monic_div(a,b):
    assert b[-1]==1
    out=[0]*max(1,len(a)-len(b)+1);work=list(a)
    for j in range(len(a)-1,len(b)-2,-1):
        c=work[j];out[j-len(b)+1]=c
        for i,v in enumerate(b):work[j-len(b)+1+i]-=c*v
    return trim(out),trim(work[:len(b)-1] or [0])

def divisors(n):return [d for d in range(1,n+1) if n%d==0]

def phi(n):return sum(gcd(k,n)==1 for k in range(1,n+1))

CYCL={}
def cyclo(n):
    if n in CYCL:return CYCL[n]
    p=[-1]+[0]*(n-1)+[1]
    for d in divisors(n)[:-1]:
        p,rem=monic_div(p,cyclo(d));assert rem==[0]
    CYCL[n]=p
    return p

def generator(orders):
    ds=sorted(set(d for r in orders for d in divisors(r)))
    a=[1]
    for d in ds:a=mul(a,cyclo(d))
    return a,ds

def inclusion_exclusion(orders):
    value=0
    for k in range(1,len(orders)+1):
        for subset in combinations(orders,k):value+=(-1)**(k+1)*gcd(*subset)
    return value

def rank(matrix):
    a=[[F(x) for x in row] for row in matrix]
    if not a:return 0
    row=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(row,len(a)) if a[i][j]),None)
        if pivot is None:continue
        a[row],a[pivot]=a[pivot],a[row]
        x=a[row][j];a[row]=[v/x for v in a[row]]
        for i in range(row+1,len(a)):
            if a[i][j]:
                y=a[i][j];a[i]=[v-y*w for v,w in zip(a[i],a[row])]
        row+=1
        if row==len(a):break
    return row

def residue_matrix(orders,N):
    return [[int(j%r==k) for j in range(N+1)] for r in orders for k in range(r)]

def rem_root(a,r):
    v=[0]*r
    for j,c in enumerate(a):v[j%r]+=c
    return trim(v)

def real_basis(a,M,h=1):
    b=power(a,h);D=len(b)-1;J=2*M-D;s=(-1)**h
    if J<0:return []
    qs=[]
    for j in range((J+1)//2):
        q=[G()]*(J+1);q[j]=G(1);q[J-j]=G(s);qs.append(q)
        q=[G()]*(J+1);q[j]=G(0,1);q[J-j]=G(0,-s);qs.append(q)
    if J%2==0:
        q=[G()]*(J+1);q[J//2]=G(1) if s==1 else G(0,1);qs.append(q)
    out=[]
    for q in qs:
        r=mul(b,q);r=r+[G()]*(2*M+1-len(r))
        out.append(r)
    return out

def run():
    rng=Random(20261102);order_sets=[[1],[2],[2,3],[2,3,6],[2,3,6,1],[2,3,6,5],[4,6],[6,10,15],[3,5,7]]
    order_sets += [rng.sample(range(1,13),rng.randrange(1,5)) for _ in range(20)]
    count=0;reconstructed=0;realchecks=0;degreechecks=0
    for orders in order_sets:
        a,ds=generator(orders);D=len(a)-1
        assert D==sum(map(phi,ds))==inclusion_exclusion(orders)
        assert list(reversed(a))==[-v for v in a]
        assert all(rem_root(a,r)==[0] for r in orders)
        degreechecks+=1
        for N in range(0,25):
            actual=rank(residue_matrix(orders,N))
            assert actual==min(D,N+1)
            count+=1
            if N>=D:
                q=[rng.randint(-5,5) for _ in range(N-D+1)]
                p=mul(a,q);recovered,rem=monic_div(p,a)
                assert rem==[0] and recovered==trim(q)
                assert all(rem_root(p,r)==[0] for r in orders)
                assert max(map(abs,p))<=sum(map(abs,a))*max(map(abs,q))
                reconstructed+=1
        for h in (1,2,3):
            for M in range(0,12):
                basis=real_basis(a,M,h)
                assert len(basis)==max(0,2*M+1-h*D)
                for p in basis:
                    assert all(p[j]==p[2*M-j].conjugate() for j in range(2*M+1))
                    # Contact multiplicity checked by exact Gaussian quotient.
                    aa=power(a,h);quo,rem=monic_div(p,aa)
                    assert rem==[0]
                    assert all(rem_root(p,r)==[0] for r in orders)
                    realchecks+=1
    # User's 133 denominator example and a positive real scalar construction.
    q=133;L=F(1);Ts=[F(1,2),F(1,3),F(1,6),F(1),F(1,5)]
    orders=[(q*T/L).denominator for T in Ts]
    assert orders==[2,3,6,1,5]
    example=[]
    for k in range(1,len(Ts)+1):
        a,ds=generator(orders[:k]);D=len(a)-1
        example.append({'T':[str(t) for t in Ts[:k]],'orders':orders[:k], 'D':D,'degree8_free_complex_coefficients':max(0,9-D)})
    a,_=generator([2,3]);assert a==[-1,-1,0,1,1]
    P=mul(a,[1,-2]);assert len(P)-1==5
    assert sum(map(abs,a))==4 and sum(map(abs,[1,-2]))==3
    assert sum(map(abs,P))<=12
    assert rem_root(P,2)==rem_root(P,3)==[0]
    out={'status':'PASS','exact_residue_matrix_rank_checks':count,
      'cyclotomic_and_inclusion_exclusion_degree_checks':degreechecks,
      'integer_kernel_reconstruction_and_height_checks':reconstructed,
      'real_fourier_basis_and_contact_checks':realchecks,
      'q133_example':example,
      'positive_modulator':{'A':a,'Q':[1,-2],'AQ':P,'epsilon':'1/1000','max_fourier_frequency':5,'uniform_g_interval':['1','143/125']},
      'scope':'Fixed-period finite Fourier submodel, with exact algebraic certificates; not a theorem about all arbitrary modulators or pi approximation.'}
    Path(__file__).with_name('checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='q133_example'},indent=2))

if __name__=='__main__':run()
