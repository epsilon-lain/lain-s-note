"""Self-contained saturated integer endpoint kernel, extracted from the earlier constructive family implementation."""
from fractions import Fraction as F
from math import comb
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
import sympy.polys.matrices.lll as lll_module
lll_module.mfloor=lambda value: value.numerator//value.denominator

def data(c,a,b):
    M=s.zeros(c+1,c+1)
    for j in range(c):
        M[j,j]=2*j+1
        M[j+1,j]=2*j-2*c+3
    for j in range(c): M[j,c]=comb(c-1,j)
    rows=[list(M[j,:]) for j in range(a)]
    rows += [[sum(comb(j,k)*M[j,v] for j in range(k,c+1)) for v in range(c+1)] for k in range(b)]
    C=s.Matrix(rows)
    normal,left,right=smith_normal_decomp(DomainMatrix.from_Matrix(C).convert_to(s.ZZ))
    T=right.to_Matrix();rank=a+b
    assert abs(T.det())==1
    B=T[:,rank:].T
    B=B.lll(delta=F(3,4))
    assert C*B.T==s.zeros(rank,c+1-rank)
    return M,B
