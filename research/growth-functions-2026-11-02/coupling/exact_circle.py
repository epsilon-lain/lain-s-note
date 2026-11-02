"""Exact rational affine-circle generator, iterate splitter and cycle detector."""
from fractions import Fraction as Q

def floor(x):return x.numerator//x.denominator

def affine_at(x,omega,kappa):
    j=floor(x);u=x-j
    if u<=Q(1,4):return 1-kappa,omega+kappa*j
    if u<=Q(3,4):return 1+kappa,omega-kappa*(j+Q(1,2))
    return 1-kappa,omega+kappa*(j+1)

def F(x,omega,kappa):
    a,b=affine_at(x,omega,kappa);return a*x+b

def iterate_pieces(omega,kappa,q):
    omega,kappa=Q(omega),Q(kappa)
    if not 0<kappa<1 or not isinstance(q,int) or q<1:
        raise ValueError("Require 0 < kappa < 1 and positive integer q")
    pieces=[(Q(0),Q(1),Q(1),Q(0))]
    for _ in range(q):
        out=[]
        for lo,hi,A,B in pieces:
            v0,v1=A*lo+B,A*hi+B
            breaks=[lo,hi]
            for j in range(floor(v0)-1,floor(v1)+2):
                for d in [Q(1,4),Q(3,4)]:
                    t=(j+d-B)/A
                    if lo<t<hi:breaks.append(t)
            breaks=sorted(set(breaks))
            for l,h in zip(breaks,breaks[1:]):
                a,b=affine_at(A*(l+h)/2+B,omega,kappa)
                out.append((l,h,a*A,a*B+b))
        pieces=out
    return pieces

def cycles(omega,kappa,p,q):
    out=[]
    for lo,hi,A,B in iterate_pieces(omega,kappa,q):
        if A==1:
            if B==p:out.append((lo,hi,None,A))
            continue
        x=(p-B)/(A-1)
        if lo<=x<=hi:out.append((lo,hi,x,A))
    return out

if __name__=='__main__':
    for omega,p,q in [(Q(51,100),1,2),(Q(17,50),1,3)]:
        print({'omega':str(omega),'kappa':'1/2','p':p,'q':q,
               'exact_cycle_candidates':[(str(a),str(b),None if x is None else str(x),str(A))
                                        for a,b,x,A in cycles(omega,Q(1,2),p,q)]})
