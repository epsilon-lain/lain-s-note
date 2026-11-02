"""Independent all-integer certificate: finite window plus positive tails.

Reconstructs the Gram score from kernel inputs, without SymPy, pi digits,
root isolation, or the proposed generator's stored winner.
"""
from fractions import Fraction as F
from math import comb, gcd
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent


def multiply(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return out


def value(a,x):
    out=F(0)
    for c in reversed(a): out=out*x+c
    return out


def shift(a,base,direction):
    # Coefficients of a(base+direction*x), exactly.
    return [sum((a[j]*comb(j,i)*base**(j-i)*direction**i
                 for j in range(i,len(a))),F(0)) for i in range(len(a))]


def verify_kernel(row,m,c):
    rr=list(map(int,row['R']));al=int(row['alpha']);hh=list(map(int,row['H']))
    n=[F(0)]*(c+1)
    for j,v in enumerate(rr):
        n[j]+=(2*j+1)*v
        n[j+1]+=(2*j-2*c+3)*v
    for j in range(c):n[j]+=al*comb(c-1,j)
    endpoint=[F(0)]*m+[F((-1)**j*comb(2*m,j)) for j in range(2*m+1)]
    assert n==multiply(endpoint,hh)


def main():
    data=json.loads((ROOT/'dynamics/inputs.json').read_text())
    m,c=data['m'],data['c'];first,second=data['rows']
    for row in (first,second):verify_kernel(row,m,c)
    scale=2**(c-3)
    raw=[-sum(map(int,first['R'])),scale*int(first['alpha']),
         -sum(map(int,second['R'])),scale*int(second['alpha'])]
    content=gcd(*raw)
    p0,a0,p1,a1=[r//content for r in raw]
    assert (p0,a0,p1,a1)==(313044,99645,-4948,-1575)
    # Unimodular coordinates give the *all k* gcd identity.
    assert 113*p1-355*a1==1
    assert 113*p0-355*a0==-3
    assert p1*a0-a1*p0==840
    D=840
    h0=list(map(int,first['H']));h1=list(map(int,second['H']))
    energy=[F(0)]*3
    for i in range(len(h0)):
        for j in range(len(h0)):
            entry=multiply([h0[i],h1[i]],[h0[j],h1[j]])
            for k,v in enumerate(entry):energy[k]+=v/F(2*i+2*j+1)
    factor=F(16*scale*scale,content*content)*F(2,27)**(2*m)
    raw_score=[factor*v for v in multiply(multiply([a0,a1],[a0,a1]),energy)]
    cap=128
    records=[]
    for k in range(-cap,cap+1):
        P,A=p0+k*p1,a0+k*a1
        if not A:continue
        g=gcd(P,A)
        assert g==gcd(k-3,D)
        score=value(raw_score,k)/g**4
        records.append((score,k,F(P,A)))
    score,k,rat=min(records)
    assert k==3 and rat==F(355,113)
    # On either infinite tail, F(k)>=best*D^4 follows coefficient by
    # coefficient after k=+cap+x or -cap-x, for every real x>=0.
    tails=[]
    for base,direction in [(cap,1),(-cap,-1)]:
        coefficients=shift(raw_score,base,direction)
        coefficients[0]-=score*D**4
        assert all(v>0 for v in coefficients)
        tails.append({'base':base,'direction':direction,
                      'positive_coefficients':[str(v) for v in coefficients]})
    assert sum(r[0]==score for r in records)==1
    # Entire maximal-divisor chain; reduced pairs have constant determinant 1.
    assert 355*(-1575)-113*(-4948)==-1
    for j in range(-100,101):
        assert gcd(355-4948*j,113-1575*j)==1
        assert F(p0+p1*(3+840*j),a0+a1*(3+840*j))==F(355-4948*j,113-1575*j)
    result={'status':'PASS','arithmetic':'standard-library Fraction only',
            'selection_uses_pi_or_CF':False,'all_integer_global_minimum':True,
            'scope':'All k in this fixed affine chain; not all approximants or a new record.',
            'score_coefficients':[str(v) for v in raw_score],
            'finite_window':[-cap,cap],'finite_cases':len(records),
            'best':{'k':k,'p':rat.numerator,'q':rat.denominator,'Q_bound_squared':str(score)},
            'infinite_tail_certificates':tails,
            'gcd_identity':'gcd(P(k),A(k))=gcd(k-3,840)',
            'maximal_gcd_chain':'k=3+840*j; (355-4948*j)/(113-1575*j)'}
    (ROOT/'independent_global_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['score_coefficients','infinite_tail_certificates']},indent=2))


if __name__=='__main__':main()
