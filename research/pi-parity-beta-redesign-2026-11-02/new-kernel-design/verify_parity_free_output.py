"""Independent standard-library audit of saved full-R candidates and periods."""
from fractions import Fraction as F
from math import comb,gcd
from functools import reduce
from pathlib import Path
import json


def atan(q,n=180):
    v=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(n)),F(0))
    return v,v+F(1,(2*n+1)*q**(2*n+1))


def mul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):r[i+j]+=x*y
    return r


def numerator(rr,alpha,c):
    n=[0]*(2*c+1)
    for k,r in enumerate(rr,1):
        n[k-1]+=k*r
        n[k+1]+=(k-2*c+2)*r
    for j in range(c):n[2*j]+=alpha*comb(c-1,j)
    return n


def main():
    root=Path(__file__).resolve().parent
    data=json.loads((root/'full_parity_search.json').read_text())
    a,b=atan(5);u,v=atan(239);lo,hi=16*a-4*v,16*b-4*u
    count=periodchecks=0
    images=[]
    for block in data['per_m']:
        m=block['m'];c=4*m+1;A=B=2*m
        basis=block['integer_kernel_basis']
        assert len(basis)==4*m+2
        for row in basis:
            n=numerator(row[:-1],row[-1],c)
            assert not any(n[:A])
            assert all(sum(comb(j,k)*n[j] for j in range(k,len(n)))==0 for k in range(B))
        periods=[(r[-1],sum(r[:-1])) for r in basis]
        hnf=block['period_image_hnf']
        d1,d2=hnf[0][0],hnf[1][1]
        assert hnf==[[d1,0],[0,d2]]
        assert all(a%d1==0 and b%d2==0 for a,b in periods)
        minors=[abs(a*d-c0*b) for i,(a,b) in enumerate(periods) for c0,d in periods[i+1:]]
        assert reduce(gcd,minors)==d1*d2
        periodchecks+=1;images.append({'m':m,'image':hnf})
        for rec in block['all_candidates']:
            rr=[int(z) for z in rec['R_coefficients']]
            alpha=int(rec['alpha']);n=[int(z) for z in rec['N_coefficients']]
            assert numerator(rr,alpha,c)==n
            assert n[2*c-1]==0
            hh=[int(z) for z in rec['H_coefficients']]
            endpoint=[0]*A+[(-1)**j*comb(B,j) for j in range(B+1)]
            expected=mul(endpoint,hh)+[0]*(len(n)-len(endpoint)-len(hh)+1)
            assert expected==n
            rat=F(-4*sum(rr),2**(c-1)*alpha)
            pp,q=int(rec['p']),int(rec['q'])
            assert rat==F(pp,q) and gcd(pp,q)==1 and q>0
            if rat<lo:L,U=(lo-rat)*q*q,(hi-rat)*q*q
            elif rat>hi:L,U=(rat-hi)*q*q,(rat-lo)*q*q
            else:raise AssertionError('Need finer Machin interval')
            assert F(rec['quality_lower'])==L and F(rec['quality_upper'])==U
            assert F(rec['rational_norm_quality_upper'])==F(22*q*q*sum(map(abs,hh)),7*abs(alpha))*F(4,121)**m
            count+=1
    result=dict(status='passed',exact_candidate_checks=count,
                period_image_minor_checks=periodchecks,period_images=images,
                scope='Finite exact output audit; no claim of family optimality or new approximation record.')
    (root/'parity_free_output_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
