"""Independent standard-library audit of saved Hilbert quotient lifts."""
from fractions import Fraction as F
from math import comb,gcd,lcm
from pathlib import Path
import json


def inner(a,b):
    return sum((F(x)*y/F(i+j+1) for i,x in enumerate(a) if x for j,y in enumerate(b) if y),F(0))


def gs(rows):
    stars=[];norms=[]
    for row in rows:
        v=list(map(F,row))
        for u,nu in zip(stars,norms):
            t=inner(v,u)/nu
            v=[x-t*y for x,y in zip(v,u)]
        nu=inner(v,v);assert nu>0
        stars.append(v);norms.append(nu)
    return stars,norms


def orthogonal(v,stars,norms):
    v=list(map(F,v))
    for u,nu in zip(stars,norms):
        t=inner(v,u)/nu
        v=[x-t*y for x,y in zip(v,u)]
    return v


def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def numerator(r,alpha,c):
    out=[0]*(2*c+1)
    for k,v in enumerate(r,1):
        out[k-1]+=k*v;out[k+1]+=(k-2*c+2)*v
    for j in range(c):out[2*j]+=alpha*comb(c-1,j)
    return out


def check_primitive(vec,h,m):
    c=4*m+1
    n=numerator(vec[:-1],vec[-1],c)
    expected=mul([0]*(2*m)+[(-1)**j*comb(2*m,j) for j in range(2*m+1)],h)
    expected +=[0]*(len(n)-len(expected))
    assert n==expected and n[2*c-1]==0


def atan(q,n=180):
    a=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(n)),F(0))
    return a,a+F(1,(2*n+1)*q**(2*n+1))


def check_file(filename):
    root=Path(__file__).resolve().parent
    data=json.loads((root/filename).read_text())
    aa,ab=atan(5);ba,bb=atan(239)
    lo,hi=16*aa-4*bb,16*ab-4*ba
    basis_count=candidate_count=0
    summaries=[]
    for block in data['per_m']:
        m=block['m'];c=4*m+1
        CP=block['period_lift_coordinates'];CH=block['period_lift_H']
        ZP=block['zero_period_lift_coordinates'];ZH=block['zero_period_H']
        periods=block['period_basis']
        original=block['unrestricted_period_basis'];T=block['period_sublattice_transform']
        force=block['forced_period_factor']
        expected=1
        if data['forced_period_congruence']:
            for k in range(1,2*m+1):expected=lcm(expected,k)
            while expected%2==0:expected//=2
        assert force==expected==block['period_sublattice_index']
        assert abs(T[0][0]*T[1][1]-T[0][1]*T[1][0])==force
        assert periods==[[sum(original[i][k]*T[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
        assert all(x%force==0 for row in periods for x in row)
        for j,(v,h) in enumerate(zip(CP,CH)):
            check_primitive(v,h,m)
            assert [v[-1],sum(v[:-1])]==[periods[0][j],periods[1][j]]
            basis_count+=1
        for v,h in zip(ZP,ZH):
            check_primitive(v,h,m)
            assert v[-1]==0 and sum(v[:-1])==0
            basis_count+=1
        assert len(ZH)==4*m
        stars,norms=gs(ZH)
        assert norms==[F(x) for x in block['orthogonal_zero_lengths']]
        projected=[orthogonal(h,stars,norms) for h in CH]
        S=[[inner(x,y) for y in projected] for x in projected]
        assert S==[[F(x) for x in row] for row in block['schur_metric']]
        assert S[0][0]>0 and S[0][0]*S[1][1]>S[0][1]**2
        radius=sum(norms,F(0))/4
        for rec in block['candidates']:
            a=rec['quotient_coordinates'];z=rec['zero_period_correction']
            vec=[sum(a[k]*CP[k][j] for k in range(2))+sum(t*r[j] for t,r in zip(z,ZP))
                 for j in range(2*c)]
            h=[sum(a[k]*CH[k][j] for k in range(2))+sum(t*r[j] for t,r in zip(z,ZH))
               for j in range(len(CH[0]))]
            assert h==list(map(int,rec['H_coefficients']))
            assert vec[:-1]==list(map(int,rec['R_coefficients']))
            assert vec[-1]==int(rec['alpha'])
            check_primitive(vec,h,m)
            alpha=vec[-1];period=[sum(periods[i][j]*a[j] for j in range(2)) for i in range(2)]
            assert period==[alpha,sum(vec[:-1])]
            assert alpha%force==0 and period[1]%force==0
            content=0
            for x in vec:content=gcd(content,abs(x))
            assert content==rec['integer_coordinate_content']
            assert rec['forced_period_factor']==force
            assert rec['period_gcd']==gcd(abs(alpha),abs(period[1]))
            assert rec['raw_fraction_gcd']==gcd(abs(4*period[1]),abs(2**(c-1)*alpha))
            proj=sum(F(a[i])*S[i][j]*a[j] for i in range(2) for j in range(2))
            errors=[inner(h,u)/nu for u,nu in zip(stars,norms)]
            assert all(abs(t)<=F(1,2) for t in errors)
            correction=sum(t*t*nu for t,nu in zip(errors,norms))
            actual=inner(h,h)
            assert actual==proj+correction and correction<=radius
            assert actual==F(rec['hilbert_norm_squared'])
            assert proj==F(rec['projected_hilbert_norm_squared'])
            assert correction==F(rec['babai_rounding_cost_squared'])
            assert radius==F(rec['uniform_lift_radius_squared'])
            rat=F(-4*sum(vec[:-1]),2**(c-1)*alpha)
            p,q=int(rec['p']),int(rec['q'])
            assert rat==F(p,q) and gcd(p,q)==1 and q>0
            if rat<lo:L,U=(lo-rat)*q*q,(hi-rat)*q*q
            elif rat>hi:L,U=(rat-hi)*q*q,(rat-lo)*q*q
            else:raise AssertionError('Need finer interval')
            assert F(rec['quality_lower'])<=L<U<=F(rec['quality_upper'])
            cert=F(72,7)*F(q**4,alpha**2)*F(4,121)**(2*m)*actual
            assert cert==F(rec['joint_quality_squared_upper'])
            assert U*U<cert
            candidate_count+=1
        summaries.append(dict(m=m,candidates=len(block['candidates']),zero_period_rank=len(ZH)))
    return dict(filename=filename,status='passed',primitive_basis_checks=basis_count,
             exact_hilbert_lift_checks=candidate_count,per_m=summaries,
             scope='Exact saved-output audit. Schur metrics, integer lift costs, periods, actual gcd and rational quality upper certificates are verified without SymPy.')


def main():
    root=Path(__file__).resolve().parent
    out=dict(status='passed',datasets=[check_file(name) for name in ['hilbert_quotient_lift.json','hilbert_forced_period_lift.json']])
    (root/'hilbert_output_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
