"""Independent standard-library audit of saved exact Hilbert period lifts."""
from fractions import Fraction as F
from math import comb,gcd,lcm
from pathlib import Path
import json
import sys


def inner(a,b):
    return sum((F(x)*F(y)/F(i+j+1) for i,x in enumerate(a)
                for j,y in enumerate(b)),F(0))


def gs(rows):
    stars=[];lengths=[]
    for row in rows:
        v=list(map(F,row))
        for s,d in zip(stars,lengths):
            coefficient=inner(v,s)/d
            v=[x-coefficient*y for x,y in zip(v,s)]
        d=inner(v,v)
        assert d>0
        stars.append(v);lengths.append(d)
    return stars,lengths


def project(v,stars,lengths):
    out=list(map(F,v))
    for s,d in zip(stars,lengths):
        coefficient=inner(out,s)/d
        out=[x-coefficient*y for x,y in zip(out,s)]
    return out


def nearest(x):return (2*x.numerator+x.denominator)//(2*x.denominator)


def numerator(row,c):
    r,alpha=row[:-1],row[-1]
    out=[0]*(2*c+1)
    for k,value in enumerate(r,1):
        out[k-1]+=k*value
        out[k+1]+=(k-2*c+2)*value
    for j in range(c):out[2*j]+=alpha*comb(c-1,j)
    return out


def quotient(n,a,b):
    assert not any(n[:a])
    data=list(n[a:]);divisor=[(-1)**j*comb(b,j) for j in range(b+1)]
    out=[0]*(len(data)-b)
    for k in range(len(out)-1,-1,-1):
        value=data[k+b]//divisor[-1]
        assert value*divisor[-1]==data[k+b]
        out[k]=value
        for j,t in enumerate(divisor):data[k+j]-=value*t
    assert not any(data)
    return out


def atan(q,n=180):
    out=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(n)),F(0))
    return out,out+F(1,(2*n+1)*q**(2*n+1))


def main():
    forced='--forced' in sys.argv
    filename='hilbert_forced_period_lift.json' if forced else 'hilbert_quotient_lift.json'
    for ancestor in Path(__file__).resolve().parents:
        source=ancestor/'period-main-kernel'/filename
        if source.exists():break
    else:raise RuntimeError('Saved Hilbert period-lift data are required.')
    data=json.loads(source.read_text())
    a,b=atan(5);u,v=atan(239)
    lo,hi=16*a-4*v,16*b-4*u
    original=json.loads(source.with_name('hilbert_quotient_lift.json').read_text()) if forced else None
    blocks=candidates=0;indices=[];primitive=0
    for block in data['per_m']:
        m,c=block['m'],block['c']
        bb=block['period_lift_coordinates'];zz=block['zero_period_lift_coordinates']
        hh=block['period_lift_H'];zh=block['zero_period_H']
        pb=block['period_basis']
        assert len(bb)==2 and len(zz)==4*m
        assert pb==[[row[-1] for row in bb],[sum(row[:-1]) for row in bb]]
        assert pb[0][0]*pb[1][1]-pb[0][1]*pb[1][0]!=0
        D=1
        if forced:
            for k in range(1,2*m+1):D=lcm(D,k)
            while D%2==0:D//=2
            assert block['forced_period_factor']==D
            T=block['period_sublattice_transform'];p0=block['unrestricted_period_basis']
            assert abs(T[0][0]*T[1][1]-T[0][1]*T[1][0])==D==block['period_sublattice_index']
            assert gcd(gcd(*p0[1]),D)==1
            assert all(x%D==0 for x in p0[0])
            assert pb==[[sum(p0[i][k]*T[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
            assert all(x%D==0 for row in pb for x in row)
            prior=next(b for b in original['per_m'] if b['m']==m)
            assert prior['period_basis']==p0
            assert zz==prior['zero_period_lift_coordinates'] and zh==prior['zero_period_H']
            for current,key in [(bb,'period_lift_coordinates'),(hh,'period_lift_H')]:
                old=prior[key]
                expected=[[sum(T[i][j]*old[i][k] for i in range(2)) for k in range(len(old[0]))] for j in range(2)]
                assert current==expected
            oldS=[list(map(F,row)) for row in prior['schur_metric']]
            expectedS=[[sum(T[k][i]*oldS[k][l]*T[l][j] for k in range(2) for l in range(2)) for j in range(2)] for i in range(2)]
            assert expectedS==[list(map(F,row)) for row in block['schur_metric']]
            indices.append(D)
        for r,h in list(zip(bb,hh))+list(zip(zz,zh)):
            assert quotient(numerator(r,c),2*m,2*m)==h
        for row in zz:assert row[-1]==0 and sum(row[:-1])==0
        stars,lengths=gs(zh)
        assert lengths==list(map(F,block['orthogonal_zero_lengths']))
        residual=[project(row,stars,lengths) for row in hh]
        metric=[[inner(x,y) for y in residual] for x in residual]
        assert metric==[list(map(F,row)) for row in block['schur_metric']]
        assert metric[0][0]>0 and metric[0][0]*metric[1][1]>metric[0][1]**2
        radius=sum(lengths,F(0))/4
        for rec in block['candidates']:
            qcoord=rec['quotient_coordinates']
            h=[qcoord[0]*hh[0][j]+qcoord[1]*hh[1][j] for j in range(len(hh[0]))]
            r=[qcoord[0]*bb[0][j]+qcoord[1]*bb[1][j] for j in range(len(bb[0]))]
            chosen=[];errors=[]
            for hrow,rrow,s,d in reversed(list(zip(zh,zz,stars,lengths))):
                coordinate=-inner(h,s)/d
                z=nearest(coordinate)
                errors.append(coordinate-z);chosen.append(z)
                h=[x+z*y for x,y in zip(h,hrow)]
                r=[x+z*y for x,y in zip(r,rrow)]
            assert list(reversed(chosen))==rec['zero_period_correction']
            assert h==list(map(int,rec['H_coefficients']))
            assert r[:-1]==list(map(int,rec['R_coefficients']))
            alpha=int(rec['alpha']);assert r[-1]==alpha
            assert [alpha,sum(r[:-1])]==[sum(pb[k][j]*qcoord[j] for j in range(2)) for k in range(2)]
            if forced:
                content=0
                for x in r:content=gcd(content,x)
                period_gcd=gcd(alpha,sum(r[:-1]))
                raw_gcd=gcd(4*sum(r[:-1]),2**(c-1)*alpha)
                assert content==rec['integer_coordinate_content']==1
                assert period_gcd==rec['period_gcd'] and period_gcd%D==0
                assert raw_gcd==rec['raw_fraction_gcd'] and raw_gcd%D==0
                assert alpha%D==0 and sum(r[:-1])%D==0
                primitive+=1
            actual=inner(h,h)
            projected=sum(qcoord[i]*metric[i][j]*qcoord[j] for i in range(2) for j in range(2))
            correction=sum(e*e*d for e,d in zip(errors,reversed(lengths)))
            assert actual==projected+correction
            assert all(abs(e)<=F(1,2) for e in errors)
            assert correction<=radius
            assert actual==F(rec['hilbert_norm_squared'])
            assert projected==F(rec['projected_hilbert_norm_squared'])
            assert correction==F(rec['babai_rounding_cost_squared'])
            assert radius==F(rec['uniform_lift_radius_squared'])
            rat=F(-4*sum(r[:-1]),2**(c-1)*alpha)
            p,q=int(rec['p']),int(rec['q'])
            assert rat==F(p,q) and gcd(p,q)==1 and q>0
            if rat<lo:L,U=(lo-rat)*q*q,(hi-rat)*q*q
            elif rat>hi:L,U=(rat-hi)*q*q,(rat-lo)*q*q
            else:raise AssertionError('Need a finer pi interval.')
            assert F(rec['quality_lower'])<=L<=U<=F(rec['quality_upper'])
            certificate=F(72*q**4,7*alpha**2)*F(4,121)**(2*m)*actual
            assert certificate==F(rec['joint_quality_squared_upper']) and U*U<certificate
            candidates+=1
        blocks+=1
    result=dict(status='passed',forced_period_congruence=forced,
                hilbert_period_blocks=blocks,integer_nearest_plane_candidates=candidates,
                exact_congruence_indices=indices,primitive_forced_candidates=primitive,
                reproduced=['Hilbert Gram and Schur metric','orthogonal zero lengths','integer lifts and unchanged period pairs','exact rounding-cost decomposition','actual reduced q and rational quality intervals','72/7 joint squared-quality certificate'],
                scope='Saved-data audit. Completeness of integer kernel bases is supplied by the generator Smith/unimodular proof; no asymptotic quality claim.')
    destination=Path(__file__).with_name('verify_forced_hilbert_lifts_independent.json') if forced else Path(__file__).with_suffix('.json')
    destination.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
