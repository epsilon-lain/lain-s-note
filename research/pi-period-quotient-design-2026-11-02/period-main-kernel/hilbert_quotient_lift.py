"""A proved Hilbert-norm objective on the rank-two period quotient.

The zero-period basis is integrally preconditioned using LLL. Quotient norms,
projection and the actual lift certificate use the exact rational Hilbert Gram.
"""
from quotient_lift import base,F,comb,gcd,lcm,Path,json,time,sp,DomainMatrix,smith_normal_decomp,nearest,gauss_two
import sys


def inner(a,b):
    return sum((F(x)*y/F(i+j+1) for i,x in enumerate(a) if x for j,y in enumerate(b) if y),F(0))


def quotient(n,A,B):
    assert not any(n[:A])
    data=list(map(int,n[A:])); divisor=[(-1)**j*comb(B,j) for j in range(B+1)]
    out=[0]*(len(data)-B)
    for k in range(len(out)-1,-1,-1):
        a=data[k+B]//divisor[-1]
        assert a*divisor[-1]==data[k+B]
        out[k]=a
        for j,t in enumerate(divisor):data[k+j]-=a*t
    assert not any(data)
    return out


def gs(rows):
    stars=[];norms=[]
    for row in rows:
        v=list(map(F,row))
        for u,nu in zip(stars,norms):
            a=inner(v,u)/nu
            v=[x-a*y for x,y in zip(v,u)]
        nu=inner(v,v);assert nu>0
        stars.append(v);norms.append(nu)
    return stars,norms


def project(v,stars,norms):
    out=list(map(F,v))
    for u,nu in zip(stars,norms):
        a=inner(out,u)/nu
        out=[x-a*y for x,y in zip(out,u)]
    return out


def bezout(a,b):
    if b==0:
        return (abs(a),1 if a>=0 else -1,0)
    g,x,y=bezout(b,a%b)
    return g,y,x-(a//b)*y


def run(maximum=3,forced=False):
    alo,ahi=base.atan(5);blo,bhi=base.atan(239)
    pi_lo,pi_hi=16*alo-4*bhi,16*ahi-4*blo
    blocks=[];start=time.monotonic()
    for m in range(1,maximum+1):
        c=4*m+1;ncoord=2*c;hcount=4*m+3
        M,C,K=base.kernel(c,2*m,2*m)
        P=sp.Matrix([[int(row[-1]) for row in K.tolist()],
                     [sum(map(int,row[:-1])) for row in K.tolist()]])
        pm=DomainMatrix.from_Matrix(P).convert_to(sp.ZZ)
        diag,left,right=smith_normal_decomp(pm)
        assert diag==left*pm*right
        U=right.to_Matrix();B=U.T*K;periods=P*U
        assert abs(U.det())==1 and periods[:,2:]==sp.zeros(2,K.rows-2)
        BB=[[int(x) for x in row] for row in B.tolist()]
        HH=[quotient([int(x) for x in M*sp.Matrix(row)],2*m,2*m) for row in BB]
        original_periods=[[int(periods[i,j]) for j in range(2)] for i in range(2)]
        force=1
        T=sp.eye(2)
        if forced:
            for k in range(1,2*m+1):force=lcm(force,k)
            while force%2==0:force//=2
            z0,z1=int(periods[1,0]),int(periods[1,1])
            g,u,v=bezout(z0,z1)
            assert gcd(g,force)==1
            T=sp.Matrix([[u*force,-z1//g],[v*force,z0//g]])
            assert abs(T.det())==force
            assert all(int(x)%force==0 for x in periods[:,:2]*T)
            oldBB=BB[:2];oldHH=HH[:2]
            BB[:2]=[[sum(int(T[i,j])*oldBB[i][k] for i in range(2)) for k in range(ncoord)] for j in range(2)]
            HH[:2]=[[sum(int(T[i,j])*oldHH[i][k] for i in range(2)) for k in range(hcount)] for j in range(2)]
            periods[:,:2]=periods[:,:2]*T
        assert all(len(h)==hcount for h in HH)
        # Integral preconditioning only. No claim that Euclidean LLL solves Hilbert CVP.
        Zembedding=sp.Matrix([HH[j]+BB[j] for j in range(2,len(HH))]).lll(delta=F(3,4))
        ZZ=[list(map(int,row)) for row in Zembedding.tolist()]
        ZH=[r[:hcount] for r in ZZ];ZB=[r[hcount:] for r in ZZ]
        for h,z in zip(ZH,ZB):
            assert z[-1]==0 and sum(z[:-1])==0
            assert quotient([int(t) for t in M*sp.Matrix(z)],2*m,2*m)==h
        stars,norms=gs(ZH)
        projected=[project(h,stars,norms) for h in HH[:2]]
        metric=[[inner(x,y) for y in projected] for x in projected]
        radius=sum(norms,F(0))/4
        seen={};nonzero=0
        for e in [0,8*m,12*m,16*m,20*m]:
            w=F(1,2**e)
            S=[[metric[i][j]+w*w*sum(int(periods[k,i])*int(periods[k,j]) for k in range(2))
                for j in range(2)] for i in range(2)]
            pair=gauss_two(S);a,b=pair
            for qcoord in [a,b,[a[0]+b[0],a[1]+b[1]],[a[0]-b[0],a[1]-b[1]]]:
                h=[qcoord[0]*HH[0][j]+qcoord[1]*HH[1][j] for j in range(hcount)]
                vec=[qcoord[0]*BB[0][j]+qcoord[1]*BB[1][j] for j in range(ncoord)]
                errors=[];co=[]
                for zh,zb,u,nu in reversed(list(zip(ZH,ZB,stars,norms))):
                    coordinate=-inner(h,u)/nu
                    z=nearest(coordinate);errors.append(coordinate-z);co.append(z)
                    h=[x+z*y for x,y in zip(h,zh)];vec=[x+z*y for x,y in zip(vec,zb)]
                alpha=vec[-1];r1=sum(vec[:-1])
                assert alpha%force==0 and r1%force==0
                assert sp.Matrix([alpha,r1])==periods[:,:2]*sp.Matrix(qcoord)
                assert C*sp.Matrix(vec)==sp.zeros(4*m,1)
                assert h==quotient([int(t) for t in M*sp.Matrix(vec)],2*m,2*m)
                projcost=sum(F(qcoord[i])*metric[i][j]*qcoord[j] for i in range(2) for j in range(2))
                actual=inner(h,h)
                correction=sum(er*er*nu for er,nu in zip(errors,reversed(norms)))
                assert all(abs(er)<=F(1,2) for er in errors)
                assert actual==projcost+correction and actual<=projcost+radius
                if not alpha:continue
                rat=F(-4*r1,2**(c-1)*alpha)
                if rat<=0:continue
                nonzero+=1;qn=rat.denominator
                if rat<pi_lo:L,Uq=(pi_lo-rat)*qn*qn,(pi_hi-rat)*qn*qn
                elif rat>pi_hi:L,Uq=(rat-pi_hi)*qn*qn,(rat-pi_lo)*qn*qn
                else:raise AssertionError('Need finer interval')
                certificate=F(72,7)*F(qn**4,alpha**2)*F(4,121)**(2*m)*actual
                assert Uq*Uq<certificate
                scale=10**24
                saved_lower=F((L*scale).numerator//(L*scale).denominator,scale)
                value=Uq*scale
                saved_upper=F(-((-value.numerator)//value.denominator),scale)
                assert saved_lower<=L<Uq<=saved_upper
                seen[(rat.numerator,qn)]={'m':m,'c':c,'cost_exponent':e,
                    'quotient_coordinates':qcoord,'zero_period_correction':list(reversed(co)),
                    'p':str(rat.numerator),'q':str(qn),'alpha':str(alpha),
                    'forced_period_factor':force,
                    'integer_coordinate_content':int(sp.igcd(*vec)),
                    'period_gcd':gcd(abs(alpha),abs(r1)),
                    'raw_fraction_gcd':gcd(abs(4*r1),abs(2**(c-1)*alpha)),
                    'R_coefficients':list(map(str,vec[:-1])),'H_coefficients':list(map(str,h)),
                    'hilbert_norm_squared':str(actual),'projected_hilbert_norm_squared':str(projcost),
                    'babai_rounding_cost_squared':str(correction),'uniform_lift_radius_squared':str(radius),
                    'joint_quality_squared_upper':str(certificate),
                    'quality_lower':str(saved_lower),'quality_upper':str(saved_upper),'displayed_Q':float(Uq)}
        candidates=sorted(seen.values(),key=lambda r:F(r['quality_upper']))
        blocks.append({'m':m,'c':c,'kernel_rank':K.rows,'zero_period_rank':len(ZH),
                       'period_basis':[[int(periods[i,j]) for j in range(2)] for i in range(2)],
                       'unrestricted_period_basis':original_periods,
                       'forced_period_factor':force,
                       'period_sublattice_transform':[[int(T[i,j]) for j in range(2)] for i in range(2)],
                       'period_sublattice_index':force,
                       'period_lift_coordinates':BB[:2],'period_lift_H':HH[:2],
                       'zero_period_lift_coordinates':ZB,'zero_period_H':ZH,
                       'schur_metric':[[str(x) for x in row] for row in metric],
                       'orthogonal_zero_lengths':[str(x) for x in norms],
                       'nonzero_positive_lifts':nonzero,'distinct_candidates':len(candidates),
                       'best':candidates[0] if candidates else None,'candidates':candidates})
    result={'status':'passed','rule':'Exact rank-two Hilbert quotient; integer zero-period nearest-plane lifts',
            'generation':'Rational Hilbert Gram plus fixed dyadic alpha/R1 costs; no pi input. Machin intervals are post-checks only.',
            'zero_period_preconditioning':'Integer Euclidean LLL on H coefficients with original-coordinate guards; exact Hilbert projection follows.',
            'elapsed_seconds':time.monotonic()-start,'per_m':blocks,
            'scope':'Bounded m=1,2,3 construction; exact lift radius and joint certificates are proved, but no asymptotic quality statement is inferred.'}
    result['forced_period_congruence']=forced
    destination=Path(__file__).with_name('hilbert_forced_period_lift.json') if forced else Path(__file__).with_suffix('.json')
    destination.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([dict(m=b['m'],q=b['best']['q'] if b['best'] else None,
                           Q=b['best']['displayed_Q'] if b['best'] else None,
                           candidates=b['distinct_candidates']) for b in blocks],indent=2))


if __name__=='__main__':run(forced='--force-period' in sys.argv)
