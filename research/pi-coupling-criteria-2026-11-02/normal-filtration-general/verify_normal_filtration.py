"""Exact standard-library checks for normal profiles and straightening.

Infinite statements are proved in the companion note. These finite checks
audit raw coefficients, changes of coordinates, support cuts, and legal packets.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
import json


def det(matrix):
    a=[list(map(F,row)) for row in matrix]
    result=F(1)
    for i in range(len(a)):
        pivot=next((j for j in range(i,len(a)) if a[j][i]),None)
        if pivot is None:return F(0)
        if pivot!=i:a[i],a[pivot]=a[pivot],a[i];result=-result
        diagonal=a[i][i];result*=diagonal
        for j in range(i+1,len(a)):
            ratio=a[j][i]/diagonal
            for k in range(i+1,len(a)):a[j][k]-=ratio*a[i][k]
    return result


def ge(a,b):
    return all(x>=y for x,y in zip(a,b))


def coefficient(h,degree,s,x):
    if degree<0:return F(0)
    return sum(F(comb(degree,k))*x**(degree-k)*F(h**(s-k),factorial(s-k))
               for k in range(min(s,degree)+1))


def raw(row,col,z):
    j,s,nu,b=row;h,kappa,a=col
    if a<b or not ge(kappa,nu):return F(0)
    factor=prod(comb(k,v) for k,v in zip(kappa,nu))*comb(a,b)
    return factor*coefficient(h,a-b+sum(kappa)-sum(nu),s,j*z)


def coordinate_checks():
    row_checks=column_checks=0
    for d in (1,2):
        normals=list(product(range(3),repeat=d))
        for nu in normals:
            for kappa in normals:
                for b in range(3):
                    for s in (0,1):
                        for j in (0,1):
                            a,h,z=3,2,F(3,2)
                            row=(j,s,nu,b)
                            value=F(0)
                            for eta in product(range(b+1),repeat=d):
                                if sum(eta)>b:continue
                                shifted=tuple(x+y for x,y in zip(nu,eta))
                                factor=prod(comb(x+y,x) for x,y in zip(nu,eta))
                                value+=factor*raw((j,s,shifted,b-sum(eta)),(h,kappa,a),z)
                            direct=F(0)
                            if ge(kappa,nu):
                                degree=a+sum(kappa)-sum(nu)
                                if degree>=b:
                                    direct=prod(comb(k,v) for k,v in zip(kappa,nu))*comb(degree,b)*coefficient(h,degree-b,s,j*z)
                            assert value==direct
                            row_checks+=1
                            value=F(0)
                            for lam in product(*(range(k+1) for k in kappa)):
                                shift=sum(kappa)-sum(lam)
                                factor=(-1)**shift*prod(comb(k,l) for k,l in zip(kappa,lam))
                                value+=factor*raw(row,(h,lam,a+shift),z)
                            direct=F(0)
                            shift=sum(kappa)-sum(nu)
                            if ge(kappa,nu) and 0<=b-shift<=a:
                                direct=(-1)**shift*prod(comb(k,v) for k,v in zip(kappa,nu))*comb(a,b-shift)*coefficient(h,a-b+shift,s,j*z)
                            assert value==direct
                            column_checks+=1
    return row_checks,column_checks


def matching_size(rows,cols):
    owner={}
    def augment(i,seen):
        for j,kappa in enumerate(cols):
            if j in seen or not ge(kappa,rows[i]):continue
            seen.add(j)
            if j not in owner or augment(owner[j],seen):
                owner[j]=i
                return True
        return False
    return sum(augment(i,set()) for i in range(len(rows)))


def support_checks():
    cases=zero_surplus=0
    labels=list(product(range(2),repeat=2))
    upper=[]
    for mask in range(1<<len(labels)):
        subset={i for i in range(len(labels)) if mask>>i&1}
        if all(j in subset for i in subset for j in range(len(labels)) if ge(labels[j],labels[i])):
            upper.append(subset)
    for r in product(range(3),repeat=len(labels)):
        for c in product(range(3),repeat=len(labels)):
            if sum(r)!=sum(c) or not 0<sum(r)<=5:continue
            rows=[x for x,count in zip(labels,r) for _ in range(count)]
            cols=[x for x,count in zip(labels,c) for _ in range(count)]
            M=len(rows)
            cuts=[M-sum(r[i] for i in U)+sum(c[i] for i in U) for U in upper]
            match=matching_size(rows,cols)
            assert match<=min(cuts)
            assert (match==M)==all(sum(r[i] for i in U)<=sum(c[i] for i in U) for U in upper)
            surplus=sum(map(sum,cols))-sum(map(sum,rows))
            if match==M:assert surplus>=0
            if surplus==0:
                assert (match==M)==(r==c)
                zero_surplus+=1
            cases+=1
    return cases,len(upper),zero_surplus


def closure_checks():
    count=0
    # Complete one-normal packets, at integer heights. Check exact finite
    # closure and its sufficient weight ordering, including active base jets.
    for H in range(2,9):
        for vT,vX in product(range(1,5),repeat=2):
            rows={(s,nu,b) for s,nu,b in product(range(H+1),repeat=3)
                  if s+vX*nu+vT*b<H}
            closed=all((s,nu+k,b-k) in rows for s,nu,b in rows for k in range(b+1))
            if vT>=vX:assert closed
            cols={(h,a,k) for h,a,k in product(range(H+1),repeat=3)
                  if h+vT*a+vX*k<=H}
            cclosed=all((h,a+k-lam,lam) in cols for h,a,k in cols for lam in range(k+1))
            if vT<=vX:assert cclosed
            count+=1
    return count


def legal_family_checks():
    cases=[];raw_checks=row_changes=0
    for N in range(2,6):
        normal=range((N-1)//2+1)
        rows=[(j,s,(nu,),b) for nu in normal for j in (0,1)
              for s in range(N-2*nu) for b in range(N-2*nu-s)]
        expected={(j,s,(nu,),b) for j in (0,1) for s,nu,b in product(range(N+1),repeat=3)
                  if s+b+2*nu<N}
        assert set(rows)==expected
        cols=[(h,(nu,),a) for nu in normal for h in range(N-2*nu)
              for a in range(2*(N-2*nu-h))]
        assert len(rows)==len(cols)
        assert all(F(h,4)+F(a,3)+F(2*kappa[0],3)<N for h,kappa,a in cols)
        assert all(sum(row[2]==(nu,) for row in rows)==sum(col[1]==(nu,) for col in cols) for nu in normal)
        exponent=sum(c[2] for c in cols)-sum(r[3] for r in rows)
        assert exponent==sum(sum(k*k for k in range(1,N-2*nu+1)) for nu in normal)
        constant=det([[raw(row,col,F(1)) for col in cols] for row in rows])
        assert constant
        for z in (F(1),F(2),F(-3,2)):
            actual=det([[raw(row,col,z) for col in cols] for row in rows])
            assert actual==constant*z**exponent
            product_blocks=F(1)
            for nu in normal:
                rr=[r for r in rows if r[2]==(nu,)]
                cc=[c for c in cols if c[1]==(nu,)]
                product_blocks*=det([[raw(row,col,z) for col in cc] for row in rr])
                for row in rr:
                    j,s,norm,b=row
                    for h,kappa,a in cc:
                        transformed=sum((-1)**k*comb(b+k,k)*raw((j,s-k,norm,b+k),(h,kappa,a),z)
                                        for k in range(s+1))
                        direct=F(h**s,factorial(s))*comb(a,b)*(j*z)**(a-b) if a>=b else F(0)
                        assert transformed==direct
                        row_changes+=1
            assert actual==product_blocks
            raw_checks+=1
        cases.append({'N':N,'matrix_size':len(rows),'normal_blocks':len(normal),
                      'active_base_derivatives':any(r[3]>0 for r in rows),
                      'monomial_constant':str(constant),'monomial_exponent':exponent})
    return cases,raw_checks,row_changes


def main():
    row,col=coordinate_checks()
    profiles,upper,zero=support_checks()
    cases,raw_count,changes=legal_family_checks()
    result={'status':'passed','row_coordinate_checks':row,'column_coordinate_checks':col,
            'normal_profile_support_checks':profiles,'upper_sets_in_test_poset':upper,
            'zero_surplus_profiles':zero,'complete_weighted_closure_cases':closure_checks(),
            'legal_complete_packets':cases,'raw_monomial_determinant_checks':raw_count,
            'time_to_base_row_change_entries':changes,
            'scope':'Exact finite checks supplement general support/closure/monomial proofs. '
                    'No universal nonvanishing assertion, no critical-density construction, '
                    'and no pi approximation claim.'}
    Path(__file__).with_name('normal_filtration_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
