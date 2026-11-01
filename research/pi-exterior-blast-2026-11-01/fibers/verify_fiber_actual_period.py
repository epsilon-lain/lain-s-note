from math import comb,factorial
from fractions import Fraction
from pathlib import Path
import json
from itertools import product
P=65537
H=20
weights=(2,5)
rows=[(j,s,(b1,b2)) for j in range(2) for b2 in range(3) for b1 in range(6) for s in range(H) if s+4*b1+10*b2<H]
columns=[(h,(a1,a2)) for a2 in range(5) for a1 in range(11) for h in range(H+1) if h+2*a1+5*a2<=H]
target=list(rows)
for old,new in [((0,7,(3,0)),(0,6,(4,0))),((0,9,(0,1)),(0,8,(1,1)))]:
    target[target.index(old)]=new

def entry(row,column,centers):
    j,s,beta=row;h,alpha=column
    if any(b>a for a,b in zip(alpha,beta)):return 0
    value=0
    for k1 in range(min(s,alpha[0]-beta[0])+1):
        for k2 in range(min(s-k1,alpha[1]-beta[1])+1):
            power=s-k1-k2
            value+=comb(alpha[0],beta[0])*comb(alpha[1],beta[1])*comb(alpha[0]-beta[0],k1)*comb(alpha[1]-beta[1],k2)*(factorial(s)//factorial(power))*pow(h,power,P)*pow(j*centers[0]%P,alpha[0]-beta[0]-k1,P)*pow(j*centers[1]%P,alpha[1]-beta[1]-k2,P)
    return value%P

def pivots(matrix):
    a=[r[:] for r in matrix];out=[];row=0
    for c in range(len(a[0])):
        pivot=next((i for i in range(row,len(a)) if a[i][c]),None)
        if pivot is None:continue
        a[row],a[pivot]=a[pivot],a[row]
        inv=pow(a[row][c],-1,P)
        a[row]=[(x*inv)%P for x in a[row]]
        for i in range(row+1,len(a)):
            factor=a[i][c]
            if factor:
                a[i]=[(x-factor*y)%P for x,y in zip(a[i],a[row])]
        out.append(c);row+=1
        if row==len(a):break
    return out


def det_mod(matrix):
    a=[r[:] for r in matrix];out=1
    for c in range(len(a)):
        pivot=next((i for i in range(c,len(a)) if a[i][c]),None)
        if pivot is None:return 0
        if pivot!=c:
            a[c],a[pivot]=a[pivot],a[c]
            out=-out
        pc=a[c][c];out=out*pc%P
        inverse=pow(pc,-1,P)
        for i in range(c+1,len(a)):
            factor=a[i][c]*inverse%P
            if factor:
                for k in range(c+1,len(a)):
                    a[i][k]=(a[i][k]-factor*a[c][k])%P
    return out%P


def det_exact(matrix):
    a=[[Fraction(x) for x in row] for row in matrix];out=Fraction(1)
    for c in range(len(a)):
        pivot=next((i for i in range(c,len(a)) if a[i][c]),None)
        if pivot is None:return Fraction(0)
        if pivot!=c:a[c],a[pivot]=a[pivot],a[c];out=-out
        value=a[c][c];out*=value
        for i in range(c+1,len(a)):
            factor=a[i][c]/value
            for k in range(c+1,len(a)):
                a[i][k]-=factor*a[c][k]
    return out


def fiber_checks():
    count=0
    for c1,c2,height in ((1,3,8),(3,9,10),(2,7,15),(5,19,27)):
        fibers=[(y,(height-c2*y+c1-1)//c1) for y in range(height) if c2*y<height]
        source=[(x,y) for y,n in fibers for x in range(n)]
        for mask in range(1<<len(fibers)):
            chosen={y for i,(y,n) in enumerate(fibers) if (mask>>i)&1}
            target=[(x if x<n-1 or y not in chosen else n,y) for y,n in fibers for x in range(n)]
            matrix=[[comb(a,b)*comb(y,z) if a>=b and y>=z else 0 for a,y in target] for b,z in source]
            expected=1
            for y,n in fibers:
                if y in chosen:expected*=n
            assert det_exact(matrix)==expected
            assert sum(a for a,y in target)-sum(b for b,z in source)==len(chosen)
            assert all(a+y<=height for a,y in target)
            count+=1
    return count


if __name__=='__main__':
    assert P==65537 and all(P%d for d in range(2,257))
    assert 256*256%P==P-1
    assert len(rows)==156 and len(columns)==242 and len(set(target))==156
    assert all(h+2*a[0]+5*a[1]<=20 for h,a in columns)
    assert all(j in(0,1) and s+4*b[0]+10*b[1]<20 for j,s,b in rows)
    assert all(2*b[0]+5*b[1]<=20 for j,s,b in target)
    period=[[entry(r,c,(1,1)) for c in columns] for r in target]
    rational_centers=(44*256*pow(7,-1,P)%P,710*256*pow(113,-1,P)%P)
    source=[[entry(r,c,rational_centers) for c in columns] for r in rows]
    selected=pivots(period)
    assert len(selected)==156
    period_residue=det_mod([[line[i] for i in selected] for line in period])
    source_residue=det_mod([[line[i] for i in selected] for line in source])
    assert period_residue==21680 and source_residue==35592
    factorial_products=[]
    for packet in (target,rows):
        value=1
        for j,s,b in packet:value=value*factorial(s)%P
        factorial_products.append(value)
    normalized_period=period_residue*pow(factorial_products[0],-1,P)%P
    normalized_source=source_residue*pow(factorial_products[1],-1,P)%P
    assert normalized_period==57507 and normalized_source==36181
    beta=[(0,0),(1,0),(2,0),(3,0),(0,1)]
    nodes=[(0,0),(1,0),(2,0),(4,0),(1,1)]
    pascal=[[comb(a,b)*comb(y,z) if a>=b and y>=z else 0 for a,y in nodes] for b,z in beta]
    assert det_exact(pascal)==4
    assert sum(a for a,y in nodes)-sum(b for b,z in beta)==2
    assert sum(y for a,y in nodes)==sum(z for b,z in beta)
    result={'status':'passed','prime':P,'imaginary_unit_image':256,'source_rows':156,'available_columns':242,
        'selected_column_indices':selected,'selected_raw_columns':[{'h':columns[i][0],'alpha':list(columns[i][1])} for i in selected],
        'rational_centers':['44i/7','710i/113'],'formal_period_evaluation':1,
        'factorial_scaled_period_target_residue':period_residue,'factorial_scaled_rational_source_residue':source_residue,
        'period_target_factorial_product_mod_prime':factorial_products[0],'rational_source_factorial_product_mod_prime':factorial_products[1],
        'normalized_period_target_residue':normalized_period,'normalized_rational_source_residue':normalized_source,
        'normalized_full_transport_coefficient':'4*zeta_1^2','factorial_scaled_full_transport_coefficient':'252*zeta_1^2',
        'fiber_factorization_cases':fiber_checks(),
        'scope':'Nonzero modular residues certify a Gaussian-rational source minor and a nonzero formal-period polynomial on the same raw columns. Known transcendence transfers the latter to 2*pi*i. This supplies no lower bound on its magnitude and no high-dimensional endpoint result.'}
    Path(__file__).with_name('fiber_actual_period_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('selected_column_indices','selected_raw_columns')},indent=2))
