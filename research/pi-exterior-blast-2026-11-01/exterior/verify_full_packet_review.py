"""Independent convolution-based check of the full-packet modular witness.

No double Taylor-allocation formula from the originating probe is imported.
Ordinary time derivatives are normalized by multiplying their rows by s!;
this does not change nonvanishing over characteristic zero or this field.
"""
from math import comb, factorial
from pathlib import Path
import json

P=65537
HEIGHT=20
ROWS=[(j,s,(b1,b2)) for j in range(2) for b2 in range(3) for b1 in range(6) for s in range(HEIGHT) if s+4*b1+10*b2<HEIGHT]
COLS=[(h,(a1,a2)) for a2 in range(5) for a1 in range(11) for h in range(HEIGHT+1) if h+2*a1+5*a2<=HEIGHT]
TARGET=list(ROWS)
CHANGES=[((0,7,(3,0)),(0,6,(4,0))),((0,9,(0,1)),(0,8,(1,1)))]
for old,new in CHANGES:TARGET[TARGET.index(old)]=new
FACT=[factorial(k)%P for k in range(HEIGHT+1)]
INVFACT=[pow(k,-1,P) for k in FACT]


def convolution(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=(out[i+j]+x*y)%P
    return out


def entry(row,column,centers):
    j,s,beta=row;h,alpha=column
    if any(b>a for b,a in zip(beta,alpha)):return 0
    poly=[1]
    for base,a,b in zip(centers,alpha,beta):
        n=a-b;base=j*base%P
        factor=[comb(a,b)*comb(n,k)*pow(base,n-k,P)%P for k in range(n+1)]
        poly=convolution(poly,factor)
    value=sum(poly[k]*pow(h,s-k,P)*INVFACT[s-k] for k in range(min(s,len(poly)-1)+1))
    return value*FACT[s]%P


def lex_pivots(matrix):
    a=[row[:] for row in matrix]
    result=[]
    k=0
    for col in range(len(a[0])):
        pivot=next((j for j in range(k,len(a)) if a[j][col]%P),None)
        if pivot is None:continue
        a[k],a[pivot]=a[pivot],a[k]
        factor=pow(a[k][col],-1,P)
        a[k]=[x*factor%P for x in a[k]]
        # Reduced-row elimination is independent of the probe's echelon-only path.
        for row in range(len(a)):
            if row==k:continue
            factor=a[row][col]
            if factor:
                a[row]=[(x-factor*y)%P for x,y in zip(a[row],a[k])]
        result.append(col);k+=1
        if k==len(a):break
    return result


def determinant_mod(matrix):
    a=[row[:] for row in matrix];det=1
    for k in range(len(a)):
        pivot=next((row for row in range(k,len(a)) if a[row][k]%P),None)
        if pivot is None:return 0
        if pivot!=k:
            a[k],a[pivot]=a[pivot],a[k];det=-det
        p=a[k][k]%P;det=det*p%P
        inv=pow(p,-1,P)
        for row in range(k+1,len(a)):
            ratio=a[row][k]*inv%P
            if ratio:
                for col in range(k+1,len(a)):a[row][col]=(a[row][col]-ratio*a[k][col])%P
    return det%P


def test_fiber_coefficient():
    # At total order 10 source first-coordinate fibers have lengths 4 and 1.
    source=[(0,0),(1,0),(2,0),(3,0),(0,1)]
    target=[(0,0),(1,0),(2,0),(4,0),(1,1)]
    coeff=determinant_mod([[comb(a1,b1)*comb(a2,b2) if a1>=b1 and a2>=b2 else 0 for a1,a2 in target] for b1,b2 in source])
    delta=tuple(sum(a[i] for a in target)-sum(b[i] for b in source) for i in range(2))
    assert coeff==4 and delta==(2,0)
    assert len(set(TARGET))==len(TARGET)
    assert all(s>=0 and s+sum(a)==10 for _,s,a in [x[1] for x in CHANGES])
    assert all(s+4*a[0]+10*a[1]>=HEIGHT for _,s,a in [x[1] for x in CHANGES])
    assert all(2*a[0]+5*a[1]<=HEIGHT for _,s,a in TARGET)
    return {'integer_coefficient':coeff,'zeta_exponent':list(delta),'replaced_rows':2}


if __name__=='__main__':
    assert 256*256%P==P-1
    period=[[entry(row,col,(1,1)) for col in COLS] for row in TARGET]
    chosen=lex_pivots(period)
    assert len(chosen)==len(ROWS)==156 and len(COLS)==242
    rational_centers=(44*256*pow(7,-1,P)%P,710*256*pow(113,-1,P)%P)
    period_minor=[[row[c] for c in chosen] for row in period]
    source_minor=[[entry(row,COLS[c],rational_centers) for c in chosen] for row in ROWS]
    det_period=determinant_mod(period_minor)
    det_source=determinant_mod(source_minor)
    assert det_period and det_source
    result={'status':'passed','prime':P,'image_of_i':256,'source_rows':len(ROWS),'ambient_columns':len(COLS),'selected_column_indices':chosen,'period_target_determinant_at_formal_period_1_mod_p':det_period,'rational_source_determinant_mod_p':det_source,'transport':test_fiber_coefficient(),'scope':'Finite complete-packet witness. Nonzero period polynomial follows from its modular specialization; its actual-period value uses known transcendence. No high-dimensional endpoint inequality or separated-weight interpolation threshold is certified.'}
    Path(__file__).with_name('full_packet_review_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='selected_column_indices'},indent=2))
