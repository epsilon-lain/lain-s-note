"""Exact raw-column determinants and certified roots close to 2*pi*i.

Standard library only. Decimal strings below propose rational root brackets;
every bracket is proved by exact sign and derivative interval checks. Pi is
enclosed independently by Machin's formula and alternating arctan series.
No numerical root finder or floating-point arithmetic is used here.
"""
from fractions import Fraction as F
from math import comb, factorial, prod
from pathlib import Path
import json


def trim(p):
    while len(p)>1 and p[-1]==0: p.pop()
    return p


def evaluate(p,x):
    y=F(0)
    for c in reversed(p): y=y*x+c
    return y


def interpolate(values):
    out=[F(0)]*len(values)
    basis=[F(1)]
    dif=list(map(F,values))
    for k in range(len(values)):
        for j,c in enumerate(basis):out[j]+=c*dif[0]
        dif=[b-a for a,b in zip(dif,dif[1:])]
        new=[F(0)]*(len(basis)+1)
        for j,c in enumerate(basis):
            new[j]-=c*k/(k+1)
            new[j+1]+=c/(k+1)
        basis=new
    return trim(out)


def det(matrix):
    a=[list(map(F,row)) for row in matrix]
    n=len(a)
    sign=1
    previous=F(1)
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:return F(0)
        if pivot!=k:
            a[k],a[pivot]=a[pivot],a[k]
            sign=-sign
        diagonal=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                a[i][j]=(a[i][j]*diagonal-a[i][k]*a[k][j])/previous
            a[i][k]=F(0)
        previous=diagonal
    return sign*a[-1][-1]


def raw_entry(j,s,h,a,z):
    # s! times the normalized Taylor coefficient.
    return sum(comb(a,k)*(factorial(s)//factorial(s-k))*h**(s-k)
               *(j*z)**(a-k) for k in range(min(a,s)+1))


def period_polynomial(n):
    rows=[(j,s) for j in (0,1) for s in range(n)]
    cols=[(h,a) for h in (0,1) for a in range(n)]
    H=2*n
    complete=[(j,s,b) for j in (0,1) for s in range(n+1)
              for b in range(2) if 2*s+2*n*b<H]
    assert complete==[(j,s,0) for j,s in rows]
    assert all(h+2*a<=H for h,a in cols)
    denominator=prod(factorial(s) for j,s in rows)
    bound=sum(a for h,a in cols)
    def at(z):
        return det([[raw_entry(j,s,h,a,z) for h,a in cols] for j,s in rows])/denominator
    p=interpolate([at(z) for z in range(bound+1)])
    for z in (-1,bound+2):assert evaluate(p,z)==at(z)
    assert all(c==0 for k,c in enumerate(p) if k<n or (k-n)%2)
    return p,denominator


def mul(A,B):
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]


def subtract(A,B):
    return [[a-b for a,b in zip(r,s)] for r,s in zip(A,B)]


def block_identity(n,p):
    E=[[F(1,factorial(s-a)) if s>=a else F(0) for a in range(n)] for s in range(n)]
    for z in (-2,-1,0,1,2):
        T=[[F(comb(a,s))*z**(a-s) if a>=s else F(0) for a in range(n)] for s in range(n)]
        assert det(subtract(mul(E,T),mul(T,E)))==evaluate(p,z)
    # Exact rank-one truncation defect for the multiplication/derivative pair.
    D=[[F(a) if a==s+1 else F(0) for a in range(n)] for s in range(n)]
    defect=subtract(subtract(mul(D,E),mul(E,D)),E)
    expected=[[F(-n,factorial(n-a)) if s==n-1 else F(0) for a in range(n)] for s in range(n)]
    assert defect==expected


def imul(a,b):
    xs=[x*y for x in a for y in b]
    return min(xs),max(xs)


def ieval(p,x):
    y=(F(0),F(0))
    for c in reversed(p):
        lo,hi=imul(y,x)
        y=lo+c,hi+c
    return y


def outward(x,digits=50):
    scale=10**digits
    return F((x[0]*scale).__floor__(),scale),F((x[1]*scale).__ceil__(),scale)


def arctan_interval(q,N=80):
    s=sum((F(1) if k%2==0 else F(-1))/((2*k+1)*q**(2*k+1)) for k in range(N))
    nxt=(F(1) if N%2==0 else F(-1))/((2*N+1)*q**(2*N+1))
    return min(s,s+nxt),max(s,s+nxt)


def pi_interval():
    a,b=arctan_interval(5),arctan_interval(239)
    return outward((16*a[0]-4*b[1],16*a[1]-4*b[0]),90)


def ss(x):return [str(a) for a in x]


LEFT={
 5:'6.295498999506229431413738468580566502705586150633543066822526322097014',
 6:'6.282871924499176268613861505823969274850092017660555222646702287818335',
 7:'6.283190917605573347253072450734694876381367776016005842822321191414658',
 8:'6.283185233359473468397202657690050424101413553835727806109722975629018',
 9:'6.283185307922172806659992299879340127680582707349139423390588393187889',
 10:'6.283185307173696882205229983212857020324607743902779418888496497111926',
}


def main():
    pi=pi_interval()
    omega=(2*pi[0],2*pi[1])
    records=[]
    for n in range(4,11):
        p,denominator=period_polynomial(n)
        block_identity(n,p)
        q=[p[k]*(-1)**((k-n)//2) if (k-n)%2==0 else F(0) for k in range(n,len(p))]
        # P(i*x)=i^n*x^n*q(x), q has rational real coefficients.
        if n==4:
            bracket=(F(6),F(6))
            assert evaluate(q,F(6))==0
        else:
            # Use 45 digits from the suggested decimal, then a one-unit enclosure.
            left=F(LEFT[n][:47])
            bracket=(left,left+F(1,10**45))
            assert evaluate(q,bracket[0])*evaluate(q,bracket[1])<0
        derivative=[k*q[k] for k in range(1,len(q))]
        slope=ieval(derivative,bracket)
        assert slope[0]>0 or slope[1]<0
        signed=(bracket[0]-omega[1],bracket[1]-omega[0])
        assert signed[0]>0 or signed[1]<0
        distance=(min(abs(signed[0]),abs(signed[1])),max(abs(signed[0]),abs(signed[1])))
        if n==10:assert distance[1]<F(6,10**12) and bracket[0]/distance[1]>10**12
        q_at=ieval(q,omega)
        abs_q=(min(abs(q_at[0]),abs(q_at[1])),max(abs(q_at[0]),abs(q_at[1])))
        assert q_at[0]>0 or q_at[1]<0
        value=(omega[0]**n*abs_q[0],omega[1]**n*abs_q[1])
        termwise=[abs(c) for c in p]
        absolute_sum=(evaluate(termwise,omega[0]),evaluate(termwise,omega[1]))
        cancellation_ratio=(value[0]/absolute_sum[1],value[1]/absolute_sum[0])
        if n==10:
            assert value[1]<F(13,10**17)
            assert cancellation_ratio[1]<F(18,10**33)
        records.append({
            'n':n,'H':2*n,'K':2,'m':1,'w0':1,'w1':2,'v0':2,'v1':2*n,'theta':str(F(1,n)),
            'source_rows':2*n,'raw_columns':[[h,a] for h in (0,1) for a in range(n)],
            'degree_surplus_D':n*(n-1),'degree':len(p)-1,
            'polynomial':{str(k):str(c) for k,c in enumerate(p) if c},
            'row_factorial_denominator':str(denominator),
            'root_imaginary_part_interval':ss(bracket),
            'signed_root_minus_two_pi':ss(outward(signed)),
            'distance_to_actual_period_interval':ss(outward(distance)),
            'root_separation_statistic_lower_bound':str((bracket[0]/distance[1]).__floor__()),
            'actual_period_determinant_absolute_interval':ss(outward(value)),
            'actual_value_over_termwise_absolute_sum_interval':ss(outward(cancellation_ratio)),
            'unique_root_in_bracket':True,'block_identity_checked':True,
            'rank_one_boundary_defect_checked':True,
        })
    result={'status':'passed','pi_interval':ss(pi),'cases':records,
            'scope':'Complete legal raw-column packets. Exact finite root proximity and actual-period '
                    'determinant cancellation. No infinite root-convergence theorem, critical half-density '
                    'high-dimensional construction, or rational approximation improvement is asserted.'}
    Path(__file__).with_name('raw_near_period_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','complete_packets':len(records),'root_certificates':len(records),
                     'commutator_evaluations':5*len(records),'boundary_defect_checks':len(records),
                     'n10_root_distance_upper':'6e-12','n10_root_statistic_lower':'1e12',
                     'n10_actual_determinant_upper':'1.3e-16','n10_cancellation_ratio_upper':'1.8e-32'},indent=2))


if __name__=='__main__':main()
