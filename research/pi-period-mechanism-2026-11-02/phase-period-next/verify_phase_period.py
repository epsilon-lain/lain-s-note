"""Exact same-column period-rank certificates for complete phase-dilated packets.

Standard library only. Gaussian rational evaluation is reduced modulo a proved
prime; nonzero residues prove characteristic-zero determinants are nonzero.
Actual 2*pi*i nonvanishing then follows from rational polynomial coefficients.
No determinant magnitude, irrationality exponent, or asymptotic rank is claimed.
"""
from fractions import Fraction as F
from itertools import product
from collections import Counter,defaultdict
from math import comb,prod,isqrt,factorial
from pathlib import Path
import json

P=1000033
assert all(P%d for d in range(2,isqrt(P)+1))
I=next(pow(a,(P-1)//4,P) for a in range(2,100)
       if pow(pow(a,(P-1)//4,P),2,P)==P-1)
R=44*pow(7,-1,P)*I%P

def ceiling(x):return -(-x.numerator//x.denominator)
def fall(a,k):return prod(range(a-k+1,a+1))

def matrix_entry(row,col,Z):
    j,t,a=row;h,alpha=col
    if any(x<y for x,y in zip(alpha,a)):return 0
    d=sum(x-y for x,y in zip(alpha,a))
    pref=prod(fall(x,y) for x,y in zip(alpha,a))%P
    return pref*sum(comb(t,k)*pow(h,t-k,P)*(fall(d,k)%P)*
                    pow(j*Z,d-k,P) for k in range(min(t,d)+1))%P

def pivot_columns(A):
    a=[row[:] for row in A];r=0;cols=[]
    for c in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][c]%P),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        inv=pow(a[r][c],-1,P)
        a[r]=[(x*inv)%P for x in a[r]]
        for i in range(r+1,len(a)):
            u=a[i][c]%P
            if u:a[i]=[(x-u*y)%P for x,y in zip(a[i],a[r])]
        cols.append(c);r+=1
        if r==len(a):break
    return cols

def determinant_mod(A):
    a=[row[:] for row in A];n=len(a);det=1
    for c in range(n):
        r=next((i for i in range(c,n) if a[i][c]%P),None)
        if r is None:return 0
        if r!=c:a[c],a[r]=a[r],a[c];det=-det
        v=a[c][c]%P;det=det*v%P;inv=pow(v,-1,P)
        for i in range(c+1,n):
            u=a[i][c]*inv%P
            if u:
                for j in range(c+1,n):a[i][j]=(a[i][j]-u*a[c][j])%P
            a[i][c]=0
    return det%P

def bareiss(A):
    a=[row[:] for row in A];n=len(a);old=1;sign=1
    for k in range(n-1):
        r=next((i for i in range(k,n) if a[i][k]),None)
        if r is None:return 0
        if r!=k:a[k],a[r]=a[r],a[k];sign=-sign
        v=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=a[i][j]*v-a[i][k]*a[k][j]
                assert num%old==0
                a[i][j]=num//old
            a[i][k]=0
        old=v
    return sign*a[-1][-1] if n else 1

def make_packet(H,w0,w,theta):
    K=2;v0=F(1);L=2;Q=1;m=len(w)
    v=tuple(x/theta for x in w);A=L*theta;S=m*(L-1)
    assert v0==2*K*theta**m*w0
    assert w0>(S+2)*(m+1)*v0
    assert all(x>v0 for x in v) and A<1
    N0=ceiling(A*H/min(w))+S
    types=[b for b in product(*(range(ceiling(H/vi)) for vi in v))
           if sum(vi*bi for vi,bi in zip(v,b))<H]
    rows=[(j,s,b) for j in range(K) for b in types
          for s in range(ceiling((H-sum(vi*bi for vi,bi in zip(v,b)))/v0))]
    targets=[];buckets=defaultdict(list)
    for j,s,b in rows:
        N=s+sum(b);phase=tuple((N//L**i)%L for i in range(m))
        a=b if N<N0 else tuple(L*bi+ri for bi,ri in zip(b,phase))
        t=N-sum(a)
        assert t>=0 and all(ai>=bi for ai,bi in zip(a,b))
        assert sum(wi*ai for wi,ai in zip(w,a))<=H
        targets.append((j,t,a))
        if j==0:buckets[N].append(b)
    assert len(set(targets))==len(rows)
    cols=[(h,a) for a in product(*(range(int(H/wi)+1) for wi in w))
          if sum(wi*ai for wi,ai in zip(w,a))<=H
          for h in range(int((H-sum(wi*ai for wi,ai in zip(w,a)))/w0)+1)]
    return rows,targets,cols,buckets,N0,v

def check(H,w0,w,theta):
    rows,targets,ambient,buckets,N0,v=make_packet(H,F(w0),w,theta)
    M=len(rows);m=len(w)
    source=[[matrix_entry(r,c,R) for c in ambient] for r in rows]
    chosen_indices=pivot_columns(source)
    assert len(chosen_indices)==M
    columns=[ambient[i] for i in chosen_indices]
    source_minor=[[row[i] for i in chosen_indices] for row in source]
    target_minor=[[matrix_entry(r,c,R) for c in columns] for r in targets]
    ds=determinant_mod(source_minor);dt=determinant_mod(target_minor)
    assert ds and dt
    gamma=[sum(a[i] for _,_,a in targets)-sum(b[i] for _,_,b in rows) for i in range(m)]
    residual=[sum(a[i] for _,a in columns)-sum(a[i] for _,_,a in targets) for i in range(m)]
    assert all(x>=0 for x in gamma+residual)
    counts=Counter(a for _,a in columns)
    max_time=int(H/w0)+1
    slacks=[]
    for r in range(1,max_time+1):
        d=sum(max(k-r,0) for k in counts.values())
        source_tail=sum(s>=r for _,s,_ in rows)
        target_tail=sum(t>=r for _,t,_ in targets)
        assert d<=target_tail<=source_tail
        slacks.append(target_tail-d)
    transport_blocks=0
    for N,bs in buckets.items():
        bs=sorted(bs)
        if N<N0:continue
        phase=tuple((N//2**i)%2 for i in range(m))
        nodes=[tuple(2*b[i]+phase[i] for i in range(m)) for b in bs]
        mat=[[prod(comb(a[i],b[i]) if a[i]>=b[i] else 0 for i in range(m)) for a in nodes] for b in bs]
        value=bareiss(mat)
        # For integer dilation 2 every diagonal Newton factor is 2^|beta|.
        assert value==2**sum(sum(b) for b in bs)
        transport_blocks+=1
    outside=sum(t+sum(vi*ai for vi,ai in zip(v,a))>=H for _,t,a in targets)
    return {'H':H,'m':m,'K':2,'v0':'1','w0':str(w0),'w':list(map(str,w)),
            'theta':str(theta),'v':list(map(str,v)),'dilation':'2','phase_period':2**m,
            'source_rows':M,'ambient_columns':len(ambient),'selected_columns':len(columns),
            'source_minor_residue':ds,'target_minor_residue':dt,
            'outside_rows':outside,'outside_fraction':str(F(outside,M)),
            'source_types':len(set(b for _,_,b in rows)),
            'target_types':len(set(a for _,_,a in targets)),
            'gamma':gamma,'residual_coordinate_degree':residual,
            'time_cuts_checked':max_time,'minimum_cut_slack':min(slacks),
            'pascal_bucket_checks':transport_blocks,
            'selected_raw_columns':[[h,list(a)] for h,a in columns]}

if __name__=='__main__':
    cases=[check(H,7,(F(1),),F(1,28)) for H in (20,28,30,32,36,40,48,56,64)]
    cases += [check(H,64,(F(1,2),F(4)),F(1,16)) for H in (12,16,20,24,28,32,40)]
    data={'status':'passed','prime':P,'image_of_i':I,'image_of_44i_over_7':R,
          'evaluation_center':'44*i/7','row_scaling':'t! times product(a_i!)',
          'cases':cases,'scope':'Exact finite same-raw-column nonvanishing for source and phase target at Gaussian rational centers. Rational polynomial nonvanishing implies actual-period nonvanishing by transcendence. Not a magnitude bound, uniform-in-H theorem, high-dimensional theorem, or arithmetic-log-weight approximation construction.'}
    Path(__file__).with_name('phase_period_verification.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(cases),'prime':P,'rows':[c['source_rows'] for c in cases],
                      'outside_rows':[c['outside_rows'] for c in cases],
                      'residues':[[c['source_minor_residue'],c['target_minor_residue']] for c in cases]},indent=2))
