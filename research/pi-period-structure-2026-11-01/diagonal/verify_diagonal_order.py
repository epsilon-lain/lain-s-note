"""Exact actual-period regrouping diagnostics; Python standard library only."""
from fractions import Fraction as Q
from math import comb, lcm, factorial
from itertools import permutations, product
from pathlib import Path
import json

# Two-variable polynomials in r1,r2.
def add(a,b):
    c=dict(a)
    for k,v in b.items():
        c[k]=c.get(k,Q(0))+v
        if not c[k]:del c[k]
    return c

def mul(a,b):
    c={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():
            key=(i+k,j+l);c[key]=c.get(key,Q(0))+v*w
    return {k:v for k,v in c.items() if v}

def scale(a,c):return {k:v*c for k,v in a.items() if v*c}

def power(a,n):
    out={(0,0):Q(1)}
    for _ in range(n):out=mul(out,a)
    return out

def evaluate(a,x,y):return sum(v*x**i*y**j for (i,j),v in a.items())

def detpoly(matrix):
    n=len(matrix);out={}
    for p in permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        term={(0,0):Q(sign)}
        for i,j in enumerate(p):term=mul(term,matrix[i][j])
        out=add(out,term)
    return out

def multiply_t(a,b,N):
    out=[{} for _ in range(N+1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=N:out[i+j]=add(out[i+j],mul(x,y))
    return out

def entry(j,s,beta,h,alpha):
    if any(a<b for a,b in zip(alpha,beta)):return {}
    out=[{(0,0):Q(comb(h,k))} if k<=h else {} for k in range(s+1)]
    factor=1
    for i,(a,b) in enumerate(zip(alpha,beta)):
        factor*=comb(a,b)
        constant={(1,0) if i==0 else (0,1):Q(j)} if j else {}
        base=[constant]+[{(0,0):Q((-1)**(k+1),k)} for k in range(1,s+1)]
        for _ in range(a-b):out=multiply_t(out,base,s)
    return scale(out[s],factor)

def packet(H,v0,weights,theta,K):
    rr=[]
    for j in range(K):
        for s in range(int(H/v0)+1):
            for beta in product(*(range(int(H*theta/w)+1) for w in weights)):
                if v0*s+sum(w*b for w,b in zip(weights,beta))/theta<H:
                    rr.append((j,s,beta))
    return rr

def matrix(rr,cols):
    return [[sum_poly(scale(entry(j,s,b,h,a),coef) for (h,a),coef in col.items()) for col in cols] for j,s,b in rr]

def sum_poly(items):
    out={}
    for p in items:out=add(out,p)
    return out

def normal_form(poly):
    # r1=z, r2=z+delta. Keys now (z degree, delta degree).
    out={}
    for (a,b),v in poly.items():
        for k in range(b+1):
            key=(a+b-k,k);out[key]=out.get(key,Q(0))+v*comb(b,k)
    return {k:v for k,v in out.items() if v}

def epsilon_expansion(poly):
    # r1=W+e1, r2=W+e2. Keys (W degree,e1 degree,e2 degree).
    out={}
    for (a,b),v in poly.items():
        for i in range(a+1):
            for j in range(b+1):
                key=(a+b-i-j,i,j)
                out[key]=out.get(key,Q(0))+v*comb(a,i)*comb(b,j)
    return {k:v for k,v in out.items() if v}

def norm_order(poly):return min(k for a,k in normal_form(poly)) if poly else None

def eps_order(poly):return min(i+j for a,i,j in epsilon_expansion(poly)) if poly else None

def show(poly):return [{'powers':list(k),'coefficient':str(v)} for k,v in sorted(poly.items())]

def check_degree_bound(rr,cols,weights):
    colcost=sum(max(sum(w*a for w,a in zip(weights,alpha)) for h,alpha in col) for col in cols)
    rowcost=sum(sum(w*b for w,b in zip(weights,beta)) for j,s,beta in rr)
    E=colcost-rowcost
    D=detpoly(matrix(rr,cols))
    if D:
        assert all(sum(w*a for w,a in zip(weights,powers))<=E for powers in D)
        normal=normal_form(D)
        assert all(weights[0]*z+weights[1]*d<=E for z,d in normal)
        assert norm_order(D)*weights[1]<=E
        assert norm_order(D)==eps_order(D)
    return D,E


def det_numeric(matrix):
    a=[[Q(x) for x in row] for row in matrix];value=Q(1)
    for k in range(len(a)):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:return Q(0)
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];value=-value
        x=a[k][k];value*=x
        for i in range(k+1,len(a)):
            c=a[i][k]/x
            for j in range(k+1,len(a)):a[i][j]-=c*a[k][j]
    return value


def check_coefficient_norm(rr,cols,D,K=2):
    B=[sum(abs(coef)*Q(3,2)**h*(4*K)**sum(a) for (h,a),coef in col.items()) for col in cols]
    checked=0
    for (j,s,beta),line in zip(rr,matrix(rr,cols)):
        for p,b in zip(line,B):
            assert sum(abs(v) for v in normal_form(p).values())<=2**s*b
            checked+=1
    bound=Q(factorial(len(rr)))*2**sum(s for j,s,beta in rr)
    for b in B:bound*=b
    assert sum(abs(v) for v in normal_form(D).values())<=bound
    return checked


def check_regrouping(rr,cols,D):
    # Explicitly sum determinant row choices, independently of expanding D.
    assert all(j in (0,1) and b==(0,0) for j,s,b in rr)
    maximum=[max(alpha[i] for col in cols for h,alpha in col) for i in range(2)]
    candidates=list(product(*(range(a+1) for a in maximum)))
    free=[i for i,(j,s,b) in enumerate(rr) if j==1]
    degree=sum(max(sum(alpha) for h,alpha in col) for col in cols)
    expanded=epsilon_expansion(D);cases=0
    for W in range(degree+1):
        grouped={}
        for selection in product(candidates,repeat=len(free)):
            rr_new=list(rr)
            for index,a in zip(free,selection):
                j,s,b=rr_new[index];rr_new[index]=(j,s,a)
            mat=matrix(rr_new,cols)
            coefficient=det_numeric([[evaluate(p,Q(W),Q(W)) for p in row] for row in mat])
            key=tuple(sum(a[i] for a in selection) for i in range(2))
            grouped[key]=grouped.get(key,Q(0))+coefficient
            cases+=1
        grouped={k:v for k,v in grouped.items() if v}
        expected={}
        for (z,i,j),v in expanded.items():expected[i,j]=expected.get((i,j),Q(0))+v*W**z
        assert grouped=={k:v for k,v in expected.items() if v}
    return cases

def main():
    raw=[{(0,(0,0)):Q(1)},{(0,(1,0)):Q(1)},
         {(1,(1,0)):Q(1)},{(0,(0,1)):Q(1)}]
    rr=packet(Q(2),Q(1),(Q(1),Q(1)),Q(1,3),2)
    assert rr==[(j,s,(0,0)) for j in range(2) for s in range(2)]
    D,E=check_degree_bound(rr,raw,(Q(1),Q(1)))
    assert D=={(2,0):Q(1),(1,1):Q(-1)}
    assert norm_order(D)==1
    assert all(h+sum(a)<=2 for col in raw for h,a in col)
    regroup_checks=check_regrouping(rr,raw,D)
    norm_checks=check_coefficient_norm(rr,raw,D)
    # Same raw family with separated weights in a valid finite packet.
    rr2=packet(Q(4),Q(2),(Q(1),Q(3)),Q(1,5),2)
    assert rr2==rr
    Ds,Es=check_degree_bound(rr2,raw,(Q(1),Q(3)))
    assert Ds==D and Es==5
    # Pure normal-factor mixed columns. Match the stated center-major row order.
    Z2={(0,(0,2)):Q(1),(0,(1,1)):Q(-2),(0,(2,0)):Q(1)}
    tZ2={}
    for (h,a),v in Z2.items():tZ2[1,a]=v;tZ2[0,a]=-v
    mixed=[{(0,(0,0)):Q(1)},{(1,(0,0)):Q(1),(0,(0,0)):Q(-1)},Z2,tZ2]
    rm=packet(Q(3),Q(2),(Q(1),Q(1)),Q(1,4),2)
    assert rm==rr
    for col in mixed:
        assert all(Q(1,4)*h+sum(a)<=3 for h,a in col)
    P,EP=check_degree_bound(rm,mixed,(Q(1),Q(1)))
    delta={(0,1):Q(1),(1,0):Q(-1)}
    assert P==power(delta,4)
    assert norm_order(P)==4 and normal_form(P)=={(0,4):Q(1)}
    regroup_checks+=check_regrouping(rm,mixed,P)
    norm_checks+=check_coefficient_norm(rm,mixed,P)
    # The degree-versus-order bound is also sharp with unequal weights.
    rm_sep=packet(Q(7),Q(4),(Q(1),Q(3)),Q(1,8),2)
    assert rm_sep==rr
    assert all(Q(1,4)*h+a[0]+3*a[1]<=7 for col in mixed for h,a in col)
    P_sep,E_sep=check_degree_bound(rm_sep,mixed,(Q(1),Q(3)))
    assert P_sep==P and E_sep==12 and norm_order(P_sep)*3==E_sep
    # Exhaust a small rational polynomial universe and certify order equivalence,
    # including higher-degree coefficient polynomials and nonhomogeneous cancellation.
    terms=[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
    count=0
    for coeffs in product((-1,0,1),repeat=len(terms)):
        p={k:Q(v) for k,v in zip(terms,coeffs) if v}
        if not p:continue
        assert norm_order(p)==eps_order(p)
        nf=normal_form(p)
        for x,y in [(Q(2),Q(3)),(Q(-1,2),Q(7,5))]:
            assert evaluate(p,x,y)==evaluate(nf,x,y-x)
        count+=1
    # Exact spacing after normal factor extraction: use the actual known fractions,
    # but no numerical assertion about their errors relative to pi is needed.
    fractions=[(3,1),(22,7),(355,113)]
    pairs=[]
    for i in range(len(fractions)):
        for j in range(i+1,len(fractions)):
            p1,q1=fractions[i];p2,q2=fractions[j];Qden=lcm(q1,q2)
            difference=Q(p2,q2)-Q(p1,q1)
            integer=Qden*difference
            assert integer.denominator==1 and integer
            absolute_det=16*difference**4
            assert absolute_det>=Q(16,Qden**4)
            pairs.append({'fractions':[[p1,q1],[p2,q2]],'lcm':Qden,
                          'cleared_difference':str(integer),'absolute_fourth_power_minor':str(absolute_det)})
    # Exact weak-anchor difference control for separated error magnitudes.
    errors=0
    for a,b in product(range(-12,13),repeat=2):
        e1=Q(a,10);e2=Q(b,20)
        if e1 and abs(e2)<=abs(e1)/2:
            delta_abs=2*abs(e2-e1)
            assert abs(e1)<=delta_abs<=3*abs(e1)
            errors+=1
    # A coefficient l1 upper bound, tested on exact rational W in place of omega.
    upper_checks=0
    for p in [D,P,{(3,0):Q(2),(2,1):Q(-4),(1,2):Q(2)}]:
        nf=normal_form(p);k=norm_order(p);degree=max(sum(t) for t in p)
        L=sum(abs(v) for v in nf.values())
        for W,e1,e2 in product((Q(-3),Q(2)),(Q(-1,8),Q(0),Q(1,8)),(Q(-1,8),Q(0),Q(1,8))):
            eta=max(abs(e1),abs(e2))
            bound=L*max(Q(1),abs(W)+1)**degree*(2*eta)**k
            assert abs(evaluate(p,W+e1,W+e2))<=bound
            upper_checks+=1
    out={'status':'passed','raw_minor':show(D),'raw_normal_form':show(normal_form(D)),
         'raw_actual_period_order':norm_order(D),'raw_weight_budget':str(E),
         'separated_weight_budget':str(Es),'normal_factor_minor':show(P),
         'normal_factor_order':norm_order(P),'normal_factor_weight_budget':str(EP),
         'entry_coefficient_norm_checks':norm_checks,
         'explicit_row_choice_regrouping_checks':regroup_checks,
         'unequal_weight_sharp_bound':str(E_sep),'polynomial_order_equivalence_cases':count,'rational_spacing_examples':pairs,
         'weak_anchor_error_bound_cases':errors,'coefficient_upper_bound_cases':upper_checks,
         'scope':'Exact actual-diagonal order and finite source-format examples; no interpolation-threshold certification or pi bad-approximability result.'}
    Path(__file__).with_name('diagonal_order_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
