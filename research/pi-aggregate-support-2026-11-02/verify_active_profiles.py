"""Exact finite checks for active spans, normal profiles, and quotient capacities.

Standard library only. Formal rational centers are never treated as periods.
The asymptotic results are proved in the note, not inferred from these checks.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial, prod, ceil
from pathlib import Path
import json

from verify_bivariate_blocks import determinant, rank, multiply, transpose


def inverse(matrix):
    n=len(matrix)
    a=[[F(x) for x in row]+[F(i==j) for j in range(n)]
       for i,row in enumerate(matrix)]
    for k in range(n):
        p=next(i for i in range(k,n) if a[i][k])
        a[k],a[p]=a[p],a[k]
        q=a[k][k];a[k]=[x/q for x in a[k]]
        for i in range(n):
            if i==k:continue
            q=a[i][k]
            if q:a[i]=[x-q*y for x,y in zip(a[i],a[k])]
    return [row[n:] for row in a]


def multiindices(m,maximum):
    return sorted((a for a in product(range(maximum+1),repeat=m)
                   if sum(a)<=maximum),key=lambda a:(sum(a),a))


@lru_cache(None)
def expand_linear(rows,exponents):
    """Coefficient dictionary of product_i (rows[i] dot x)^exponents[i]."""
    dimension=len(rows[0]);zero=(0,)*dimension
    out={zero:F(1)}
    for linear,power in zip(rows,exponents):
        for _ in range(power):
            nxt={}
            for a,value in out.items():
                for i,c in enumerate(linear):
                    if not c:continue
                    b=list(a);b[i]+=1;b=tuple(b)
                    nxt[b]=nxt.get(b,F(0))+value*c
            out={a:x for a,x in nxt.items() if x}
    return out


def poly_multiply(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def jet_entry(j,s,beta,h,alpha,centers,flow):
    """Time-factorial-scaled ordinary Taylor coefficient, exactly."""
    if any(b>a for b,a in zip(beta,alpha)):return F(0)
    polynomial=[F(1)]
    prefactor=1
    for b,a,c,f in zip(beta,alpha,centers[j],flow):
        prefactor*=comb(a,b)
        degree=a-b
        polynomial=poly_multiply(polynomial,[comb(degree,k)*c**(degree-k)*f**k
                                             for k in range(degree+1)])
    return prefactor*sum(value*F(factorial(s),factorial(s-k))*h**(s-k)
                         for k,value in enumerate(polynomial) if k<=s)


def raw_columns(m,w0,w,H):
    return [(h,a) for a in multiindices(m,H//w)
            for h in range((H-w*sum(a))//w0+1)]


def source_rows(m,K,v0,v,H):
    return [(j,s,b) for j in range(K)
            for b in multiindices(m,ceil(F(H)/v)-1)
            if v*sum(b)<H
            for s in range(ceil((H-v*sum(b))/v0))]


def coordinate_model(m,d):
    """Old X=L*(active T, normal U), with integral unimodular L."""
    if d==1:
        L=[[1]+[int(i==k+1) for k in range(m-1)] for i in range(m)]
        r=[2]*m
    else:
        L=[[1,i]+[int(i==k+2) for k in range(m-2)] for i in range(m)]
        r=list(range(m))
    return tuple(tuple(F(x) for x in row) for row in L),tuple(F(x) for x in r)


def exact_profile_checks():
    models=[(2,1,1,1,1,2), (3,1,2,F(1,2),2,3),
            (3,2,2,F(1,2),2,3), (4,2,2,F(1,4),2,3),
            (3,2,2,2,2,3)]
    records=[];entries=0;gram_checks=0
    for m,d,K,v0,v,H in models:
        w0=w=1;v0=F(v0);v=F(v)
        rows=source_rows(m,K,v0,v,H);cols=raw_columns(m,w0,w,H)
        L,r=coordinate_model(m,d);Linv=tuple(map(tuple,inverse(L)))
        centers=[tuple(j*x for x in r) for j in range(K)]
        flow=(F(1),)*m
        newcenters=[tuple(sum(x*y for x,y in zip(row,c)) for row in Linv)
                    for c in centers]
        newflow=tuple(sum(row) for row in Linv)
        assert all(not x for c in newcenters for x in c[d:])
        assert all(not x for x in newflow[d:])
        raw=[[jet_entry(j,s,b,h,a,centers,flow) for h,a in cols]
             for j,s,b in rows]
        adapted=[[jet_entry(j,s,b,h,a,newcenters,newflow) for h,a in cols]
                 for j,s,b in rows]
        # New coefficient rows are obtained by substituting u=L*eta.
        row_change=[[expand_linear(L,b).get(g,F(0)) if (j,s)==(jj,ss) else F(0)
                     for jj,ss,b in rows] for j,s,g in rows]
        # Adapted coordinate columns use Q=L^-1 X.
        column_change=[[expand_linear(Linv,a).get(b,F(0)) if h==hh else F(0)
                        for h,a in cols] for hh,b in cols]
        assert determinant(row_change)==1
        assert determinant(column_change)==1
        assert multiply(multiply(row_change,raw),column_change)==adapted
        entries+=len(rows)*len(cols)
        profiles=sorted({b[d:] for _,_,b in rows})
        block_ranks=[];capacity=0;visible_columns=0;block_grams=[]
        for nu in profiles:
            ri=[i for i,(_,_,b) in enumerate(rows) if b[d:]==nu]
            ci=[i for i,(_,a) in enumerate(cols) if a[d:]==nu]
            block=[[adapted[i][j] for j in ci] for i in ri]
            assert all(not adapted[i][j] for i in ri for j,(_,a) in enumerate(cols)
                       if a[d:]!=nu)
            block_ranks.append(rank(block));capacity+=min(len(ri),len(ci))
            visible_columns+=len(ci)
            block_grams.append(multiply(block,transpose(block)))
        actual=rank(raw)
        assert actual==rank(adapted)==sum(block_ranks)
        assert actual<=capacity<=min(len(rows),visible_columns)
        if len(rows)<=24:
            assert determinant(multiply(adapted,transpose(adapted)))==prod(
                determinant(g) for g in block_grams)
            gram_checks+=1
        records.append({'m':m,'active_dimension':d,'K':K,'H':H,
                        'v0':str(v0),'v':str(v),'source_rows':len(rows),
                        'all_columns':len(cols),'exact_rank':actual,
                        'sum_block_capacities':capacity,'visible_columns':visible_columns})
    return records,entries,gram_checks


def simplex_count(d,time_weight,transverse_weight,height,strict):
    if height<0 or (strict and height<=0):return 0
    count=0
    maximum=ceil(height/transverse_weight)-1 if strict else int(height//transverse_weight)
    for k in range(maximum+1):
        residual=height-transverse_weight*k
        times=ceil(residual/time_weight) if strict else int(residual//time_weight)+1
        count+=comb(k+d-1,d-1)*max(0,times)
    return count


def aggregated_counts(m,d,K,w0,w,v0,v,H):
    n=m-d;rows=capacity=visible=0
    for k in range(ceil(F(H)/v)):
        multiplicity=comb(k+n-1,n-1)
        r=K*simplex_count(d,v0,v,H-v*k,True)
        c=simplex_count(d,w0,w,H-w*k,False)
        rows+=multiplicity*r;capacity+=multiplicity*min(r,c);visible+=multiplicity*c
    assert rows==K*simplex_count(m,v0,v,F(H),True)
    return rows,capacity,visible


def beta_cdf(n,b,theta):
    normalization=F(factorial(n-1)*factorial(b-1),factorial(n+b-1))
    integral=sum((-1)**j*comb(b-1,j)*theta**(n+j)/F(n+j) for j in range(b))
    result=integral/normalization
    bernoulli=sum(comb(n+b-1,j)*theta**j*(1-theta)**(n+b-1-j)
                  for j in range(n,n+b))
    assert result==bernoulli
    return result


def asymptotic_formula_checks():
    count=0
    for m in range(3,26):
        for theta in (F(1,3),F(1,2),F(3,4)):
            one=theta**(m-1)*((1-theta)**2*m*m+(1-theta*theta)*m+2*theta)
            two=2*theta**(m-2)*(comb(m+1,3)*(1-theta)**3+
                 comb(m+1,2)*theta*(1-theta)**2+(m+1)*theta**2*(1-theta)+theta**3)
            assert one==2*beta_cdf(m-1,3,theta)
            assert two==2*beta_cdf(m-2,4,theta)
            if theta==F(1,2):
                assert one==F(m*m+3*m+4,2**(m+1))
                assert two==F(m**3+3*m*m+8*m+12,6*2**m)
            count+=2
    assert 2*beta_cdf(9,3,F(1,2))==F(67,1024)
    assert 2*beta_cdf(8,4,F(1,2))==F(29,128)
    records=[]
    for m,d in ((3,1),(4,2),(6,1),(6,2),(10,2)):
        K=2;theta=F(1,2);v0=2*K*theta**m
        for H in (8,16,32):
            M,cap,visible=aggregated_counts(m,d,K,F(1),F(1),v0,F(2),H)
            N=simplex_count(m,F(1),F(1),F(H),False)
            records.append({'m':m,'d':d,'H':H,'source_rows':M,'all_columns':N,
                            'finite_capacity_fraction':str(F(cap,M)),
                            'finite_visible_fraction':str(F(visible,M)),
                            'proved_visible_limit':str(2*beta_cdf(m-d,d+2,theta))})
    return count,records


def greedy_basis(forms,weights,reverse=False):
    selected=[]
    for i in sorted(range(len(weights)),key=lambda i:(weights[i],i),reverse=reverse):
        if rank([forms[j] for j in selected+[i]])>len(selected):selected.append(i)
    return selected


def greedy_checks():
    instances=[((1,2,5,9),(0,1,2,3)),((1,2,5,9),(0,0,1,1)),
               ((5,1,7,3),(0,2,2,0)),((1,4,4,9),(2,2,3,5)),
               ((1,3,7),(2,2,2))]
    records=[];support_checks=0;image_checks=0
    for weights,r in instances:
        d=1 if len(set(r))==1 else 2
        forms=tuple((F(1),) if d==1 else (F(1),F(x)) for x in r)
        low=greedy_basis(forms,weights);high=greedy_basis(forms,weights,True)
        low_inverse=inverse([forms[i] for i in low])
        high_inverse=inverse([forms[i] for i in high])
        low_coeff=multiply(forms,low_inverse);high_coeff=multiply(forms,high_inverse)
        for i in range(len(weights)):
            for q in range(d):
                if low_coeff[i][q]:assert weights[low[q]]<=weights[i]
                if high_coeff[i][q]:assert weights[i]<=weights[high[q]]
                support_checks+=2
        # Restrict the actual weighted monomials, and compute the image rank.
        for H in (3,5,7):
            raw=[(h,a) for a in multiindices(len(weights),H//min(weights))
                 if sum(w*x for w,x in zip(weights,a))<=H
                 for h in range(H-sum(w*x for w,x in zip(weights,a))+1)]
            image=[(h,a) for a in multiindices(d,H//min(weights))
                   for h in range(H-sum(a)+1)]
            matrix=[[expand_linear(forms,a).get(b,F(0)) if hh==h else F(0)
                     for hh,a in raw] for h,b in image]
            expected=sum(H-sum(weights[i]*x for i,x in zip(low,a))+1
                         for a in multiindices(d,H//min(weights))
                         if sum(weights[i]*x for i,x in zip(low,a))<=H)
            assert rank(matrix)==expected
            image_checks+=1
        records.append({'weights':weights,'r':r,'active_dimension':d,
                        'cheapest_pivots':low,'most_expensive_pivots':high,
                        'cheapest_weights':[weights[i] for i in low],
                        'most_expensive_weights':[weights[i] for i in high]})
    return records,support_checks,image_checks


if __name__=='__main__':
    models,entries,grams=exact_profile_checks()
    formulas,counts=asymptotic_formula_checks()
    greedy,support,images=greedy_checks()
    result={'status':'passed','exact_coordinate_entries':entries,
            'exact_profile_models':models,'adapted_gram_products':grams,
            'rational_beta_formula_checks':formulas,'finite_profile_counts':counts,
            'greedy_models':greedy,'greedy_filtration_checks':support,
            'exact_quotient_image_ranks':images,
            'scope':'Exact finite formal-center identities and ranks. All-profile splitting '
                    'uses equal transverse weights; unequal weights use only the proved quotient test. '
                    'No sufficient interpolation threshold or pi approximation claim.'}
    Path(__file__).with_name('active_profile_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:len(v) if isinstance(v,list) else v for k,v in result.items()},indent=2))
