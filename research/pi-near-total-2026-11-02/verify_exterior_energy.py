"""Exact certificates for exterior column energy and the thermal gap model.

All arithmetic is integral or rational; the standard library suffices.
The full-jet model has theta=1 and zero boundary cost, so it is not a
counterexample to the separated-weight pi endpoint construction.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
import json


def determinant(matrix):
    if not matrix:
        return F(1)
    a=[[F(x) for x in row] for row in matrix]
    answer=F(1)
    for c in range(len(a)):
        pivot=next((i for i in range(c,len(a)) if a[i][c]),None)
        if pivot is None:
            return F(0)
        if pivot!=c:
            a[c],a[pivot]=a[pivot],a[c]
            answer=-answer
        diagonal=a[c][c]
        answer*=diagonal
        for i in range(c+1,len(a)):
            quotient=a[i][c]/diagonal
            for j in range(c+1,len(a)):
                a[i][j]-=quotient*a[c][j]
    return answer


def transpose(a):
    return [list(row) for row in zip(*a)]


def multiply(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def indices(d,height,strict=False):
    result=[a for a in product(range(height+1),repeat=d)
            if sum(a)<height if strict]
    if not strict:
        result=[a for a in product(range(height+1),repeat=d) if sum(a)<=height]
    return sorted(result,key=lambda a:(sum(a),a))


def direct_entry(row,column):
    s,beta=row[0],row[1:]
    h,alpha=column[0],column[1:]
    return comb(h,s) if beta==alpha and h>=s else 0


def ordinary_entry(row,column):
    s,beta=row[0],row[1:]
    h,alpha=column[0],column[1:]
    if any(a<b for a,b in zip(alpha,beta)):
        return F(0)
    degree=sum(a-b for a,b in zip(alpha,beta))
    if degree>s:
        return F(0)
    return F(prod(comb(a,b) for a,b in zip(alpha,beta))*h**(s-degree),
             factorial(s-degree))


def change_entry(old,new):
    s,beta=old[0],old[1:]
    r,gamma=new[0],new[1:]
    if any(g<b for g,b in zip(gamma,beta)):
        return F(0)
    degree=sum(g-b for g,b in zip(gamma,beta))
    k=s-degree
    if k<r or k<0:
        return F(0)
    exponential=sum((-1)**(r-j)*comb(r,j)*j**k for j in range(r+1))
    return F(prod(comb(g,b) for g,b in zip(gamma,beta))*exponential,
             factorial(k))


def geometric_pascal(s,t,q):
    # Products of two binomial polynomials are counted by overlap size l.
    answer=F(0)
    for l in range(min(s,t)+1):
        k=s+t-l
        multiplier=F(factorial(k),factorial(l)*factorial(s-l)*factorial(t-l))
        answer+=multiplier*q**k/(1-q)**(k+1)
    return answer


def direct_infinite_entry(row,other,q):
    if row[1:]!=other[1:]:
        return F(0)
    return q**sum(row[1:])*geometric_pascal(row[0],other[0],q)


def full_model_checks():
    records=[]
    changed_entries=0
    for d in (1,2,3):
        for n in (1,2,3,4):
            rows=indices(d,n,True)
            cols=indices(d,n)
            M=len(rows)
            assert M==comb(n+d-1,d)
            direct=[[F(direct_entry(row,col)) for col in cols] for row in rows]
            ordinary=[[ordinary_entry(row,col) for col in cols] for row in rows]
            change=[[change_entry(row,new) for new in rows] for row in rows]
            assert determinant(change)==1
            assert multiply(change,direct)==ordinary
            changed_entries+=M*len(cols)
            gram=multiply(ordinary,transpose(ordinary))
            finite=determinant(gram)
            beta_types=indices(d-1,n,True) if d>1 else [()]
            lengths=[n-sum(beta) for beta in beta_types]
            expected=prod(comb(2*l,l) for l in lengths)
            assert finite==expected
            assert finite<=4**M
            E0=sum(sum(row) for row in rows)
            T=sum(l*l for l in lengths)
            assert E0==F(d*(n-1),d+1)*M
            assert T==F(2*n+d-1,d+1)*M
            Emax=sum(sorted(map(sum,cols),reverse=True)[:M])
            assert Emax<=E0+d*M and Emax<=M*n
            q=F(1,n+1)
            infinite_direct=[[direct_infinite_entry(row,other,q) for other in rows]
                             for row in rows]
            infinite_ordinary=multiply(multiply(change,infinite_direct),transpose(change))
            infinite=determinant(infinite_ordinary)
            assert infinite==q**E0/(1-q)**T
            repaired=q**(-Emax)*infinite
            crude=q**(-M*n)*infinite
            assert finite<=repaired<=crude
            q_star=F(n+d,3*n+2*d-1)
            A=M*n-E0
            assert q_star==F(A,A+T)
            crude_opt=q_star**(-A)/(1-q_star)**T
            assert crude_opt<=crude
            records.append({'d':d,'H':n,'source_rows':M,'raw_columns':len(cols),
                            'finite_full_gram_determinant':str(finite),
                            'source_total_degree':E0,'sum_time_block_squares':T,
                            'maximum_allowed_exterior_energy':Emax,
                            'crude_thermal_optimal_q':str(q_star),
                            'crude_thermal_optimal_squared_bound':str(crude_opt),
                            'repaired_q':str(q),'repaired_squared_bound':str(repaired)})
    return records,changed_entries


def geometric_checks():
    count=0
    for l in range(1,7):
        for q in (F(1,4),F(1,3),F(2,5)):
            gram=[[geometric_pascal(s,t,q) for t in range(l)] for s in range(l)]
            assert determinant(gram)==q**(l*(l-1)//2)/(1-q)**(l*l)
            count+=1
    return count


def weighted_energy_checks():
    records=[]
    for weights in ((1,),(1,3),(1,2,5)):
        d=len(weights)
        for H in (12,24,48):
            costs=sorted(sum(w*a for w,a in zip(weights,nu))
                         for nu in product(*(range(H//w+1) for w in weights))
                         if sum(w*a for w,a in zip(weights,nu))<=H)
            for cut in (F(1,2),F(3,4)):
                selected=[c for c in costs if c>=cut*H]
                M=len(selected)
                assert sum(selected)==sum(costs[-M:])
                predicted=F(d,d+1)*(1-cut**(d+1))/(1-cut**d)
                records.append({'weights':weights,'H':H,'columns':len(costs),
                                'selected_columns':M,'cost_cutoff_fraction':str(cut),
                                'exact_maximum_energy':sum(selected),
                                'finite_normalized_maximum':str(F(sum(selected),M*H)),
                                'proved_limit_for_tail_fraction':str(predicted)})
    return records


if __name__=='__main__':
    model,entries=full_model_checks()
    result={'status':'passed','full_packet_models':model,
            'local_coordinate_entries':entries,'geometric_gram_blocks':geometric_checks(),
            'weighted_simplex_energy_cases':weighted_energy_checks(),
            'scope':'Exact finite full-jet polynomial kernels and exterior energy. '
                    'The model has equal source/raw weights, theta=1 and no '
                    'boundary escape; it is not the critical pi parameter family. '
                    'Infinite losses and repairs follow from the written formulas.'}
    Path(__file__).with_name('exterior_energy_verification.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:(len(v) if isinstance(v,list) else v)
                      for k,v in result.items()},indent=2))
