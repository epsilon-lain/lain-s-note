"""Pell-synchronised Euler arctangent tails with a shared integer norm.

All candidate choices are integer / rational; pi interval used only for scoring.
"""
from fractions import Fraction as F
from math import factorial, gcd, lcm
from pathlib import Path
from decimal import Decimal, localcontext
import json

def coeff(n):
    return F(4**n*factorial(n)**2,factorial(2*n+1))

def euler(x,n):
    z=x*x/(1+x*x)
    h=x/(1+x*x)
    return h*sum((coeff(k)*z**k for k in range(n)),F())

def atan_interval(k,n):
    a=sum((F((-1)**j,(2*j+1)*k**(2*j+1)) for j in range(n)),F())
    b=a+F((-1)**n,(2*n+1)*k**(2*n+1))
    return min(a,b),max(a,b)

def pi_interval():
    a,b=atan_interval(5,90),atan_interval(239,30)
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]

def fmt(x):
    with localcontext() as c:
        c.prec=17
        return str(Decimal(x.numerator)/Decimal(x.denominator))

def quality(r,bounds):
    lo,hi=bounds
    e0=F(0) if lo<=r<=hi else min(abs(lo-r),abs(hi-r))
    e1=max(abs(lo-r),abs(hi-r))
    return e0*r.denominator**2,e1*r.denominator**2

def generate():
    rows=[]
    u,b=1,1
    for j in range(1,8):
        a=u+b
        D=a*a+b*b
        m=b*b+2*a*b-a*a
        k=a*a+2*a*b-b*b
        assert abs(m)==1 and m*m+k*k==2*D*D
        r,s,d=F(b,a),F(a-b,a+b),F(m,k)
        zr,zs,zd=(x*x/(1+x*x) for x in (r,s,d))
        assert zr-zs==F(m,2*D) and zd==F(1,2*D*D)
        for n in range(1,9):
            A=2**(n+1)*a*b**(2*n+1)
            B=(a*a-b*b)*(a-b)**(2*n)
            assert A!=B
            if n==1:
                assert B-A==-m*D
            if n%4==1:
                assert (B-A)%D==0
            wa,wb=A//gcd(A,B),B//gcd(A,B)
            sr,ss=euler(r,n),euler(s,n)
            for ell in range(1,7):
                sd=euler(d,ell)
                rr,rs=8*sr-4*sd,8*ss+4*sd
                combo=(B*rr-A*rs)/(B-A)
                naive=4*(sr+ss)
                den=lcm(rr.denominator,rs.denominator)
                p1=int(rr*den);p2=int(rs*den)
                raw_p=B*p1-A*p2
                raw_q=(B-A)*den
                gg=gcd(abs(raw_p),abs(raw_q))
                assert combo.denominator==abs(raw_q)//gg
                # All denominators are supported on 2, D, and coefficient costs.
                coeff_den=lcm(*(coeff(h).denominator for h in range(max(n,ell))))
                budget=coeff_den*(2*D)**max(n,2*ell)
                assert budget%rr.denominator==0 and budget%rs.denominator==0
                rows.append({'pell_index':j,'a':a,'b':b,'D':D,'m':m,'k':k,
                    'n':n,'delta_terms':ell,'A':str(A),'B':str(B),
                    'primitive_weights':[str(wa),str(wb)],
                    'fraction':str(combo),'p':str(combo.numerator),'q':str(combo.denominator),
                    'naive_fraction':str(naive),'first_fraction':str(rr),'second_fraction':str(rs),
                    'raw_p':str(raw_p),'raw_q':str(raw_q),'gcd':str(gg),
                    'q_bits':combo.denominator.bit_length(),
                    'common_denominator_budget_bits':budget.bit_length()})
        u,b=u+2*b,u+b
    return rows

def main():
    rows=generate()
    bounds=pi_interval()
    for row in rows:
        r=F(row['fraction'])
        q0,q1=quality(r,bounds)
        row.update(quality_interval=[str(q0),str(q1)],quality_display=fmt((q0+q1)/2))
        qs=quality(F(row['naive_fraction']),bounds)
        row['naive_quality_display']=fmt(sum(qs)/2)
        row['strictly_better_than_naive']=q1<qs[0]
    ranked=sorted(rows,key=lambda x:F(x['quality_interval'][1]))
    compact_rows=[]
    for row in rows:
        compact_rows.append({key:value for key,value in row.items() if key not in
            ['quality_interval','raw_p','raw_q','A','B','primitive_weights','first_fraction','second_fraction']})
    result={'scope':'Finite correlated tail cancellation, no new pi record or asymptotic quality claim.',
        'candidate_count':len(rows),'pi_interval':[str(x) for x in bounds],
        'all_generation_without_pi_digits':True,
        'strict_quality_improvements_over_naive':sum(x['strictly_better_than_naive'] for x in rows),
        'best':ranked[:10],'rows':compact_rows}
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['rows','pi_interval','best']},indent=2))
    print(json.dumps([{k:r[k] for k in ['pell_index','a','b','n','delta_terms','fraction','quality_display','naive_quality_display','strictly_better_than_naive']} for r in ranked[:8]],indent=2))

if __name__=='__main__':
    main()
