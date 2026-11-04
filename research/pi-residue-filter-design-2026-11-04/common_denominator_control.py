"""Try the same residue-filter target with a shared type-II Pade denominator.

This is a separate control: the ordinary-arctan residue theorem is not assumed
to hold for these simultaneous approximants.
"""
from fractions import Fraction as F
from pathlib import Path
from decimal import Decimal, localcontext
import json
import sympy as sp

ROOT=Path(__file__).resolve().parent

def main():
    cs=[F(2,5),F(3,7),F(1)]
    weights=[35,35,-29]
    rows=[]
    for n in range(1,6):
        N=3*n
        fs=[[F((-1)**j,2*j+1)*c**(2*j+1) for j in range(N+n)] for c in cs]
        mat=[];rhs=[]
        for f in fs:
            for j in range(N,N+n):
                mat.append([sp.Rational(f[j-i]) for i in range(1,N+1)])
                rhs.append(-sp.Rational(f[j]))
        q=[F(1)]+[F(x) for x in sp.Matrix(mat).inv()*sp.Matrix(rhs)]
        Ps=[[sum((q[i]*f[j-i] for i in range(j+1)),F()) for j in range(N)] for f in fs]
        assert sum(q,F())!=0
        r=F(2,3)*sum((w*sum(P,F()) for w,P in zip(weights,Ps)),F())/sum(q,F())
        rows.append({'n':n,'N':N,'Q_coefficients':list(map(str,q)),
                     'P_coefficients':[list(map(str,P)) for P in Ps],
                     'p':str(r.numerator),'q':str(r.denominator),
                     'actual_q_bits':r.denominator.bit_length()})
    # Read a saved independent pi certificate only after generating every row.
    data=json.loads((ROOT/'results.json').read_text())
    lo,hi=map(F,data['pi_interval'])
    for row in rows:
        r=F(int(row['p']),int(row['q']));q=r.denominator
        assert r<lo or r>hi
        a,b=min(abs(lo-r),abs(hi-r))*q*q,max(abs(lo-r),abs(hi-r))*q*q
        with localcontext() as c:
            c.prec=16
            row['quality_display']=str(Decimal((a+b).numerator)/Decimal((a+b).denominator)/2)
        row['quality_strictly_above_one']=a>1
    (ROOT/'common_denominator_results.json').write_text(json.dumps({'rows':rows,
        'generation_reads_pi':False,'scope':'Five finite shared-denominator controls, no transferred residue theorem.'},indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ['n','actual_q_bits','quality_display']} for r in rows],indent=2))

if __name__=='__main__':main()
