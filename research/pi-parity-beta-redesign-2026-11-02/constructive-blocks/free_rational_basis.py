"""Small saturated-kernel control for the rational-basis construction."""
from rational_basis_truncation import *
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
from sympy import ZZ
import sympy.polys.matrices.lll as lll_module
lll_module.mfloor=lambda x:x.numerator//x.denominator

def run():
    lo,hi=rational_pi();rows=[]
    for m in range(1,4):
      c=4*m+1
      for K in [2*m,c]:
       for d in sorted(set([m+2,2*m+1,3*m+1])):
        mu=[]
        for j in range(d+1):
          aa=bb=F(0)
          for k in range(K+1):
            a,b=rational_moment(2*m+j+k,2*m+k,c+k)
            f=comb(c+k-1,k);aa+=f*a;bb+=f*b
          mu.append((aa,bb))
        rr=[]
        for col in range(2):
          v=[x[col] for x in mu];L=lcm(*(x.denominator for x in v))
          row=[int(x*L) for x in v];g=reduce(gcd,map(abs,row));rr.append([x//g for x in row])
        mat=DomainMatrix([[ZZ(x) for x in r] for r in rr],(2,d+1),ZZ)
        diag,left,right=smith_normal_decomp(mat)
        assert diag==left*mat*right
        rank=sum(diag.to_Matrix()[i,i]!=0 for i in range(2))
        kernel=right.to_Matrix()[:,rank:].T
        assert sp.Matrix(rr)*kernel.T==sp.zeros(2,d+1-rank)
        reduced=kernel.lll()
        base=[0]*(2*m)+[(-1)**j*comb(2*m,j) for j in range(2*m+1)]
        best=None;zero=0
        for row in reduced.tolist():
          h=[int(x) for x in row]
          assert all(sum(F(x)*mu[j][k] for j,x in enumerate(h))==0 for k in range(2))
          beta,alpha=pi_period(mul(base,h),c)
          if not alpha:zero+=1;continue
          pq=-4*beta/alpha;q=pq.denominator
          up=max(abs(lo-pq),abs(hi-pq))*q*q
          v=dict(m=m,K=K,d=d,H=[str(x) for x in h],p=str(pq.numerator),q=str(q),
                 Q_upper=str(up),displayed_Q=float(up),alpha=str(alpha),zero_rows=zero)
          if best is None or up<F(best['Q_upper']):best=v
        rows.append(best or dict(m=m,K=K,d=d,zero_rows=zero,status='all rows zero period'))
    out=dict(status='completed bounded saturated controls',samples=rows,
             scope='No quality theorem inferred; fixed generation contains no pi value and exact intervals are used afterwards.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{k:x.get(k) for k in ['m','K','d','status','displayed_Q','p','q','zero_rows']} for x in rows],indent=2))

if __name__=='__main__':run()
