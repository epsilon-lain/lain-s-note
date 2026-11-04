"""An exact Gaussian-orbit control; no pi value enters its linear system."""
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as sp

ROOT=Path(__file__).resolve().parent

def main():
    data=[];A,B,w=1,0,0
    for k in range(1,8):
        prev=A;A,B=2*A-B,A+2*B
        if prev*A<0:w+=1
        if k not in [1,6,7]:continue
        h,z=F(A*B,A*A+B*B),F(B*B,A*A+B*B)
        b=F(1);fs=[]
        for j in range(4):
            fs.append(h*b*z**j)
            b*=F(2*(j+1),2*j+3)
        data.append({'k':k,'A':A,'B':B,'w':w,'coefficients':fs})
    rows=[];rhs=[]
    for n in range(1,4):
        rows.append([sp.Rational(row['coefficients'][n-ell]) for row in data for ell in range(2)])
        rhs.append(0)
    for key,target in [('k',0),('w',-1),(None,0)]:
        rows.append([row[key] if key else 1 for row in data for ell in range(2)])
        rhs.append(target)
    v=list(map(F,sp.Matrix(rows).inv()*sp.Matrix(rhs)))
    value=sum((v[2*j]*row['coefficients'][0] for j,row in enumerate(data)),F())
    assert value==F(23040475932864512,7333971781828125)
    factors={3:1,5:6,7:1,29833:1,749209:1}
    product=1
    for p,e in factors.items():
        assert sp.isprime(p);product*=p**e
    assert product==value.denominator
    for row in data:row['coefficients']=list(map(str,row['coefficients']))
    saved=json.loads((ROOT/'results.json').read_text())
    lo,hi=map(F,saved['pi_interval'])
    qlo=min(abs(lo-value),abs(hi-value))*value.denominator**2
    qhi=max(abs(lo-value),abs(hi-value))*value.denominator**2
    assert 10**26<qlo<qhi<10**27
    out={'angles':data,'polynomial_weights':[list(map(str,v[2*j:2*j+2])) for j in range(3)],
         'candidate':str(value),'denominator_factors':factors,
         'generation_reads_pi':False,'quality_between':['1e26','1e27'],
         'scope':'A failed bounded control; common norm support is not a reduced-denominator theorem.'}
    (ROOT/'orbit_control_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'candidate':str(value),'denominator_factors':factors,'status':'PASS'},indent=2))

if __name__=='__main__':main()
