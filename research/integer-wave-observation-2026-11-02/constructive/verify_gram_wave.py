"""Independent exact checks of the stored Gram construction and fixed-fiber projection."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import json,sys
import sympy as s
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'support'))
from integral_certificate import certificate
t=s.Symbol('t'); checked=0

def check_row(m,row,fiber=False):
    global checked
    c=4*m+1; R=list(map(int,row['R'])); alpha=int(row['alpha'])
    N=[0]*(c+1)
    for j in range(c+1):
        if j<c: N[j]+=(2*j+1)*R[j]+alpha*comb(c-1,j)
        if j: N[j]+=(2*j-2*c+1)*R[j-1]
    pol=sum(z*t**j for j,z in enumerate(N))
    H,rem=s.div(pol,t**m*(1-t)**(2*m)); assert rem==0
    h=[s.expand(H).coeff(t,j) for j in range(m+2)]
    if not fiber: assert list(map(int,h))==row['H']
    energy=sum((h[i]*h[j]/s.Integer(2*i+2*j+1) for i in range(m+2) for j in range(m+2)),s.Rational(0))
    energy=F(int(energy.p),int(energy.q))
    if fiber: assert energy/F(alpha*alpha)==F(row['optimal_normalized_energy'])
    else: assert energy==F(row['energy'])
    rat=F(-sum(R),2**(c-3)*alpha)
    assert (str(rat.numerator),str(rat.denominator))==(row['p'],row['q'])
    prefactor=F(4*rat.denominator**2,abs(alpha))*F(2,27)**m
    assert prefactor**2*energy==F(row['gram_Q_bound_squared'])
    cert=certificate(c,R,alpha,K=160)
    assert F(cert['quality_upper'])**2<F(row['gram_Q_bound_squared'])
    checked+=1

z=json.loads((HERE/'gram_wave_checks.json').read_text())
for item in z['search']:
    for row in item['top_by_gram_certificate']:check_row(item['m'],row)
f=json.loads((HERE/'fiber_projection_short_checks.json').read_text())
for row in f['results']:check_row(row['m'],row,True)
assert s.expand((1+t)**4-16*t*(1-t)**2-(t*t-6*t+1)**2)==0
print(json.dumps({'status':'PASS','independent_rows_checked':checked,'sharp_kernel_SOS_identity':True,'pi_input':False}))
