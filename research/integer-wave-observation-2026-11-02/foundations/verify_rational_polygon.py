"""Exact Fraction polygon winding and rigorous margins for rational circle paths."""
from fractions import Fraction as F
from pathlib import Path
import json


def cross(a,b): return a[0]*b[1]-a[1]*b[0]
def dot(a,b): return a[0]*b[0]+a[1]*b[1]

def winding(vertices):
    """Signed ray-crossing count; rejects any segment containing the origin."""
    assert vertices[0]==vertices[-1]
    result=0
    for a,b in zip(vertices,vertices[1:]):
        if cross(a,b)==0 and dot(a,b)<=0:
            raise ValueError('boundary meets origin')
        c=cross(a,b)
        if a[1]<=0<b[1] and c>0: result+=1
        if b[1]<=0<a[1] and c<0: result-=1
    return result

def circle(t):
    """Unit circle, four rationally parametrized quarter arcs, period 1."""
    t=t%1
    q=int(4*t)
    u=4*t-q
    a,b=(1-u*u)/(1+u*u),2*u/(1+u*u)
    return [(a,b),(-b,a),(-a,-b),(b,-a)][q]

rows=[]
for k in [-9,-3,-1,0,1,2,7,13]:
    for shift in [F(0),F(1,3),F(3,4),F(5,4),F(2)]:
        for eta in [F(0),F(1,100)]:
            r=F(1)
            m=abs(shift-r)
            L=8*r*abs(k)
            N=max(16,16*abs(k))
            while eta+L/(2*N)>=m: N*=2
            points=[]
            for j in range(N):
                a,b=circle(F(k*j,N))
                # Measured x error varies deterministically; norm<=eta exactly.
                noise=eta*F((j*7)%11-5,5)
                points.append((r*a+shift+noise,r*b))
            points.append(points[0])
            observed=winding(points)
            expected=k if shift<r else 0
            assert observed==expected
            assert eta+L/(2*N)<m
            rows.append(dict(k=k,shift=str(shift),eta=str(eta),N=N,
                             m=str(m),L=str(L),max_error=str(eta+L/(2*N)),
                             winding=observed,expected=expected))
# Boundary rejection and orientation/addition tests.
try:
    winding([(F(-1),F(0)),(F(1),F(0)),(F(-1),F(0))])
except ValueError: pass
else: raise AssertionError('origin crossing accepted')
square=[(F(1),F(1)),(F(-1),F(1)),(F(-1),F(-1)),(F(1),F(-1)),(F(1),F(1))]
assert winding(square)==1 and winding(square[::-1])==-1
assert winding(square[:-1]+square)==2
out={'status':'PASS','arithmetic':'fractions.Fraction only','margin_certified_examples':len(rows),'boundary_rejection':True,'orientation_and_additivity':True,'examples':rows}
Path(__file__).with_name('rational_polygon_checks.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='examples'}))
