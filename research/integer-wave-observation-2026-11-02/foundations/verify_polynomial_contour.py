"""Exact, data-driven contour nonvanishing and polynomial root-count certificates."""
from fractions import Fraction as F
from pathlib import Path
import json


def add(a,b): return a[0]+b[0],a[1]+b[1]
def mul(a,b): return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def sub(a,b): return a[0]-b[0],a[1]-b[1]
def scale(a,t): return a[0]*t,a[1]*t

def dot(a,b): return a[0]*b[0]+a[1]*b[1]
def cross(a,b): return a[0]*b[1]-a[1]*b[0]

def segment_distance_squared(a,b):
    v=sub(b,a)
    vv=dot(v,v)
    if vv==0: return dot(a,a)
    t=max(F(0),min(F(1),-dot(a,v)/vv))
    p=add(a,scale(v,t))
    return dot(p,p)

def polygon_distance_squared(points):
    return min(segment_distance_squared(a,b) for a,b in zip(points,points[1:]))

def winding(points):
    assert points[0]==points[-1]
    result=0
    for a,b in zip(points,points[1:]):
        if cross(a,b)==0 and dot(a,b)<=0:
            raise ValueError('image polygon meets origin')
        c=cross(a,b)
        if a[1]<=0<b[1] and c>0: result+=1
        if b[1]<=0<a[1] and c<0: result-=1
    return result

def square(t):
    t=t%1
    side=int(4*t)
    u=4*t-side
    return [(2*u-1,F(-1)),(F(1),2*u-1),(1-2*u,F(1)),(F(-1),1-2*u)][side]

def evaluate(coeff,z):
    p=(F(0),F(0))
    for a in reversed(coeff): p=add(mul(p,z),(a,F(0)))
    return p

def certificate(coeff, eta=F(0),maxN=4096):
    # Square contour speed 8, |z|<=2; coefficients are real Fractions.
    L=8*sum(j*abs(a)*2**(j-1) for j,a in enumerate(coeff) if j)
    N=8
    while N<=maxN:
        points=[]
        for j in range(N):
            value=evaluate(coeff,square(F(j,N)))
            # Measured error is known only through its rational bound eta.
            noise=eta*F((j*7)%11-5,5)
            points.append((value[0]+noise,value[1]))
        points.append(points[0])
        delta=eta+L/(2*N)
        d2=polygon_distance_squared(points)
        if d2>delta*delta:
            return {'certified':True,'N':N,'eta':str(eta),'L':str(L),'delta':str(delta),
                    'polygon_distance_squared':str(d2),'strict_margin_squared':str(d2-delta*delta),
                    'inside_zero_count':winding(points),'coefficients':[str(a) for a in coeff]}
        if d2==0:
            # This polygon currently fails. Exact zero samples prove true boundary zeros
            # only in the noiseless case; do not infer that from a crossing edge.
            if eta==0 and any(p==(0,0) for p in points):
                return {'certified':False,'reason':'exact polynomial zero on contour',
                        'N':N,'coefficients':[str(a) for a in coeff]}
        N*=2
    return {'certified':False,'reason':'budget exhausted'}

cases=[('two_inside',[F(1,4),F(0),F(1)],2),
       ('zero_inside',[F(4),F(0),F(1)],0),
       ('triple_inside',[F(-1,8),F(3,4),F(-3,2),F(1)],3)]
rows=[]
for name,coeff,expected in cases:
    for eta in [F(0),F(1,1000)]:
        row=certificate(coeff,eta)
        assert row['certified'] and row['inside_zero_count']==expected,(name,row)
        row.update(name=name,expected=expected)
        rows.append(row)
# A boundary zero must not be reported as an interior count.
boundary=certificate([F(-1),F(1)])
assert not boundary['certified'] and boundary['reason']=='exact polynomial zero on contour'
# Distance projections: closest point can lie in the interior of an edge.
assert segment_distance_squared((F(-2),F(1)),(F(3),F(1)))==1
assert segment_distance_squared((F(2),F(1)),(F(3),F(1)))==5
out={'status':'PASS','arithmetic':'fractions.Fraction only','data_driven_root_certificates':len(rows),
     'requires_true_boundary_margin':False,'boundary_zero_rejected':True,'examples':rows,'boundary_case':boundary}
Path(__file__).with_name('polynomial_contour_checks.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k not in ['examples','boundary_case']}))
print(json.dumps(rows,indent=2))
