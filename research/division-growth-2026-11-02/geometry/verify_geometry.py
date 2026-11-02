"""Exact geometry of a quotient that removes a shared growth factor.
SymPy handles polynomial/rational identities; Fraction bounds atan endpoint error.
This is not an irrationality-measure computation or a new pi approximation.
"""
from pathlib import Path
from fractions import Fraction as Q
import json
import sympy as s

t=s.symbols('t',real=True)

def angular(X,Y):
    return (X*s.diff(Y,t)-Y*s.diff(X,t))/(X**2+Y**2)

def atan_interval(x,N=50):
    x=Q(x)
    if x<0:
        lo,hi=atan_interval(-x,N);return -hi,-lo
    assert 0<=x<1
    a=sum(((-1)**j*x**(2*j+1)/Q(2*j+1) for j in range(N)),Q(0))
    b=a+(-1)**N*x**(2*N+1)/Q(2*N+1)
    return min(a,b),max(a,b)

def run():
    X,Y,a=[s.Function(name)(t) for name in ('X','Y','a')]
    assert s.factor(angular(a*X,a*Y)-angular(X,Y))==0
    x=(1-t*t)/(1+t*t);y=2*t/(1+t*t)
    omega=2/(1+t*t)
    assert s.factor(x*x+y*y-1)==0
    assert s.factor(angular(x,y)-omega)==0
    assert (x.subs(t,0),y.subs(t,0),x.subs(t,1),y.subs(t,1))==(1,0,0,1)
    # Common scale cancels from angular speed but not normalized Euclidean speed.
    V=s.Matrix([a*x,a*y])
    euclidean_squared=s.factor(V.diff(t).dot(V.diff(t))/a**2)
    assert s.factor(euclidean_squared-(s.diff(a,t)/a)**2-omega**2)==0
    u=s.Function('u')(t);v=s.Function('v')(t)
    Rx=(1+u)*x-v*y;Ry=(1+u)*y+v*x
    phase_derivative=((1+u)*s.diff(v,t)-v*s.diff(u,t))/((1+u)**2+v**2)
    assert s.factor(angular(Rx,Ry)-omega-phase_derivative)==0
    assert s.simplify(s.diff(s.atan(v/(1+u)),t)-phase_derivative)==0
    # A noncircular rational curve, whose integrated angular quotient is unchanged.
    A=1+t/2;U=1+t/3;V=t*(1-t)/5
    X=A*(U*x-V*y);Y=A*(U*y+V*x)
    exact_extra=15*(3-6*t-t*t)/(9*t**4-18*t**3+34*t*t+150*t+225)
    assert s.factor(angular(X,Y)-omega-exact_extra)==0
    assert s.factor(exact_extra-s.diff(s.atan(V/U),t))==0
    assert V.subs(t,0)==V.subs(t,1)==0
    assert (X.subs(t,0),Y.subs(t,0),X.subs(t,1),Y.subs(t,1))==(1,0,0,2)
    # Four points fail the exact concyclicity determinant: this is not any circle.
    rows=[]
    for q in (s.Rational(0),s.Rational(1,3),s.Rational(2,3),s.Rational(1)):
        xx=X.subs(t,q);yy=Y.subs(t,q)
        rows.append([xx*xx+yy*yy,xx,yy,1])
    circle_det=s.factor(s.Matrix(rows).det())
    assert circle_det==-s.Rational(25330898,159924375)
    rational_cases=0
    for ac in (s.Rational(1,4),s.Rational(1,2),s.Rational(1)):
        for uc in (s.Rational(-1,2),s.Rational(0),s.Rational(1,3)):
            for vc in (s.Rational(-1,3),s.Rational(1,5),s.Rational(2,5)):
                A=1+ac*t;U=1+uc*t;V=vc*t*(1-t)
                X=A*(U*x-V*y);Y=A*(U*y+V*x)
                assert s.factor(angular(X,Y)-omega-s.diff(s.atan(V/U),t))==0
                assert U.subs(t,0)>0 and U.subs(t,1)>0
                rational_cases+=1
    # First-order endpoint phase error cancels; verify the rigorous second-order bound.
    endpoint_cases=0
    for pa in (Q(-1),Q(0),Q(1)):
        for pb in (Q(-1),Q(0),Q(1)):
            for q0 in (Q(-1,2),Q(1,3),Q(1)):
                for eps in (Q(1,20),Q(1,10),Q(1,4)):
                    eta=max(abs(eps*pa),abs(eps*pb));assert eta<1
                    if pa==pb:actual=(Q(0),Q(0))
                    else:
                        aa=atan_interval(eps*q0/(1+eps*pa))
                        bb=atan_interval(eps*q0/(1+eps*pb))
                        actual=(bb[0]-aa[1],bb[1]-aa[0])
                    bound=eps**2*abs(q0)*abs(pb-pa)/(1-eta)**2
                    assert max(abs(z) for z in actual)<=bound
                    endpoint_cases+=1
    # Higher-dimensional angular bivector removes the same scalar factor.
    z=s.Matrix([2*t,2*t*t,1-t*t-t**4])/(1+t*t+t**4)
    assert s.factor(z.dot(z)-1)==0
    a=1+t/2
    az=a*z;sum_w_squared=0
    for i in range(3):
        for j in range(i+1,3):
            W=z[i]*s.diff(z[j],t)-z[j]*s.diff(z[i],t)
            Aij=(az[i]*s.diff(az[j],t)-az[j]*s.diff(az[i],t))/az.dot(az)
            assert s.factor(Aij-W)==0
            sum_w_squared+=W*W
    assert s.factor(sum_w_squared-z.diff(t).dot(z.diff(t)))==0
    # At a common zero the observed quotient is undefined, even if algebra cancels.
    rawX=t*x;rawY=t*y
    assert (rawX**2+rawY**2).subs(t,0)==0
    assert s.factor(angular(rawX,rawY)-omega)==0
    report={'status':'PASS','symbolic_scalar_gauge_identity':True,'rational_circle_identity':True,
       'noncircular_rational_deformation_cases':rational_cases,
       'second_order_endpoint_error_interval_checks':endpoint_cases,
       'higher_dimensional_bivector_checks':3,'common_zero_domain_condition_checked':True,
       'noncircular_four_point_determinant':str(circle_det),
       'explicit_nontrivial_extra_integrand':str(exact_extra),
       'scope':'Exact quotient and endpoint-phase identities; no new rational approximation to pi or theorem about its bad approximability.',
       'dependencies':{'python':'3','sympy':s.__version__}}
    Path(__file__).with_name('checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':run()
