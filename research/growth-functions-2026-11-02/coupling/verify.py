"""Exact rational controls of the separately proved locking and smoothing certificates."""
from fractions import Fraction as Q
from math import gcd,lcm
from hashlib import sha256
from pathlib import Path
import json
from exact_circle import F,cycles,iterate_pieces
HERE=Path(__file__).resolve().parent

def iterate(x,omega,kappa,q):
    for _ in range(q):x=F(x,omega,kappa)
    return x

def interval_certificate(omega,kappa,p,q,left,right):
    assert gcd(p,q)==1 and left<right
    # F is strictly increasing: endpoint evaluation encloses the exact interval image.
    lo=iterate(left,omega,kappa,q)-p
    hi=iterate(right,omega,kappa,q)-p
    assert left<lo<=hi<right
    margin=min(lo-left,right-hi)
    return {'p':p,'q':q,'interval':[str(left),str(right)],
            'iterate_minus_p_image':[str(lo),str(hi)],'strict_margin':str(margin)}

# General q-step splitter uses rational affine pieces only.
piece_checks=0;point_checks=0;max_pieces=0
for omega in [Q(0),Q(1,7),Q(17,50),Q(51,100),Q(99,100)]:
    for kappa in [Q(1,4),Q(1,2),Q(3,4)]:
        for q in range(1,17):
            pieces=iterate_pieces(omega,kappa,q)
            assert len(pieces)<=2*q+1
            max_pieces=max(max_pieces,len(pieces))
            assert pieces[0][0]==0 and pieces[-1][1]==1
            for (lo,hi,A,B) in pieces:
                assert A>0 and lo<hi
                for x in [lo,(lo+hi)/2,hi]:
                    assert iterate(x,omega,kappa,q)==A*x+B
                    point_checks+=1
            for old,new in zip(pieces,pieces[1:]):
                assert old[1]==new[0]
                x=old[1]
                assert old[2]*x+old[3]==new[2]*x+new[3]
            assert iterate(Q(1),omega,kappa,q)==iterate(Q(0),omega,kappa,q)+1
            piece_checks+=1

# The full closed locking band of the explicit orbit formula is certified;
# this code does not assert it is the maximal tongue.
band_checks=0
for kappa in [Q(1,4),Q(1,2),Q(3,4)]:
    width=kappa*kappa/(4*(2+kappa))
    for j in range(-50,51):
        delta=width*Q(j,50);omega=Q(1,2)+delta
        x0=(2+kappa)*delta/(kappa*kappa)
        x1=Q(1,2)+(2-kappa)*delta/(kappa*kappa)
        assert -Q(1,4)<=x0<=Q(1,4)
        assert Q(1,4)<x1<Q(3,4)
        assert F(x0,omega,kappa)==x1
        assert F(x1,omega,kappa)==x0+1
        band_checks+=1

examples=[]
for omega,p,q,x in [(Q(51,100),1,2,Q(1,10)),(Q(17,50),1,3,Q(28,125))]:
    kappa=Q(1,2)
    discovered=cycles(omega,kappa,p,q)
    assert any(row[2]==x for row in discovered)
    orbit=[];z=x
    for j in range(q):orbit.append(str(z));z=F(z,omega,kappa)
    assert z==x+p
    for j in range(1,q):assert (iterate(x,omega,kappa,j)-x).denominator!=1
    cert=interval_certificate(omega,kappa,p,q,x-Q(1,100),x+Q(1,100))
    cert.update({'omega':str(omega),'kappa':str(kappa),'orbit':orbit,
                 'last_lift':str(z),'piece_count':len(iterate_pieces(omega,kappa,q)),
                 'minimum_circle_period':q,'candidate_affine_roots':len(discovered)})
    examples.append(cert)

# Fejer smoothing error: H_1024<8, pi^2>9.
N=1024
H=sum((Q(1,n) for n in range(1,N+1)),Q(0))
assert H<8
smooth_error=Q(1,9)*(Q(8,N+1)+Q(1,N)) # 1/2 times the psi error
assert smooth_error<Q(1,1000)
L=Q(3,2)
for cert in examples:
    q=cert['q'];amplification=sum((L**j for j in range(q)),Q(0))
    assert amplification*smooth_error<Q(cert['strict_margin'])
    cert['smooth_iterate_perturbation_bound']=str(amplification*smooth_error)
    cert['smooth_lock_persists']=True

# Charge the common denominator of the rational part of all Fourier coefficients.
coefficients=[(n,Q((1 if n%4==1 else -1)*(N+1-n),(N+1)*n*n))
              for n in range(1,N+1,2)]
D=lcm(*(r.denominator for n,r in coefficients))
ints=[(n,int(r*D)) for n,r in coefficients]
assert all(r*D==b for (n,r),(nn,b) in zip(coefficients,ints))
assert gcd(D,*[abs(b) for n,b in ints])==1
payload='\n'.join(f'{n}:{b}' for n,b in ints).encode()

# Analytic invariant-tube bound for the separate anchored 2D coupling.
r=Q(1,10);eps=Q(1,10)
contraction=Q(2,3)*Q(22,7)**2*r*r+2*eps
assert contraction==Q(977,3675) and contraction<1

out={'status':'PASS','scope':'Exact implementation controls for separately proved general theorems; no pi approximation theorem.',
     'piece_partition_cases':piece_checks,'exact_point_evaluation_checks':point_checks,
     'maximum_piece_count':max_pieces,'closed_band_formula_checks':band_checks,
     'locking_examples':examples,
     'fejer':{'N':N,'nonzero_sine_terms':len(coefficients),'F_sup_error_bound':str(smooth_error),
               'rational_common_denominator_bits':D.bit_length(),
               'max_integer_coefficient_bits':max(abs(b).bit_length() for n,b in ints),
               'integer_coefficients_sha256':sha256(payload).hexdigest(),
               'common_analytic_prefactor':'2/pi^2; retained, not treated as a rational coefficient'},
     'anchored_2D_coupling':{'radius':str(r),'epsilon':str(eps),'contraction_bound':str(contraction)},
     'sources':['https://www.math.stonybrook.edu/~jack/DYNOTES/dn14.pdf']}
(HERE/'checks.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(out,ensure_ascii=False,indent=2))
