"""Bounded constructive family a=m,b=2m,c=4m+1; no pi in LLL features."""
from fractions import Fraction as F
from pathlib import Path
from math import lcm,comb
import json,signal,time
import sympy as s
from antiderivative_search_saturated import data,samples,atan

def timeout_handler(signum,frame):raise TimeoutError('30-second per-m budget exceeded')
signal.signal(signal.SIGALRM,timeout_handler)
start=time.monotonic();l,h=atan(5,300);u,v=atan(239,300);lo,hi=16*l-4*v,16*h-4*u
out=[];stopped=None
for m in range(1,7):
 tick=time.monotonic();signal.alarm(30);a=m;b=2*m;c=4*m+1;rows=[];seen=set();zero=0
 try:
  M,B=data(c,a,b);S=samples(c)
  for budget in [0,8*m,12*m,16*m,20*m]:
   cost=F(1,2**budget)
   feats=S+[[cost if j==c else F(0) for j in range(c+1)],[cost if j<c else F(0) for j in range(c+1)]]
   D=lcm(*(v.denominator for row in feats for v in row))
   feature=s.Matrix([[int(v*D) for v in row] for row in feats]).T
   extended=(B*feature).row_join(B)
   reduced=extended.lll(delta=F(3,4))
   for row in reduced.tolist():
    R=list(map(int,row[-(c+1):]));alpha=R[-1];r1=sum(R[:-1])
    if not alpha:zero+=1;continue
    rat=F(-4*r1,alpha*2**(c-1));p,q=rat.numerator,rat.denominator
    if p<=0 or (p,q) in seen:continue
    seen.add((p,q))
    if rat<lo:L,H=(lo-rat)*q*q,(hi-rat)*q*q
    elif rat>hi:L,H=(rat-hi)*q*q,(rat-lo)*q*q
    else:raise AssertionError('narrower pi enclosure needed')
    t=s.symbols('t');N=list(map(int,M*s.Matrix(R)));P=sum(z*t**j for j,z in enumerate(N));Hp,rem=s.div(P,t**a*(1-t)**b);assert rem==0
    hc=[int(s.expand(Hp).coeff(t,j)) for j in range(c-a-b+1)];hnorm=sum(abs(z) for z in hc)
    shape=F(4*q*q*hnorm,abs(alpha))*F(2,27)**m
    scale=10**25
    rows.append({'cost_weight_power_of_two':-budget,'alpha':str(alpha),'R_coefficients':[str(z) for z in R[:-1]],'R_at_one':str(r1),'N_coefficients':[str(z) for z in N],'H_coefficients':[str(z) for z in hc],'H_l1':str(hnorm),'p':str(p),'q':str(q),'quality_display':float(H),'quality_rational_bounds':[str(F(L.numerator*scale//L.denominator,scale)),str(F(-((-H.numerator*scale)//H.denominator),scale))],'shape_quality_upper':str(shape),'shape_quality_upper_display':float(shape)})
  best=min(rows,key=lambda z:F(z['quality_rational_bounds'][1])) if rows else None
  out.append({'m':m,'a':a,'b':b,'c':c,'H_degree_bound':m+1,'elapsed_seconds':time.monotonic()-tick,'candidate_count':len(rows),'zero_alpha_rows':zero,'best_observed':best,'candidates':rows,'scope':'Finite candidates of a predeclared family; quality used only after LLL selection.'})
 except TimeoutError as e:
  stopped={'m':m,'reason':str(e),'elapsed_seconds':time.monotonic()-tick};break
 finally:signal.alarm(0)
result={'status':'PASS' if stopped is None else 'BOUNDED_PARTIAL','family':'a=m,b=2m,c=4m+1,H_degree<=m+1','generation':'Full integer endpoint kernel and rational sampled LLL features; no pi or target fraction in candidate generation.','quality_certification':'Machin 300-term rational interval after selection; root may separately certify directly from the integer integrand.','completed_m':[r['m'] for r in out],'per_m':out,'stopped':stopped,'elapsed_seconds':time.monotonic()-start,'scope':'Constructive finite family experiment, not new rational approximation records or an infinite quality-zero theorem.'}
Path(__file__).with_name('new_family_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='per_m'},indent=2))
for r in out:print(json.dumps({'m':r['m'],'candidate_count':r['candidate_count'],'best':{k:v for k,v in (r['best_observed'] or {}).items() if k in ['p','q','quality_display','shape_quality_upper_display','alpha']}}))
