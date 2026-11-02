"""Exact rational checks for the modular subsystem occurring in the S6 period data."""
from fractions import Fraction as F
from math import comb,factorial,lcm
from pathlib import Path
import json

def mul(a,b,N):
 out=[F(0)]*(N+1)
 for i,x in enumerate(a[:N+1]):
  if not x:continue
  for j,y in enumerate(b[:N+1-i]):
   if y:out[i+j]+=x*y
 return out

def power(a,k,N):
 out=[F(1)]+[F(0)]*N
 for _ in range(k):out=mul(out,a,N)
 return out

def binomial_power(a,r,N):
 assert a[0]==1
 v=list(a[:N+1])+[F(0)]*max(0,N+1-len(a));v[0]=0
 term=[F(1)]+[F(0)]*N;out=term[:];coef=F(1)
 for k in range(1,N+1):
  term=mul(term,v,N);coef*=F(r-k+1,k)
  for j,t in enumerate(term):out[j]+=coef*t
 return out

def compose(a,z,N):
 term=[F(1)]+[F(0)]*N;out=[F(0)]*(N+1)
 for coef in a[:N+1]:
  for j,t in enumerate(term):out[j]+=coef*t
  term=mul(term,z,N)
 return out

def main():
 N=40;a,b=F(1,12),F(5,12)
 f=[F(1)];h=[F(0)];S=F(0);M=[1];L=1
 for n in range(1,N+1):
  k=n-1;f.append(f[-1]*(k+a)*(k+b)/(n*n))
  S+=1/(k+a)+1/(k+b)-F(2,n);h.append(f[-1]*S)
  assert (1728**n*f[-1]).denominator==1
  L=lcm(L,*range(12*k+1,12*n+1))
  assert (L*1728**n*h[-1]).denominator==1
  M.append(24*(2*k+1)*(6*k+1)*(6*k+5)*M[-1]//n**3)
  assert M[-1]==factorial(6*n)//(factorial(3*n)*factorial(n)**3)
  assert h[n]*n*n-(n-1+a)*(n-1+b)*h[n-1] == -2*n*f[n]+(2*(n-1)+a+b)*f[n-1]
 G=mul(f,f,N)
 assert G==[F(m,1728**n) for n,m in enumerate(M)]
 # z W(f,g)=f^2+z(f h'-f' h)=(1-z)^(-1/2)
 cross=mul(f,[F((n+1))*h[n+1] for n in range(N)],N-1)
 cross2=mul([F((n+1))*f[n+1] for n in range(N)],h,N-1)
 W=[G[0]]+[G[n]+cross[n-1]-cross2[n-1] for n in range(1,N+1)]
 assert W==[F(comb(2*n,n),4**n) for n in range(N+1)]
 # Explicit rational pullback into the lambda-based Ramanujan series.
 P=30;z=[F(0),F(0)]+[F(27*comb(k+2,2),64*4**k) for k in range(P-1)]
 H=[F(comb(2*n,n)**3,64**n) for n in range(P+1)]
 sq=binomial_power([F(1),F(-1,4)],F(1,2),P)
 assert compose(G,z,P)==mul(sq,H,P)
 # Independent q-series checks of E4^(1/4), j=E4^3/Delta and cusp normalization.
 Q=12;T=Q+1
 E4=[F(1)]+[F(240*sum(d**3 for d in range(1,n+1) if n%d==0)) for n in range(1,T+1)]
 Delta=[F(0),F(1)]+[F(0)]*(T-1)
 for n in range(1,T+1):
  factor=[F(0)]*(T+1)
  for j in range(0,min(24,T//n)+1):factor[n*j]=F((-1)**j*comb(24,j))
  Delta=mul(Delta,factor,T)
 zq=[1728*t for t in mul(Delta,binomial_power(E4,F(-3),T),T)]
 fz=compose(f,zq,Q)
 assert fz==binomial_power(E4,F(1,4),Q)
 unit=[zq[n+1]/1728 for n in range(Q+1)]
 v=unit[:];v[0]-=1;term=[F(1)]+[F(0)]*Q;logunit=[F(0)]*(Q+1)
 for k in range(1,Q+1):
  term=mul(term,v,Q)
  for j,x in enumerate(term):logunit[j]+=F((-1)**(k+1),k)*x
 ratio=mul(compose(h,zq,Q),binomial_power(fz,F(-1),Q),Q)
 assert all(x+y==0 for x,y in zip(logunit,ratio))
 out={'status':'PASS','Frobenius_ODE_and_integrality_coefficients':N,'Wronskian_and_Clausen_coefficients':N,'rational_pullback_coefficients':P,'independent_modular_q_series_coefficients':Q,'cusp_log_normalization':'2*pi*i*tau=g/f-log(1728)','first_f_coefficients':[str(x) for x in f[:6]],'first_1728_scaled_f_coefficients':[str(f[n]*1728**n) for n in range(6)],'first_h_coefficients':[str(x) for x in h[:4]],'first_square_integer_coefficients':[str(x) for x in M[:6]],'scope':'Exact finite audits supporting the all-parameter derivations. No new pi approximation record.'}
 Path(__file__).with_name('modular_subsystem_checks.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
