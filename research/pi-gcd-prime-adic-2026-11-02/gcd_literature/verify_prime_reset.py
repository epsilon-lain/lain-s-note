import json,math
from pathlib import Path
N=2000
A,B=[1,1],[0,1]
for n in range(1,N):
 A.append((2*n+1)*A[-1]+n*n*A[-2]);B.append((2*n+1)*B[-1]+n*n*B[-2])
def primes(n):
 s=bytearray(b'\1')*(n+1);s[:2]=b'\0\0'
 for p in range(2,math.isqrt(n)+1):
  if s[p]:s[p*p:n+1:p]=b'\0'*(((n-p*p)//p)+1)
 return [p for p in range(2,n+1) if s[p]]
def v(n,p):
 ans=0
 while n%p==0:ans+=1;n//=p
 return ans
def vf(n,p):
 ans=0
 while n:n//=p;ans+=n
 return ans
checks=0;block_examples={}
for p in primes(500):
 if p==2:continue
 chi=(-1)**((p+1)//2)
 assert A[p-1]%p==chi%p
 assert B[p]%p==chi%p
 assert all(A[n]%p==0 for n in range(p,N+1))
 for r in range(p+1):
  assert B[p+r]%p==chi*A[r]%p; checks+=1
 for n in range(1,min(3*p,N)+1):
  expected=False if n<p else (A[n-p]%p==0 if n<2*p else True)
  assert (math.gcd(A[n],4*B[n])%p==0)==expected;checks+=1
 if p<50:block_examples[p]=[r for r in range(p) if A[r]%p==0]
for k in range((N-3)//4+1):
 assert v(A[4*k],2)==2*k
 assert v(A[4*k+1],2)==2*k
 assert v(A[4*k+2],2)==2*k+2
 assert v(A[4*k+3],2)>=2*k+3
excess=[]
f=1
for n in range(1,N+1):
 f*=n;g=math.gcd(A[n],4*B[n]);oddg=g//(g&-g);extra=oddg//math.gcd(oddg,f)
 if extra>1:
  fac=[]
  for p in primes(n):
   if extra%p==0:fac.append([p,v(g,p),vf(n,p)])
  excess.append([n,fac])
result={'max_n':N,'prime_upper':500,'prime_reset_checks':checks,'first_block_zero_residues':block_examples,'first_higher_valuation_excess':excess[:8],'higher_valuation_excess_count':len(excess),'partial_two_adic_theorem_checked':True,'note':'Finite checks audit proved statements; stronger full 2-adic formula remains unproved.'}
Path(__file__).with_name('prime_reset_checks.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
