"""Exact early-stop sieve for the fixed Legendre arctangent family.

Main sieve constructs U_n and short B prefixes, never full B_n. Optional
--audit independently builds the full integer A,B sequences to check every
bound. All decisions use rational arithmetic, not a floating pi value.
"""
from fractions import Fraction as F
from math import gcd, isqrt, lcm
from pathlib import Path
import argparse
import json
from time import perf_counter


class GMPGCD:
    """Optional exact gcd accelerator using the system GMP shared library.

    The rest of the sieve remains native Python integer arithmetic. No
    third-party Python package is required; omit --gmp for stdlib only.
    """
    def __init__(self):
        import ctypes as ct
        from ctypes.util import find_library
        import atexit
        import platform
        if platform.system()!='Linux' or ct.sizeof(ct.c_void_p)!=8 or ct.sizeof(ct.c_ulong)!=8:
            raise RuntimeError('--gmp currently supports the conventional 64-bit Linux GMP ABI only')
        class Z(ct.Structure):
            _fields_=[('alloc',ct.c_int),('size',ct.c_int),
                      ('limbs',ct.POINTER(ct.c_ulong))]
        self.ct=ct
        self.lib=ct.CDLL(find_library('gmp') or 'libgmp.so.10')
        if ct.c_int.in_dll(self.lib,'__gmp_bits_per_limb').value!=64:
            raise RuntimeError('--gmp requires 64-bit GMP limbs')
        zp=ct.POINTER(Z)
        self.init=getattr(self.lib,'__gmpz_init'); self.init.argtypes=[zp]
        self.clear=getattr(self.lib,'__gmpz_clear'); self.clear.argtypes=[zp]
        self.imp=getattr(self.lib,'__gmpz_import')
        self.imp.argtypes=[zp,ct.c_size_t,ct.c_int,ct.c_size_t,
                           ct.c_int,ct.c_size_t,ct.c_void_p]
        self.exp=getattr(self.lib,'__gmpz_export')
        self.exp.argtypes=[ct.c_void_p,ct.POINTER(ct.c_size_t),ct.c_int,
                           ct.c_size_t,ct.c_int,ct.c_size_t,zp]
        self.op=getattr(self.lib,'__gmpz_gcd'); self.op.argtypes=[zp,zp,zp]
        self.z=[Z(),Z(),Z()]
        for z in self.z:
            self.init(ct.byref(z))
        atexit.register(self.close)

    def close(self):
        for z in self.z:
            self.clear(self.ct.byref(z))
        self.z=[]

    def __call__(self,a,b):
        from math import gcd as native_gcd
        if min(a.bit_length(),b.bit_length())<8192:
            return native_gcd(a,b)
        ct=self.ct
        for z,value in zip(self.z,(abs(a),abs(b))):
            data=value.to_bytes(max(1,(value.bit_length()+7)//8),'little')
            self.imp(ct.byref(z),len(data),-1,1,0,0,data)
        self.op(ct.byref(self.z[2]),ct.byref(self.z[0]),ct.byref(self.z[1]))
        out=ct.create_string_buffer(max(1,(min(a.bit_length(),b.bit_length())+7)//8))
        count=ct.c_size_t()
        self.exp(out,ct.byref(count),-1,1,0,0,ct.byref(self.z[2]))
        return int.from_bytes(out.raw[:count.value],'little')


def v2(x):
    return (x & -x).bit_length()-1


def odd(x):
    return x >> v2(x)


def primes_to(N):
    flags = bytearray(b'\1')*(N+1)
    flags[:2] = b'\0\0'
    for p in range(2, isqrt(N)+1):
        if flags[p]:
            flags[p*p:N+1:p] = b'\0'*((N-p*p)//p+1)
    return [p for p in range(3,N+1,2) if flags[p]]


def quality_floor(n, d, G, D):
    num, den = 2*D*D, (2*n+1)*G*G
    if d >= 0:
        den <<= 2*d
    else:
        num <<= -2*d
    return F(num, den)


def prune_overvaluation(n, C, G, primes):
    """Remove whole prime powers that cannot be positive surplus.

    G is the old gcd(C,G0). If its p-exponent is less than the universal
    cap k, then m=v_p(C) is already that exponent and is safe. At equality
    test C modulo p^(k+1); a zero means m>k and removes the full factor.
    Only the small cap G is factored, not the full integer C.
    """
    remaining=G
    factors=[]
    for p in primes:
        if p*p>remaining:
            break
        if remaining%p==0:
            value=1
            while remaining%p==0:
                remaining//=p
                value*=p
            factors.append((p,value))
    if remaining>1:
        factors.append((remaining,remaining))
    removed=[]
    result=G
    for p,value in factors:
        power=p
        k=0
        while power<=n-1:
            k+=1
            power*=p
        u=n
        while u%p==0:
            u//=p
            k-=1
        assert k>=1
        if value==p**k and C%(p**(k+1))==0:
            result//=value
            removed.append({'p':p,'cap':k,'removed_factor':str(value)})
    return result,removed


def run(N, target=F(1), audit=False, output=None, use_gmp=False):
    global gcd
    from math import gcd as native_gcd
    gcd=GMPGCD() if use_gmp else native_gcd
    started=perf_counter()
    primes = primes_to(N)
    Bsmall = [0,1]
    Uprev, U = 1,2
    L = 1
    Aprev=A=1
    Bprev,B=0,1
    fac=1
    survivors=[]
    traces=[]
    max_tests=max_prefix=0
    certified=baseline=0
    old_baseline=pruned_indices=removed_factor_count=total_tests=0
    for n in range(1,N+1):
        if n > 1:
            k=n-1
            value=2*(2*k+1)*U+4*k*Uprev
            assert value%n==0
            Uprev,U=U,value//n
            L=lcm(L,n-1)
            if audit:
                Aprev,A=A,(2*k+1)*A+k*k*Aprev
                Bprev,B=B,(2*k+1)*B+k*k*Bprev
                fac*=n
        C=odd(U)
        G0=odd(L//gcd(L,n))
        Gold=gcd(C,G0)
        G,removed=prune_overvaluation(n,C,Gold,primes)
        if removed:
            pruned_indices+=1
            removed_factor_count+=len(removed)
        # d is an upper bound for v2(g_n/n!), sufficient even at n=3.
        d=v2(U)-n
        if quality_floor(n,d,Gold,1)>target:
            old_baseline+=1
        D=1
        bound=quality_floor(n,d,G,D)
        used=[]
        tests=0
        if bound>target:
            baseline+=1
        else:
            for p in reversed(primes):
                if p>n:
                    continue
                power=p
                while power*p<=n:
                    power*=p
                digit=n//power
                assert 1<=digit<p
                while len(Bsmall)<=digit:
                    j=len(Bsmall)-1
                    Bsmall.append((2*j+1)*Bsmall[-1]+j*j*Bsmall[-2])
                max_prefix=max(max_prefix,digit)
                tests+=1
                if C%p and Bsmall[digit]%p:
                    D*=power
                    used.append({'p':p,'prime_power':power,'leading_digit':digit})
                    bound=quality_floor(n,d,G,D)
                if bound>target:
                    break
        max_tests=max(max_tests,tests)
        total_tests+=tests
        if audit:
            assert U*fac==(1<<n)*A
            actual=F(gcd(A,4*B),fac)
            upper=F(2)**d*F(G,D)
            assert actual<=upper
            surplus=gcd(A,4*B)//gcd(gcd(A,4*B),fac)
            surplus=odd(surplus)
            assert G%surplus==0
            # The exact Legendre lower bound is 2/((2n+1)*actual**2).
            assert bound<=F(2,2*n+1)/actual**2
        if bound>target:
            certified+=1
        else:
            survivors.append(n)
        if n in (4,7,21,107,N):
            bound_record = str(bound) if bound.numerator.bit_length()<512 else {
                'numerator_bits':bound.numerator.bit_length(),
                'denominator':str(bound.denominator),
                'exact_bound_in_formula':'2*Dhat^2*2^(-2d)/((2n+1)*G^2)'}
            traces.append({'n':n,'G_old':str(Gold),'G':str(G),
                           'removed_by_initial_data':removed,'Dhat':str(D),'d':d,
                           'bound':bound_record,'prime_tests':tests,
                           'used_prime_certificates':used})
    result={'status':'PASS','max_index':N,'target':str(target),
            'gcd_backend':'system GMP' if use_gmp else 'Python stdlib',
            'full_integer_audit':audit,'certified_rejections':certified,
            'surviving_indices':survivors,'baseline_rejections':baseline,
            'old_baseline_rejections':old_baseline,
            'indices_with_stronger_G':pruned_indices,
            'removed_prime_power_count':removed_factor_count,
            'total_deficit_prime_tests':total_tests,
            'U_recurrence_steps':max(0,N-1),
            'lcm_update_steps':max(0,N-1),
            'max_prime_tests_per_index':max_tests,
            'largest_B_prefix_index_requested':max_prefix,'traces':traces,
            'elapsed_seconds':perf_counter()-started,
            'scope':'Finite exact early-stop certificates. A survivor is not '
                    'a successful approximation. No infinite quality bound.'}
    output_path=Path(output) if output else Path(__file__).with_name('adaptive_sieve_checks.json')
    output_path.write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--N',type=int,default=5000)
    ap.add_argument('--target',default='1')
    ap.add_argument('--audit',action='store_true')
    ap.add_argument('--output',default=None)
    ap.add_argument('--gmp',action='store_true',help='use optional system GMP gcd')
    args=ap.parse_args()
    run(args.N,F(args.target),args.audit,args.output,args.gmp)
