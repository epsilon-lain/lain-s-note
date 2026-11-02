"""Independent manifest, prime-square lift, and 2-adic saturation checks."""
from math import gcd
from pathlib import Path
import argparse
import hashlib
import json


def v2(x):
    assert x>0
    return (x & -x).bit_length()-1


def valuation(x,p):
    assert x>0
    count=0
    while x%p==0:
        x//=p
        count+=1
    return count


def primes_upto(limit):
    flags=[True]*(limit+1)
    flags[:2]=[False,False]
    for p in range(2,limit+1):
        if flags[p]:
            for k in range(p*p,limit+1,p):
                flags[k]=False
    return [p for p in range(3,limit+1) if flags[p]]


def f(n):
    m,r=divmod(n,4)
    if r in (0,1):
        return 2*m
    if r==2:
        return 2*m+2
    return 2*m+3+v2(m+1)


def difference_formula(n):
    m,r=divmod(n,4)
    ones=m.bit_count()
    if r==0:
        return 4*m-3-2*ones-3*v2(m)
    if r==1:
        return 4*m+2-2*ones
    if r==2:
        return 4*m-2*ones
    return 4*m-1-2*ones-v2(m+1)


def manifest_checks(source):
    manifest=source/'GCD_RESEARCH_MANIFEST.json'
    contents=manifest.read_bytes()
    data=json.loads(contents)
    for entry in data['files']:
        content=(source/entry['path']).read_bytes()
        assert len(content)==entry['bytes']
        assert hashlib.sha256(content).hexdigest()==entry['sha256']
    return {'manifest_files_verified':len(data['files']),
            'manifest_status':'original input package integrity verified',
            'original_package_directory':str(source.resolve()),
            'original_manifest_sha256':hashlib.sha256(contents).hexdigest(),
            'original_package_declared_baseline':data.get('baseline_commit')}


def recurrence_checks(maximum=5000,prime_limit=199):
    aa,bb,dd=[1,1],[0,1],[0,2]
    for n in range(1,maximum):
        aa.append((2*n+1)*aa[n]+n*n*aa[n-1])
        bb.append((2*n+1)*bb[n]+n*n*bb[n-1])
        if n<2*prime_limit:
            dd.append((2*n+1)*dd[n]+n*n*dd[n-1]+2*aa[n]+2*n*aa[n-1])
    negative=[]
    for n in range(maximum+1):
        assert v2(aa[n])==f(n)
        g=gcd(aa[n],4*bb[n])
        assert v2(g)==f(n)-(1 if n==3 else 0)
        if n:
            actual=2+2*((n-1)-(n-1).bit_count())-f(n)-f(n-1)
            assert actual==difference_formula(n)
            if actual<0:
                negative.append(n)
            if n>=5:
                assert actual>=0
    assert negative==[3,4]
    lift_checks=transfer_checks=0
    examples=[]
    for p in primes_upto(prime_limit):
        modulus=p*p
        shifted=[1,(2*p+1)%modulus]
        for r in range(1,p-1):
            shifted.append(((2*p+2*r+1)*shifted[r]+(p+r)**2*shifted[r-1])%modulus)
        assert bb[p]%p
        assert aa[p]%p==0
        for r in range(1,p):
            n=p+r
            assert shifted[r]==(aa[r]+p*dd[r])%modulus
            assert aa[n]%modulus==shifted[r]*aa[p]%modulus
            assert bb[n]%modulus==shifted[r]*bb[p]%modulus
            transfer_checks+=3
            if aa[r]%p:
                predicted=0
            else:
                predicted=2 if (aa[r]//p+dd[r])%p==0 else 1
            assert valuation(gcd(aa[n],4*bb[n]),p)==predicted
            lift_checks+=1
            if (p,r) in ((17,4),(73,34),(131,35)):
                examples.append({'p':p,'r':r,'n':n,
                                 'A_r_div_p_mod_p':aa[r]//p%p,
                                 'D_r_mod_p':dd[r]%p,
                                 'exact_gcd_valuation':predicted})
    return {'maximum_saturation_index':maximum,
            'A_valuation_checks':maximum+1,
            'actual_gcd_saturation_checks':maximum+1,
            'difference_valuation_checks':maximum,
            'negative_adjacent_difference_indices':negative,
            'odd_prime_limit':prime_limit,
            'shifted_mod_p_squared_equalities':transfer_checks,
            'exact_top_half_lift_classifications':lift_checks,
            'examples':examples}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,
                        help='Optional original input directory containing '
                             'GCD_RESEARCH_MANIFEST.json; no implicit lookup.')
    parser.add_argument('--output',type=Path,
                        default=Path(__file__).with_name(
                            'prime_lift_and_saturation_verification.json'),
                        help='Output JSON path; defaults beside this script.')
    args=parser.parse_args()
    if args.source is not None:
        if not (args.source/'GCD_RESEARCH_MANIFEST.json').is_file():
            parser.error('--source must contain GCD_RESEARCH_MANIFEST.json')
        manifest_result=manifest_checks(args.source)
    else:
        manifest_result={
            'manifest_files_verified':None,
            'manifest_status':'not requested; original package integrity '
                              'audited in this checkpoint'}
    result={'status':'passed',
            **manifest_result,
            **recurrence_checks(),
            'scope':'Finite exact audits, with optional source-manifest '
                    'integrity verification. General '
                    'prime-square lifting and 2-adic saturation follow from '
                    'the proofs. The n<=50000 partial-prime scan was not rerun; '
                    'no aggregate odd-prime or global pi conclusion.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
