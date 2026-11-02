"""Reproducible independent optional-backend checks, not needed by the sieve."""
import contextlib
import io
import json
from math import gcd
from pathlib import Path
from random import Random
from tempfile import TemporaryDirectory
import adaptive_exact_sieve as sieve


accelerated=sieve.GMPGCD()
pairs=[(0,0),(0,1<<100000),(1<<100000,0),
       (1<<100000,1<<120000),(-(1<<120000),1<<100000),
       ((1<<120000)-1,(1<<100000)-1)]
r=Random(20261102)
for _ in range(500):
    a=r.getrandbits(r.randrange(8192,60000))
    b=r.getrandbits(r.randrange(8192,60000))
    common=r.getrandbits(r.randrange(1,500))
    pairs.append((a*common,b*common))
for a,b in pairs:
    assert accelerated(a,b)==gcd(a,b)

with TemporaryDirectory() as work:
    with contextlib.redirect_stdout(io.StringIO()):
        native=sieve.run(10000,output=Path(work)/'native.json')
        gmp=sieve.run(10000,output=Path(work)/'gmp.json',use_gmp=True)
for key in native:
    if key not in ('gcd_backend','elapsed_seconds'):
        assert native[key]==gmp[key],key
result={'status':'PASS','gcd_pair_checks':len(pairs),
        'edge_cases':'zero operands, signed operands, powers of 2, powers minus 1',
        'random_seed':20261102,'backend_comparison_max_index':10000,
        'survivors':native['surviving_indices'],
        'native_seconds':native['elapsed_seconds'],
        'gmp_seconds':gmp['elapsed_seconds'],
        'gmp_abi':'conventional 64-bit Linux, 64-bit GMP limbs; checked at runtime'}
Path(__file__).with_name('gmp_backend_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
