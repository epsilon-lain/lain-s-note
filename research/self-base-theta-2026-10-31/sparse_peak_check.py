"""Exact finite checks for the sparse continued-fraction construction.
The infinite conclusions are proved in SPARSE_PEAK_OBSTRUCTION.md.
"""
import json
def digit(k):
    return k.bit_length() if k>=2 and k&(k-1)==0 else 1
qprev,q=0,1
total=0
rows=[]
for n in range(1,4097):
    a=digit(n)
    qprev,q=q,a*q+qprev
    total+=a
    J=n.bit_length()-1
    assert total==n+J*(J+1)//2
    assert total<=2*n
    assert q>=2**(n//2)
    if n in [16,64,256,1024,4096]:
        rows.append(dict(n=n,partial_quotient_sum=total,max_partial_quotient=J+1,
                         denominator_bits=q.bit_length()))
print(json.dumps(dict(status="Finite integer checks; infinite theorem proved in note",
                      checked_prefixes=4096,rows=rows),indent=2))
