"""Numerical consistency checks, not proof of the sampling theorem."""
from cmath import exp, phase
from math import pi, ceil
import random, json
from pathlib import Path
rng=random.Random(20261102)
rows=[]
for k in [-11,-5,-1,0,1,3,8,17]:
    for shift in [0., .2, .75, 1.25, 2.]:
        for eta in [0.,.005,.02]:
            margin=abs(1-shift)
            L=2*(22/7)*abs(k)
            N=max(8,ceil(L/(2*(margin-eta)))+1)
            samples=[]
            for j in range(N):
                noise=eta*exp(2j*pi*rng.random())
                samples.append(exp(2j*pi*k*j/N)+shift+noise)
            samples.append(samples[0])
            estimate=sum(phase(samples[j+1]/samples[j]) for j in range(N))/(2*pi)
            expected=k if shift<1 else 0
            assert abs(estimate-expected)<1e-9,(k,shift,eta,N,estimate)
            assert eta+L/(2*N)<margin
            rows.append(dict(k=k,shift=shift,noise_bound=eta,N=N,winding=round(estimate),expected=expected))
# Aliasing: every sampled value is 1, but winding is any multiple of N.
for N in [8,16,32]:
    for ell in [-2,-1,0,1,3]:
        assert max(abs(exp(2j*pi*(ell*N)*j/N)-1) for j in range(N))<1e-11
out={'status':'PASS','certified_margin_examples':len(rows),'aliasing_examples':15,'note':'Analytic theorem is proved in NOTE.md; floats only verify examples.','examples':rows}
Path(__file__).with_name('checks.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='examples'}))
