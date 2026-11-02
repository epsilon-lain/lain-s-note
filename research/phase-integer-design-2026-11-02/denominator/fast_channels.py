#!/usr/bin/env python3
"""Known-error-channel cancellation before mixing: finite Machin comparison."""
from fractions import Fraction as F
from pathlib import Path
import json
from run_design import pi_interval,absolute_interval,display

def partial(k,n):
    return sum((F((-1)**j,(2*j+1)*k**(2*j+1)) for j in range(n)),F())

def mode(k,n):
    value=partial(k,n)
    term=F((-1)**n,(2*n+1)*k**(2*n+1))
    ratio=F(-(2*n+1),(2*n+3)*k*k)
    return value+term/(1-ratio)

def blind(n):
    a,b,c=[16*partial(5,n+j)-4*partial(239,n+j) for j in range(3)]
    return (a*c-b*b)/(a-2*b+c)

def main():
    rows=[]
    for n in range(10):
        separate=16*mode(5,n)-4*mode(239,n)
        mixed=blind(n)
        rows.append({'n':n,'separate_fraction':str(separate),'mixed_fraction':str(mixed),
                     'separate_q_bits':separate.denominator.bit_length(),
                     'mixed_q_bits':mixed.denominator.bit_length()})
    bounds=pi_interval()
    for row in rows:
        separate,mixed=F(row['separate_fraction']),F(row['mixed_fraction'])
        se=absolute_interval(separate,bounds);me=absolute_interval(mixed,bounds)
        sq=[separate.denominator**2*x for x in se]
        mq=[mixed.denominator**2*x for x in me]
        lower,upper=mq[0]/sq[1],mq[1]/sq[0]
        assert lower>1
        row.update({'quality_gain_interval':[str(lower),str(upper)],
                    'quality_gain_display':display((lower+upper)/2),
                    'separate_Q_display':display(sum(sq)/2),
                    'mixed_Q_display':display(sum(mq)/2)})
    result={'scope':'Finite classical Machin/Aitken control experiment; all absolute scores poor; no new pi record.',
            'candidate_generation_uses_pi':False,'rows':rows,
            'status':'PASS','strict_joint_improvements':len(rows)}
    Path(__file__).with_name('fast_channels.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','cases':len(rows),'last_case':{k:v for k,v in rows[-1].items() if 'fraction' not in k and 'interval' not in k}},indent=2))

if __name__=='__main__':main()
