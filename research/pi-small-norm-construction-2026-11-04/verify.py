"""Independent exact checks: backward fractions, tree boundaries, Taylor bounds.

Does not import construct.py or use a decimal value / precomputed CF of pi.
All decisions use integers or fractions. Decimal enclosures are checked as
fractions, so printed digits cannot overstate the certificate's accuracy.
"""
from fractions import Fraction as F
from functools import cache
from math import gcd
from pathlib import Path
import argparse,hashlib,json,sys,time

sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parent

def read_pair(pair):return F(int(pair[0]),int(pair[1]))

@cache
def backward(t,n):
    tail=F((2*n-1)*t)
    for j in range(n-1,0,-1):
        tail=(2*j-1)*t+F(j*j)/tail
    return 1/tail

def atan_taylor(t,bits):
    """Directed integer rounding of alternating Taylor terms on a dyadic grid.

    Before the first omitted term, every term is enclosed to one grid unit;
    the entire remaining alternating tail is strictly less than one unit.
    """
    scale=1<<bits
    power=t;low=high=0;j=0
    while True:
        down,rem=divmod(scale,(2*j+1)*power)
        if not down:
            if j%2:low-=1
            else:high+=1
            break
        up=down+bool(rem)
        if j%2:low-=up;high-=down
        else:low+=down;high+=up
        power*=t*t;j+=1
    return F(low,scale),F(high,scale),j

def gaussian_multiply(x,y):
    return x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]

def check_machin():
    z=(1,0)
    for _ in range(4):z=gaussian_multiply(z,(5,1))
    assert z==gaussian_multiply((2,2),(239,1))==(476,480)
    # 0 < 4 atan(1/5) - atan(1/239) < 4/5 < pi/2.
    # Positivity uses atan x > x-x^3/3 and atan x < x.
    assert 4*(F(1,5)-F(1,3*5**3))-F(1,239)>0

def tree_boundary_check(node,direction,k,lo,hi):
    """Check the k-th and (k+1)-st mediants, not a floor-ratio formula."""
    a,b,c,d=node
    assert b*c-a*d==1 and F(a,b)<lo<hi<F(c,d)
    assert k>=1
    if direction=='R':
        v=F(a+k*c,b+k*d);w=F(a+(k+1)*c,b+(k+1)*d)
        assert v<=lo<hi<=w
        p,q=c,d
        bounds=(q*(p-q*hi),q*(p-q*lo))
        new=(a+k*c,b+k*d,c,d)
    else:
        assert direction=='L'
        v=F(c+k*a,d+k*b);w=F(c+(k+1)*a,d+(k+1)*b)
        assert w<=lo<hi<=v
        p,q=a,b
        bounds=(q*(q*lo-p),q*(q*hi-p))
        new=(a,b,c+k*a,d+k*b)
    assert gcd(p,q)==1 and bounds[0]>0
    return new,(p,q),bounds

def check_decimal(item):
    lo,hi=map(F,item['quality_interval'])
    dl,du=map(F,item['quality_decimal_enclosure'])
    assert dl<=lo<hi<=du

def verify(path):
    start=time.monotonic()
    data=json.loads(path.read_text())
    check_machin()
    levels=data['source_levels'];trace=data['bulk_trace']
    intervals=[];source_checks=0;gap_bits=0
    for lev in levels:
        for t in (5,239):
            n=lev['order'+str(t)]
            pair=sorted([backward(t,n),backward(t,n+1)])
            assert pair==list(map(read_pair,lev['atan'+str(t)]))
            gap=pair[1]-pair[0]
            gap_bits=max(gap_bits,gap.denominator.bit_length()-gap.numerator.bit_length()+2)
            source_checks+=2
        l5,u5=map(read_pair,lev['atan5'])
        l239,u239=map(read_pair,lev['atan239'])
        interval=(16*l5-4*u239,16*u5-4*l239)
        assert interval==tuple(map(read_pair,lev['pi_interval']))
        assert 3<interval[0]<interval[1]<4
        intervals.append(interval)
    node=(3,1,4,1);steps=[];level_before=-1
    for i,block in enumerate(trace):
        level=block['source_level'];assert level_before<=level<len(levels)
        if i:assert block['direction']!=trace[i-1]['direction']
        old=node
        node,pq,bounds=tree_boundary_check(old,block['direction'],int(block['run']),*intervals[level])
        steps.append((old,node,pq,bounds))
        level_before=level
    # A different formula and rounding scheme supplies one very narrow interval.
    bits=max(gap_bits+256,2*max(node[1].bit_length(),node[3].bit_length())+256)
    tl5,tu5,terms5=atan_taylor(5,bits)
    tl239,tu239,terms239=atan_taylor(239,bits)
    tlo,thi=16*tl5-4*tu239,16*tu5-4*tl239
    for lev in levels:
        for t,tl,tu in [(5,tl5,tu5),(239,tl239,tu239)]:
            lower,upper=map(read_pair,lev['atan'+str(t)])
            assert lower<tl<tu<upper
    # Recheck every emitted block against Taylor, independently of source CFs.
    best=None;record_steps=[];first_hits={}
    for i,(old,new,pq,bounds) in enumerate(steps):
        block=trace[i]
        new2,pq2,tbounds=tree_boundary_check(old,block['direction'],int(block['run']),tlo,thi)
        assert new2==new and pq2==pq
        if best is None or tbounds[1]<best:
            best=tbounds[1];record_steps.append(i)
        for target in data['targets']:
            if tbounds[1]<F(target) and target not in first_hits:first_hits[target]=i
    assert record_steps==[r['step'] for r in data['records']]
    for record in data['records']:check_decimal(record)
    verified_hits=[]
    for hit in data['hits']:
        i=hit['step'];old,new,(p,q),bounds=steps[i]
        assert list(map(int,hit['before_node']))==list(old)
        assert (p,q)==(int(hit['p']),int(hit['q']))
        assert list(bounds)==list(map(F,hit['quality_interval']))
        assert bounds[1]<F(hit['target']) and first_hits[hit['target']]==i
        assert hit['primitive_gcd']==gcd(p,q)==1
        assert hit['q_digits']==len(str(q)) and hit['q_bits']==q.bit_length()
        assert hit['run']==trace[i]['run'] and hit['source_level']==trace[i]['source_level']
        assert hit['direction']==trace[i]['direction']
        ratio=F(old[1],old[3]) if hit['direction']=='R' else F(old[3],old[1])
        assert F(hit['branch_quality_upper'])==1/(int(hit['run'])+ratio)
        assert bounds[1]<=F(hit['branch_quality_upper'])
        check_decimal(hit)
        verified_hits.append({'target':hit['target'],'step':i,'q_digits':len(str(q)),
                              'quality_decimal_enclosure':hit['quality_decimal_enclosure']})
    assert set(first_hits)=={h['target'] for h in data['hits']}
    assert data['summary']['bulk_blocks']==len(trace)
    assert data['summary']['targets_reached']==len(data['hits'])
    assert data['summary']['targets_requested']==len(data['targets'])
    return {'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'checks':'PASS','backward_fraction_checks':source_checks,
            'tree_blocks_source_and_taylor':len(trace),'source_levels':len(levels),
            'independent_taylor_grid_bits':bits,'independent_taylor_terms':[terms5,terms239],
            'last_node_denominator_digits':[len(str(node[1])),len(str(node[3]))],
            'largest_certified_run':max(int(z['run']) for z in trace),
            'record_steps':record_steps,'verified_hits':verified_hits,
            'unreached_targets':[t for t in data['targets'] if t not in first_hits],
            'elapsed_seconds_display':round(time.monotonic()-start,6)}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files',nargs='*',default=['results.json','deep_scan.json'])
    parser.add_argument('--output',default='verification.json')
    args=parser.parse_args()
    report={'reports':[verify(ROOT/name) for name in args.files],
            'scope':'Independent implementation, not an external or multi-agent review; finite certificates only.'}
    (ROOT/args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
