#!/usr/bin/env python3
"""Independent exact verification of the finite-base observation experiments.

No imports from experiment.py or smooth_scan.py. The independent pi interval
uses 4*(atan(1/2)+atan(1/3)) at 7800 bits, instead of the main generator's
Machin identity. Continued-fraction digits are obtained by the ordinary
Euclidean algorithms of both rational endpoints and their common prefix.
All finite conclusions retain the denominator bounds stated in README.md.
"""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import hashlib
import json

BITS = 7800
TERMS = 4100
SCALE = 1 << BITS


def arctan_enclosure(inverse):
    lower = upper = 0
    power = inverse
    for k in range(TERMS):
        d = (2*k+1)*power
        if k % 2:
            lower -= (SCALE+d-1)//d
            upper -= SCALE//d
        else:
            lower += SCALE//d
            upper += (SCALE+d-1)//d
        power *= inverse*inverse
    d = (2*TERMS+1)*power
    remainder_ceiling = (SCALE+d-1)//d
    if TERMS % 2:
        lower -= remainder_ceiling
    else:
        upper += remainder_ceiling
    return lower, upper


def independent_pi():
    l2,u2 = arctan_enclosure(2)
    l3,u3 = arctan_enclosure(3)
    return 4*(l2+l3),4*(u2+u3)


def euclidean_digits(numerator, denominator):
    out = []
    while denominator:
        quotient,remainder = divmod(numerator,denominator)
        out.append(quotient)
        numerator,denominator = denominator,remainder
    return out


def common_prefix(first,second):
    n=0
    while n<min(len(first),len(second)) and first[n]==second[n]:
        n+=1
    return first[:n]


def qinterval(p,q,lower,upper,scale):
    a=q*q*lower-p*q*scale
    b=q*q*upper-p*q*scale
    assert a*b>0,('Unresolved side',p,q)
    return min(abs(a),abs(b)),max(abs(a),abs(b))


def below_one(bounds,scale):
    if bounds[1]<scale:
        return True
    assert bounds[0]>scale,'Unresolved threshold'
    return False


def printed_interval(record,bounds,scale,key='quality_decimal',digits=24):
    lo,hi=map(F,record[key])
    assert lo<=F(bounds[0],scale)<F(bounds[1],scale)<=hi
    assert hi-lo==F(1,10**digits)


def terminates(q,base):
    # Repeated gcd removal, rather than the generator's prime-factor loop.
    while True:
        divisor=gcd(q,base)
        if divisor==1:
            return q==1
        q//=divisor


def verify_cf(data,pi_low,pi_high,source_pi):
    lo,hi,scale=source_pi
    low_digits=euclidean_digits(pi_low,SCALE)
    high_digits=euclidean_digits(pi_high,SCALE)
    prefix=common_prefix(low_digits,high_digits)
    source=data['complete_subunit_scan']
    limit=int(source['denominator_limit'])
    used=[]
    convergents=[(0,1),(1,0)]
    for a in prefix:
        p0,q0=convergents[-2]
        p1,q1=convergents[-1]
        pair=(a*p1+p0,a*q1+q0)
        convergents.append(pair)
        used.append(a)
        if pair[1]>limit:
            break
    else:
        raise AssertionError('Independent precision did not reach the bound')
    assert used==source['certified_cf_digits']
    assert len(used)==source['cf_digit_count']
    assert convergents[-2][1]<=limit<convergents[-1][1]
    assert str(convergents[-1][1])==source['first_convergent_denominator_above_limit']
    count=hits=0
    integer_hits=[]
    selected=[]
    by_base={str(b):[] for b in list(range(2,37))+[113]}
    digest=hashlib.sha256()

    def check(p,q):
        nonlocal count,hits
        assert 0<q<=limit and gcd(p,q)==1
        bound=qinterval(p,q,pi_low,pi_high,SCALE)
        source_bound=qinterval(p,q,lo,hi,scale)
        is_hit=below_one(bound,SCALE)
        assert is_hit==below_one(source_bound,scale)
        digest.update(f'{p},{q},{source_bound[0]},{source_bound[1]}\n'.encode())
        count+=1
        if not is_hit:
            return
        hits+=1
        if q==1:
            integer_hits.append((str(p),str(q)))
            return
        bases=[int(b) for b in by_base if terminates(q,int(b))]
        if bases:
            selected.append((str(p),str(q),bases))
            for b in bases:
                by_base[str(b)].append({'p':str(p),'q':str(q)})

    check(*convergents[2])
    for i in range(1,len(used)):
        # For digit i, convergents[i+1] is the immediately preceding one;
        # convergents[i] is its predecessor, including the 1/0 sentinel.
        pp,qq=convergents[i+1]
        rp,rq=convergents[i]
        stop=min(used[i],(limit-rq)//qq)
        for multiplier in range(1,stop+1):
            check(multiplier*pp+rp,multiplier*qq+rq)
    assert count==source['candidate_count']
    assert hits==source['subunit_candidate_count_all_denominators']
    assert digest.hexdigest()==source['candidate_stream_sha256']
    assert integer_hits==[(r['p'],r['q']) for r in source['integer_hits']]
    assert selected==[(r['p'],r['q'],r['qualifying_bases']) for r in source['noninteger_hits_for_selected_bases']]
    assert by_base==source['noninteger_hits_by_base']
    assert by_base['10']==source['decimal_noninteger_subunit_hits']==[]
    for r in source['integer_hits']+source['noninteger_hits_for_selected_bases']+data['explicit_examples']:
        printed_interval(r,qinterval(int(r['p']),int(r['q']),pi_low,pi_high,SCALE),SCALE)
    return {'independent_common_prefix_length':len(prefix),'used_cf_digits':len(used),
            'candidates':count,'subunit_candidates_all_denominators':hits,
            'noninteger_decimal_hits':0,
            'selected_hits':[{'p':p,'q':q,'bases':b} for p,q,b in selected],
            'integer_hits':[list(pair) for pair in integer_hits]}


def verify_grids(data,pi_low,pi_high,source_pi):
    lo,hi,scale=source_pi
    counters={'grid_scales':0,'grid_candidates':0,'noninteger_candidates':0}
    stream=hashlib.sha256()
    for entry in data['power_grids']['bases']:
        base=entry['base']
        assert entry['digit_range']==[1,1000 if base==10 else 120]
        best={}
        subunit=[]
        for n in range(1,entry['digit_range'][1]+1):
            B=base**n
            floor_low=pi_low*B//SCALE
            assert floor_low==pi_high*B//SCALE
            closest=(pi_low*B+SCALE//2)//SCALE
            assert closest==(pi_high*B+SCALE//2)//SCALE
            counters['grid_scales']+=1
            for mode,a in [('truncation',floor_low),('ceiling',floor_low+1),('nearest',closest)]:
                g=gcd(a,B)
                p,q=a//g,B//g
                bounds=qinterval(p,q,pi_low,pi_high,SCALE)
                source_bounds=qinterval(p,q,lo,hi,scale)
                stream.update(f'{base},{n},{mode},{a},{B},{g},{source_bounds[0]},{source_bounds[1]}\n'.encode())
                counters['grid_candidates']+=1
                if q==1:
                    continue
                counters['noninteger_candidates']+=1
                if mode not in best or bounds[1]<best[mode]['bounds'][0]:
                    best[mode]={'p':str(p),'q':str(q),'n':n,'raw_p':str(a),'raw_denominator':str(B),'gcd':str(g),'bounds':bounds}
                else:
                    assert bounds[0]>best[mode]['bounds'][1] or (str(p),str(q))==(best[mode]['p'],best[mode]['q'])
                if below_one(bounds,SCALE):
                    subunit.append({'n':n,'mode':mode,'p':str(p),'q':str(q)})
        assert subunit==entry['all_noninteger_subunit_occurrences']
        assert best.keys()==entry['best_noninteger_by_mode'].keys()
        for mode,item in best.items():
            saved=entry['best_noninteger_by_mode'][mode]
            for key in ('p','q','n','raw_p','raw_denominator','gcd'):
                assert item[key]==saved[key]
            printed_interval(saved,item['bounds'],SCALE)
    assert counters==data['power_grids']['counts']
    assert stream.hexdigest()==data['power_grids']['candidate_stream_sha256']
    return counters


def verify_smooth(data,pi_low,pi_high):
    limit=int(data['denominator_bound'])
    stored_lo=F(data['pi_certificate']['lower'])
    stored_hi=F(data['pi_certificate']['upper'])
    assert stored_lo<F(pi_low,SCALE)<F(pi_high,SCALE)<stored_hi
    # Enumerate powers of 5 first, opposite to the producer's nesting.
    numbers=[]
    five,c=1,0
    while five<=limit:
        q,a=five,0
        while q<=limit:
            if q>1:
                numbers.append((q,a,c))
            q*=2
            a+=1
        five*=5
        c+=1
    numbers.sort()
    assert len(numbers)==len(set(q for q,_,_ in numbers))
    stream=hashlib.sha256()
    best=None
    all_winners=[]
    max_down=max_up=0
    for q,a,c in numbers:
        middle=q*pi_low//SCALE
        assert middle==q*pi_high//SCALE
        # A period of admissible residues divides 10. This complete window
        # finds the nearest coprime integer on each side without the producer's
        # unbounded while loops.
        down=max(p for p in range(middle-9,middle+1) if gcd(p,q)==1)
        up=min(p for p in range(middle+1,middle+11) if gcd(p,q)==1)
        max_down=max(max_down,middle-down)
        max_up=max(max_up,up-middle-1)
        bd=qinterval(down,q,pi_low,pi_high,SCALE)
        bu=qinterval(up,q,pi_low,pi_high,SCALE)
        if bd[1]<bu[0]:
            p,bounds=down,bd
        else:
            assert bu[1]<bd[0]
            p,bounds=up,bu
        assert bounds[0]>SCALE
        stream.update(f'{q}:{a}:{c}:{down}:{up}:{p}\n'.encode())
        row=(bounds,p,q,a,c)
        all_winners.append(row)
        if best is None or bounds[1]<best[0][0]:
            best=row
        else:
            assert bounds[0]>best[0][1]
    assert len(numbers)==data['counts']['denominators']
    assert 2*len(numbers)==data['counts']['coprime_side_candidates']
    assert data['counts']['below_one']==data['counts']['ambiguous_pi_floors']==0
    assert max_down==data['counts']['max_downward_coprime_steps']
    assert max_up==data['counts']['max_upward_coprime_steps']
    assert stream.hexdigest()==data['ordered_candidate_digest_sha256']
    bounds,p,q,a,c=best
    saved=data['minimum']
    assert (p,q)==(saved['numerator'],saved['denominator'])==(25,8)
    assert [a,c]==saved['denominator_powers_2_5']
    assert F(saved['quality_lower'])<F(bounds[0],SCALE)<F(bounds[1],SCALE)<F(saved['quality_upper'])
    printed_interval(saved,bounds,SCALE,'quality_interval_decimal',30)
    ordered=sorted(all_winners,key=lambda row:row[0][1])
    for left,right in zip(ordered[:10],ordered[1:11]):
        assert left[0][1]<right[0][0]
    for row,saved in zip(ordered[:10],data['top_10']):
        bounds,p,q,_,_=row
        assert (p,q)==(saved['numerator'],saved['denominator'])
        printed_interval(saved,bounds,SCALE,'quality_interval_decimal',30)
    return {'denominators':len(numbers),'side_candidates':2*len(numbers),
            'minimum':[25,8],'all_noninteger_qualities_above_one':True,
            'ordered_candidate_digest_sha256':stream.hexdigest()}


def run():
    folder=Path(__file__).resolve().parent
    main=json.loads((folder/'results.json').read_text())
    smooth=json.loads((folder/'smooth_results.json').read_text())
    assert main['status']==smooth['status']=='PASS'
    assert main['script_sha256']==hashlib.sha256((folder/'experiment.py').read_bytes()).hexdigest()
    assert smooth['program_sha256']==hashlib.sha256((folder/'smooth_scan.py').read_bytes()).hexdigest()
    pi_low,pi_high=independent_pi()
    src=main['pi_interval']
    source_pi=tuple(int(src[k]) for k in ('lower_numerator','upper_numerator','denominator'))
    lo,hi,scale=source_pi
    assert hi-lo==int(src['width_numerator'])
    assert F(lo,scale)<F(pi_low,SCALE)<F(pi_high,SCALE)<F(hi,scale)
    cf=verify_cf(main,pi_low,pi_high,source_pi)
    grids=verify_grids(main,pi_low,pi_high,source_pi)
    smooth_check=verify_smooth(smooth,pi_low,pi_high)
    return {'status':'PASS','arithmetic':'Exact integers and Fraction; no floating point',
            'independent_pi':{'identity':'4*(atan(1/2)+atan(1/3))','bits':BITS,'terms_per_arctan':TERMS,
                              'lower_numerator':str(pi_low),'upper_numerator':str(pi_high),'denominator':str(SCALE),
                              'width_numerator':str(pi_high-pi_low),
                              'containment':'Strictly inside both stored pi intervals'},
            'cf':cf,'power_grids':grids,'smooth_denominators':smooth_check,
            'independence':{'code':'No generator imports or execution',
                            'cf':'Euclidean CF of both independent rational pi endpoints, followed by common-prefix comparison',
                            'smooth':'Prime-power enumeration with reversed nesting; fixed complete residue windows for nearest coprime numerators',
                            'stream_checks':'All main CF/grid streams and direct smooth candidate digests match'},
            'theory_review':{'metric_claim':'For fixed finite S, sum over S-smooth q of 1/q is finite. The circle measure of q^2|x-p/q|<M for some p is at most 2M/q. Borel-Cantelli, then a countable intersection over integer M, gives restricted quality tending to infinity for almost every x.',
                             'pi_limitation':'The almost-everywhere result does not certify membership for pi.',
                             'finite_completeness':'CF enumeration uses the classical convergent-or-intermediate criterion; no new proof of that theorem is claimed.'},
            'scope':['No infinite-denominator conclusion for pi','No external formal proof checker','Integer q=1 exceptions 3 and 4 remain included in the CF report'],
            'files_sha256':{name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in ['experiment.py','results.json','smooth_scan.py','smooth_results.json','verify.py']}}


if __name__=='__main__':
    output=run()
    Path(__file__).with_name('verification.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':output['status'],'cf':output['cf'],'power_grids':output['power_grids'],'smooth':output['smooth_denominators']},ensure_ascii=False))
