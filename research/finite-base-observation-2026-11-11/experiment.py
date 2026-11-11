#!/usr/bin/env python3
"""Exact finite-base pi experiments. Standard library; no floating point.

The power-grid scan is not complete for reduced denominator quality. The
separate continued-fraction/semiconvergent scan is complete for Q < 1 up to
the stated denominator bound, by the classical intermediate-convergent test.
Finite experiments do not prove a statement at all denominator heights.
"""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import hashlib
import json

BITS = 7200
SCALE = 1 << BITS
DENOMINATOR_LIMIT = 10**1000
BASES = list(range(2, 37))


def dyadic_arctan(inv):
    lo = hi = k = 0
    power = inv
    while True:
        denominator = (2*k+1)*power
        if denominator.bit_length() > BITS+32:
            break
        numerator = SCALE if k % 2 == 0 else -SCALE
        lo += numerator // denominator
        hi += -((-numerator)//denominator)
        power *= inv*inv
        k += 1
    numerator = SCALE if k % 2 == 0 else -SCALE
    lo += min(0, numerator//denominator)
    hi += max(0, -((-numerator)//denominator))
    return lo, hi, k


def pi_interval():
    lo5, hi5, n5 = dyadic_arctan(5)
    lo239, hi239, n239 = dyadic_arctan(239)
    return 16*lo5-4*hi239, 16*hi5-4*lo239, {"atan_1_over_5":n5,"atan_1_over_239":n239}


def quality(p, q, lo, hi):
    """Return integer bounds divided by SCALE for reduced-denominator Q."""
    assert q>0 and gcd(p,q)==1
    pscale = p*SCALE
    if pscale < q*lo:
        return q*(q*lo-pscale), q*(q*hi-pscale)
    if pscale > q*hi:
        return q*(pscale-q*hi), q*(pscale-q*lo)
    raise AssertionError("Precision does not separate pi from a candidate")


def dec_interval(bounds, digits=24):
    factor = 10**digits
    a = bounds[0]*factor//SCALE
    b = -((-bounds[1]*factor)//SCALE)
    def out(n):
        assert n >= 0
        return str(n//factor)+"."+str(n%factor).zfill(digits)
    return [out(a),out(b)]


def record(p,q,lo,hi,**kwargs):
    bounds=quality(p,q,lo,hi)
    return {"p":str(p),"q":str(q),"quality_decimal":dec_interval(bounds),**kwargs}


def power_grids(lo,hi):
    summaries=[]
    checksum=hashlib.sha256()
    counts={"grid_scales":0,"grid_candidates":0,"noninteger_candidates":0}
    for base in BASES:
        max_n = 1000 if base==10 else 120
        B=1
        best={}
        subunit=[]
        for n in range(1,max_n+1):
            B*=base
            p=lo*B//SCALE
            assert hi*B//SCALE==p
            # Strict irrationality lets us use p+1 as the ceiling.
            near=(2*lo*B+SCALE)//(2*SCALE)
            assert (2*hi*B+SCALE)//(2*SCALE)==near
            counts["grid_scales"]+=1
            for mode,raw_p in (("truncation",p),("ceiling",p+1),("nearest",near)):
                g=gcd(raw_p,B)
                P,Q=raw_p//g,B//g
                lower,upper=quality(P,Q,lo,hi)
                checksum.update(f"{base},{n},{mode},{raw_p},{B},{g},{lower},{upper}\n".encode())
                counts["grid_candidates"]+=1
                if Q==1:
                    continue
                counts["noninteger_candidates"]+=1
                previous=best.get(mode)
                if previous is None or upper<previous[1]["_lower"]:
                    best[mode]=(upper,record(P,Q,lo,hi,n=n,raw_p=str(raw_p),raw_denominator=str(B),gcd=str(g)))
                else:
                    # This check ensures no unresolved overlap affects the winner.
                    assert lower>previous[0] or (P,Q)==(int(previous[1]["p"]),int(previous[1]["q"]))
                if best[mode][1].get("_lower") is None:
                    best[mode][1]["_lower"]=lower
                if upper<SCALE:
                    subunit.append({"n":n,"mode":mode,"p":str(P),"q":str(Q)})
                else:
                    assert lower>SCALE
        for _,item in best.values():
            item.pop("_lower",None)
        summaries.append({"base":base,"digit_range":[1,max_n],
                          "best_noninteger_by_mode":{m:item for m,(_,item) in best.items()},
                          "all_noninteger_subunit_occurrences":subunit})
    return {"scope":"Only floor, ceiling, and nearest numerators on the power grids; not all finite-base rationals",
            "counts":counts,"candidate_stream_sha256":checksum.hexdigest(),"bases":summaries}


def prime_support(n):
    support=[]
    p=2
    while p*p<=n:
        if n%p==0:
            support.append(p)
            while n%p==0:n//=p
        p+=1
    if n>1:support.append(n)
    return support


def allowed_denominator(q, base):
    for p in prime_support(base):
        while q%p==0:q//=p
    return q==1


def complete_subunit_scan(lo,hi):
    """Every Q<1 rational is a convergent or intermediate convergent.

    The interval CF algorithm certifies all digits needed through the first
    convergent exceeding DENOMINATOR_LIMIT. All intermediate denominators
    within the limit are generated; thresholds are separated exactly.
    """
    lower,upper=F(lo,SCALE),F(hi,SCALE)
    a0=lower.numerator//lower.denominator
    assert upper.numerator//upper.denominator==a0
    digits=[a0]
    previous_p,previous_q=1,0
    p,q=a0,1
    hits=[]
    all_hits=0
    candidates=0
    stream=hashlib.sha256()
    base_hits={str(b):[] for b in BASES+[113]}
    integer_hits=[]

    def check(P,Q):
        nonlocal all_hits,candidates
        candidates+=1
        l,u=quality(P,Q,lo,hi)
        stream.update(f"{P},{Q},{l},{u}\n".encode())
        if u<SCALE:
            all_hits+=1
            qualifying=[b for b in BASES+[113] if allowed_denominator(Q,b)]
            if Q==1:
                integer_hits.append(record(P,Q,lo,hi))
                return
            if qualifying:
                item=record(P,Q,lo,hi,qualifying_bases=qualifying)
                hits.append(item)
                for b in qualifying:base_hits[str(b)].append({"p":str(P),"q":str(Q)})
        else:
            assert l>SCALE, ("Threshold unresolved",P,Q)

    check(p,q)
    while q<=DENOMINATOR_LIMIT:
        a=digits[-1]
        lower,upper=1/(upper-a),1/(lower-a)
        a=lower.numerator//lower.denominator
        assert upper.numerator//upper.denominator==a, "Need more precision"
        digits.append(a)
        max_j=min(a,(DENOMINATOR_LIMIT-previous_q)//q)
        for j in range(1,max_j+1):
            check(j*p+previous_p,j*q+previous_q)
        previous_p,p=p,a*p+previous_p
        previous_q,q=q,a*q+previous_q
    return {"scope":"Complete for all reduced rationals with 1<=q<=10^1000 and q^2*abs(pi-p/q)<1, via the classical convergent-or-intermediate criterion",
            "denominator_limit":str(DENOMINATOR_LIMIT),
            "cf_digit_count":len(digits),"certified_cf_digits":digits,
            "first_convergent_denominator_above_limit":str(q),
            "candidate_count":candidates,"subunit_candidate_count_all_denominators":all_hits,
            "candidate_stream_sha256":stream.hexdigest(),
            "integer_hits":integer_hits,"noninteger_hits_for_selected_bases":hits,
            "noninteger_hits_by_base":base_hits,
            "decimal_noninteger_subunit_hits":base_hits["10"]}


def run():
    lo,hi,terms=pi_interval()
    assert 3*SCALE<lo<hi<4*SCALE
    result={"status":"PASS","arithmetic":"Exact integers and fractions.Fraction; no floating point",
            "scope":"Finite certified experiments; no all-n theorem for pi",
            "pi_interval":{"method":"16 atan(1/5)-4 atan(1/239), directed dyadic rounding plus alternating remainder",
                           "bits":BITS,"terms":terms,"lower_numerator":str(lo),"upper_numerator":str(hi),
                           "denominator":str(SCALE),"width_numerator":str(hi-lo)},
            "power_grids":power_grids(lo,hi),"complete_subunit_scan":complete_subunit_scan(lo,hi),
            "explicit_examples":[record(19,6,lo,hi,finite_in_base=6),record(22,7,lo,hi,finite_in_base=7),
                                 record(355,113,lo,hi,finite_in_base=113),record(25,8,lo,hi,finite_decimal="3.125"),
                                 record(157,50,lo,hi,finite_decimal="3.14")],
            "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    return result


if __name__=="__main__":
    result=run()
    Path(__file__).with_name('results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({"status":result["status"],"power_counts":result["power_grids"]["counts"],
                      "cf_candidates":result["complete_subunit_scan"]["candidate_count"],
                      "cf_digits":result["complete_subunit_scan"]["cf_digit_count"],
                      "noninteger_hits":result["complete_subunit_scan"]["noninteger_hits_for_selected_bases"]},ensure_ascii=False))
