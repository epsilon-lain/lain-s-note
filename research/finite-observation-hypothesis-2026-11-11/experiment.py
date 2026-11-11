#!/usr/bin/env python3
"""Exact global rounding-gcd maxima on a finite integer resolution range.

No infinite conclusion for pi. Standard library only; no floating point.
Legendre's criterion makes the compressed candidate list complete above G=1.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import json

BITS = 11200
SCALE = 1 << BITS
LIMIT = 10**1000
ROOT_LIMIT = isqrt(LIMIT)
DIRECT_LIMIT = 100000


def atan_interval(inv):
    lo = hi = k = 0
    power = inv
    while True:
        denominator = (2*k+1)*power
        if denominator.bit_length() > BITS+40:
            break
        numerator = SCALE if k % 2 == 0 else -SCALE
        lo += numerator//denominator
        hi += -((-numerator)//denominator)
        power *= inv*inv
        k += 1
    numerator = SCALE if k % 2 == 0 else -SCALE
    lo += min(0, numerator//denominator)
    hi += max(0, -((-numerator)//denominator))
    return lo, hi, k


def pi_interval():
    l5,u5,n5 = atan_interval(5)
    l239,u239,n239 = atan_interval(239)
    return 16*l5-4*u239,16*u5-4*l239,[n5,n239]


def round_readout(B,lo,hi):
    a = (2*B*lo+SCALE)//(2*SCALE)
    b = (2*B*hi+SCALE)//(2*SCALE)
    assert a == b, "Increase precision to decide rounding"
    return a


def quality(p,q,lo,hi):
    assert gcd(p,q)==1
    a,b = q*lo-p*SCALE,q*hi-p*SCALE
    if a>0:
        return q*a,q*b,"below"
    assert b<0
    return -q*b,-q*a,"above"


def decimal_ratio(num,den,places=30):
    s=10**places
    a=num*s//den
    b=-((-num*s)//den)
    def display(v):
        return str(v//s)+"."+str(v%s).zfill(places)
    return [display(a),display(b)]


def decimal_interval(lo,hi,places=30):
    return [decimal_ratio(lo,SCALE,places)[0],
            decimal_ratio(hi,SCALE,places)[1]]


def convergents(lo,hi):
    lower,upper=F(lo,SCALE),F(hi,SCALE)
    pm2,pm1,qm2,qm1=0,1,1,0
    digits,rows=[],[]
    while True:
        a=lower.numerator//lower.denominator
        assert a==upper.numerator//upper.denominator
        digits.append(a)
        p,q=a*pm1+pm2,a*qm1+qm2
        if q>ROOT_LIMIT:
            return digits,rows,(p,q)
        rows.append((p,q))
        pm2,pm1,qm2,qm1=pm1,p,qm1,q
        lower,upper=1/(upper-a),1/(lower-a)


def main():
    lo,hi,terms=pi_interval()
    digits,rows,next_row=convergents(lo,hi)
    centers=[]
    exact=[]
    stream=sha256()
    for index,(p,q) in enumerate(rows):
        l,u,side=quality(p,q,lo,hi)
        k_lo=q*SCALE//(2*u)
        k_hi=q*SCALE//(2*l)
        assert k_lo==k_hi, "Increase precision for maximal multiplier"
        K=k_lo
        assert 2*K*u<q*SCALE
        assert 2*(K+1)*l>q*SCALE
        if K:
            B=K*q
            assert round_readout(B,lo,hi)==K*p
            assert gcd(K*p,B)==K
        assert round_readout((K+1)*q,lo,hi)!=(K+1)*p
        centers.append({"index":index,"p":str(p),"q":str(q),
                        "full_rounding_multiplier_K":str(K),
                        "quality_decimal_interval":decimal_interval(l,u),
                        "side_of_fraction":side})
        exact.append((p,q,K,l,u))
        stream.update(f"{index},{p},{q},{K},{l},{u}\n".encode())

    def best_at(H):
        # B=1 has G=1. Any improvement must have q^2<B<=H,
        # and quality<1/2, so it appears among these convergents.
        best_num,best_den,best_index,best_k=1,1,None,None
        for index,(p,q,K,_,_) in enumerate(exact):
            if q*q>=H:
                break
            k=min(K,H//q)
            if k*best_den>best_num*q:
                best_num,best_den,best_index,best_k=k,q,index,k
        return best_num,best_den,best_index,best_k

    def winner_record(H):
        num,den,index,k=best_at(H)
        if index is None:
            return {"H":str(H),"G_ratio":["1","1"],"B":"1",
                    "A":"3","gcd":"1","index":None}
        p,q,K,l,u=exact[index]
        B=k*q
        A=round_readout(B,lo,hi)
        assert A==k*p and gcd(A,B)==k
        value=F(num,den)
        return {"H":str(H),"index":index,"p":str(p),"q":str(q),
                "full_K":str(K),"capped_K":str(k),"B":str(B),
                "A":str(A),"gcd":str(k),
                "G_ratio":[str(value.numerator),str(value.denominator)],
                "G_decimal_interval":decimal_ratio(num,den),
                "quality_decimal_interval":decimal_interval(l,u)}

    winner=winner_record(LIMIT)
    changes=[]
    previous=None
    for n in range(1,1001):
        H=10**n
        best=best_at(H)
        ratio=F(best[0],best[1])
        key=(ratio,best[2])
        if key!=previous:
            item=winner_record(H)
            item["decimal_exponent"]=n
            changes.append(item)
            previous=key

    example_indices=[0,1,3]
    assert rows[3]==(355,113)
    if winner["index"] not in example_indices:
        example_indices.append(winner["index"])
    examples=[]
    for index in example_indices:
        p,q,K,l,u=exact[index]
        assert K>0
        at=round_readout(K*q,lo,hi)
        after=round_readout((K+1)*q,lo,hi)
        examples.append({"index":index,"p":str(p),"q":str(q),
                         "full_K":str(K),"last_resolution":str(K*q),
                         "last_readout":str(at),"last_gcd":str(gcd(at,K*q)),
                         "next_multiple_resolution":str((K+1)*q),
                         "next_actual_readout":str(after),
                         "next_center_numerator":str((K+1)*p),
                         "next_multiple_exits_center":after!=(K+1)*p,
                         "G_at_last_multiple":decimal_ratio(K,q),
                         "quality_decimal_interval":decimal_interval(l,u)})

    direct_hash=sha256()
    direct_best=F(0)
    direct_records=[]
    reduced_cf={(p,q) for p,q in rows}
    improvements_over_one=0
    for B in range(1,DIRECT_LIMIT+1):
        A=round_readout(B,lo,hi)
        g=gcd(A,B)
        G=F(g*g,B)
        if G>1:
            improvements_over_one+=1
            assert (A//g,B//g) in reduced_cf
            assert (B//g)**2<B
        if G>direct_best:
            direct_best=G
            direct_records.append({"B":B,"A":A,"gcd":g,
                                   "G_ratio":[str(G.numerator),str(G.denominator)]})
        num,den,_,_=best_at(B)
        assert direct_best==F(num,den)
        direct_hash.update(f"{B},{A},{g}\n".encode())

    # A separate nested observation grid: all reduced denominators dividing D_n.
    # The baseline 3/1 has quality below 1/2, so every possible improvement is
    # a Legendre convergent; all such denominators are at most D_1000<10^500.
    D=1
    lcm_records=[]
    lcm_stream=sha256()
    last_best=None
    divisibility_candidates=0
    for n in range(1,1001):
        D=D*n//gcd(D,n)
        assert D<=ROOT_LIMIT
        best_index=0
        for index,(p,q,K,l,u) in enumerate(exact):
            if q>D:
                break
            if D%q:
                continue
            divisibility_candidates+=1
            if index==best_index:
                continue
            best_l,best_u=exact[best_index][3:5]
            if u<best_l:
                best_index=index
            else:
                assert l>best_u, "Increase precision to compare LCM-grid minima"
        p,q,K,l,u=exact[best_index]
        grid_numerator=p*(D//q)
        assert gcd(grid_numerator,D)==D//q
        lcm_stream.update(f"{n},{D},{best_index},{p},{q}\n".encode())
        if best_index!=last_best:
            lcm_records.append({"first_n":n,"D_n":str(D),"index":best_index,
                                "p":str(p),"q":str(q),
                                "grid_numerator":str(grid_numerator),
                                "grid_gcd":str(D//q),
                                "quality_decimal_interval":decimal_interval(l,u)})
            last_best=best_index
    lcm_final={"n":1000,"D_n":str(D),"D_n_decimal_digits":len(str(D)),
               "index":best_index,"p":str(p),"q":str(q),
               "quality_decimal_interval":decimal_interval(l,u)}

    result={"status":"PASS",
            "scope":"Complete maximum of G(B) for integer 1<=B<=10^1000; no infinite pi claim",
            "definitions":{"G":"gcd(round(B*pi),B)^2/B",
                           "round":"nearest integer; no ties for irrational pi",
                           "full_K":"largest positive k with round(k*q*pi)=k*p; zero if none",
                           "capped_K":"min(full_K,floor(H/q))"},
            "pi_interval":{"method":"16 atan(1/5)-4 atan(1/239), directed dyadic terms and alternating remainder",
                           "bits":BITS,"terms":terms,"lower_numerator":str(lo),
                           "upper_numerator":str(hi),"width_numerator":str(hi-lo)},
            "resolution_limit":str(LIMIT),"convergent_denominator_limit":str(ROOT_LIMIT),
            "certified_cf_digits":digits,"convergents_within_root_limit":len(rows),
            "first_convergent_above_root_limit":{"p":str(next_row[0]),"q":str(next_row[1])},
            "completeness":"G>1 forces a Legendre convergent with q^2<B<=H; all others bounded by baseline G(1)=1",
            "center_stream_sha256":stream.hexdigest(),"centers":centers,
            "global_maximum":winner,
            "power_of_ten_changes":changes,
            "examples":examples,
            "lcm_nested_grid":{"definition":"D_n=lcm(1,...,n), minimum reduced-denominator quality over all A/D_n",
                               "scope":"All layers 1<=n<=1000; all q dividing D_n, not only q<=n",
                               "record_changes":lcm_records,"final":lcm_final,
                               "divisibility_candidates":divisibility_candidates,
                               "all_layers_sha256":lcm_stream.hexdigest()},
            "direct_crosscheck":{"resolutions":[1,DIRECT_LIMIT],"count":DIRECT_LIMIT,
                                  "all_prefix_maxima_match":True,
                                  "G_above_one_count":improvements_over_one,
                                  "all_above_one_have_convergent_centers":True,
                                  "candidate_sha256":direct_hash.hexdigest(),
                                  "record_count":len(direct_records),"records":direct_records},
            "script_sha256":sha256(Path(__file__).read_bytes()).hexdigest()}
    target=Path(__file__).with_name("results.json")
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"PASS","convergents":len(rows),
                      "winner_index":winner["index"],"winner_q_digits":len(winner["q"]),
                      "G_decimal":winner["G_decimal_interval"],
                      "power_changes":len(changes),"direct_count":DIRECT_LIMIT,
                      "example_355_113_K":examples[2]["full_K"],
                      "lcm_record_layers":[r["first_n"] for r in lcm_records],
                      "bytes":target.stat().st_size},ensure_ascii=False))


if __name__=="__main__":
    main()
