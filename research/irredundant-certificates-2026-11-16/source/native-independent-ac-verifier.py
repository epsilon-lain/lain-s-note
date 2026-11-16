"""Independent umbrella verification: raw action trees, string cancellation.
Assumes frozen archived 14-action specification; not official-core validation.
"""
from collections import defaultdict, Counter
from pathlib import Path
import hashlib
import itertools
import json

START=("BA","BBAAbba")
TARGET=("a","b")
INVERSE=(0,1,3,2,5,4,7,6,9,8,11,10,13,12)
OLD=(7,4,12,1,12,2,5,5)
TRANSLATE=str.maketrans("aAbB","AaBb")
INTEGER={"a":1,"A":-1,"b":2,"B":-2}

def inverse(w):
    return w.translate(TRANSLATE)[::-1]

def normal(w):
    while True:
        nxt=w.replace("aA","").replace("Aa","").replace("bB","").replace("Bb","")
        if nxt==w:
            return w
        w=nxt

def step(s,m):
    a,b=s
    ia,ib=inverse(a),inverse(b)
    options=[
        (ia,b),(a,ib),(a+b,b),(a+ib,b),(a,b+a),(a,b+ia),
        ("a"+a+"A",b),("A"+a+"a",b),
        ("b"+a+"B",b),("B"+a+"b",b),
        (a,"a"+b+"A"),(a,"A"+b+"a"),
        (a,"b"+b+"B"),(a,"B"+b+"b")]
    x,y=options[m]
    return normal(x),normal(y)

def raw_tree(source):
    layers=[[(source,())]]
    seen={source}
    balls=[{source}]
    for d in range(4):
        layer=[]
        for state,word in layers[-1]:
            for m in range(14):
                nxt=step(state,m)
                assert step(nxt,INVERSE[m])==state
                layer.append((nxt,word+(m,)))
        layers.append(layer)
        seen=seen|{s for s,w in layer}
        balls.append(seen)
    exact4=defaultdict(list)
    for state,word in layers[4]:
        exact4[state].append(word)
    return balls,exact4,[len(layer) for layer in layers]

def state_hash(states):
    numeric=[tuple(tuple(INTEGER[c] for c in w) for w in s) for s in states]
    raw=json.dumps(sorted(numeric),separators=(",",":")).encode("ascii")
    return hashlib.sha256(raw).hexdigest()

def replay(moves):
    states=[START]
    for m in moves:
        states.append(step(states[-1],m))
    assert states[-1]==TARGET
    sizes=[sum(map(len,s)) for s in states]
    return states,sizes

def sig(states):
    sizes=[sum(map(len,s)) for s in states]
    return len(states)-1,max(sizes),sum(sizes[1:])

def product(a,b):
    return a[0]+b[0],max(a[1],b[1]),a[2]+b[2]

def tokens(word,macro):
    n=len(word)
    dp=[n+1]*(n+1)
    dp[0]=0
    for i in range(n):
        dp[i+1]=min(dp[i+1],dp[i]+1)
        j=i+len(macro)
        if j<=n and word[i:j]==macro:
            dp[j]=min(dp[j],dp[i]+1)
    return dp[n]

def main():
    f,front,raw_f=raw_tree(START)
    b,back,raw_b=raw_tree(TARGET)
    assert not f[3]&b[4]
    assert len(f[3])==1329 and len(f[4])==11985 and len(b[4])==5115
    assert state_hash(f[3])=="e1f6dd0863a66abdb4bae8809fa1e39c6e66ee274d95549822f75c0a906bedae"
    assert state_hash(b[4])=="3737904cb75b8f52f58027c33bd8f5c22a8f5a3379ac844961cf6d6db3158bfc"
    middle=sorted(front.keys() & back.keys())
    paths=sorted(p+tuple(INVERSE[m] for m in reversed(q))
        for s in middle for p in front[s] for q in back[s])
    assert len(paths)==len(set(paths))==8 and OLD in paths and len(middle)==2
    records=[]
    hist=Counter()
    for word in paths:
        states,sizes=replay(word)
        peak,work=max(sizes),sum(sizes)
        hist[(peak,work)]+=1
        for cut in range(9):
            assert product(sig(states[:cut+1]),sig(states[cut:]))==sig(states)
            left=sum(sum(map(len,s)) for s in states[:cut+1])
            right=sum(sum(map(len,s)) for s in states[cut:])
            assert left+right-sizes[cut]==work
        records.append(dict(moves=word,lengths=sizes,peak=peak,work=work))
    assert hist=={(9,48):2,(9,49):4,(9,51):2}
    assert next(x["work"] for x in records if tuple(x["moves"])==OLD)==49
    candidates=sorted({p[i:i+k] for p in paths for k in (2,3,4)
                       for i in range(9-k)})
    costs=[(len(m)+sum(tokens(p,m) for p in paths),m) for m in candidates]
    best=min(costs)[0]
    winners=[m for cost,m in costs if cost==best]
    single=[min([8]+[len(m)+tokens(p,m) for m in candidates]) for p in paths]
    assert best==55 and (2,5,5) in winners and single==[8]*8
    result=dict(format="umbrella-independent-ac-audit/1",status="passed",
        method="All 14^4 labelled words each side, string pair cancellation, prefix token DP; not BFS parent enumeration",
        challenge="ac-08491",raw_layer_counts=[raw_f,raw_b],
        ball_counts=dict(forward3=len(f[3]),forward4=len(f[4]),backward4=len(b[4])),
        hashes=dict(forward3=state_hash(f[3]),backward4=state_hash(b[4])),
        intersections=dict(radius3_4=0,radius4_4=len(middle)),
        all_length8_paths=records,
        distribution=[dict(peak=k[0],work=k[1],count=v) for k,v in sorted(hist.items())],
        flat_macro=dict(candidates=len(candidates),baseline=64,best_paid_cost=best,
            optimal_macros=winners,single_path_paid_optima=single,
            language="zero/one flat macro, body length2..4, each primitive/ref1, definition body length, header0"),
        adopted_math=[dict(id="M-5606eba96406db45",revision=1),
                      dict(id="external/ac-08491-distance-8",revision=1)],
        limitations="Shared frozen move-table assumption; no official CLI, no platform upload, no optimization over longer paths or nested grammars")
    Path("/output/umbrella-verify.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,sort_keys=True))

if __name__=="__main__":
    main()
