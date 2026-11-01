from fractions import Fraction as F
from itertools import product
from pathlib import Path
from collections import Counter
import sys,json,random,time
import matrix_core as ex

P=ex.PRIME

def ceilf(x):
 x=F(x);return -((-x.numerator)//x.denominator)

def rows(H,weights,K,theta,v0):
 out=[]
 for j in range(K):
  for b in product(*(range(int(H*theta/w)+1) for w in weights)):
   ns=max(0,ceilf((H-sum(w*x for w,x in zip(weights,b))/theta)/v0))
   out.extend((j,s,b) for s in range(ns))
 return out

def params():
 for weights in ((1,3),(1,4),(1,5),(2,5),(1,2,5)):
  for H in (8,10,12,16):
   for theta in (F(1,3),F(1,4)):
    yield H,weights,2,theta,F(1,2),F(1)

def main():
 results=[]
 for H,weights,K,theta,v0,w0 in params():
  m=len(weights)
  assert v0<w0 and K*theta**m<1 and K*w0/v0*theta**m<1
  rr=rows(H,weights,K,theta,v0)
  cc=ex.columns(m,H,weights,w0)
  if len(rr)>len(cc) or len(rr)>220:continue
  base=[r[2] for r in rr]
  rank,piv=ex.rank_pivots(ex.vectors(cc,rr,base,W=2))
  record={'H':H,'weights':weights,'K':K,'theta':str(theta),'v0':str(v0),'w0':str(w0),
          'rows':len(rr),'columns':len(cc),'base_rank':rank}
  if rank==len(rr):
   rng=random.Random(2718+H+sum(weights));best=base; cp=ex.cost_parts(rr,base,H,weights)
   moves=0
   for _ in range(150):
    idx=rng.randrange(len(rr));j,s,beta=rr[idx]
    if j==0:continue
    a=list(best[idx]);i=rng.randrange(m);a[i]+=1
    if sum(w*x for w,x in zip(weights,a))>H:continue
    cand=best[:];cand[idx]=tuple(a)
    if len(set((r[0],r[1],a) for r,a in zip(rr,cand)))<len(rr):continue
    ncp=ex.cost_parts(rr,cand,H,weights)
    if ncp[0]+F(1,6)*ncp[1] >=cp[0]+F(1,6)*cp[1]:continue
    nr,np=ex.rank_pivots(ex.vectors(cc,rr,cand,W=2))
    if nr==len(rr):best,cp,piv=cand,ncp,np;moves+=1
   record.update({'accepted_moves':moves,'cost_parts':[str(x) for x in cp],
                  'row_data':rr,'column_data':cc,'assignment':best,'minor_indices':piv})
  results.append(record)
  print({k:v for k,v in record.items() if k not in ('row_data','column_data','assignment','minor_indices')},flush=True)
 Path(__file__).with_name('results.json').write_text(json.dumps(results,indent=2))
 selected_keys={(8,(1,5),'1/4'),(12,(2,5),'1/4'),(8,(1,2,5),'1/3'),(16,(1,5),'1/3')}
 selected=[{k:v for k,v in x.items() if k not in ('row_data','column_data')} for x in results if (x['H'],tuple(x['weights']),x['theta']) in selected_keys]
 assert len(selected)==4
 Path(__file__).with_name('fixtures.json').write_text(json.dumps(selected,indent=2))
 summary=[{k:v for k,v in x.items() if k not in ('row_data','column_data','assignment','minor_indices')} for x in results]
 Path(__file__).with_name('scan_summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__':main()
