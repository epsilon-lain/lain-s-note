"""Connect an exact integer pi source to bulk Farey-tree branch decisions.

No decimal pi values or precomputed continued-fraction terms are inputs.
Classical arctan convergents provide rigorous rational intervals.  All emitted
branch lengths and ball hits are certified by integer inequalities.
"""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
import argparse,json,sys,time

sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parent
TARGETS=[F(1,300),F(1,1000),F(1,10000)]
MAX_ORDER=2048

class ArctanSource:
    def __init__(self,t):
        self.t=t;self.U=[1,t];self.V=[0,1]
    def extend(self,n):
        while len(self.U)<=n:
            j=len(self.U)-1
            self.U.append((2*j+1)*self.t*self.U[-1]+j*j*self.U[-2])
            self.V.append((2*j+1)*self.t*self.V[-1]+j*j*self.V[-2])
    def at(self,n):
        self.extend(n)
        return F(self.V[n],self.U[n])
    def interval(self,n):
        x,y=self.at(n),self.at(n+1)
        lo,hi=(x,y) if n%2==0 else (y,x)
        assert lo<hi
        return lo,hi

def pair(x):return [str(x.numerator),str(x.denominator)]
def floor_ratio(a,b):
    assert b>0
    return a//b
def ceil_ratio(a,b):
    assert b>0
    return -((-a)//b)

def display_interval(lo,hi):
    with localcontext() as c:
        c.prec=20
        c.rounding=ROUND_FLOOR
        lower=str(Decimal(lo.numerator)/Decimal(lo.denominator))
        c.rounding=ROUND_CEILING
        upper=str(Decimal(hi.numerator)/Decimal(hi.denominator))
        return [lower,upper]

def run(targets=TARGETS,max_order=MAX_ORDER,output='results.json'):
    t0=time.monotonic();sources={5:ArctanSource(5),239:ArctanSource(239)}
    levels=[];trace=[];records=[];hits=[]
    a,b,c,d=3,1,4,1
    order=0;lo=hi=None
    def refine():
        nonlocal order,lo,hi
        order=4 if order==0 else 2*order
        if order>max_order:return False
        m=max(2,(2*order+4)//5)
        l5,u5=sources[5].interval(order);l239,u239=sources[239].interval(m)
        lo,hi=16*l5-4*u239,16*u5-4*l239
        assert 3<lo<hi<4
        levels.append({'order5':order,'order239':m,'atan5':[pair(l5),pair(u5)],
                       'atan239':[pair(l239),pair(u239)],'pi_interval':[pair(lo),pair(hi)],
                       'source5_raw_denominator_bits':sources[5].U[order+1].bit_length(),
                       'source239_raw_denominator_bits':sources[239].U[m+1].bit_length()})
        print(json.dumps({'source_level':len(levels)-1,'order5':order,
                          'completed_bulk_blocks':len(trace)},ensure_ascii=False),flush=True)
        return True
    assert refine()
    best=None;next_target=0
    while next_target<len(targets):
        # The source interval must lie strictly inside the current Farey node.
        if not(a*lo.denominator<b*lo.numerator and hi.numerator*d<c*hi.denominator):
            if not refine():break
            continue
        med=F(a+c,b+d)
        if lo>=med:direction='R'
        elif hi<=med:direction='L'
        else:
            if not refine():break
            continue
        if direction=='R':
            low_num=b*lo.numerator-a*lo.denominator
            low_den=c*lo.denominator-d*lo.numerator
            high_num=b*hi.numerator-a*hi.denominator
            high_den=c*hi.denominator-d*hi.numerator
            k0=floor_ratio(low_num,low_den)
            k1=ceil_ratio(high_num,high_den)-1
            p,q=c,d # Stationary upper endpoint, not the updated semiconvergent.
            ratio=F(b,d)
            qlo=q*(p-q*hi);qhi=q*(p-q*lo)
        else:
            low_num=c*hi.denominator-d*hi.numerator
            low_den=b*hi.numerator-a*hi.denominator
            high_num=c*lo.denominator-d*lo.numerator
            high_den=b*lo.numerator-a*lo.denominator
            k0=floor_ratio(low_num,low_den)
            k1=ceil_ratio(high_num,high_den)-1
            p,q=a,b
            ratio=F(d,b)
            qlo=q*(q*lo-p);qhi=q*(q*hi-p)
        if k0!=k1 or k0<1:
            if not refine():break
            continue
        k=k0
        assert b*c-a*d==1 and gcd(p,q)==1 and qlo>0
        branch_bound=1/(k+ratio)
        assert qhi<=branch_bound # Follows even from the coarser integer branch data.
        step=len(trace)
        trace.append({'direction':direction,'run':str(k),'source_level':len(levels)-1})
        if best is None or qhi<best:
            best=qhi
            records.append({'step':step,'p':str(p),'q':str(q),'q_digits':len(str(q)),
                            'run':str(k),'quality_interval':[str(qlo),str(qhi)],
                            'quality_decimal_enclosure':display_interval(qlo,qhi)})
        while next_target<len(targets) and qhi<targets[next_target]:
            hit={'target':str(targets[next_target]),'step':step,'direction':direction,
                 'before_node':list(map(str,[a,b,c,d])),'p':str(p),'q':str(q),
                 'q_digits':len(str(q)),'q_bits':q.bit_length(),'run':str(k),
                 'source_level':len(levels)-1,'quality_interval':[str(qlo),str(qhi)],
                 'quality_decimal_enclosure':display_interval(qlo,qhi),'branch_quality_upper':str(branch_bound),
                 'primitive_gcd':gcd(p,q)}
            hits.append(hit)
            print(json.dumps({'hit':str(targets[next_target]),'step':step,
                              'q_digits':hit['q_digits'],'run':str(k),
                              'quality_enclosure':hit['quality_decimal_enclosure']},ensure_ascii=False),flush=True)
            next_target+=1
        if direction=='R':a,b=a+k*c,b+k*d
        else:c,d=c+k*a,d+k*b
    out={'schema':1,'source':'Integer arctan recurrence with 5 and 239; no input decimal pi or CF table.',
         'targets':list(map(str,targets)),'max_source_order':max_order,
         'source_levels':levels,'bulk_trace':trace,'records':records,'hits':hits,
         'summary':{'targets_reached':len(hits),'targets_requested':len(targets),
                    'bulk_blocks':len(trace),'last_source_order5':levels[-1]['order5'],
                    'elapsed_seconds_display':round(time.monotonic()-t0,6)},
         'scope':'Finite constructive recovery of good fractions via certified covering-tree branches; no new pi record or infinite unbounded-run theorem.'}
    (ROOT/output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out['summary'],indent=2),flush=True)
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--targets',nargs='+',default=list(map(str,TARGETS)))
    parser.add_argument('--max-order',type=int,default=MAX_ORDER)
    parser.add_argument('--output',default='results.json')
    args=parser.parse_args()
    targets=list(map(F,args.targets))
    assert all(t>0 for t in targets)
    assert all(a>b for a,b in zip(targets,targets[1:]))
    assert args.max_order>=4
    run(targets,args.max_order,args.output)
