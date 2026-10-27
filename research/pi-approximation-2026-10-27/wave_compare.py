"""Exact interval digit comparisons; optional plot uses matplotlib.
Run in the same folder as verify.py. No floating point decides any digit.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
from decimal import Decimal, localcontext
import json
from verify import pi_bounds

ROOT=Path(__file__).parent

def digits(x,b,count):
    out=[]
    for _ in range(count):
        x*=b
        d=x.__floor__()
        out.append(d)
        x-=d
    return out

def dec(x):
    with localcontext() as ctx:
        ctx.prec=22
        return str(Decimal(x.numerator)/Decimal(x.denominator))

def interval_constant():
    lo=sum((F(1,10**factorial(k)) for k in range(1,7)),F())
    return lo,lo+F(2,10**factorial(7))

def analyze(label,lo,hi,q,n):
    b=q-1
    assert 1<=n<=b
    r=F(n,q)
    target=digits(lo,b,10)
    assert target==digits(hi,b,10)
    rational=digits(r,b,10)
    assert rational==([n-1,b-n]*5)
    J=next(i+1 for i in range(10) if target[i]!=rational[i])
    jump=target[J-1]-rational[J-1]
    first=F(jump,b**J)
    normalized_lo,normalized_hi=q*q*(lo-r),q*q*(hi-r)
    # The entire certified error interval respects first-break bounds.
    sign=1 if normalized_lo>0 else -1
    abs_lo,abs_hi=sorted((sign*normalized_lo,sign*normalized_hi))
    assert abs_lo>0
    assert min(r,1-r)*F(q*q,b**J)<=abs_lo
    assert F(q*q*max(0,abs(jump)-1),b**J)<=abs_lo
    assert abs_hi<=F(q*q*(abs(jump)+1),b**J)
    return {'label':label,'base':b,'q':q,'numerator_fractional':n,
        'target_digits':target,'rational_digits':rational,
        'common_prefix':J-1,'first_different_position':J,'first_jump':jump,
        'relative_jump':dec(F(abs(jump),b-1)),
        'normalized_error_lower_approx':dec(abs_lo),
        'normalized_error_upper_approx':dec(abs_hi),
        'first_difference_contribution':dec(first),
        'remaining_signed_contribution_approx':dec((lo+hi)/2-r-first)}

if __name__=='__main__':
    a,c=pi_bounds();a-=3;c-=3
    rows=[analyze('pi / 22-7',a,c,7,1),
          analyze('pi / 355-113',a,c,113,16),
          analyze('pi / next convergent',a,c,33102,4687),
          analyze('pi / later convergent',a,c,1725033,244252),
          analyze('pi / large convergent',a,c,6701487259,948881364)]
    lo,hi=interval_constant()
    for k in range(1,6):
        q=10**factorial(k)
        r=sum((F(1,10**factorial(j)) for j in range(1,k+1)),F())
        rows.append(analyze('Liouville k='+str(k),lo,hi,q,int(q*r)))
    (ROOT/'wave_results.json').write_text(json.dumps(rows,indent=2))
    for row in rows:
        print(row['label'],'base',row['base'],'prefix',row['common_prefix'],
              'jump',row['first_jump'],'q2_error',row['normalized_error_lower_approx'])
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib unavailable; exact results are complete.')
    else:
        fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
        panels=[rows[0],rows[1],rows[6],rows[7]]
        titles=['Pi vs 22/7 | base 6','Pi vs 355/113 | base 112',
                'Liouville constant vs truncation | base 99',
                'Liouville constant vs truncation | base 999,999']
        for ax,row,title in zip(axes.flat,panels,titles):
            b=row['base'];t=row['target_digits'];r=row['rational_digits']
            ax.plot(range(1,11),[v/(b-1) for v in r],color='#ea8737',ls='--',lw=2,
                    marker='s',ms=5,label='Rational: repeating two-digit wave')
            ax.plot(range(1,11),[v/(b-1) for v in t],color='#3264b0',lw=1.7,
                    marker='o',ms=4,label='Target number')
            j=row['first_different_position']
            ax.axvline(j,color='#b8324d',alpha=.45,ls=':')
            ax.set_title(title,fontsize=12,fontweight='bold')
            ax.set_xlabel('Fractional digit position')
            ax.set_ylabel('Digit value / (base - 1)')
            ax.set_xticks(range(1,11));ax.set_ylim(-.04,1.06)
            ax.grid(alpha=.15)
            ax.text(.02,.04,'First difference: position '+str(j)+
                    '; digit jump = '+str(row['first_jump'])+
                    '\nq^2 x distance = '+row['normalized_error_lower_approx'],
                    transform=ax.transAxes,fontsize=9,
                    bbox={'facecolor':'white','alpha':.92,'edgecolor':'#dddddd'})
        axes[0,0].legend(loc='upper right',fontsize=8)
        fig.suptitle('Digit waves: initial overlap and the first departure',fontsize=16)
        fig.savefig(ROOT/'wave_comparison.png',dpi=160)
        print('Plot saved.')
