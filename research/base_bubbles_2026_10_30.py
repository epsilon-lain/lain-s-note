"""Parameter-space cylinders for fixed pi-3; exact endpoints plus a plot."""
from fractions import Fraction as F
from pathlib import Path
from dyadic_rectangles_2026_10_29 import pi_bounds

def finite_value(word,b):
    return sum((F(d)*b**(-j) for j,d in enumerate(word,1)),F())

def periodic_value(prefix,cycle,b):
    return finite_value(prefix,b)+b**(-len(prefix))*finite_value(cycle,b)/(1-b**(-len(cycle)))

def decreasing_root(fn,xlo,xhi,steps=130):
    left,right=F(2),F(16)
    assert fn(left)>xhi and fn(right)<xlo
    for _ in range(steps):
        mid=(left+right)/2; value=fn(mid)
        if value>xhi:left=mid
        elif value<xlo:right=mid
        else:raise ArithmeticError('Increase pi precision')
    assert fn(left)>xhi and fn(right)<xlo
    return left,right

def midpoint(interval):return sum(interval)/2

def cylinder(word,xlo,xhi):
    # For these examples digits <=2 and all endpoints >3.  The final
    # remainder in [0,1) then guarantees all earlier remainders in [0,1).
    lower=decreasing_root(lambda b:finite_value(word,b),xlo,xhi)
    upper=decreasing_root(lambda b:finite_value(word,b)+b**(-len(word)),xlo,xhi)
    assert lower[0]>3 and lower[1]<upper[0]
    return lower,upper

def greedy(x,b,n):
    out=[]
    for _ in range(n):
        x*=b;d=x.__floor__();out.append(d);x-=d
    return out

def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plo,phi=pi_bounds();xlo,xhi=plo-3,phi-3;x=(xlo+xhi)/2
    B=(1+1/xhi,1+1/xlo);beta=midpoint(B)
    families=[([], [1], '1 repeated'),([], [1,0], '10 repeated'),([], [1,2], '12 repeated'),([1],[0,2],'1 + 02 repeated')]
    data=[]
    for prefix,cycle,label in families:
        centre=decreasing_root(lambda b:periodic_value(prefix,cycle,b),xlo,xhi)
        rows=[]
        for k in range(2,7):
            word=(prefix+cycle*(k+1))[:k];lower,upper=cylinder(word,xlo,xhi)
            assert lower[1]<centre[0]<centre[1]<upper[0]
            m=(lower[1]+upper[0])/2
            assert greedy(xlo,m,len(word))==word==greedy(xhi,m,len(word))
            if rows:assert lower[0]>=rows[-1][1][0] and upper[1]<=rows[-1][2][1]
            rows.append((k,lower,upper))
        data.append((label,centre,rows))
        print(label,'centre',f'{float(midpoint(centre)):.15f}')
    print('all-1 bubbles: N, left endpoint, right endpoint, width')
    all_one=[]
    for n in range(2,13):
        lower,upper=cylinder([1]*n,xlo,xhi)
        l,r=midpoint(lower),midpoint(upper)
        all_one.append((n,l,r))
        if n in (2,3,4,6,8,10,12):print(n,f'{float(l):.15f}',f'{float(r):.15f}',f'{float(r-l):.12g}')
    assert greedy(xlo,F(8),3)==greedy(xhi,F(8),3)==[1,1,0]
    assert greedy(xlo,F(129,16),7)==greedy(xhi,F(129,16),7)==[1]*6+[0]
    print('Rational-base checks: base 8 has 2 leading ones; base 129/16 has 6.')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(13,5.2),gridspec_kw={'width_ratios':[1.25,1]})
    colours=['#166e90','#d4742c','#4c956c','#9263a5']
    for idx,(label,centre,rows) in enumerate(data):
        c=colours[idx];bc=float(midpoint(centre))
        axes[0].axvline(bc,color=c,alpha=.28,lw=1)
        for k,lower,upper in rows:
            y=k+(idx-1.5)*.15
            axes[0].plot([float(midpoint(lower)),float(midpoint(upper))],[y,y],color=c,lw=5,solid_capstyle='round',label=label if k==2 else None)
            axes[0].plot(bc,y,'|',color='#202b38',ms=8)
    axes[0].set(xlabel='Base b (continuous)',ylabel='Number of matched digits N',title='Finite repeats form nested intervals',xlim=(6.95,9.35),ylim=(1.45,6.65))
    axes[0].legend(loc='upper right',frameon=False,fontsize=9)
    for n,l,r in all_one:
        if n<3:continue
        zl=float(beta**n*(l-beta));zr=float(beta**n*(r-beta))
        axes[1].plot([zl,zr],[n,n],color=colours[0],lw=7,alpha=.7,solid_capstyle='round')
        axes[1].plot(0,n,'|',color='#202b38',ms=12)
    axes[1].axvline(-float(beta-1),color='#d4742c',ls='--',lw=1.2,label='Limiting left edge')
    axes[1].axvline(float((beta-1)*(beta-2)),color='#4c956c',ls='--',lw=1.2,label='Limiting right edge')
    axes[1].set(xlabel='Magnified offset: B^N (b - B)',ylabel='Number of leading 1s, N',title='All-1 bubble: exponential shrinking',xlim=(-10,47),ylim=(2.4,13.8))
    axes[1].legend(loc='upper right',frameon=False,fontsize=9)
    fig.suptitle('Fixed target x = pi - 3; only the base changes',fontsize=15)
    fig.tight_layout()
    folder=Path(__file__).parent
    fig.savefig(folder/'base_bubbles_2026_10_30.png',dpi=160)
    fig.savefig(folder/'base_bubbles_2026_10_30.svg')
    print('Certified endpoint brackets have width 14 / 2^130; exact checks passed.')
if __name__=='__main__':main()
