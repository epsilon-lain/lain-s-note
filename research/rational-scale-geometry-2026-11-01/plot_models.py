#!/usr/bin/env python3
"""Scientific displays of the specified finite models; no inferred pi tail."""
from pathlib import Path
from fractions import Fraction as F
from math import gcd,ceil,floor
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

HERE=Path(__file__).resolve().parent
DATA=json.loads((HERE/'model_results.json').read_text())
PI=float((F(DATA['pi_lower'])+F(DATA['pi_upper']))/2)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})


def main():
    fig,axs=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    ax=axs[0]
    for q in range(1,33):
        r=1/(2*q*q)
        for p in range(floor((3.10-r)*q),ceil((3.20+r)*q)+1):
            if gcd(p,q)!=1:continue
            ax.add_patch(Circle((p/q,r),r,facecolor='#b9d8f0',edgecolor='#28669a',alpha=.5,lw=.7))
    ax.axvline(PI,color='#c03558',lw=1.6,label=r'$x=\pi$')
    ax.axvline(22/7,color='#29804d',lw=1.3,ls='--',label=r'$x=22/7$')
    ax.set(xlim=(3.10,3.20),ylim=(0,.044),xlabel='Rational point position p/q',ylabel='Height (denominator encoded)',
           title='Denominator-labelled circles: radius = 1/(2q²)')
    ax.set_aspect('equal',adjustable='box');ax.legend(loc='upper left',fontsize=9)
    ax.annotate('22/7 circle',(22/7,.0102),xytext=(3.167,.022),arrowprops={'arrowstyle':'->','color':'#28669a'},fontsize=9)
    ax.annotate('Integer 3: lower arc',(3.185,.0355),xytext=(3.145,.039),arrowprops={'arrowstyle':'->'},fontsize=8)
    ax=axs[1]
    d=np.array([r['n'] for r in DATA['plotting_midpoints_not_interval_endpoints']])
    value=np.array([r['pi_running_envelope'] for r in DATA['plotting_midpoints_not_interval_endpoints']])
    ax.step(d,value,where='post',color='#c03558',lw=2,label=r'$E_D(\pi)$')
    for D,p in ((1,3),(7,22),(113,355)):
        y=value[D-1]
        ax.scatter([D],[y],color='#c03558',s=22,zorder=3)
        ax.annotate(f'{p}/{D}',(D,y),xytext=(6,8),textcoords='offset points',fontsize=9)
    ax.text(.03,.08,'Integer 3: E_D = 0 from D = 1\nRational 22/7: E_D = 0 from D = 7',transform=ax.transAxes,
            bbox={'facecolor':'#edf5ef','edgecolor':'none','pad':6},fontsize=9)
    ax.set(yscale='log',xlim=(1,256),ylim=(.0015,.3),xlabel='Largest allowed denominator D',ylabel='Best denominator-normalized error',
           title=r'$E_D(x)=\min_{1\leq d\leq D,\ p\in\mathbb{Z}}d^2|x-p/d|$')
    ax.grid(alpha=.2);ax.legend(loc='upper right')
    fig.suptitle('Finite geometry: exact rational hits versus the pi observation path',fontsize=13)
    fig.savefig(HERE/'scale_geometry.png',dpi=160)
    plt.close(fig)

    fig,axs=plt.subplots(1,3,figsize=(12.5,4),layout='constrained')
    b=np.linspace(1.25,3.35,501)
    for ax,x,label in zip(axs[:2],(3.,PI),('Integer target x = 3','Target x = pi')):
        delta=np.min(np.abs(x-np.arange(4)[:,None]*b[None,:]),axis=0)
        V=(delta[:,None]**2+delta[None,:]**2)/2
        ax.contourf(b,b,V,levels=[0,.003,.015,.04,.10],colors=['#194a69','#447e99','#83b4bd','#d6e8df'],extend='max')
        for k in (1,2):
            for l in (1,2):ax.plot(x/k,x/l,'x',color='#be374c',ms=6)
        ax.set(title=label,xlabel='Scale b1 (real extension)',ylabel='Scale b2 (real extension)',xlim=(1.25,3.35),ylim=(1.25,3.35))
        ax.set_aspect('equal')
    ax=axs[2]
    t=np.linspace(1/3,1/2,500)
    terms=np.array([9*t-3,t,2-4*t])
    for term,label,color in zip(terms,('source 1/3','integer source 0','source 1/2'),('#28669a','#29804d','#bd4964')):
        ax.plot(t,term,ls='--',alpha=.7,label=label,color=color)
    ax.plot(t,np.min(terms,axis=0),color='#202833',lw=2.3,label='first arrival E3')
    ax.axvspan(3/8,2/5,color='#29804d',alpha=.10)
    ax.set(xlim=(1/3,1/2),ylim=(0,.5),xlabel='Target coordinate t',ylabel='Weighted arrival time',title='Farther integer can win')
    ax.legend(fontsize=8,loc='lower center');ax.grid(alpha=.15)
    fig.suptitle('Dark = smaller error; crosses = ideal alignments in the real extension',fontsize=12)
    fig.savefig(HERE/'basins_and_envelope.png',dpi=160)
    print(json.dumps({'figures':['scale_geometry.png','basins_and_envelope.png'],
                      'scope':'Rendering midpoints and continuous lens fields; rational pi bounds are in model_results.json'}))


if __name__=='__main__':main()
