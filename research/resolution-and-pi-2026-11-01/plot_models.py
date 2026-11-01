#!/usr/bin/env python3
"""Plot the controlled quantization comparison; numerical renderings only."""
from pathlib import Path
from fractions import Fraction as F
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent


def main():
    data=json.loads((HERE/'mechanism_results.json').read_text())
    fig,axs=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
    ax=axs[0]
    t=np.linspace(0,np.pi/2,2001)
    points=np.c_[np.cos(t),np.sin(t)]
    rounded=np.floor(points*7+.5)/7
    ax.plot(points[:,0],points[:,1],color='#252f3f',lw=2,label='Exact unit circle')
    ax.plot(rounded[:,0],rounded[:,1],color='#8295b1',lw=1.5,marker='.',ms=2,label='Observed on grid h = 1/7')
    for true,seen,color,label in [((3/5,4/5),(4/7,6/7),'#267b51','Outward: radius squared = 52/49'),
                                 ((5/13,12/13),(3/7,6/7),'#be4963','Inward: radius squared = 45/49')]:
        ax.scatter(*true,color='#252f3f',s=22,zorder=4)
        ax.scatter(*seen,color=color,s=45,zorder=5,label=label)
        ax.annotate('',xy=seen,xytext=true,arrowprops={'arrowstyle':'->','color':color,'lw':1.5})
    ax.set_aspect('equal');ax.set(xlim=(-.03,1.08),ylim=(-.03,1.08),xlabel='Position coordinate 1',ylabel='Position coordinate 2',
                                title='Apparent inward/outward motion from resolution')
    ax.legend(loc='lower left',fontsize=8);ax.grid(alpha=.15)
    ax=axs[1]
    row=next(r for r in data['pi_observations'] if r['q']==113)
    powers=[r['power'] for r in row['observations']]
    seen=[float(F(r['measured_kappa'])) for r in row['observations']]
    exact=row['true_kappa_display']
    ax.axhline(exact,color='#252f3f',ls='--',lw=1.8,label='True positive normalized error')
    ax.plot(powers,seen,'o-',color='#be4963',lw=1.8,label='Error inferred from quantized pi')
    ax.set(xlim=(.7,6.4),ylim=(-.00035,.0042),xticks=powers,
           xlabel='Resolution exponent k: grid step = 113^(-k)\nLarger k means finer observation',
           ylabel='113² times the distance to 355/113',title='An observed zero can hide a real positive error')
    ax.text(1.3,.00035,'Reported zero at coarse resolutions',fontsize=9,color='#9b304c')
    ax.text(3.2,exact+.00018,'True value about 0.00340631',fontsize=9)
    ax.legend(fontsize=8,loc='upper left');ax.grid(alpha=.2)
    fig.suptitle('Specified measurement models, not a physical field or a pi theorem',fontsize=12)
    fig.savefig(HERE/'resolution_controls.png',dpi=160)
    print(json.dumps(dict(output='resolution_controls.png',scope='Finite controlled displays; exact evidence is in the JSON and proofs.')))


if __name__=='__main__':main()
