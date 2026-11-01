#!/usr/bin/env python3
"""Plots of proved formulas; display arithmetic is not the certification."""
from pathlib import Path
from fractions import Fraction as F
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent


def main():
    data=json.loads((HERE/'pi_transfer_results.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    e=np.logspace(-3,np.log10(.4),120)
    volume=e**3/(4*np.sqrt((1-e*e/4)*(1-e*e)*(1-e*e/2)**2))
    axes[0].loglog(e,volume,color='#a53d59',lw=2.2,label='All infinitely many coordinates')
    axes[0].loglog(e,e**3/4,color='#326d95',lw=1.5,ls='--',label='First 3 coordinates')
    axes[0].set(xlabel='Node spacing epsilon (nodes 0, epsilon, 2 epsilon)',
                ylabel='Independence volume = sqrt(Gram determinant)',
                title='Cubic compression survives the infinite-coordinate limit')
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    rows=data['pi_candidates'];x=np.arange(len(rows));width=.32
    a=[r['scalar_kappa_display'] for r in rows]
    b=[r['hilbert_kappa_display'] for r in rows]
    axes[1].bar(x-width/2,a,width,color='#326d95',label='Original scalar q² error')
    axes[1].bar(x+width/2,b,width,color='#a53d59',label='Hilbert distance times same q²')
    axes[1].set(xticks=x,xticklabels=[f"{r['p']}/{r['q']}" for r in rows],yscale='log',
                xlabel='The same rational candidates (no new fractions)',ylabel='Normalized distance',
                title='Infinite embedding changes a fixed factor, not the exponent')
    for i,r in enumerate(rows):
        axes[1].text(i,max(a[i],b[i])*1.12,f"ratio = {r['ratio_display']:.5f}",ha='center',fontsize=9)
    axes[1].set_ylim(.001,.15);axes[1].legend(fontsize=8);axes[1].grid(axis='y',alpha=.2)
    fig.suptitle('One coherent infinite object with certified finite interfaces',fontsize=13)
    fig.savefig(HERE/'infinite_controls.png',dpi=160)
    print(json.dumps(dict(output='infinite_controls.png',scope='Finite plots of full-kernel formulas; no pi tail extrapolation.')))


if __name__=='__main__':main()
