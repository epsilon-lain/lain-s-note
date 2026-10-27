"""Optional visual of the proved real-base cells; requires matplotlib."""
from pathlib import Path
from fractions import Fraction as F
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).parent
data=json.loads((ROOT/'base_space_results.json').read_text())
rows={r['prefix_length']:r for r in data['rows'] if r['digit']==1}
pi_fraction=np.pi-3
beta_zero=1+1/pi_fraction

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13.2,5.4),layout='constrained',
                            gridspec_kw={'width_ratios':[1.1,1]})
colors=['#2864a3','#118375','#a35180']
for n,color in zip((1,2,3),colors):
    b=np.linspace(n+1.001,24,600)
    ax1.plot(b,n/(b-1),lw=2.3,color=color,label=f'Horizontal digit {n}')
    intersection=1+n/pi_fraction
    ax1.scatter([intersection],[pi_fraction],s=55,color=color,zorder=4)
ax1.axhline(pi_fraction,color='#dc7a22',ls='--',lw=1.8,label='Target: pi − 3')
ax1.axvline(8,color='#8c9199',ls=':',lw=1.4)
ax1.text(8.15,.86,'integer base 8',color='#686f77',fontsize=9)
ax1.set_xlim(2,24);ax1.set_ylim(0,1)
ax1.set_xlabel('Real base β');ax1.set_ylabel('Fractional value')
ax1.set_title('Exact horizontal-wave curves',fontweight='bold')
ax1.grid(alpha=.17);ax1.legend(loc='upper right',fontsize=9)

for i,k in enumerate([2,3,4,5]):
    row=rows[k]
    left=float((F(row['left_root_interval_exact'][0])+F(row['left_root_interval_exact'][1]))/2)
    right=float((F(row['right_root_interval_exact'][0])+F(row['right_root_interval_exact'][1]))/2)
    ax2.barh(i,right-left,left=left,height=.55,color=colors[0],alpha=.75)
    if k==5:
        ax2.plot([left,right],[i,i],color='#2864a3',lw=4)
ax2.axvline(8,color='#cf6d26',lw=2,label='Integer base 8')
ax2.axvline(beta_zero,color='#183c62',ls='--',lw=1.8,label='Infinite horizontal wave')
ax2.set_xlim(7.90,8.72);ax2.set_yticks(range(4),['K = 2','K = 3','K = 4','K = 5'])
ax2.invert_yaxis();ax2.set_xlabel('Real base β')
ax2.set_title('First K digits all equal 1',fontweight='bold')
ax2.grid(axis='x',alpha=.17);ax2.legend(loc='lower right',fontsize=9)
fig.suptitle('Base space for the fractional part of π',fontsize=15,fontweight='bold')
fig.savefig(ROOT/'base_space_diagram.png',dpi=150)
print('Saved base_space_diagram.png')
