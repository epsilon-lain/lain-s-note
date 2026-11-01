"""Render rational error intervals for two explicit non-pi toy numbers."""
from pathlib import Path
from fractions import Fraction
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
data=json.loads((root/'results_core_models.json').read_text())
plt.rcParams['svg.hashsalt']='pi-construction-models-v1'
fig,ax=plt.subplots(figsize=(7.8,4.3),layout='constrained')
for model,color,label in zip(data['continued_fractions'],['#4073b8','#ba4863'],
                             ['bounded digits [1; 1, 1, ...]','growing digits [0; 1, 2, 3, ...]']):
    rows=model['rows'];n=[r['n'] for r in rows]
    lo=[float(Fraction(r['normalized_error_lower'])) for r in rows]
    hi=[float(Fraction(r['normalized_error_upper'])) for r in rows]
    ax.fill_between(n,lo,hi,color=color,alpha=.2)
    ax.plot(n,lo,color=color,lw=1.6,label=label)
    ax.plot(n,hi,color=color,lw=1.6)
ax.set(xlabel='convergent index n',ylabel=r'$q_n^2 |x-p_n/q_n|$',
       title='Both toy numbers have irrationality exponent 2',yscale='log')
ax.text(.02,.03,'Shaded bands are rational bounds for an unknown infinite tail.\nThese models are not pi.',
        transform=ax.transAxes,fontsize=9)
ax.grid(alpha=.18);ax.legend(frameon=False,loc='upper right',fontsize=9)
fig.savefig(root/'continued_fraction_models.svg',metadata={'Date':None})
fig.savefig('/workspace/scratch/66870cb2f937/source_context/pi_models_preview.png',dpi=150)
