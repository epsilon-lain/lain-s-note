"""Illustration only. Exact claims are proved in NOTE.md and checked separately."""
from pathlib import Path
from math import factorial
from fractions import Fraction
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
HERE=Path(__file__).resolve().parent
fontpath=Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
if fontpath.exists():
    font=FontProperties(fname=str(fontpath))
    plt.rcParams['font.family']=font.get_name()
plt.rcParams.update({'svg.fonttype':'path','axes.spines.top':False,'axes.spines.right':False,'font.size':10})
def H(x,r):
    y=x.copy()
    for j in range(1,r+1):
        c=float(Fraction(2*(-1)**(j+1)*factorial(r)**2,factorial(r-j)*factorial(r+j)))
        y-=c*np.sin(2*np.pi*j*x)/(2*np.pi*j)
    return y
fig,ax=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
x=np.linspace(-1.1,1.1,1601)
ax[0].plot(x,x,color='#9ca3af',ls='--',label='reference y=x')
for T,col in [(1,'#236a8d'),(2,'#d9793d')]:ax[0].plot(x,H(T*x,1)/T,color=col,label=f'T={T}: x g(Tx)')
ax[0].scatter([-1,0,1],[-1,0,1],color='#111827',s=23,zorder=4)
ax[0].set(title='One shared g, two scales',xlabel='x',ylabel='f_T(x)');ax[0].legend(fontsize=8)
u=np.linspace(-.5,.5,1201)
for r,col in [(1,'#246a8d'),(2,'#3a8b73'),(4,'#d68631'),(8,'#8a5a9e')]:ax[1].plot(u,H(u,r),color=col,label=f'r={r}, order {2*r+1}')
ax[1].plot(u,u,color='#9ca3af',ls='--',lw=.8)
ax[1].set(title='More flatness constraints, higher correction order',xlabel='normalized error u',ylabel='next error H_r(u)');ax[1].legend(fontsize=8)
for a in ax:a.grid(alpha=.15)
fig.savefig(HERE/'shared_shapes.svg')
# QA raster is an intermediate file, not a separate deliverable.
qa=HERE.parent/'tmp';qa.mkdir(exist_ok=True)
fig.savefig(qa/'shared_shapes.png',dpi=130)
