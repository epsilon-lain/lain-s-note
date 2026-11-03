#!/usr/bin/env python3
"""Draw exact-node polygons with proved whole-domain error bands."""
from fractions import Fraction
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

HERE = Path(__file__).resolve().parent
font = Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
if font.exists():
    font_manager.fontManager.addfont(str(font))
    plt.rcParams['font.family'] = 'Noto Sans CJK JP'
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,
                     'font.size':11,'svg.fonttype':'none'})
data = json.loads((HERE/'results.json').read_text())
fig, axs = plt.subplots(2,1,figsize=(11.5,7.3),sharex=True)
fig.patch.set_facecolor('#f7f8fc')
for ax,curve in zip(axs,data['curves']):
    b = curve['base']
    n = curve['grid_denominator']
    x = np.arange(n+1)/n
    y = np.array(curve['value_numerators'])/curve['value_denominator']
    error = float(Fraction(curve['uniform_error_upper_exact']))
    top = float(Fraction(curve['maximum_exact']))
    avg = float(Fraction(curve['mean_exact']))
    ax.set_facecolor('white')
    ax.plot(x,y,color='#176a82',lw=1.5,label='精确节点折线（全局误差见右下）')
    ax.axhline(avg,color='#9d6b24',lw=1,ls='--',label=f'全区间平均值 {curve["mean_exact"]}')
    if b == 2:
        points = np.array([1/4,3/4])
        ax.scatter(points,[top]*len(points),s=32,color='#b33d6a',zorder=4,label='全部最大点')
        subtitle = '最大值 1/4；最大点为 1/4、3/4'
    else:
        points = [float(Fraction(z)) for z in data['ternary_maximum_witnesses_exact']]
        ax.scatter(points,[top]*len(points),marker='|',s=60,color='#b33d6a',zorder=4,
                   label='最大点集的 81 个精确见证')
        subtitle = '最大值 3/13；所有偶数位为 1 的三进制地址都达到它'
    ax.set_title(f'固定进制 b={b}：{subtitle}',loc='left',fontsize=12,pad=11)
    ax.set_ylabel('几何读数 $F_b(x)$')
    ax.set_ylim(-.012,.275)
    ax.grid(axis='y',alpha=.16)
    ax.text(.99,.05,f'{n:,} 个子区间；全局误差 ≤ {error:.3g}',
            transform=ax.transAxes,ha='right',color='#46526b',fontsize=10)
    ax.legend(loc='upper right',fontsize=8,framealpha=.9)
axs[-1].set_xlabel('面积 / 实数地址 x（横轴遍历整个 [0,1]，进制保持固定）')
axs[-1].set_xlim(0,1)
fig.suptitle('一个递推方程，同时确定整个实数区间的曲线',x=.08,ha='left',fontsize=17,y=.98)
fig.text(.08,.925,'原始矩形积分模型 · 全区间定理 + 精确有理数验证',color='#46526b')
fig.subplots_adjust(left=.08,right=.98,bottom=.09,top=.85,hspace=.4)
fig.savefig(HERE/'whole_interval.png',dpi=170,facecolor=fig.get_facecolor())
fig.savefig(HERE/'whole_interval.svg',facecolor=fig.get_facecolor())
print('Saved whole_interval.png and whole_interval.svg')
