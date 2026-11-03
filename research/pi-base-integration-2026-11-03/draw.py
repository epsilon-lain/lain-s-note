#!/usr/bin/env python3
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE=Path(__file__).resolve().parent


def main():
    font=Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
    if font.exists():
        font_manager.fontManager.addfont(str(font))
        plt.rcParams['font.family']=font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none'})
    data=json.loads((HERE/'results.json').read_text())
    rows=data['scan'];curve={r[0]:r for r in data['constructor']['curve']}
    fig,axes=plt.subplots(1,2,figsize=(13,5.4))
    fig.subplots_adjust(left=.07,right=.98,bottom=.23,top=.78,wspace=.24)
    fig.suptitle('旧有理构造可以编译成整段进制曲线的模板',fontsize=18,weight='bold',y=.96)
    ax=axes[0];short=[r for r in rows if r[0]<=192]
    ax.plot([r[0] for r in short],[r[4] for r in short],color='#be7c24',lw=1.2,label='有理模板 16/113')
    ax.scatter([r[0] for r in short],[r[1] for r in short],s=10,color='#2873a7',label='π 的几何函数',zorder=3)
    ax.set(xlim=(0,194),ylim=(-.007,.26),xlabel='整数进制 b',ylabel='J 函数值',title='整个区间内，两条曲线几乎重合')
    ax.legend(fontsize=10,loc='lower right',framealpha=.95)
    ax=axes[1]
    ax.axhspan(-1,1,color='#b6ccdd',alpha=.3,label='已有构造证书给出的范围')
    xs=[r[0] for r in rows]
    ys=[(r[1]-r[4])/curve[r[0]][2] for r in rows]
    ax.scatter(xs,ys,s=13,color='#2873a7',alpha=.8)
    ax.axhline(0,color='#7f8b96',lw=.7)
    ax.set_xscale('log');ax.set_ylim(-1.09,1.09)
    ax.set(xlabel='整数进制 b（对数轴）',ylabel='实际曲线差 ÷ 证书允许误差',title='放大差异后，仍在同一误差带内')
    ax.legend(fontsize=9,loc='lower left',framealpha=.95)
    for ax in axes:
        ax.grid(alpha=.16);ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    fig.text(.5,.075,'模板及误差带仅用旧构造器的 355/113 和质量上界 0.00366 编译；π 只用于事后验证。\n实际检查 516 个底数：2…512 全部，以及 1024、2048、4096、8192、12769；全区间结论由证明给出。',ha='center',va='center',fontsize=10,color='#52616f')
    fig.savefig(HERE/'integration.png',dpi=170)
    fig.savefig(HERE/'integration.svg')
    plt.close(fig)


if __name__=='__main__': main()
