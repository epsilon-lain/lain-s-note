#!/usr/bin/env python3
"""Plot discrete geometric base functions and their exact second-layer correction."""
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
    rows=data['objects']['pi_minus_3']['rows']
    fig,axes=plt.subplots(1,2,figsize=(13,5.4))
    fig.subplots_adjust(left=.07,right=.98,bottom=.24,top=.78,wspace=.22)
    fig.suptitle('让进制成为自变量：面积图的边界怎样随换底变化',fontsize=18,weight='bold',y=.96)
    b=[r[0] for r in rows];j=[r[1] for r in rows]
    ax=axes[0]
    palette=['#2873a7','#b8791c','#389074','#7959a0','#bf5961','#63848a','#8c713a']
    for residue,c in enumerate(palette):
        branch=[r for r in rows if r[0]%7==residue]
        ax.plot([r[0] for r in branch],[r[1] for r in branch],lw=.8,color=c,alpha=.7)
        ax.scatter([r[0] for r in branch],[r[1] for r in branch],s=15,color=c,zorder=3,label=f'余 {residue}')
    row=next(r for r in rows if r[0]==113)
    ax.scatter([113],[row[1]],s=65,color='#b8791c',zorder=4)
    ax.annotate('113：边界接近铺整齐',(113,row[1]),xytext=(75,.065),
                arrowprops={'arrowstyle':'->','color':'#b8791c'},color='#986415',fontsize=10)
    ax.set(xlabel='整数进制 b',ylabel='边界混合程度 Jπ−3(b)',xlim=(0,194),ylim=(-.007,.26))
    ax.set_title('原几何函数：按进制除以 7 的余数分组',fontsize=12,pad=18)
    ax.legend(ncol=4,fontsize=8,loc='upper center',bbox_to_anchor=(.5,1.01),framealpha=.9)
    small=[r for r in rows if 16<=r[0]<=96]
    ax=axes[1]
    ax.plot([r[0] for r in small],[r[4] for r in small],lw=.9,ls='--',color='#bf7919',label='第二层抛物线波')
    ax.scatter([r[0] for r in small],[r[3] for r in small],s=17,color='#2873a7',label='剥离第一层后，放大 b² 倍',zorder=3)
    ax.set(xlabel='整数进制 b',ylabel='放大后的第二层变化',xlim=(14,98),ylim=(-.007,.26))
    ax.set_title('再放大一层：仍能看见波，而非一条平均线',fontsize=12,pad=18)
    ax.legend(loc='upper left',bbox_to_anchor=(.02,1.01),fontsize=9,framealpha=.92)
    for ax in axes:
        ax.grid(alpha=.18)
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    fig.text(.5,.085,'只研究整数进制；连线仅帮助读图。右图蓝点 = b²[ψ(R1) − J]，虚线 = ψ(R2)，ψ(t)=t(1−t)。\n两者差不超过 1/(4b)。所有点由数位矩形积分与多层公式两条独立计算路径夹逼。',ha='center',va='center',fontsize=10,color='#52616f')
    fig.savefig(HERE/'base_function.png',dpi=170)
    fig.savefig(HERE/'base_function.svg')
    plt.close(fig)
    # A separate diagnostic makes the universal leading shape visible.
    fig,ax=plt.subplots(figsize=(8.4,5.4))
    for name,label,c,marker in [('pi_minus_3','π − 3','#2873a7','o'),
                               ('sqrt2_minus_1','√2 − 1','#7959a0','x'),
                               ('one_third','1/3','#b8791c','s'),
                               ('liouville','刘维尔数','#389074','+')]:
        rr=[r for r in data['objects'][name]['rows'] if r[0]>=16]
        ax.scatter([r[5] for r in rr],[r[1] for r in rr],s=16,alpha=.55,label=label,color=c,marker=marker)
    tt=[i/500 for i in range(501)]
    ax.plot(tt,[t*(1-t) for t in tt],color='#243846',lw=1,label='共同主形状 r(1−r)')
    ax.set(xlabel='第一层剩余比例 r = {bα}',ylabel='Jα(b)',title='不同数共享主抛物线；随进制变化的相位轨迹不同')
    ax.grid(alpha=.2);ax.legend(fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE/'phase_comparison.svg')
    plt.close(fig)


if __name__=='__main__':
    main()
