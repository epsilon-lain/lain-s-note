#!/usr/bin/env python3
"""Scientific comparison: common cell normalization and actual radix neighborhood."""
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
    data=json.loads((HERE/'comparison.json').read_text())
    near=json.loads((HERE/'neighborhood_results.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(13,5.4),gridspec_kw={'width_ratios':[1,1.2]})
    fig.subplots_adjust(left=.16,right=.98,bottom=.23,top=.76,wspace=.28)
    fig.suptitle('纪录迁移：极小占比可以重现，附近整数进制没有一起变好',fontsize=17,weight='bold',y=.96)
    ax=axes[0]
    vals=[float(data[k]['C_display_bounds'][0]) for k in ['355_over_113','local_record','reported_record']]
    colors=['#b77318','#2873a7','#7959a0']
    for y,(v,c) in enumerate(zip(vals,colors)):
        ax.scatter([v],[y],s=85,color=c,facecolors='none' if y==2 else c,zorder=3)
        ax.annotate(f'{v:.5g}',(v,y),xytext=(0,13),textcoords='offset points',ha='center',color=c,fontsize=11)
    ax.set_yticks([0,1,2],['355/113\n内部缺口','本地纪录\n外部新增面积','巨型纪录\n依报告推导'])
    ax.set_xscale('log');ax.set_xlim(1e-13,1e-1);ax.set_ylim(2.6,-.6)
    ax.set_xticks([1e-12,1e-9,1e-6,1e-3])
    ax.set_xlabel('差异面积 ÷ 所在小块面积（越小越好）',labelpad=12)
    ax.set_title('各自放大面积为 1/q² 的小块',fontsize=12,pad=20)
    ax.grid(axis='x',alpha=.23)
    ax=axes[1]
    rows=near['neighbors_q_plus_h']
    hs=[r['h'] for r in rows]
    ys=[float(r['distance_to_integer_display_bounds'][0]) for r in rows]
    ax.plot(hs,ys,color='#a8b6c2',lw=1,zorder=1)
    ax.scatter(hs,ys,color='#2873a7',s=20,zorder=2)
    ax.scatter([0],[0],color='#b77318',s=65,zorder=3)
    ax.annotate('q：残差极小\n此比例尺下贴近 0',xy=(0,0),xytext=(8,.095),
                arrowprops={'arrowstyle':'->','color':'#b77318'},color='#985f13',fontsize=10)
    ax.set_xlim(-34,34);ax.set_ylim(-.018,.53)
    ax.set_xticks([-32,-16,0,16,32]);ax.set_xlabel('分块数 B = q + h，横轴为 h',labelpad=12)
    ax.set_ylabel('Bπ 距离最近整数的差')
    ax.set_title('本地纪录的 65 个整数进制',fontsize=12,pad=20)
    ax.grid(alpha=.2)
    for ax in axes:
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    fig.text(.5,.075,'左：统一比较相对面积，未按真实几何大小绘制。巨型纪录未重建完整 p、q。\n右：实际计算本地 233,325 位分母的邻域；所有 64 个非零偏移在约分后仍有 C > 1。仅覆盖所列有限邻域。',
             ha='center',va='center',fontsize=10,color='#52616f')
    fig.savefig(HERE/'record_migration.png',dpi=170)
    fig.savefig(HERE/'record_migration.svg')
    plt.close(fig)


if __name__=='__main__':
    main()
