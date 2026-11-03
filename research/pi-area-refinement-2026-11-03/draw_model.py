#!/usr/bin/env python3
"""Draw equal-area nested subdivisions, with each active cell shown enlarged."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import font_manager

HERE=Path(__file__).resolve().parent


def main():
    data=json.loads((HERE/'results.json').read_text())
    font=Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
    if font.exists():
        font_manager.fontManager.addfont(str(font))
        plt.rcParams['font.family']=font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,4,figsize=(14,6.1))
    fig.subplots_adjust(left=.015,right=.99,bottom=.28,top=.78,wspace=.14)
    fig.suptitle('更自然的面积图：只细分尚未填完的那一块',fontsize=19,weight='bold',y=.97)
    blue,orange,dark,gray='#2873a7','#e7a74a','#ab5d12','#677482'
    for i,ax in enumerate(axes.flat,1):
        digit=data['pi_fractional_digits'][i-1]
        vertical=i%2==1
        template_children=16 if i==1 else 113
        if vertical:
            ax.add_patch(Rectangle((0,0),template_children/113,1,facecolor=orange,alpha=.4,lw=0))
            ax.add_patch(Rectangle((0,0),digit/113,1,facecolor=blue,lw=0))
            active=(digit/113,0,1/113,1)
        else:
            ax.add_patch(Rectangle((0,0),1,template_children/113,facecolor=orange,alpha=.4,lw=0))
            ax.add_patch(Rectangle((0,0),1,digit/113,facecolor=blue,lw=0))
            active=(0,digit/113,1,1/113)
        for j in range(114):
            if vertical: ax.plot([j/113,j/113],[0,1],color='white',lw=.23,alpha=.65)
            else: ax.plot([0,1],[j/113,j/113],color='white',lw=.23,alpha=.65)
        ax.add_patch(Rectangle((active[0],active[1]),active[2],active[3],
                              facecolor='none',edgecolor=dark,lw=1.5))
        ax.add_patch(Rectangle((0,0),1,1,facecolor='none',edgecolor=gray,lw=1))
        ax.set_aspect('equal')
        ax.set_xlim(-.015,1.015);ax.set_ylim(-.015,1.015)
        ax.set_xticks([]);ax.set_yticks([])
        for spine in ax.spines.values(): spine.set_visible(False)
        ax.set_title(f'第 {i} 层：填 {digit} 份',fontsize=13,pad=32)
        if i==1:
            desc='第 16 份尚未填完\n下一幅放大这一份'
        elif i<4:
            desc='填满 112 份，只留最后一份\n下一幅继续放大它'
        else:
            desc='填 69 份；第 70 份待细分\n其后 43 份全部为空'
        ax.set_xlabel(desc,fontsize=10,labelpad=13)
        label='1' if i==1 else f'1/113^{i-1}'
        ax.text(.5,1.02,'当前块实际面积：'+label,transform=ax.transAxes,
                ha='center',va='bottom',fontsize=9,color=gray)
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(facecolor=blue,label='π 已填满的面积'),
                        Patch(facecolor=orange,alpha=.4,label='355/113 模板与当前前缀间的空白'),
                        Patch(facecolor='none',edgecolor=dark,label='下一层继续细分的子块')],
               loc='lower center',bbox_to_anchor=(.5,.12),ncol=3,fontsize=10,frameon=False)
    fig.text(.5,.062,'整数部分 3 已略去。每幅的局部横纵坐标分别缩放，未保持统一长宽比；实际面积见上方标注，各幅面积不相加。\n先竖切、再横切，始终按同一规则细分；填满全部 113 个子块，就自然合成父块。',
             fontsize=10,color=gray,ha='center',va='center')
    fig.savefig(HERE/'nested_area.png',dpi=170)
    fig.savefig(HERE/'nested_area.svg')
    plt.close(fig)


if __name__=='__main__':
    main()
