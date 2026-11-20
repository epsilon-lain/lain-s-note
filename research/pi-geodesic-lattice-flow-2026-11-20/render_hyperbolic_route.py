"""Render the half-plane route and a side view of one raised connector."""
from pathlib import Path
from math import pi
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import fontManager
from matplotlib.lines import Line2D

from verify_hyperbolic_route import farey_order, optimized_route, stack_endpoints


def arc(left, right):
    x, y = map(float, left)
    u, v = map(float, right)
    center = (u*u + v*v - x*x - y*y) / (2*(u-x))
    radius_squared = (x-center)**2 + y*y
    points = np.linspace(x, u, 120)
    return points, np.sqrt(radius_squared - (points-center)**2)


def main():
    fontManager.addfont('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
    plt.rcParams.update({'font.family': ['Noto Sans CJK JP', 'DejaVu Sans'],
                         'svg.fonttype': 'none', 'axes.unicode_minus': False})
    dark, blue, rose, teal = '#122c39', '#274bad', '#b44259', '#168374'
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.8, 4.7),
                              gridspec_kw={'width_ratios': [1.23, 1]})
    fig.patch.set_facecolor('#f7fafb')
    fig.subplots_adjust(left=.085, right=.97, bottom=.17, top=.76, wspace=.34)
    for ax in (a, b):
        ax.set_facecolor('#ffffff')
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['left', 'bottom']].set_color('#a5b6bd')
        ax.tick_params(colors=dark, labelsize=10)
        ax.grid(color='#dbe5e9', linewidth=.8, zorder=0)
        ax.set_axisbelow(True)

    cap = 6
    fractions = list(farey_order(cap))
    stacks = [stack_endpoints(p, q, cap) for p, q in fractions]
    _, route = optimized_route(stacks)
    for before, after in zip(route, route[1:]):
        x, y = arc(before[1], after[0])
        a.plot(x+3, y, color=blue, linewidth=1.8, zorder=2)
    for (p, q), (lo, hi) in zip(fractions, route):
        x = 3 + p/q
        a.plot([x, x], [float(lo[1]), float(hi[1])], color=teal, linewidth=2.0)
        a.scatter([x]*(cap//q), [1/(q*j)**2 for j in range(1, cap//q+1)],
                  color=teal, edgecolor='white', linewidth=.5, s=23, zorder=4)
    a.axvline(pi, linewidth=1.55, color=rose, linestyle='--')
    a.set(yscale='log', ylim=(.022, 1.2), xlim=(2.975, 4.025),
          xticks=[3, 3.25, 3.5, 3.75, 4], xlabel='数值坐标  x',
          ylabel='分母高度  y  （对数刻度）')
    handles = [
        Line2D([], [], color=blue, lw=2, label='双曲最短弧'),
        Line2D([], [], color=teal, lw=1.6, marker='o', ms=4, label='同值分数'),
        Line2D([], [], color=rose, lw=1.6, ls='--', label='目标竖线'),
    ]
    a.legend(handles=handles, loc='upper center', ncol=3, fontsize=8.4, frameon=False)

    b.set(xlim=(2.992, 3.178), ylim=(-.13, 1.22),
          xticks=[3, pi, 3+1/6], yticks=[0, 1],
          xlabel='数值坐标  x', ylabel='新增坐标  h')
    b.set_xticklabels(['3', 'π', '19/6'])
    b.plot([3, 3+1/6], [0, 0], color='#9fadb2', lw=1.3, ls='--')
    b.plot([pi, pi], [0, 1.1], color=rose, alpha=.35, lw=1.2, ls=':')
    b.plot([3, 3, 3+1/6, 3+1/6], [0, 1, 1, 0], lw=2.4, color=blue)
    b.scatter([3, 3+1/6], [0, 0], c=teal, s=35, zorder=4)
    b.scatter([pi], [0], c=rose, s=40, zorder=4)
    b.annotate('目标只在 h=0', xy=(pi, 0), xytext=(3.055, .31), color=rose,
               fontsize=10, arrowprops={'arrowstyle': '->', 'color': rose})
    b.text(3.024, .90, '抬高后经过 x=π', color=blue, fontsize=10)

    a.set_title('在双曲平面里选走法', fontsize=14, color=dark, loc='left', y=1.045)
    b.set_title('加一个方向，绕过人工零点', fontsize=14, color=dark, loc='left', y=1.045)
    fig.text(.085, .925, '单条曲线怎样经过所有分数', fontsize=23, weight='bold', color=dark)
    fig.text(.085, .87, '左图为 3≤p/q≤4、q≤6 的路线。右图投影到 (x,h)，省去 y 坐标。',
             fontsize=11, color='#456371')
    target = Path(sys.argv[1])
    fig.savefig(target.with_suffix('.svg'))
    fig.savefig(target.with_suffix('.png'), dpi=150)


if __name__ == '__main__':
    main()
