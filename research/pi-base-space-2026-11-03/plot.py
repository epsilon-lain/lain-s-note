#!/usr/bin/env python3
"""Scientific figure: observed zero gaps, denominator cost, and periodic templates."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE = Path(__file__).resolve().parent


def main():
    data = json.loads((HERE/'truncations/plot_input.json').read_text())
    periods = json.loads((HERE/'periodic/results.json').read_text())
    font = Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
    if font.exists():
        font_manager.fontManager.addfont(str(font))
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.fonttype': 'none'})
    blue, orange, green, gray = '#185f9e', '#b84d19', '#31705c', '#626877'
    fig, axes = plt.subplots(1, 2, figsize=(13.7, 6.1), layout='constrained',
                             gridspec_kw={'width_ratios': [1.1, 1]})
    fig.suptitle('进制里的“直”：平直段有多长，分母要花多少？', fontsize=19, weight='bold')
    ax = axes[0]
    pi = data['pi_decimal_nonzero_positions']
    ax.scatter([r['m'] for r in pi], [r['next_nonzero_gap_z_plus_1'] for r in pi],
               s=9, alpha=.45, color=blue, label='π：十进制实测')
    for key, label, color, marker in [
        ('classical_liouville_nonzero_positions', '经典刘维尔常数：阶乘位置放 1', orange, 'o'),
        ('square_gap_control_nonzero_positions', '对照数：平方位置放 1', green, 's')]:
        rows = data[key]
        ax.plot([r['m'] for r in rows], [r['next_nonzero_gap_z_plus_1'] for r in rows],
                color=color, marker=marker, markersize=4, lw=1.4, label=label)
    m = np.geomspace(1, 800, 150)
    ax.plot(m, m+1, color=gray, lw=1.2, ls='--', label='平方门槛：零段长度 = 前缀长度（未计约分）')
    ax.set(xscale='log', yscale='log', xlabel='截断位置 m（十进制位数）',
           ylabel='到下一个非零数字的距离（零段长度 + 1）',
           title='① 两种“越来越长的直线”，收益不同')
    ax.set_ylim(.75, 9000)
    ax.set_xlim(.9, 850)
    ax.grid(alpha=.17, which='major')
    ax.legend(loc='upper left', fontsize=8.8)
    ax.annotate('前缀 120 位 → 后续 599 个零\n'+r'分母 $10^{120}$，误差约 $10^{-720}$',
                xy=(120, 600), xytext=(30, 1800), fontsize=9, color=orange,
                arrowprops={'arrowstyle': '->', 'color': orange, 'lw': 1})
    ax = axes[1]
    stream = periods['special_bases']['114']['first_12_fractional_digits']
    positions = np.arange(1, len(stream)+1)
    ax.axvspan(.6, 2.4, color=orange, alpha=.08)
    ax.plot(positions, stream, marker='o', color=blue, lw=1.5, label='π 的真实数位')
    ax.axhline(16, color=orange, lw=1.5, ls='--', label='355/113：将数位 16 无限重复')
    ax.set(xticks=positions, xlabel='小数位位置', ylabel='114 进制的数位值',
           title='② π 在 114 进制的短暂平直，生成已知好逼近',
           ylim=(-3, 116), xlim=(.6, 12.4))
    ax.grid(alpha=.17)
    ax.legend(loc='upper right', fontsize=9)
    ax.text(.04, .9, '每个圆点是一个数位；数位可大于 9。', transform=ax.transAxes,
            fontsize=9, color=gray)
    fig.supxlabel('左：π 取前 720 位；两种稀疏数由定义计算，虚线未计额外约分。右：扫描底数 2–512、周期 1–4，共 2044 个模板。\n114 底的常数模板只匹配 π 前两位；最佳结果仍为已知的 355/113。有限实验不能判定 π 的整体逼近性质。',
                  fontsize=9.5, color=gray)
    fig.savefig(HERE/'base_space_curves.png', dpi=170)
    fig.savefig(HERE/'base_space_curves.svg')
    plt.close(fig)


if __name__ == '__main__':
    main()
