"""Exact finite example for properties versus whole structures; no model calls.

Python >= 3.10; the default verification uses only the standard library.
Run with --plot to regenerate the figure (requires matplotlib).
"""
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import json
import math
import argparse

ROOT = Path(__file__).resolve().parent
N = 8
PAIRS = list(combinations(range(N), 2))
BITS = {p: 1 << i for i, p in enumerate(PAIRS)}


def mask_of_cycles(order, sizes):
    answer, offset = 0, 0
    for size in sizes:
        cycle = order[offset:offset + size]
        offset += size
        for i, a in enumerate(cycle):
            b = cycle[(i + 1) % size]
            answer |= BITS[tuple(sorted((a, b)))]
    return answer


def adjacency(mask):
    adj = [set() for _ in range(N)]
    for i, (a, b) in enumerate(PAIRS):
        if mask & (1 << i):
            adj[a].add(b)
            adj[b].add(a)
    return adj


def component_sizes(adj):
    unseen, sizes = set(range(len(adj))), []
    while unseen:
        stack = [unseen.pop()]
        count = 0
        while stack:
            node = stack.pop()
            count += 1
            new = adj[node] & unseen
            unseen.difference_update(new)
            stack.extend(new)
        sizes.append(count)
    return tuple(sorted(sizes))


def triangles(adj):
    return sum(b in adj[a] and c in adj[a] and c in adj[b]
               for a, b, c in combinations(range(len(adj)), 3))


def radius_one_is_p3(adj):
    return all(len(adj[v]) == 2 and
               not any(b in adj[a] for a, b in combinations(adj[v], 2))
               for v in range(len(adj)))


def describe(mask):
    adj = adjacency(mask)
    return {
        "vertices": N,
        "edges": mask.bit_count(),
        "degrees": [len(a) for a in adj],
        "triangles": triangles(adj),
        "all_induced_radius_one_neighborhoods_are_p3": radius_one_is_p3(adj),
        "component_sizes": list(component_sizes(adj)),
        "hom_c4": sum(b in adj[a] and c in adj[b] and d in adj[c] and a in adj[d]
                      for a, b, c, d in product(range(N), repeat=4)),
        "edges_list": [list(e) for i, e in enumerate(PAIRS) if mask & (1 << i)],
    }


def main(make_plot=False, font_path=None):
    # Every finite simple degree-two graph is a disjoint union of cycles >= 3.
    # The only partitions of 8 into such sizes are 8, 5+3, and 4+4.
    all_graphs = set()
    for order in permutations(range(N)):
        for sizes in ((8,), (5, 3), (4, 4)):
            all_graphs.add(mask_of_cycles(order, sizes))
    distribution = Counter()
    triangle_free = []
    for mask in all_graphs:
        adj = adjacency(mask)
        assert all(len(a) == 2 for a in adj)
        assert mask.bit_count() == 8
        distribution[component_sizes(adj)] += 1
        if triangles(adj) == 0:
            assert radius_one_is_p3(adj)
            triangle_free.append(mask)

    expected = {(8,): math.factorial(8) // (2 * 8),
                (3, 5): math.factorial(8) // ((2 * 3) * (2 * 5)),
                (4, 4): math.factorial(8) // ((2 * 4)**2 * 2)}
    assert dict(distribution) == expected
    assert len(all_graphs) == 3507
    assert len(triangle_free) == 2835

    left = mask_of_cycles(tuple(range(8)), (8,))
    right = mask_of_cycles(tuple(range(8)), (4, 4))
    facts = {"scope": "all simple undirected degree-two graphs on labels 0..7",
             "exhaustiveness_reason": "cycle decomposition and all vertex permutations",
             "labeled_graph_count": len(all_graphs),
             "cycle_type_counts": {str(k): v for k, v in sorted(distribution.items())},
             "triangle_free_labeled_count": len(triangle_free),
             "triangle_free_isomorphism_types": [[8], [4, 4]],
             "C8": describe(left), "two_C4": describe(right),
             "interpretation": "Labels are repeated representations, not prior probabilities."}
    assert facts['C8']['hom_c4'] == 48
    assert facts['two_C4']['hom_c4'] == 64
    assert facts['C8']['component_sizes'] == [8]
    assert facts['two_C4']['component_sizes'] == [4, 4]
    (ROOT / 'exact_results.json').write_text(
        json.dumps(facts, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: facts[k] for k in ('labeled_graph_count', 'cycle_type_counts',
                                          'triangle_free_labeled_count',
                                          'triangle_free_isomorphism_types')}, ensure_ascii=False))

    if make_plot:
        plot_example(left, right, font_path)


def plot_example(left, right, font_path=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.patches import Circle
    # No fontconfig command or Linux-specific path is required. A font may be
    # supplied explicitly; otherwise use an installed CJK font if available.
    if font_path:
        font_file = str(Path(font_path).expanduser().resolve())
        font_manager.fontManager.addfont(font_file)
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=font_file).get_name()
    else:
        installed = {font.name for font in font_manager.fontManager.ttflist}
        cjk_fonts = ('Noto Sans CJK SC', 'Source Han Sans SC', 'Microsoft YaHei',
                     'PingFang SC', 'SimHei', 'WenQuanYi Zen Hei', 'Noto Sans CJK JP')
        selected = next((name for name in cjk_fonts if name in installed), None)
        if selected:
            plt.rcParams['font.family'] = selected
        else:
            raise SystemExit('Exact checks passed. Plotting Chinese labels needs an installed '
                             'CJK font, or --font /path/to/font.ttf.')
    plt.rcParams['axes.unicode_minus'] = False
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.7), dpi=180)
    fig.patch.set_facecolor('#faf9f6')
    colors = ['#255b73', '#756141']
    titles = ['八个顶点连成一圈', '两个互不相连的四点方环']
    positions = [{i: (1.25*math.cos(math.pi/8+2*math.pi*i/8),
                      1.25*math.sin(math.pi/8+2*math.pi*i/8)) for i in range(8)},
                 {i: ((-0.95 if i < 4 else 0.95)+0.58*math.cos(math.pi/4+2*math.pi*(i%4)/4),
                      0.58*math.sin(math.pi/4+2*math.pi*(i%4)/4)) for i in range(8)}]
    for ax, mask, pos, color, title, count in zip(axes, [left, right], positions, colors, titles, [1, 2]):
        ax.set_facecolor('#faf9f6')
        for a, b in describe(mask)['edges_list']:
            ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color=color, linewidth=2.5)
        for x, y in pos.values():
            ax.add_patch(Circle((x, y), 0.075, facecolor='#faf9f6', edgecolor=color, linewidth=2.2, zorder=3))
        ax.set_title(title, fontsize=15, color='#222a30', pad=10)
        ax.text(0, -1.64, f'连通部分数：{count}', ha='center', fontsize=13, color=color)
        ax.set(xlim=(-1.9, 1.9), ylim=(-1.85, 1.55), aspect='equal')
        ax.axis('off')
    fig.suptitle('相同的已知性质与局部连接，可以对应不同整体', fontsize=17, y=.97, color='#222a30')
    fig.text(.5, .095, '两边均为 8 个顶点、8 条边；每个顶点连接 2 个邻居；没有三角形。',
             ha='center', fontsize=11, color='#414a50')
    fig.text(.5, .035, '此处圆点是图的顶点；上述各条性质才对应原始示意图中的“白点”。',
             ha='center', fontsize=10, color='#647078')
    fig.subplots_adjust(left=.035, right=.965, top=.84, bottom=.18, wspace=.10)
    fig.savefig(ROOT/'figures'/'local_same_global_different.png', dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot', action='store_true', help='also regenerate the figure')
    parser.add_argument('--font', help='CJK font path for the optional plot')
    args = parser.parse_args()
    main(make_plot=args.plot, font_path=args.font)
