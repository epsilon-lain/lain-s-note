# OpenAI 数学语料接入与第二轮 P2848 工具分析

2026-11-01。上游：[`openai/math`](https://github.com/openai/math)，固定提交 [`4eaee87548250b2682affda7d05d947c05554bae`](https://github.com/openai/math/commit/4eaee87548250b2682affda7d05d947c05554bae)。

## 接入内容与数量核对

当前官方 `CONTENTS.md` 和 `overview.tex` 均为 **721 篇稿件、373 个研究主题组**，不是对话中提到的 772 篇。逐条解析得到 721 个不同 PDF 路径；同版本 `preprints` 树共有 751 个目录，其中 30 个未列入当前目录。未列入不自动表示撤稿或无效；本次没有逐个推断其状态。

| 内容 | 入口 | 本次状态 |
|---|---|---|
| 全部上游源码、PDF、Lean 与历史材料 | 仓库根目录 [external/openai-math](../../external/openai-math) 子模块 | Gitlink 固定到上述提交；包含上游整个仓库，而非仅本轮精选论文 |
| 721 篇可浏览论文索引 | [CATALOG.md](CATALOG.md) | 每篇都有固定版本 PDF 链接、主题组与学科 |
| 标题、摘要、分类与来源元数据 | [catalog.json](catalog.json)、[families.json](families.json) | 可离线搜索；原文结果描述是上游主张，不代表独立验证 |
| 30 个目录外文件夹 | [uncatalogued.json](uncatalogued.json) | 保留精确目录名与上游 tree SHA，不混入 721 篇计数 |
| 快照与许可 | [snapshots](snapshots)、[manifest.json](manifest.json) | 官方目录、总览、形式化清单、历史、Apache-2.0 许可与本地 SHA-256 |
| 第二轮定向分析与新工具 | 下文及三个分支报告 | 全量元数据检索，7 篇重点稿件的指定定理／证明段落阅读；不是 721 篇证明的逐行审计 |

形式化信息保留两个不同字段：175 篇稿件的路径列在上游 `formalization.yaml` 的 sources 中；245 个主题组的目录说明带有 Lean 范围页链接。它们不等于本次运行并通过的 Lean 证明数，也不能直接互相换算为形式化完成率。本轮没有运行这些论文的 Lean。

子模块接入避免复制整套大型 PDF，同时保存精确版本关系。第一次在本地克隆时使用：

```bash
git clone --recurse-submodules https://github.com/epsilon-lain/lain-s-note.git
```

已有工作副本则在仓库根目录运行：

```bash
git submodule update --init external/openai-math
```

这一步实际下载完整上游内容。此提交已把 GitHub 上的来源关联和固定版本保存好；本轮工作环境没有下载全部 PDF 或逐篇读取全文。只运行下面的元数据搜索不需要下载子模块。

## 立即检索

Python 3 标准库即可搜索标题、摘要、主题说明与分类：

```bash
python3 research/openai-math-corpus-2026-11-01/corpus.py search 'matroid'
python3 research/openai-math-corpus-2026-11-01/corpus.py search --family 017
python3 research/openai-math-corpus-2026-11-01/corpus.py search --tag determinant_rank --limit 30
```

`--json` 输出结构化结果。查询按空格分词并要求每个词都出现，是确定性的文本检索，不是语义模型；例如可以搜索 `identity`，或直接查 family 116。

子模块下载后，安装 ripgrep，即可检索全部 TeX、Markdown、Lean 源文：

```bash
python3 research/openai-math-corpus-2026-11-01/corpus.py fulltext 'Wronskian'
```

此命令先核对本地子模块 HEAD 与清单版本一致，再运行 `rg`。它搜索源文本，不把 PDF 二进制当文本。无 TeX 的其他材料仍通过完整子模块与 PDF 链接访问。

从保存的官方快照重建索引：

```bash
python3 research/openai-math-corpus-2026-11-01/corpus.py build
python3 research/openai-math-corpus-2026-11-01/verify_corpus.py
```

没有设置定时更新或自动移动上游版本。今后升级必须同时更新 gitlink、快照、固定 SHA、索引与分析来源，避免把旧证明笔记混成新版本。

## 第二轮分析的范围

对全部 721 篇标题、摘要、主题说明进行可复现关键词标记：丢番图方向 56 篇、行列式／秩方向 65 篇、拟阵／优化方向 32 篇、动力系统方向 14 篇、算术求和方向 28 篇。类别重叠；关键词命中既不等于相关定理成立，也不保证对本项目有用。不同学科对“rank”等词含义不同，后续按对象与量词筛除误匹配。

定向阅读以下 **7 篇不同稿件**，含主定理或相关证明段落：005、022、108、113、114、116、153。三个分支之间的交叉阅读只计一次。017 的 π 论文已在本项目此前审读，本轮未把它重复计作新阅读。精确阅读范围、未复算部分和固定版本原文链接见：

- [数论与动力系统分析](analysis_number_theory.md)：005、022、153；另筛读数论 001–031 与动力系统 143–154 的主题概述。
- [代数与消失阶分析](analysis_algebra.md)：116、108；与组合分支交叉核对 114。
- [组合与算法分析](analysis_combinatorics.md)：114、113；34 个关键词候选主题组的筛选与对象匹配。

## 真正接到当前缺口的工具

当前目标仍是：对**同一个**非零行列式，控制列成本、最低非零展开代价、清分母成本和解析上界。我们尚未证明 π 坏逼近或非坏逼近。

| 优先级 | 来源 | 能提供的具体工具 | 仍缺的连接 |
|---|---|---|---|
| 最高 | 116 微分系统与有限零检测 | 为同类周期展开行设置精确容量；计算充分 Taylor 深度 | 容量在 H 增长时是否持续提供足够收益 |
| 高 | 114 两拟阵共同基 | 每个原始行选一项且向量独立的有限组合模型；可探索计数／抽样 | 当前平方占用成本非加性；计数不控制带符号相消 |
| 高 | 005 Catalan 行列式 | 固定模素数块证书推出增长矩阵族非零；算子分解避免逆 Vandermonde | 尚未为我们的实际 π 矩阵构造相同块结构与范数估计 |
| 中 | 153 有界系数紧性 | 反证法得到受控系数族的统一有限阶非消失 | 行列式系数高度随 H 增长；常数项归一化也有成本 |
| 辅助 | 113 匹配熵 | 正权匹配和／permanent 上界 | 对 signed determinant 没有相应下界；仅改进 M! 计数属低阶 |
| 暂缓 | 108 行列式表示复杂度 | 一般参数几何、极向交点工具 | 表示所需矩阵大小不等于指定点的值或算术小量 |
| 不可直接套用 | 022 Duffin–Schaeffer | 几乎处处逼近与重合分类 | 几乎处处不能指定到 π |

### 已经落成的新约束

固定选列 P=YʰX^α 和展开指标 a。可用列的周期函数都是
\[
\binom\alpha a e^{hz}z^{|\alpha|-|a|}.
\]
按函数类型 \((h,|\alpha|-|a|)\) 去重，得到数量 rₐ。若某个展开组合用了 nₐ 个这一类行，那么非零行列式必须满足
\[
\boxed{n_a\le r_a.}
\]
理由是这些行经过的函数空间维数至多 rₐ；超过该维数必定线性相关。完整证明独立写在代数分析中，不需要相信 family 116 的整个主定理。

在已有四列 `1, X₂, X₂², X₁` 中，不移动时只有三种函数类型；四行超过容量 3，故这一展开组合严格为零。这将上一轮的精确秩现象变成能直接加入占用优化的必要约束。

另一工具给出有限 Taylor 深度：对每个 h 补齐次数 0,…,dₕ 后，取 D=Σₕ(dₕ+1)，前 D 个系数足以判定该指数多项式空间的零性。证明用 D 阶常系数微分方程，附合流 Vandermonde 公式。

[capacity_tool.py](capacity_tool.py) 与 [capacity_results.json](capacity_results.json) 保存实现和六组精确验证，包括独立求导取系数、符号秩和行列式公式。它们要求 SymPy，本次使用 1.14.0。容量只是必要条件，不能替代不同 a 之间的完整独立性检查。

## 本轮后的研究顺序

1. 对同一组已认证非零的列计算 rₐ，给现有占用优化加入 nₐ≤rₐ，同时保留二次成本 cΣnₐ²。
2. 用完整符号秩校验筛过容量的候选，比较 H 增长时最低存活代价与原松弛的差。
3. 若差不能达到所需尺度，重点尝试 005 的固定块非零证书和抗节点合流的算子估计；每次完整追踪分母。

共同基数量多、有限样例通过、某篇标题声称大定理，都不能直接替代这个量化比较。此次尚未闭合 P2848 的反证。

## 来源与许可

上游内容由 OpenAI 发布，许可证为 Apache-2.0；本目录快照保持来源标识，许可原文见 [snapshots/LICENSE](snapshots/LICENSE)。生成的索引与本轮分析是对此固定快照的整理与研究注记。上游仓库说明其内容处于不同核验阶段，视作研究数据；本接入没有赋予所有稿件新的独立正确性认证。
