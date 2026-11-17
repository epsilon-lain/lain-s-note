# 证明重构：目标结构、完整费用与片段组合

2026-11-17。已实际完成云端理论/应用两方向辅助研究；原生 Yohaku 旧阶段未改动。本目录保存数学推导与有限核验，完整实验、实现修正、独立证书重放和协作记录见[过程包](https://github.com/epsilon-lain/yohaku/blob/de48e5389ae0a322345127bc816cd95df3a24528/research/proof-reconstruction-round2-2026-11-17/README.zh.md)。

## 数学内容

- [理论全文](theory/README.zh.md)：若所需共轭词为 h v^±1 h^-1，利用现成关系词 v，给出 2|h|+4 个官方原语的实现上界。证明 donor 不变、执行顺序和 rank 条件；不声称这个上界总是最短。
- 同文给出任意精确片段替换后的完整 length/work/peak 公式，并证明固定候选图上的 Pareto 动态规划可以选出预算内最优组合。有限证明图中，单次最大收益的贪心给3步，兼容组合给2步；不是官方计分实例。
- [目标结构短注](GOAL-STRUCTURE-NOTE.zh.md)：在实际训练证书的精确入口 (y,x)，原“先整理为(x,y)再删除”的7步可直接改成2步。更一般地，正置换单字母关系词组到 Stable 空目标，rank 为 k 时最短恰 k 步；下界来自每步最多把 rank 减1。

这些是标准群恒等式、费用计数和有限图方法的自含推导，不主张数学首创、IUT 推论或一般 AC 最短性。

## 已实际执行的结果

理论作者核验112条 rank2…8 合成编译路线，并检查一份真正到 Stable 空目标的合成证书：35→31步，work521→437，peak31不变。[代码](theory/check_theory.py) · [完整核验](theory/checks.json)。结果记录实际执行脚本哈希；当前脚本只调整输出格式，沿革见[方法记录](theory/METHOD-LOG.md)。

应用固定21例，普通清理与增加整体词库两臂的最终动作列表完全相同：各将总步数345降为310，七个公开训练实例各省5步，全部来自同一个已有单原语替换机制。新增整体词库的额外收益为0；真实计分题 sac-07445 仍15步/work128。训练例不计比赛分数。

主代理用独立动作解释器复核63条原/最终证书，并通过官方 CLI 入口验收7条训练改善与旧真实题路线；未独立重跑优化搜索。[完整独立回执](https://github.com/epsilon-lain/yohaku/blob/de48e5389ae0a322345127bc816cd95df3a24528/research/proof-reconstruction-round2-2026-11-17/INDEPENDENT-REVIEW.json) · [CLI回执](https://github.com/epsilon-lain/yohaku/blob/de48e5389ae0a322345127bc816cd95df3a24528/research/proof-reconstruction-round2-2026-11-17/CLI-ACCEPTANCE.json) · [原始应用结果](https://github.com/epsilon-lain/yohaku/blob/de48e5389ae0a322345127bc816cd95df3a24528/research/proof-reconstruction-round2-2026-11-17/application/results.json)。

七例改善都是绕过同一多余中间目标，不能归功于未命中的新宏，也不能当作七种独立机制。固定候选图的最优组合算法尚未用于本轮应用比较。没有平台上传、排行榜查询或速度/token收益声明。

[固定协议](PROTOCOL.json) · [完整可复现代码与过程](https://github.com/epsilon-lain/yohaku/blob/de48e5389ae0a322345127bc816cd95df3a24528/research/proof-reconstruction-round2-2026-11-17/README.zh.md)。
