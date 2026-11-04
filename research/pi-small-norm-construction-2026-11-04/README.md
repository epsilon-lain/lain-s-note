# π² 小范数构造与预置公因子：有限尝试

2026-11-04。目标是构造分母无界、且 `Q=q²|π-p/q|→0` 的分数。

**本轮没有得到新的 π 逼近纪录，也没有得到无限成功族。** 实际尝试了一个移动二次型构造，并加入可主动安排公因子的整数权重。它确实降低了辅助 π² 分数的算术成本，但尚未产生随规模延续的小范数解。

## 实际做了什么

令两列数 `u_n,v_n` 满足经典 Apéry 递推

\[
(n+1)^2X_{n+1}=(11n^2+11n+3)X_n+n^2X_{n-1},
\]

初值分别为 `(1,3)`、`(0,5)`。相邻的 `6v_n/u_n` 严格夹住 π²。将这些源或它们的正整数加权组合约成 `A/B`，再寻找整数解

\[
Bp^2-Aq^2=h.
\]

这次同时收费：既限制 `h`，也限制 `q²` 乘以 π² 源的误差上界。具体充分条件、完整推导及搜索覆盖范围见 [THEORY.md](THEORY.md)。

| 试验 | 参数或形式数 | 检查的根式收敛分数 | 输出的不同分数 |
|---|---:|---:|---|
| 原源及相邻权重 `k=0,...,8` | 1152 个二次型 | 82075 | `3/1,22/7,355/113` |
| 按模核预置公因子，再选择正权重 | 128 个二次型 | 12231 | `3/1,22/7,355/113` |

两组都取 `1≤n≤128`，各测试 `|h|≤min(25,B)` 和 `|h|≤floor(sqrt(B))` 两种预算。共 94306 个根式收敛分数；重复命中并不计作新分数。该计数不包括任意更大分母、其他范数预算或其他构造。

## 有用的正向变化：先安排约分

设相邻的两个 π² 源已经约成 `A0/B0,A1/B1`。以

\[
G=|A_0B_1-A_1B_0|
\]

为模数，求一个整数 `k`，使两坐标同时满足

\[
A_0+kA_1\equiv B_0+kB_1\equiv0\pmod G.
\]

然后从这个模核中选择非负整数权重，形成新的源。权重的大小、完整公因子和最终既约分母均写入证书；没有把乘上权重的成本藏掉。

第 128 项的实测结果：

| 同一 π² 源区间的收费项 | 原源 | 加权并完全约分后 |
|---|---:|---:|
| 既约分母位数 | 238 | 134 |
| `B² × 区间宽度` | 约 `2.544×10^208` | 约 `2.365` |

这是**辅助 π² 源的有限高度改善**。它使用经典有理连分数／二维整数格机制；不构成新的 π 定理，也没有证明这些权重全局最优。把源压到这个尺度后，仍需要小范数解；两件事必须一起成立。

一个精确的有限种子是

\[
1112\cdot355^2-10975\cdot113^2=25.
\]

它来自第 4 项，并通过了余项预算。但测试中没有找到能让这种种子随 `n` 无界延续的规则。

## 后续构造的具体对象

应把正权重与小范数解**一起设计**。若 `g` 是实际公因子，要同时满足

\[
\begin{aligned}
xA_0+yA_1&=gA,\\
xB_0+yB_1&=gB,\\
x(B_0p^2-A_0q^2)+y(B_1p^2-A_1q^2)&=gh,
\end{aligned}
\]

并守住 `q² × 源误差` 的预算。当前实现先压缩源的分母，再求小范数解；尚未构造出使两步同时成功的无限族。这是待攻的构造问题，不是已经证明的存在定理。

## 覆盖树程序的定位

附带完成了“整数 arctan 源 → 严格有理区间 → 整段 Farey 分支 → 覆盖命中”的程序。它从精确整数递推起算，不输入 π 小数或连分数表，认证了 3988 段分支。重现的最佳分数分母为 216 位，质量严格落在

\[
0.000048130561246289404628
<Q<
0.000048130561246289404629.
\]

其后继部分商是 20776。**这是已知结果的校验样例，不能作为研究突破。** [OEIS A033089](https://oeis.org/A033089) 早已列出该旧纪录，并记载 Max Frank 于 2025-11-01 报告的 653520000000 个连分数项计算；本轮没有下载或复核那个大型数据集。继续扩大这里的小型扫描不是当前构造目标。

## 文件与复现

- [conic_attempt.py](conic_attempt.py)、[conic_results.json](conic_results.json)：1152 个原始／相邻混合二次型。
- [content_lift.py](content_lift.py)、[content_lift_results.json](content_lift_results.json)：128 个预置公因子的源及二次型。
- [verify_conic.py](verify_conic.py)、[conic_verification.json](conic_verification.json)：二项式和、Casoratian 求和、独立根式区间算法复核。
- [verify_content.py](verify_content.py)、[content_verification.json](content_verification.json)：整数权重、公因子、既约分母与全部有限解预算的独立复核。
- [construct.py](construct.py)、[results.json](results.json)、[deep_scan.json](deep_scan.json)：经典覆盖树校验程序与有限轨迹。
- [verify.py](verify.py)、[verification.json](verification.json)：反向广义连分数、树边界、另一种 Taylor 误差公式交叉检查。

仅使用 Python 3 标准库。在本目录运行：

```sh
python3 conic_attempt.py
python3 content_lift.py
python3 verify_conic.py
python3 verify_content.py
```

可选的覆盖程序校验：

```sh
python3 construct.py
python3 construct.py --targets 1/300 1/1000 1/10000 1/100000 1/1000000 --max-order 2048 --output deep_scan.json
python3 verify.py
```

三份验证报告均为 `PASS`，其中覆盖验证含短、长两份数据。这里的“独立”指不同公式与实现，不是新增外部审稿或多代理复核。耗时字段因机器而变；数学数据与整数不等式才是证书。

## 参考来源

1. L. Jain and P. Tzermias, [Beukers’ integrals and Apéry’s recurrences](https://cs.uwaterloo.ca/journals/JIS/VOL8/Tzermias/tzermias5.pdf), *Journal of Integer Sequences* 8 (2005), Article 05.1.1：π² 源、递推与积分的经典依据。
2. B. Lynn, [Continued Fractions — Convergence](https://crypto.stanford.edu/pbc/notes/contfrac/converge.html)：Legendre 判据。
3. NIST DLMF, [§4.25, equation 4.25.4](https://dlmf.nist.gov/4.25.E4)：arctan 的经典广义连分数。
4. [OEIS A033089](https://oeis.org/A033089)、[A001203](https://oeis.org/A001203)：既有数值计算的规模与有限旧纪录；不用于候选生成。
