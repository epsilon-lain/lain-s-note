# 找回的历史双波与 Machin 相位线索

2026-11-03 整理。两份来源都是历史聊天导出，不是本轮找到的 Lean 源码；不复制对话中的个人信息。以下独立核对了列出的代数恒等式与圆周稠密性推论，没有审计原对话的其他“已证”标签。

## 1. 双波的精确代数部分

来自 [2026-09-15 历史材料](https://chatgpt.com/library?file_id=libfile_1aeb4f289c84819192aff601e92bc560&account_id=personal) 中约第 2800–3100 行。原聊天数学格式部分损坏；尚未确认它是否就是旧恢复清单中要求的带 \(q,s\asymp H\) 量词的完整 Twin-Wave 版本。

取整数向量 \((a,b),(A,B)\)，假设
\[
aB-bA=1,\qquad b\ne0,\qquad e=a+\pi b\ne0.
\]
定义
\[
\alpha=\frac{A+\pi B}{e},\qquad \beta=\frac{B}{b},\qquad M=be.
\]
直接由行列式得到 \(\beta-\alpha=1/M\)。对整数 \(h,k\)，设
\[
(P,Q)=h(A,B)-k(a,b),\qquad \theta=h\alpha-k.
\]
则
\[
\boxed{Q(P+\pi Q)=M\theta\left(\theta+\frac{h}{M}\right).}
\]
证明仅用 \(P+\pi Q=e\theta\) 与 \(Q=b(\theta+h/M)\)。这是真实的两因子恒等式，但两因子互相关联，不能只凭它们各自可能小就宣布同时命中。

如果另外构造一列满足 \(|Q|\to\infty\)、两因子乘积与 \(|M|\) 的乘积趋零，则对应有理数 \(-P/Q\) 的归一化误差趋零；约分只会把该质量再除以 gcd 的平方。原聊天中 Fejér 捕获条件、最近整数一致、有效参数范围与原始 Twin-Wave 全部量词尚需恢复和审校，不能从本代数恒等式直接宣布解决 P2848。

## 2. Machin 双链的固定宽度相位窗口

来自 [2026-09-01 历史材料](https://chatgpt.com/library?file_id=libfile_9ccea9a08a5c8191aaf5ce786f3aa505&account_id=personal) 的相位讨论。原稿中具体增长率、误差渐近与无理性作为前置条件；本次未独立核验原稿的整条 L1–L8 链，不能把这些前提当作已证。

仅保留可独立验证的普遍事实：对于无理 \(\theta\)，圆周旋转 \(x\mapsto x+\theta\pmod1\) 对任何非空开区间的返回具有统一有界间隔。证明是：所有逆向平移的目标区间构成圆周开覆盖，再用紧性选取有限子覆盖。

这与[当前双 arctan 构造](../coupled-arctan-2026-11-03/JOINT_DIVISIBILITY.md)中在长窗口内寻找相位对齐的机制相连。但固定相位窗口带来的是固定程度的相消；索引返回有界本身不控制已约分分母、其对数尺度间距，或临界 \(q^2\) 质量。若要求误差相对比例趋零，还要让窗口缩小并另计返回时间。

## 3. 查找边界

旧 Twin-Wave 的完整量词与尺度约束、plateau routing 原始定义、PellLambert/RiccatiJet Lean 源码和编译日志，本次仍未完整恢复。历史聊天中有关证明完成或形式化通过的描述不能代替这些材料；本记录也不声称发现 π 的新逼近纪录。
