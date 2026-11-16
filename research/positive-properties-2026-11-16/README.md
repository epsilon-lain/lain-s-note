# 最少白点怎样识别整体：正性质的一个完整图模型

2026-11-16。来自本轮 Sol/Max 理论队的数学研究，由 Sol/High 外部观察员独立复核。执行面是 ChatGPT Work，原 Yohaku 宿主未在这里迁移。结论限定于明确对象库、正边性质和单位证据成本，不作文献首创声明。

**English abstract.** In a fixed family of labeled simple 2-regular graphs on 2m vertices with girth at least m, positive edge observations distinguish a single cycle from two m-cycles with minimum certificate size m. We classify all minimum certificates, separate isomorphism-type identification from labeled identification, and give consistent closure-exchange counterexamples. General proofs are paired with exhaustive finite checks. The results are scoped mathematical deductions, not claims of universal novelty or competition performance.

完整论证见 [THEORY.md](THEORY.md)。固定整数 `m>=3`，背景已知对象是一个 `C_(2m)` 或两个 `C_m`，每条已经确立的带标签正边性质成本1，未知边不算非边：

| 研究目标 | 本轮结论 |
|---|---|
| 识别同构类型 | 两种类型各恰需 `m` 条正边性质 |
| 单环的全部最小证书 | `m` 条连续环边；奇数 `m` 时另有两个完美匹配 |
| 两环的全部最小证书 | 完整观察其中一个 `m` 环 |
| 原八点例 | 同构类型识别最少4条，精确标号识别最少6条 |
| “不可再删”是否等于“最少” | 不等于：同一八点目标存在4条与5条的不可删证书 |
| 闭包能否直接套用线性独立／拟阵结论 | 不能：在非空候选集合上给出交换律失败的显式反例 |

精确标号的一般公式及三角形例外在完整证明中分开；不存在把有限枚举外推为一般定理的步骤。证明先将候选完成化为路径分量大小的子集和，再证明适用于所有参数的整数引理。程序完整核验 `m=3,4` 的候选图与观察子集，并在更大有限参数范围核对临界模式和引理；这些检查不代替一般证明。

在本目录复现，Python标准库即可：

```sh
python3 verify_white_points.py
python3 verify_white_points.py --check
```

结果见 [RESULTS.json](RESULTS.json)，独立复核范围见 [REVIEW_RECEIPT.json](REVIEW_RECEIPT.json)。`THEORY.md` 保留原始运行正文，其中 `theory/...` 命令和工作流路径相对于完整两队运行根目录；本目录应使用上面的命令。[SOURCE.json](SOURCE.json)给出四个逐字节副本的原始哈希。

[完整两队研究与工作流观察](https://github.com/epsilon-lain/yohaku/tree/codex/initial-release/research/theory-application-sol-2026-11-16)另存共享数学版本、方法变化、应用队证书与可恢复检查点。两队读取了对方的数学，但当前对象范围不同，没有把图识别定理强行当成 AC 搜索或13步结果的原因。
