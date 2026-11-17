# 来源与独立审阅范围

日期：2026-11-17。本文是成熟方法的建模、收费与重现，不作数学优先权声明。

## 一手资料

1. Gilbert Strang，MIT OpenCourseWare：[Lecture 22: Diagonalization and Powers of A](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-22-diagonalization-and-powers-of-a/)。对角化用来简化矩阵幂；本文自己推导2×2形状的变换和精确费用，没有挪用任意矩阵均可同样处理的结论。
2. J. M. Pollard，[The Fast Fourier Transform in a Finite Field](https://luca-giuzzi.unibs.it/corsi/Support/papers-cryptography/2004932.pdf)，Mathematics of Computation 25(114), 365–374 (1971)，DOI 10.1090/S0025-5718-1971-0301966-0。上述链接是大学站点保存的原论文副本。本文借鉴有限域傅里叶变换与卷积的经典关系；具体常数、程序及收费是本轮自行实现与推导。AMS原PDF获取返回403，没有把它记成已读。
3. J. W. Cooley 与 J. W. Tukey，[原论文的IBM条目](https://research.ibm.com/publications/an-algorithm-for-the-machine-calculation-of-complex-fourier-series)及[大学保存的原论文PDF](https://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Cooley65.pdf)，An Algorithm for the Machine Calculation of Complex Fourier Series (1965)。本文使用二进制蝶形框架，按自己的标量费用计价。
4. SymPy [Discrete官方文档](https://docs.sympy.org/latest/modules/discrete.html)对NTT、逆NTT、模数与卷积提供实现参考。本轮检查器不依赖SymPy，不把其内部自动填零的线性卷积接口误当成本文的n项循环卷积。

没有重新发布完整论文或长段原文，也未调用IUT定理。

## 本轮实际执行

- 主代理实现 `check_matrix.py` 并实际运行，穷举F17全部83,521输入四元组；每组分别验证结果和精确收费。
- 辅助任务 `/root/net_gain_fourier_audit`（请求模型Sol/high）独立推导并实现NTT完整账本，实际运行583对输入，原脚本与结果未经主代理改写复制到ntt目录。n=16含全部256个标准基对，其他规模为结构例子和固定种子样本。
- `/root/net_gain_matrix_audit`（Sol/high）独立解析推导矩阵例子，并在主代理加计常数2构造后静态审阅源码／收费，未再次执行枚举。变化矩阵链的4m+7／6m公式获解析审阅，没有另作大规模数值实验。
- `/root/net_gain_information_audit`（Sol/high）静态审阅NTT源码／结果，确认单位根、递归根索引、逆归一化和账本；未独立重跑。审阅特别要求区分场内算术与内存／学习等未计费用。
- 主代理阅读两种NTT算法及结果并核对代数推导；没有把同一作者的多种测试称作第二份独立程序复算。

## 未采用的候选

也考察了前缀和表示。它适合说明付费转换把稠密读出变成稀疏读出，但其逐查询原始读取下界需要禁止派生缓存；允许跨查询复用后，某些滑动区间批次可用滚动和。为避免把表示限制误称全体算法下界，它没有被拿来作为本文主要数值优势证据。

## 原生状态

此轮只在当前云端会话研究与复现。未接续WSL原负责人，未增加原生模型turn，未修改shared_stop、旧阶段、预算或截止，未运行比赛搜索或平台提交。既有原生成果和失败均保留。
