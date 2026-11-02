# 输入来源、核查范围与一处补充限定

2026-11-02。输入是用户提供的 `pi_next_work_20261102.zip`，归档 SHA-256 为 `45f0f7c529c2002b46214b6158821bdded09995771690f1de5ba47e7efddd5f5`。`imported-next-work/` 保存其全部原始文件；原 `RESEARCH_MANIFEST.json` 中的 19 项文件字节数与 SHA-256 已逐项匹配。没有把输入自带的独立审查当成本轮审查的替代。

## 本轮重新核查的数学内容

- 近根的全部混合端点系数、Laguerre 阶乘积分、全阶几何上界、移动 Z 时的一致性，以及 Rouché 定理和对称性带来的唯一简单虚轴根。
- 完整反射因子的常数项、符号及系数高度／Mahler 测度主阶。有限因子的系数由次数上界加足够多的精确行列式值认证；其模素数证书由另一套标准库有限域程序认证。
- 指数 Padé 对照的归一化、余项、周期小值、最近根尺度、完整有理求值分母及阈值。有限数值根图例没有用来证明渐近结果。
- 本轮新建的反正切／正积分整数线性形式、真正约分后的分母、误差与 gcd 的渐近等价、邻近整除约束，以及固定指标差的质量乘积界。

各自的假设及原仓库输入仍见对应证明和 `reviews/`。原稿的完整高维插值定理及 μ(π)=2 的宣称不是本轮联合算术论证的输入。

## 需要补在输入中的范围

输入 `literature/PADE_LINEAR_DEGREE_BENCHMARK.md` 的式 (13) 直接使用式 (3) 的非零有理求值下界，因而该推导应明确限制为 **b≥2**，等价于既约 p/q 的 q≥3。q=1、2 应另行处理。这个有限分母限定不影响文中 m/b 渐近阈值或大分母结论。为保留输入字节指纹，原文件没有改写；本说明与独立审查同时记录该限定。

## 实际打开核对的主要来源

1. Kathy A. Driver 与 Nico M. Temme，*Zero and pole distribution of diagonal Padé approximants to the exponential function*，1999。[作者机构归档全文](https://ir.cwi.nl/pub/1229/1229D.pdf)。核对其引言中的显式指数 Padé 系数、积分表示及余项，公式 (1.2)–(1.5)；其规范化不同，输入中的常数另由精确幂级数计算核查。
2. Doron Zeilberger 与 Wadim Zudilin，*The irrationality measure of π is at most 7.103205334137…*，2019／2020。[作者全文](https://sites.math.rutgers.edu/~zeilberg/mamarim/mamarimPDF/pimeas.pdf)。读取方法说明及 Part II 的积分、算术清分母与额外共同因子结构。它作为“联合核算线性形式与共同约数”的历史对照；此处不把标题数值称为当前最佳界。
3. George Csordas、Marios Charalambides 与 Fabian Waleffe，*A new property of a class of Jacobi polynomials*。[作者全文](https://people.math.wisc.edu/~fwaleffe/submitted/Jacobi_CCW2.pdf)。核对文内 Hermite–Biehler 表述及 Corollary 3 的 Bessel 稳定性，用于输入 Padé 对照的可选全根几何说明。
4. John Cullinan 与 Nick Scheel，*On the arithmetic of Padé approximants to the exponential function*。[预印本页面](https://arxiv.org/abs/2007.01329)。核对作者、标题、摘要及研究范围。没有把其完整 Padé 多项式的不可约性移用于奇部或本项目反射因子。

原输入还包含其他文献和 2026 年论文现状的介绍。这些现状段落作为输入原文保留，未被本轮采作已认证事实或证明依赖；本轮没有进行当前最佳无理性测度的全面文献核验。

## 关于复核与有限计算

独立复核是其他研究代理重新推导并运行检查，不是外部同行评审。全 n 结论来自解析证明及精确恒等式；有限 n 表格不用于推断一般不可约性、未来 gcd 增长或 π 的坏逼近性。新增六个标准库复现入口见 README，实际执行记录及校验和见 CHECKS.json。
