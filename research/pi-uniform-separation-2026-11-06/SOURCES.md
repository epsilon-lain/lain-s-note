# 来源、版本及实际阅读范围

检索/读取日期：2026-11-06。以下记录本轮外部原文，数学推导见 [DETERMINANT_AUDIT.md](DETERMINANT_AUDIT.md)。没有把搜索摘要当作已核验定理，也不声称读过列出范围之外的全部内容。

## 1. Zudilin：正矩与 Hankel 机制

- Wadim Zudilin, *A determinantal approach to irrationality*，arXiv:1507.05697v2，2016-02-02；期刊 Constr. Approx. 45(2) (2017), 301–310，DOI 10.1007/s00365-016-9333-7。
- [版本页](https://arxiv.org/abs/1507.05697v2)，[实际读取的 PDF](https://arxiv.org/pdf/1507.05697v2)。读取 §§1–2（PDF pp.1–3）及 §4（p.5），并读取 §5 的开头。
- 使用位置：Proposition 2 的条件为正矩表示、共同清分母序列 δn、整除链 δn|δn+1、指数增长及 εΔ^(3/2)/4<1。证明用 Heine/Vandermonde 正性及次数 n 的整数多项式，在假设 ξ=p/q 下付 q^n。Theorem 2 针对 π。
- 本轮只迁移正矩、正性及算术审查方法；不从无理性结论推导坏逼近性。我们的 μk 序列另有 rank1，不能照抄一般次数成本。
- 版本页显示 v1→v2，无撤稿提示。本轮另以题名加 erratum/correction 检索，未发现独立勘误；这是检索结果而非不存在勘误的证明。PDF 页眉 v2/2016-02-02 与首页正文的修订年份笔误须区分。

## 2. Fischler–Rivoal：指数极限的适用范围

- Stéphane Fischler and Tanguy Rivoal, *Irrationality exponent and rational approximations with prescribed growth*，arXiv:0910.4448v1，2009-10-23。
- [版本页](https://arxiv.org/abs/0910.4448v1)，[PDF](https://arxiv.org/pdf/0910.4448v1)，[作者站点 PDF](https://rivoal.perso.math.cnrs.fr/articles/iera.pdf)。读取 arXiv PDF pp.1–2 的 Lemma 1、Theorem 1 和 §4.1（p.5）及 §4.2 开头。
- Lemma 1 假定整数形式误差 n 次根趋于 0<α<1、系数 n 次根上极限≤β>1，推出 μ(ξ)≤1−logβ/logα。Theorem 1 在严格指数不等式下给逆向构造；§4.1 Proposition 1 仍允许 o(1) 指数。
- 用于辨认统一常数缺口：n 次根极限不能控制次指数因子，因此不能直接替代本轮所需的常数级窗口。本文不把该观察当作对任意构造的不可能定理。
- 版本页只有 v1；按题名加 correction/erratum 检索未见独立修正。此处不依赖版本页期刊页码元数据。

## 3. Zeilberger–Zudilin：必须使用修正后的算术引理

- Doron Zeilberger and Wadim Zudilin, *The irrationality measure of π is at most 7.103205334137…*，arXiv:1912.06345v2，2020-01-08；Moscow J. Comb. Number Th. 9 (2020), 407–419。
- [版本页](https://arxiv.org/abs/1912.06345v2)，[PDF](https://arxiv.org/pdf/1912.06345v2)，[作者页及代码链接](https://sites.math.rutgers.edu/~zeilberg/mamarim/mamarimhtml/pimeas.html)。实际读取引言、PDF p.2 修正说明、Part II 的 Proposition 1–2（PDF p.12）及其前后的整数化表述；未逐行重新证明 Part II 所有算术引理，也未运行 Maple。
- 原文明确说明 v1 的 Lemma 2 有误，v2 修复；因此本轮使用 v2 作为机制参照。作者页另有 2020-01-13 对 p.7 恒等式的机器证明链接。
- 这里只引用有限指数上界的构造范式，既不宣称它是截至本日的最新纪录，也不从其标题数字推出临界指数 2。

## 4. Legendre 基本关系的核对来源

- NIST Digital Library of Mathematical Functions，[§14.10，式14.10.3](https://dlmf.nist.gov/14.10.E3) 和 [§18.3，Table18.3.1](https://dlmf.nist.gov/18.3)。实际读取递推公式、正交多项式定义及表格说明。
- 对整数 n，令一般阶数 μ=0 的递推指标平移，即得本轮使用的 Legendre 三项递推。它是多项式恒等式，所以在 i 处代入不涉及一般非整数 Legendre 函数的分支。
- 本轮证明另从多项式差商积分导出 Wn 递推，从奇偶性和正交性导出全部矩根等式，不依赖不明版本的数值表。DLMF §18.3 标示 1.2.0 的说明文字修订，本轮未使用该新说明作独立定理。

## 5. 仅作检索线索、未认证的项目

检索出现 2026 年 π 无理性指数改进预印本及 OpenAI 数学语料的相关主张。本轮外部机制审查没有读取其全部证明、核验其算术引理或认证其结论，故不把它们当作已证前提。仓库已有语料版本和旧审查的实际状态，应在本轮 README 的继承检查点中单独记录。

本轮没有复制外部论文或其完整代码；这里保存的是阅读范围、条件及独立推导。

## 6. Minkowski 体积保证的精确版本

- Jan-Hendrik Evertse, *Notes on Diophantine Approximation*, Chapter 6, *Geometry of numbers*，PDF 首页日期2007-11-06；[作者托管的原讲义](https://pub.math.leidenuniv.nl/~evertsejh/Minkowski.pdf)。本轮实际读取 §6.1 的凸体定义、§6.2 Theorem 6.2.1、Lemma 6.2.2 及两者完整证明（印刷页46–48，PDF页1–3），另读 Corollary 6.2.3。
- 定理要求中心对称、紧、凸、内部含0，格满秩；体积至少 `2^n det(L)` 可得非零格点。本轮使用严格大于版本，(a,b,c)格的协体积为1，盒体积10>8；所有条件在 THEORY §2 逐项验证。讲义并未给出“格点在某个有理数处取值非零”，这正是本轮检验而失败的额外步骤。
- 本轮不依赖讲义 §6.3 的第二定理，也未声称重证其一般版本。该讲义未显示修订/勘误记录；本轮使用的第一定理另已在笔记给出标准自足证明梗概，数学成本是显式行列式1。

## 7. 有限范数菜单命题使用的经典超越性输入

- F. Lindemann, *Ueber die Zahl π*, Math. Ann.20 (1882), 213–225，DOI10.1007/BF01446522。[期刊书目](https://link.springer.com/article/10.1007/BF01446522)已读，但正文要求订阅；旧 Göttingen/Digizeitschriften PDF链接实际返回404。本轮**没有读取或重新认证1882年原文**。
- 为核对实际使用的定理，另读取 R. Steinberg and R. M. Redheffer, *Analytic proof of the Lindemann theorem*, Pacific J. Math.2 (1952),231–242，[期刊原文 PDF](https://msp.org/pjm/1952/2-2/pjm-v2-n2-p11-p.pdf)，DOI10.2140/pjm.1952.2.231。读取首页 Lindemann Theorem 的完整陈述、§3 Lemmas1–8、§4 Theorem1及证明、§5 Lemmas9–12与Theorem2及推导。纸本pp231–242；该期刊PDF另含封面/目录页。
- 使用条件：代数的互异指数的指数函数值在代数数上无非平凡线性关系。若 π 代数，互异的0和iπ皆代数，而 `1+exp(iπ)=0` 与该定理矛盾。由此仅调用 π 超越性于 COVER_AUDIT 命题1；**不会推出坏逼近性**。主文件的平方退化定理适用于任意实数，完全不依赖这个外部输入。
- 该1952原文将一般证明分为解析部分及代数对称化，特别强调不能由乘积等于零就跳过非零系数问题；本轮读到这一非零性处理，不把讲义或搜索摘要当作1882原文。本轮未重建全文的所有外部代数引理，调用的是经典 Lindemann 定理。题名检索未发现单独勘误，不能据此断言不存在后续修正。
