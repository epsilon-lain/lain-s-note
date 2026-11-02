# 本轮来源核对

2026-11-02。以下均为实际打开的主要来源；访问覆盖范围不等于全文审查或采用全部结论。

1. [NIST DLMF 4.25.4](https://dlmf.nist.gov/4.25#E4)：核对经典反正切的广义连分数。当前整数递推对应 x=1 的这族截断；重建这个族不构成新的 π 逼近族。
2. Peter Bala，2024，[Arctanh(z) and the Legendre polynomials](https://oeis.org/A003881/a003881.pdf)：核对 n!P_n(z) 作为分母的递推及 Legendre 表示。它帮助辨认经典结构；本轮 gcd 赋值结论有单独证明。
3. Tony D. Noe，2006，[On the Divisibility of Generalized Central Trinomial Coefficients](https://cs.uwaterloo.ca/journals/JIS/VOL9/Noe/noe35.pdf)：读取生成函数、递推、Legendre 对应、Schur／Holt 同余和素数整除章节。本轮 ell 的 Lucas 数位分解属于这里讨论的经典 Schur 型现象，未声称该同余首创。
4. Frits Beukers，2000，[A rational approach to π](https://www.nieuwarchief.nl/serie5/pdf/naw5-2000-01-4-372.pdf)：核对第 377–378 页的同一个正积分、首批分数及连分数联系。该文还介绍另一类具有更好清分母估计的积分；本轮没有把后者静默换成当前固定族。
5. Shalosh B. Ekhad 与 Doron Zeilberger，2014，[Searching for Apéry-Style Miracles](https://arxiv.org/pdf/1405.4445)：读取积分与反正切构造部分。它展示解析小值、整数清分母和经典质量估计的配合。此处不以其中给出的质量指数数值，代替对实际约分后稀疏子序列的共同约数控制。

OEIS [A012244](https://oeis.org/A012244) 和 [A006139](https://oeis.org/A006139) 用于寻找这些原始文献及核对初值。A012244 是本轮 A_n；A006139 是整数 U_n=2^n ell_n，其中心系数表达式直接证明了筛法所需的整性。证明不依赖 OEIS 表格外推。

本轮没有进行 π 当前最佳无理性测度的全面文献检索，也没有将用户原稿的全部高维结论视为独立认证结果。新增的反例、赋值证明及必要筛选条件以本目录中的完整论证和可运行证书为准。
