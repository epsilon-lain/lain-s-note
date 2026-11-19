# 六点行列式的模素数非零证书：Vandermonde 降维

日期：2026-11-19  
状态：**严格代数引理**。它给 exact phase / invariant 输出一个可直接检查的模素数非零接口；是否适用于当前六点矩阵，取决于该矩阵能否在模某个素数下化为六个相位点的多项式评价矩阵。

## 1. 评价矩阵引理

令 \(\ell\) 为素数，令 \(x_1,\ldots,x_6\in\mathbb F_\ell\) 两两不同。对 \(j=0,\ldots,5\)，取多项式 \(f_j(X)\in\mathbb F_\ell[X]\)，次数至多 5。定义

\[
A_{ij}=f_{j-1}(x_i),\qquad 1\le i\le6,\quad 1\le j\le6.
\]

将六个多项式按幂基展开：

\[
f_j(X)=\sum_{r=0}^5 c_{rj}X^r,
\qquad C=(c_{rj})_{0\le r,j\le5}.
\]

令 \(V=(x_i^r)_{1\le i\le6,\,0\le r\le5}\)。则 \(A=VC\)，所以

\[
\boxed{
\det A
=\det C\prod_{1\le i<j\le6}(x_j-x_i).
}
\tag{1}
\]

因此只要 \(\det C\ne0\)，便有 \(\det A\ne0\) 于 \(\mathbb F_\ell\)。若 \(A\) 是某个整数矩阵 \(\widetilde A\) 的模 \(\ell\) 约化，则

\[
\boxed{\det\widetilde A\not\equiv0\pmod\ell
\quad\Longrightarrow\quad
\det\widetilde A\ne0\text{ 于 }\mathbb Z.}
\tag{2}
\]

### 证明

按定义，
\[
(VC)_{ij}=\sum_{r=0}^5x_i^r c_{r,j-1}=f_{j-1}(x_i)=A_{ij}.
\]
故 \(A=VC\)。取行列式并使用 Vandermonde 公式即得 (1)。若 \(\det A\ne0\) 模 \(\ell\)，整数行列式不可能为零，得到 (2)。证毕。

## 2. 对六点相位构造意味着什么

假设 exact phase 能为六个截断点给出模 \(\ell\) 的相位坐标
\[
x_1,\ldots,x_6\in\mathbb F_\ell
\]
且它们两两不同；又假设六点整数矩阵 \(A\) 的列，在模 \(\ell\) 下是六个次数至多 5 的相位多项式 \(f_0,\ldots,f_5\) 的取值。那么非零性不需要直接展开一个巨大的 \(6\times6\) 行列式，只需证明：

1. 六个相位剩余 \(x_i\) 两两不同；
2. 系数矩阵 \(C\) 可逆。

这两项分别是六个一次同余检查和一个固定大小的模素数行列式检查。若六个 \(f_j\) 恰好构成 \(\mathbb F_\ell[X]_{\le5}\) 的一组基，则 \(\det C\ne0\) 自动成立。

这提供了一个具体的“exact phase \(\Rightarrow\) modular nonvanishing”目标，而不再只是笼统地说“找个素数让行列式不为零”。

## 3. 与细带上界拼接

若同一个整数矩阵 \(A\) 还满足前一份笔记中的秩五近似条件，则

\[
|\det A|\le1296\varepsilon(H+\varepsilon)^5.
\]

若右侧严格小于 1，则 \(\det A=0\)。另一方面，若第 1 节的相位评价结构在某个素数 \(\ell\) 下成立，则 \(\det A\ne0\)。两者矛盾。

因此完整的六点爆破解可以被拆成两个互不循环的证书：

- **解析证书**：给出秩至多 5 的矩阵 \(B\)，以及 \(H,\varepsilon\)，使 \(1296\varepsilon(H+\varepsilon)^5<1\)。
- **算术证书**：给出素数 \(\ell\)、六个两两不同的相位剩余、系数矩阵 \(C\)，并验证 \(\det C\not\equiv0\pmod\ell\)。

目前第二个证书仍需从当前实际六点构造中提取；本引理本身不声称已经找到这样的 \(\ell\) 或 \(x_i\)。

## 4. 必须防止的伪证书

- 六个整数点彼此不同，不代表六个相位模 \(\ell\) 后仍然不同。
- 六个相位两两不同，只能保证 Vandermonde 部分非零；若 \(\det C=0\)，整个行列式仍然为零。
- 只计算浮点行列式或浮点相位不能证明模素数非零。
- 若清分母时乘过的整数因子含有 \(\ell\)，必须在清分母前后明确追踪该因子；不可把一个被乘成零的模数证书当成原始矩阵的证书。

## 5. 下一步的实际代入表

取得真实六点矩阵后，按以下顺序填写：

| 对象 | 必须提供的精确数据 |
|---|---|
| 六个相位 | \(x_1,\ldots,x_6\pmod\ell\) |
| 相位分离 | 六个差 \(x_j-x_i\not\equiv0\pmod\ell\) |
| 六个列函数 | \(f_0,\ldots,f_5\pmod\ell\) |
| 基变换 | 系数矩阵 \(C\) |
| 非零检查 | \(\det C\pmod\ell\) |
| 原矩阵关系 | \(A\equiv VC\pmod\ell\) |
| 解析端 | \(H,\varepsilon\) 及严格不等式 \(1296\varepsilon(H+\varepsilon)^5<1\) |

只有所有行都由精确整数/有理运算验证后，才可宣称六点行列式矛盾闭合。
