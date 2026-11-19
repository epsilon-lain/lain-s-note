# 六点细带的行列式门槛：一个可直接调用的严格引理

日期：2026-11-19  
状态：**已证明的线性代数引理**；尚未代入最新 Padé/相位构造的具体六点矩阵。  
用途：把“六个整数点落在同一条极窄细带”转成精确的行列式估计，并明确指出非零性必须独立证明。

## 1. 精确接口

令 \(A\in M_6(\mathbb Z)\)。假设存在一个实矩阵 \(B\in M_6(\mathbb R)\)，满足

\[
\operatorname{rank}(B)\le 5,\qquad
\|B\|_{\max}\le H,\qquad
\|A-B\|_{\max}\le\varepsilon,
\]

其中
\[
\|M\|_{\max}:=\max_{1\le i,j\le6}|M_{ij}|,
\qquad H,\varepsilon\ge0.
\]

则

\[
\boxed{
|\det A|\le1296\,\varepsilon(H+\varepsilon)^5.
}
\tag{1}
\]

### 证明

记 \(A_j,B_j\) 为 \(A,B\) 的第 \(j\) 列。由于 \(\operatorname{rank}(B)\le5\)，有 \(\det B=0\)。逐列替换给出望远镜恒等式

\[
\det A-\det B
=\sum_{j=1}^6
\det(A_1,\ldots,A_{j-1},A_j-B_j,B_{j+1},\ldots,B_6).
\]

每个行列式中：

- 特殊列 \(A_j-B_j\) 的 Euclidean 范数至多 \(\sqrt6\,\varepsilon\)；
- 其余五列的每个坐标绝对值至多 \(H+\varepsilon\)，故每列 Euclidean 范数至多 \(\sqrt6(H+\varepsilon)\)。

由 Hadamard 不等式，每一项绝对值至多

\[
(\sqrt6\,\varepsilon)
\bigl(\sqrt6(H+\varepsilon)\bigr)^5
=216\,\varepsilon(H+\varepsilon)^5.
\]

求和六项即得 (1)。证毕。

## 2. 整数性给出的结论，以及它没有给出的结论

因为 \(A\) 是整数矩阵，

\[
\det A\in\mathbb Z.
\]

若

\[
1296\,\varepsilon(H+\varepsilon)^5<1,
\tag{2}
\]

则 (1) 只能推出

\[
\boxed{\det A=0.}
\tag{3}
\]

**注意：这本身不是矛盾。** 要得到爆破解，还必须从独立的算术不变量证明

\[
\boxed{\det A\ne0.}
\tag{4}
\]

最直接的证书是找到一个素数 \(\ell\)，并精确证明

\[
\det(A\bmod\ell)\ne0\quad\text{于 }\mathbb F_\ell.
\]

这会立即推出 (4)。也可以用精确整数消元、非零子式或其他独立方法，但不能从“六点互异”或“六点处于细带”直接推断行列式非零。

因此六点路线的逻辑闭环必须是：

\[
\underbrace{\det A\in\mathbb Z}_{\text{清分母后的整数矩阵}}
,\qquad
\underbrace{\det A\ne0}_{\text{独立非零证书}}
,\qquad
\underbrace{|\det A|<1}_{\text{细带误差界}},
\]

三者同时成立才产生矛盾。

## 3. 细带宽度的定量预算

若 \(0<\varepsilon\le1\) 且 \(H\ge1\)，则 \(H+\varepsilon\le2H\)，所以

\[
|\det A|\le 1296\cdot32\,\varepsilon H^5
=41472\,\varepsilon H^5.
\]

一个简单的充分条件是

\[
\boxed{\varepsilon<\frac1{41472H^5}.}
\tag{5}
\]

也就是说，若六点矩阵的坐标规模为 \(H\)，则要仅靠一个秩至多五的近似矩阵把整数行列式压到 1 以下，细带误差至少需要达到 \(H^{-5}\) 量级（常数如上）。若矩阵条目实际具有不同尺度，应先逐行/逐列做**精确有理归一化**，再用对应的行列式缩放因子恢复整数性；不能把归一化后的非整数行列式直接称为整数。

## 4. 对当前两组反正切路线的可证伪检查

一旦 Padé endpoint、exact phase、invariant 和 certificate phase 的产物可读，必须逐项填入：

1. 六个整数点如何组成 \(A\in M_6(\mathbb Z)\)；给出完整矩阵，不能只给行列式数值。
2. 由 exact phase 得到的秩五矩阵 \(B\) 是什么；明确写出其五维行/列空间或精确线性关系。
3. 六点每个坐标与 \(B\) 的偏差，给出统一有理上界 \(\varepsilon\)。
4. \(H=\|B\|_{\max}\) 的精确上界。
5. 用 (1) 或更强的专用估计核验 \(|\det A|<1\)。
6. 用独立的模素数秩证书或精确消元核验 \(\det A\ne0\)。

**红队测试：** 对同一矩阵构造，分别替换成有理目标、黄金比例和一个无理性指数为 2 但非坏逼近的稀疏连分数目标。若非零证书对所有目标都声称成立，而细带上界又只依赖一般余项估计，则必须检查是否遗漏了目标特有的相位假设；否则很可能把一般不可能成立的结论藏进了“非零性”步骤。

## 5. 当前研究状态

- 本文件中的行列式扰动引理及其证明：**已完成**。
- 当前六点构造的具体 \(A,B,H,\varepsilon\)：本轮可见材料中没有完整条目，故尚不能代入。
- 当前 exact phase / invariant / certificate phase 的实际输出：只知道状态标签为 finished，尚没有可核验的矩阵、模数、余数或有理区间内容。
- 因此，本文件不宣称已证明 \(\det A\ne0\)、\(|\det A|<1\)，也不宣称已经证明 \(\pi\) 是坏逼近数。

下一步不是继续抽象讨论“六点”，而是把最新六点矩阵原样接入第 1 节，先独立验证模素数非零，再把 Padé 余项界代入第 3 节的 \(H^{-5}\) 预算。若任一环节失败，保留精确反例并修改构造。
