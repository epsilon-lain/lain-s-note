# 从实际初值筛到可证书化搜索：本轮推进与下一目标

日期：2026-11-02。基于仓库检查点 `05e6e140c6a41a5ed62f38d6d583487d0455bc7c` 的多区间亏损定理。

**这轮没有构造出质量趋零的子序列。** 实质推进是：利用实际初值，严格缩小可提供额外约分的素数集合；把成功指标间隔定量化；实现可提前拒绝指标的精确搜索；另提出一个不必逐指标控制异常聚集的充分目标。

固定记号为
\[
X_{n+1}=(2n+1)X_n+n^2X_{n-1},\quad
(A_0,A_1)=(1,1),\quad(B_0,B_1)=(0,1),
\]
\[
a_n=A_n/n!,\quad b_n=B_n/n!,\quad
g_n=\gcd(A_n,4B_n),\quad
p_n=4B_n/g_n,\quad q_n=A_n/g_n,
\]
\[
Q_n=q_n^2|\pi-p_n/q_n|,\quad U_n=2^na_n,\quad C_n=\operatorname{odd}(U_n).
\]
对奇素数令 \(s_p(n)=v_p(g_n)-v_p(n!)\)，正、负部分乘积分别为 \(E_n,D_n\)，\(J_n=\log(E_n/D_n)\)。

## 1. 真正用到初值的筛：高次整除可能完全无益

已有离散 Wronskian 与分母界是
\[
n(a_nb_{n-1}-a_{n-1}b_n)=(-1)^n,\qquad
a_j\in\mathbb Z[1/2],\quad
\operatorname{den}(b_j)\mid\operatorname{lcm}(1,\ldots,j)
\quad\text{于奇素数部分}.
\]
令
\[
m_p=v_p(a_n)=v_p(C_n),\qquad
k_p=\lfloor\log_p(n-1)\rfloor-v_p(n).
\]
对所有 \(n\ge2\) 和奇素数 \(p\)，若 \(m_p>k_p\)，Wronskian 第一项严格被 \(p\) 整除，故第二项必须为单位。因此
\[
\boxed{m_p>k_p\implies
s_p(n)=v_p(b_n)=-v_p(n)-v_p(a_{n-1})\le0.}
\]
特别地，正超额必须满足
\[
\boxed{s_p(n)>0\implies1\le v_p(C_n)\le k_p.}
\]
这比仅控制 \(\min(v_p(a_n),v_p(b_n))\) 更强，包含 \(p\mid n\) 和素数幂指标；不需要任何扩域。

因此可把旧的正超额倍数完整剥除这些无益素数：
\[
\boxed{E_n\mid G_n,\qquad
G_n=\prod_{\substack{p<n\text{ 奇}\\1\le m_p\le k_p}}p^{m_p}.}
\]
例如 \(n=68\) 时，\(v_{37}(C_{68})=2>k_{37}=1\)，旧倍数中的整个 \(37\) 应删除，而不是把它当作可能的收益。固定亏损证书时，质量下界提高 \(37^2\) 倍。

本轮也复核了 Legendre 微分 Bézout 恒等式，但它在 \(i\) 处等价于上述 Wronskian，**不把它宣传成另一条独立新恒等式**。新内容是更强的赋值用途。完整证明见 [initial-data/ACTUAL_INITIAL_BEZOUT_FILTER.md](initial-data/ACTUAL_INITIAL_BEZOUT_FILTER.md)。

## 2. 从首层异常到“恰一次异常”的新必要条件

保留前轮所选区间 \(t_n<p\le n\)，其中 \(t_n\asymp n^{2/3}\)。定义
\[
W(n)=\sum_{\substack{t_n<p\le n\\p\mid C_n}}\log p,
\qquad
W_1(n)=\sum_{\substack{t_n<p\le n\\v_p(C_n)=1}}\log p.
\]
选定大素数若提供正超额，必须在 \(C_n\) 中**恰出现一次**。结合前轮亏损证书，得到
\[
\boxed{J_n\le W(n)+W_1(n)-n+o(n),}
\]
\[
\boxed{\log Q_n\ge(2+\log2)n-2W(n)-2W_1(n)+o(n).}
\]
因此，任何固定 \(M>0\) 的 \(Q_n\le M\) 子序列都要求
\[
W(n)+W_1(n)\ge(1+\tfrac12\log2)n-o(n).
\]
又因 \(2W-W_1\le\log C_n\) 且已有
\(\log C_n=\log(2+\sqrt2)n+O(\log n)\)，得到
\[
\boxed{W_1(n)\ge
\left(\frac{2+\log2-\log(2+\sqrt2)}3-o(1)\right)n
=(0.4884000010868\ldots-o(1))n.}
\]
这是对“恰一次出现的大素数”的必要规模，不能把它误称为实际正超额的质量下界。它仍然没有排除聚集。

## 3. 成功指标至少间隔约 \(0.3466n/\log n\)

对任意固定 \(M>0\)，两个充分大的 \(Q_n,Q_m\le M\) 指标 \(m>n\) 满足
\[
\boxed{m-n\ge\left(\frac{\log2}{2}-o(1)\right)\frac n{\log n}.}
\]
所以 \([1,N]\) 内满足 \(Q_n\le M\) 的指标数至多为
\[
\boxed{\left(\frac1{\log2}+o(1)\right)(\log N)^2.}
\]
证明来自有限整数整除证书。设 \(P_n\) 是 \(W(n)\) 对应的平方自由素数乘积，\(B=\lfloor n^{1/3}\rfloor\)，\(h=m-n\)。递推转移给出正整数多项式 \(G_h(n)\)，使
\[
P_nP_m\mid
\left(\prod_{\substack{p\le m\\p\text{ 奇}}}p\right)
\left(\prod_{b=1}^{B}U_b\right)G_h(n),
\qquad0<G_h(n)\le(5m)^h.
\]
这一步不需要在模素数意义下除以任何递推分母。共同大素数质量于是至多为 \(h\log(5m)+O(n^{2/3})\)，和前轮单指标必要质量结合即可。

完整证明见 [spacing/quantitative-spacing.md](spacing/quantitative-spacing.md)；独立审计见 [spacing-audit/INDEPENDENT_AUDIT.md](spacing-audit/INDEPENDENT_AUDIT.md)。这加强了“间隔趋于无穷”，仍允许无限稀疏子序列；渐近间隔没有有效起点，不能据此跳过有限搜索区间。

## 4. 已落地的精确搜索

直接逐项构造阶乘尺度的 \(A_n,B_n\) 很浪费。筛选器只递推 \(O(n)\) 位的整数
\[
(n+1)U_{n+1}=2(2n+1)U_n+4nU_{n-1},\quad U_0=1,\ U_1=2.
\]
先计算并剥除 \(G_n\)，再按需加入可证明的亏损乘积 \(\widehat D_n\mid D_n\)。令 \(d_n=v_2(U_n)-n\)，已有严格误差下界给出
\[
\boxed{Q_n>
\frac{2}{2n+1}\,2^{-2d_n}\frac{\widehat D_n^2}{G_n^2}.}
\]
右侧与任意正有理目标做精确整数比较；一旦超过目标便停止处理该指标。未被排除的指标只是候选，才需要完整 \(B_n\) 和误差计算。无需分解整个 \(C_n\)。

本轮标准库版本已精确证明
\[
\boxed{Q_n>1\qquad(4\le n\le50000).}
\]
这一有限结论不是数值近似：主筛不构造完整 \(B_n\)，累计只做 8 次额外亏损素数测试；每项最多 2 次，只请求了 \(B_1=1\)。新初值筛在 298 个指标上完整删除了 300 个无益素数幂。另对 \(n\le10000\) 独立构造原始 \(A_n,B_n\)，逐项核对了精确 gcd 上界，全部通过。

随后用可选的系统 GMP 精确 gcd 后端把证书延伸到 **100000**，仍证明所有 \(4\le n\le100000\) 的 \(Q_n>1\)，仅留下指标 1、2、3。后端本身通过 506 组精确 gcd 对照，并与标准库在整个 \(n\le10000\) 范围逐字段比较数学结果；该扩展记录单独保存在 [search/adaptive_sieve_100000_gmp.json](search/adaptive_sieve_100000_gmp.json)。它使用同一整数判据，不依赖浮点 π。这个后端目前限于常规 64 位 Linux/GMP ABI；标准库版本不依赖它。

在本目录运行：

```bash
python search/adaptive_exact_sieve.py --N 10000 --audit --output audit-rerun.json
python search/adaptive_exact_sieve.py --N 50000 --output sieve-rerun.json
# 可选后端；先运行其独立对照检查：
python search/verify_gmp_backend.py
python search/adaptive_exact_sieve.py --N 100000 --gmp --output sieve-gmp-rerun.json
```

公式、边界、操作计数和精确结果见 [search/ADAPTIVE_EXACT_SIEVE.md](search/ADAPTIVE_EXACT_SIEVE.md) 与 [search/adaptive_sieve_50000.json](search/adaptive_sieve_50000.json)。性能统计仅描述这次运行，不是复杂度保证。

## 5. 下一刀：先限制“可用素数池”的容量

这条路线改变要证明的对象。利用广义中央三项式身份
\[
U_n=[x^n](1+2x+2x^2)^n,
\]
定义一个与指标 \(n\) 无关的素数集合
\[
\mathcal Z=\{p\text{ 奇素数}:\exists\,0\le r<p,\quad U_r\equiv0\pmod p\},
\qquad\vartheta_{\mathcal Z}(x)=\sum_{\substack{p\le x\\p\in\mathcal Z}}\log p.
\]
Lucas 数字乘法说明：若 \(p\mid U_n\)，必有 \(p\in\mathcal Z\)。所以对**所有指标**直接有
\[
W(n)\le\vartheta_{\mathcal Z}(n).
\]
由前轮定理立即得到一个充分研究目标：如果能对某个固定 \(\eta>0\) 证明
\[
\boxed{\vartheta_{\mathcal Z}(x)\le
\left(\frac{2+\log2}{4}-\eta\right)x
\quad\text{对所有充分大的 }x,}
\]
那么 \(\log Q_n\ge4\eta n+o(n)\)，本族所有质量有界的稀疏子序列都被排除。这样不必再逐指标证明上百上千个异常无法同步。

精确分类所有 \(p\le20000\) 的 2261 个奇素数，915 个属于 \(\mathcal Z\)。

| \(x\) | \(\vartheta_{\mathcal Z}(x)/x\) |
| ---: | ---: |
| 5000 | 0.40716859 |
| 10000 | 0.40738777 |
| 20000 | 0.39969607 |

这些结果有助于选择理论目标，**不构成渐近密度定理**。模 8、16 的每个奇剩余类都观察到大量活跃与不活跃素数，不能凭简单剩余类猜测宣布完成。文献中的固定素数数字规律，也不能当成跨素数的密度上界。

优先推进这一素数池上界及其平均零点数版本；搜索则用精确筛保留真正未排除的候选，不继续把单个高次整除当作突破。若目标是构造新的 \(Q_n\to0\) 子序列，这轮的搜索仍没有给出候选；上述素数池路线若成功，结论是当前族失败，届时需要改变构造来继续研究 π。

主文献为 Noe (2006) 的广义 Lucas/Holt 同余和 Kohen (2026) 的数字乘法表述。它们支持序列识别与同余，**不提供这里所需的跨素数上界**。详细来源与扫描见 [active-primes/ACTIVE_PRIME_CAPACITY.md](active-primes/ACTIVE_PRIME_CAPACITY.md)。

## 6. 验证范围和仍未完成的部分

- 初值赋值定理：精确多项式恒等式 0–25 次、全序列至 1000，含所有奇素数及素数幂边界；另有独立审阅。
- 间隔证书：原推导 4000 组转移/高度核验；独立程序 12000 组整除证书核验。
- 搜索器：原始整数全审计至 10000；标准库精确证书至 50000；通过后端独立对照的 GMP 扩展证书至 100000。
- 素数池：分类至 20000，另直接核对小素数的完整 Holt 同余 11066 次。

尚未证明上述素数池密度目标、尚未证明 \(S(n)\le n/6\)、尚未排除本族全部无限稀疏成功子序列，也未证明 π 的坏逼近性或非坏逼近性。有限验证、渐近必要条件、待证充分目标在本记录中分别标明。
