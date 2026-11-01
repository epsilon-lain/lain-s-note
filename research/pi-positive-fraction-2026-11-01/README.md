# 正比例逃逸、周期拼接反例与完整外幂和

2026-11-01。延续 [上一轮外幂研究](../pi-exterior-blast-2026-11-01/README.md)。这轮把少数行的逃逸推进到正比例行，并检验局部非零能否保证完整周期独立性。

**没有证明 π 是坏逼近数，也没有构造新的 π 有理逼近纪录。** 以下结果是运输矩阵的渐近定理、完整有限包的精确证书和一个中间命题的反例。它们不等于原稿的端点矛盾。坏逼近要求一个对所有分母统一的正下界；无理测度为 2 本身并不提供这个常数。

## 1. 正比例逃逸确实有成本，但高维仍会分摊它

设源行满足

\[
v_0s+\sum_{i=1}^m v_i\beta_i<H,\qquad
c_i=v_i-v_0>0,\qquad
\mathcal A=\prod_i\frac{v_i}{c_i}.
\]

每次运输保留中心和总阶数 N=s+|beta|。共有 M_H 行，若其中逃出原源包的比例趋于 rho，所有目标选择都满足

\[
\boxed{\liminf_{H\to\infty}
\frac{c\cdot\Gamma}{M_HH}
\ge\frac{\rho^2}{(m+1)\mathcal A}.}
\]

这里 Gamma 是全部目标类型与源类型之差。固定维数、权重后先取 H 极限；不能把这些参数一起变大而沿用误差项。

我们同时构造了完整包、非零 Pascal 运输系数的例子。在 H 极限之后，再使 c_*/v_0 增大，逃逸比例可趋于 rho，而归一化位移成本趋于

\[
\frac{2\rho^2}{m+2}.
\]

因此“必须移动很多行”仍不能恢复维度无关的平方惩罚。成本约为 1/m 的分摊可以发生在非零运输系数上。**这还没有证明对应的周期子式非零或足够小。** 详见 [POSITIVE_FRACTION_LAW.md](POSITIVE_FRACTION_LAW.md)。

## 2. 局部打包非零，不能自动推出整体周期独立

找到了一个精确反例：完整源包有 10 行，每个局部 Pascal 行列式都非零，低阶行全部保留，目标中也没有零行，完整运输系数为 zeta^6。但目标周期矩阵的秩只有 8。

这个结论对每个非零形式周期 Z 都成立：源子式为 -9216 Z^5，目标有两个明确的多项式行关系，同时一个 8 阶子式为 16 Z^4。因此它也适用于实际周期 2*pi*i。

原因可以直接看出来：列预算中某个横向类型只剩常数与 e^t 两种时间形状，它们的一阶、二阶导数提供不了两个独立方向。**“目标类型存在”与“它有足够的独立时间方向”是不同要求。**

这是局部拼接推理的反例，不是 OpenAI 原稿的反例，也不排除别的目标选择成功。详见 [PERIOD_GLUING_OBSTRUCTION.md](PERIOD_GLUING_OBSTRUCTION.md)。

## 3. 也找到了完整的大批运输成功例子

另一个完整包包含 156 行、242 个候选原始列。保留全部低总阶数行，在 N=8,...,19 的整层使用 a=(2 beta_1+j,beta_2)：66 行位于这些层，51 个行标签改变，33 行逃出源包。

| 精确检查 | 结果 |
|---|---:|
| 完整运输系数，普通 Taylor 归一化 | 2^60 (-zeta_1)^93 |
| 同一组 156 列上，目标周期子式在 Z=1 模 65537 的值 | 9333 |
| 同一组列上，有理源子式模 65537 的值 | 44706 |

表中两个模值使用阶乘缩放行；归一化因子已完整记录。两套独立实现得到同一组选列和同一结果。目标周期多项式非零，再用已知超越性，得到其在实际周期的值非零。

这是一个完整有限包的成功拼接证书，尚不是渐近族。它只证明一个目标集合对应的项非零；同次的其他项仍可能与它抵消。详见 [BULK_PERIOD_CERTIFICATE.md](BULK_PERIOD_CERTIFICATE.md)。

## 4. 完整余项不能只按低阶的小秩处理

固定误差次数 D 时，余项只作用在 O(D H^m) 条边界行上，可以用行列式引理精确压缩。但完整余项未必如此：一维情形下，带真实列类型预算的运输边界算子满足

\[
\frac{\operatorname{rank}R(H)}{M_H}
\longrightarrow\frac{v_1-v_0}{v_1+v_0}.
\]

这个比例随权重比增加可以接近 1。它是运输系数矩阵的秩，不自动等于乘上周期矩阵后的秩。详见 [BOUNDARY_RANK_FILTRATION.md](BOUNDARY_RANK_FILTRATION.md)。

## 5. 下一步所需的估计已经明确

在每个固定参数选择下，所有可用目标行只有 O(H^(m+1)) 条。完整外幂求和的项数，取对数后除以 M_H H，会趋于零。因而用绝对值逐项控制的路线可以归结为

\[
\frac{\log|\det(TB)|}{M_HH}
\le \frac{\log2}{w_*}
+\max_J\left\{
\frac{\log|\det B_J|}{M_HH}
-\frac{\tau\rho_J^2}{(m+1)\mathcal A}\right\}+o(1),
\]

其中 rho_J 是逃逸比例，tau=min_i(-log|zeta_i|/c_i)。这个上界已经包含整个目标集合求和；它不提供下界，也不声称最大项不会被抵消。原源包的恒等项还可以单独处理。

当前缺口是 **实际周期子式随逃逸比例变化的统一增长上界**，尤其必须覆盖逃逸比例趋于零的选择。证明目标子式非零不能替代这个量级估计。若走利用相消的路线，则需要完整求和的相消定理。详见 [EXTERIOR_SUM_ENVELOPE.md](EXTERIOR_SUM_ENVELOPE.md)。

## 复现与审查

在本目录运行以下命令；全部只用 Python 标准库。

```sh
python verify_positive_fraction.py
python verify_rank_filtration.py
python verify_bulk_transport.py
python verify_bulk_review.py
python verify_period_obstruction.py
```

精确检查覆盖 175 个成本选择、12 个非零纤维行列式、12 个完整包预算实例、864 个 Cauchy 秩块、348 个边界支持检查、6 个截断行列式恒等式、两套 156 行证书，以及拼接反例的完整多项式恒等式。无限极限结论来自文中的证明，不由有限实验外推。

独立复核记录为 [INDEPENDENT_REVIEW.md](INDEPENDENT_REVIEW.md) 与 [INDEPENDENT_RANK_REVIEW.md](INDEPENDENT_RANK_REVIEW.md)；最终执行记录与文件哈希见 [CHECKS.json](CHECKS.json)。

研究背景是用户提供的 OpenAI 手稿 *The irrationality exponent of pi is 2*（2026-09-24；输入文件 SHA256：`9da33449f4d2aed948401fa2e3aa9e01e8bc753708891b37710762cb56c6d1cd`）。本轮没有重新独立验证其 Theorem 2.1，也没有把它升级为无条件的坏逼近结论。本轮的有限秩证书均直接计算，不依赖该定理。
