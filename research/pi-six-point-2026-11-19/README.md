# π 六点整数行列式研究

本目录保存六点细带路线的严格接口、模素数非零方法和一个从两条反正切级数重新构造的有限候选。

- [SIX_POINT_DETERMINANT_GATE.md](SIX_POINT_DETERMINANT_GATE.md)：秩五近似导致整数行列式为零的定量引理。
- [MODULAR_VANDERMONDE_CERTIFICATE.md](MODULAR_VANDERMONDE_CERTIFICATE.md)：把相位评价矩阵的非零性化为模素数 Vandermonde 因子与系数矩阵检查。
- [TWO_ARCTAN_SIX_POINT_CANDIDATE.md](TWO_ARCTAN_SIX_POINT_CANDIDATE.md)：由 arctan(2/5)+arctan(3/7)=π/4 构造六个整数截断点，得到 det=64、模 5 非零证书，并审计细带误差预算。
- [verify_two_arctan_six_point.py](verify_two_arctan_six_point.py)：只使用 Python 标准库和 Fraction，精确验证尾项比值包围与整数行列式。

**状态：** 六点候选及有限证书已完成；尚未得到随高度趋零的误差界，也没有证明 π 是坏逼近数。有限证书不能外推为无限族。

- [移位收敛分数六点引理](SHIFTED_CONVERGENT_SIX_POINT_LEMMA.md)：沿尾项斜率 \(\lambda=\log(2/5)/\log(3/7)\) 的连分数收敛分数构造六个整数点，使两组反正切余项比值同时趋于 1，并得到有界非零整数行列式 \(\det A=\pm4\kappa\)。这是一条严格的辅助几何构造，尚未给出 π 坏逼近结论。
- [移位六点精确复算](verify_shifted_six_point.py)：用有理交错余项区间和 Bareiss 消元复核 \(93/86\) 的有限样本。

- [移位六点定量审计](SHIFTED_SIX_POINT_QUANTITATIVE_AUDIT.md)：把六个反正切尾项比值收敛加强为统一 \(O(q_j^{-1})\) 误差，并精确说明有界非零行列式为何尚未产生 π 的有理逼近矛盾。

- [精确尾项列行列式](EXACT_TAIL_COLUMN_DETERMINANT.md)：把 \(f_{m,n}=\pi/4-T_m(2/5)-T_n(3/7)\) 作为第六列，严格证明清分母后的行列式为整数，并给出 \(216LH^6\delta<1\) 的充分矛盾门槛；同时审计原始 Taylor 公分母成本。当前估计尚不足以满足门槛，也没有证明该路线不可能。
- [尾项列有限证书复算](verify_exact_tail_column.py)：以 Fraction 计算六个有理截断、精确清分母并用 Bareiss 与模 5 检查非零性。

- [奇偶兼容的移位六点族](PARITY_COMPATIBLE_SIX_POINT_FAMILY.md)：修复此前六点族未确保每个 \((m,n)\) 奇偶相反的缺口；利用相邻收敛分数在 \(\mathbb F_2^2\) 的基性质，构造无限多个 \(p\) 偶、\(q\) 奇的近似，再得到六个奇偶相反的整数点。给出尾项比值 \(1+O(q^{-1})\)、有界非零二次行列式，以及一个实际尾项列行列式模 5 非零的有限证书。
- [奇偶兼容尾项列复算](verify_parity_tail_column.py)：精确计算六个有理截断，清除公分母并以 Bareiss 和模 5 复核非零性。

- [BBP 接力与既有 Padé 工作分流](RESEARCH_HANDOFF_BBP_VS_PADE.md)：记录 11 月 19 日 Yohaku `pi-math` 稀疏例外接力的结论，区分 BBP–有限积分二进赋值路线与已经完成的 Niven-Padé、指数 Padé、OpenAI 插值行列式路线，避免重复研究。
