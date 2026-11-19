# π 六点整数行列式研究

本目录保存六点细带路线的严格接口、模素数非零方法和一个从两条反正切级数重新构造的有限候选。

- [SIX_POINT_DETERMINANT_GATE.md](SIX_POINT_DETERMINANT_GATE.md)：秩五近似导致整数行列式为零的定量引理。
- [MODULAR_VANDERMONDE_CERTIFICATE.md](MODULAR_VANDERMONDE_CERTIFICATE.md)：把相位评价矩阵的非零性化为模素数 Vandermonde 因子与系数矩阵检查。
- [TWO_ARCTAN_SIX_POINT_CANDIDATE.md](TWO_ARCTAN_SIX_POINT_CANDIDATE.md)：由 arctan(2/5)+arctan(3/7)=π/4 构造六个整数截断点，得到 det=64、模 5 非零证书，并审计细带误差预算。
- [verify_two_arctan_six_point.py](verify_two_arctan_six_point.py)：只使用 Python 标准库和 Fraction，精确验证尾项比值包围与整数行列式。

**状态：** 六点候选及有限证书已完成；尚未得到随高度趋零的误差界，也没有证明 π 是坏逼近数。有限证书不能外推为无限族。

- [移位收敛分数六点引理](SHIFTED_CONVERGENT_SIX_POINT_LEMMA.md)：沿尾项斜率 \(\lambda=\log(2/5)/\log(3/7)\) 的连分数收敛分数构造六个整数点，使两组反正切余项比值同时趋于 1，并得到有界非零整数行列式 \(\det A=\pm4\kappa\)。这是一条严格的辅助几何构造，尚未给出 π 坏逼近结论。
- [移位六点精确复算](verify_shifted_six_point.py)：用有理交错余项区间和 Bareiss 消元复核 \(93/86\) 的有限样本。
