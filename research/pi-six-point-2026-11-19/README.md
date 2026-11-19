# π 六点整数行列式研究

本目录保存六点细带路线的严格接口、模素数非零方法和一个从两条反正切级数重新构造的有限候选。

- [SIX_POINT_DETERMINANT_GATE.md](SIX_POINT_DETERMINANT_GATE.md)：秩五近似导致整数行列式为零的定量引理。
- [MODULAR_VANDERMONDE_CERTIFICATE.md](MODULAR_VANDERMONDE_CERTIFICATE.md)：把相位评价矩阵的非零性化为模素数 Vandermonde 因子与系数矩阵检查。
- [TWO_ARCTAN_SIX_POINT_CANDIDATE.md](TWO_ARCTAN_SIX_POINT_CANDIDATE.md)：由 arctan(2/5)+arctan(3/7)=π/4 构造六个整数截断点，得到 det=64、模 5 非零证书，并审计细带误差预算。
- [verify_two_arctan_six_point.py](verify_two_arctan_six_point.py)：只使用 Python 标准库和 Fraction，精确验证尾项比值包围与整数行列式。

**状态：** 六点候选及有限证书已完成；尚未得到随高度趋零的误差界，也没有证明 π 是坏逼近数。有限证书不能外推为无限族。
