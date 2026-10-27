# π 的变进制重复数位与有理逼近：研究记录

日期：2026-10-27。结论性质：初等恒等式、带证明的等价改写、反例与有限精确实验；没有证明 P2848，也未声称所列改写具有文献上的新颖性。未运行 Lean。

## 1. 结果概览

原始构想是把分母 q 放到 b=q+1 进制中，因为 n/q=0.\overline n_b。研究结果：

1. 可以仅根据前三个数位，对距离给出一个严格区间；乘以 q² 后，区间宽度小于 1/q。
2. P2848 等价于：存在进制 b 趋于无穷，使 π 的前两位小数相同，而且第三位与它们之差除以 b 趋于零。第三位不必完全相同。
3. “前两位相同”本身不够；sqrt(2)-1 有无穷多个这样的进制，却是坏逼近数。
4. “前三位完全相同无穷多次”过强，会产生三次幂精度，不能直接当作 P2848 的等价目标。甚至“第三位与前两位仅差一个固定常数”仍然过强。
5. 精确扫描 4≤b≤1,000,000：10 个进制前两位相同，没有前三位全同的进制。这只是有限结论。

## 2. 对象与约分

令 x=π−3∈(0,1)。令 b≥4 为整数进制，q=b−1。

用贪心展开定义数位 d₁,d₂,d₃,…：

    x = d₁/b + d₂/b² + d₃/b³ + …,
    0≤dⱼ≤b−1.

每个 dⱼ 是一个数位的数值，不是它的十进制字符串。例如 114 进制中，数位 16 是一个位置的值。x 无理，因此不存在有限展开与尾随最大数位之间的二义性。

令 n=d₁，则 0≤n≤q，且

    R=n/q=0.\overline n_b.

两个端点 n=0,q 分别表示 0,1。对固定 x∈(0,1)，当 b 足够大时 n/q→x，故最终远离端点。

逼近 π 的数是 3+n/q=(3q+n)/q。它的约简后分母为 q/gcd(n,q)，可能小于 q。本记录用 q² 归一化误差；若约分，约简后的归一化误差只会更小。沿误差趋于零的无穷序列，约简后分母也必趋于无穷，否则无理数 x 会成为一组有限有理数的极限。

## 3. 有限数位的精确保守距离界

令 K≥1，A_K 是前 K 位在 b 进制下组成的整数：

    A_K=d₁ b^(K−1)+…+d_K.

令 t_K=b^K x−A_K∈[0,1)，这是所有剩余数位的总贡献在平移缩放后的值。令

    M_K=A_K−n(b^K−1)/(b−1).

M_K 是整数，表示实际前缀与 K 个 n 构成的前缀之差。定义有符号归一化误差

    E=q²(x−n/q).

直接代数化简得到精确恒等式

    E=(q²/b^K)(M_K+t_K−n/q).                 (1)

因此，仅知道有限位，也有严格双侧区间

    (q²/b^K)(M_K−n/q)
       ≤ E <
    (q²/b^K)(M_K+1−n/q).                    (2)

绝对误差的上下界可以从这个区间直接计算：若区间跨过零，则下界为零；否则下界取端点绝对值的较小值，上界取较大值。

不使用 n/q 细节时，还能简化为

    (q²/b^K) max(0,|M_K|−1)
       ≤ |E| <
    (q²/b^K)(|M_K|+1).                      (3)

区间 (2) 的宽度恰为 q²/b^K。K=2 时宽度趋于 1；K=3 时宽度为 q²/(q+1)³<1/q，趋于零。这解释了为什么三个数位足以越来越准确地判定 q² 尺度上的误差，同时不要求三个数位相同。

## 4. 前 K 位相同的准确条件

给定数位 n∈{0,…,b−1}，记 R=n/q。前 K 位全为 n 当且仅当

    −R/b^K ≤ x−R < (1−R)/b^K.              (4)

证明：K 个 n 的有限和为 R(1−b^(−K))，具有这个前缀的数恰位于从该有限和开始、宽度 b^(−K) 的半开区间。减去 R 即得。

这同时准确处理了区间边界、借位以及正负误差的不对称性。对固定 x，好的 R 会趋于 x，故这两个容许宽度的比例最终接近 x 与 1−x，而不是退化成零。

## 5. P2848 的三数位等价表述

定理。对任意固定无理数 x∈(0,1)，以下两条等价：

(A) 存在 qⱼ→∞ 和整数 nⱼ，使

    qⱼ² |x−nⱼ/qⱼ| → 0.

(B) 存在 bⱼ→∞，使 x 在 bⱼ 进制下满足

    d₁=d₂=n，且 |d₃−n|/bⱼ → 0.

这里 n 随进制变化，不要求不同进制的重复数位相同。

证明 B⇒A：当 d₁=d₂=n 时，M₃=d₃−n。由 (1)，

    E=(q²/b³)(d₃−n+t₃−n/q).

因 |t₃−n/q|≤1，有

    ||E|−(q/b)² |d₃−n|/b| ≤ q²/b³.         (5)

右端趋于零，q/b→1，故 B 推出 A。还可得到便于使用的界

    ||E|−|d₃−n|/b| < 3/b.                 (6)

证明 A⇒B：令 δ=x−n/q，E=q²δ→0，b=q+1，R=n/q→x。对足够大的 q，有 0<R<1。计算

    bx=n+R+bδ,
    b²x=n(b+1)+R+b²δ.

因为 bδ→0、b²δ→0，而 R→x∈(0,1)，这两个式子的尾部最终都在 (0,1) 中。因此前两位确为 n。此后应用 (5)，可知 |d₃−n|/b→0。

应用于 x=π−3，即得 P2848 的等价数位条件。

量词必须保留：等价地说，对每个 η>0 和每个进制下限 B₀，都要存在 b≥B₀，使前两位相同且 |d₃−d₁|<ηb。不能用某个固定 η 或有限样本替代。

## 6. 为什么对应原始的 sin 问题

P2848 的原始形式是 liminf_{m→∞} m|sin m|=0，其中 m 为正整数。

若 q²|π−a/q|→0，令 D=|a−πq|，则 a≈πq，且

    a|sin a|≤aD≈πqD→0.

反过来，若 m|sin m|→0，取最近的 π 的整数倍 πq，令 D=|m−πq|≤π/2。因 sin D≥2D/π，得

    qD≤(m/π+1/2)(π/2)|sin m|→0.

这正是 q²|π−m/q|→0。足够好的逼近最终位于 (3,4)，可写成 3+n/q。

## 7. 两位完全相同仍然不足：严格对照

取 x=sqrt(2)−1。对任何 n/q∈[0,1]，

    q²|x−n/q|
      = |2q²−(q+n)²|/(sqrt(2)+1+n/q)
      > 1/4.                              (7)

分子是非零整数，分母小于 4。所以 x 不满足 A。

另一方面，取 Pell 解 A_k+Q_k sqrt(2)=(1+sqrt(2))^k，并令 N_k=A_k−Q_k。对于 k 足够大，0<N_k<Q_k，且 A_k²−2Q_k²=(−1)^k。于是

    Q_k²(x−N_k/Q_k)=−(−1)^k/(sqrt(2)+A_k/Q_k).

绝对值趋于 1/(2sqrt(2))≈0.353553。重复两位的条件 (4) 在 K=2 下，两侧极限界分别是 −x≈−0.414214 与 1−x≈0.585786；上述正负极限都严格落在里面。因此所有充分大的 Pell 分母 Q_k，在进制 Q_k+1 下，前两位都相同。

这证明：无穷多次两位相同，并不推出 P2848 类型的性质。第三位的相对偏差是必要的附加控制。

## 8. 三位完全相同为什么过强

若 d₁=d₂=d₃=n，则 |x−n/q|<b^(−3)<q^(−3)。无穷多次这样的事件，要求无理性指数至少为 3；P2848 只要求二次幂误差前的常数可以趋于零，两者不同。

更一般地，若 |d₃−n|≤C，其中 C 是不随 b 增长的固定常数，则由 (1)

    |x−n/q|<(C+1)/b³.

无穷多次这样的事件同样要求无理性指数至少为 3。允许“第三位只差 1、2、3”并没有回到 P2848 的正确强度。

一个明确的严格反例族：y=[0;1,2,3,4,…]。设第 k 个渐近分数为 P_k/Q_k。由连分数尾项公式，

    1/(k+3)<Q_k²|y−P_k/Q_k|<1/(k+1).

故 y 满足 A。但 Q_k>k+3 对 k≥3 成立，因而这些渐近分数的误差大于 1/Q_k³。若任意充分大约简分数误差小于 1/s³，则由 Legendre 定理它必为渐近分数，产生矛盾。于是 y 只有有限多个三次幂精度的逼近，更不可能在无穷多个进制下前三位相同。

此处仅调用经典 Legendre 定理：约简分数 r/s 若满足 |y−r/s|<1/(2s²)，则是渐近分数。参考教学原始来源：

https://sites.millersville.edu/bikenaga/number-theory/approximation-by-rationals/approximation-by-rationals.html

本记录没有把 OpenAI 关于 μ(π)=2 的预印本当作已独立核验的输入。若该结论成立，则它还意味着 π 的前三位全同只能发生有限多次，与 P2848 仍可能成立并不冲突。预印本原文：

https://github.com/openai/math/blob/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/main.tex

## 9. 精确数值实验

用 Machin 恒等式 π=16 arctan(1/5)−4 arctan(1/239)，及反正切交错级数的下一项余项界，得到宽度恰为 10^(−200) 的有理包围区间。用整数运算检查上下端点是否有相同数位前缀；不同则立即报错，绝不把不确定结果当作一个确定数位。

穷举 4≤b≤1,000,000，共 999,997 个进制。前两位相同的进制恰为：

    8, 114, 227, 340, 453, 566, 679, 33103, 66318, 265382.

其中没有前三位完全相同的进制。注意 b=227,…,679 的一部分例子来自 355/113 的未约分倍数，不是新的约简有理数；脚本记录真正分母。

代表性结果（中括号内的每一项是一个数位的十进制数值）：

| 进制 b | 分母 q=b−1 | π 的前三位小数 | q²×有理逼近距离 |
|---:|---:|:---|---:|
| 8 | 7 | [1,1,0] | 0.0619599741001313… |
| 114 | 113 | [16,16,15] | 0.00340631193013807… |
| 1,725,034 | 1,725,033 | [244252,244252,130577] | 0.0658970801129740… |
| 6,701,487,260 | 6,701,487,259 | [948881364,948881364,1027625562] | 0.0117502571361641… |

后两行由认证的连分数递推选出，并非穷举到这些巨大进制。所有数字判定都用上述包围区间核验。

例如 b=114 时，重复数位 16 表示

    3+16/113=355/113.

只用前三位 [16,16,15]、不使用后续任何数位，公式 (2) 就严格给出：

    −0.009840 < 113²(π−355/113) < −0.001220.

实际数值约为 −0.00340631193，位于这个保守区间中。严格端点分数见附录 JSON。

脚本另在 12,000 个固定种子随机有理数/进制样本上，分别检查 K=1,2,3,4 的恒等式和区间，共 48,000 次精确检查，全部通过。这些检查辅助发现实现错误，数学有效性来自上面的证明。

## 10. 接下来真正需要证明的事

经过这轮推导，原始目标可以保留为下面这个精确问题：

    对每个 η>0，是否有任意大的进制 b，
    使 π 的前两位小数相同，
    且第三位与第一位的差小于 ηb？

这与 P2848 等价，尚未解决。初等恒等式没有提供 π 特有的新下界或上界，也没有证明合适进制的无穷性。

值得继续探索的是不同进制之间的约束及上述数位差的下极限；不能仅假设各进制的数位相互独立，也不能把它们随机化后得到的平均结论直接用于指定的 π。

从经典连分数角度，足够小的约简归一化误差最终来自渐近分数，因此当前数位判据提供了另一种表示与有限位证书，没有绕开连分数部分商无界的核心困难。

## 附录 A：复现方法

把附录 B 代码保存为 verify.py，用 Python 3.10 或更新版本运行：

    python3 verify.py

只用 Python 标准库，不需要付费模型、联网、Lean 或第三方包。运行后会在脚本旁生成 results.json。

附录 C 为本轮实际输出的完整 JSON。JSON 中 normalized_error_lower/upper 是方便阅读的有限位近似字符串，不把这些舍入字符串当作严谨端点；精确端点在 normalized_error_exact_interval 中以整数分子、共同整数分母保存。其他表中小数也是显示近似值。


## 附录 B：完整验证程序

```python
"""Variable-base rational-approximation experiment; Python standard library only.

All floor decisions use exact rational enclosures, never float rounding.
The output is a finite certificate, not a proof of infinitely many hits.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from math import isqrt, gcd
from pathlib import Path
import json
import random


def atan_bounds(d, terms):
    total = sum((F((-1)**k, (2*k+1)*d**(2*k+1))
                 for k in range(terms)), F(0))
    following = F((-1)**terms, (2*terms+1)*d**(2*terms+1))
    return min(total, total+following), max(total, total+following)


def pi_bounds():
    a, b = atan_bounds(5, 160)
    c, d = atan_bounds(239, 50)
    return 16*a-4*d, 16*b-4*c


def rounded_enclosure(lo, hi, places=200):
    scale = 10**places
    lower = lo.numerator*scale // lo.denominator
    upper = -(-hi.numerator*scale // hi.denominator)
    return lower, upper, scale


def certified_digits(lower, upper, scale, base, length=3):
    """Both interval endpoints must belong to the same length-digit cell."""
    factor = base**length
    prefix = lower*factor // scale
    if prefix != upper*factor // scale:
        raise ArithmeticError('Insufficient enclosure precision')
    result = []
    for _ in range(length):
        result.append(prefix % base)
        prefix //= base
    assert prefix == 0
    return list(reversed(result))


def decimal_string(value, places=16):
    with localcontext() as context:
        context.prec = places
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def test_exact_identities():
    rng = random.Random(2848)
    for _ in range(12000):
        base = rng.randrange(2, 1000)
        denominator = rng.randrange(2, 10000)
        alpha = F(rng.randrange(1, denominator), denominator)
        q = base-1
        first = (base*alpha).__floor__()
        actual = q*q*(alpha-F(first,q))
        for length in (1,2,3,4):
            power = base**length
            prefix = (power*alpha).__floor__()
            tail = power*alpha-prefix
            repeated_prefix = first*(power-1)//q
            difference = prefix-repeated_prefix
            factor = F(q*q,power)
            reconstructed = factor*(difference+tail-F(first,q))
            assert actual == reconstructed
            lower = factor*(difference-F(first,q))
            upper = factor*(difference+1-F(first,q))
            assert lower <= actual < upper
            assert factor*max(0,abs(difference)-1) <= abs(actual)
            assert abs(actual) < factor*(abs(difference)+1)
    return 48000


def analyze_pi(limit=1000000):
    lo, hi = pi_bounds()
    lower, upper, scale = rounded_enclosure(lo,hi)
    alpha_lower, alpha_upper = lower-3*scale, upper-3*scale
    repeated2, repeated3, record_rows, repeated_rows = 0, [], [], []
    best = None
    for base in range(4,limit+1):
        q = base-1
        d1,d2,d3 = certified_digits(alpha_lower,alpha_upper,scale,base)
        if d1 == d2:
            repeated2 += 1
            factor=F(q*q,base**3)
            cert_lo=factor*(d3-d1-F(d1,q))
            cert_hi=factor*(d3-d1+1-F(d1,q))
            actual_lo=F(q*(q*alpha_lower-d1*scale),scale)
            actual_hi=F(q*(q*alpha_upper-d1*scale),scale)
            assert cert_lo <= actual_lo <= actual_hi < cert_hi
            repeated_rows.append({'base':base,'q':q,'digits':[d1,d2,d3],
                'signed_three_digit_bound_exact':[str(cert_lo),str(cert_hi)]})
            if d2 == d3:
                repeated3.append(base)
        # q times the distance of q*pi to its nearest integer.
        nearest = (2*q*lower+scale)//(2*scale)
        assert nearest == (2*q*upper+scale)//(2*scale)
        gap_lower, gap_upper = q*lower-nearest*scale,q*upper-nearest*scale
        assert gap_lower*gap_upper > 0
        candidate_lower = q*min(abs(gap_lower),abs(gap_upper))
        candidate_upper = q*max(abs(gap_lower),abs(gap_upper))
        if best is None or candidate_upper < best[0]:
            best = candidate_lower,candidate_upper
            record_rows.append({
                'base':base,'q':q,'numerator':nearest,
                'reduced_denominator':q//gcd(nearest,q),
                'digits':[d1,d2,d3],
                'normalized_error_lower':decimal_string(F(candidate_lower,scale),25),
                'normalized_error_upper':decimal_string(F(candidate_upper,scale),25),
                'normalized_error_exact_interval':{
                    'lower_numerator':str(candidate_lower),
                    'upper_numerator':str(candidate_upper),'denominator':str(scale)},
                'relative_third_gap':decimal_string(F(abs(d3-d1),base)),
            })
        elif best[0] <= candidate_upper and candidate_lower <= best[1]:
            raise ArithmeticError('Unresolved comparison of record errors')
    # Certify selected convergents beyond the exhaustive scan.
    cf_lo,cf_hi=F(lower,scale),F(upper,scale)
    p_prev2,p_prev,q_prev2,q_prev=0,1,1,0
    convergents=[]
    for index in range(40):
        a=cf_lo.__floor__()
        assert a == cf_hi.__floor__()
        p=a*p_prev+p_prev2
        q=a*q_prev+q_prev2
        if q>2:
            base=q+1
            digits=certified_digits(alpha_lower,alpha_upper,scale,base)
            error1=F(q*(q*lower-p*scale),scale)
            error2=F(q*(q*upper-p*scale),scale)
            assert error1*error2 > 0
            convergents.append({'index':index,'numerator':p,'q':q,'base':base,
                'digits':digits,'normalized_error':decimal_string(abs((error1+error2)/2),20),
                'relative_third_gap':decimal_string(F(abs(digits[2]-digits[0]),base),20)})
        cf_lo,cf_hi=1/(cf_hi-a),1/(cf_lo-a)
        p_prev2,p_prev,q_prev2,q_prev=p_prev,p,q_prev,q
    return {'base_min':4,'base_max':limit,'tested_bases':limit-3,
            'pi_enclosure_places':200,'pi_enclosure_width_units':upper-lower,
            'repeated_two_count':repeated2,'repeated_three_bases':repeated3,
            'repeated_two_rows':repeated_rows,
            'records':record_rows,'convergents':convergents}


def sqrt2_control():
    scale=10**200
    root=isqrt(2*scale*scale)
    lower,upper=root-scale,root+1-scale
    # sqrt(2)-1=[0;2,2,...].
    p_prev2,p_prev,q_prev2,q_prev=0,1,1,0
    rows=[]
    for index in range(21):
        a=0 if index==0 else 2
        p,q=a*p_prev+p_prev2,a*q_prev+q_prev2
        if q>=5:
            digits=certified_digits(lower,upper,scale,q+1)
            assert digits[0]==digits[1]
            error=abs(F(q*(q*(lower+upper)-2*p*scale),2*scale))
            assert error>F(1,4)
            rows.append({'q':q,'base':q+1,'digits':digits,
                'normalized_error':decimal_string(error,20),
                'relative_third_gap':decimal_string(F(abs(digits[2]-digits[0]),q+1),20)})
        p_prev2,p_prev,q_prev2,q_prev=p_prev,p,q_prev,q
    return rows


if __name__ == '__main__':
    result={'identity_checks':test_exact_identities(),
            'pi':analyze_pi(),'sqrt2_control':sqrt2_control()}
    destination=Path(__file__).with_name('results.json')
    destination.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('Exact identity checks:',result['identity_checks'])
    print('Pi bases tested:',result['pi']['tested_bases'])
    print('First two equal:',result['pi']['repeated_two_count'])
    print('First three equal:',result['pi']['repeated_three_bases'])
    print('First two equal bases:',[row['base'] for row in result['pi']['repeated_two_rows']])
    print('Record rows:')
    for row in result['pi']['records']:
        print({key:value for key,value in row.items() if key!='normalized_error_exact_interval'})
    print('Selected further convergents:')
    for row in result['pi']['convergents']:
        if row['digits'][0]==row['digits'][1] and float(row['normalized_error'])<0.1:
            print(row)
    print('Sqrt2 controls:',result['sqrt2_control'][:4])

```

## 附录 C：本轮实际结果

```json
{
  "identity_checks": 48000,
  "pi": {
    "base_min": 4,
    "base_max": 1000000,
    "tested_bases": 999997,
    "pi_enclosure_places": 200,
    "pi_enclosure_width_units": 1,
    "repeated_two_count": 10,
    "repeated_three_bases": [],
    "repeated_two_rows": [
      {
        "base": 8,
        "q": 7,
        "digits": [
          1,
          1,
          0
        ],
        "signed_three_digit_bound_exact": [
          "-7/64",
          "-7/512"
        ]
      },
      {
        "base": 114,
        "q": 113,
        "digits": [
          16,
          16,
          15
        ],
        "signed_three_digit_bound_exact": [
          "-4859/493848",
          "-226/185193"
        ]
      },
      {
        "base": 227,
        "q": 226,
        "digits": [
          32,
          32,
          29
        ],
        "signed_three_digit_bound_exact": [
          "-160460/11697083",
          "-109384/11697083"
        ]
      },
      {
        "base": 340,
        "q": 339,
        "digits": [
          48,
          48,
          37
        ],
        "signed_three_digit_bound_exact": [
          "-1280403/39304000",
          "-582741/19652000"
        ]
      },
      {
        "base": 453,
        "q": 452,
        "digits": [
          64,
          64,
          39
        ],
        "signed_three_digit_bound_exact": [
          "-1712176/30986559",
          "-4932224/92959677"
        ]
      },
      {
        "base": 566,
        "q": 565,
        "digits": [
          80,
          80,
          31
        ],
        "signed_three_digit_bound_exact": [
          "-15687225/181321496",
          "-1921000/22665187"
        ]
      },
      {
        "base": 679,
        "q": 678,
        "digits": [
          96,
          96,
          12
        ],
        "signed_three_digit_bound_exact": [
          "-38678544/313046839",
          "-38218860/313046839"
        ]
      },
      {
        "base": 33103,
        "q": 33102,
        "digits": [
          4687,
          4687,
          25649
        ],
        "signed_three_digit_bound_exact": [
          "22968797123574/36274552383727",
          "22969892865978/36274552383727"
        ]
      },
      {
        "base": 66318,
        "q": 66317,
        "digits": [
          9390,
          9390,
          45078
        ],
        "signed_three_digit_bound_exact": [
          "26158870034467/48611946451572",
          "156957618151291/291671678709432"
        ]
      },
      {
        "base": 265382,
        "q": 265381,
        "digits": [
          37576,
          37576,
          200470
        ],
        "signed_three_digit_bound_exact": [
          "5736069004659739/9345109457661484",
          "11472208436394639/18690218915322968"
        ]
      }
    ],
    "records": [
      {
        "base": 4,
        "q": 3,
        "numerator": 9,
        "reduced_denominator": 1,
        "digits": [
          0,
          2,
          1
        ],
        "normalized_error_lower": "1.274333882308139146163790",
        "normalized_error_upper": "1.274333882308139146163790",
        "normalized_error_exact_interval": {
          "lower_numerator": "127433388230813914616379044951552595777452459437595238877450133077034765657588098765231342807905361183933277861954075982384460148595524008552823467315633005705255692431744668995003680160654059437343764",
          "upper_numerator": "127433388230813914616379044951552595777452459437595238877450133077034765657588098765231342807905361183933277861954075982384460148595524008552823467315633005705255692431744668995003680160654059437343773",
          "denominator": "100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        },
        "relative_third_gap": "0.25"
      },
      {
        "base": 7,
        "q": 6,
        "numerator": 19,
        "reduced_denominator": 6,
        "digits": [
          0,
          6,
          6
        ],
        "normalized_error_lower": "0.9026644707674434153448382",
        "normalized_error_upper": "0.9026644707674434153448382",
        "normalized_error_exact_interval": {
          "lower_numerator": "90266447076744341534483820193789616890190162249619044490199467691860937369647604939074628768378555264266888552183696070462159405617903965788706130737467977178977230273021324019985279357383762250624908",
          "upper_numerator": "90266447076744341534483820193789616890190162249619044490199467691860937369647604939074628768378555264266888552183696070462159405617903965788706130737467977178977230273021324019985279357383762250624944",
          "denominator": "100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        },
        "relative_third_gap": "0.8571428571428571"
      },
      {
        "base": 8,
        "q": 7,
        "numerator": 22,
        "reduced_denominator": 7,
        "digits": [
          1,
          1,
          0
        ],
        "normalized_error_lower": "0.06195997410013131533047422",
        "normalized_error_upper": "0.06195997410013131533047422",
        "normalized_error_exact_interval": {
          "lower_numerator": "6195997410013131533047421930435867433869943061981477222771497691699609197575906722629355823626366887474376084916697429240161413202147064545738900170442524493607896760501246582757741347550120841128347",
          "upper_numerator": "6195997410013131533047421930435867433869943061981477222771497691699609197575906722629355823626366887474376084916697429240161413202147064545738900170442524493607896760501246582757741347550120841128396",
          "denominator": "100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        },
        "relative_third_gap": "0.125"
      },
      {
        "base": 114,
        "q": 113,
        "numerator": 355,
        "reduced_denominator": 113,
        "digits": [
          16,
          16,
          15
        ],
        "normalized_error_lower": "0.003406311930138070506638904",
        "normalized_error_upper": "0.003406311930138070506638904",
        "normalized_error_exact_interval": {
          "lower_numerator": "340631193013807050663890402767168634393937927377197093250082149230813139729651862331520650715893595108331189822642325869818064861548309888571760742461127732229259894702400310889780956479449395262507",
          "upper_numerator": "340631193013807050663890402767168634393937927377197093250082149230813139729651862331520650715893595108331189822642325869818064861548309888571760742461127732229259894702400310889780956479449395275276",
          "denominator": "100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        },
        "relative_third_gap": "0.008771929824561404"
      }
    ],
    "convergents": [
      {
        "index": 1,
        "numerator": 22,
        "q": 7,
        "base": 8,
        "digits": [
          1,
          1,
          0
        ],
        "normalized_error": "0.061959974100131315330",
        "relative_third_gap": "0.125"
      },
      {
        "index": 2,
        "numerator": 333,
        "q": 106,
        "base": 107,
        "digits": [
          15,
          16,
          10
        ],
        "normalized_error": "0.93505573491682736626",
        "relative_third_gap": "0.046728971962616822430"
      },
      {
        "index": 3,
        "numerator": 355,
        "q": 113,
        "base": 114,
        "digits": [
          16,
          16,
          15
        ],
        "normalized_error": "0.0034063119301380705066",
        "relative_third_gap": "0.0087719298245614035088"
      },
      {
        "index": 4,
        "numerator": 103993,
        "q": 33102,
        "base": 33103,
        "digits": [
          4687,
          4687,
          25649
        ],
        "normalized_error": "0.63321927297600212499",
        "relative_third_gap": "0.63323565839954082712"
      },
      {
        "index": 5,
        "numerator": 104348,
        "q": 33215,
        "base": 33216,
        "digits": [
          4703,
          4702,
          25765
        ],
        "normalized_error": "0.36586380906794851031",
        "relative_third_gap": "0.63409200385356454721"
      },
      {
        "index": 6,
        "numerator": 208341,
        "q": 66317,
        "base": 66318,
        "digits": [
          9390,
          9390,
          45078
        ],
        "normalized_error": "0.53811723974624529987",
        "relative_third_gap": "0.53813444313760969872"
      },
      {
        "index": 7,
        "numerator": 312689,
        "q": 99532,
        "base": 99533,
        "digits": [
          14093,
          14092,
          84889
        ],
        "normalized_error": "0.28871241161940854588",
        "relative_third_gap": "0.71128168547115027177"
      },
      {
        "index": 8,
        "numerator": 833719,
        "q": 265381,
        "base": 265382,
        "digits": [
          37576,
          37576,
          200470
        ],
        "normalized_error": "0.61380486765818164496",
        "relative_third_gap": "0.61380952739824102614"
      },
      {
        "index": 9,
        "numerator": 1146408,
        "q": 364913,
        "base": 364914,
        "digits": [
          51669,
          51668,
          338312
        ],
        "normalized_error": "0.21448855324407582011",
        "relative_third_gap": "0.78550836635481236675"
      },
      {
        "index": 10,
        "numerator": 4272943,
        "q": 1360120,
        "base": 1360121,
        "digits": [
          192583,
          192583,
          1209267
        ],
        "normalized_error": "0.74749406656204237614",
        "relative_third_gap": "0.74749525961293149654"
      },
      {
        "index": 11,
        "numerator": 5419351,
        "q": 1725033,
        "base": 1725034,
        "digits": [
          244252,
          244252,
          130577
        ],
        "normalized_error": "0.065897080112973993924",
        "relative_third_gap": "0.065897251880252795017"
      },
      {
        "index": 12,
        "numerator": 80143857,
        "q": 25510582,
        "base": 25510583,
        "digits": [
          3612111,
          3612111,
          13226130
        ],
        "normalized_error": "0.37686392012315839953",
        "relative_third_gap": "0.37686394701367663765"
      },
      {
        "index": 13,
        "numerator": 165707065,
        "q": 52746197,
        "base": 52746198,
        "digits": [
          7468474,
          7468473,
          36135672
        ],
        "normalized_error": "0.45650680656059565093",
        "relative_third_gap": "0.54349316324183214115"
      },
      {
        "index": 14,
        "numerator": 245850922,
        "q": 78256779,
        "base": 78256780,
        "digits": [
          11080585,
          11080585,
          48548374
        ],
        "normalized_error": "0.47878009058506871915",
        "relative_third_gap": "0.47878010058681177529"
      },
      {
        "index": 15,
        "numerator": 411557987,
        "q": 131002976,
        "base": 131002977,
        "digits": [
          18549059,
          18549058,
          106017473
        ],
        "normalized_error": "0.33231735207421226309",
        "relative_third_gap": "0.66768264357839745886"
      },
      {
        "index": 16,
        "numerator": 1068966896,
        "q": 340262731,
        "base": 340262732,
        "digits": [
          48178703,
          48178703,
          169125014
        ],
        "normalized_error": "0.35544977243112388075",
        "relative_third_gap": "0.35544977344154163789"
      },
      {
        "index": 17,
        "numerator": 2549491779,
        "q": 811528438,
        "base": 811528439,
        "digits": [
          114906465,
          114906464,
          631754328
        ],
        "normalized_error": "0.36311798880051063087",
        "relative_third_gap": "0.63688200950404401047"
      },
      {
        "index": 18,
        "numerator": 6167950454,
        "q": 1963319607,
        "base": 1963319608,
        "digits": [
          277991633,
          277991633,
          855161009
        ],
        "normalized_error": "0.29397627022727851246",
        "relative_third_gap": "0.29397627041883035072"
      },
      {
        "index": 19,
        "numerator": 14885392687,
        "q": 4738167652,
        "base": 4738167653,
        "digits": [
          670889731,
          670889730,
          2086851496
        ],
        "normalized_error": "0.70115836529728447280",
        "relative_third_gap": "0.29884163429790735604"
      },
      {
        "index": 20,
        "numerator": 21053343141,
        "q": 6701487259,
        "base": 6701487260,
        "digits": [
          948881364,
          948881364,
          1027625562
        ],
        "normalized_error": "0.011750257136164143625",
        "relative_third_gap": "0.011750257061594413894"
      },
      {
        "index": 21,
        "numerator": 1783366216531,
        "q": 567663097408,
        "base": 567663097409,
        "digits": [
          80376924307,
          80376924306,
          423454450312
        ],
        "normalized_error": "0.39563179714858975461",
        "relative_third_gap": "0.60436820284939784457"
      },
      {
        "index": 22,
        "numerator": 3587785776203,
        "q": 1142027682075,
        "base": 1142027682076,
        "digits": [
          161702729978,
          161702729978,
          630551584865
        ],
        "normalized_error": "0.41054070951615969982",
        "relative_third_gap": "0.41054070951653070882"
      },
      {
        "index": 23,
        "numerator": 5371151992734,
        "q": 1709690779483,
        "base": 1709690779484,
        "digits": [
          242079654285,
          242079654284,
          965348172153
        ],
        "normalized_error": "0.57695945573961178590",
        "relative_third_gap": "0.42304054425928934403"
      },
      {
        "index": 24,
        "numerator": 8958937768937,
        "q": 2851718461558,
        "base": 2851718461559,
        "digits": [
          403782384263,
          403782384263,
          582854062261
        ],
        "normalized_error": "0.062794304701685582439",
        "relative_third_gap": "0.062794304701489949667"
      },
      {
        "index": 25,
        "numerator": 139755218526789,
        "q": 44485467702853,
        "base": 44485467702854,
        "digits": [
          6298815418230,
          6298815418229,
          36601034561395
        ],
        "normalized_error": "0.31882880617164428636",
        "relative_third_gap": "0.68117119382833727997"
      },
      {
        "index": 26,
        "numerator": 428224593349304,
        "q": 136308121570117,
        "base": 136308121570118,
        "digits": [
          19300228638953,
          19300228638953,
          28937879468775
        ],
        "normalized_error": "0.070704890646331714857",
        "relative_third_gap": "0.070704890646331110048"
      },
      {
        "index": 27,
        "numerator": 5706674932067741,
        "q": 1816491048114374,
        "base": 1816491048114375,
        "digits": [
          257201787724619,
          257201787724618,
          675455169534828
        ],
        "normalized_error": "0.76974652187557840304",
        "relative_third_gap": "0.23025347812442054955"
      },
      {
        "index": 28,
        "numerator": 6134899525417045,
        "q": 1952799169684491,
        "base": 1952799169684492,
        "digits": [
          276502016363572,
          276502016363572,
          638621222738417
        ],
        "normalized_error": "0.18543596904199374983",
        "relative_third_gap": "0.18543596904199397612"
      },
      {
        "index": 29,
        "numerator": 30246273033735921,
        "q": 9627687726852338,
        "base": 9627687726852339,
        "digits": [
          1363209853178907,
          1363209853178906,
          6919995506422308
        ],
        "normalized_error": "0.42283279112334387802",
        "relative_third_gap": "0.57716720887665596461"
      },
      {
        "index": 30,
        "numerator": 66627445592888887,
        "q": 21208174623389167,
        "base": 21208174623389168,
        "digits": [
          3002921722721386,
          3002921722721386,
          6206404246788051
        ],
        "normalized_error": "0.15104942226068549891",
        "relative_third_gap": "0.15104942226068549466"
      },
      {
        "index": 31,
        "numerator": 430010946591069243,
        "q": 136876735467187340,
        "base": 136876735467187341,
        "digits": [
          19380740189507223,
          19380740189507222,
          134055148864601337
        ],
        "normalized_error": "0.16220672356271720634",
        "relative_third_gap": "0.83779327643728279224"
      },
      {
        "index": 32,
        "numerator": 2646693125139304345,
        "q": 842468587426513207,
        "base": 842468587426513208,
        "digits": [
          119287362859764724,
          119287362859764724,
          127722700020716135
        ],
        "normalized_error": "0.010012642948170703053",
        "relative_third_gap": "0.010012642948170703171"
      },
      {
        "index": 33,
        "numerator": 262452630335382199398,
        "q": 83541266890691994833,
        "base": 83541266890691994834,
        "digits": [
          11828829663306214899,
          11828829663306214898,
          36373768201596807548
        ],
        "normalized_error": "0.70619384345217127808",
        "relative_third_gap": "0.29380615654782872190"
      },
      {
        "index": 34,
        "numerator": 265099323460521503743,
        "q": 84383735478118508040,
        "base": 84383735478118508041,
        "digits": [
          11948117026165979623,
          11948117026165979623,
          36383590555406659173
        ],
        "normalized_error": "0.28957563197207625885",
        "relative_third_gap": "0.28957563197207625885"
      },
      {
        "index": 35,
        "numerator": 792651277256425206884,
        "q": 252308737846929010913,
        "base": 252308737846929010914,
        "digits": [
          35725063715638174145,
          35725063715638174144,
          186819388287578353356
        ],
        "normalized_error": "0.40115302442039768724",
        "relative_third_gap": "0.59884697557960231276"
      },
      {
        "index": 36,
        "numerator": 1850401877973371917511,
        "q": 589001211171976529866,
        "base": 589001211171976529867,
        "digits": [
          83398244457442327913,
          83398244457442327913,
          170751177115603052886
        ],
        "normalized_error": "0.14830688121056413615",
        "relative_third_gap": "0.14830688121056413615"
      },
      {
        "index": 37,
        "numerator": 11895062545096656711950,
        "q": 3786316004878788190109,
        "base": 3786316004878788190110,
        "digits": [
          536114530460292141623,
          536114530460292141622,
          3187490467959945991265
        ],
        "normalized_error": "0.29974784616939740068",
        "relative_third_gap": "0.70025215383060259932"
      },
      {
        "index": 38,
        "numerator": 37535589513263342053361,
        "q": 11947949225808341100193,
        "base": 11947949225808341100194,
        "digits": [
          1691741835838318752782,
          1691741835838318752782,
          3732484021424870258076
        ],
        "normalized_error": "0.17080271660164212402",
        "relative_third_gap": "0.17080271660164212402"
      },
      {
        "index": 39,
        "numerator": 199573010111413366978755,
        "q": 63526062133920493691074,
        "base": 63526062133920493691075,
        "digits": [
          8994823709651885905533,
          8994823709651885905532,
          41494467120955066703194
        ],
        "normalized_error": "0.48840456468417533059",
        "relative_third_gap": "0.51159543531582466941"
      }
    ]
  },
  "sqrt2_control": [
    {
      "q": 5,
      "base": 6,
      "digits": [
        2,
        2,
        5
      ],
      "normalized_error": "0.35533905932737622004",
      "relative_third_gap": "0.5"
    },
    {
      "q": 12,
      "base": 13,
      "digits": [
        5,
        5,
        0
      ],
      "normalized_error": "0.35324701827431297256",
      "relative_third_gap": "0.38461538461538461538"
    },
    {
      "q": 29,
      "base": 30,
      "digits": [
        12,
        12,
        23
      ],
      "normalized_error": "0.35360595577293604222",
      "relative_third_gap": "0.36666666666666666667"
    },
    {
      "q": 70,
      "base": 71,
      "digits": [
        29,
        29,
        3
      ],
      "normalized_error": "0.35354437183426087173",
      "relative_third_gap": "0.36619718309859154930"
    },
    {
      "q": 169,
      "base": 170,
      "digits": [
        70,
        70,
        131
      ],
      "normalized_error": "0.35355493796768882503",
      "relative_third_gap": "0.35882352941176470588"
    },
    {
      "q": 408,
      "base": 409,
      "digits": [
        169,
        169,
        24
      ],
      "normalized_error": "0.35355312510579627569",
      "relative_third_gap": "0.35452322738386308068"
    },
    {
      "q": 985,
      "base": 986,
      "digits": [
        408,
        408,
        757
      ],
      "normalized_error": "0.35355343614372361844",
      "relative_third_gap": "0.35395537525354969574"
    },
    {
      "q": 2378,
      "base": 2379,
      "digits": [
        985,
        985,
        143
      ],
      "normalized_error": "0.35355338277805211126",
      "relative_third_gap": "0.35393022278268179908"
    },
    {
      "q": 5741,
      "base": 5742,
      "digits": [
        2378,
        2378,
        4409
      ],
      "normalized_error": "0.35355339193415381160",
      "relative_third_gap": "0.35370950888192267503"
    },
    {
      "q": 13860,
      "base": 13861,
      "digits": [
        5741,
        5741,
        840
      ],
      "normalized_error": "0.35355339036321511675",
      "relative_third_gap": "0.35358199264122357694"
    },
    {
      "q": 33461,
      "base": 33462,
      "digits": [
        13860,
        13860,
        25691
      ],
      "normalized_error": "0.35355339063274558548",
      "relative_third_gap": "0.35356523818062279601"
    },
    {
      "q": 80782,
      "base": 80783,
      "digits": [
        33461,
        33461,
        4899
      ],
      "normalized_error": "0.35355339058650146799",
      "relative_third_gap": "0.35356448757783196960"
    },
    {
      "q": 195025,
      "base": 195026,
      "digits": [
        80782,
        80782,
        149735
      ],
      "normalized_error": "0.35355339059443570419",
      "relative_third_gap": "0.35355798714017618164"
    },
    {
      "q": 470832,
      "base": 470833,
      "digits": [
        195025,
        195025,
        28560
      ],
      "normalized_error": "0.35355339059307440447",
      "relative_third_gap": "0.35355423260476644585"
    },
    {
      "q": 1136689,
      "base": 1136690,
      "digits": [
        470832,
        470832,
        872713
      ],
      "normalized_error": "0.35355339059330796658",
      "relative_third_gap": "0.35355373936605406927"
    },
    {
      "q": 2744210,
      "base": 2744211,
      "digits": [
        1136689,
        1136689,
        166463
      ],
      "normalized_error": "0.35355339059326789366",
      "relative_third_gap": "0.35355371726153710484"
    },
    {
      "q": 6625109,
      "base": 6625110,
      "digits": [
        2744210,
        2744210,
        5086541
      ],
      "normalized_error": "0.35355339059327476908",
      "relative_third_gap": "0.35355352590372084388"
    },
    {
      "q": 15994428,
      "base": 15994429,
      "digits": [
        6625109,
        6625109,
        970224
      ],
      "normalized_error": "0.35355339059327358945",
      "relative_third_gap": "0.35355341537981756023"
    },
    {
      "q": 38613965,
      "base": 38613966,
      "digits": [
        15994428,
        15994428,
        29646527
      ],
      "normalized_error": "0.35355339059327379184",
      "relative_third_gap": "0.35355340086019654132"
    }
  ]
}
```

## 后续补充：π 的统一界与刘维尔数的任意长重复前缀

本节回应“π 可能在任何进制下都不会前 1000 位相同；刘维尔数又怎样”的问题。这里的重复均从小数点后第一位开始，进制取整数 b≥2；不讨论从很后面某一位开始的重复块。

### A. π：经典 Mahler 下界排除前 42 位相同

Mahler《On the approximation of π》(1953)，原文定理 1（原始第 33 页；2019 重印第 561 页）证明：对正整数分子 p 和 q≥2，

    |π−p/q|>q^(−42).

此处是对所有 q≥2 的无未知常数显式不等式，不是只对充分大分母的渐近指数陈述。核对来源：

https://ems.press/books/dms/252/4986

https://content.ems.press/assets/public/full-texts/books/252/chapters/online-pdf/252-chapter-4986.pdf

若在某个 b≥3 进制下，π 的前 42 位小数全是 n，则令 q=b−1，由共同前缀的距离界，

    |π−(3+n/q)|<b^(−42)<q^(−42),

与 Mahler 矛盾。此处 p=3q+n 为正整数，允许未约分；Mahler 的原始不等式直接适用。

二进制单独处理：1/8<π−3<1/4，故前三位小数是 001，不能全同。

所以在任意整数进制下，π 从小数点后第一位开始的同数位连续重复长度至多 41。这个上界并不声称最优。特别地，前 1000 位全同不可能。它不排除从后面某一位置开始的长重复块。

### B. 经典刘维尔常数：每个指定长度都能构造

定义

    L=Σ_{j≥1}10^(−j!)=0.110001000000000000000001… .

十进制的第 1!,2!,3!,… 位是 1，其他位是 0。这个经典常数是刘维尔数的一例；“刘维尔数”是一个数类，不只这一个常数。

定义与背景核对来源：

https://math.hmc.edu/funfacts/liouville-numbers/

https://isa-afp.org/entries/Liouville_Numbers.html

给定正整数 k，令

    q_k=10^(k!),
    b_k=q_k+1,
    R_k=Σ_{j=1}^k 10^(−j!),
    N_k=q_k R_k=Σ_{j=1}^k10^(k!−j!).

N_k 是整数，且 0<N_k<q_k，故是 b_k 进制下的一个合法数位。R_k=N_k/q_k 在 b_k 进制下恰为 0.\overline{N_k}。

尾项满足

    q_k^(−(k+1))<L−R_k<2q_k^(−(k+1)).       (A)

证明下界只需保留第一项 10^(−(k+1)!)；上界可把后续稀疏项放大为全部十进制尾项，用等比级数求和。

又因 R_k<0.12、q_k≥10、k≤q_k，有

    b_k^k (L−R_k)
      <(2/q_k)(1+1/q_k)^k
      <6/q_k≤0.6<0.88<1−R_k.

代入本记录公式 (4)，证明 L 在 b_k 进制下前 k 位小数全部为 N_k。

另一方面，由 (A)，

    L−R_k>q_k^(−(k+1))>(1−R_k)/b_k^(k+1),

所以第 k+1 位必定不再等于 N_k：初始重复长度恰为 k。

因此，取 k=1000，进制 b=10^(1000!)+1，即严格得到前 1000 位全为同一个数位 N_1000。这个进制极大，但结论是证明，不依赖实际打印其展开。

小规模例子已用有理区间实际核验：

| k | 进制 b_k | 重复数位 N_k | 已核验的开头 |
|---:|---:|---:|:---|
| 1 | 11 | 1 | [1][2]… |
| 2 | 101 | 11 | [11][11][12]… |
| 3 | 1,000,001 | 110001 | [110001][110001][110001][110002]… |
| 4 | 10^24+1 | 一个大整数 | 初始重复长度恰为 4 |
| 5 | 10^120+1 | 一个大整数 | 初始重复长度恰为 5 |

表中的 [11]、[110001] 都分别是一个数位。计算时取前六项的精确有理和 L_6，以及上界 L_6+2×10^(−7!)；检查整个有理区间的前缀取整值一致，且下一位不能保持同一数位，五个例子全部通过。

### C. 一般等价：任意长重复前缀恰好刻画刘维尔数

定理。对无理数 x∈(0,1)，下述两条等价：

1. x 是刘维尔数：对每个正整数 M，存在无穷多个有理数 a/q，使 0<|x−a/q|<q^(−M)。
2. 对每个正整数 K，存在一个整数进制，使 x 的前 K 位小数全部相同。

证明 1⇒2：固定 K，选择指数 K+1 的无穷多个逼近 a/q，令 b=q+1、R=a/q。q→∞ 时 R→x，且

    b^K|x−R|<(1+1/q)^K/q→0.

由于 x∈(0,1)，R 最终在 (0,1) 内且远离两端，故公式 (4) 的正负误差条件最终都满足。这样得到的不仅是一个进制，而且是无穷多个进制。

证明 2⇒1：选择长度 K 趋于无穷的重复前缀。对应进制必须趋于无穷：若始终来自某个有限集合，则其中一个固定进制、固定数位会拥有任意长重复前缀，使 x 等于该数位的无限循环小数，成为有理数。令 q=b−1，对应有理数 a/q 的误差小于 b^(−K)<q^(−K)。对每个固定 M，选取 K≥M 的无穷多个例子即得到所需逼近；约分不会使界失效，约简分母也必然无界。

### D. 研究意义与边界

因此，用户提出的对比确实有严格的数学区别：π 的初始重复长度在所有整数进制中统一有界；刘维尔数的初始重复长度无界。

但“统一有界”仅刻画非刘维尔性，并不等同于坏逼近性。P2848 只要求 q² 归一化误差前的常数能任意小；刘维尔数则要求逼近指数本身无界。π 不可能前 1000 位全同，与 P2848 可能成立完全相容。

这些等价改写没有被宣称为新的文献成果；Mahler 定理是明确调用的外部输入，重复前缀与刘维尔性的对应是本记录给出的初等推导。
