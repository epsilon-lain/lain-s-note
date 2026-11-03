#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""WSL 下的 π 连分数纪录搜索器（Python 3.10+，依赖 gmpy2）。

安装与启动：
  python3 -m venv ~/pi-search-venv
  ~/pi-search-venv/bin/pip install gmpy2
  ~/pi-search-venv/bin/python pi_record_search.py \
      --digits 1000000 --max-digits 16000000 --workdir ~/pi-search

Ctrl+C 安全暂停；再次执行相同命令自动续跑。长计算中的暂停需等当前
GMP/MPFR 调用返回。--max-digits 是精度上限，不是已经搜索过的项数。

默认从 π 的开头出发，逐步复现已知范围；本程序没有 Max Frank 的外部
检查点，不能直接跳到第 6535.2 亿项。百万位规模用于观察本机速度。
批量小矩阵减少完整大整数运算次数，但本实现不是外存 half-GCD 引擎，
不声称适用于数千亿项规模，也不把小规模测速线性外推到该规模。

数学目标：C(p,q)=q²|π-p/q|，p/q 既约，q>0，C 越小越好。
公开基线 A=845188380288；若公开项数据正确，它对应的 C>1/(A+2)。
默认要求候选 C<1/(A+2)，这是充分的保守门槛；只用 A 无法判断极小幅
刷新同一大部分商下的精确纪录。可用 --target-quality 分子/分母 改门槛。

认证方法：MPFR 向下/向上舍入分别计算 π，形成有理区间；每批候选
连分数前缀都对整个区间验证。无法确认时提高精度，达到上限则保存并
停止，不继续猜测。上下界、整数矩阵、候选分数和断点均使用精确整数。
这依赖 GMP/MPFR 的正确性，未提供 Lean 证明。

输出：checkpoint.json（带 SHA-256 的断点）、history.json（摘要）、
records/*.json（完整分子分母的十六进制表示及误差证书）。
验证某条结果：python pi_record_search.py --verify records/文件名.json
先试一个可快速命中的自定义目标：
  python pi_record_search.py --digits 10000 --max-digits 10000 \
      --target-quality 1/10000 --workdir demo-pi

参考：
  https://oeis.org/A033089
  https://oeis.org/A033090
  https://gmpy2.readthedocs.io/en/stable/
  https://www.mpfr.org/mpfr-current/mpfr.html
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import sys
import time

try:
    import gmpy2 as g
except ImportError:
    raise SystemExit("缺少 gmpy2。请先在虚拟环境中执行：python -m pip install gmpy2")

VERSION = "1.0.1"
SCHEMA = "pi-cf-search-v1"
REFERENCE_A = 845188380288
WORLD_TARGET = g.mpq(1, REFERENCE_A + 2)
IDENTITY = (1, 0, 0, 1)
STOP = False


def message(text):
    print(text, flush=True)


def on_signal(signum, frame):
    global STOP
    STOP = True


def zhex(n):
    return g.mpz(n).digits(16)


def rat_data(x):
    return [zhex(x.numerator), zhex(x.denominator)]


def read_rat(data):
    return g.mpq(g.mpz(data[0], 16), g.mpz(data[1], 16))


def display(x):
    with g.context(precision=80):
        return format(g.mpfr(x), ".12g")


def canonical(data):
    return json.dumps(data, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":")).encode("utf-8")


def atomic_json(path, data):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as stream:
        stream.write(canonical(data))
        stream.write(b"\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def pi_bounds(bits):
    # MPFR's const_pi is correctly rounded in the requested direction.
    with g.context(precision=bits, round=g.RoundDown):
        lower = g.mpq(g.const_pi())
    with g.context(precision=bits, round=g.RoundUp):
        upper = g.mpq(g.const_pi())
    if not 3 < lower < upper < 4:
        raise ArithmeticError("π 区间异常")
    return lower, upper


def precision_bits(digits):
    # 3.32193 > log_2(10), plus guard bits.
    return (digits * 332193 + 99999) // 100000 + 96


def inverse_endpoint(matrix, numerator, denominator):
    a, b, c, d = matrix
    n = b * denominator - d * numerator
    den = c * numerator - a * denominator
    if den < 0:
        n, den = -n, -den
    if den <= 0 or n <= den:
        return None
    return n, den


def inverse_interval(matrix, tail):
    if matrix == IDENTITY:
        return tail
    ln, ld, un, ud = tail
    # Both endpoints must be on the same side of the pole a/c.
    a, b, c, d = matrix
    dl = c * ln - a * ld
    du = c * un - a * ud
    if dl == 0 or du == 0 or (dl < 0) != (du < 0):
        return None
    x = inverse_endpoint(matrix, ln, ld)
    y = inverse_endpoint(matrix, un, ud)
    if x is None or y is None:
        return None
    # The inverse Mobius derivative has the sign of det(matrix).
    # This avoids cross-multiplying two enormous rational endpoints per block.
    if a * d - b * c < 0:
        x, y = y, x
    return x[0], x[1], y[0], y[1]


def multiply(left, right):
    a, b, c, d = left
    e, f, h, k = right
    return (a * e + b * h, a * f + b * k,
            c * e + d * h, c * f + d * k)


def proposed_block(tail, cap=62):
    """Predict from high bits; certify the WHOLE prefix against exact bounds."""
    ln, ld, un, ud = tail
    shift = max(0, max(ln.bit_length(), ld.bit_length()) - 192)
    n, den = int(ln >> shift), int(ld >> shift)
    matrix = IDENTITY
    plan = []
    while den:
        a, remainder = divmod(n, den)
        next_matrix = multiply(matrix, (a, 1, 1, 0))
        if max(abs(v).bit_length() for v in next_matrix) > cap:
            break
        plan.append((a, next_matrix))
        matrix = next_matrix
        n, den = den, remainder
    while plan:
        if inverse_interval(plan[-1][1], tail) is not None:
            return plan
        plan.pop()
    # Very large a, or close to a boundary: use an exact one-term check.
    a = ln // ld
    if a != un // ud:
        return []
    matrix = (a, 1, 1, 0)
    if inverse_interval(matrix, tail) is None:
        return []
    return [(a, matrix)]


def quality_from_tail(tail, q, q_previous):
    ln, ld, un, ud = tail
    return (g.mpq(ud * q, un * q + ud * q_previous),
            g.mpq(ld * q, ln * q + ld * q_previous))


def quality_from_pi(p, q, lower, upper):
    left = q * lower - p
    right = q * upper - p
    if left <= 0 <= right:
        low = g.mpq(0)
    else:
        low = q * min(abs(left), abs(right))
    high = q * max(abs(left), abs(right))
    return low, high


def encode_best(best):
    if best is None:
        return None
    return {"p_hex": zhex(best["p"]), "q_hex": zhex(best["q"]),
            "index": best["index"], "following_a": str(best["a"]),
            "lower": rat_data(best["lower"]), "upper": rat_data(best["upper"])}


def decode_best(data):
    if data is None:
        return None
    return {"p": g.mpz(data["p_hex"], 16), "q": g.mpz(data["q_hex"], 16),
            "index": data["index"], "a": g.mpz(data["following_a"]),
            "lower": read_rat(data["lower"]), "upper": read_rat(data["upper"])}


def save_checkpoint(state, folder):
    payload = {
        "schema": SCHEMA, "version": VERSION,
        "digits": state["digits"], "bits": state["bits"],
        "terms": state["terms"], "status": state["status"],
        "target": rat_data(state["target"]),
        "matrix": [zhex(v) for v in state["matrix"]],
        "tail": [zhex(v) for v in state["tail"]],
        "pi": [rat_data(v) for v in state["pi"]],
        "best": encode_best(state["best"]),
        "max_a": str(state["max_a"]), "history": state["history"],
        "elapsed": state["elapsed"],
    }
    digest = hashlib.sha256(canonical(payload)).hexdigest()
    atomic_json(folder / "checkpoint.json", {"sha256": digest, "payload": payload})
    atomic_json(folder / "history.json", state["history"])


def load_checkpoint(path):
    envelope = json.loads(path.read_text(encoding="utf-8"))
    payload = envelope["payload"]
    if hashlib.sha256(canonical(payload)).hexdigest() != envelope["sha256"]:
        raise ValueError("断点校验失败，文件可能不完整；未覆盖此文件")
    if payload["schema"] != SCHEMA:
        raise ValueError("不支持此断点格式")
    state = dict(payload)
    state["matrix"] = tuple(g.mpz(v, 16) for v in payload["matrix"])
    state["tail"] = tuple(g.mpz(v, 16) for v in payload["tail"])
    state["target"] = read_rat(payload["target"])
    state["pi"] = tuple(read_rat(v) for v in payload["pi"])
    state["best"] = decode_best(payload["best"])
    state["max_a"] = g.mpz(payload["max_a"])
    a, b, c, d = state["matrix"]
    if a * d - b * c != (-1 if state["terms"] % 2 else 1):
        raise ValueError("断点矩阵不满足连分数行列式恒等式")
    ln, ld, un, ud = state["tail"]
    if min(ln, ld, un, ud) <= 0 or ln * ud >= un * ld:
        raise ValueError("断点区间异常")
    return state


def refine(state, digits):
    bits = precision_bits(digits)
    message(f"[精度] 计算 π 的有理上下界：{digits:,} 位，{bits:,} bit …")
    lower, upper = pi_bounds(bits)
    tail = (lower.numerator, lower.denominator, upper.numerator, upper.denominator)
    if state["terms"]:
        tail = inverse_interval(state["matrix"], tail)
        if tail is None:
            return False
    state["pi"] = (lower, upper)
    state["tail"] = tail
    state["digits"], state["bits"] = digits, bits
    if state["best"]:
        best = dict(state["best"])
        best["lower"], best["upper"] = quality_from_pi(best["p"], best["q"], lower, upper)
        state["best"] = best
    message("[精度] 上下界就绪；从已确认的连分数位置接续。")
    return True


def record_document(best, state, target):
    lo, hi = best["lower"], best["upper"]
    return {
        "schema": SCHEMA, "version": VERSION,
        "meaning": "C=q^2*abs(pi-p/q); lower <= C <= upper",
        "convergent_index_zero_based": best["index"],
        "following_partial_quotient": str(best["a"]),
        "p_hex": zhex(best["p"]), "q_hex": zhex(best["q"]),
        "q_bits": best["q"].bit_length(),
        "quality_lower_hex": rat_data(lo), "quality_upper_hex": rat_data(hi),
        "quality_display": [display(lo), display(hi)],
        "pi_precision_bits": state["bits"],
        "target_quality_hex": rat_data(target),
        "target_met": bool(hi < target),
        "conservative_reference_threshold_met": bool(hi < WORLD_TARGET),
        "reference_A": str(REFERENCE_A),
        "reference_caveat": "The public record is assumed correct; this threshold is sufficient, not necessary.",
        "verification": "Recompute directed MPFR pi bounds with --verify; no Lean certificate.",
    }


def process_block(state, plan, folder, keep_going):
    """Commit a whole verified prefix, or stop before an ambiguous comparison."""
    before = IDENTITY
    accepted = 0
    best = state["best"]
    cutoff = (int(best["upper"].denominator // best["upper"].numerator) - 1
              if best else -1)
    max_a = state["max_a"]
    events = []
    need_precision = False
    hit = False
    for a, after in plan:
        i = state["terms"] + accepted
        # Usually skip all large products for small partial quotients.
        possible = i > 0 and a >= cutoff
        if possible:
            matrix = multiply(state["matrix"], before)
            local_tail = inverse_interval(before, state["tail"])
            if local_tail is None:
                raise ArithmeticError("内部前缀认证失败")
            lo, hi = quality_from_tail(local_tail, matrix[2], matrix[3])
            if best is None or hi < best["lower"]:
                best = {"p": matrix[0], "q": matrix[2], "index": i - 1,
                        "a": a, "lower": lo, "upper": hi}
                cutoff = int(hi.denominator // hi.numerator) - 1
                events.append(best)
                hit = hit or hi < state["target"]
            elif lo < best["upper"]:
                need_precision = True
                break
        max_a = max(max_a, a)
        accepted += 1
        before = after
        if hit and not keep_going:
            break
    if accepted:
        new_tail = inverse_interval(before, state["tail"])
        if new_tail is None:
            raise ArithmeticError("批次尾项认证失败")
        # Write evidence before changing the checkpoint state.
        summaries = []
        for event in events:
            filename = f"convergent_{event['index']:012d}.json"
            doc = record_document(event, state, state["target"])
            atomic_json(folder / "records" / filename, doc)
            summary = {"index": event["index"], "following_a": str(event["a"]),
                       "q_bits": event["q"].bit_length(),
                       "quality": doc["quality_display"],
                       "file": "records/" + filename, "target_met": doc["target_met"]}
            summaries.append(summary)
            kind = "达到目标" if doc["target_met"] else "本地纪录"
            message(f"[{kind}] 收敛下标 {event['index']:,}；后项 {event['a']}；"
                    f"q={event['q'].bit_length():,} bit；C∈[{display(event['lower'])}, {display(event['upper'])}]")
        state["tail"] = new_tail
        state["matrix"] = multiply(state["matrix"], before)
        state["terms"] += accepted
        state["best"] = best
        state["max_a"] = max_a
        state["history"].extend(summaries)
    return accepted, need_precision, hit


def verify_record(filename):
    doc = json.loads(Path(filename).read_text(encoding="utf-8"))
    p, q = g.mpz(doc["p_hex"], 16), g.mpz(doc["q_hex"], 16)
    if q <= 0 or g.gcd(p, q) != 1:
        raise ValueError("分数不是正分母的既约分数")
    bits = max(int(doc["pi_precision_bits"]) + 64, 2 * q.bit_length() + 256)
    target = read_rat(doc["target_quality_hex"])
    old_lo, old_hi = read_rat(doc["quality_lower_hex"]), read_rat(doc["quality_upper_hex"])
    for _ in range(4):
        lower, upper = pi_bounds(bits)
        lo, hi = quality_from_pi(p, q, lower, upper)
        if hi < old_lo or lo > old_hi:
            raise ValueError("重新计算的质量区间与记录不相容")
        if lo < old_lo or hi > old_hi:
            bits *= 2
            continue
        message(f"[验证通过] 既约；重新计算区间包含在证书区间内。C∈[{display(lo)}, {display(hi)}]")
        message(f"[目标] C < {display(target)}：{'已确认' if hi < target else '未确认'}")
        message(f"[公开基线保守门槛] {'已超过' if hi < WORLD_TARGET else '尚未超过'}")
        return 0
    raise ArithmeticError("精度提高后仍不能复核证书；未宣称通过")


def acquire_lock(folder):
    # fcntl is available in WSL. The kernel releases it even after a crash.
    import fcntl
    stream = (folder / "run.lock").open("a+")
    try:
        fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        stream.close()
        raise ValueError("这个工作目录已有搜索进程，请使用另一目录")
    return stream


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", type=Path, default=Path("pi-search"), help="结果与断点目录")
    parser.add_argument("--digits", type=int, default=1000000, help="起始 π 十进制精度，默认 1000000")
    parser.add_argument("--max-digits", type=int, default=16000000, help="自动提高精度的上限，默认 16000000")
    parser.add_argument("--max-terms", type=int, default=0, help="本次最多新增确认多少项，0 表示不限")
    parser.add_argument("--target-quality", default=None, help="自定义 C 门槛，如 1/10000；默认公开纪录保守门槛")
    parser.add_argument("--keep-going", action="store_true", help="达到目标后继续寻找更小的本地 C")
    parser.add_argument("--checkpoint-seconds", type=float, default=30, help="定时保存间隔，默认 30 秒")
    parser.add_argument("--progress-seconds", type=float, default=2, help="进度输出间隔，默认 2 秒")
    parser.add_argument("--verify", type=Path, help="独立复核一条 records JSON 后退出")
    args = parser.parse_args(argv)
    if args.verify:
        return verify_record(args.verify)
    if min(args.digits, args.max_digits) < 30 or args.max_terms < 0:
        parser.error("精度至少 30 位，--max-terms 不得为负")
    if args.checkpoint_seconds <= 0 or args.progress_seconds <= 0:
        parser.error("保存和进度间隔必须为正数")
    target = g.mpq(args.target_quality) if args.target_quality else None
    if target is not None and not 0 < target < g.mpq(1, 2):
        parser.error("目标质量必须在 0 与 1/2 之间")
    folder = args.workdir.expanduser().resolve()
    folder.mkdir(parents=True, exist_ok=True)
    lock = acquire_lock(folder)
    try:
        (folder / "records").mkdir(exist_ok=True)
        checkpoint = folder / "checkpoint.json"
        if checkpoint.is_file():
            state = load_checkpoint(checkpoint)
            if target is not None and target != state["target"]:
                raise ValueError("目标与已有断点不一致；请用另一个 --workdir")
            message(f"[续跑] 已确认 {state['terms']:,} 项；π 精度 {state['digits']:,} 位。")
        else:
            state = {"terms": 0, "matrix": tuple(g.mpz(v) for v in IDENTITY),
                     "tail": (g.mpz(3), g.mpz(1), g.mpz(4), g.mpz(1)),
                     "pi": (g.mpq(3), g.mpq(4)), "best": None,
                     "max_a": g.mpz(0), "history": [], "digits": 0, "bits": 0,
                     "target": target if target is not None else WORLD_TARGET,
                     "status": "new", "elapsed": 0.0}
        message(f"[目标] C=q²|π-p/q| < {display(state['target'])}")
        message("[说明] 本地纪录是当前程序复现范围内的纪录；默认从 π 开头扫描。")
        if state["best"] and state["best"]["upper"] < state["target"] and not args.keep_going:
            message("[完成] 此工作目录已经达到目标；--keep-going 可继续搜索。")
            return 0
        signal.signal(signal.SIGINT, on_signal)
        signal.signal(signal.SIGTERM, on_signal)
        start = time.monotonic()
        elapsed_before = state["elapsed"]
        base_terms = state["terms"]
        term_limit = base_terms + args.max_terms if args.max_terms else None
        cap = max(args.max_digits, args.digits, state["digits"])
        last_save = last_progress = start
        state["status"] = "running"
        try:
            if args.digits > state["digits"]:
                if not refine(state, args.digits):
                    raise ArithmeticError("要求的精度不足以接续断点")
                state["elapsed"] = elapsed_before + time.monotonic() - start
                save_checkpoint(state, folder)
            while not STOP:
                if term_limit is not None and state["terms"] >= term_limit:
                    state["status"] = "term_budget"
                    break
                plan = proposed_block(state["tail"])
                if term_limit is not None:
                    plan = plan[:term_limit - state["terms"]]
                accepted, needs_precision, hit = (0, True, False)
                if plan:
                    accepted, needs_precision, hit = process_block(state, plan, folder, args.keep_going)
                if hit and not args.keep_going:
                    state["status"] = "target_hit"
                    break
                if not accepted or needs_precision:
                    if state["digits"] >= cap:
                        state["status"] = "precision_limit"
                        break
                    next_digits = min(cap, max(state["digits"] + 100, state["digits"] * 2))
                    if not refine(state, next_digits):
                        raise ArithmeticError("提高精度后仍无法认证断点尾项")
                    state["elapsed"] = elapsed_before + time.monotonic() - start
                    save_checkpoint(state, folder)
                    last_save = time.monotonic()
                now = time.monotonic()
                if now - last_progress >= args.progress_seconds:
                    rate = (state["terms"] - base_terms) / max(now - start, 1e-9)
                    current = display(state["best"]["upper"]) if state["best"] else "暂无"
                    message(f"[进度] {state['terms']:,} 项；本轮约 {rate:,.0f} 项/秒；"
                            f"最大部分商 {state['max_a']}；本地 C≤{current}")
                    last_progress = now
                if now - last_save >= args.checkpoint_seconds:
                    state["elapsed"] = elapsed_before + now - start
                    save_checkpoint(state, folder)
                    last_save = time.monotonic()
            if STOP:
                state["status"] = "paused"
        except Exception:
            state["status"] = "error"
            raise
        finally:
            state["elapsed"] = elapsed_before + time.monotonic() - start
            if state["digits"]:
                save_checkpoint(state, folder)
                message(f"[已保存] {state['terms']:,} 项；状态 {state['status']}；目录 {folder}")
        if state["status"] == "precision_limit":
            message("[精度上限] 有效位数已用完或比较无法认证；提高 --max-digits 后运行同一命令即可续跑。")
        return 0
    finally:
        lock.close()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, ArithmeticError, OSError, KeyError, TypeError) as error:
        message(f"[错误] {error}")
        raise SystemExit(1)
