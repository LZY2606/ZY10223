"""生成固定验收 fixture：短读重复时刻、一次仪器换档跳变、两个流量阶跃恰落在样本时刻。

确定性（LCG 噪声 + 固定舍入），可随时重新生成比对。
单位：时间 h，流量 m3/h，压力 kPa，内部一致的抽象油藏单位。
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

P0 = 30000.0
STAGES = [{"start": 0.0, "rate": 0.0},
          {"start": 10.0, "rate": 15.0},
          {"start": 600.0, "rate": 25.0}]
# 开井前短读（q=0）；开井 t=10 与第二次阶跃 t=600 各有一条短读重复。
# 开井初期加密以覆盖井筒储集单元斜率段。
TIMES = [0.01, 0.05, 0.2, 1.0, 4.0, 9.0,
         10.0, 10.0, 10.02, 10.05, 10.1, 10.2, 10.4, 10.8, 11.5,
         13.0, 16.0, 22.0, 32.0, 46.0, 66.0, 95.0, 135.0, 190.0,
         270.0, 380.0, 520.0, 600.0, 600.0, 600.2, 600.6,
         620.0, 680.0, 780.0, 900.0, 1030.0, 1180.0, 1360.0,
         1560.0, 1800.0, 2080.0, 2400.0, 2760.0, 3000.0]
SHIFT_INDEX = 37          # t=1360 处仪器换档
SHIFT_OFFSET = 18.0       # 观测值上跳 18 kPa
STORAGE_TAU = 1.5
BOUNDARY_TB = 1400.0
NOISE_SIGMA = 0.02
SEED = 20260922



def _e1_scalar(x: float) -> float:
    """指数积分 E1(x)（无 scipy 依赖）。"""
    if x <= 0.0:
        return float("inf")
    gamma = 0.5772156649015329
    if x <= 4.0:
        # E1(x) = -gamma - ln x - sum_{k>=1} (-x)^k / (k * k!)
        total = 0.0
        term = 1.0
        for k in range(1, 100):
            term *= -x / k
            total += term / k
            if abs(term / k) < 1e-15:
                break
        return -gamma - math.log(x) - total
    tiny = 1e-30
    b = x + 1.0
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 200):
        a2 = -i * i
        b += 2.0
        d = a2 * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + a2 / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-14:
            break
    return math.exp(-x) * h


_E1 = np.frompyfunc(_e1_scalar, 1, 1)


def _line_source(u):
    u = np.asarray(u, dtype=float)
    arg = np.divide(1.0, 4.0 * u, out=np.full_like(u, 1e30), where=u > 0)
    return 0.5 * np.asarray(_E1(np.minimum(arg, 700.0)), dtype=float)


def _wbs(u):
    """井筒储集核：u*exp(-u/tau) 叠加线源；小 u 单元斜率，单调过渡到径向平台。"""
    u = np.asarray(u, dtype=float)
    return _line_source(u) + u * np.exp(-u / STORAGE_TAU)


def _boundary(u):
    """封闭边界晚期项：u>tb 后线性抬升，使对数导数持续上抬（pseudosteady 趋势）。"""
    u = np.asarray(u, dtype=float)
    x = u / BOUNDARY_TB
    return np.where(u > BOUNDARY_TB, (x - 1.0) * 0.5, 0.0)


def ideal_pressure(t: float) -> float:
    prev_q = 0.0
    val = P0
    for st in STAGES:
        if t + 1e-12 < st["start"]:
            break
        dq = st["rate"] - prev_q
        u = t - st["start"]
        if u >= -1e-12 and abs(dq) > 0:
            val -= dq * (_wbs(np.array(u)).item()
                         + _boundary(np.array(u)).item())
        prev_q = st["rate"]
    return val


def lcg_normal(seed: int, n: int) -> np.ndarray:
    out = []
    state = seed
    for _ in range(n):
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        u1 = (state % 100000) / 100000.0 + 1e-6
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        u2 = (state % 100000) / 100000.0 + 1e-6
        out.append(math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2))
    return np.array(out[:n])


def main() -> None:
    noise = lcg_normal(SEED, len(TIMES)) * NOISE_SIGMA
    pressures = []
    for i, t in enumerate(TIMES):
        p = ideal_pressure(t) + noise[i]
        if i >= SHIFT_INDEX:
            p += SHIFT_OFFSET
        pressures.append(round(p, 3))

    fixture = {
        "fixture_version": "pwgsb-fixture-1.0",
        "description": "变流量压降：两次阶跃恰落在样本时刻(含短读重复)，t=1360 仪器换档上跳 18 kPa",
        "units": {"time": "h", "rate": "m3/h", "pressure": "kPa"},
        "reservoir": {"viscosity": 1.0, "formation_volume_factor": 1.0,
                      "pay_thickness_m": 10.0},
        "stages": STAGES,
        "times": TIMES,
        "pressures": pressures,
        "reference_pressure": P0,
        "known_shifts": [{"index": SHIFT_INDEX, "offset": SHIFT_OFFSET,
                          "time": TIMES[SHIFT_INDEX]}],
    }
    out = Path(__file__).resolve().parents[1] / "fixtures" / "well_test_fixture.json"
    out.write_text(json.dumps(fixture, ensure_ascii=False, indent=2) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
