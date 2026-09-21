"""压降译井台 - 核心计算模块。

数据口径（详见 README）:
- 时间为自开井起的实际经过时间 (h)，流量阶段按右连续 (right-continuous) 规则取值；
- 多流量等效时间采用标准变流量叠加时间:
      teq(t) = exp( 1/qn * SUM_{k<=n} (qk-q_{k-1}) * ln(t - t_{k-1}) )  (t 在第 n 阶段内)
  在阶段起点样本上采用明确的右连续规则，teq 记 0 并标记 teq_nonpositive；
- 压力变化 dp = p_ref - p_corrected（压降为正）；
- 对数导数在实际、非等间隔的等效时间坐标上用 Bourdet 三点法计算，
  重复时刻与 teq<=0 的样本保留但不参与导数；
- 仪器换档按分段边界处理，导数绝不跨越换档点，平滑不能跨段。
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

# 单位制（内部一致的抽象油藏单位，见 README）
VISCOSITY = 1.0        # mPa*s
FORMATION_FACTOR = 1.0  # rb/stb
PAY_THICKNESS = 10.0    # m

# 边界距离导出的数量级常数（d 与 r_inv 同量纲，用于候选提示，非严格反演）
BOUNDARY_RADIUS_FACTOR = 1.0

DEFAULT_SMOOTH_NEIGHBORS = 2
DEFAULT_SMOOTH_L = 1.2
EPS = 1e-9


def _to_float_list(values: Sequence[Any]) -> List[float]:
    return [float(v) for v in values]


def normalize_stages(stages: Sequence[Dict[str, Any]]) -> List[Dict[str, float]]:
    """归一化流量阶段，按开始时刻排序并校验。"""
    norm = [
        {"start": float(s["start"]), "rate": float(s["rate"])}
        for s in stages
    ]
    norm.sort(key=lambda s: s["start"])
    if not norm:
        raise ValueError("至少需要一个流量阶段")
    if any(norm[i]["start"] == norm[i - 1]["start"] for i in range(1, len(norm))):
        raise ValueError("流量阶段开始时刻不能重复")
    return norm


def stage_index_at(stages: Sequence[Dict[str, float]], t: float) -> int:
    """右连续：返回 t 所属阶段索引（start <= t 的最后一个阶段）。"""
    idx = 0
    for i, s in enumerate(stages):
        if t + EPS >= s["start"]:
            idx = i
        else:
            break
    return idx


def rate_at(stages: Sequence[Dict[str, float]], t: float) -> float:
    return stages[stage_index_at(stages, t)]["rate"]


def equivalent_time(stages: Sequence[Dict[str, float]], t: float) -> float:
    """单个时刻的多流量等效时间；阶段起点（右连续）记 0。"""
    n = stage_index_at(stages, t)
    # 明确的右连续规则：样本恰好落在阶段起点时，等效时间归零（叠加奇异点）。
    if n > 0 and abs(t - stages[n]["start"]) <= EPS:
        return 0.0
    if n == 0 and t <= stages[0]["start"] + EPS:
        return 0.0
    prev_rate = 0.0
    weighted = 0.0
    current = stages[n]["rate"]
    for k in range(n + 1):
        qk = stages[k]["rate"]
        tk = stages[k]["start"]
        dt = t - tk
        if dt <= EPS:
            continue
        weighted += (qk - prev_rate) * math.log(dt)
        prev_rate = qk
    if abs(current) < EPS:
        return 0.0
    teq = math.exp(weighted / current)
    if teq <= 0.0 or not math.isfinite(teq):
        return 0.0
    return teq


def _estimate_shift_offset(times: np.ndarray, pressure: np.ndarray, at: int,
                           window: int = 4) -> float:
    """用换档点前后各 window 个点做局部线性外推，估计换档跳变幅度。"""
    lo = max(0, at - window)
    hi = min(len(times), at + 1 + window)
    left_t = times[lo:at]
    left_p = pressure[lo:at]
    right_t = times[at:hi]
    right_p = pressure[at:hi]
    if len(left_t) < 2 or len(right_t) < 2:
        return float(pressure[at] - pressure[at - 1]) if at > 0 else 0.0
    left_fit = np.polyfit(left_t, left_p, 1)
    right_fit = np.polyfit(right_t, right_p, 1)
    tc = times[at]
    return float(np.polyval(right_fit, tc) - np.polyval(left_fit, tc))


def apply_shifts(times: np.ndarray, pressure: np.ndarray,
                 shifts: Sequence[Dict[str, Any]]) -> Tuple[np.ndarray, np.ndarray, List[Dict[str, float]]]:
    """解析换档并做分段校正。

    换档定义为 index 处仪器记录发生跳变；index 及之后样本属于新段。
    offset 为观测值相对地层真值的跳变；校正后 pressure_corr = pressure - cum_offset。
    返回校正压力、段编号数组和归一化换档列表。
    """
    n = len(times)
    resolved: List[Dict[str, float]] = []
    for sh in shifts:
        idx = int(sh["index"])
        if idx <= 0 or idx >= n:
            continue
        offset = sh.get("offset")
        if offset is None:
            offset = _estimate_shift_offset(times, pressure, idx)
        resolved.append({"index": idx, "offset": float(offset),
                         "time": float(times[idx])})
    resolved.sort(key=lambda s: s["index"])

    corrected = pressure.astype(float).copy()
    segment = np.zeros(n, dtype=int)
    cum = 0.0
    seg_id = 0
    next_idx = resolved[0]["index"] if resolved else n + 1
    ri = 0
    for i in range(n):
        if ri < len(resolved) and i == resolved[ri]["index"]:
            cum += resolved[ri]["offset"]
            seg_id += 1
            ri += 1
            next_idx = resolved[ri]["index"] if ri < len(resolved) else n + 1
        corrected[i] = pressure[i] - cum
        segment[i] = seg_id
    return corrected, segment, resolved


def _segment_groups(segment: np.ndarray) -> List[Tuple[int, int]]:
    groups = []
    start = 0
    for i in range(1, len(segment)):
        if segment[i] != segment[start]:
            groups.append((start, i))
            start = i
    groups.append((start, len(segment)))
    return groups


def _bourdet(x: np.ndarray, y: np.ndarray, seg: np.ndarray,
             neighbors: int = 2, l_window: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
    """非等间隔 Bourdet 三点对数导数，严格限制在同一换档段内。

    d = dy/d(ln x)。选点规则（二者取更靠近的约束）:
    - neighbors: 左/右各取这么多个有效邻居（自适应非等间隔，默认 2）；
    - l_window: 若给出 ln 半窗，则不超过该窗。
    导数绝不跨越换档段；段边缘退化为单侧并标记 one_sided。
    """
    n = len(x)
    d = np.full(n, np.nan)
    one_sided = np.zeros(n, dtype=bool)
    groups = _segment_groups(seg)

    for g0, g1 in groups:
        gx, gy = x[g0:g1], y[g0:g1]
        m = len(gx)
        if m < 2:
            continue
        lx = np.log(gx)
        for j in range(m):
            left = j
            count = 0
            while left - 1 >= 0 and count < neighbors:
                if l_window is not None and (lx[j] - lx[left - 1]) > l_window + EPS:
                    break
                left -= 1
                count += 1
            right = j
            count = 0
            while right + 1 < m and count < neighbors:
                if l_window is not None and (lx[right + 1] - lx[j]) > l_window + EPS:
                    break
                right += 1
                count += 1

            dxl = lx[j] - lx[left] if j > left else 0.0
            dxr = lx[right] - lx[j] if right > j else 0.0
            ml = (gy[j] - gy[left]) / dxl if dxl > EPS else np.nan
            mr = (gy[right] - gy[j]) / dxr if dxr > EPS else np.nan

            if dxl > EPS and dxr > EPS:
                val = (dxr * ml + dxl * mr) / (dxl + dxr)
            elif dxl > EPS:
                val, one_sided[g0 + j] = ml, True
            elif dxr > EPS:
                val, one_sided[g0 + j] = mr, True
            else:
                val = np.nan
            d[g0 + j] = val
    return d, one_sided


@dataclass
class Sample:
    index: int
    time: float
    pressure_raw: float
    pressure: float
    dp: float
    teq: float
    segment: int
    rate: float
    flags: List[str] = field(default_factory=list)
    derivative: float = math.nan
    derivative_real: float = math.nan
    one_sided: bool = False


def build_samples(stages: Sequence[Dict[str, Any]], times: Sequence[float],
                  pressures: Sequence[float],
                  shifts: Optional[Sequence[Dict[str, Any]]] = None,
                  smooth_l: float = DEFAULT_SMOOTH_L,
                  smooth_neighbors: int = DEFAULT_SMOOTH_NEIGHBORS,
                  reference_pressure: Optional[float] = None) -> Dict[str, Any]:
    """构建逐样本诊断：等效时间、压力变化、换档分段与对数导数。"""
    stages = normalize_stages(stages)
    t = np.asarray(_to_float_list(times), dtype=float)
    p = np.asarray(_to_float_list(pressures), dtype=float)

    order = np.argsort(t, kind="stable")
    t, p = t[order], p[order]
    original_index = order

    corrected, segment, resolved = apply_shifts(t, p, shifts or [])
    if reference_pressure is None:
        reference_pressure = float(corrected[0])
    reference_pressure = float(reference_pressure)

    samples: List[Sample] = []
    seen_time: Dict[float, int] = {}
    for i in range(len(t)):
        flags: List[str] = []
        teq = equivalent_time(stages, t[i])
        rate = rate_at(stages, t[i])
        n_stage = stage_index_at(stages, t[i])
        if t[i] in seen_time:
            flags.append("duplicate_time")
        seen_time[t[i]] = i
        if teq <= 0.0:
            flags.append("teq_nonpositive")
        if n_stage > 0 and t[i] <= stages[n_stage]["start"] + EPS:
            flags.append("rate_boundary_sample")
        samples.append(Sample(
            index=int(original_index[i]), time=float(t[i]),
            pressure_raw=float(p[i]), pressure=float(corrected[i]),
            dp=reference_pressure - float(corrected[i]),
            teq=float(teq), segment=int(segment[i]), rate=float(rate),
            flags=flags,
        ))

    # 有效导数点：teq>0 且无重复时刻；重复/非正等效时间保留诊断但不参与。
    valid = np.array(
        [s.teq > 0.0 and "duplicate_time" not in s.flags for s in samples],
        dtype=bool,
    )

    def _run_derivative(xcoord: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        d = np.full(len(t), np.nan)
        one = np.zeros(len(t), dtype=bool)
        groups = _segment_groups(segment)
        for g0, g1 in groups:
            local_valid = np.where(valid[g0:g1])[0]
            if len(local_valid) < 2:
                continue
            xv = xcoord[g0:g1][local_valid]
            yv = np.array([samples[g0 + k].dp for k in local_valid])
            sv = segment[g0:g1][local_valid]
            dd, oo = _bourdet(xv, yv, sv, neighbors=smooth_neighbors, l_window=None)
            for k, gk in enumerate(local_valid):
                d[g0 + gk] = dd[k]
                one[g0 + gk] = oo[k]
        return d, one

    teq_arr = np.array([s.teq for s in samples], dtype=float)
    d_teq, one_teq = _run_derivative(teq_arr)
    d_real, one_real = _run_derivative(t)

    out = []
    for i, s in enumerate(samples):
        s.derivative = float(d_teq[i]) if not math.isnan(d_teq[i]) else math.nan
        s.derivative_real = float(d_real[i]) if not math.isnan(d_real[i]) else math.nan
        s.one_sided = bool(one_teq[i])
        out.append({
            "index": s.index,
            "time": s.time,
            "rate": s.rate,
            "pressure_raw": s.pressure_raw,
            "pressure": s.pressure,
            "dp": s.dp,
            "teq": s.teq,
            "segment": s.segment,
            "derivative": None if math.isnan(s.derivative) else s.derivative,
            "derivative_real": None if math.isnan(s.derivative_real) else s.derivative_real,
            "one_sided": s.one_sided,
            "flags": s.flags,
        })

    return {
        "reference_pressure": reference_pressure,
        "resolved_shifts": resolved,
        "samples": out,
    }


def _linfit(x: np.ndarray, y: np.ndarray,
            fixed_slope: Optional[float] = None) -> Tuple[float, float]:
    if fixed_slope is None:
        slope, intercept = np.polyfit(x, y, 1)
        return float(slope), float(intercept)
    intercept = float(np.mean(y - fixed_slope * x))
    return float(fixed_slope), intercept


def _valid_points(samples: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    pts = []
    for s in samples:
        if (s["teq"] > 0 and s["derivative"] is not None
                and s["dp"] > 0 and s["derivative"] > 0):
            pts.append(s)
    return sorted(pts, key=lambda s: s["teq"])


def detect_candidates(samples: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """识别早期井筒储集、径向流、边界效应三个候选区段（不跨换档段）。"""
    pts = _valid_points(samples)
    if len(pts) < 5:
        return []
    lx = np.array([math.log10(s["teq"]) for s in pts])
    lder = np.array([math.log10(s["derivative"]) for s in pts])

    # 有效点虽按 teq 排序，但不同流量阶段在 teq 轴上会交错。
    # 为检测建立“同换档段+同流量”块内按真实时间排列的邻接表。
    block_members: Dict[Tuple[int, float], List[int]] = {}
    for i, sp in enumerate(pts):
        block_members.setdefault((sp["segment"], round(sp["rate"], 9)), []).append(i)
    for members in block_members.values():
        members.sort(key=lambda k: pts[k]["time"])
    pos_in_block = {}
    for key, members in block_members.items():
        for pos, idx in enumerate(members):
            pos_in_block[idx] = (key, pos)

    def same_block(a: int, b: int) -> bool:
        return (pts[a]["segment"] == pts[b]["segment"]
                and abs(pts[a]["rate"] - pts[b]["rate"]) <= EPS)

    def local_slope(i: int, half: int = 3) -> float:
        # 在同换档段、同流量块内按真实时间邻接取窗（绝不跨阶跃或换档点）。
        key, pos = pos_in_block[i]
        members = block_members[key]
        a = max(0, pos - half)
        b = min(len(members), pos + half + 1)
        idxs = members[a:b]
        if len(idxs) < 3:
            return np.nan
        return float(np.polyfit(lx[idxs], lder[idxs], 1)[0])

    slopes = np.array([local_slope(i) for i in range(len(pts))])
    ratio = np.array([s["derivative"] / s["dp"] for s in pts])

    # WBS：最早一段 der/dp ≈ 1（单元斜率）
    wbs_pred = [0.7 <= ratio[i] <= 1.85 for i in range(len(pts))]
    wbs_run = (-1, -1)
    i = 0
    while i < len(wbs_pred):
        if wbs_pred[i]:
            j = i
            while j + 1 < len(wbs_pred) and wbs_pred[j + 1] and same_block(j, j + 1):
                j += 1
            if j - i + 1 >= 3:
                wbs_run = (i, j + 1)
                break
            i = j + 1
        else:
            i += 1

    # 径向流：在“同换档段+同流量”块内枚举平坦区间（|局部斜率|<=0.22）。
    # 平台基准取每个块导数的“高位聚簇”：对块导数做直方图式稳健估计
    # （取上四分位与最大值之间的中位数），再用 ±15% 门限排除阶跃瞬态与早期过渡。
    plateau_ref = {}
    for key in {(sp["segment"], round(sp["rate"], 9)) for sp in pts}:
        vals = np.array(sorted(
            sp["derivative"] for sp in pts
            if sp["segment"] == key[0] and abs(sp["rate"] - key[1]) <= EPS))
        q3 = float(np.quantile(vals, 0.6))
        cluster = vals[vals >= q3]
        plateau_ref[key] = float(np.median(cluster))

    def near_plateau(i: int) -> bool:
        key = (pts[i]["segment"], round(pts[i]["rate"], 9))
        ref = plateau_ref[key]
        return ref > 0 and abs(pts[i]["derivative"] / ref - 1.0) <= 0.15

    def radial_ok(i: int) -> bool:
        return (not np.isnan(slopes[i])) and abs(slopes[i]) <= 0.22 and near_plateau(i)

    rad_members: List[int] = []
    for members in block_members.values():
        k = 0
        while k < len(members):
            idx = members[k]
            if radial_ok(idx):
                k2 = k
                while k2 + 1 < len(members) and radial_ok(members[k2 + 1]):
                    k2 += 1
                if k2 - k + 1 >= 4:
                    cand = members[k:k2 + 1]
                    vals = np.array([pts[m]["derivative"] for m in cand])
                    med = float(np.median(vals))
                    stable = float(np.mean(np.abs(vals / med - 1.0))) <= 0.06
                    if stable and len(cand) > len(rad_members):
                        rad_members = cand
                k = k2 + 1
            else:
                k += 1

    # 边界：在最末换档段内按真实时间顺序找导数单调上升的尾部。
    bnd_members: List[int] = []
    last_seg = pts[-1]["segment"]
    tail = [i for i, sp in enumerate(pts) if sp["segment"] == last_seg]
    tail.sort(key=lambda k: pts[k]["time"])
    if len(tail) >= 4:
        end = len(tail)
        start = end - 1
        while start - 1 >= 0 and (lder[tail[start]] - lder[tail[start - 1]]) >= -0.05:
            start -= 1
        if (end - start) >= 4 and (lder[tail[end - 1]] - lder[tail[start]]) >= math.log10(1.2):
            bnd_members = tail[start:end]

    groups = [("wbs", ("r", wbs_run)), ("radial", ("m", rad_members)),
              ("boundary", ("m", bnd_members))]
    out = []
    for kind, run in groups:
        if run[0] == "m":
            idxs = sorted(run[1], key=lambda k: pts[k]["teq"])
        else:
            a, b = run[1]
            if a < 0:
                continue
            idxs = list(range(a, b))
        if not idxs:
            continue
        grp = [pts[k] for k in idxs]
        out.append({
            "regime": kind,
            "left": grp[0]["teq"],
            "right": grp[-1]["teq"],
            "left_open": False,
            "right_open": False,
            "sample_count": len(grp),
            "segments": sorted({sp["segment"] for sp in grp}),
            "rates": sorted({sp["rate"] for sp in grp}),
            "sample_indices": [sp["index"] for sp in grp],
            "locked_slope": None,
        })
    return out


def interval_points(samples: Sequence[Dict[str, Any]], interval: Dict[str, Any]) -> List[Dict[str, Any]]:
    """按开闭端点与换档段筛选有效样本。"""
    left, right = float(interval["left"]), float(interval["right"])
    segs = interval.get("segments")
    rates = interval.get("rates")
    allowed = set(interval.get("sample_indices")) if interval.get("sample_indices") else None
    res = []
    for sp in _valid_points(samples):
        ok_left = sp["teq"] > left if interval.get("left_open") else sp["teq"] >= left
        ok_right = sp["teq"] < right if interval.get("right_open") else sp["teq"] <= right
        ok_seg = (segs is None) or (sp["segment"] in segs)
        ok_rate = (rates is None) or any(abs(sp["rate"] - r) <= EPS for r in rates)
        ok_idx = (allowed is None) or (sp["index"] in allowed)
        if ok_left and ok_right and ok_seg and ok_rate and ok_idx:
            res.append(sp)
    return res


def derive_interval(samples: Sequence[Dict[str, Any]], interval: Dict[str, Any]) -> Dict[str, Any]:
    """计算区段拟合与导出参数。"""
    pts = interval_points(samples, interval)
    regime = interval["regime"]
    result: Dict[str, Any] = {
        "regime": regime,
        "effective_samples": len(pts),
        "smooth_l": float(interval.get("smooth_l") or DEFAULT_SMOOTH_L),
    }
    if len(pts) < 2:
        result["status"] = "insufficient_samples"
        return result

    x = np.array([math.log(s["teq"]) for s in pts])
    yp = np.array([s["dp"] for s in pts])
    x10 = np.array([math.log10(s["teq"]) for s in pts])
    yd10 = np.array([math.log10(s["derivative"]) for s in pts])

    if regime == "wbs":
        slope_dp, intercept_dp = _linfit(x, yp)
        result.update({
            "status": "ok",
            "slope_dp_dlnt": slope_dp,
            "intercept_dp": intercept_dp,
            "storage_coefficient": slope_dp / float(np.mean([s["rate"] for s in pts])),
            "note": "C = (d dp / d teq)/q （teq 为等效时间）",
        })
    elif regime == "radial":
        locked = interval.get("locked_slope")
        slope, intercept = _linfit(x10, yd10,
                                   fixed_slope=float(locked) if locked is not None else None)
        if locked is not None:
            plateau = 10.0 ** intercept  # 锁定斜率时由截距给出 teq=1 处的平台
        else:
            plateau = float(np.mean([10.0 ** v for v in yd10]))
        q_mean = float(np.mean([s["rate"] for s in pts]))
        permeability = None
        if abs(q_mean) > EPS:
            permeability = 162.6 * q_mean * VISCOSITY * FORMATION_FACTOR / (PAY_THICKNESS * plateau)
        result.update({
            "status": "ok",
            "slope_loglog": slope,
            "plateau_derivative": float(plateau),
            "permeability": permeability,
            "rate_used": q_mean,
            "locked": locked is not None,
            "note": "m = 162.6 q mu B / (h * 导数平台)",
        })
    else:  # boundary
        slope, intercept = _linfit(x10, yd10)
        teq_intercept = 10.0 ** (-intercept / slope) if abs(slope) > EPS else None
        result.update({
            "status": "ok",
            "slope_loglog": slope,
            "intercept_loglog": intercept,
            "teq_unit_line_intercept": teq_intercept,
            "boundary_radius_estimate": (BOUNDARY_RADIUS_FACTOR * teq_intercept) if teq_intercept else None,
            "note": "斜率上抬指示边界；r_b 为数量级估计",
        })
    return result


ANALYSIS_VERSION = "pwgsb-analysis-1.0"


def canonical_json(obj: Any) -> str:
    import json

    return json.dumps(obj, sort_keys=True, ensure_ascii=False,
                      separators=(",", ":"), default=str)


def run_fingerprint(payload: Dict[str, Any]) -> str:
    """运行指纹：包含流量修订版本、换档标注、平滑尺度、分析器版本。"""
    import hashlib

    material = {
        "analysis_version": ANALYSIS_VERSION,
        "rate_revision": payload.get("rate_revision"),
        "stages": payload.get("stages"),
        "shifts": payload.get("shifts"),
        "smooth_l": payload.get("smooth_l", DEFAULT_SMOOTH_L),
        "reference_pressure": payload.get("reference_pressure"),
        "intervals": payload.get("intervals"),
        "data_digest": payload.get("data_digest"),
    }
    return hashlib.sha256(canonical_json(material).encode("utf-8")).hexdigest()[:16]


def next_rate_revision(revisions: Sequence[Dict[str, Any]]) -> int:
    if not revisions:
        return 1
    return max(int(r.get("version", 0)) for r in revisions) + 1
