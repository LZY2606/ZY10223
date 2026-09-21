"""核心计算验收：右连续边界、非等间隔导数、换档分段、重复时刻、候选与指纹。"""
import math

import analysis
import numpy as np


def test_fixture_shape(fixture):
    times = fixture["times"]
    # 两个流量阶跃恰落在样本时刻
    assert 10.0 in times and 600.0 in times
    # 每个阶跃时刻都有短读重复
    assert times.count(10.0) == 2
    assert times.count(600.0) == 2
    # 一次仪器换档
    sh = fixture["known_shifts"][0]
    assert times[sh["index"]] == sh["time"] == 1360.0
    assert len(times) == len(fixture["pressures"])


def test_right_continuous_boundary(fixture):
    """阶段起点样本采用右连续规则：属于新阶段，teq=0 并标记。"""
    res = analysis.build_samples(
        fixture["stages"], fixture["times"], fixture["pressures"],
        shifts=[], reference_pressure=fixture["reference_pressure"])
    at10 = [s for s in res["samples"] if s["time"] == 10.0]
    at600 = [s for s in res["samples"] if s["time"] == 600.0]
    assert len(at10) == 2 and len(at600) == 2
    for s in at10:
        assert s["rate"] == 15.0          # 右连续：取阶跃后流量
        assert s["teq"] == 0.0
        assert "rate_boundary_sample" in s["flags"]
        assert "teq_nonpositive" in s["flags"]
    for s in at600:
        assert s["rate"] == 25.0
        assert s["teq"] == 0.0
        assert "rate_boundary_sample" in s["flags"]
    # 阶跃前最后一个样本仍取旧流量
    pre9 = [s for s in res["samples"] if s["time"] == 9.0][0]
    pre520 = [s for s in res["samples"] if s["time"] == 520.0][0]
    assert pre9["rate"] == 0.0 and pre520["rate"] == 15.0


def test_duplicate_times_retained_not_derivative(fixture):
    res = analysis.build_samples(
        fixture["stages"], fixture["times"], fixture["pressures"],
        shifts=[], reference_pressure=fixture["reference_pressure"])
    dups = [s for s in res["samples"] if "duplicate_time" in s["flags"]]
    assert len(dups) == 2
    for s in dups:
        assert s["derivative"] is None       # 保留诊断但不参与导数
        assert s["derivative_real"] is None


def test_nonpositive_teq_retained(fixture):
    res = analysis.build_samples(
        fixture["stages"], fixture["times"], fixture["pressures"],
        shifts=[], reference_pressure=fixture["reference_pressure"])
    flagged = [s for s in res["samples"] if "teq_nonpositive" in s["flags"]]
    # 开井前 6 个短读 + 两个阶段起点(各含重复共 4 个) = 10
    assert len(flagged) >= 10
    assert all(s["derivative"] is None for s in flagged)


def test_equivalent_time_multirate_formula():
    stages = [{"start": 0, "rate": 0}, {"start": 10, "rate": 15}, {"start": 600, "rate": 25}]
    # 第二阶段内（t=135）只有 15 的阶跃：teq = t-10
    assert abs(analysis.equivalent_time(stages, 135.0) - 125.0) < 1e-9
    # 第三阶段内（t=900）叠加：15*ln(t-10) + 10*ln(t-600)，除以 25
    t = 900.0
    expected = math.exp((15 * math.log(t - 10) + 10 * math.log(t - 600)) / 25)
    assert abs(analysis.equivalent_time(stages, t) - expected) < 1e-6


def test_nonuniform_derivative_on_teq(fixture):
    """非等间隔 Bourdet 导数在实际（非重采样）等效时间坐标上计算。"""
    res = analysis.build_samples(
        fixture["stages"], fixture["times"], fixture["pressures"],
        shifts=[{"index": fixture["known_shifts"][0]["index"]}],
        reference_pressure=fixture["reference_pressure"])
    radial = [s for s in res["samples"] if 56 <= s["teq"] <= 260
              and s["segment"] == 0 and s["derivative"] is not None]
    assert len(radial) >= 5
    # 径向平台 ≈ 0.5*q = 7.5（单位制核），容差含噪声与平滑
    plateau = np.median([s["derivative"] for s in radial])
    assert 6.8 < plateau < 8.2
    # 实际时间坐标导数也存在且有限
    real = [s for s in radial if s["derivative_real"] is not None]
    assert real and all(math.isfinite(s["derivative_real"]) for s in real)


def test_shift_segmentation_not_smoothed(fixture):
    """换档点成为分段边界，导数绝不跨段；跳变不被平滑成地层响应。"""
    idx = fixture["known_shifts"][0]["index"]
    res = analysis.build_samples(
        fixture["stages"], fixture["times"], fixture["pressures"],
        shifts=[{"index": idx}],
        reference_pressure=fixture["reference_pressure"])
    by_idx = {s["index"]: s for s in res["samples"]}
    before = by_idx[idx - 1]
    at = by_idx[idx]
    after = by_idx[idx + 1]
    assert before["segment"] == 0 and at["segment"] == 1 and after["segment"] == 1
    # 估计跳变接近已知 18 kPa（局部线性外推，容差给足）
    resolved = res["resolved_shifts"][0]
    assert abs(resolved["offset"] - fixture["known_shifts"][0]["offset"]) < 4.0
    # 换档两侧压力经校正后连续（跳变被移除）
    assert abs(after["pressure"] - before["pressure"]) < 8.0
    # 换档点的导数只能来自新段单侧（不使用旧段点）
    assert at["one_sided"] is True
    # 原始跳变仍然保留在 pressure_raw（不覆盖原始读数）
    assert abs(at["pressure_raw"] - before["pressure_raw"] -
               fixture["known_shifts"][0]["offset"]) < 4.0


def test_shift_derivative_does_not_cross_segment():
    """构造纯线性+单点跳变，验证导数不被跳变污染。"""
    times = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    pressure = [1000.0 - 5.0 * t for t in times]
    for i in range(6, 10):
        pressure[i] += 100.0   # 仪器换档后读数整体上抬，非地层响应
    res = analysis.build_samples(
        [{"start": 0, "rate": 10}], times, pressure,
        shifts=[{"index": 6, "offset": 100.0}], reference_pressure=1000.0)
    by_t = {s["time"]: s for s in res["samples"]}
    assert by_t[4.0]["segment"] == 0 and by_t[5.0]["segment"] == 1
    # 换档边界点只能用新段单侧导数
    assert by_t[5.0]["one_sided"] is True
    for s in res["samples"]:
        if s["derivative"] is not None and not s["one_sided"]:
            # 校正后 Δp=5t 线性；导数不出现跳变量级（100）的尖峰
            assert s["derivative"] < 60.0
            assert abs(s["derivative"] - s["dp"]) / s["dp"] < 0.30


def test_candidates_three_regimes(fixture):
    res = analysis.build_samples(
        fixture["stages"], fixture["times"], fixture["pressures"],
        shifts=[{"index": fixture["known_shifts"][0]["index"]}],
        reference_pressure=fixture["reference_pressure"])
    cands = {c["regime"]: c for c in analysis.detect_candidates(res["samples"])}
    assert set(["wbs", "radial", "boundary"]).issubset(cands)
    # 径向平台导出渗透率合理（q=15, m=7.5 → k≈32.5）
    radial = analysis.derive_interval(res["samples"], cands["radial"])
    assert 6.9 < radial["plateau_derivative"] < 8.1
    assert 28 < radial["permeability"] < 38
    assert radial["rate_used"] == 15.0
    # WBS 是最早簇
    assert cands["wbs"]["left"] < cands["radial"]["left"] < cands["boundary"]["left"]
    # 边界在换档后的最末段、导数上抬
    assert cands["boundary"]["segments"] == [1]


def test_open_endpoint_filtering(fixture):
    res = analysis.build_samples(
        fixture["stages"], fixture["times"], fixture["pressures"],
        shifts=[{"index": fixture["known_shifts"][0]["index"]}],
        reference_pressure=fixture["reference_pressure"])
    cands = {c["regime"]: c for c in analysis.detect_candidates(res["samples"])}
    rad = dict(cands["radial"])
    closed = len(analysis.interval_points(res["samples"], rad))
    rad["left_open"] = rad["right_open"] = True
    opened = analysis.interval_points(res["samples"], rad)
    # 两端各去掉一个样本
    assert len(opened) == closed - 2


def test_locked_plateau_slope(fixture):
    res = analysis.build_samples(
        fixture["stages"], fixture["times"], fixture["pressures"],
        shifts=[{"index": fixture["known_shifts"][0]["index"]}],
        reference_pressure=fixture["reference_pressure"])
    rad = analysis.detect_candidates(res["samples"])[1]
    rad["locked_slope"] = 0.0
    out = analysis.derive_interval(res["samples"], rad)
    assert out["locked"] is True
    assert abs(out["slope_loglog"]) < 1e-9


def test_fingerprint_changes_with_rate_revision(fixture):
    stages = fixture["stages"]
    shifts = [{"index": fixture["known_shifts"][0]["index"]}]
    base = analysis.run_fingerprint({
        "rate_revision": 0, "stages": stages, "shifts": shifts,
        "smooth_l": 1.2, "reference_pressure": fixture["reference_pressure"]})
    revised_stages = [dict(s) for s in stages]
    revised_stages[1]["rate"] = 16.0
    rev = analysis.run_fingerprint({
        "rate_revision": 1, "stages": revised_stages, "shifts": shifts,
        "smooth_l": 1.2, "reference_pressure": fixture["reference_pressure"]})
    assert base != rev
