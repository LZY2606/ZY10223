"""Pressure drawdown interpretation engine.

Conventions:

* time     hours since test start (t = 0)
* rate     m^3/day, positive production, right-continuous at step times
* pressure kPa absolute raw gauge readings

Variable-rate equivalent time (geometric superposition)::

    t_eq = exp(sum_k (dq_k / q_N) * ln(t - t_k))      t_k < t

which satisfies d ln t_eq / dt = S(t) / q_N with
S(t) = sum_{t_k<t} dq_k / (t - t_k), so the radial plateau stays constant
across rate changes.  Step instants, zero/negative current rate, duplicate
timestamps and gauge shifts are retained in the diagnostic output but are
never used as derivative neighbours.

The Bourdet log derivative is computed in actual-time coordinates: neighbour
search uses |ln(delta t)| on the real time axis, while the secants are taken
against ln t_eq.  Neighbours never cross a gauge-shift segment, therefore a
shift jump cannot be smeared into a "formation" response.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np

HOUR_TO_S = 3600.0
KPA_TO_PA = 1000.0
M3D_TO_M3S = 1.0 / 86400.0

# Reservoir constants shared with the fixture (documented in README).
VISCOSITY = 1.0e-3            # Pa.s
FVF = 1.0
TOTAL_COMPRESSIBILITY = 4.0e-10   # 1/Pa
POROSITY = 0.2
WELL_RADIUS = 0.1             # m
FORMATION_THICKNESS = 20.0    # m


@dataclass
class RateStep:
    time: float
    rate: float


@dataclass
class Shift:
    time: float
    offset: float = 0.0        # cumulative kPa correction applied from here


@dataclass
class Sample:
    idx: int
    time: float
    pressure_raw: float


@dataclass
class ComputedSample:
    idx: int
    time: float
    pressure_raw: float
    pressure_corrected: float
    delta_p: Optional[float]
    segment: int
    current_rate: float
    teq: Optional[float]
    teq_note: str
    derivative: Optional[float]
    derivative_note: str
    duplicate: bool


@dataclass
class FitResult:
    regime: str
    n_points: int
    indices: list[int]
    slope: Optional[float] = None
    value: Optional[float] = None
    permeability_md: Optional[float] = None
    storage_si: Optional[float] = None
    storage_field: Optional[float] = None
    skin: Optional[float] = None
    boundary_distance_m: Optional[float] = None
    doubling_ratio: Optional[float] = None
    note: str = ""


def canonical_steps(rows: list[dict]) -> list[RateStep]:
    """Collapse (time, rate) revision rows to right-continuous steps."""
    latest: dict[float, float] = {}
    for row in rows:
        latest[float(row["time"])] = float(row["rate"])
    steps = [RateStep(t, q) for t, q in latest.items()]
    steps.sort(key=lambda s: s.time)
    if not steps or steps[0].time != 0.0:
        steps.insert(0, RateStep(0.0, 0.0))
    return steps


def rate_at(steps: list[RateStep], t: float) -> float:
    """Right-continuous lookup: a step at exactly t already applies."""
    current = steps[0].rate
    for step in steps:
        if step.time <= t:
            current = step.rate
        else:
            break
    return current


def equivalent_time(steps: list[RateStep], t: float):
    """Return (t_eq, note).  t_eq is None when it is not diagnostically valid."""
    qn = rate_at(steps, t)
    if qn <= 0.0:
        return None, "nonpositive_rate"
    prior = [s for s in steps if s.time < t]
    aligned = any(s.time == t for s in steps if s.time > 0.0)
    if aligned:
        # Singular instant of a rate step that lands exactly on a sample.
        return 0.0, "rate_step_aligned"
    if not prior:
        return 0.0, "initial"
    log_sum = 0.0
    q_prev = 0.0
    for step in prior:
        dq = step.rate - q_prev
        dt = t - step.time
        if dt <= 0.0:
            return None, "non_positive_elapsed"
        log_sum += (dq / qn) * math.log(dt)
        q_prev = step.rate
    teq = math.exp(log_sum)
    if not (teq > 0.0) or not math.isfinite(teq):
        return None, "non_positive_teq"
    return teq, "ok"


def assign_segments(times: list[float], shift_times: list[float]) -> list[int]:
    """Right-continuous: a shift at exactly a sample time owns that sample."""
    shifts = sorted(shift_times)
    return [sum(1 for ts in shifts if ts <= t) for t in times]


def estimate_shifts(samples: list[Sample], shift_times: list[float]) -> list[Shift]:
    """Estimate the constant kPa jump of each marked gauge shift.

    The local pre-shift trend is extrapolated to the first post-shift
    sample, so a smooth formation response is not mistaken for the jump.
    """
    times = [s.time for s in samples]
    raw = [s.pressure_raw for s in samples]
    segments = assign_segments(times, shift_times)
    counts: dict[float, int] = {}
    active_index: dict[float, int] = {}
    usable = [True] * len(samples)
    for i, t in enumerate(times):
        counts[t] = counts.get(t, 0) + 1
        if t in active_index:
            usable[active_index[t]] = False
        active_index[t] = i
    shifts: list[Shift] = []
    cumulative = 0.0
    for seg_index, st in enumerate(sorted(shift_times), start=1):
        left = [i for i, g in enumerate(segments) if g == seg_index - 1
                and usable[i]]
        right = [i for i, g in enumerate(segments) if g == seg_index
                 and usable[i]]
        if not left or not right:
            shifts.append(Shift(float(st), cumulative))
            continue
        i_r = right[0]
        i_l = left[-1]
        pre = left[-3:]
        expected = raw[i_l]
        if len(pre) >= 2:
            # Local log-time trend: a radial/near-radial pressure is linear
            # in ln t, and a short window keeps the extrapolation unbiased.
            xs = [math.log(max(times[j], 1.0e-12)) for j in pre]
            ys = [raw[j] for j in pre]
            xm = sum(xs) / len(xs)
            ym = sum(ys) / len(ys)
            den = sum((x - xm) ** 2 for x in xs)
            slope = (sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / den
                     if den else 0.0)
            expected = ym + slope * (math.log(max(times[i_r], 1.0e-12)) - xm)
        cumulative += raw[i_r] - expected
        shifts.append(Shift(float(st), cumulative))
    return shifts


def bourdet_derivatives(cs: list[ComputedSample], L: float = 0.0) -> None:
    """Fill .derivative / .derivative_note in place.

    Neighbours are the nearest eligible samples in actual-time log distance
    within the same gauge segment, no farther than exp(L) in time ratio.
    Secants use x = ln t_eq; a non-monotonic x is reported, never bridged.
    """
    n = len(cs)
    eligible = [
        i
        for i, c in enumerate(cs)
        if (not c.duplicate and c.teq is not None and c.teq > 0.0
            and c.teq_note == "ok" and c.delta_p is not None)
    ]
    pool = set(eligible)
    for i in range(n):
        c = cs[i]
        if i not in pool:
            if c.derivative_note == "":
                c.derivative_note = (
                    "duplicate_timestamp" if c.duplicate
                    else ("no_positive_teq" if c.teq_note != "ok" else "ineligible")
                )
            continue
        seg = c.segment
        ldist = [(abs(math.log(cs[j].time / c.time)), j)
                 for j in eligible
                 if j < i and cs[j].segment == seg and cs[j].time > 0.0]
        rdist = [(abs(math.log(cs[j].time / c.time)), j)
                 for j in eligible
                 if j > i and cs[j].segment == seg]
        if L > 0.0:
            ldist = [x for x in ldist if x[0] <= L]
            rdist = [x for x in rdist if x[0] <= L]
        if not ldist or not rdist:
            c.derivative = None
            c.derivative_note = "one_sided_segment_boundary"
            continue
        j = min(ldist)[1]
        k = min(rdist)[1]
        x_l, x_i, x_r = math.log(cs[j].teq), math.log(c.teq), math.log(cs[k].teq)
        if not (x_l < x_i < x_r):
            c.derivative = None
            c.derivative_note = "non_monotonic_teq"
            continue
        y_l, y_i, y_r = cs[j].delta_p, c.delta_p, cs[k].delta_p
        dl, dr = x_i - x_l, x_r - x_i
        m_l = (y_i - y_l) / dl
        m_r = (y_r - y_i) / dr
        c.derivative = (dr / (dl + dr)) * m_l + (dl / (dl + dr)) * m_r
        c.derivative_note = "ok"


def compute(
    samples: list[Sample],
    steps: list[RateStep],
    shifts: Optional[list[Shift]] = None,
    smoothing: float = 0.0,
) -> dict:
    """Run the full diagnostic computation."""
    steps = sorted(steps, key=lambda s: s.time)
    shift_times = [s.time for s in (shifts or [])]
    segments = assign_segments([s.time for s in samples], shift_times)
    if shifts is None or all(abs(s.offset) < 1.0e-12 for s in shifts):
        shifts = estimate_shifts(samples, shift_times)
    seg_offset = {0: 0.0}
    for si, sh in enumerate(shifts, start=1):
        seg_offset[si] = sh.offset

    # Duplicate timestamps: the last reading at a time is the active one.
    seen: dict[float, int] = {}
    duplicate = [False] * len(samples)
    for i, s in enumerate(samples):
        if s.time in seen:
            duplicate[seen[s.time]] = False
            duplicate[i] = False
            duplicate[seen[s.time]] = True
        seen[s.time] = i

    p0 = samples[0].pressure_raw
    computed: list[ComputedSample] = []
    for s, seg, dup in zip(samples, segments, duplicate):
        teq, note = equivalent_time(steps, s.time)
        corrected = s.pressure_raw - seg_offset.get(seg, 0.0)
        computed.append(ComputedSample(
            idx=s.idx,
            time=s.time,
            pressure_raw=s.pressure_raw,
            pressure_corrected=corrected,
            delta_p=None if dup else p0 - corrected,
            segment=seg,
            current_rate=rate_at(steps, s.time),
            teq=teq,
            teq_note=note,
            derivative=None,
            derivative_note="",
            duplicate=dup,
        ))
    bourdet_derivatives(computed, L=smoothing)
    return {
        "samples": computed,
        "steps": steps,
        "shifts": shifts,
        "diagnostics": {
            "duplicates": sum(duplicate),
            "step_aligned": sum(1 for c in computed if c.teq_note == "rate_step_aligned"),
            "nonpositive_teq": sum(
                1 for c in computed if c.teq_note in ("nonpositive_rate", "non_positive_teq")
            ),
            "shift_segments": len(set(segments)),
            "one_sided": sum(
                1 for c in computed if c.derivative_note == "one_sided_segment_boundary"
            ),
        },
    }


def _interval_members(
    computed: list[ComputedSample],
    start: float,
    end: float,
    start_open: bool,
    end_open: bool,
    require_derivative: bool,
) -> list[ComputedSample]:
    out = []
    for c in computed:
        if c.duplicate:
            continue
        lo_ok = (c.time > start) if start_open else (c.time >= start)
        hi_ok = (c.time < end) if end_open else (c.time <= end)
        if not (lo_ok and hi_ok):
            continue
        if require_derivative and (c.derivative is None or c.delta_p is None):
            continue
        if not require_derivative and c.delta_p is None:
            continue
        out.append(c)
    return out


def _fit_slope(x: list[float], y: list[float]) -> float:
    xm, ym = sum(x) / len(x), sum(y) / len(y)
    den = sum((xi - xm) ** 2 for xi in x)
    return sum((xi - xm) * (yi - ym) for xi, yi in zip(x, y)) / den if den else float("nan")


def fit_interval(
    computed: list[ComputedSample],
    regime: str,
    start: float,
    end: float,
    *,
    start_open: bool = False,
    end_open: bool = False,
    smoothing: float = 0.0,
    plateau_override: Optional[float] = None,
    net_thickness: float = FORMATION_THICKNESS,
) -> FitResult:
    """Fit one regime interval.  Regime in storage | radial | boundary."""
    members = _interval_members(
        computed, start, end, start_open, end_open,
        require_derivative=(regime != "storage"),
    )
    indices = [c.idx for c in members]
    if regime == "storage":
        # Wellbore storage: Np = C * dp, with cumulative surface volume
        # Np = integral of rate dt (constant rate per interval member).
        # Use elapsed time from the earliest member (the storage restart),
        # fit C = sum(Np*dp)/sum(dp^2) through the origin, all SI units.
        t0 = min(c.time for c in members)
        np_pa = []
        dp_pa = []
        for c in sorted(members, key=lambda z: z.time):
            elapsed_s = (c.time - t0) * HOUR_TO_S
            np_pa.append(c.current_rate * M3D_TO_M3S * elapsed_s)
            dp_pa.append(c.delta_p * KPA_TO_PA)
        den = sum(d * d for d in dp_pa)
        c_storage = sum(v * d for v, d in zip(np_pa, dp_pa)) / den if den else float("nan")
        m_kpa_per_h = None
        if members and abs(members[0].current_rate) > 1e-12:
            m_kpa_per_h = (members[0].current_rate * M3D_TO_M3S
                           / (c_storage * KPA_TO_PA)) * HOUR_TO_S
        log_slope = (
            _fit_slope([math.log(max(c.time - t0, 1e-12)) for c in members],
                       [math.log(c.delta_p) for c in members])
            if len(members) >= 2 else float("nan")
        )
        return FitResult(
            regime=regime, n_points=len(members), indices=indices,
            slope=log_slope, value=m_kpa_per_h,
            storage_si=c_storage,
            storage_field=c_storage * KPA_TO_PA,
            note="C in m^3/Pa; field column kPa per (m^3 -> kPa)",
        )

    if not members:
        return FitResult(regime=regime, n_points=0, indices=indices,
                         note="no eligible derivative samples in interval")
    x = [math.log(c.teq) for c in members]
    y = [math.log(c.derivative) for c in members]
    log_slope = _fit_slope(x, y) if len(members) >= 2 else float("nan")
    plateau = plateau_override
    if plateau is None:
        plateau = float(math.exp(sum(y) / len(y)))      # geometric mean

    if regime == "radial":
        # Last (highest-rate) rate in the interval governs the plateau.
        qn = max(c.current_rate for c in members) * M3D_TO_M3S
        m_kpa = plateau
        k_si = VISCOSITY * FVF * qn / (4.0 * math.pi * net_thickness * (m_kpa * KPA_TO_PA))
        permeability_md = k_si * 1.01325e15
        mid = members[len(members) // 2]
        te = mid.teq * HOUR_TO_S
        # Standard line-source drawdown skin (t_eq in seconds).
        skin = 0.5 * (
            (mid.delta_p * KPA_TO_PA) / (m_kpa * KPA_TO_PA)
            - math.log(4.0 * k_si * te
                       / (2.2458 * POROSITY * VISCOSITY
                          * TOTAL_COMPRESSIBILITY * WELL_RADIUS ** 2))
        )
        return FitResult(
            regime=regime, n_points=len(members), indices=indices,
            slope=log_slope, value=plateau,
            permeability_md=permeability_md, skin=skin,
            note="permeability from locked/fitted plateau; skin at interval midpoint",
        )

    # boundary ----------------------------------------------------------------
    # Doubling ratio is measured against the radial plateau of the same
    # constant-rate segment: the minimum derivative of this rate's run is the
    # closest available pre-boundary level.
    same_rate = [
        c for c in computed
        if (not c.duplicate and c.derivative is not None
            and c.segment == members[0].segment
            and abs(c.current_rate - members[0].current_rate) < 1.0e-9
            and c.derivative > 0)
    ]
    early_level = min(c.derivative for c in same_rate) if same_rate else None
    if early_level is None:
        pre = [c.derivative for c in members][: max(1, len(members) // 3)]
        early_level = math.exp(sum(math.log(v) for v in pre) / len(pre))
    late_level = max(c.derivative for c in members)
    ratio = late_level / early_level
    qn = max(c.current_rate for c in members) * M3D_TO_M3S
    k_si = VISCOSITY * FVF * qn / (4.0 * math.pi * net_thickness * (early_level * KPA_TO_PA))
    # Sealing-fault approximation from the boundary onset (actual elapsed
    # time since test start for the first flow; a multi-rate correction
    # would use rate-normalized time, documented in README).
    t_first = members[0].time * HOUR_TO_S
    distance = math.sqrt(
        k_si * t_first
        / (0.25 * POROSITY * VISCOSITY * TOTAL_COMPRESSIBILITY)
    )
    return FitResult(
        regime=regime, n_points=len(members), indices=indices,
        slope=log_slope, value=late_level,
        permeability_md=k_si * 1.01325e15,
        boundary_distance_m=distance, doubling_ratio=ratio,
        note="sealing fault approximation (doubling ratio 1 -> 2)",
    )


def suggest_candidates(
    computed: list[ComputedSample],
    smoothing: float = 0.0,
) -> dict:
    """Return up to one candidate interval per regime.

    Radial and boundary runs are searched inside contiguous groups that share
    a gauge segment and a (constant) current rate, so a step restart or a gauge
    jump can never be joined into one candidate.
    """
    good = [
        c for c in computed
        if c.derivative is not None and c.delta_p is not None and c.derivative > 0
    ]
    result: dict[str, tuple[float, float]] = {}
    if not good:
        return result

    groups: list[list[ComputedSample]] = []
    current = [good[0]]
    for prev, nxt in zip(good, good[1:]):
        same = (prev.segment == nxt.segment
                and abs(prev.current_rate - nxt.current_rate) < 1.0e-9)
        if same:
            current.append(nxt)
        else:
            groups.append(current)
            current = [nxt]
    groups.append(current)

    # Storage: within a constant-rate run, the earliest points whose dp grows
    # as t^1.  Uses all positive-dp points (a derivative is not required,
    # storage diagnostics read directly off dp vs actual time).
    usable_dp = [
        c for c in computed
        if not c.duplicate and c.delta_p is not None and c.delta_p > 0
    ]
    rate_groups: list[list[ComputedSample]] = []
    rg = [usable_dp[0]] if usable_dp else []
    for prev, nxt in zip(usable_dp, usable_dp[1:]):
        if (prev.segment == nxt.segment
                and abs(prev.current_rate - nxt.current_rate) < 1.0e-9):
            rg.append(nxt)
        else:
            rate_groups.append(rg)
            rg = [nxt]
    if rg:
        rate_groups.append(rg)
    for group in rate_groups:
        early = group[:7]
        t0 = early[0].time
        pts = [c for c in early if c.time - t0 > 0.0]
        if len(pts) < 3:
            continue
        slope = _fit_slope(
            [math.log(c.time - t0) for c in pts],
            [math.log(c.delta_p) for c in pts],
        )
        if 0.6 <= slope <= 2.6:
            result["storage"] = (pts[0].time, pts[-1].time)
            break

    # Radial: longest flat run (pairwise |slope| <= 0.25) whose points stay
    # within +/-15% of their geometric mean.  The tightness rejects a WBS
    # hump or a boundary-doubling segment.
    best: list[ComputedSample] = []
    best_group_index = -1
    for gi, group in enumerate(groups):
        runs: list[list[ComputedSample]] = []
        run = [group[0]]
        for a, b in zip(group, group[1:]):
            dx = math.log(b.teq) - math.log(a.teq)
            dy = math.log(b.derivative) - math.log(a.derivative)
            flat = dx > 0 and abs(dy / dx) <= 0.25
            if flat:
                run.append(b)
            else:
                runs.append(run)
                run = [b]
        runs.append(run)
        for cand in runs:
            if len(cand) < 3:
                continue
            # A genuine radial run sits beyond the storage transition.
            if min(c.teq for c in cand if c.teq) < 6.0:
                continue
            gm = math.exp(sum(math.log(c.derivative) for c in cand) / len(cand))
            tight = max(c.derivative / gm for c in cand) < 1.15
            tight = tight and min(c.derivative / gm for c in cand) > 0.87
            if tight and len(cand) > len(best):
                best = cand
                best_group_index = gi
    if best:
        result["radial"] = (best[0].time, best[-1].time)

    # Boundary: a sustained rise AFTER the identified radial plateau.  The
    # plateau may live in the same group (the boundary then bends the tail of
    # that run up) or in a later constant-rate group.
    start_index = best_group_index if best_group_index >= 0 else 0
    radial_end_time = best[-1].time if best else -1.0
    for group in reversed(groups[start_index:]):
        if len(group) < 4:
            continue
        if group[-1].time <= radial_end_time:
            continue
        levels = [c.derivative for c in group]
        pre = [c for c in group if c.time <= radial_end_time] or group[:3]
        base = min(c.derivative for c in pre[:max(3, len(pre))])
        rising = [c for c in group
                  if c.time > radial_end_time and c.derivative / base >= 1.18]
        xs = [math.log(c.teq) for c in group]
        ys = [math.log(c.derivative) for c in group]
        trend = _fit_slope(xs, ys)
        if len(rising) >= 2 and trend > 0.05 and max(levels) / base >= 1.3:
            result["boundary"] = (rising[0].time, rising[-1].time)
            break
    return result
