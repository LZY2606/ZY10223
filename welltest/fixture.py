"""Fixed synthetic fixture for the drawdown interpretation platform.

Deterministic (no RNG).  Deliberately contains:

* two short duplicate readings (0.020 h and 150.0 h),
* one gauge shift (+8.000 kPa) starting exactly at the t = 48 h sample,
* two rate steps (10 h: 0 -> 60 m3/d and 200 h: 60 -> 200 m3/d) landing
  exactly on pressure samples, so the right-continuous rule is testable,
* an extra rate change 200 -> 80 at 2400 h which does NOT land on a sample,
* a single sealing fault at 0.8 m.

Formation response: the line-source (E1) infinite-acting well convolved
with wellbore storage through C_D dpD/dtD = pD_rad(tD) - pD(tD) (exact unit
slope early, exact radial plateau late), plus a sealing-fault image well.
"""

from __future__ import annotations

import math

import numpy as np

from .engine import (
    FVF,
    FORMATION_THICKNESS,
    POROSITY,
    TOTAL_COMPRESSIBILITY,
    VISCOSITY,
    WELL_RADIUS,
)

INITIAL_PRESSURE_KPA = 25000.0
GAUGE_SHIFT_KPA = 8.0
GAUGE_SHIFT_TIME_H = 48.0

TRUE = {
    "permeability_md": 1000.0,
    "thickness_m": FORMATION_THICKNESS,
    "viscosity_pas": VISCOSITY,
    "porosity": POROSITY,
    "ct_pa_inv": TOTAL_COMPRESSIBILITY,
    "well_radius_m": WELL_RADIUS,
    "cd": 5.0,
    "skin": 0.0,
    "fault_distance_m": 1.5,
}

RATE_STEPS = [
    {"time": 0.0, "rate": 0.0},
    {"time": 10.0, "rate": 60.0},
    {"time": 200.0, "rate": 200.0},
    {"time": 2400.0, "rate": 80.0},
]

_GRID = np.array(
    [0.001, 0.002, 0.005, 0.010, 0.020, 0.050, 0.100, 0.200, 0.500,
     1.0, 2.0, 5.0, 10.0, 10.5, 11.0, 12.0, 13.0, 15.0, 18.0, 22.0,
     28.0, 34.0, 40.0, 48.0, 60.0, 70.0, 80.0, 92.0, 105.0, 120.0,
     138.0, 150.0, 158.0, 178.0, 200.0, 260.0,
     340.0, 380.0, 430.0, 470.0, 510.0, 555.0, 600.0, 660.0, 720.0,
     860.0, 1000.0, 1200.0, 1450.0, 1750.0, 2100.0, 2300.0, 2500.0,
     2700.0, 3000.0]
)

_DUPLICATES = {0.020: 0.11, 150.0: -0.07}


def _e1(z: float) -> float:
    """Exponential integral E1(z) = -Ei(-z), z > 0."""
    if z >= 6.0:
        # E1(z) ~ e^-z/z * [1 - 1/z + 2!/z^2 - 3!/z^3 + ...] (truncated
        # at the least term, which is the optimal asymptotic stopping point).
        total = 1.0
        term = 1.0
        for k in range(1, 100):
            nxt = -term * k / z
            if abs(nxt) >= abs(term):
                break
            term = nxt
            total += term
        return math.exp(-z) / z * total
    gamma = 0.5772156649015329
    total = -gamma - math.log(z)
    power = -z          # (-z)^1
    factorial = 1.0
    for k in range(1, 200):
        total -= power / (k * factorial)
        power *= -z
        factorial *= k + 1
        if abs(power / ((k + 1) * factorial)) < 1.0e-16:
            break
    return total


def _radial_pd(td: float) -> float:
    return 0.5 * _e1(0.25 / td)


def _dim_pressure(td: float, cd: float) -> float:
    """Wellbore-storage convolution at a single dimensionless time.

    Pre-computed on a fixed log grid in ``build_fixture``; this scalar form
    exists for tests/inspection and solves the ODE on a local deterministic
    grid each call (fixture generation uses the vector version for speed).
    """
    return float(_wbs_curve(td, cd))


def _wbs_grid(td_max: float, cd: float, n: int = 4000) -> tuple[np.ndarray, np.ndarray]:
    """Integrate CD dpD/dtD = pD_rad(tD) - pD on a log-spaced grid."""
    tg = np.logspace(-6.0, math.log10(max(td_max, 1.0e-5)), n)
    p = np.zeros(n)
    for i in range(1, n):
        dt = tg[i] - tg[i - 1]
        r = _radial_pd(float(tg[i]))
        p[i] = p[i - 1] + (r - p[i - 1]) * (1.0 - math.exp(-dt / cd))
    return tg, p


def _wbs_curve(td: float, cd: float) -> float:
    tg, p = _wbs_grid(max(td * 1.01, 1.0e-5), cd)
    return float(np.interp(td, tg, p))


def build_fixture() -> dict:
    k_si = TRUE["permeability_md"] / 1.01325e15
    phi, mu = POROSITY, VISCOSITY
    ct, rw, h = TOTAL_COMPRESSIBILITY, WELL_RADIUS, FORMATION_THICKNESS
    b = FVF
    cd = TRUE["cd"]
    fault_l = TRUE["fault_distance_m"]

    alpha = k_si / (phi * mu * ct)                  # m^2/s
    tcoef = alpha / rw ** 2 / 3600.0               # dimensionless t per hour
    fault_td = (2.0 * fault_l / (2.0 * rw)) ** 2
    # p = mu B q_d / (2 pi k h) * pD with q_d = q[m3/d]/86400; the Pa result
    # is divided by 1000 to give kPa per (m3/d) per unit pD.
    m0 = mu * b / (2.0 * math.pi * k_si * h) / 86400.0 / 1000.0

    steps = [(0.0, 0.0), (10.0, 60.0), (200.0, 200.0), (2400.0, 80.0)]
    td_max = tcoef * float(_GRID.max())
    tg_wbs, pg_rad = _wbs_grid(td_max, cd, n=3000)

    # WBS convolution of the total (radial + sealing-fault image) kernel.
    p_tot = np.zeros_like(tg_wbs)
    for i in range(1, len(tg_wbs)):
        dt = tg_wbs[i] - tg_wbs[i - 1]
        kernel = _radial_pd(float(tg_wbs[i])) + 0.5 * _e1(fault_td / tg_wbs[i])
        p_tot[i] = p_tot[i - 1] + (kernel - p_tot[i - 1]) * (
            1.0 - math.exp(-dt / cd))
    pg_fault = p_tot - pg_rad

    def drawdown(t: float) -> float:
        total = 0.0
        q_prev = 0.0
        for tk, qk in steps:
            if t <= tk:
                break
            dq = qk - q_prev
            td = tcoef * (t - tk)
            if td > 0.0 and dq != 0.0:
                pd_r = float(np.interp(td, tg_wbs, pg_rad))
                pd_f = float(np.interp(td, tg_wbs, pg_fault))
                pd = pd_r + pd_f
                total += dq * m0 * pd
            q_prev = qk
        return total

    rows = []
    idx = 0
    for t in _GRID:
        dp = drawdown(float(t))
        raw = INITIAL_PRESSURE_KPA - dp
        if t >= GAUGE_SHIFT_TIME_H:
            raw += GAUGE_SHIFT_KPA
        rows.append({"idx": idx, "time": float(t),
                     "pressure_raw": round(raw, 3)})
        idx += 1
        if float(t) in _DUPLICATES:
            rows.append({"idx": idx, "time": float(t),
                         "pressure_raw": round(raw + _DUPLICATES[float(t)], 3),
                         "duplicate": True})
            idx += 1

    return {
        "name": "ZY10223-drawdown",
        "description": (
            "Synthetic multi-rate drawdown: duplicate short readings, "
            "gauge shift at 48 h, rate steps at 10/200 h (on samples) and "
            "2400 h (off grid), sealing fault at 0.8 m."
        ),
        "units": {"time": "h", "rate": "m3/d", "pressure": "kPa"},
        "initial_pressure_kpa": INITIAL_PRESSURE_KPA,
        "rate_steps": RATE_STEPS,
        "gauge_shifts": [{"time": GAUGE_SHIFT_TIME_H,
                          "offset_kpa": GAUGE_SHIFT_KPA}],
        "truth": TRUE,
        "samples": rows,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(build_fixture(), ensure_ascii=False, indent=2))
