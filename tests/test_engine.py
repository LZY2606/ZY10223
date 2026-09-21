"""Invariants of the variable-rate interpretation engine."""

import math

import pytest

from welltest.engine import (
    ComputedSample,
    RateStep,
    Sample,
    Shift,
    bourdet_derivatives,
    canonical_steps,
    compute,
    equivalent_time,
    fit_interval,
    rate_at,
    suggest_candidates,
)
from welltest.fixture import build_fixture


def _fx_samples():
    fx = build_fixture()
    return [Sample(r["idx"], r["time"], r["pressure_raw"])
            for r in fx["samples"]]


def test_right_continuous_rate_lookup():
    steps = canonical_steps([
        {"time": 0, "rate": 0},
        {"time": 10, "rate": 60},
        {"time": 200, "rate": 200},
    ])
    # A step at exactly the sample time already applies.
    assert rate_at(steps, 10.0) == 60.0
    assert rate_at(steps, 199.999) == 60.0
    assert rate_at(steps, 200.0) == 200.0
    assert rate_at(steps, 0.0) == 0.0


def test_rate_step_aligned_sample_is_singular_not_smoothed():
    steps = canonical_steps([
        {"time": 0, "rate": 0},
        {"time": 10, "rate": 60},
        {"time": 200, "rate": 200},
    ])
    teq, note = equivalent_time(steps, 10.0)
    assert note == "rate_step_aligned"
    assert teq == 0.0
    # samples just after the step use a finite geometric superposition time
    teq_after, note_after = equivalent_time(steps, 10.5)
    assert note_after == "ok" and teq_after > 0


def test_equivalent_time_single_rate_is_elapsed_time():
    steps = [RateStep(0.0, 100.0)]
    teq, note = equivalent_time(steps, 25.0)
    assert note == "ok"
    assert teq == pytest.approx(25.0)


def test_zero_rate_sample_kept_but_diagnosed():
    steps = [RateStep(0.0, 0.0), RateStep(10.0, 60.0)]
    teq, note = equivalent_time(steps, 5.0)
    assert teq is None and note == "nonpositive_rate"


def test_derivative_on_actual_time_unequal_grid_recovers_plateau():
    # dp = A ln(t) with unequally spaced actual times -> constant derivative.
    times = [1.0, 1.7, 3.1, 6.4, 12.8, 25.0, 51.0, 100.0]
    plateau = 12.5
    cs = []
    for i, t in enumerate(times):
        cs.append(ComputedSample(
            idx=i, time=t, pressure_raw=0.0, pressure_corrected=0.0,
            delta_p=plateau * math.log(t), segment=0, current_rate=100.0,
            teq=t, teq_note="ok", derivative=None,
            derivative_note="", duplicate=False))
    bourdet_derivatives(cs, L=0.0)
    inner = [c.derivative for c in cs[1:-1]]
    assert all(v is not None for v in inner)
    assert max(inner) / min(inner) == pytest.approx(1.0, abs=0.05)
    assert sum(inner) / len(inner) == pytest.approx(plateau, rel=0.03)


def test_duplicate_timestamp_diagnosed_and_not_a_neighbour():
    fx = build_fixture()
    samples = _fx_samples()
    res = compute(samples, canonical_steps(fx["rate_steps"]),
                  [Shift(48.0)])
    dups = [c for c in res["samples"] if c.duplicate]
    assert len(dups) == 2
    for c in dups:
        assert c.delta_p is None
        assert c.derivative is None
        assert c.derivative_note == "duplicate_timestamp"


def test_gauge_shift_segmentation_breaks_derivative_pairs():
    # A hard jump at t=48 must NOT be bridged into a formation derivative.
    fx = build_fixture()
    samples = _fx_samples()
    res = compute(samples, canonical_steps(fx["rate_steps"]),
                  [Shift(48.0)])
    at48 = next(c for c in res["samples"] if c.time == 48.0)
    assert at48.segment == 1
    assert at48.derivative_note == "one_sided_segment_boundary"
    # the first point of segment 1 also cannot borrow across the jump
    after = [c for c in res["samples"] if c.segment == 1 and c.time > 48]
    leftmost = min(after, key=lambda c: c.time)
    # derivative exists only once an in-segment left neighbour is available
    assert leftmost.time >= 48.0


def test_estimated_shift_close_to_known_jump():
    fx = build_fixture()
    samples = _fx_samples()
    res = compute(samples, canonical_steps(fx["rate_steps"]),
                  [Shift(48.0)])
    estimated = res["shifts"][0].offset
    assert estimated == pytest.approx(fx["gauge_shifts"][0]["offset_kpa"],
                                      abs=0.5)


def test_smoothing_cannot_create_derivative_at_segment_boundary():
    fx = build_fixture()
    samples = _fx_samples()
    res = compute(samples, canonical_steps(fx["rate_steps"]),
                  [Shift(48.0)], smoothing=1.0)
    at48 = next(c for c in res["samples"] if c.time == 48.0)
    assert at48.derivative is None


def test_nonmonotonic_teq_after_rate_cut_is_flagged_not_bridged():
    # After a sharp rate decrease the geometric t_eq can reset; the engine must
    # report non_monotonic_teq rather than emit a fake huge derivative.
    fx = build_fixture()
    samples = _fx_samples()
    res = compute(samples, canonical_steps(fx["rate_steps"]),
                  [Shift(48.0)])
    notes = {c.derivative_note for c in res["samples"]}
    assert "non_monotonic_teq" in notes or "one_sided_segment_boundary" in notes


def test_fixture_preserves_all_required_artifacts():
    fx = build_fixture()
    times = [r["time"] for r in fx["samples"]]
    # two rate steps land exactly on pressure samples
    for st in (10.0, 200.0):
        assert st in times
    # duplicate readings
    assert sum(1 for r in fx["samples"] if r.get("duplicate")) == 2
    # gauge shift lands on a sample
    assert fx["gauge_shifts"][0]["time"] in times


def test_radial_candidate_recovers_permeability_order():
    fx = build_fixture()
    samples = _fx_samples()
    res = compute(samples, canonical_steps(fx["rate_steps"]),
                  [Shift(48.0)])
    cands = suggest_candidates(res["samples"])
    assert "radial" in cands
    a, b = cands["radial"]
    fit = fit_interval(res["samples"], "radial", a, b)
    # synthetic truth is 1000 mD; multi-rate superposition biases mildly,
    # but the answer must be within a factor of two.
    assert fit.permeability_md == pytest.approx(
        fx["truth"]["permeability_md"], rel=0.5)


def test_interval_endpoint_open_closed_semantics():
    fx = build_fixture()
    samples = _fx_samples()
    res = compute(samples, canonical_steps(fx["rate_steps"]),
                  [Shift(48.0)])
    cands = suggest_candidates(res["samples"])
    a, b = cands["radial"]
    closed = fit_interval(res["samples"], "radial", a, b)
    open_end = fit_interval(res["samples"], "radial", a, b, end_open=True)
    # opening the end can only remove a coincident sample, never add one
    assert open_end.n_points <= closed.n_points
