"""FastAPI application for the drawdown interpretation platform.

Run:
    .venv/bin/uvicorn app:app --host 127.0.0.1 --port 5563

Database path defaults to data/drawdown.db and can be overridden with the
DRAWDOWN_DB environment variable (used by the test suite).
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict
from typing import Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from welltest import engine, store

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
STATIC_DIR = os.path.join(BASE_DIR, "static")
DB_PATH = os.environ.get(
    "DRAWDOWN_DB", os.path.join(BASE_DIR, "data", "drawdown.db"))

app = FastAPI(title="drawdown-platform", version="1.0.0")


def get_conn():
    os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
    conn = store.connect(DB_PATH)
    store.init_db(conn)
    if not store.is_seeded(conn):
        store.seed_fixture(conn)
    return conn


def _sample_dict(c: engine.ComputedSample) -> dict:
    return {
        "idx": c.idx,
        "time": c.time,
        "pressure_raw": c.pressure_raw,
        "pressure_corrected": c.pressure_corrected,
        "delta_p": c.delta_p,
        "segment": c.segment,
        "current_rate": c.current_rate,
        "teq": c.teq,
        "teq_note": c.teq_note,
        "derivative": c.derivative,
        "derivative_note": c.derivative_note,
        "duplicate": c.duplicate,
    }


def _computed_from_dicts(rows):
    return [
        engine.ComputedSample(
            idx=c["idx"], time=c["time"], pressure_raw=c["pressure_raw"],
            pressure_corrected=c["pressure_corrected"], delta_p=c["delta_p"],
            segment=c["segment"], current_rate=c["current_rate"],
            teq=c["teq"], teq_note=c["teq_note"], derivative=c["derivative"],
            derivative_note=c["derivative_note"], duplicate=c["duplicate"])
        for c in rows
    ]


def compute_state(smoothing: float = 0.0) -> dict:
    conn = get_conn()
    try:
        state = store.get_state(conn)
        samples = [engine.Sample(int(r["idx"]), float(r["time"]),
                                 float(r["pressure_raw"]))
                   for r in state["samples"]]
        steps = engine.canonical_steps(state["steps"])
        shifts_in = [
            engine.Shift(float(s["time"]),
                         0.0 if s.get("offset_kpa") is None
                         else float(s["offset_kpa"]))
            for s in state["shifts"]
        ]
        result = engine.compute(samples, steps, shifts_in or None,
                                smoothing=smoothing)
        return {
            "rate_version": state["rate_version"],
            "smoothing": smoothing,
            "diagnostics": result["diagnostics"],
            "steps": [{"time": s.time, "rate": s.rate}
                      for s in result["steps"]],
            "shifts": [{"time": s.time, "offset_kpa": s.offset}
                       for s in result["shifts"]],
            "samples": [_sample_dict(c) for c in result["samples"]],
            "candidates": engine.suggest_candidates(result["samples"],
                                                    smoothing),
            "schemes": store.list_schemes(conn),
        }
    finally:
        conn.close()


class RateRevisionIn(BaseModel):
    time: float
    rate: float
    note: Optional[str] = "user"


class ShiftIn(BaseModel):
    time: float
    offset_kpa: Optional[float] = None
    note: Optional[str] = "user"


class SchemeIn(BaseModel):
    name: str
    smoothing: float = 0.0


class IntervalIn(BaseModel):
    regime: str
    start_time: float
    end_time: float
    start_open: bool = False
    end_open: bool = False
    smoothing: float = 0.0
    locked_plateau: Optional[float] = None


class RunIn(BaseModel):
    label: str


@app.get("/")
def index():
    return FileResponse(os.path.join(STATIC_DIR, "index.html"))


@app.get("/api/state")
def api_state(smoothing: float = 0.0):
    return compute_state(smoothing)


@app.post("/api/rates")
def api_revise_rate(body: RateRevisionIn):
    conn = get_conn()
    try:
        version = store.add_rate_revision(
            conn, body.time, body.rate, body.note or "user")
        return {"ok": True, "rate_version": version}
    finally:
        conn.close()


@app.post("/api/shifts")
def api_mark_shift(body: ShiftIn):
    conn = get_conn()
    try:
        version = store.add_shift(
            conn, body.time, body.offset_kpa, body.note or "user")
        return {"ok": True, "shift_version": version}
    finally:
        conn.close()


@app.post("/api/schemes")
def api_create_scheme(body: SchemeIn):
    conn = get_conn()
    try:
        sid = store.create_scheme(conn, body.name, body.smoothing)
        return {"ok": True, "scheme_id": sid}
    finally:
        conn.close()


@app.get("/api/schemes")
def api_schemes():
    conn = get_conn()
    try:
        return {"schemes": store.list_schemes(conn)}
    finally:
        conn.close()


@app.post("/api/schemes/{scheme_id}/intervals")
def api_add_interval(scheme_id: int, body: IntervalIn,
                     smoothing: float = 0.0):
    state = compute_state(smoothing)
    cs = _computed_from_dicts(state["samples"])
    fit = engine.fit_interval(
        cs, body.regime, body.start_time, body.end_time,
        start_open=body.start_open, end_open=body.end_open,
        smoothing=body.smoothing, plateau_override=body.locked_plateau)
    payload = body.model_dump()
    payload["fit"] = asdict(fit)
    conn = get_conn()
    try:
        exists = conn.execute(
            "SELECT id FROM schemes WHERE id=?", (scheme_id,)).fetchone()
        if exists is None:
            raise HTTPException(404, "scheme not found")
        iid = store.add_interval(conn, scheme_id, payload)
        return {"ok": True, "interval_id": iid, "fit": payload["fit"]}
    finally:
        conn.close()


@app.post("/api/runs")
def api_save_run(body: RunIn, smoothing: float = 0.0):
    state = compute_state(smoothing)
    conn = get_conn()
    try:
        return store.save_run(conn, body.label, state)
    finally:
        conn.close()


@app.get("/api/runs")
def api_runs():
    conn = get_conn()
    try:
        return {"runs": store.list_runs(conn)}
    finally:
        conn.close()


@app.get("/api/runs/{run_id}")
def api_run(run_id: int):
    conn = get_conn()
    try:
        run = store.get_run(conn, run_id)
        if run is None:
            raise HTTPException(404, "run not found")
        return run
    finally:
        conn.close()


@app.get("/api/export")
def api_export():
    conn = get_conn()
    try:
        return store.export_bundle(conn)
    finally:
        conn.close()


@app.post("/api/import")
async def api_import(request: Request):
    try:
        bundle = await request.json()
    except Exception:
        raise HTTPException(400, "invalid JSON bundle")
    conn = get_conn()
    try:
        state = store.import_bundle(conn, bundle)
        return {"ok": True, "fingerprint": store.fingerprint(state),
                "rate_version": state["rate_version"]}
    except ValueError as exc:
        raise HTTPException(409, str(exc))
    finally:
        conn.close()


@app.post("/api/reset")
def api_reset():
    conn = get_conn()
    try:
        store.clear_all(conn)
        store.seed_fixture(conn)
        return {"ok": True}
    finally:
        conn.close()


app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
