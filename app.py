"""压降译井台 - FastAPI 本地服务。"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, Field

import analysis
from db import DEFAULT_DB, Database

BASE_DIR = Path(__file__).resolve().parent
FIXTURE_PATH = BASE_DIR / "fixtures" / "well_test_fixture.json"

app = FastAPI(title="压降译井台", version="1.0.0")
db = Database(DEFAULT_DB)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _load_fixture() -> Dict[str, Any]:
    return json.loads(FIXTURE_PATH.read_text())


def _ensure_raw() -> Dict[str, Any]:
    raw = db.get_raw()
    if raw is None:
        raise HTTPException(409, "尚未导入压力数据；请先调用 /api/import")
    return raw


def _current_stages() -> List[Dict[str, Any]]:
    revs = db.list_revisions()
    if revs:
        return revs[-1]["stages"]
    return _ensure_raw()["stages"]


def _data_digest(times: List[float], pressures: List[float]) -> str:
    h = hashlib.sha256()
    h.update(json.dumps([times, pressures], separators=(",", ":")).encode())
    return h.hexdigest()[:16]


class ShiftIn(BaseModel):
    index: int
    offset: Optional[float] = None


class CreateRunIn(BaseModel):
    name: str = "未命名方案"
    shifts: List[ShiftIn] = Field(default_factory=list)
    smooth_neighbors: int = analysis.DEFAULT_SMOOTH_NEIGHBORS
    smooth_l: float = analysis.DEFAULT_SMOOTH_L
    auto_candidates: bool = True


class IntervalIn(BaseModel):
    regime: str
    left: float
    right: float
    left_open: bool = False
    right_open: bool = False
    segments: Optional[List[int]] = None
    rates: Optional[List[float]] = None
    sample_indices: Optional[List[int]] = None
    locked_slope: Optional[float] = None
    smooth_l: Optional[float] = None


class RevisionIn(BaseModel):
    stages: List[Dict[str, Any]]
    note: str = ""


def _serialize_run(run: Dict[str, Any], include_samples: bool = True) -> Dict[str, Any]:
    out = {k: v for k, v in run.items() if include_samples or k != "samples"}
    out["interval_results"] = [
        analysis.derive_interval(run["samples"], iv) for iv in run["intervals"]
    ]
    return out


@app.get("/")
def index():
    return FileResponse(BASE_DIR / "static" / "index.html")


@app.get("/api/health")
def health():
    return {"status": "ok", "has_data": db.has_raw(),
            "analysis_version": analysis.ANALYSIS_VERSION}


@app.post("/api/import")
def import_fixture():
    """从固定 fixture 重新导入原始读数（清空旧库后写入，便于复核）。"""
    fx = _load_fixture()
    db.reset()
    db.import_raw(fx, _now())
    return {"ok": True, "n": len(fx["times"]),
            "fixture_version": fx.get("fixture_version"),
            "known_shifts": fx.get("known_shifts", [])}


@app.post("/api/reset")
def reset_database():
    db.reset()
    return {"ok": True}


@app.get("/api/raw")
def get_raw():
    raw = _ensure_raw()
    fx = _load_fixture()
    return {
        **raw,
        "units": fx.get("units"),
        "reservoir": fx.get("reservoir"),
        "known_shifts": fx.get("known_shifts", []),
        "revisions": db.list_revisions(),
        "current_stages": _current_stages(),
    }


@app.post("/api/revisions")
def create_revision(body: RevisionIn):
    _ensure_raw()
    stages = analysis.normalize_stages(body.stages)
    version = analysis.next_rate_revision(db.list_revisions())
    db.add_revision(version, stages, body.note, _now())
    return {"version": version, "stages": stages}


@app.get("/api/revisions")
def list_revisions():
    return {"revisions": db.list_revisions()}


@app.get("/api/runs")
def list_runs():
    raw = _ensure_raw()
    stages = _current_stages()
    revs = db.list_revisions()
    return {"runs": db.list_runs(),
            "rate_revision": revs[-1]["version"] if revs else 0,
            "stages": stages}


@app.post("/api/runs")
def create_run(body: CreateRunIn):
    raw = _ensure_raw()
    revs = db.list_revisions()
    revision = revs[-1]["version"] if revs else 0
    stages = revs[-1]["stages"] if revs else raw["stages"]
    shifts = [s.model_dump() for s in body.shifts]

    computed = analysis.build_samples(
        stages, raw["times"], raw["pressures"], shifts=shifts,
        smooth_l=body.smooth_l, smooth_neighbors=body.smooth_neighbors,
        reference_pressure=raw["reference_pressure"])

    intervals = analysis.detect_candidates(computed["samples"]) if body.auto_candidates else []

    fingerprint = analysis.run_fingerprint({
        "rate_revision": revision,
        "stages": stages,
        "shifts": shifts,
        "smooth_l": body.smooth_l,
        "smooth_neighbors": body.smooth_neighbors,
        "reference_pressure": computed["reference_pressure"],
        "intervals": intervals,
        "data_digest": _data_digest(raw["times"], raw["pressures"]),
    })
    run = {
        "name": body.name, "rate_revision": revision, "shifts": shifts,
        "smooth_neighbors": body.smooth_neighbors, "smooth_l": body.smooth_l,
        "intervals": intervals, "samples": computed["samples"],
        "resolved_shifts": computed["resolved_shifts"],
        "fingerprint": fingerprint, "created_at": _now(),
    }
    run_id = db.insert_run(run)
    stored = db.get_run(run_id)
    return _serialize_run(stored)


@app.get("/api/runs/{run_id}")
def get_run(run_id: int):
    run = db.get_run(run_id)
    if not run:
        raise HTTPException(404, "运行方案不存在")
    return _serialize_run(run)


@app.put("/api/runs/{run_id}/intervals")
def update_intervals(run_id: int, intervals: List[IntervalIn]):
    run = db.get_run(run_id)
    if not run:
        raise HTTPException(404, "运行方案不存在")
    normalized = []
    for iv in intervals:
        d = iv.model_dump()
        d.setdefault("smooth_l", run["smooth_l"])
        normalized.append(d)
    run["intervals"] = normalized
    fingerprint = analysis.run_fingerprint({
        "rate_revision": run["rate_revision"],
        "stages": _current_stages(),
        "shifts": run["shifts"],
        "smooth_l": run["smooth_l"],
        "smooth_neighbors": run["smooth_neighbors"],
        "reference_pressure": run["samples"][0]["pressure"] if run["samples"] else None,
        "intervals": normalized,
        "data_digest": _data_digest(_ensure_raw()["times"], _ensure_raw()["pressures"]),
    })
    # 直接更新该行
    with db.connect() as conn:
        conn.execute(
            "UPDATE runs SET intervals_json=?, fingerprint=? WHERE id=?",
            (json.dumps(normalized), fingerprint, run_id))
    stored = db.get_run(run_id)
    return _serialize_run(stored)


@app.get("/api/runs/{run_id}/export")
def export_run(run_id: int):
    run = db.get_run(run_id)
    if not run:
        raise HTTPException(404, "运行方案不存在")
    raw = _ensure_raw()
    payload = {
        "format": "pwgsb-run-export-1.0",
        "exported_at": _now(),
        "raw": raw,
        "rate_revisions": db.list_revisions(),
        "run": _serialize_run(run),
    }
    return JSONResponse(payload)


@app.post("/api/runs/replay")
def replay_run(payload: Dict[str, Any]):
    """从导出包重放：用其中的流量版本/换档/平滑重新计算并核对指纹。

    不写库；返回重算结果与指纹一致性。也会用当前库的原始读数（若已导入）交叉核对。
    """
    if payload.get("format") != "pwgsb-run-export-1.0":
        raise HTTPException(400, "导出包格式不支持")
    src_run = payload["run"]
    rev_version = src_run["rate_revision"]
    revisions = payload.get("rate_revisions", [])
    match = [r for r in revisions if r["version"] == rev_version]
    stages = match[0]["stages"] if match else payload["raw"]["stages"]
    raw = payload["raw"]

    recomputed = analysis.build_samples(
        stages, raw["times"], raw["pressures"], shifts=src_run["shifts"],
        smooth_l=src_run["smooth_l"], smooth_neighbors=src_run["smooth_neighbors"],
        reference_pressure=raw["reference_pressure"])
    fingerprint = analysis.run_fingerprint({
        "rate_revision": rev_version,
        "stages": stages,
        "shifts": src_run["shifts"],
        "smooth_l": src_run["smooth_l"],
        "smooth_neighbors": src_run["smooth_neighbors"],
        "reference_pressure": raw["reference_pressure"],
        "intervals": src_run["intervals"],
        "data_digest": _data_digest(raw["times"], raw["pressures"]),
    })
    return {
        "fingerprint": fingerprint,
        "original_fingerprint": src_run["fingerprint"],
        "matches": fingerprint == src_run["fingerprint"],
        "resolved_shifts": recomputed["resolved_shifts"],
        "samples": recomputed["samples"],
        "interval_results": [
            analysis.derive_interval(recomputed["samples"], iv)
            for iv in src_run["intervals"]
        ],
    }
