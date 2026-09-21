"""SQLite persistence for the drawdown interpretation platform.

Tables keep the *raw* readings immutable; every user action is an append-only
revision:

* ``pressure_samples`` - imported raw readings (never overwritten),
* ``rate_revisions``  - ordered (time, rate) corrections with a version,
* ``shift_revisions`` - marked/edited gauge shifts,
* ``schemes`` / ``intervals`` - side-by-side interpretation candidates,
* ``runs`` - frozen fingerprints + serialized state for replay/export.
"""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import time
from typing import Any, Optional

from .fixture import build_fixture

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS pressure_samples (
    idx INTEGER PRIMARY KEY,
    time REAL NOT NULL,
    pressure_raw REAL NOT NULL,
    duplicate INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE IF NOT EXISTS rate_revisions (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    time REAL NOT NULL,
    rate REAL NOT NULL,
    version INTEGER NOT NULL,
    note TEXT NOT NULL DEFAULT '',
    created_at REAL NOT NULL
);
CREATE TABLE IF NOT EXISTS shift_revisions (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    time REAL NOT NULL,
    offset_kpa REAL,
    version INTEGER NOT NULL,
    note TEXT NOT NULL DEFAULT '',
    created_at REAL NOT NULL
);
CREATE TABLE IF NOT EXISTS schemes (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    smoothing REAL NOT NULL DEFAULT 0.0,
    created_at REAL NOT NULL
);
CREATE TABLE IF NOT EXISTS intervals (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    scheme_id INTEGER NOT NULL REFERENCES schemes(id) ON DELETE CASCADE,
    regime TEXT NOT NULL,
    start_time REAL NOT NULL,
    end_time REAL NOT NULL,
    start_open INTEGER NOT NULL DEFAULT 0,
    end_open INTEGER NOT NULL DEFAULT 0,
    smoothing REAL NOT NULL DEFAULT 0.0,
    locked_plateau REAL,
    payload TEXT NOT NULL,
    created_at REAL NOT NULL
);
CREATE TABLE IF NOT EXISTS runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    label TEXT NOT NULL,
    fingerprint TEXT NOT NULL,
    rate_version INTEGER NOT NULL,
    payload TEXT NOT NULL,
    created_at REAL NOT NULL
);
"""


def connect(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA)
    conn.commit()


def is_seeded(conn: sqlite3.Connection) -> bool:
    row = conn.execute("SELECT COUNT(*) AS n FROM pressure_samples").fetchone()
    return row["n"] > 0


def seed_fixture(conn: sqlite3.Connection) -> dict:
    """Import the fixed fixture idempotently; never mutates existing rows."""
    init_db(conn)
    if is_seeded(conn):
        return get_state(conn)
    fx = build_fixture()
    rows = [
        (int(r["idx"]), float(r["time"]), float(r["pressure_raw"]),
         1 if r.get("duplicate") else 0)
        for r in fx["samples"]
    ]
    conn.executemany(
        "INSERT INTO pressure_samples(idx,time,pressure_raw,duplicate)"
        " VALUES (?,?,?,?)", rows)
    version = 1
    for st in fx["rate_steps"]:
        conn.execute(
            "INSERT INTO rate_revisions(time,rate,version,note,created_at)"
            " VALUES (?,?,?,'fixture',?)",
            (float(st["time"]), float(st["rate"]), version, time.time()))
    for sh in fx["gauge_shifts"]:
        conn.execute(
            "INSERT INTO shift_revisions(time,offset_kpa,version,note,created_at)"
            " VALUES (?,?,?,'fixture',?)",
            (float(sh["time"]), float(sh["offset_kpa"]), version, time.time()))
    conn.execute(
        "INSERT OR REPLACE INTO meta(key,value) VALUES('fixture_name',?)",
        (fx["name"],))
    conn.commit()
    return get_state(conn)


def rate_version(conn: sqlite3.Connection) -> int:
    row = conn.execute("SELECT COALESCE(MAX(version),0) AS v FROM rate_revisions").fetchone()
    return int(row["v"])


def current_rate_steps(conn: sqlite3.Connection) -> list[dict]:
    """Collapse revisions to the latest rate at each step time."""
    rows = conn.execute(
        "SELECT time, rate FROM rate_revisions r WHERE id = "
        "(SELECT id FROM rate_revisions WHERE time = r.time "
        " ORDER BY version DESC, id DESC LIMIT 1) ORDER BY time"
    ).fetchall()
    out = [dict(r) for r in rows]
    if not out or out[0]["time"] != 0.0:
        out.insert(0, {"time": 0.0, "rate": 0.0})
    return out


def current_shifts(conn: sqlite3.Connection) -> list[dict]:
    rows = conn.execute(
        "SELECT time, offset_kpa FROM shift_revisions r WHERE id = "
        "(SELECT id FROM shift_revisions WHERE time = r.time "
        " ORDER BY version DESC, id DESC LIMIT 1) ORDER BY time"
    ).fetchall()
    return [{"time": r["time"], "offset_kpa": r["offset_kpa"]} for r in rows]


def add_rate_revision(conn: sqlite3.Connection, t: float, rate: float,
                      note: str = "user") -> int:
    version = rate_version(conn) + 1
    conn.execute(
        "INSERT INTO rate_revisions(time,rate,version,note,created_at)"
        " VALUES (?,?,?,?,?)",
        (float(t), float(rate), version, note, time.time()))
    conn.commit()
    return version


def add_shift(conn: sqlite3.Connection, t: float,
              offset_kpa: Optional[float], note: str = "user") -> int:
    row = conn.execute("SELECT COALESCE(MAX(version),0) AS v FROM shift_revisions").fetchone()
    version = int(row["v"]) + 1
    conn.execute(
        "INSERT INTO shift_revisions(time,offset_kpa,version,note,created_at)"
        " VALUES (?,?,?,?,?)",
        (float(t), offset_kpa, version, note, time.time()))
    conn.commit()
    return version


def samples(conn: sqlite3.Connection) -> list[dict]:
    rows = conn.execute(
        "SELECT idx,time,pressure_raw,duplicate FROM pressure_samples"
        " ORDER BY time,idx").fetchall()
    return [dict(r) for r in rows]


def create_scheme(conn: sqlite3.Connection, name: str,
                  smoothing: float = 0.0) -> int:
    cur = conn.execute(
        "INSERT INTO schemes(name,smoothing,created_at) VALUES(?,?,?)",
        (name, float(smoothing), time.time()))
    conn.commit()
    return int(cur.lastrowid)


def add_interval(conn: sqlite3.Connection, scheme_id: int, payload: dict) -> int:
    cur = conn.execute(
        "INSERT INTO intervals(scheme_id,regime,start_time,end_time,start_open,"
        "end_open,smoothing,locked_plateau,payload,created_at)"
        " VALUES (?,?,?,?,?,?,?,?,?,?)",
        (scheme_id, payload["regime"], float(payload["start_time"]),
         float(payload["end_time"]), 1 if payload.get("start_open") else 0,
         1 if payload.get("end_open") else 0,
         float(payload.get("smoothing", 0.0)),
         payload.get("locked_plateau"), json.dumps(payload, ensure_ascii=False),
         time.time()))
    conn.commit()
    return int(cur.lastrowid)


def list_schemes(conn: sqlite3.Connection) -> list[dict]:
    schemes = [dict(r) for r in conn.execute(
        "SELECT * FROM schemes ORDER BY id").fetchall()]
    for s in schemes:
        iv = [dict(r) for r in conn.execute(
            "SELECT * FROM intervals WHERE scheme_id=? ORDER BY id",
            (s["id"],)).fetchall()]
        for i in iv:
            i["data"] = json.loads(i.pop("payload"))
        s["intervals"] = iv
    return schemes


def fingerprint(state: dict) -> str:
    """Stable hash including raw readings and the rate revision version."""
    canon = {
        "samples": [
            [s["idx"], round(float(s["time"]), 9),
             round(float(s["pressure_raw"]), 6)]
            for s in state["samples"]
        ],
        "steps": [
            [round(float(s["time"]), 9), round(float(s["rate"]), 9)]
            for s in state["steps"]
        ],
        "shifts": [
            [round(float(s["time"]), 9),
             None if s.get("offset_kpa") is None else round(float(s["offset_kpa"]), 6)]
            for s in state.get("shifts", [])
        ],
        "rate_version": state["rate_version"],
    }
    blob = json.dumps(canon, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


def save_run(conn: sqlite3.Connection, label: str, state: dict) -> dict:
    fp = fingerprint(state)
    cur = conn.execute(
        "INSERT INTO runs(label,fingerprint,rate_version,payload,created_at)"
        " VALUES (?,?,?,?,?)",
        (label, fp, int(state["rate_version"]),
         json.dumps(state, ensure_ascii=False), time.time()))
    conn.commit()
    return {"id": int(cur.lastrowid), "fingerprint": fp,
            "rate_version": state["rate_version"], "label": label}


def list_runs(conn: sqlite3.Connection) -> list[dict]:
    rows = conn.execute(
        "SELECT id,label,fingerprint,rate_version,created_at FROM runs"
        " ORDER BY id").fetchall()
    return [dict(r) for r in rows]


def get_run(conn: sqlite3.Connection, run_id: int) -> Optional[dict]:
    row = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
    if row is None:
        return None
    data = dict(row)
    data["state"] = json.loads(data.pop("payload"))
    return data


def get_state(conn: sqlite3.Connection) -> dict:
    return {
        "samples": samples(conn),
        "steps": current_rate_steps(conn),
        "shifts": current_shifts(conn),
        "rate_version": rate_version(conn),
    }


def clear_all(conn: sqlite3.Connection) -> None:
    for table in ("runs", "intervals", "schemes", "shift_revisions",
                  "rate_revisions", "pressure_samples", "meta"):
        conn.execute(f"DELETE FROM {table}")
    conn.commit()


def export_bundle(conn: sqlite3.Connection) -> dict:
    state = get_state(conn)
    return {
        "format": "drawdown-platform-v1",
        "exported_at": time.time(),
        "state": state,
        "fingerprint": fingerprint(state),
        "schemes": list_schemes(conn),
        "runs": list_runs(conn),
    }


def verify_bundle(bundle: dict) -> bool:
    state = bundle.get("state", {})
    return bundle.get("fingerprint") == fingerprint(state)


def import_bundle(conn: sqlite3.Connection, bundle: dict,
                  verify: bool = True) -> dict:
    if verify and not verify_bundle(bundle):
        raise ValueError("fingerprint mismatch; bundle is not verifiable")
    clear_all(conn)
    init_db(conn)
    state = bundle["state"]
    conn.executemany(
        "INSERT INTO pressure_samples(idx,time,pressure_raw,duplicate)"
        " VALUES (?,?,?,?)",
        [(int(s["idx"]), float(s["time"]), float(s["pressure_raw"]),
          int(s.get("duplicate", 0))) for s in state["samples"]])
    # Preserve the exported rate revision version so the fingerprint replays
    # bit-for-bit after a clear/import cycle.
    version = 0
    target_version = int(state.get("rate_version", len(state["steps"])))
    for st in state["steps"]:
        version += 1
        conn.execute(
            "INSERT INTO rate_revisions(time,rate,version,note,created_at)"
            " VALUES (?,?,?,'import',?)",
            (float(st["time"]), float(st["rate"]),
             target_version - len(state["steps"]) + version, time.time()))
    sh_version = 0
    for sh in state.get("shifts", []):
        sh_version += 1
        conn.execute(
            "INSERT INTO shift_revisions(time,offset_kpa,version,note,created_at)"
            " VALUES (?,?,?,'import',?)",
            (float(sh["time"]), sh.get("offset_kpa"), sh_version, time.time()))
    conn.commit()
    return get_state(conn)


DEFAULT_DB = os.environ.get(
    "DRAWDOWN_DB",
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                 "data", "drawdown.db"))
