"""SQLite 持久层：原始读数永不覆盖，流量修订与运行方案独立版本化。"""
from __future__ import annotations

import json
import sqlite3
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, List, Optional

DEFAULT_DB = "data/pressure_transient.db"

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS raw_readings (
    id INTEGER PRIMARY KEY CHECK (id = 1),
    times_json TEXT NOT NULL,
    pressures_json TEXT NOT NULL,
    stages_json TEXT NOT NULL,
    reference_pressure REAL,
    imported_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS rate_revisions (
    version INTEGER PRIMARY KEY,
    stages_json TEXT NOT NULL,
    note TEXT,
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    rate_revision INTEGER NOT NULL,
    shifts_json TEXT NOT NULL,
    smooth_neighbors INTEGER NOT NULL,
    smooth_l REAL NOT NULL,
    intervals_json TEXT NOT NULL,
    samples_json TEXT NOT NULL,
    resolved_shifts_json TEXT NOT NULL,
    fingerprint TEXT NOT NULL,
    created_at TEXT NOT NULL
);
"""


class Database:
    def __init__(self, path: str = DEFAULT_DB):
        self.path = path
        if ":memory:" not in path:
            Path(path).parent.mkdir(parents=True, exist_ok=True)
        self.init_schema()

    @contextmanager
    def connect(self):
        conn = sqlite3.connect(self.path)
        conn.row_factory = sqlite3.Row
        try:
            yield conn
            conn.commit()
        finally:
            conn.close()

    def init_schema(self) -> None:
        with self.connect() as conn:
            conn.executescript(SCHEMA)

    def reset(self) -> None:
        with self.connect() as conn:
            conn.executescript(
                "DELETE FROM runs; DELETE FROM rate_revisions; "
                "DELETE FROM raw_readings; DELETE FROM meta;")

    # --- 原始数据 ---
    def has_raw(self) -> bool:
        with self.connect() as conn:
            return conn.execute("SELECT 1 FROM raw_readings WHERE id=1").fetchone() is not None

    def import_raw(self, fixture: Dict[str, Any], imported_at: str) -> None:
        with self.connect() as conn:
            conn.execute(
                "INSERT OR REPLACE INTO raw_readings "
                "(id, times_json, pressures_json, stages_json, reference_pressure, imported_at) "
                "VALUES (1,?,?,?,?,?)",
                (json.dumps(fixture["times"]), json.dumps(fixture["pressures"]),
                 json.dumps(fixture["stages"]), fixture.get("reference_pressure"), imported_at))
            conn.execute(
                "INSERT OR REPLACE INTO meta (key,value) VALUES ('fixture_version',?)",
                (fixture.get("fixture_version", ""),))

    def get_raw(self) -> Optional[Dict[str, Any]]:
        with self.connect() as conn:
            row = conn.execute("SELECT * FROM raw_readings WHERE id=1").fetchone()
        if not row:
            return None
        return {
            "times": json.loads(row["times_json"]),
            "pressures": json.loads(row["pressures_json"]),
            "stages": json.loads(row["stages_json"]),
            "reference_pressure": row["reference_pressure"],
        }

    # --- 流量修订 ---
    def list_revisions(self) -> List[Dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute("SELECT * FROM rate_revisions ORDER BY version").fetchall()
        return [{"version": r["version"], "stages": json.loads(r["stages_json"]),
                 "note": r["note"], "created_at": r["created_at"]} for r in rows]

    def add_revision(self, version: int, stages: List[Dict[str, Any]],
                     note: str, created_at: str) -> None:
        with self.connect() as conn:
            conn.execute(
                "INSERT INTO rate_revisions (version, stages_json, note, created_at) "
                "VALUES (?,?,?,?)",
                (version, json.dumps(stages), note, created_at))

    # --- 运行 ---
    def insert_run(self, run: Dict[str, Any]) -> int:
        with self.connect() as conn:
            cur = conn.execute(
                "INSERT INTO runs (name, rate_revision, shifts_json, smooth_neighbors, "
                "smooth_l, intervals_json, samples_json, resolved_shifts_json, fingerprint, "
                "created_at) VALUES (?,?,?,?,?,?,?,?,?,?)",
                (run["name"], run["rate_revision"], json.dumps(run["shifts"]),
                 run["smooth_neighbors"], run["smooth_l"], json.dumps(run["intervals"]),
                 json.dumps(run["samples"]), json.dumps(run["resolved_shifts"]),
                 run["fingerprint"], run["created_at"]))
            return int(cur.lastrowid)

    def list_runs(self) -> List[Dict[str, Any]]:
        with self.connect() as conn:
            rows = conn.execute(
                "SELECT id,name,rate_revision,fingerprint,created_at FROM runs ORDER BY id"
            ).fetchall()
        return [dict(r) for r in rows]

    def get_run(self, run_id: int) -> Optional[Dict[str, Any]]:
        with self.connect() as conn:
            row = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
        if not row:
            return None
        return {
            "id": row["id"], "name": row["name"],
            "rate_revision": row["rate_revision"],
            "shifts": json.loads(row["shifts_json"]),
            "smooth_neighbors": row["smooth_neighbors"],
            "smooth_l": row["smooth_l"],
            "intervals": json.loads(row["intervals_json"]),
            "samples": json.loads(row["samples_json"]),
            "resolved_shifts": json.loads(row["resolved_shifts_json"]),
            "fingerprint": row["fingerprint"], "created_at": row["created_at"],
        }
