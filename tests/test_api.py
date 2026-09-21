"""End-to-end FastAPI tests: revisions, schemes, runs, export/import/reset."""

from fastapi.testclient import TestClient

import app as app_module
from fastapi.testclient import TestClient as _TC

from welltest import store


def _client(temp_db):
    app_module.DB_PATH = temp_db
    return TestClient(app_module.app)


def test_index_shows_platform_title(temp_db):
    client = _client(temp_db)
    html = client.get("/").text
    assert "压降译井台" in html


def test_state_seeds_fixture_and_reports_diagnostics(temp_db):
    client = _client(temp_db)
    data = client.get("/api/state").json()
    assert data["rate_version"] >= 1
    assert data["diagnostics"]["duplicates"] == 2
    assert data["diagnostics"]["step_aligned"] == 2
    assert data["diagnostics"]["shift_segments"] == 2


def test_rate_revision_increments_version_and_changes_plateau(temp_db):
    client = _client(temp_db)
    before = client.get("/api/state").json()
    r = client.post("/api/rates", json={"time": 200, "rate": 120}).json()
    assert r["rate_version"] == before["rate_version"] + 1
    after = client.get("/api/state").json()
    step = next(s for s in after["steps"] if s["time"] == 200)
    assert step["rate"] == 120


def test_run_fingerprint_contains_rate_version(temp_db):
    client = _client(temp_db)
    run1 = client.post("/api/runs", json={"label": "v1"}).json()
    client.post("/api/rates", json={"time": 200, "rate": 120})
    run2 = client.post("/api/runs", json={"label": "v2"}).json()
    assert run1["fingerprint"] != run2["fingerprint"]
    assert run2["rate_version"] > run1["rate_version"]


def test_scheme_and_intervals_side_by_side(temp_db):
    client = _client(temp_db)
    s1 = client.post("/api/schemes", json={"name": "A"}).json()["scheme_id"]
    s2 = client.post("/api/schemes", json={"name": "B", "smoothing": 0.2}) \
        .json()["scheme_id"]
    body = {"regime": "radial", "start_time": 60, "end_time": 150}
    r = client.post(f"/api/schemes/{s1}/intervals", json=body).json()
    assert r["ok"] and r["fit"]["n_points"] >= 1
    schemes = client.get("/api/schemes").json()["schemes"]
    assert {s["id"] for s in schemes} == {s1, s2}


def test_export_import_roundtrip_verifies_and_replays(temp_db):
    client = _client(temp_db)
    client.post("/api/rates", json={"time": 200, "rate": 137})
    client.post("/api/runs", json={"label": "checkpoint"})
    bundle = client.get("/api/export").json()
    assert bundle["fingerprint"]
    # wipe via reset then re-import the exported bundle
    client.post("/api/reset")
    re = client.post("/api/import", json=bundle).json()
    assert re["fingerprint"] == bundle["fingerprint"]
    state = client.get("/api/state").json()
    step = next(s for s in state["steps"] if s["time"] == 200)
    assert step["rate"] == 137


def test_import_rejects_tampered_bundle(temp_db):
    client = _client(temp_db)
    bundle = client.get("/api/export").json()
    bundle["state"]["steps"][0]["rate"] = -999
    res = client.post("/api/import", json=bundle)
    assert res.status_code == 409


def test_raw_readings_immutable_after_revision(temp_db):
    client = _client(temp_db)
    original = client.get("/api/state").json()["samples"]
    client.post("/api/rates", json={"time": 200, "rate": 120})
    after = client.get("/api/state").json()["samples"]
    assert len(after) == len(original)
    assert all(a["pressure_raw"] == b["pressure_raw"]
               for a, b in zip(after, original))


def test_reset_reseeds_fixed_fixture(temp_db):
    client = _client(temp_db)
    client.post("/api/rates", json={"time": 200, "rate": 120})
    client.post("/api/reset")
    state = client.get("/api/state").json()
    step = next(s for s in state["steps"] if s["time"] == 200)
    assert step["rate"] == 200
