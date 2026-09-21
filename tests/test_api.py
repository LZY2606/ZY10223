"""服务层验收：导入/清空/重导、流量修订进入指纹、导出-重放一致、页面标题。"""
import json


def test_index_title(client):
    r = client.get("/")
    assert r.status_code == 200
    assert "压降译井台" in r.text


def test_health_and_raw(client, fixture):
    h = client.get("/api/health").json()
    assert h["has_data"] is True
    raw = client.get("/api/raw").json()
    assert len(raw["times"]) == len(fixture["times"])
    assert raw["current_stages"][1]["rate"] == 15.0


def test_reset_and_reimport(client):
    assert client.post("/api/reset").json()["ok"] is True
    assert client.get("/api/health").json()["has_data"] is False
    # 清空后重新导入可复核
    d = client.post("/api/import").json()
    assert d["n"] > 0
    raw = client.get("/api/raw").json()
    assert len(raw["times"]) == d["n"]


def test_run_contains_candidates_and_fingerprint(client, fixture):
    idx = fixture["known_shifts"][0]["index"]
    r = client.post("/api/runs", json={"name": "基准", "shifts": [{"index": idx}]})
    assert r.status_code == 200
    run = r.json()
    assert run["fingerprint"]
    assert f"(nothing here)" not in run["fingerprint"]
    regimes = {iv["regime"] for iv in run["intervals"]}
    assert {"wbs", "radial", "boundary"}.issubset(regimes)
    # 换档被解析成段
    assert run["resolved_shifts"]
    # 样本含诊断标记
    flags = {f for s in run["samples"] for f in s["flags"]}
    assert "duplicate_time" in flags and "rate_boundary_sample" in flags


def test_rate_revision_changes_fingerprint(client, fixture):
    idx = fixture["known_shifts"][0]["index"]
    r1 = client.post("/api/runs", json={"shifts": [{"index": idx}]}).json()
    fp1 = r1["fingerprint"]

    # 修改一个流量阶段并保存为新版本
    stages = client.get("/api/raw").json()["current_stages"]
    stages[1]["rate"] = 16.0
    rev = client.post("/api/revisions", json={"stages": stages, "note": "调参"}).json()
    assert rev["version"] >= 1

    r2 = client.post("/api/runs", json={"shifts": [{"index": idx}]}).json()
    assert r2["rate_revision"] != r1["rate_revision"]
    assert r2["fingerprint"] != fp1
    # 原始读数未被覆盖
    raw = client.get("/api/raw").json()
    fx = fixture
    assert raw["times"] == fx["times"] and raw["pressures"] == fx["pressures"]


def test_interval_lock_and_derived(client, fixture):
    idx = fixture["known_shifts"][0]["index"]
    run = client.post("/api/runs", json={"shifts": [{"index": idx}]}).json()
    radial = next(iv for iv in run["intervals"] if iv["regime"] == "radial")
    radial["locked_slope"] = 0.0
    r = client.put(f"/api/runs/{run['id']}/intervals", json=run["intervals"])
    assert r.status_code == 200
    updated = r.json()
    res = next(x for x in updated["interval_results"] if x["regime"] == "radial")
    assert res["locked"] is True
    assert abs(res["slope_loglog"]) < 1e-9


def test_export_replay_matches(client, fixture):
    idx = fixture["known_shifts"][0]["index"]
    run = client.post("/api/runs", json={"shifts": [{"index": idx}]}).json()
    export = client.get(f"/api/runs/{run['id']}/export").json()
    replay = client.post("/api/runs/replay", json=export).json()
    assert replay["matches"] is True
    assert replay["fingerprint"] == run["fingerprint"]
    # 导出包携带流量修订记录（指纹含修订版本）
    assert "rate_revisions" in export


def test_multiple_runs_side_by_side(client, fixture):
    idx = fixture["known_shifts"][0]["index"]
    a = client.post("/api/runs", json={"name": "A", "shifts": [{"index": idx}]}).json()
    b = client.post("/api/runs", json={"name": "B", "shifts": []}).json()
    assert a["id"] != b["id"]
    listing = client.get("/api/runs").json()["runs"]
    assert len(listing) == 2
    # 两个方案互不覆盖
    got_a = client.get(f"/api/runs/{a['id']}").json()
    got_b = client.get(f"/api/runs/{b['id']}").json()
    assert got_a["shifts"] and not got_b["shifts"]


def test_409_without_data(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    import app as app_module
    from db import Database
    monkeypatch.setattr(app_module, "db", Database(str(tmp_path / "empty.db")))
    with TestClient(app_module.app) as c:
        assert c.get("/api/raw").status_code == 409
