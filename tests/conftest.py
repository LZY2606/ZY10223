import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from db import Database
import app as app_module


@pytest.fixture
def fixture():
    return json.loads((ROOT / "fixtures" / "well_test_fixture.json").read_text())


@pytest.fixture
def client(monkeypatch, tmp_path, fixture):
    from fastapi.testclient import TestClient

    db_path = str(tmp_path / "test.db")
    monkeypatch.setattr(app_module, "db", Database(db_path))
    with TestClient(app_module.app) as c:
        c.post("/api/import")
        yield c
