import os
import tempfile

import pytest


@pytest.fixture()
def temp_db(monkeypatch):
    fd, path = tempfile.mkstemp(suffix=".db")
    os.close(fd)
    monkeypatch.setenv("DRAWDOWN_DB", path)
    yield path
    if os.path.exists(path):
        os.remove(path)
