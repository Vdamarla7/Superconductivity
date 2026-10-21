import pytest


@pytest.fixture(autouse=True)
def no_network_sharing(monkeypatch, tmp_path):
    """Tests never contact GitHub; test_share.py re-enables sharing with a fake API."""
    monkeypatch.setenv("SCFLOW_SHARE", "0")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    monkeypatch.delenv("GITHUB_TOKEN", raising=False)
    monkeypatch.delenv("GH_TOKEN", raising=False)
