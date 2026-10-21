import json
from pathlib import Path

import pytest
import yaml

from scflow import db, share
from scflow.validation import validate

ROOT = Path(__file__).resolve().parents[1]
FIX = Path(__file__).parent / "fixtures"


@pytest.fixture
def record(tmp_path):
    cfg = yaml.safe_load((ROOT / "compounds/Nb_0GPa.yaml").read_text())
    cfg.setdefault("parallel", {})["nproc"] = 1
    con = db.connect(tmp_path / "t.db")
    rid = db.record_run(con, cfg, FIX / "Nb", validate(cfg, FIX / "Nb"), host="my-laptop.local")
    return db.export_record(con, rid)


class FakeGitHub:
    def __init__(self, login="someone", fail_on=None):
        self.login, self.fail_on, self.calls = login, fail_on, []

    def __call__(self, method, path, token, body=None):
        self.calls.append((method, path, body))
        if self.fail_on and self.fail_on in path:
            raise share.ShareError("boom")
        if path == "/user":
            return {"login": self.login}
        if path.endswith("/forks"):
            return {"full_name": f"{self.login}/scflow-results"}
        if "/git/ref/heads/main" in path:
            return {"object": {"sha": "abc123"}}
        if path.endswith("/pulls"):
            return {"html_url": "https://github.com/Vdamarla7/scflow-results/pull/1"}
        return {}


def test_submission_drops_host_and_paths(record):
    sub = share.make_submission(record)
    text = json.dumps(sub)
    assert "my-laptop" not in text and "run_dir" not in sub["run"] and "host" not in sub["run"]
    assert "/home/" not in text and "/Users/" not in text
    assert sub["submission"]["verified"] is False
    assert share.repo_path(sub) == f"records/Nb/0GPa/{sub['run']['id']}.json"


def test_opt_out_precedence(monkeypatch, tmp_path):
    monkeypatch.delenv("SCFLOW_SHARE")
    assert share.share_enabled() is True
    cfg = share.config_dir() / "config.yaml"
    cfg.parent.mkdir(parents=True)
    cfg.write_text("share: false\n")
    assert share.share_enabled() is False
    monkeypatch.setenv("SCFLOW_SHARE", "1")
    assert share.share_enabled() is True          # env beats config file
    assert share.share_enabled(no_share=True) is False  # flag beats env


def test_no_credentials_queues_in_outbox(record, tmp_path, capsys):
    share.share_record(record, tmp_path / "runs")
    files = share.pending(tmp_path / "runs")
    assert len(files) == 1
    assert "scflow share" in capsys.readouterr().out


def test_disabled_writes_nothing(record, tmp_path):
    share.share_record(record, tmp_path / "runs", enabled=False)
    assert share.pending(tmp_path / "runs") == []


def test_send_opens_pr_from_fork(record, tmp_path):
    share.write_outbox(share.make_submission(record), tmp_path)
    gh = FakeGitHub()
    url = share.send(share.pending(tmp_path), token="t", api=gh)
    assert url.endswith("/pull/1")
    methods = [(m, p.split("/")[-1]) for m, p, _ in gh.calls]
    assert ("POST", "forks") in methods
    put = [c for c in gh.calls if c[0] == "PUT"][0]
    assert put[1].startswith("/repos/someone/scflow-results/contents/records/Nb/0GPa/")
    pr = gh.calls[-1][2]
    assert pr["head"].startswith("someone:scflow-") and pr["base"] == "main"
    assert share.pending(tmp_path) == [] and len(list((tmp_path / "outbox/sent").glob("*.json"))) == 1


def test_owner_pushes_branch_without_fork(record, tmp_path):
    share.write_outbox(share.make_submission(record), tmp_path)
    gh = FakeGitHub(login="Vdamarla7")
    share.send(share.pending(tmp_path), token="t", repo="Vdamarla7/scflow-results", api=gh)
    assert not any(p.endswith("/forks") for _, p, _ in gh.calls)


def test_failed_send_keeps_outbox(record, tmp_path):
    share.write_outbox(share.make_submission(record), tmp_path)
    with pytest.raises(share.ShareError):
        share.send(share.pending(tmp_path), token="t", api=FakeGitHub(fail_on="/pulls"))
    assert len(share.pending(tmp_path)) == 1


def test_record_matches_published_schema(record):
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((ROOT / "deploy/results-repo/schema/record.v1.json").read_text())
    jsonschema.validate(share.make_submission(record), schema)
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(record | {"submission": {"verified": False}}, schema)  # host still present


def test_cli_dry_run_writes_nothing(record, tmp_path, monkeypatch):
    from scflow.cli import main
    monkeypatch.delenv("SCFLOW_SHARE")
    con_db = tmp_path / "t.db"
    assert main(["share", "--all", "--dry-run", "--db", str(con_db), "--runs-dir", str(tmp_path / "runs")]) == 0
    assert not (tmp_path / "runs" / "outbox").exists()


def test_cli_no_share_flag_keeps_outbox_empty(tmp_path, monkeypatch):
    from scflow.cli import main
    monkeypatch.delenv("SCFLOW_SHARE")
    rc = main(["validate", str(ROOT / "compounds/Nb_0GPa.yaml"), "--run-dir", str(FIX / "Nb"), "--record",
               "--db", str(tmp_path / "runs/scflow.db"), "--no-share"])
    assert rc == 0 and share.pending(tmp_path / "runs") == []


def test_absolute_paths_in_check_messages_are_scrubbed(record):
    record["checks"].append({"name": "x", "status": "skip", "message": "/home/alice/runs/Nb/2_phonons/ph.out not found"})
    msg = share.make_submission(record)["checks"][-1]["message"]
    assert msg == "ph.out not found"
