"""Share finished runs with the public scflow results database, a GitHub repository.

Sharing is on by default. Each finished run's record (structure, settings, results, check
outcomes, code versions, file hashes — no hostnames or paths) is queued in runs/outbox/ and
sent as a pull request when GitHub credentials are available (GITHUB_TOKEN / GH_TOKEN or a
logged-in `gh` CLI). Without credentials it stays queued until `scflow share` is run.

Disable with --no-share, SCFLOW_SHARE=0, or `share: false` in ~/.config/scflow/config.yaml.
"""

import base64
import json
import os
import re
import shutil
import subprocess
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import yaml

from . import __version__

ABS_PATH = re.compile(r"(?:/[\w.-]+){2,}")
RESULTS_REPO = os.environ.get("SCFLOW_RESULTS_REPO", "Vdamarla7/scflow-results")
API = "https://api.github.com"


class ShareError(RuntimeError):
    pass


def config_dir():
    return Path(os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config") / "scflow"


def share_enabled(no_share=False):
    if no_share:
        return False
    env = os.environ.get("SCFLOW_SHARE")
    if env is not None:
        return env.strip().lower() not in ("0", "false", "no", "off")
    cfg = config_dir() / "config.yaml"
    if cfg.exists():
        try:
            return bool((yaml.safe_load(cfg.read_text()) or {}).get("share", True))
        except yaml.YAMLError:
            return True
    return True


def print_notice():
    print(f"  Results will be shared publicly at github.com/{RESULTS_REPO} (disable: --no-share)")
    marker = config_dir() / "share-notice-shown"
    if marker.exists():
        return
    print(f"""
  ── Data sharing ──────────────────────────────────────────────────────────────
  scflow contributes finished runs to a public results database:
    github.com/{RESULTS_REPO}
  Sent: the crystal structure and all calculation settings, the results (λ, ω_log, Tc),
  every check outcome, scflow/QE/EPW versions, and SHA-256 hashes of input/output files.
  Not sent: hostname, file paths, raw output files, credentials, anything else.
  Records are public and labeled "self-reported, unverified".
  Turn it off: --no-share on a command, SCFLOW_SHARE=0, or `share: false` in
  {config_dir() / 'config.yaml'}
  ──────────────────────────────────────────────────────────────────────────────
""")
    try:
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(datetime.now(timezone.utc).isoformat() + "\n")
    except OSError:
        pass


def make_submission(record):
    """The record as it will be published: no hostname or run directory."""
    sub = json.loads(json.dumps(record, default=str))
    sub["run"].pop("host", None)
    sub["run"].pop("run_dir", None)
    sub["artifacts"] = [{"path": a["path"], "sha256": a["sha256"], "bytes": a["bytes"]}
                        for a in sub.get("artifacts", []) if not Path(a["path"]).is_absolute()]
    for c in sub.get("checks", []):
        if c.get("message"):
            c["message"] = ABS_PATH.sub(lambda m: Path(m.group(0)).name, c["message"])
    sub["submission"] = {
        "scflow_version": __version__,
        "submitted_at": datetime.now(timezone.utc).isoformat(),
        "verified": False,
    }
    return sub


def repo_path(sub):
    run = sub["run"]
    return f"records/{run['formula']}/{float(run['pressure_gpa']):g}GPa/{run['id']}.json"


def outbox_dir(runs_dir):
    return Path(runs_dir) / "outbox"


def write_outbox(sub, runs_dir):
    path = outbox_dir(runs_dir) / f"{sub['run']['id']}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(sub, indent=2) + "\n")
    return path


def pending(runs_dir):
    d = outbox_dir(runs_dir)
    return sorted(d.glob("*.json")) if d.exists() else []


def github_token():
    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if token:
        return token
    if shutil.which("gh"):
        try:
            out = subprocess.run(["gh", "auth", "token"], capture_output=True, text=True, timeout=10)
            if out.returncode == 0 and out.stdout.strip():
                return out.stdout.strip()
        except (OSError, subprocess.TimeoutExpired):
            pass
    return None


def _api(method, path, token, body=None):
    req = urllib.request.Request(API + path, method=method,
                                 data=json.dumps(body).encode() if body is not None else None,
                                 headers={"Authorization": f"Bearer {token}",
                                          "Accept": "application/vnd.github+json",
                                          "User-Agent": f"scflow/{__version__}"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.loads(r.read() or b"{}")
    except urllib.error.HTTPError as e:
        raise ShareError(f"GitHub API {method} {path}: HTTP {e.code} {e.read()[:200]!r}") from e
    except (urllib.error.URLError, TimeoutError) as e:
        raise ShareError(f"GitHub API unreachable: {e}") from e


def send(paths, token=None, repo=RESULTS_REPO, api=_api):
    """Open one pull request adding these outbox files to the results repo. Returns the PR URL."""
    paths = [Path(p) for p in paths]
    if not paths:
        raise ShareError("nothing to send")
    token = token or github_token()
    if not token:
        raise ShareError("no GitHub credentials (set GITHUB_TOKEN or run `gh auth login`)")
    owner = repo.split("/")[0]
    login = api("GET", "/user", token)["login"]
    head_repo = repo if login.lower() == owner.lower() else api("POST", f"/repos/{repo}/forks", token)["full_name"]
    base_sha = api("GET", f"/repos/{repo}/git/ref/heads/main", token)["object"]["sha"]
    branch = f"scflow-{datetime.now(timezone.utc).strftime('%Y%m%d-%H%M%S')}-{os.getpid()}"
    for attempt in range(6):  # a new fork can take a few seconds to become writable
        try:
            api("POST", f"/repos/{head_repo}/git/refs", token, {"ref": f"refs/heads/{branch}", "sha": base_sha})
            break
        except ShareError:
            if attempt == 5:
                raise
            time.sleep(5)
    subs = [json.loads(p.read_text()) for p in paths]
    for sub in subs:
        api("PUT", f"/repos/{head_repo}/contents/{repo_path(sub)}", token,
            {"message": f"Add {sub['run']['formula']} {sub['run']['pressure_gpa']:g} GPa run {sub['run']['id'][:8]}",
             "content": base64.b64encode((json.dumps(sub, indent=2) + "\n").encode()).decode(),
             "branch": branch})
    lines = [f"- {s['run']['formula']} @ {s['run']['pressure_gpa']:g} GPa: "
             f"λ={s['results'].get('lambda')}, Tc(Eliashberg)={s['results'].get('tc_eliashberg_iso_k')} K, "
             f"self-reported tier `{s['confidence']['tier']}`" for s in subs]
    pr = api("POST", f"/repos/{repo}/pulls", token, {
        "title": f"scflow results: {', '.join(sorted({s['run']['formula'] for s in subs}))}",
        "head": f"{login}:{branch}" if head_repo != repo else branch,
        "base": "main",
        "body": "Submitted automatically by scflow " + __version__ + ".\n\n" + "\n".join(lines)})
    sent = paths[0].parent / "sent"
    sent.mkdir(exist_ok=True)
    for p in paths:
        p.rename(sent / p.name)
    return pr["html_url"]


def share_record(record, runs_dir, enabled=True):
    """Queue a finished run and try to send everything pending. Never raises."""
    if not enabled:
        return
    try:
        write_outbox(make_submission(record), runs_dir)
    except OSError as e:
        print(f"  ⚠ could not queue the result for sharing: {e}")
        return
    try:
        url = send(pending(runs_dir))
        print(f"  Shared with the public results database: {url}")
    except ShareError as e:
        print(f"  Result queued in {outbox_dir(runs_dir)} ({e}). Send it later with `scflow share`.")
