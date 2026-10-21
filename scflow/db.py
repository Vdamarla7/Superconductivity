"""Results database: SQLite locally, with a JSON record format for submitting to a shared database.

Confidence tiers (computed, never hand-set):
  invalid      — at least one validation check failed
  preliminary  — checks pass, but Tc convergence w.r.t. the fine grid has not been shown
  converged    — ≥2 runs of the same physics on different fine grids agree on λ within 10%
  benchmark    — converged, and Tc is within 25% of a reference value at the same pressure (±10%)
"""

import hashlib
import json
import sqlite3
import uuid
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .config import formula, physics_fingerprint, structure_fingerprint

SCHEMA_VERSION = 1
CONVERGENCE_TOL = 0.10

SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    id              TEXT PRIMARY KEY,
    created_at      TEXT NOT NULL,
    prefix          TEXT NOT NULL,
    formula         TEXT NOT NULL,
    pressure_gpa    REAL NOT NULL,
    spacegroup      TEXT,
    structure_hash  TEXT NOT NULL,
    physics_hash    TEXT NOT NULL,
    method_hash     TEXT NOT NULL,
    phonon_method   TEXT NOT NULL,
    tc_path         TEXT NOT NULL,
    ecutwfc_ry      REAL,
    coarse_q        TEXT,
    coarse_k        TEXT,
    fine_k          TEXT,
    fine_k_points   INTEGER,   -- fine k-points × fine q-points
    muc             REAL,
    scflow_version  TEXT,
    qe_version      TEXT,
    epw_version     TEXT,
    host            TEXT,
    backend         TEXT,
    run_dir         TEXT,
    worst_check     TEXT NOT NULL,
    config_json     TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS results (
    run_id                  TEXT PRIMARY KEY REFERENCES runs(id) ON DELETE CASCADE,
    lambda                  REAL,
    omega_log_mev           REAL,
    omega2_mev              REAL,
    omega_max_mev           REAL,
    tc_allen_dynes_k        REAL,
    tc_allen_dynes_f1f2_k   REAL,
    tc_eliashberg_iso_k     REAL,
    tc_epw_k                REAL,
    scf_pressure_gpa        REAL,
    phonon_min_cm1          REAL,
    tc_reference_k          REAL
);
CREATE TABLE IF NOT EXISTS checks (
    run_id   TEXT NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    name     TEXT NOT NULL,
    status   TEXT NOT NULL,
    message  TEXT,
    value    TEXT,
    PRIMARY KEY (run_id, name)
);
CREATE TABLE IF NOT EXISTS artifacts (
    run_id  TEXT NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    path    TEXT NOT NULL,
    sha256  TEXT NOT NULL,
    bytes   INTEGER NOT NULL,
    PRIMARY KEY (run_id, path)
);
CREATE INDEX IF NOT EXISTS runs_method ON runs(method_hash);
"""

ARTIFACT_GLOBS = ["*/*.in", "*/scf.out", "*/ph.out", "*/nscf.out", "*/epw_pass*.out", "*/*.wout",
                  "*/*.a2f*", "manifest.json", "relaxed_structure.yaml"]


def connect(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA foreign_keys = ON")
    con.executescript(SCHEMA)
    return con


def method_fingerprint(cfg):
    """Everything that defines the physics except the convergence parameter being tested (fine grid)."""
    c = json.loads(json.dumps(cfg, default=str))
    c.get("epw", {}).pop("fine_kgrid", None)
    c.get("epw", {}).pop("fine_qgrid", None)
    return physics_fingerprint(c)


def fine_sampling(cfg):
    """(label, number of fine k-points × fine q-points) — the convergence parameter."""
    k = cfg["epw"]["fine_kgrid"]
    q = cfg["epw"].get("fine_qgrid", k)
    label = json.dumps(k) if q == k else f"k{json.dumps(k)} q{json.dumps(q)}"
    return label, k[0] * k[1] * k[2] * q[0] * q[1] * q[2]


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def record_run(con, cfg, run_dir, report, backend="local", host=None):
    run_dir = Path(run_dir).resolve()
    m = report.metrics
    fine_label, fine_n = fine_sampling(cfg)
    run_id = str(uuid.uuid4())
    con.execute(
        "INSERT INTO runs VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (run_id, datetime.now(timezone.utc).isoformat(), cfg["compound"]["prefix"], formula(cfg),
         float(cfg["compound"]["pressure_GPa"]), m.get("spacegroup"), structure_fingerprint(cfg),
         physics_fingerprint(cfg), method_fingerprint(cfg),
         cfg.get("phonons", {}).get("method", "harmonic"),
         "a2f" if cfg["epw"].get("use_a2f") else "eliashberg",
         cfg["dft"]["ecutwfc"], json.dumps(cfg["dft"]["coarse_kgrid"]),
         json.dumps(cfg["dft"].get("epw_coarse_kgrid", cfg["dft"]["coarse_kgrid"])),
         fine_label, fine_n, cfg["eliashberg"]["muc"],
         __version__, m.get("qe_version"), m.get("epw_version"), host, backend, str(run_dir),
         report.worst, json.dumps({k: v for k, v in cfg.items() if k != "paths"}, default=str)))
    con.execute(
        "INSERT INTO results VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
        (run_id, m.get("lambda"), m.get("omega_log_meV"), m.get("omega2_meV"), m.get("omega_max_meV"),
         m.get("tc_allen_dynes_K"), m.get("tc_allen_dynes_f1f2_K"), m.get("tc_eliashberg_iso_K"),
         m.get("tc_allen_dynes_epw_K"), m.get("scf_pressure_GPa"), m.get("phonon_min_cm1"),
         m.get("tc_reference_K")))
    con.executemany("INSERT INTO checks VALUES (?,?,?,?,?)",
                    [(run_id, c.name, c.status, c.message, json.dumps(c.value, default=str)) for c in report.checks])
    files = sorted({p for g in ARTIFACT_GLOBS for p in run_dir.glob(g) if p.is_file()})
    con.executemany("INSERT INTO artifacts VALUES (?,?,?,?)",
                    [(run_id, str(p.relative_to(run_dir)), _sha256(p), p.stat().st_size) for p in files])
    con.commit()
    return run_id


def convergence(con, method_hash):
    """Fine-grid convergence for one method: (status, message, [(fine_k_points, λ), ...])."""
    rows = con.execute(
        """SELECT r.fine_k_points AS n, res.lambda AS lam FROM runs r JOIN results res ON res.run_id = r.id
           WHERE r.method_hash = ? AND r.worst_check != 'fail' AND res.lambda IS NOT NULL
           ORDER BY r.fine_k_points""", (method_hash,)).fetchall()
    by_grid = {}
    for r in rows:
        by_grid[r["n"]] = r["lam"]
    series = sorted(by_grid.items())
    if len(series) < 2:
        return "not_assessed", "needs ≥2 runs on different fine grids", series
    (_, l1), (_, l2) = series[-2], series[-1]
    rel = abs(l2 - l1) / l2
    if rel <= CONVERGENCE_TOL:
        return "converged", f"λ changes {rel * 100:.1f}% between the two densest fine grids", series
    return "unconverged", f"λ changes {rel * 100:.1f}% between the two densest fine grids (>{CONVERGENCE_TOL * 100:.0f}%)", series


def confidence_tier(con, run_id):
    run = con.execute("SELECT * FROM runs WHERE id = ?", (run_id,)).fetchone()
    if run["worst_check"] == "fail":
        return "invalid", "failed validation checks"
    status, msg, series = convergence(con, run["method_hash"])
    densest = series[-1][0] if series else None
    if status != "converged" or run["fine_k_points"] != densest:
        return "preliminary", msg if status != "converged" else "a denser fine grid exists for this method"
    ref = con.execute("SELECT status FROM checks WHERE run_id = ? AND name = 'reference.tc'", (run_id,)).fetchone()
    if ref and ref["status"] == "pass":
        return "benchmark", msg + "; matches reference"
    return "converged", msg


def list_runs(con):
    rows = con.execute(
        """SELECT r.id, r.prefix, r.formula, r.pressure_gpa, r.fine_k, r.muc, r.worst_check, r.created_at,
                  res.lambda, res.tc_eliashberg_iso_k, res.tc_allen_dynes_k, res.tc_reference_k
           FROM runs r LEFT JOIN results res ON res.run_id = r.id ORDER BY r.created_at""").fetchall()
    return [dict(r) | {"tier": confidence_tier(con, r["id"])[0]} for r in rows]


def export_record(con, run_id):
    """Self-contained JSON record for one run — the unit submitted to a shared database."""
    run = dict(con.execute("SELECT * FROM runs WHERE id = ?", (run_id,)).fetchone())
    res = con.execute("SELECT * FROM results WHERE run_id = ?", (run_id,)).fetchone()
    checks = [dict(r) for r in con.execute("SELECT name, status, message, value FROM checks WHERE run_id = ?", (run_id,))]
    arts = [dict(r) for r in con.execute("SELECT path, sha256, bytes FROM artifacts WHERE run_id = ?", (run_id,))]
    tier, why = confidence_tier(con, run_id)
    run["config"] = json.loads(run.pop("config_json"))
    run.pop("run_dir")
    return {
        "schema_version": SCHEMA_VERSION,
        "run": run,
        "results": {k: v for k, v in dict(res).items() if k != "run_id"} if res else {},
        "checks": checks,
        "artifacts": arts,
        "confidence": {"tier": tier, "reason": why},
    }
