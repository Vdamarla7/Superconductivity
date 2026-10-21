"""Physics-level tests: Eliashberg integrals, the validation gate, and the database tiers.

Fixtures are trimmed outputs of real QE 7.4.1 / EPW 5.9 runs.
"""

import copy
from pathlib import Path

import numpy as np
import pytest
import yaml

from scflow import db
from scflow.eliashberg import A2F, allen_dynes, eliashberg_tc, summarize
from scflow.postprocess import parse_relaxed_structure
from scflow.validation import Report, check_structure, min_distances, validate

ROOT = Path(__file__).resolve().parents[1]
FIX = Path(__file__).parent / "fixtures"


def load(name):
    cfg = yaml.safe_load((ROOT / "compounds" / name).read_text())
    cfg.setdefault("parallel", {})["nproc"] = 1
    return cfg


def einstein(lam=1.0, w_e=50.0, width=1.0):
    w = np.linspace(0.05, 2 * w_e, 4000)
    f = np.exp(-(w - w_e) ** 2 / (2 * width ** 2))
    f *= lam / (2 * np.trapezoid(f / w, w))
    return A2F(w, f)


# ── Eliashberg theory ────────────────────────────────────────────────────────

def test_einstein_moments():
    lam, wlog, w2 = einstein(lam=1.3).moments()
    assert lam == pytest.approx(1.3, rel=1e-6)
    assert wlog == pytest.approx(50.0, rel=1e-3)
    assert w2 == pytest.approx(50.0, rel=1e-3)


def test_allen_dynes_matches_closed_form():
    tc = allen_dynes(1.0, 50.0, 0.1)
    expected = 50.0 * 11.604518 / 1.2 * np.exp(-1.04 * 2.0 / (1.0 - 0.1 * 1.62))
    assert tc == pytest.approx(expected, rel=1e-12)
    assert allen_dynes(0.1, 50.0, 0.13) == 0.0


def test_eliashberg_close_to_allen_dynes_f1f2_without_coulomb():
    # Allen & Dynes fitted their formula to Eliashberg solutions to within a few % for μ*=0.
    spec = einstein(lam=1.0)
    tc_el = eliashberg_tc(spec, 0.0, wc_mev=2000.0)
    tc_ad = allen_dynes(1.0, 50.0, 0.0, 50.0)
    assert tc_el == pytest.approx(tc_ad, rel=0.05)


def test_eliashberg_converged_in_matsubara_cutoff():
    spec = einstein(lam=1.0)
    assert eliashberg_tc(spec, 0.0, 1000.0) == pytest.approx(eliashberg_tc(spec, 0.0, 2000.0), rel=0.01)


def test_eliashberg_monotonic_in_muc_and_lambda():
    spec = einstein(lam=1.5)
    assert eliashberg_tc(spec, 0.10, 500.0) > eliashberg_tc(spec, 0.15, 500.0)
    assert eliashberg_tc(einstein(lam=2.0), 0.1, 500.0) > eliashberg_tc(spec, 0.1, 500.0)


def test_nb_a2f_reproduces_epw_numbers():
    spec = A2F.from_epw(FIX / "Nb/3_epw/Nb.a2f.01.300.000")
    s = summarize(spec, 0.10, 300.0)
    assert s["lambda"] == pytest.approx(1.8365505, rel=2e-3)
    assert s["omega_log_meV"] == pytest.approx(0.0006630 * 13605.693, rel=2e-3)
    assert s["tc_allen_dynes_K"] == pytest.approx(14.166529, rel=5e-3)


# ── Structure / validation gate ──────────────────────────────────────────────

def test_cah6_config_has_correct_12d_site():
    cfg = load("CaH6_200GPa.yaml")
    d = min_distances(cfg)
    assert d[("Ca", "H")] > 1.8
    report = Report()
    check_structure(cfg, report)
    assert report.worst != "fail"
    assert report.metrics["spacegroup"] == "Im-3m"


def test_wrong_wyckoff_site_is_caught_before_any_compute():
    cfg = load("CaH6_200GPa.yaml")
    cfg["structure"]["positions"] = [["Ca", 0, 0, 0], ["H", 0, .25, .25], ["H", 0, .75, .75], ["H", .25, 0, .25],
                                     ["H", .75, 0, .75], ["H", .25, .25, 0], ["H", .75, .75, 0]]
    report = Report()
    check_structure(cfg, report)
    by = {c.name: c for c in report.checks}
    assert by["structure.distance.Ca-H"].status == "fail"
    # Same space group — the symmetry check alone cannot catch this.
    assert by["structure.spacegroup"].status == "pass"


def test_wrong_site_cah6_run_fails_pressure_and_stability():
    cfg = load("CaH6_200GPa.yaml")
    report = validate(cfg, FIX / "CaH6_wrong_site")
    by = {c.name: c for c in report.checks}
    assert by["scf.pressure"].status == "fail"
    assert report.metrics["scf_pressure_GPa"] == pytest.approx(426.8, abs=0.1)
    assert by["phonons.stability"].status == "fail"


def test_nb_run_validation():
    cfg = load("Nb_0GPa.yaml")
    report = validate(cfg, FIX / "Nb")
    by = {c.name: c.status for c in report.checks}
    assert by["scf.pressure"] == "pass"
    assert by["scf.symmetry"] == "pass"
    assert by["phonons.stability"] == "pass"
    assert by["eliashberg.cross_check"] == "pass"
    assert by["wannier.disentanglement"] == "warn"
    assert by["epw.fine_grid_phonons"] == "warn"
    assert by["reference.tc"] == "warn"  # λ=1.84 → Tc far above experiment: not converged
    assert report.worst == "warn"


def test_expected_spacegroup_mismatch_fails():
    cfg = load("Nb_0GPa.yaml")
    cfg["structure"]["expected_spacegroup"] = "Fm-3m"
    report = Report()
    check_structure(cfg, report)
    assert {c.name: c.status for c in report.checks}["structure.spacegroup"] == "fail"


# ── vc-relax geometry propagation ────────────────────────────────────────────

def test_relaxed_structure_ibrav3_alat_units():
    cell, pos = parse_relaxed_structure(FIX / "vcrelax/Nb_ibrav3.out")
    a = 2 * max(abs(x) for x in cell[0])
    assert 3.2 < a < 3.4
    assert pos[0][0] == "Nb"


def test_relaxed_structure_ibrav0_angstrom():
    cell, pos = parse_relaxed_structure(FIX / "vcrelax/CaH6_ibrav0.out")
    assert cell[0][0] == pytest.approx(-1.699096493)
    assert len(pos) == 7 and pos[1] == ["H", 0.5, 0.75, 0.25]


# ── Database tiers ───────────────────────────────────────────────────────────

def test_confidence_tiers(tmp_path):
    con = db.connect(tmp_path / "t.db")
    cfg = load("Nb_0GPa.yaml")
    report = validate(cfg, FIX / "Nb")
    rid = db.record_run(con, cfg, FIX / "Nb", report)
    assert db.confidence_tier(con, rid)[0] == "preliminary"

    denser = copy.deepcopy(cfg)
    denser["epw"]["fine_kgrid"] = [20, 20, 20]
    rep2 = copy.deepcopy(report)
    rep2.metrics["lambda"] *= 1.05
    rid2 = db.record_run(con, denser, FIX / "Nb", rep2)
    assert db.confidence_tier(con, rid2)[0] == "converged"
    assert db.confidence_tier(con, rid)[0] == "preliminary"  # superseded by the denser grid

    rec = db.export_record(con, rid2)
    assert rec["schema_version"] == db.SCHEMA_VERSION
    assert "paths" not in rec["run"]["config"]
    assert any(a["path"].endswith("scf.out") for a in rec["artifacts"])


def test_failed_run_is_invalid_and_excluded_from_convergence(tmp_path):
    con = db.connect(tmp_path / "t.db")
    cfg = load("Nb_0GPa.yaml")
    bad = validate(cfg, FIX / "Nb")
    bad.add("synthetic.check", "fail", "synthetic failure")
    rid = db.record_run(con, cfg, FIX / "Nb", bad)
    assert db.confidence_tier(con, rid)[0] == "invalid"


def test_excluded_band_count():
    from scflow.postprocess import excluded_band_count
    assert excluded_band_count("exclude_bands = 1-4") == 4
    assert excluded_band_count("exclude_bands = 1-34") == 34
    assert excluded_band_count("exclude_bands = 1-4,7,9-10") == 7
    assert excluded_band_count(None) == 0


def test_frozen_window_clearance():
    from scflow.postprocess import frozen_window_clearance
    eigs = np.array([[1.0, 2.0, 2.0, 5.0], [1.5, 2.49, 3.0, 6.0]])
    assert frozen_window_clearance(eigs, 2.5) == pytest.approx(0.01)
    assert frozen_window_clearance(eigs, 4.0) == pytest.approx(1.0)
