"""Every shipped compound config loads, passes pre-flight, and generates QE/EPW inputs."""

from pathlib import Path

import pytest

from scflow.cli import main
from scflow.config import ConfigError, load_config
from scflow.generators import gen_epw_pass1, gen_epw_pass2, gen_nscf, gen_phonon, gen_scf, gen_vcrelax
from scflow.validation import Report, check_structure

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = sorted(p for p in (ROOT / "compounds").glob("*.yaml") if p.name != "template.yaml")


@pytest.fixture(autouse=True)
def fake_paths(monkeypatch, tmp_path):
    monkeypatch.setenv("SCFLOW_QE_BIN", str(tmp_path / "bin"))
    monkeypatch.setenv("SCFLOW_PSEUDO_DIR", str(tmp_path / "pseudo"))


@pytest.mark.parametrize("path", CONFIGS, ids=lambda p: p.stem)
def test_config_generates_inputs(path):
    cfg = load_config(path)
    for gen in (gen_vcrelax, gen_scf, gen_phonon, gen_nscf, gen_epw_pass1, gen_epw_pass2):
        text = gen(cfg)
        assert cfg["compound"]["prefix"] in text or gen is gen_phonon
    assert "/home/" not in gen_scf(cfg).replace(cfg["paths"]["pseudo_dir"], "")


@pytest.mark.parametrize("path", CONFIGS, ids=lambda p: p.stem)
def test_config_passes_preflight(path):
    report = Report()
    check_structure(load_config(path), report)
    assert report.worst != "fail", [c.message for c in report.checks if c.status == "fail"]


def test_nscf_has_full_unreduced_grid():
    cfg = load_config(ROOT / "compounds/Nb_0GPa.yaml")
    n = cfg["dft"]["coarse_kgrid"]
    text = gen_nscf(cfg)
    assert f"K_POINTS {{crystal}}\n{n[0] * n[1] * n[2]}" in text


def test_missing_paths_is_a_clear_error(monkeypatch):
    monkeypatch.delenv("SCFLOW_QE_BIN")
    with pytest.raises(ConfigError, match="SCFLOW_QE_BIN"):
        load_config(ROOT / "compounds/Nb_0GPa.yaml")


def test_dry_run_writes_all_inputs(tmp_path):
    out = tmp_path / "nb"
    assert main(["run", str(ROOT / "compounds/Nb_0GPa.yaml"), "--outdir", str(out), "--dry-run"]) == 0
    assert (out / "2_phonons/scf.in").exists()


def test_slurm_script(tmp_path):
    assert main(["submit", str(ROOT / "compounds/Nb_0GPa.yaml"), "--backend", "slurm", "--runs-dir",
                 str(tmp_path), "--nproc", "32", "--dry-run"]) == 0
    script = (tmp_path / "Nb.sbatch").read_text()
    assert "#SBATCH --ntasks=32" in script and "apptainer run" in script


def test_nscf_prints_eigenvalues_for_window_scan():
    cfg = load_config(ROOT / "compounds/CaH6_200GPa.yaml")
    assert "verbosity     = 'high'" in gen_nscf(cfg)


def test_changed_input_is_not_silently_cached(tmp_path):
    from scflow.executors import StepFailed
    from scflow.workflow import Workflow

    cfg = load_config(ROOT / "compounds/Nb_0GPa.yaml")
    wf = Workflow(cfg, tmp_path / "nb", executor=None, dry_run=True, log=lambda *a: None)
    wf.scf()
    (tmp_path / "nb/2_phonons/scf.out").write_text("JOB DONE\n")
    wf.cfg["dft"]["ecutwfc"] = 90.0
    with pytest.raises(StepFailed, match="new --outdir"):
        wf.scf()
