"""Load and validate compound YAML configs, resolving machine-specific paths."""

import copy
import hashlib
import json
import os
from pathlib import Path

import yaml

REQUIRED = {
    "compound": ["prefix", "pressure_GPa"],
    "structure": ["ibrav", "nat", "ntyp", "species", "positions"],
    "wannier": ["nbndsub", "projections", "dis_win_min", "dis_win_max", "dis_froz_min", "dis_froz_max"],
    "eliashberg": ["muc"],
    "dft": ["ecutwfc", "ecutrho", "degauss", "scf_kgrid", "coarse_kgrid"],
    "epw": ["fine_kgrid", "fsthick", "wscut"],
}


class ConfigError(ValueError):
    pass


def load_config(path, qe_bin=None, pseudo_dir=None, nproc=None):
    """Load a compound YAML. Paths resolve as: argument > SCFLOW_* env var > YAML `paths:`."""
    with open(path) as f:
        cfg = yaml.safe_load(f)
    _check_required(cfg, path)

    paths = cfg.setdefault("paths", {})
    paths["qe_bin"] = qe_bin or os.environ.get("SCFLOW_QE_BIN") or paths.get("qe_bin")
    paths["pseudo_dir"] = pseudo_dir or os.environ.get("SCFLOW_PSEUDO_DIR") or paths.get("pseudo_dir")
    for key in ("qe_bin", "pseudo_dir"):
        if not paths[key]:
            raise ConfigError(
                f"{path}: no {key}. Set SCFLOW_{key.upper()}, pass --{key.replace('_', '-')}, "
                f"or add paths.{key} to the YAML."
            )
    paths["pseudo_dir"] = str(paths["pseudo_dir"]).rstrip("/") + "/"

    par = cfg.setdefault("parallel", {})
    env_nproc = os.environ.get("SCFLOW_NPROC")
    par["nproc"] = int(nproc or env_nproc or par.get("nproc") or os.cpu_count() or 1)

    s = cfg["structure"]
    if s["ibrav"] in (2, 3) and s.get("A") is None:
        raise ConfigError(f"{path}: structure.A is required for ibrav={s['ibrav']}")
    if s["ibrav"] == 0 and "cell_parameters_angstrom" not in s:
        raise ConfigError(f"{path}: structure.cell_parameters_angstrom is required for ibrav=0")
    if len(s["positions"]) != s["nat"]:
        raise ConfigError(f"{path}: nat={s['nat']} but {len(s['positions'])} positions given")
    if len(s["species"]) != s["ntyp"]:
        raise ConfigError(f"{path}: ntyp={s['ntyp']} but {len(s['species'])} species given")

    nq = cfg["dft"]["coarse_kgrid"]
    nk = cfg["dft"].get("epw_coarse_kgrid", nq)
    if any(k % q for k, q in zip(nk, nq)):
        raise ConfigError(f"{path}: epw_coarse_kgrid {nk} must be a multiple of coarse_kgrid {nq}")
    return cfg


def _check_required(cfg, path):
    missing = []
    for section, keys in REQUIRED.items():
        if section not in cfg:
            missing.append(section)
            continue
        missing += [f"{section}.{k}" for k in keys if k not in cfg[section]]
    if missing:
        raise ConfigError(f"{path}: missing required keys: {', '.join(missing)}")


def physics_fingerprint(cfg):
    """Hash of everything that determines the physics, excluding machine paths and parallelism."""
    c = copy.deepcopy(cfg)
    c.pop("paths", None)
    c.pop("parallel", None)
    c.pop("reference", None)
    blob = json.dumps(c, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def structure_fingerprint(cfg):
    """Hash of the crystal structure + pressure only, used to group runs of the same material."""
    s = cfg["structure"]
    key = {
        "pressure_GPa": float(cfg["compound"]["pressure_GPa"]),
        "ibrav": s["ibrav"],
        "A": s.get("A"),
        "cell": s.get("cell_parameters_angstrom"),
        "positions": [[p[0]] + [round(float(x), 6) for x in p[1:]] for p in s["positions"]],
    }
    return hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:16]


def formula(cfg):
    counts = {}
    for p in cfg["structure"]["positions"]:
        counts[p[0]] = counts.get(p[0], 0) + 1
    return "".join(f"{el}{n if n > 1 else ''}" for el, n in counts.items())


def default_outdir(cfg):
    return Path(os.environ.get("SCFLOW_RUNS_DIR", "runs")) / cfg["compound"]["prefix"]
