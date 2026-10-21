"""Compare MACE-MP-0 phonons with scflow DFPT phonons at the same geometry and q-points.

  python ml/validate_mlip.py CONFIG.yaml PH_OUT [CONFIG.yaml PH_OUT ...]

Answers: does the universal potential reach the same dynamic-stability verdict as DFT,
and how large are its frequency errors? (Phonons are computed at the config geometry,
which is what the DFPT run used, so only the force constants are being compared.)
"""

import json
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from mlip_stability import calculator, force_constants, frequencies_at, to_atoms  # noqa: E402
from predict import load_structure  # noqa: E402

from scflow.postprocess import parse_phonon_qpoints  # noqa: E402
from scflow.validation import lattice_vectors  # noqa: E402


def qe_alat(cfg):
    s = cfg["structure"]
    return float(s["A"]) if s["ibrav"] != 0 else float(np.linalg.norm(lattice_vectors(cfg)[0]))


def compare(config, ph_out, calc):
    cfg = yaml.safe_load(open(config))
    structure, _ = load_structure(config)
    ph = force_constants(to_atoms(structure), calc)
    A = lattice_vectors(cfg) / qe_alat(cfg)
    dft, mlip = [], []
    for q_cart, freqs in parse_phonon_qpoints(ph_out):
        q_frac = A @ np.array(q_cart)
        dft.append(sorted(freqs))
        mlip.append(sorted(frequencies_at(ph, [q_frac])[0]))
    dft, mlip = np.array(dft), np.array(mlip)
    optical = np.abs(dft) > 50
    return {
        "config": str(config),
        "n_qpoints": len(dft),
        "dft_min_cm1": float(dft.min()), "mlip_min_cm1": float(mlip.min()),
        "dft_max_cm1": float(dft.max()), "mlip_max_cm1": float(mlip.max()),
        "mae_cm1": float(np.abs(dft - mlip).mean()),
        "rel_err_optical": float(np.median(np.abs(dft - mlip)[optical] / np.abs(dft)[optical])) if optical.any() else None,
        "stability_agrees": bool((dft.min() > -20) == (mlip.min() > -20)),
    }


def main():
    args = sys.argv[1:]
    if not args or len(args) % 2:
        raise SystemExit(__doc__)
    calc = calculator()
    rows = [compare(args[i], args[i + 1], calc) for i in range(0, len(args), 2)]
    print(json.dumps(rows, indent=2))
    print(f"\n  {'config':<40}{'DFT min':>9}{'MACE min':>10}{'DFT max':>9}{'MACE max':>10}{'MAE':>7}  verdict")
    for r in rows:
        print(f"  {Path(r['config']).name:<40}{r['dft_min_cm1']:>9.0f}{r['mlip_min_cm1']:>10.0f}{r['dft_max_cm1']:>9.0f}"
              f"{r['mlip_max_cm1']:>10.0f}{r['mae_cm1']:>7.0f}  {'agrees' if r['stability_agrees'] else 'DISAGREES'}")


if __name__ == "__main__":
    main()
