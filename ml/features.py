"""Features computed from a crystal structure alone (no DFT outputs), so models apply to new structures."""

import warnings
from functools import lru_cache

import numpy as np
import pandas as pd
from pymatgen.core import Composition, Structure

warnings.filterwarnings("ignore", message=".*fractional.*")
warnings.filterwarnings("ignore", category=UserWarning, module="matminer")

LIGHT = {"H", "Li", "Be", "B", "C", "N", "O"}


@lru_cache(maxsize=1)
def _featurizers():
    from matminer.featurizers.composition import ElementProperty, Stoichiometry, ValenceOrbital

    return ElementProperty.from_preset("magpie", impute_nan=True), Stoichiometry(), ValenceOrbital(impute_nan=True)


def composition_features(comp: Composition):
    ep, st, vo = _featurizers()
    vals = ep.featurize(comp) + st.featurize(comp) + vo.featurize(comp)
    names = ep.feature_labels() + st.feature_labels() + vo.feature_labels()
    out = dict(zip(names, vals))
    frac = comp.fractional_composition
    out["frac_H"] = frac["H"] if "H" in frac else 0.0
    out["frac_light"] = sum(frac[e] for e in frac.elements if e.symbol in LIGHT)
    out["n_elements"] = len(comp.elements)
    return out


def structure_features(s: Structure):
    comp = s.composition
    natoms = comp.num_atoms
    vpa = s.volume / natoms
    a, b, c = s.lattice.abc
    radii = []
    for site in s:
        for el, occ in site.species.items():
            r = getattr(el, "atomic_radius", None) or getattr(el, "atomic_radius_calculated", None)
            if r:
                radii.append((float(r), occ))
    packing = sum(4 / 3 * np.pi * r ** 3 * occ for r, occ in radii) / s.volume if radii else np.nan
    try:
        sg = s.get_space_group_info(symprec=0.1)[1]
    except Exception:
        sg = 0
    return {
        "volume_per_atom": vpa,
        "density": float(s.density),
        "packing_fraction": packing,
        "spacegroup": sg,
        "crystal_system": _crystal_system(sg),
        "abc_max_over_min": max(a, b, c) / min(a, b, c),
        "nsites": len(s),
        "is_ordered": float(s.is_ordered),
    }


def _crystal_system(sg):
    for upper, idx in ((2, 0), (15, 1), (74, 2), (142, 3), (167, 4), (194, 5), (230, 6)):
        if sg <= upper:
            return idx if sg > 0 else -1
    return -1


def featurize(structure: Structure):
    return composition_features(structure.composition) | structure_features(structure)


def featurize_many(structures, n_jobs=8):
    from concurrent.futures import ProcessPoolExecutor

    with ProcessPoolExecutor(n_jobs) as ex:
        rows = list(ex.map(_safe_featurize, structures, chunksize=64))
    return pd.DataFrame(rows)


def _safe_featurize(s):
    try:
        return featurize(s)
    except Exception:
        return {}


def chemsys(comp: Composition):
    return "-".join(sorted(e.symbol for e in comp.elements))


def structure_from_scflow_yaml(path):
    """Build a pymatgen Structure from an scflow compound YAML (lattice from ibrav/A or explicit cell)."""
    import yaml

    from scflow.validation import lattice_vectors

    cfg = yaml.safe_load(open(path))
    s = cfg["structure"]
    return Structure(lattice_vectors(cfg), [p[0] for p in s["positions"]],
                     [[float(x) for x in p[1:4]] for p in s["positions"]]), cfg
