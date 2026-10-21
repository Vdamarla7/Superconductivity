"""Physics-based stability screen with a universal machine-learned interatomic potential (MACE-MP-0).

  python ml/mlip_stability.py compounds/CaH6_200GPa.yaml
  python ml/mlip_stability.py structure.cif --pressure 0 --no-hull

For one structure at a given pressure:
  1. relax cell + positions with MACE-MP-0 at that pressure
  2. dynamic stability: finite-displacement phonons on a supercell, frequencies on a q-mesh
  3. thermodynamic stability: enthalpy above the convex hull built from every Materials Project
     polymorph in the chemical system, all re-relaxed with the same potential at the same
     pressure (a self-consistent hull, valid at high pressure where MP's 0 GPa hull is not)

MACE-MP-0 is trained on near-ambient DFT data (MPtrj). At megabar pressures it is
extrapolating; scflow's DFPT is the arbiter.
"""

import argparse
import json
import math
import sys
import tempfile
from itertools import combinations
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE
from ase.phonons import Phonons
from ase.units import GPa

sys.path.insert(0, str(Path(__file__).parent))

CM1_PER_EV = 8065.544


def calculator(device=None):
    import torch
    from mace.calculators import mace_mp

    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    return mace_mp(model="medium", device=device, default_dtype="float64")


def to_atoms(structure):
    return Atoms([str(s.specie) for s in structure], cell=structure.lattice.matrix,
                 scaled_positions=structure.frac_coords, pbc=True)


def relax(atoms, calc, pressure_gpa=0.0, fmax=0.005, steps=1000):
    atoms = atoms.copy()
    atoms.calc = calc
    opt = FIRE(FrechetCellFilter(atoms, scalar_pressure=pressure_gpa * GPa), logfile=None)
    converged = opt.run(fmax=fmax, steps=steps)
    return atoms, bool(converged)


def enthalpy_per_atom(atoms, pressure_gpa):
    return (atoms.get_potential_energy() + pressure_gpa * GPa * atoms.get_volume()) / len(atoms)


def supercell_for(atoms, min_length=9.0):
    lengths = atoms.cell.lengths()
    return tuple(max(2, math.ceil(min_length / L)) for L in lengths)


def _signed(w):
    """ASE returns imaginary frequencies as complex numbers; report them as negative reals."""
    w = np.asarray(w)
    if np.iscomplexobj(w):
        return np.where(np.abs(w.imag) > 1e-8, -np.abs(w.imag), w.real)
    return w


def force_constants(atoms, calc, supercell=None, delta=0.01):
    supercell = supercell or supercell_for(atoms)
    tmp = tempfile.mkdtemp(prefix="scflow-ph-")
    ph = Phonons(atoms, calc, supercell=supercell, delta=delta, name=str(Path(tmp) / "ph"))
    ph.run()
    ph.read(acoustic=True)
    ph.clean()
    return ph


def frequencies_at(ph, qpoints_frac):
    """Frequencies (cm⁻¹, imaginary negative) at q-points given in reciprocal-lattice (crystal) units."""
    return _signed(ph.band_structure(np.asarray(qpoints_frac, dtype=float))) * CM1_PER_EV


def phonon_frequencies(atoms, calc, supercell=None, mesh=6):
    ph = force_constants(atoms, calc, supercell)
    q = np.array([[i / mesh, j / mesh, k / mesh] for i in range(mesh) for j in range(mesh) for k in range(mesh)])
    bands = frequencies_at(ph, q)
    return {"supercell": list(ph.supercell if hasattr(ph, "supercell") else supercell_for(atoms)), "q_mesh": mesh,
            "gamma_cm1": sorted(float(x) for x in bands[0]),
            "min_cm1_off_gamma": float(bands[1:].min()), "max_cm1": float(bands.max())}


def hull_distance(atoms, calc, pressure_gpa, max_sites=48, ehull_cutoff=0.2):
    """Enthalpy above a MACE hull built from MP polymorphs re-relaxed at the same pressure (eV/atom)."""
    from mp_api.client import MPRester
    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
    from pymatgen.core import Composition

    from fetch_data import api_key

    elements = sorted(set(atoms.get_chemical_symbols()))
    systems = ["-".join(c) for i in range(1, len(elements) + 1) for c in combinations(elements, i)]
    docs = []
    with MPRester(api_key()) as mpr:
        for cs in systems:
            docs += mpr.materials.summary.search(chemsys=cs, deprecated=False,
                                                 fields=["material_id", "structure", "energy_above_hull"])
    entries, used = [], []
    for d in docs:
        if d.energy_above_hull is None or d.energy_above_hull > ehull_cutoff or len(d.structure) > max_sites:
            continue
        a, ok = relax(to_atoms(d.structure), calc, pressure_gpa, fmax=0.01, steps=600)
        h = enthalpy_per_atom(a, pressure_gpa)
        entries.append(PDEntry(Composition(a.get_chemical_formula()), h * len(a), name=str(d.material_id)))
        used.append(str(d.material_id))
    target = PDEntry(Composition(atoms.get_chemical_formula()), enthalpy_per_atom(atoms, pressure_gpa) * len(atoms),
                     name="candidate")
    pd_ = PhaseDiagram(entries + [target])
    e_hull = pd_.get_e_above_hull(target)
    decomp = pd_.get_decomposition(target.composition)
    return {"e_above_hull_eV_per_atom": float(e_hull), "n_reference_phases": len(entries),
            "decomposes_to": sorted(e.name for e in decomp) if e_hull > 1e-6 else []}


def screen(path, pressure=None, do_hull=True, device=None):
    from predict import load_structure

    structure, p_cfg = load_structure(path)
    pressure = p_cfg if pressure is None else pressure
    calc = calculator(device)
    atoms, converged = relax(to_atoms(structure), calc, pressure)
    a0 = to_atoms(structure)
    result = {
        "input": str(path),
        "formula": atoms.get_chemical_formula(mode="reduce"),
        "pressure_GPa": pressure,
        "relax_converged": converged,
        "volume_per_atom_input": a0.get_volume() / len(a0),
        "volume_per_atom_relaxed": atoms.get_volume() / len(atoms),
        "cell_lengths_relaxed": [float(x) for x in atoms.cell.lengths()],
    }
    ph = phonon_frequencies(atoms, calc)
    result["phonons"] = ph
    result["dynamically_stable"] = ph["min_cm1_off_gamma"] > -20.0
    if do_hull:
        result["hull"] = hull_distance(atoms, calc, pressure)
    if pressure > 20:
        result["warning"] = ("MACE-MP-0 is trained on near-ambient DFT; at this pressure it extrapolates. "
                             "Treat as a screen and verify with scflow DFPT.")
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structures", nargs="+")
    ap.add_argument("--pressure", type=float, help="GPa (default: from scflow YAML, else 0)")
    ap.add_argument("--no-hull", action="store_true")
    ap.add_argument("--device")
    args = ap.parse_args()
    for p in args.structures:
        print(json.dumps(screen(p, args.pressure, not args.no_hull, args.device), indent=2))


if __name__ == "__main__":
    main()
