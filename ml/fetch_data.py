"""Download training data into ml/data/.

  mp.parquet    Materials Project: structure, energy above hull, formation energy (stability labels)
  3dsc_mp.csv   3DSC_MP (Sommer et al., Sci. Data 10, 816 (2023), CC BY 4.0): SuperCon Tc
                matched to Materials Project structures, including measured non-superconductors

Needs MP_API_KEY in the environment or in a .env file at the repo root.
"""

import argparse
import os
import subprocess
import urllib.request
from pathlib import Path

import pandas as pd

DATA = Path(__file__).parent / "data"
DSC_URL = "https://raw.githubusercontent.com/aimat-lab/3DSC/main/superconductors_3D/data/final/MP/3DSC_MP.csv"
MP_FIELDS = ["material_id", "formula_pretty", "structure", "energy_above_hull", "formation_energy_per_atom",
             "is_stable", "theoretical", "nsites", "symmetry", "volume", "density", "band_gap", "is_metal"]


def api_key():
    key = os.environ.get("MP_API_KEY")
    env = Path(__file__).resolve().parents[1] / ".env"
    if not key and env.exists():
        for line in env.read_text().splitlines():
            if line.startswith("MP_API_KEY="):
                key = line.split("=", 1)[1].strip()
    if not key:
        raise SystemExit("MP_API_KEY not set (environment or .env)")
    return key


def fetch_3dsc():
    out = DATA / "3dsc_mp.csv"
    if not out.exists():
        urllib.request.urlretrieve(DSC_URL, out)
    repo = DATA / "3dsc_repo"
    if not repo.exists():
        # The CSV references one CIF file per entry; fetch only that folder.
        subprocess.run(["git", "clone", "-q", "--depth", "1", "--filter=blob:none", "--sparse",
                        "https://github.com/aimat-lab/3DSC.git", str(repo)], check=True)
        subprocess.run(["git", "sparse-checkout", "set", "superconductors_3D/data/final/MP/cifs"], cwd=repo, check=True)
    print(f"3DSC: {out} ({out.stat().st_size / 1e6:.1f} MB) + CIFs in {repo}")


def fetch_mp(limit=None):
    from mp_api.client import MPRester

    out = DATA / "mp.parquet"
    with MPRester(api_key()) as mpr:
        docs = mpr.materials.summary.search(deprecated=False, fields=MP_FIELDS,
                                            num_chunks=None if limit is None else 1,
                                            chunk_size=1000 if limit is None else limit)
    rows = []
    for d in docs:
        rows.append({
            "material_id": str(d.material_id),
            "formula": d.formula_pretty,
            "structure": d.structure.to_json(),
            "e_above_hull": d.energy_above_hull,
            "e_form": d.formation_energy_per_atom,
            "is_stable": d.is_stable,
            "theoretical": d.theoretical,
            "nsites": d.nsites,
            "spacegroup": d.symmetry.number if d.symmetry else None,
            "volume": d.volume,
            "density": d.density,
            "band_gap": d.band_gap,
            "is_metal": d.is_metal,
        })
    pd.DataFrame(rows).to_parquet(out)
    print(f"MP: {len(rows)} materials → {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, help="fetch only this many MP entries (for testing)")
    ap.add_argument("--skip-mp", action="store_true")
    args = ap.parse_args()
    DATA.mkdir(parents=True, exist_ok=True)
    fetch_3dsc()
    if not args.skip_mp:
        fetch_mp(args.limit)
