"""CI check for submitted records: schema, location, duplicates, no machine-identifying data.

  python scripts/check_records.py            # every record
  python scripts/check_records.py FILE ...   # only these
"""

import json
import re
import sys
from pathlib import Path

import jsonschema

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = json.loads((ROOT / "schema/record.v1.json").read_text())
LOCAL_PATH = re.compile(r"(/home/|/Users/|[A-Za-z]:\\\\|/root/|/scratch/)")


def check(path):
    errors = []
    try:
        rec = json.loads(path.read_text())
    except json.JSONDecodeError as e:
        return [f"not valid JSON: {e}"]
    try:
        jsonschema.validate(rec, SCHEMA)
    except jsonschema.ValidationError as e:
        errors.append(f"schema: {e.message} at {'/'.join(map(str, e.absolute_path))}")
        return errors
    run = rec["run"]
    expected = ROOT / "records" / run["formula"] / f"{float(run['pressure_gpa']):g}GPa" / f"{run['id']}.json"
    if path.resolve() != expected.resolve():
        errors.append(f"should be at {expected.relative_to(ROOT)}")
    if LOCAL_PATH.search(path.read_text()):
        errors.append("contains a local file path")
    return errors


def main():
    files = [Path(f) for f in sys.argv[1:]] or sorted((ROOT / "records").rglob("*.json"))
    ids, failed = {}, False
    for f in sorted((ROOT / "records").rglob("*.json")):
        rid = f.stem
        if rid in ids:
            print(f"✗ duplicate run id {rid}: {ids[rid]} and {f}")
            failed = True
        ids[rid] = f
    summary = []
    for f in files:
        errs = check(f)
        if errs:
            failed = True
            for e in errs:
                print(f"✗ {f}: {e}")
        else:
            rec = json.loads(f.read_text())
            summary.append(f"| {rec['run']['formula']} | {rec['run']['pressure_gpa']:g} | "
                           f"{rec['results'].get('lambda')} | {rec['results'].get('tc_eliashberg_iso_k')} | "
                           f"{rec['confidence']['tier']} (self-reported) |")
            print(f"✓ {f}")
    if summary:
        print("\n| formula | P (GPa) | λ | Tc Eliashberg (K) | tier |\n|---|---|---|---|---|")
        print("\n".join(summary))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
