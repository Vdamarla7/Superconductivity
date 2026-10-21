"""Build the static results page (site/index.html) from records/ and verified.txt."""

import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    verified = set()
    vfile = ROOT / "verified.txt"
    if vfile.exists():
        verified = {line.split()[0] for line in vfile.read_text().splitlines() if line.strip() and not line.startswith("#")}
    rows = []
    for f in sorted((ROOT / "records").rglob("*.json")):
        rec = json.loads(f.read_text())
        run, res = rec["run"], rec["results"]
        status = "verified" if run["id"] in verified else "self-reported, unverified"
        rows.append((run["formula"], run["pressure_gpa"], res.get("lambda"), res.get("tc_eliashberg_iso_k"),
                     res.get("tc_allen_dynes_k"), rec["confidence"]["tier"], status, run.get("scflow_version"),
                     f.relative_to(ROOT).as_posix()))
    rows.sort(key=lambda r: (r[0], r[1]))

    def cell(v, fmt="{:.2f}"):
        return "—" if v is None else (fmt.format(v) if isinstance(v, float) else html.escape(str(v)))

    body = "\n".join(
        f"<tr><td>{cell(r[0])}</td><td>{cell(r[1], '{:g}')}</td><td>{cell(r[2])}</td><td>{cell(r[3], '{:.1f}')}</td>"
        f"<td>{cell(r[4], '{:.1f}')}</td><td>{cell(r[5])}</td><td>{cell(r[6])}</td><td>{cell(r[7])}</td>"
        f"<td><a href='https://github.com/Vdamarla7/scflow-results/blob/main/{html.escape(r[8])}'>json</a></td></tr>"
        for r in rows)
    out = ROOT / "site"
    out.mkdir(exist_ok=True)
    (out / "index.html").write_text(f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>scflow results</title>
<style>
body {{ font: 15px/1.5 system-ui, sans-serif; margin: 2rem auto; max-width: 1100px; padding: 0 1rem; color: #222; }}
table {{ border-collapse: collapse; width: 100%; }} th, td {{ border-bottom: 1px solid #ddd; padding: .4rem .6rem; text-align: left; }}
th {{ background: #f4f4f4; }} p.note {{ color: #555; }}
</style></head><body>
<h1>scflow results</h1>
<p class="note">Superconducting Tc from first principles (Quantum ESPRESSO + EPW), submitted by people running
<a href="https://github.com/Vdamarla7/Superconductivity">scflow</a>. Records are self-reported: the tier comes from
scflow's checks on the submitter's machine and has not been reproduced unless marked verified.
{len(rows)} records. Data: CC BY 4.0.</p>
<table><thead><tr><th>Formula</th><th>P (GPa)</th><th>λ</th><th>Tc Eliashberg (K)</th><th>Tc Allen-Dynes (K)</th>
<th>Tier</th><th>Status</th><th>scflow</th><th></th></tr></thead>
<tbody>
{body}
</tbody></table></body></html>
""")
    print(f"wrote {out / 'index.html'} ({len(rows)} records)")


if __name__ == "__main__":
    main()
