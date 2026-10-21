"""scflow command line: run / check / validate / submit / db / doctor."""

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

from . import __version__
from . import backends, db
from .config import ConfigError, default_outdir, load_config
from .executors import LocalExecutor, StepFailed
from .validation import Report, check_structure, format_report, validate
from .workflow import STEPS, Workflow

DEFAULT_DB = os.environ.get("SCFLOW_DB", "runs/scflow.db")


def _load(args, need_paths=True):
    path = args.config
    if path.startswith("s3://"):
        path = backends.s3_download(path, tempfile.mkdtemp(prefix="scflow-"))
    if not need_paths:
        os.environ.setdefault("SCFLOW_QE_BIN", "unused")
        os.environ.setdefault("SCFLOW_PSEUDO_DIR", "unused")
    return load_config(path, qe_bin=getattr(args, "qe_bin", None), pseudo_dir=getattr(args, "pseudo_dir", None),
                       nproc=getattr(args, "nproc", None))


def cmd_check(args):
    cfg = _load(args, need_paths=False)
    report = Report()
    check_structure(cfg, report)
    print(format_report(report))
    return 1 if report.worst == "fail" else 0


def cmd_run(args):
    cfg = _load(args)
    pre = Report()
    check_structure(cfg, pre)
    print("\n  Pre-flight structure checks")
    print(format_report(pre))
    if pre.worst == "fail" and not args.force:
        print("\n  ✗ Refusing to start: a structure check failed. Fix the config or pass --force.")
        return 2

    outdir = Path(args.outdir) if args.outdir else default_outdir(cfg)
    ex = LocalExecutor(cfg["paths"]["qe_bin"])
    missing = ex.missing_binaries()
    if missing and not args.dry_run:
        print(f"  ✗ QE binaries not found in {cfg['paths']['qe_bin']}: {', '.join(missing)}")
        return 2
    wf = Workflow(cfg, outdir, ex, dry_run=args.dry_run)
    code = 0
    try:
        wf.run(start_from=args.start_from, stop_after=args.stop_after, reset_epw=args.reset_epw)
    except StepFailed as e:
        print(f"\n  ✗ {e}")
        code = 1
    if not args.dry_run:
        code = max(code, _validate_and_record(cfg, outdir, args))
    if args.s3_prefix:
        backends.s3_upload_dir(outdir, args.s3_prefix)
        backends.s3_upload_dir(Path(args.db).parent, args.s3_prefix, patterns=("scflow.db", "*.record.json"))
    return code


def _validate_and_record(cfg, run_dir, args):
    report = validate(cfg, run_dir)
    print("\n  Validation")
    print(format_report(report))
    finished = any(c.name == "epw.completed" and c.status == "pass" for c in report.checks)
    if not finished or getattr(args, "no_db", False):
        return 1 if report.worst == "fail" else 0
    con = db.connect(args.db)
    run_id = db.record_run(con, cfg, run_dir, report, backend=args.backend_label, host=os.uname().nodename)
    tier, why = db.confidence_tier(con, run_id)
    record_path = Path(args.db).parent / f"{cfg['compound']['prefix']}-{run_id[:8]}.record.json"
    record_path.write_text(json.dumps(db.export_record(con, run_id), indent=2, default=str))
    print(f"\n  Recorded run {run_id[:8]} in {args.db}  →  confidence: {tier.upper()} ({why})")
    return 1 if report.worst == "fail" else 0


def cmd_validate(args):
    cfg = _load(args, need_paths=False)
    run_dir = Path(args.run_dir) if args.run_dir else default_outdir(cfg)
    if args.json:
        print(json.dumps(validate(cfg, run_dir).to_dict(), indent=2, default=str))
        return 0
    if args.record:
        return _validate_and_record(cfg, run_dir, args)
    report = validate(cfg, run_dir)
    print(format_report(report))
    return 1 if report.worst == "fail" else 0


def cmd_submit(args):
    cfg = _load(args, need_paths=False)
    prefix = cfg["compound"]["prefix"]
    if args.backend == "docker":
        extra = ["--reset-epw"] if args.reset_epw else []
        result = backends.docker(args.config, prefix, args.runs_dir, image=args.image, nproc=args.nproc,
                                 start_from=args.start_from, stop_after=args.stop_after, extra=extra,
                                 dry_run=args.dry_run)
    elif args.backend == "slurm":
        result = backends.slurm(args.config, prefix, args.runs_dir, nproc=args.nproc or cfg["parallel"]["nproc"],
                                time=args.time, partition=args.partition, account=args.account, image=args.image,
                                native=args.native, start_from=args.start_from, stop_after=args.stop_after,
                                extra=["--reset-epw"] if args.reset_epw else [], dry_run=args.dry_run)
    else:
        if not (args.bucket and args.job_queue and args.job_definition):
            print("aws-batch needs --bucket, --job-queue and --job-definition (see deploy/aws/README.md)")
            return 2
        result = backends.aws_batch(args.config, prefix, args.bucket, args.job_queue, args.job_definition,
                                    nproc=args.nproc, start_from=args.start_from, region=args.region,
                                    dry_run=args.dry_run)
    if isinstance(result, int):
        return result
    print(json.dumps(result, indent=2) if isinstance(result, dict) else result)
    return 0


def cmd_db(args):
    con = db.connect(args.db)
    if args.db_cmd == "list":
        rows = db.list_runs(con)
        if not rows:
            print(f"  no runs in {args.db}")
            return 0
        print(f"  {'id':<9}{'compound':<12}{'P/GPa':>7}  {'fine k / q':<28}{'μ*':>5}{'λ':>7}{'Tc_El/K':>9}{'Tc_AD/K':>9}"
              f"{'ref/K':>7}  tier")
        for r in rows:
            f = lambda x, p=1: "—" if x is None else f"{x:.{p}f}"
            print(f"  {r['id'][:8]:<9}{r['prefix']:<12}{r['pressure_gpa']:>7.0f}  {r['fine_k']:<28}{r['muc']:>5.2f}"
                  f"{f(r['lambda'], 2):>7}{f(r['tc_eliashberg_iso_k']):>9}{f(r['tc_allen_dynes_k']):>9}"
                  f"{f(r['tc_reference_k'], 1):>7}  {r['tier']}")
    elif args.db_cmd == "show":
        row = con.execute("SELECT id FROM runs WHERE id LIKE ?", (args.run_id + "%",)).fetchone()
        if not row:
            print(f"  no run matching {args.run_id}")
            return 1
        print(json.dumps(db.export_record(con, row["id"]), indent=2, default=str))
    elif args.db_cmd == "export":
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        ids = [r["id"] for r in con.execute("SELECT id FROM runs")]
        for rid in ids:
            rec = db.export_record(con, rid)
            (out / f"{rec['run']['prefix']}-{rid[:8]}.json").write_text(json.dumps(rec, indent=2, default=str))
        print(f"  exported {len(ids)} records to {out}/")
    return 0


def cmd_doctor(args):
    ok = True
    qe_bin = args.qe_bin or os.environ.get("SCFLOW_QE_BIN")
    pseudo = args.pseudo_dir or os.environ.get("SCFLOW_PSEUDO_DIR")
    print(f"  scflow {__version__}, python {sys.version.split()[0]}")
    if qe_bin:
        missing = LocalExecutor(qe_bin).missing_binaries()
        print(f"  {'✓' if not missing else '✗'} QE binaries in {qe_bin}" + (f" (missing {missing})" if missing else ""))
        ok &= not missing
        pp_py = Path(qe_bin).parent / "EPW" / "bin" / "pp.py"
        print(f"  {'✓' if pp_py.exists() else '✗'} EPW pp.py at {pp_py}")
        ok &= pp_py.exists()
    else:
        print("  ✗ SCFLOW_QE_BIN not set")
        ok = False
    launcher = os.environ.get("SCFLOW_MPIRUN", "mpirun")
    print(f"  {'✓' if shutil.which(launcher) else '✗'} MPI launcher '{launcher}'")
    ok &= bool(shutil.which(launcher))
    if args.config:
        args.nproc = None
        try:
            cfg = _load(args)
        except ConfigError as e:
            print(f"  ✗ {e}")
            return 1
        for sp in cfg["structure"]["species"]:
            p = Path(cfg["paths"]["pseudo_dir"]) / sp["pseudo"]
            print(f"  {'✓' if p.exists() else '✗'} pseudopotential {sp['pseudo']}")
            ok &= p.exists()
    elif pseudo:
        n = len(list(Path(pseudo).glob("*.UPF")))
        print(f"  {'✓' if n else '✗'} {n} pseudopotentials in {pseudo}")
        ok &= n > 0
    for mod in ("spglib", "ase", "numpy", "yaml"):
        try:
            __import__(mod)
            print(f"  ✓ python module {mod}")
        except ImportError:
            print(f"  ✗ python module {mod}")
            ok = False
    return 0 if ok else 1


def build_parser():
    p = argparse.ArgumentParser(prog="scflow", description="First-principles superconducting Tc: QE + EPW, validated.")
    p.add_argument("--version", action="version", version=f"scflow {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    def paths(sp):
        sp.add_argument("--qe-bin", help="directory with pw.x, ph.x, epw.x (default: $SCFLOW_QE_BIN)")
        sp.add_argument("--pseudo-dir", help="pseudopotential directory (default: $SCFLOW_PSEUDO_DIR)")

    r = sub.add_parser("run", help="run the workflow on this machine")
    r.add_argument("config")
    r.add_argument("--outdir")
    r.add_argument("--start-from", choices=STEPS, default="vcrelax")
    r.add_argument("--stop-after", choices=STEPS)
    r.add_argument("--reset-epw", action="store_true", help="delete stale EPW files (after changing Wannier windows)")
    r.add_argument("--dry-run", action="store_true", help="write inputs only")
    r.add_argument("--nproc", type=int)
    r.add_argument("--force", action="store_true", help="run even if a pre-flight structure check fails")
    r.add_argument("--db", default=DEFAULT_DB)
    r.add_argument("--no-db", action="store_true")
    r.add_argument("--backend-label", default="local", help=argparse.SUPPRESS)
    r.add_argument("--s3-prefix", help=argparse.SUPPRESS)
    paths(r)
    r.set_defaults(func=cmd_run)

    c = sub.add_parser("check", help="pre-flight structure checks only (no QE needed)")
    c.add_argument("config")
    c.set_defaults(func=cmd_check)

    v = sub.add_parser("validate", help="run the physics checks on an existing run directory")
    v.add_argument("config")
    v.add_argument("--run-dir")
    v.add_argument("--record", action="store_true", help="also record into the database")
    v.add_argument("--json", action="store_true")
    v.add_argument("--db", default=DEFAULT_DB)
    v.add_argument("--backend-label", default="local", help=argparse.SUPPRESS)
    v.set_defaults(func=cmd_validate)

    s = sub.add_parser("submit", help="launch a run in a container, on Slurm, or on AWS Batch")
    s.add_argument("config")
    s.add_argument("--backend", choices=["docker", "slurm", "aws-batch"], required=True)
    s.add_argument("--runs-dir", default="runs")
    s.add_argument("--image", default=backends.DEFAULT_IMAGE)
    s.add_argument("--nproc", type=int)
    s.add_argument("--start-from", choices=STEPS)
    s.add_argument("--stop-after", choices=STEPS)
    s.add_argument("--reset-epw", action="store_true")
    s.add_argument("--dry-run", action="store_true", help="print/write what would be submitted")
    s.add_argument("--time", default="24:00:00", help="slurm wall time")
    s.add_argument("--partition")
    s.add_argument("--account")
    s.add_argument("--native", action="store_true", help="slurm: use a native scflow+QE install instead of Apptainer")
    s.add_argument("--bucket")
    s.add_argument("--job-queue")
    s.add_argument("--job-definition")
    s.add_argument("--region")
    s.set_defaults(func=cmd_submit)

    d = sub.add_parser("db", help="query the results database")
    d.add_argument("--db", default=DEFAULT_DB)
    dsub = d.add_subparsers(dest="db_cmd", required=True)
    dsub.add_parser("list")
    show = dsub.add_parser("show")
    show.add_argument("run_id")
    exp = dsub.add_parser("export")
    exp.add_argument("--out", default="records")
    d.set_defaults(func=cmd_db)

    doc = sub.add_parser("doctor", help="check that this machine can run scflow")
    doc.add_argument("config", nargs="?")
    paths(doc)
    doc.set_defaults(func=cmd_doctor)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except ConfigError as e:
        print(f"  ✗ {e}")
        return 2


if __name__ == "__main__":
    sys.exit(main())
