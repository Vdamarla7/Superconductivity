"""Launch a whole `scflow run` somewhere: a local container, a Slurm cluster, or AWS Batch.

Each backend runs the same container image and the same `scflow run` command, so
checkpointing, validation and the database record behave identically everywhere.
"""

import os
import shlex
import shutil
import subprocess
from pathlib import Path

DEFAULT_IMAGE = os.environ.get("SCFLOW_IMAGE", "scflow:latest")


def _run_args(config_name, prefix, nproc, start_from=None, stop_after=None, extra=()):
    args = ["run", f"/input/{config_name}", "--outdir", f"/work/runs/{prefix}", "--db", "/work/runs/scflow.db",
            "--backend-label", "{label}"]
    if nproc:
        args += ["--nproc", str(nproc)]
    if start_from:
        args += ["--start-from", start_from]
    if stop_after:
        args += ["--stop-after", stop_after]
    return args + list(extra)


def docker(config_path, prefix, runs_dir, image=DEFAULT_IMAGE, nproc=None, start_from=None, stop_after=None,
           extra=(), dry_run=False):
    config_path, runs_dir = Path(config_path).resolve(), Path(runs_dir).resolve()
    runs_dir.mkdir(parents=True, exist_ok=True)
    args = [a.format(label="docker") for a in _run_args(config_path.name, prefix, nproc, start_from, stop_after, extra)]
    cmd = ["docker", "run", "--rm", "--shm-size=4g",
           "--user", f"{os.getuid()}:{os.getgid()}",
           "-v", f"{config_path.parent}:/input:ro", "-v", f"{runs_dir}:/work/runs", image] + args
    if dry_run:
        return shlex.join(cmd)
    return subprocess.call(cmd)


SLURM_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=scflow-{prefix}
#SBATCH --nodes=1
#SBATCH --ntasks={nproc}
#SBATCH --time={time}
#SBATCH --output={runs_dir}/{prefix}.slurm-%j.log
{partition_line}{account_line}
set -euo pipefail
mkdir -p {runs_dir}
{body}
"""


def slurm(config_path, prefix, runs_dir, nproc, time="24:00:00", partition=None, account=None,
          image=DEFAULT_IMAGE, native=False, start_from=None, stop_after=None, extra=(), dry_run=False):
    """Write (and unless dry_run, sbatch) a job script.

    native=True uses an `scflow` + QE install already on the cluster (MPI via srun);
    otherwise the container image runs under Apptainer, which most HPC centres provide.
    """
    config_path, runs_dir = Path(config_path).resolve(), Path(runs_dir).resolve()
    if native:
        args = [a.format(label="slurm") for a in _run_args(config_path.name, prefix, nproc, start_from, stop_after, extra)]
        args = [a.replace("/input/", f"{config_path.parent}/").replace("/work/runs", str(runs_dir)) for a in args]
        body = "export SCFLOW_MPIRUN=srun\nscflow " + shlex.join(args)
    else:
        sif = f"{runs_dir}/scflow.sif"
        args = [a.format(label="slurm-apptainer")
                for a in _run_args(config_path.name, prefix, nproc, start_from, stop_after, extra)]
        image_ref = image if "://" in image else f"docker-daemon://{image}"
        body = (f"[ -f {sif} ] || apptainer build {sif} {image_ref}\n"
                f"apptainer run --bind {config_path.parent}:/input:ro --bind {runs_dir}:/work/runs {sif} "
                + shlex.join(args))
    script = SLURM_TEMPLATE.format(
        prefix=prefix, nproc=nproc, time=time, runs_dir=runs_dir, body=body,
        partition_line=f"#SBATCH --partition={partition}\n" if partition else "",
        account_line=f"#SBATCH --account={account}\n" if account else "")
    runs_dir.mkdir(parents=True, exist_ok=True)
    path = runs_dir / f"{prefix}.sbatch"
    path.write_text(script)
    if dry_run or not shutil.which("sbatch"):
        return f"wrote {path}" + ("" if dry_run else " (sbatch not found — submit it on the cluster)")
    out = subprocess.run(["sbatch", str(path)], capture_output=True, text=True, check=True)
    return out.stdout.strip()


def aws_batch(config_path, prefix, bucket, job_queue, job_definition, nproc=None, start_from=None,
              region=None, dry_run=False):
    """Upload the config to S3 and submit an AWS Batch job that runs scflow and syncs results back.

    Infrastructure (compute environment, queue, job definition, bucket) comes from
    deploy/aws/batch.yaml. The container's `run --s3-prefix` uploads the run directory,
    the SQLite DB and the JSON record when the job ends.
    """
    config_path = Path(config_path).resolve()
    key = f"scflow/{prefix}/input/{config_path.name}"
    s3_prefix = f"s3://{bucket}/scflow/{prefix}"
    command = ["run", f"s3://{bucket}/{key}", "--outdir", f"/work/runs/{prefix}", "--db", "/work/runs/scflow.db",
               "--backend-label", "aws-batch", "--s3-prefix", s3_prefix]
    if nproc:
        command += ["--nproc", str(nproc)]
    if start_from:
        command += ["--start-from", start_from]
    request = {"jobName": f"scflow-{prefix}", "jobQueue": job_queue, "jobDefinition": job_definition,
               "containerOverrides": {"command": command}}
    if dry_run:
        return {"upload": f"{config_path} -> s3://{bucket}/{key}", "submit_job": request}
    import boto3

    boto3.client("s3", region_name=region).upload_file(str(config_path), bucket, key)
    resp = boto3.client("batch", region_name=region).submit_job(**request)
    return {"jobId": resp["jobId"], "results": s3_prefix}


def s3_download(uri, dest_dir):
    import boto3

    bucket, key = uri[5:].split("/", 1)
    dest = Path(dest_dir) / Path(key).name
    boto3.client("s3").download_file(bucket, key, str(dest))
    return dest


def s3_upload_dir(local_dir, s3_prefix, patterns=("*.in", "*.out", "*.wout", "*.a2f*", "*.json", "*.yaml", "*.db")):
    import boto3

    bucket, prefix = s3_prefix[5:].split("/", 1)
    s3 = boto3.client("s3")
    local_dir = Path(local_dir)
    for pat in patterns:
        for p in local_dir.rglob(pat):
            if p.is_file() and "tmp" not in p.relative_to(local_dir).parts:
                s3.upload_file(str(p), bucket, f"{prefix}/{p.relative_to(local_dir)}")
