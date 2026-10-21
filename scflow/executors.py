"""How QE/EPW binaries are launched on the machine the workflow runs on."""

import os
import shlex
import subprocess
from pathlib import Path

from .postprocess import job_done


class StepFailed(RuntimeError):
    pass


class LocalExecutor:
    """Runs QE binaries from `qe_bin`, wrapped in MPI when nproc > 1.

    SCFLOW_MPIRUN overrides the launcher (e.g. "srun" on Slurm, or a full path),
    SCFLOW_MPIRUN_ARGS adds launcher flags (e.g. "--bind-to none").
    """

    def __init__(self, qe_bin, mpirun=None, mpirun_args=None):
        self.qe_bin = Path(qe_bin)
        self.mpirun = mpirun or os.environ.get("SCFLOW_MPIRUN", "mpirun")
        self.mpirun_args = shlex.split(mpirun_args if mpirun_args is not None
                                       else os.environ.get("SCFLOW_MPIRUN_ARGS", ""))

    def command(self, binary, infile, nproc=1, extra_args=None):
        cmd = [str(self.qe_bin / binary)] + list(extra_args or []) + ["-in", str(Path(infile).resolve())]
        if nproc > 1:
            cmd = [self.mpirun, *self.mpirun_args, "-n", str(nproc)] + cmd
        return cmd

    def run_qe(self, binary, infile, outfile, cwd, nproc=1, extra_args=None):
        cmd = self.command(binary, infile, nproc, extra_args)
        # One BLAS/OpenMP thread per MPI rank unless the user asks otherwise; threaded
        # OpenBLAS inside every rank oversubscribes the cores badly.
        env = os.environ.copy()
        if nproc > 1:
            env.setdefault("OMP_NUM_THREADS", "1")
            env.setdefault("OPENBLAS_NUM_THREADS", "1")
        with open(outfile, "w") as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=cwd, env=env)
        # EPW can exit non-zero from IEEE floating-point exception reporting after
        # finishing correctly, so completion is judged from the output file.
        if not job_done(outfile):
            raise StepFailed(f"{binary} did not complete (exit {result.returncode}); see {outfile}")
        return result.returncode

    def run_pp_py(self, prefix, cwd):
        pp_py = self.qe_bin.parent / "EPW" / "bin" / "pp.py"
        if not pp_py.exists():
            raise StepFailed(f"EPW pp.py not found at {pp_py}")
        return subprocess.run(["python3", str(pp_py)], input=f"{prefix}\n", text=True,
                              capture_output=True, cwd=cwd)

    def missing_binaries(self):
        needed = ["pw.x", "ph.x", "q2r.x", "matdyn.x", "epw.x"]
        return [b for b in needed if not (self.qe_bin / b).exists()]
