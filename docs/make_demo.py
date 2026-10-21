"""Record docs/demo.cast by typing and running real scflow commands, then render it with agg.

    python docs/make_demo.py               # writes docs/demo.cast
    agg --idle-time-limit 1.5 docs/demo.cast docs/demo.gif

Every line of output in the cast is genuine output of the command shown.
"""

import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WIDTH, HEIGHT = 150, 40

SCENES = [
    "# scflow: crystal structure → validated superconducting Tc (Quantum ESPRESSO + EPW)",
    "# 1. A structure that passes symmetry checks but is wrong: H on Wyckoff 12e instead of 12d",
    "scflow check examples/CaH6_wrong_site.yaml",
    "# Rejected in under a second, before any compute. The corrected structure:",
    "scflow check compounds/CaH6_200GPa.yaml",
    "# 2. Live DFT in the container: vc-relax + SCF for CaH6 at 200 GPa, then the physics gate",
    "scflow submit compounds/CaH6_200GPa.yaml --backend docker --runs-dir runs/demo --nproc 16 --stop-after scf",
    "# 3. A full Nb run (computed earlier in the same container) validated and recorded in the DB",
    "scflow validate compounds/Nb_0GPa.yaml --run-dir runs/Nb --record --db runs/demo/scflow.db",
    "scflow db --db runs/demo/scflow.db list",
    "# Nb: Eliashberg Tc far above experiment → PRELIMINARY until the grids are converged",
]


class Cast:
    def __init__(self):
        self.t0 = time.time()
        self.events = []

    def out(self, text, at=None):
        t = (time.time() - self.t0) if at is None else at
        self.events.append([round(t, 4), "o", text.replace("\n", "\r\n")])

    def type(self, text):
        self.out("\x1b[1;32m$\x1b[0m ")
        for ch in text:
            self.out(ch)
            time.sleep(random.uniform(0.012, 0.035))
        time.sleep(0.4)
        self.out("\n")

    def comment(self, text):
        self.out("\x1b[2m")
        for ch in text:
            self.out(ch)
            time.sleep(0.008)
        self.out("\x1b[0m\n")
        time.sleep(1.2)

    def run(self, cmd):
        proc = subprocess.Popen(cmd, shell=True, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                env=os.environ | {"COLUMNS": str(WIDTH)}, bufsize=0)
        for line in iter(proc.stdout.readline, b""):
            self.out(line.decode(errors="replace"))
        proc.wait()
        self.out("\n")
        time.sleep(1.5)

    def save(self, path):
        header = {"version": 2, "width": WIDTH, "height": HEIGHT, "title": "scflow demo",
                  "env": {"TERM": "xterm-256color", "SHELL": "/bin/bash"}}
        with open(path, "w") as f:
            f.write(json.dumps(header) + "\n")
            for e in self.events:
                f.write(json.dumps(e) + "\n")


def main():
    random.seed(0)
    cast = Cast()
    for cmd in SCENES:
        if cmd.startswith("#"):
            cast.comment(cmd)
        else:
            cast.type(cmd)
            cast.run(cmd)
    out = ROOT / "docs" / "demo.cast"
    cast.save(out)
    print(f"wrote {out} ({len(cast.events)} events, {cast.events[-1][0]:.0f} s)", file=sys.stderr)


if __name__ == "__main__":
    main()
