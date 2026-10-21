#!/usr/bin/env bash
# Native (non-container) install of QE 7.4.1 + EPW + Wannier90 + pslibrary + scflow on Ubuntu/WSL2.
# The container (docker/Dockerfile) is the easier route; use this for bare-metal performance.
# Run once from a terminal where sudo is available:
#   bash scripts/setup_native.sh
#
# Takes ~30-60 min (dominated by QE compilation).
set -euo pipefail
JOBS=$(nproc)
QE_DIR="$HOME/work/quantum/qe-7.4.1"
W90_DIR="$HOME/work/quantum/wannier90-3.1.0"
PSEUDO_DIR="$HOME/work/quantum/pslibrary/pbe/PSEUDOPOTENTIALS"
PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== [1/7] Installing system packages ==="
sudo apt-get update -qq
sudo apt-get install -y \
  gfortran gcc g++ make wget git \
  python3-venv \
  libopenblas-dev libfftw3-dev libfftw3-mpi-dev \
  libopenmpi-dev openmpi-bin \
  libscalapack-openmpi-dev

echo "=== [2/7] Setting up Python venv ==="
python3 -m venv "$PIPELINE_DIR/.venv"
source "$PIPELINE_DIR/.venv/bin/activate"
pip install --upgrade pip --quiet
pip install -e "$PIPELINE_DIR[sscha]"
echo "Python packages installed. Activate with: source $PIPELINE_DIR/.venv/bin/activate"

echo "=== [3/7] Building Wannier90 3.1.0 (library) ==="
mkdir -p "$HOME/work/quantum"
cd "$HOME/work/quantum"
if [ ! -f "wannier90-3.1.0.tar.gz" ]; then
  wget -q https://github.com/wannier-developers/wannier90/archive/v3.1.0.tar.gz \
       -O wannier90-3.1.0.tar.gz
fi
tar xf wannier90-3.1.0.tar.gz
cd wannier90-3.1.0
cp config/make.inc.gfort make.inc
make -j"$JOBS" libs
echo "Wannier90 library built: $W90_DIR/libwannier.a"

echo "=== [4/7] Downloading QE 7.4.1 ==="
cd "$HOME/work/quantum"
if [ ! -f "q-e-qe-7.4.1.tar.gz" ]; then
  wget -q https://gitlab.com/QEF/q-e/-/archive/qe-7.4.1/q-e-qe-7.4.1.tar.gz
fi
tar xf q-e-qe-7.4.1.tar.gz
mv q-e-qe-7.4.1 qe-7.4.1

echo "=== [5/7] Configuring and building QE ==="
cd "$QE_DIR"
./configure \
  --with-scalapack=yes \
  WANNIER90_DIR="$W90_DIR"
make -j"$JOBS" all
echo "Core QE binaries done."

echo "=== [6/7] Building EPW ==="
cd "$QE_DIR"
make -j"$JOBS" epw
echo "epw.x built."
ls -lh bin/pw.x bin/ph.x bin/q2r.x bin/matdyn.x bin/epw.x

echo "=== [7/7] Downloading pslibrary pseudopotentials ==="
cd "$HOME/work/quantum"
if [ ! -d "pslibrary" ]; then
  git clone --depth=1 https://github.com/dalcorso/pslibrary
fi
cd pslibrary
# Generate PBE USPP pseudopotentials (requires pw.x in PATH)
export PATH="$QE_DIR/bin:$PATH"
./make_all_ps 2>&1 | tee pslibrary_build.log
mkdir -p "$PSEUDO_DIR"
find . -name "*.UPF" -exec cp {} "$PSEUDO_DIR/" \;
echo "Pseudopotentials copied to $PSEUDO_DIR"

echo ""
echo "=== Setup complete ==="
echo "Add to your shell profile:"
echo "  source $PIPELINE_DIR/.venv/bin/activate"
echo "  export SCFLOW_QE_BIN=$QE_DIR/bin"
echo "  export SCFLOW_PSEUDO_DIR=$PSEUDO_DIR"
echo ""
echo "Then check the install:  scflow doctor compounds/Nb_0GPa.yaml"
