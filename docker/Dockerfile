# syntax=docker/dockerfile:1
# scflow container: Quantum ESPRESSO 7.4.1 + EPW + Wannier90 3.1.0 + pslibrary (PBE) + scflow.
#
#   docker build -f docker/Dockerfile -t scflow:latest .
#
# Stages are ordered so the slow parts (QE compile, pseudopotential generation)
# stay cached while the Python package changes.

ARG UBUNTU=24.04
ARG QE_VERSION=7.4.1
ARG W90_VERSION=3.1.0

# ── Stage 1: compile QE + EPW ────────────────────────────────────────────────
FROM ubuntu:${UBUNTU} AS qe-build
ARG QE_VERSION
ARG W90_VERSION
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        gfortran gcc g++ make wget git ca-certificates python3 \
        libopenblas-dev libfftw3-dev libfftw3-mpi-dev \
        libopenmpi-dev openmpi-bin \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
RUN wget -q https://github.com/wannier-developers/wannier90/archive/v${W90_VERSION}.tar.gz -O w90.tar.gz \
 && tar xf w90.tar.gz && mv wannier90-${W90_VERSION} /opt/wannier90 \
 && cd /opt/wannier90 && cp config/make.inc.gfort make.inc && make -j"$(nproc)" libs

RUN wget -q https://gitlab.com/QEF/q-e/-/archive/qe-${QE_VERSION}/q-e-qe-${QE_VERSION}.tar.gz -O qe.tar.gz \
 && tar xf qe.tar.gz && mv q-e-qe-${QE_VERSION} /opt/qe

WORKDIR /opt/qe
RUN (./configure --with-scalapack=no WANNIER90_DIR=/opt/wannier90 > configure.log 2>&1; \
     test -f make.inc || { tail -60 configure.log; exit 1; }) \
 && make -j"$(nproc)" pw ph pp ld1 \
 && make -j"$(nproc)" epw \
 && mkdir -p /opt/qe-dist/bin && cp -L bin/*.x /opt/qe-dist/bin/ && strip /opt/qe-dist/bin/*.x \
 && ls -l /opt/qe-dist/bin/pw.x /opt/qe-dist/bin/ph.x /opt/qe-dist/bin/q2r.x /opt/qe-dist/bin/matdyn.x /opt/qe-dist/bin/epw.x

# ── Stage 2: generate pslibrary PBE pseudopotentials with ld1.x ──────────────
FROM qe-build AS pseudos
RUN git clone --depth=1 https://github.com/dalcorso/pslibrary /opt/pslibrary \
 && cd /opt/pslibrary && printf "#!/bin/bash\nPWDIR='/opt/qe/'\n" > QE_path \
 && cd pbe && bash -c ". ../make_ps" > make_ps.log 2>&1; \
    n=$(ls /opt/pslibrary/pbe/PSEUDOPOTENTIALS/*.UPF | wc -l); \
    echo "generated $n PBE pseudopotentials"; test "$n" -gt 100

# ── Stage 3: runtime ─────────────────────────────────────────────────────────
FROM ubuntu:${UBUNTU} AS runtime
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgfortran5 libopenblas0 libfftw3-double3 libfftw3-mpi3 \
        openmpi-bin libopenmpi3t64 \
        python3 python3-venv ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=qe-build /opt/qe-dist/bin /opt/qe/bin
COPY --from=qe-build /opt/qe/EPW/bin/pp.py /opt/qe/EPW/bin/pp.py
COPY --from=pseudos /opt/pslibrary/pbe/PSEUDOPOTENTIALS /opt/pseudo/pbe

RUN python3 -m venv /opt/venv
ENV PATH=/opt/venv/bin:/opt/qe/bin:$PATH \
    SCFLOW_QE_BIN=/opt/qe/bin \
    SCFLOW_PSEUDO_DIR=/opt/pseudo/pbe \
    PYTHONUNBUFFERED=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    OMPI_MCA_rmaps_base_oversubscribe=1 \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

WORKDIR /opt/scflow
COPY pyproject.toml README.md LICENSE ./
COPY scflow ./scflow
COPY compounds ./compounds
RUN pip install --no-cache-dir .

RUN userdel -r ubuntu 2>/dev/null; useradd -m -u 1000 scflow && mkdir -p /work && chown scflow /work
USER scflow
WORKDIR /work
ENTRYPOINT ["scflow"]
CMD ["--help"]
