"""
QE and EPW input file generators.
Each function takes a config dict (loaded from compound YAML) and returns a string.
"""

from kpoints import gen_kpoints_crystal, gen_kpoints_automatic


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _cell_section(cfg):
    """Return the ibrav/A or ibrav=0 lines for the &SYSTEM block."""
    s = cfg["structure"]
    if s["ibrav"] in (2, 3):
        if s["A"] is None:
            raise ValueError(
                "structure.A must be set (provide an initial guess for vc-relax; "
                "the pipeline will update it with the relaxed value after step 1)"
            )
        return f"  ibrav    = {s['ibrav']}\n  A        = {s['A']:.6f}"
    else:
        return "  ibrav    = 0"


def _cell_parameters_block(cfg):
    """Return CELL_PARAMETERS {angstrom} block. Only needed for ibrav=0."""
    s = cfg["structure"]
    if s["ibrav"] != 0:
        return ""
    cp = s["cell_parameters_angstrom"]
    lines = ["CELL_PARAMETERS {angstrom}"]
    for row in cp:
        lines.append(f"  {row[0]:14.9f}  {row[1]:14.9f}  {row[2]:14.9f}")
    return "\n".join(lines) + "\n"


def _atomic_species_block(cfg):
    lines = ["ATOMIC_SPECIES"]
    for sp in cfg["structure"]["species"]:
        lines.append(f"  {sp['symbol']:<4s}  {sp['mass']:.5f}  {sp['pseudo']}")
    return "\n".join(lines)


def _atomic_positions_block(cfg):
    lines = ["ATOMIC_POSITIONS {crystal}"]
    for pos in cfg["structure"]["positions"]:
        lines.append(f"  {pos[0]:<4s}  {float(pos[1]):.6f}  {float(pos[2]):.6f}  {float(pos[3]):.6f}")
    return "\n".join(lines)


def _projections_block(cfg):
    lines = []
    for i, proj in enumerate(cfg["wannier"]["projections"], 1):
        lines.append(f"  proj({i})       = '{proj}'")
    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Step 1: vc-relax
# ─────────────────────────────────────────────────────────────────────────────

def gen_vcrelax(cfg):
    """Variable-cell relaxation to find equilibrium geometry at target pressure."""
    p = cfg["compound"]
    s = cfg["structure"]
    d = cfg["dft"]
    press_kbar = p["pressure_GPa"] * 10.0

    # For vc-relax, use the initial guess A (or cell params) from the config
    cell_lines = _cell_section(cfg)
    cell_params = _cell_parameters_block(cfg)

    return f"""\
&CONTROL
  calculation   = 'vc-relax'
  restart_mode  = 'from_scratch'
  prefix        = '{p['prefix']}'
  outdir        = './tmp/'
  pseudo_dir    = '{cfg['paths']['pseudo_dir']}'
  tprnfor       = .true.
  tstress       = .true.
  etot_conv_thr = 1.0d-6
  forc_conv_thr = 1.0d-5
/

&SYSTEM
{cell_lines}
  nat      = {s['nat']}
  ntyp     = {s['ntyp']}
  input_dft = 'PBE'
  ecutwfc  = {d['ecutwfc']}
  ecutrho  = {d['ecutrho']}
  occupations = 'smearing'
  smearing    = 'mp'
  degauss     = {d['degauss']}
/

&ELECTRONS
  conv_thr    = 1.0d-9
  mixing_beta = 0.4
  diagonalization = 'david'
/

&IONS
  ion_dynamics = 'bfgs'
/

&CELL
  cell_dynamics  = 'bfgs'
  press          = {press_kbar:.1f}
  press_conv_thr = 0.5
/

{_atomic_species_block(cfg)}

{cell_params}{_atomic_positions_block(cfg)}

{gen_kpoints_automatic(*d['scf_kgrid'])}
"""


# ─────────────────────────────────────────────────────────────────────────────
# Step 2: SCF
# ─────────────────────────────────────────────────────────────────────────────

def gen_scf(cfg):
    """Ground-state SCF at the relaxed geometry."""
    p = cfg["compound"]
    s = cfg["structure"]
    d = cfg["dft"]
    cell_lines = _cell_section(cfg)
    cell_params = _cell_parameters_block(cfg)

    return f"""\
&CONTROL
  calculation   = 'scf'
  restart_mode  = 'from_scratch'
  prefix        = '{p['prefix']}'
  outdir        = './tmp/'
  pseudo_dir    = '{cfg['paths']['pseudo_dir']}'
  tprnfor       = .true.
  tstress       = .true.
/

&SYSTEM
{cell_lines}
  nat      = {s['nat']}
  ntyp     = {s['ntyp']}
  input_dft = 'PBE'
  ecutwfc  = {d['ecutwfc']}
  ecutrho  = {d['ecutrho']}
  occupations = 'smearing'
  smearing    = 'mp'
  degauss     = {d['degauss']}
/

&ELECTRONS
  conv_thr    = 1.0d-12
  mixing_beta = 0.4
  diagonalization = 'david'
/

{_atomic_species_block(cfg)}

{cell_params}{_atomic_positions_block(cfg)}

{gen_kpoints_automatic(*d['scf_kgrid'])}
"""


# ─────────────────────────────────────────────────────────────────────────────
# Step 3: ph.x phonons
# ─────────────────────────────────────────────────────────────────────────────

def gen_phonon(cfg):
    """Harmonic DFPT phonons on a uniform q-grid."""
    p = cfg["compound"]
    d = cfg["dft"]
    q = d["coarse_kgrid"]

    return f"""\
Phonon calculation for {p['prefix']} at {p['pressure_GPa']} GPa
&INPUTPH
  prefix        = '{p['prefix']}'
  outdir        = './tmp/'
  fildyn        = '{p['prefix']}.dyn'
  fildvscf      = 'dvscf'
  ldisp         = .true.
  trans         = .true.
  nq1           = {q[0]}
  nq2           = {q[1]}
  nq3           = {q[2]}
  tr2_ph        = 1.0d-14
  alpha_mix(1)  = 0.3
  nmix_ph       = 4
/
"""


# ─────────────────────────────────────────────────────────────────────────────
# Step 4a: q2r
# ─────────────────────────────────────────────────────────────────────────────

def gen_q2r(cfg):
    p = cfg["compound"]
    return f"""\
&INPUT
  fildyn  = '{p['prefix']}.dyn'
  zasr    = 'simple'
  flfrc   = '{p['prefix']}.fc'
/
"""


# ─────────────────────────────────────────────────────────────────────────────
# Step 4b: matdyn (phonon DOS)
# ─────────────────────────────────────────────────────────────────────────────

def gen_matdyn(cfg):
    p = cfg["compound"]
    return f"""\
&INPUT
  asr    = 'simple'
  flfrc  = '{p['prefix']}.fc'
  flfrq  = '{p['prefix']}.freq'
  fldos  = '{p['prefix']}.dos'
  dos    = .true.
  nk1    = 20
  nk2    = 20
  nk3    = 20
  ndos   = 500
/
"""


# ─────────────────────────────────────────────────────────────────────────────
# Step 5: NSCF for EPW coarse k-grid
# ─────────────────────────────────────────────────────────────────────────────

def gen_nscf(cfg):
    """
    NSCF with fully explicit k-list — required by EPW.
    K_POINTS {automatic} gives only symmetry-reduced points; EPW needs all N³.
    """
    p = cfg["compound"]
    s = cfg["structure"]
    d = cfg["dft"]
    w = cfg["wannier"]
    cell_lines = _cell_section(cfg)
    cell_params = _cell_parameters_block(cfg)

    n1, n2, n3 = d["coarse_kgrid"]
    kpoints = gen_kpoints_crystal(n1, n2, n3)

    return f"""\
&CONTROL
  calculation   = 'nscf'
  restart_mode  = 'from_scratch'
  prefix        = '{p['prefix']}'
  outdir        = './tmp/'
  pseudo_dir    = '{cfg['paths']['pseudo_dir']}'
/

&SYSTEM
{cell_lines}
  nat      = {s['nat']}
  ntyp     = {s['ntyp']}
  input_dft = 'PBE'
  ecutwfc  = {d['ecutwfc']}
  ecutrho  = {d['ecutrho']}
  occupations = 'smearing'
  smearing    = 'mp'
  degauss     = {d['degauss']}
  nbnd        = {w.get('nbnd', 20)}
/

&ELECTRONS
  conv_thr    = 1.0d-12
  mixing_beta = 0.4
  diagonalization = 'david'
/

{_atomic_species_block(cfg)}

{cell_params}{_atomic_positions_block(cfg)}

{kpoints}
"""


# ─────────────────────────────────────────────────────────────────────────────
# Step 6: EPW pass 1 — Wannierization only, writes .epmatwp
# ─────────────────────────────────────────────────────────────────────────────

def gen_epw_pass1(cfg):
    """
    Wannierize and write .epmatwp (coarse e-ph matrices in Wannier basis).
    Uses coarse fine grid to finish quickly — real interpolation done in pass 2.
    """
    p = cfg["compound"]
    s = cfg["structure"]
    d = cfg["dft"]
    w = cfg["wannier"]
    e = cfg["epw"]

    nk = d["coarse_kgrid"]
    bands_skipped_line = ""
    if w.get("bands_skipped"):
        bands_skipped_line = f"  bands_skipped = '{w['bands_skipped']}'\n"

    return f"""\
&inputepw
  prefix        = '{p['prefix']}'
{_amass_block(cfg)}
  outdir        = './tmp/'

  ! Pass 1: Wannierize + write .epmatwp only
  ep_coupling   = .true.
  elph          = .true.
  epbwrite      = .true.
  epbread       = .false.
  epwwrite      = .true.
  epwread       = .false.
  ephwrite      = .false.
  kmaps         = .false.

  wannierize    = .true.
  num_iter      = 500
  nbndsub       = {w['nbndsub']}
{bands_skipped_line}{_projections_block(cfg)}

  dis_win_min   = {w['dis_win_min']}
  dis_win_max   = {w['dis_win_max']}
  dis_froz_min  = {w['dis_froz_min']}
  dis_froz_max  = {w['dis_froz_max']}

  phonselfen    = .false.
  a2f           = .false.
  eliashberg    = .false.
  fsthick       = {e['fsthick']}

  nk1 = {nk[0]}  nk2 = {nk[1]}  nk3 = {nk[2]}
  nq1 = {nk[0]}  nq2 = {nk[1]}  nq3 = {nk[2]}
  mp_mesh_k     = .false.

  ! Coarse fine grid for pass 1 — just need Wannierization done
  nkf1 = {nk[0]}  nkf2 = {nk[1]}  nkf3 = {nk[2]}
  nqf1 = {nk[0]}  nqf2 = {nk[1]}  nqf3 = {nk[2]}

  dvscf_dir     = '../2_phonons/save'
  iverbosity    = 0
/
"""


# ─────────────────────────────────────────────────────────────────────────────
# Step 7: EPW pass 2 — fine grid interpolation + Eliashberg
# ─────────────────────────────────────────────────────────────────────────────

def gen_epw_pass2(cfg):
    """
    Read .epmatwp from pass 1. Interpolate onto fine k/q grid.
    Compute phonon self-energy and solve Eliashberg equations.

    Critical flags (learned from debugging):
      liso = .true.   → correct Eliashberg code path that reads ephmat
      ephwrite = .true. → creates ephmat/ dir and writes freq/kqmap/eigenvalues/ephmat
      a2f = .false.   → allows write_freq() to be called (gated by IF(.NOT.a2f))
    """
    p = cfg["compound"]
    d = cfg["dft"]
    w = cfg["wannier"]
    e = cfg["epw"]
    el = cfg["eliashberg"]

    nk = d["coarse_kgrid"]
    nkf = e["fine_kgrid"]
    temps_str = "  ".join(str(t) for t in el["temperatures_K"])
    nstemp = len(el["temperatures_K"])

    bands_skipped_line = ""
    if w.get("bands_skipped"):
        bands_skipped_line = f"  bands_skipped = '{w['bands_skipped']}'\n"

    return f"""\
&inputepw
  prefix        = '{p['prefix']}'
{_amass_block(cfg)}
  outdir        = './tmp/'

  ! Pass 2: read .epmatwp, fine grid interpolation + Eliashberg
  ep_coupling   = .true.
  elph          = .true.
  epbwrite      = .false.
  epbread       = .false.
  epwwrite      = .false.
  epwread       = .true.
  ephwrite      = .true.
  kmaps         = .true.

  wannierize    = .false.
  num_iter      = 500
  nbndsub       = {w['nbndsub']}
{bands_skipped_line}{_projections_block(cfg)}

  dis_win_min   = {w['dis_win_min']}
  dis_win_max   = {w['dis_win_max']}
  dis_froz_min  = {w['dis_froz_min']}
  dis_froz_max  = {w['dis_froz_max']}

  phonselfen    = .true.
  a2f           = .false.
  fsthick       = {e['fsthick']}
  max_memlt     = 16.0

  eliashberg    = .true.
  liso          = .true.
  laniso        = .false.
  limag         = .true.
  lpade         = .true.

  muc           = {el['muc']}
  nsiter        = 500
  conv_thr_iaxis = 1.0d-4
  wscut         = {e['wscut']}

  nstemp        = {nstemp}
  temps         = {temps_str}

  nk1 = {nk[0]}  nk2 = {nk[1]}  nk3 = {nk[2]}
  nq1 = {nk[0]}  nq2 = {nk[1]}  nq3 = {nk[2]}
  mp_mesh_k     = .false.

  nkf1 = {nkf[0]}  nkf2 = {nkf[1]}  nkf3 = {nkf[2]}
  nqf1 = {nkf[0]}  nqf2 = {nkf[1]}  nqf3 = {nkf[2]}

  dvscf_dir     = '../2_phonons/save'
  iverbosity    = 0
/
"""


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _amass_block(cfg):
    """Return amass(1)=..., amass(2)=... lines for EPW."""
    lines = []
    for i, sp in enumerate(cfg["structure"]["species"], 1):
        lines.append(f"  amass({i})      = {sp['mass']:.5f}")
    return "\n".join(lines)
