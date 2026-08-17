# ht-vasp: High-Throughput VASP Workflows for Multi-Principal Element Alloys

## Description

`ht-vasp` is a Python framework for building, running, and managing high-throughput first-principles calculations using VASP (Vienna Ab initio Simulation Package) with a focus on Multi-Principal Element Alloys (MPEAs) and CALPHAD thermodynamic modeling. It wraps atomate2/jobflow makers into reusable `Worker` classes that handle VASP job execution locally or via Slurm, with built-in resume support and result persistence.

The framework supports five main workflow types:
- **Relax**: Two-step structural relaxation (R7 volume + R3 full)
- **Static**: Relaxation followed by static self-consistent calculation
- **NSCF**: Non-self-consistent DOS + band structure calculations
- **QHA**: Quasi-Harmonic Approximation for thermodynamic properties
- **OJ** (OstravaJ): Magnetic exchange coupling J and Curie temperature calculations

It also includes post-processing scripts for CALPHAD thermodynamic fitting (Birch-Murnaghan EOS, SGTE polynomial Gibbs free energy fitting) and electronic structure visualization.

## Directory Structure

```
ht-vasp/                          # Project root
├── pyproject.toml                # Package configuration, dependencies
├── CLAUDE.md                     # Current file — AI assistant guidance
├── .env                          # Environment variables (MP_API_KEY)
├── .gitignore
├── LICENSE                       # MIT License
├── README.md                     # (not present as file, referenced in pyproject.toml)
│
├── htvasp/                       # Main Python package
│   ├── __init__.py
│   ├── workflows/                # Worker implementations (the main abstraction)
│   │   ├── __init__.py           # Exports all workers
│   │   ├── base.py               # Worker base class — core framework
│   │   ├── relax.py              # RelaxWorker — structural relaxation
│   │   ├── static.py             # StaticWorker — static calculation
│   │   ├── nscf.py               # NscfWorker — DOS + band structure
│   │   ├── qha.py                # QhaWorker — quasi-harmonic approximation
│   │   └── oj.py                 # OJWorker — magnetic exchange
│   ├── utils/
│   │   ├── __init__.py
│   │   ├── local.py              # run_locally_custom — custom jobflow local runner with resume
│   │   └── local_old.py          # Older simplified run_locally implementation
│   ├── slurm/
│   │   ├── __init__.py
│   │   └── manager.py            # SlurmConfig + SlurmJobManager — sbatch wrapper
│   ├── oj/                       # OstravaJ magnetic exchange sub-package
│   │   ├── __init__.py
│   │   ├── input_set.py          # OJInputSetGenerator — VASP input for OJ
│   │   ├── jobs.py               # oj_generate, create_flip_jobs, oj_solve
│   │   ├── maker.py              # OJMaker — atomate2-style Maker
│   │   ├── task_doc.py           # OJResult — Pydantic output schema
│   │   ├── scan_supercell.py     # Supercell parameter space scanner
│   │   └── supercell_suggest.py  # Supercell size suggestion heuristic
│   └── model/
│       ├── __init__.py
│       └── endmember.py          # Endmember — fetch/build structures from Materials Project
│
├── ostravaj/                     # OstravaJ magnetic exchange tool (MIT submodule)
│   ├── ostravaj.sh               # Shell wrapper for OstravaJ binary
│   ├── setup.py
│   ├── requirements.txt
│   ├── README.md
│   ├── sync_upstream.sh
│   ├── jmixer/                   # J-mixer python package
│   └── examples/NiO/             # Example input files
│
├── tests/                        # Test scripts (no test runner — use --unit flag)
│   ├── test_relax.py
│   ├── test_static.py
│   ├── test_nscf.py
│   ├── test_qha.py
│   └── test_oj.py
│
├── data/                         # VASP calculation output data
│   ├── poscars/                  # Input POSCAR files for structures
│   ├── endmembers/               # CALPHAD endmember calculations
│   ├── CoNiHO/                   # (Co,Ni)(OH)2 NSCF data
│   ├── CoNiHOS-Co/               # S-doped variants
│   ├── CoNiHOS-Co-H/
│   ├── CoNiHOS-Ni/
│   ├── CoNiHOS-Ni-H/
│   └── modeling-sqsgen/          # SQS generation data
│
├── CMCH_at_CoNi-LDH-S/          # Research project: CMCH catalyst on CoNi-LDH-S
│   ├── CMCH_at_CoNi-LDH-S.md     # Project documentation
│   ├── build_bilayer_ldh.py      # Structure building scripts
│   ├── build_heterostructure.py
│   ├── build_oh_adsorption.py
│   ├── bader_charge_analysis.py
│   ├── comp_hetero.py
│   ├── comp_s_doping.py
│   ├── convert_to_orthorhombic.py
│   ├── conv_test.py
│   ├── gap_scan.py
│   ├── gen_conioh_sqs.py
│   ├── nscf_bulk.py / nscf_hetero.py / nscf_slab.py
│   ├── oh_gas_ref.py / oh_obsorption.py
│   ├── patch_dipole.py
│   └── data/                     # VASP calculation output per system
│
├── docs/                         # Chinese-language documentation
│   ├── 项目架构说明.md            # Project architecture overview
│   ├── 自定义atomate2工作流.md    # Custom atomate2 workflow guide
│   ├── OJ模块使用指南.md           # OJ module user guide
│   ├── OstravaJ论文解读.md        # OstravaJ paper analysis
│   ├── OstravaJ论文原文PDF/       # Original paper
│   ├── htvasp_extension_plan.md  # Extension plan
│   ├── CALPHAD端基数据库.md       # CALPHAD endmember database
│   ├── magmom-guide.md            # MAGMOM guide
│   ├── sqsgenerator参数配置指南.md # sqsgenerator config guide
│   └── AI观点方法评估.md           # AI methods evaluation
│
├── backup-dos/                   # Old DOS workflow scripts
├── backup-old/                   # Old QHA/relax scripts
├── backup-test/                  # Test job scripts (job_cpu.sh, job_gpu.sh, vasprun/)
│
├── output/                       # Analysis JSON output data
├── logs/                         # Slurm job log output
├── tries/                        # Try-out/experimental scripts
├── .claude/                      # Claude Code configuration + skills
│   └── skills/
│       ├── htvasp/               # HT-VASP specific skill
│       └── vasp-skills/          # VASP general skill
│
└── Top-level scripts:
    ├── main_oj.py                # OJ workflow run script (main entry point for magnetic exchange)
    ├── main_qha.py               # QHA workflow run script (main entry point for thermodynamics)
    ├── conioh_nscf.py            # CoNiHO NSCF workflow run script
    ├── nscf_plot.py              # DOS + band structure plotting
    ├── temp_v_e_fit.py           # EOS fitting (Birch-Murnaghan 3rd order)
    ├── temp_gibbs_fit.py         # SGTE Gibbs free energy polynomial fitting
    ├── tmp_gen_psk.py            # SQS structure generation (pymatgen mcsqs wrapper)
    └── check_store.py            # Interactive job store debug/query tool
```

## Build & Run

### Installation

The project uses setuptools with a `pyproject.toml` config:

```bash
# Install the package in development mode
pip install -e .

# Install with dev dependencies
pip install -e ".[dev]"
```

The Python version requirement is >= 3.10.

### Running Unit Tests (no VASP needed)

Each test script has `--unit`, `--local`, and `--slurm` flags:

```bash
python tests/test_relax.py --unit
python tests/test_static.py --unit
python tests/test_nscf.py --unit
python tests/test_qha.py --unit
python tests/test_oj.py --unit
```

There is no formal test runner (pytest is an optional dev dependency). Tests use simple procedural functions with argparse.

### Running Workflows Locally (requires VASP)

```bash
# A single structure locally
python tests/test_relax.py --local
python main_oj.py --oj SER-Fe
python main_oj.py --relax SER-Fe
python main_qha.py -n SER-Al -j static
python main_qha.py -n SER-Al -j qha
python conioh_nscf.py --tick CoNiHO
```

### Submitting to Slurm

```bash
# Submit a single job to Slurm via sbatch --wrap
python tests/test_relax.py --slurm
python main_oj.py --slurm
python main_qha.py -s -j static
python conioh_nscf.py --slurm

# Batch all structures locally
python main_oj.py --batch
python conioh_nscf.py --batch
```

### Environment Setup

The project relies on:
- A conda environment (typically `htvasp` or `.conda/` within the project)
- Environment variable `MP_API_KEY` (for Materials Project API access, stored in `.env`)
- VASP executables: `vasp_std` and `vasp_gam` (typically loaded via `module load vasp-cpu`)
- OstravaJ binary: accessed through `ostravaj/ostravaj.sh`
- Fast scratch directory: `/nfs_ssd/tmp` (configurable per script)

## Configuration

### Worker INCAR Configuration

Workers use a layered INCAR override system:
1. **`global_incar`**: Base defaults set in `Worker.__init__()` (ENCUT=400, ISPIN=2, ISMEAR=1, SIGMA=0.2, etc.)
2. **Per-step INCAR overrides**: e.g., `relax_incar`, `static_incar`, `nscf_dos_incar`, `nscf_band_incar`
3. **`user_incar_settings`**: Passed directly to `VaspInputGenerator` via `user_incar_settings`

Key VASP settings for transition metals / magnetic systems:
- `MAGMOM`: Dict per-element (e.g., `{"Co": 5.0, "Ni": 2.0}`)
- `LDAU` / `LDAUU` / `LDAUJ`: DFT+U parameters — enabled via `LDAU: True, LDAUTYPE: 2`
- `AMIX` / `BMIX` / `AMIX_MAG` / `BMIX_MAG`: Mixing parameters tuned for magnetic convergence
- `IVDW`: van der Waals correction (e.g., 12 = DFT-D3)
- `KPAR` / `NCORE`: Parallelization settings (tune per cluster)

### Slurm Configuration

`SlurmConfig` dataclass in `htvasp/slurm/manager.py`:
- Presets for CPU (`partCPU`) and GPU (`partGPU`) partitions
- Key fields: `ntasks`, `memory`, `partition`, `time_limit`, `cpus_per_task`, `gpus_per_task`
- Environment: `conda_prefix`, `conda_env`, `module_name`
- VASP commands: `vasp_cmd`, `vasp_gam_cmd` (can be custom shell pipelines)

### Job Store

Results are stored in `JSONStore` files (named `store.json`) inside store directories. The `Worker` uses:
- **flow_dir**: Fast storage (e.g., `/nfs_ssd/tmp`) for active computation
- **store_dir**: Persistent storage (e.g., `data/<name>/`) for finished results
- On completion, results are copied from flow_dir to store_dir and flow_dir is cleaned up

## Key Files

### Core Framework

| File | Purpose |
|------|---------|
| `htvasp/workflows/base.py` | `Worker` base class — assembles flows, runs via `run_locally_custom`, manages store lifecycle |
| `htvasp/utils/local.py` | `run_locally_custom()` — custom jobflow local runner with UUID restoration and resume support |
| `htvasp/slurm/manager.py` | `SlurmConfig` dataclass + `SlurmJobManager` for `sbatch --wrap` submission |
| `pyproject.toml` | Package metadata and pinned dependencies |

### Workflow Implementations

| File | Worker | Flow |
|------|--------|------|
| `htvasp/workflows/relax.py` | `RelaxWorker` | R7 relax -> R3 relax (DoubleRelaxMaker) |
| `htvasp/workflows/static.py` | `StaticWorker` | relax -> static |
| `htvasp/workflows/nscf.py` | `NscfWorker` | Static -> NSCF-DOS (uniform) || NSCF-Band (line) |
| `htvasp/workflows/qha.py` | `QhaWorker` | Relax -> EOS -> Phonon -> QHA |
| `htvasp/workflows/oj.py` | `OJWorker` | Generate -> N flip jobs (parallel) -> Solve |

### OstravaJ Sub-package

| File | Purpose |
|------|---------|
| `htvasp/oj/maker.py` | `OJMaker` — atomate2-style Maker assembling the three-step OJ flow |
| `htvasp/oj/jobs.py` | `@job`-decorated functions: `oj_generate`, `create_flip_jobs`, `oj_solve` |
| `htvasp/oj/input_set.py` | `OJInputSetGenerator` — extends `VaspInputGenerator` with OJ.conf generation |
| `htvasp/oj/task_doc.py` | `OJResult` — Pydantic output model with J pairs, Tc, diagnostics |
| `htvasp/oj/scan_supercell.py` | CLI tool to scan supercell parameter space (a, dc) |
| `htvasp/oj/supercell_suggest.py` | Heuristic for suggesting supercell size from structure + dist_cutoff |

### Model / Structure

| File | Purpose |
|------|---------|
| `htvasp/model/endmember.py` | `Endmember` class — fetches structures from Materials Project API and generates alloy supercells (BCC/FCC/HCP/SER templates) |

### Run Scripts

| File | Purpose |
|------|---------|
| `main_oj.py` | Entry point for OJ workflows — supports `--oj`, `--relax`, `--batch`, `--slurm`, `--dry-run` |
| `main_qha.py` | Entry point for QHA/static workflows — supports `-n/-j`, `--slurm`, `--check` |
| `conioh_nscf.py` | NSCF runs for (Co,Ni)(OH)2 and S-doped variants with DFT+U |
| `nscf_plot.py` | DOS + band structure plotting using pymatgen (`DosPlotter`, `BSPlotter`, `BSDOSPlotter`) |
| `temp_v_e_fit.py` | Birch-Murnaghan 3rd-order EOS fitting (single/batch modes) |
| `temp_gibbs_fit.py` | SGTE polynomial Gibbs free energy fitting (single/batch modes) |
| `tmp_gen_psk.py` | SQS structure generation via pymatgen's mcsqs wrapper |
| `check_store.py` | Interactive store query tool for debugging job outputs |

### Research Project: CMCH_at_CoNi-LDH-S

| File | Purpose |
|------|---------|
| `CMCH_at_CoNi-LDH-S/CMCH_at_CoNi-LDH-S.md` | Project documentation (composition, methods, results) |
| `CMCH_at_CoNi-LDH-S/build_bilayer_ldh.py` | Build LDH bilayer structures |
| `CMCH_at_CoNi-LDH-S/build_heterostructure.py` | Build heterostructures |
| `CMCH_at_CoNi-LDH-S/bader_charge_analysis.py` | Bader charge analysis |
| `CMCH_at_CoNi-LDH-S/gap_scan.py` | Band gap scanning workflow |
| `CMCH_at_CoNi-LDH-S/gen_conioh_sqs.py` | SQS generation for CoNiHO |

## Dependencies

### Core Dependencies (pinned versions in pyproject.toml)

| Package | Version | Purpose |
|---------|---------|---------|
| `atomate2[phonons]` | 0.1.4 | VASP workflow framework — makers, set generators, phonon workflow |
| `pymatgen` | 2026.5.4 | Materials structure manipulation (Structure, Lattice, DOS, BandStructure) |
| `pymatgen-core` | 2026.5.17 | Core pymatgen types |
| `emmet-core` | 0.86.4 | Materials data models (task docs) |
| `jobflow` | (atomate2 dep) | Workflow DAG engine (Flow, Job, Maker, JobStore) |
| `maggma` | (atomate2 dep) | Storage backends (JSONStore, MemoryStore) |
| `monty` | 2026.5.18 | Python utilities (serialization, MontyDecoder) |
| `pydantic` | 2.13.4 | Data validation (OJResult model) |
| `phonopy` | 2.48.0 | Phonon calculations (via atomate2 PhononMaker) |
| `seekpath` | 2.2.1 | k-point path for band structure |
| `mp-api` | 0.46.1 | Materials Project API (MPRester) |
| `lxml` | 6.1.1 | XML/vasprun.xml parsing |
| `more-itertools` | 11.1.0 | Iteration utilities |

### Optional/Dev Dependencies

| Package | Purpose |
|---------|---------|
| `pytest >= 6.0` | Test runner (not used directly — tests use argparse-based scripts) |
| `black` | Code formatter |
| `python-dotenv` | .env file loading (for MP_API_KEY) |
| `custodian` | VASP error handling (VaspErrorHandler) |
| `matplotlib` | Plotting (nscf_plot.py, temp scripts) |
| `scipy` | Curve fitting (EOS, Gibbs fitting) |
| `pandas` | Data handling (gibbs fitting) |
| `tqdm` | Progress bars (gibbs fitting) |
| `OstravaJ` | External magnetic exchange binary (via ostravaj.sh) |

### External Tools (not Python packages)

| Tool | Purpose |
|------|---------|
| **VASP** (`vasp_std`, `vasp_gam`) | DFT engine — must be available via module or PATH |
| **OstravaJ** (binary) | Magnetic exchange calculation — accessed via `ostravaj/ostravaj.sh` |
| **mcsqs** (ATAT) | Special Quasirandom Structure generation (optional, via pymatgen wrapper) |
| **Slurm** (`sbatch`, `squeue`, `scancel`, `scontrol`) | Job scheduling on HPC clusters |
| **Materials Project API** | Structure database — requires `MP_API_KEY` in `.env` |

## Conventions

### Code Style and Patterns

- **Worker Pattern**: All workflows inherit from `Worker` (htvasp/workflows/base.py). Each worker defines `_make_flow()` and overrides `get_result()` with the correct output job name regex.
- **Run Script Pattern**: Top-level scripts follow: define `VASP_ARGS`, `GLOBAL_INCAR`, structure names, `run_tick(name, force)` that creates Worker -> calls `run_flow()` -> `get_result()` -> `write_result()`. CLI via argparse with `--tick`, `--batch`, `--slurm`, `--force`/`--rerun` flags.
- **Layer Architecture**: INCAR settings are layered: global defaults in `Worker.__init__` -> `global_incar` passed by caller -> per-step INCAR overrides -> `user_incar_settings`.
- **Resume Mechanism**: The custom `run_locally_custom()` in `htvasp/utils/local.py` restores old UUIDs from the store to make OutputReference links work across runs. Each workflow defaults to `resume=True`.
- **Store Management**: flow results are written to `/nfs_ssd/tmp/<name>-<uuid[:8]>` during computation, then moved to persistent `store_dir` on completion. The store is a `JobStore` backed by `JSONStore`.

### Naming Conventions

- **Workers**: `<Feature>Worker` — PascalCase, e.g., `RelaxWorker`, `NscfWorker`, `QhaWorker`, `OJWorker`
- **Flow Makers**: Standard atomate2 naming — `RelaxMaker`, `StaticMaker`, `NonSCFMaker`, `PhononMaker`
- **Workflow files**: Short, descriptive: `relax.py`, `static.py`, `nscf.py`, `qha.py`, `oj.py`
- **Run scripts**: `main_<workflow>.py` or `<project>_<workflow>.py`, e.g., `main_oj.py`, `conioh_nscf.py`
- **Structure names**: `{BCC|FCC|HCP|SER}-{Elem1}[-{Elem2}]`, e.g., `SER-Fe`, `BCC-Al-Nb`
- **Store directories**: `data/<name>/<!flow!>/store.json` where `<!flow!>` is `ojflow`, `qhaflow`, `staticflow`, etc.
- **Store JSON output files**: `<name>-<workflow>.json`, e.g., `SER-Fe-oj.json`, `SER-Al-qha.json`
- **OJ directory naming**: `{prefix}-{job_name}` or `{prefix}-{job_name}-{counter}`, e.g., `1-init_relax`, `2-eos_relax`
- **INCAR variable naming**: `GLOBAL_INCAR`, `RELAX_INCAR`, `STATIC_INCAR`, `NSCF_INCAR` — all uppercase dicts
- **Test scripts**: `test_<workflow>.py` with `--unit`, `--local`, `--slurm` flags

### Preferred APIs

- **Structure handling**: Use `pymatgen.core.Structure` methods (`.from_file()`, `.to()`, `.replace_species()`, etc.) rather than manual POSCAR parsing.
- **VASP input generation**: Use atomate2 set generators (`RelaxSetGenerator`, `StaticSetGenerator`, `NonSCFSetGenerator`) rather than manually constructing INCAR/KPOINTS/POTCAR.
- **VASP output reading**: Use atomate2/pymatgen task docs. Job outputs from `store.get_output()` are structured dicts with fields like `output["energy"]`, `output["structure"]`, etc.
- **Workflow building**: Compose atomate2 makers (`RelaxMaker`, `StaticMaker`, `DoubleRelaxMaker`, `NonSCFMaker`, `PhononMaker`) rather than writing custom VASP job functions.
- **Custodian error handling**: Use `VaspErrorHandler` in `run_vasp_kwargs["handlers"]` — already configured in `Worker`.
- **Job store queries**: Use `store.query_one()` / `store.query()` with criteria dicts and `get_output()` for structured results.
- **MAGMOM values**: Store as per-element dict (e.g., `MAGMOM: {"Co": 5.0, "Ni": 2.0}`) in INCAR dicts. atomate2 resolves these to per-site arrays.

## Common Tasks

### Adding a New Structure for OJ Calculation

1. Add the POSCAR file to `data/poscars/` (or generate via `Endmember.get_poscar()`)
2. Add structure name to `STRUCTURE_NAMES` list in `main_oj.py`
3. Define `MAGNETIC_SPECIES` and `BASE_SPIN_MAP` entries if new elements
4. Run: `python main_oj.py --relax <name>` then `python main_oj.py --oj <name>`

### Adding a New Structure for QHA Calculation

1. Add POSCAR to `data/poscars/`
2. Add to `STRUCTURE_NAMES` list in `main_qha.py`
3. Run: `python main_qha.py -n <name> -j static` then `python main_qha.py -n <name> -j qha`

### Running NSCF for a Material System

1. Add POSCAR path to `structs` dict in `conioh_nscf.py`
2. Run: `python conioh_nscf.py --tick <name>`

### Creating a New Worker

1. Create `htvasp/workflows/<name>.py`
2. Inherit from `Worker` (htvasp/workflows/base.py)
3. Implement `__init__()` — assemble atomate2 makers, store in `self.flow_maker` or `self.flow_makers`
4. Implement `_make_flow()` — compose makers into a `jobflow.Flow`
5. Override `get_result()` with the correct output job name regex
6. Export in `htvasp/workflows/__init__.py`

### Submitting Batch Jobs to Slurm

```bash
python main_oj.py --slurm             # Submit all OJ structures
python main_qha.py -s -j static       # Submit all static calculations
python main_qha.py -s -j qha          # Submit all QHA calculations
python conioh_nscf.py --slurm         # Submit all NSCF structures
```

### Debugging Store / Job Outputs

```bash
python check_store.py     # Interactive tool to query a store by UUID
```

### Checking Calculation Results

```bash
python main_qha.py -c -j static       # Check all static job energies and magnetizations
python nscf_plot.py                   # Plot DOS and band structure (edit paths inside)
```

### Fitting EOS or Gibbs Free Energy

```bash
python temp_v_e_fit.py --single <path> --phase BCC --elements Al Nb --atom-num 2
python temp_v_e_fit.py --batch <root_dir>
python temp_gibbs_fit.py --single <path> --phase BCC --elements Al Nb --atom-num 2
python temp_gibbs_fit.py --batch <root_dir>
```

### Resuming an Interrupted Calculation

All `run_flow()` calls default to `resume=True`. To force a fresh re-run:
```bash
python main_oj.py --oj SER-Fe --rerun
python main_qha.py -n SER-Al -j qha --rerun
python conioh_nscf.py --tick CoNiHO --force
```

### Testing OJ Configuration (without VASP)

```bash
python main_oj.py --dry-run SER-Fe   # Test config generation only
python main_oj.py --dry-batch         # Test all structures
```
