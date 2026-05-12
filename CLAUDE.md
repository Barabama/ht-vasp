# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ht-vasp is a high-throughput VASP workflow framework for Multi-Principal Element Alloys (MPEA). It wraps atomate2/jobflow makers into reusable `Worker` classes that handle VASP job execution locally or via Slurm, with built-in resume support and result persistence.

## Commands

```bash
# Run unit tests (all test scripts, no VASP needed)
python tests/test_relax.py --unit
python tests/test_static.py --unit
python tests/test_nscf.py --unit

# Run a workflow locally (requires VASP in PATH)
python tests/test_relax.py --local

# Submit a workflow to Slurm
python tests/test_relax.py --slurm
```

No formal test runner is configured. Test scripts use argparse with `--unit`/`--local`/`--slurm` flags. There is no lint/format command configured beyond the optional `black` dev dependency.

## Architecture

### Worker Pattern

All workflows inherit from `Worker` ([htvasp/workflows/base.py](htvasp/workflows/base.py)). Each Worker:

1. **`__init__`**: Assembles atomate2 makers (e.g., `RelaxMaker`, `StaticMaker`, `NonSCFMaker`) wrapping them with `DoubleRelaxMaker` where needed. Accepts `global_incar` plus per-step INCAR overrides (e.g., `relax_incar`, `static_incar`).
2. **`_make_flow(structure, prev_dir)`**: Composes makers into a `jobflow.Flow`, connecting jobs via their `.output` references.
3. **`run_flow(name, structure, ...)`**: Executes the flow via `run_locally_custom()`, running in a temp `flow_dir` (fast storage like `/tmp`) and moving results to `store_dir` (persistent) on completion. Supports `resume` from previous incomplete runs.
4. **`get_result(job_name)`**: Queries the `JobStore` for a job's output by name regex.
5. **`write_result(data, json_path)`**: Serializes result dict to JSON.

### Available Workers

| Worker | Flow | Purpose |
|--------|------|---------|
| `RelaxWorker` | R7→R3 double relax | Volume + full structural optimization |
| `StaticWorker` | Relax → Static | Final energy with charge density |
| `NscfWorker` | Relax → Static → DOS ∥ Band | Electronic structure (DOS + band) |
| `QhaWorker` | Relax → EOS → Phonon → QHA | Thermodynamic properties via QHA |
| `OJWorker` | Generate → N×Flip (parallel) → Solve | Magnetic exchange coupling J & Tc |

### Slurm Integration

`SlurmJobManager` ([htvasp/slurm/manager.py](htvasp/slurm/manager.py)) submits jobs via `sbatch --wrap` without script files. `SlurmConfig` is a dataclass with CPU/GPU partition presets. The manager supports conda env activation and module loading within the job.

### OJ (OstravaJ) Sub-package

`htvasp/oj/` contains the OstravaJ magnetic exchange workflow:
- `input_set.py`: `OJInputSetGenerator` extends `VaspInputGenerator`, generates `OJ.conf` for the OstravaJ binary
- `jobs.py`: `@job`-decorated functions (`oj_generate`, `create_flip_jobs`, `oj_solve`) that shell out to `ostravaj/ostravaj.sh`
- `maker.py`: `OJMaker` assembles the three steps into a `Flow`

### Resume Mechanism

`run_locally_custom()` ([htvasp/utils/local.py](htvasp/utils/local.py)) extends jobflow's local runner with:
- Sequential directory naming: `{prefix}-{job_name}` (e.g., `1-init_relax`, `2-eos_relax`)
- UUID restoration: matches jobs by name+index to previous store entries, restores old UUIDs so `OutputReference` links work across runs
- Per-job completion check via `_validate_job_output()`

## Key Libraries and Conventions

**Always prefer using atomate2/pymatgen APIs** for VASP data extraction and processing rather than parsing raw files:

- **Structure handling**: Use `pymatgen.core.Structure` and its methods (`.from_file()`, `.to()`, `.replace_species()`, etc.). Never parse POSCAR/vasprun.xml manually.
- **VASP input generation**: Use atomate2 set generators (`RelaxSetGenerator`, `StaticSetGenerator`, `NonSCFSetGenerator`) rather than manually constructing INCAR/POTCAR/KPOINTS.
- **VASP output reading**: Use atomate2/pymatgen task docs and analysis utilities. Job outputs (accessible via `store.get_output()`) are task doc dicts with structured fields like `output["energy"]`, `output["structure"]`, etc.
- **Workflow building**: Use atomate2 makers (`RelaxMaker`, `StaticMaker`, `DoubleRelaxMaker`, `NonSCFMaker`, `PhononMaker`) from `atomate2.vasp.jobs.core` and `atomate2.vasp.flows.*`. Compose them rather than writing custom VASP job functions.
- **Custodian error handling**: Use `custodian.vasp.handlers` (e.g., `VaspErrorHandler`) in `run_vasp_kwargs["handlers"]` — already configured in the `Worker` base class.
- **Job store queries**: Use `JobStore.query()` / `query_one()` with criteria dicts and `get_output()` to retrieve structured results. The store is a `JSONStore` backed by `store.json`.

## Run Script Pattern

Top-level run scripts (e.g., `conioh_nscf.py`) follow a consistent pattern:
1. Define `VASP_ARGS` dict with `vasp_cmd` and `vasp_gamma_cmd`
2. Define `GLOBAL_INCAR` dict and per-step overrides (`RELAX_INCAR`, `STATIC_INCAR`, etc.)
3. Create structs dict mapping names to POSCAR file paths
4. Define a `run_tick(name, force)` function that creates the Worker, calls `run_flow()`, then `get_result()` + `write_result()`
5. CLI via argparse: `--tick` for single, `--batch` for all local, `--slurm` for batch submission

## Important Files

- `data/poscars/` — input POSCAR files for production runs
- `logs/` — Slurm job log output directory
- `ostravaj/ostravaj.sh` — OstravaJ binary wrapper for magnetic exchange calculations
- Store output is written to project-local directories (not `/tmp`), typically under `data/<name>/store.json`
