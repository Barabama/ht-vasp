#!/usr/bin/env python3
"""K-point convergence test — 24 independent Slurm jobs across 3 structures × 8 densities."""

import json
import logging
import argparse
from pathlib import Path

from pymatgen.core import Structure

from htvasp.workflows.base import Worker
from htvasp.slurm import SlurmJobManager

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

# ── 0. paths ──────────────────────────────────────────────────────────────────

root_dir = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S")
data_dir = root_dir / "data"

# ── 1. VASP args (CPU and GPU) ────────────────────────────────────────────────

VASP_ARGS_CPU = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
}
VASP_ARGS_GPU = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-gpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-gpu && srun vasp_gam'",
}

# ── 2. global INCAR ───────────────────────────────────────────────────────────

GLOBAL_INCAR = {
    "ENCUT": 520,
    "PREC": "Accurate",
    "ALGO": "Fast",
    "NELM": 200,
    "EDIFF": 1e-6,
    "ISPIN": 2,
    "GGA": "PE",
    "IVDW": 12,
    "LREAL": "Auto",
    "LDAU": True,
    "LDAUTYPE": 2,
    "LDAUPRINT": 1,
    "LASPH": True,
    "LMAXMIX": 4,
    "AMIX": 0.1,
    "BMIX": 1e-4,
    "AMIX_MAG": 0.4,
    "BMIX_MAG": 1e-4,
    "KPAR": 1,
    "NCORE": 2,
}

# ── 3. static INCAR overrides ─────────────────────────────────────────────────

STATIC_INCAR = {
    "ISMEAR": -5,
    "NSW": 0,
    "IBRION": -1,
    "NELM": 200,
}

# ── 4. densities ──────────────────────────────────────────────────────────────

DENSITIES = [32, 50, 64, 80, 100, 120, 150, 200]

GPU_ATOM_THRESHOLD = 36

# ── 5. structures ─────────────────────────────────────────────────────────────

structs = {
    "CMCH-bulk": {
        "poscar": str(data_dir / "CoMnH2CO5" / "3-static" / "CONTCAR.gz"),
        "incar": {
            "MAGMOM": {"Co": 5.0, "Mn": 5.0, "C": 0.6, "H": 0.6, "O": 0.6},
            "LDAUL": {"Co": 2, "Mn": 2, "C": -1, "H": -1, "O": -1},
            "LDAUU": {"Co": 3.32, "Mn": 5.00},
            "LDAUJ": {"Co": 0.0, "Mn": 0.0},
        },
    },
    "hetero_intrinsic": {
        "poscar": str(data_dir / "hetero_intrinsic" / "3-static" / "CONTCAR.gz"),
        "incar": {
            "MAGMOM": {"Co": 5.0, "Mn": 5.0, "Ni": 2.0, "C": 0.6, "H": 0.6, "O": 0.6},
            "LDAUL": {"Co": 2, "Mn": 2, "Ni": 2, "C": -1, "H": -1, "O": -1},
            "LDAUU": {"Co": 3.32, "Mn": 5.00, "Ni": 6.20},
            "LDAUJ": {"Co": 0.0, "Mn": 0.0, "Ni": 0.0},
        },
    },
    "LDH-slab": {
        "poscar": str(data_dir / "CoNiOH2-slab" / "3-static" / "CONTCAR.gz"),
        "incar": {
            "MAGMOM": {"Co": 5.0, "Ni": 2.0, "H": 0.6, "O": 0.6},
            "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1},
            "LDAUU": {"Co": 3.32, "Ni": 6.20},
            "LDAUJ": {"Co": 0.0, "Ni": 0.0},
        },
    },
}

# ── 6. single-shot StaticWorker ───────────────────────────────────────────────

class SingleStaticWorker(Worker):
    """Run one static calculation at a given reciprocal_density."""

    def __init__(self, vasp_args, global_incar, static_incar, reciprocal_density):
        super().__init__("conv_test", vasp_args, global_incar)
        from atomate2.vasp.jobs.core import StaticMaker
        from atomate2.vasp.sets.core import StaticSetGenerator

        self.flow_maker = StaticMaker(
            name="static",
            input_set_generator=StaticSetGenerator(
                user_kpoints_settings={"reciprocal_density": reciprocal_density},
                user_incar_settings={
                    **self.global_incar,
                    "ISTART": 0,
                    "ICHARG": 2,
                    "NSW": 0,
                    "IBRION": -1,
                    "NELM": 200,
                    "LWAVE": False,
                    "LCHARG": False,
                    **static_incar,
                },
            ),
            run_vasp_kwargs=self.run_vasp_kwargs,
            stop_children_kwargs={"handle_unsuccessful": False},
        )

    def _make_flow(self, structure, prev_dir=None):
        return self.flow_maker.make(structure, prev_dir)


# ── 7. run_tick ───────────────────────────────────────────────────────────────

def run_tick(name: str, density: int, device: str = "cpu", rerun: bool = False):
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = data_dir / f"{name}_conv_k{density}"
    json_path = store_dir / "energy.json"

    if not rerun and json_path.exists():
        log.info("%s k%d already done", name, density)
        return

    poscar = structs[name]["poscar"]
    if not Path(poscar).exists():
        log.warning("CONTCAR missing, skip: %s", poscar)
        return

    structure = Structure.from_file(poscar)
    natoms = len(structure)
    log.info("%s k%d: %d atoms, device=%s", name, density, natoms, device)

    worker = SingleStaticWorker(
        vasp_args=VASP_ARGS_GPU if device == "gpu" else VASP_ARGS_CPU,
        global_incar={**GLOBAL_INCAR, **structs[name]["incar"]},
        static_incar=STATIC_INCAR,
        reciprocal_density=density,
    )
    worker.run_flow(
        name=f"{name}_k{density}",
        structure=structure,
        flow_dir=flow_dir,
        store_dir=store_dir,
        resume=not rerun,
    )

    output = worker.get_result("static")
    energy = output.get("output", {}).get("energy") if output else None
    if energy is not None:
        worker.write_result(
            {"energy": energy, "natoms": natoms, "density": density},
            json_path,
        )
        log.info("%s k%d done: energy=%.6f eV", name, density, energy)
    else:
        log.error("%s k%d: no energy in output", name, density)


# ── 8. submit jobs (24 total) ─────────────────────────────────────────────────

def submit_jobs(rerun: bool = False):
    manager = SlurmJobManager()
    conda_env = str(root_dir.parent / ".conda")

    for name, cfg in structs.items():
        poscar = cfg["poscar"]
        if not Path(poscar).exists():
            log.warning("Skipping %s (no CONTCAR)", name)
            continue

        structure = Structure.from_file(poscar)
        natoms = len(structure)
        use_gpu = natoms >= GPU_ATOM_THRESHOLD

        for d in DENSITIES:
            job_name = f"{name}_conv_k{d}"
            command = (
                f"python {__file__} --tick {name} --density {d}"
                f" --device {'gpu' if use_gpu else 'cpu'}"
                f" {'--rerun' if rerun else ''}"
            )
            output_log = str(root_dir / "logs" / f"{job_name}.log")

            if use_gpu:
                config = manager.get_gpu_config(
                    job_name=job_name,
                    output_log=output_log,
                    nodes=1,
                    ntasks=1,
                    memory="16G",
                    conda_env=conda_env,
                    module_name="vasp-gpu",
                )
            else:
                config = manager.get_cpu_config(
                    job_name=job_name,
                    output_log=output_log,
                    nodes=1,
                    ntasks=32,
                    memory="64G",
                    conda_env=conda_env,
                    module_name="vasp-cpu",
                )

            jid = manager.submit_command(command=command, config=config, workdir=str(root_dir))
            tag = "GPU" if use_gpu else "CPU"
            log.info("[%s] %s (%d atoms) → %s", tag, job_name, natoms, jid or "FAILED")


# ── 9. summary ────────────────────────────────────────────────────────────────

def print_summary():
    print(f"\n{'Structure':<22s} {'Natoms':>6s} {'Device':>6s}", end="")
    for d in DENSITIES:
        print(f"  k{d:>3d}", end="")
    print(f"\n{'-' * (35 + 8 * len(DENSITIES))}")

    for name in structs:
        energies = {}
        natoms = "?"
        device = "?"

        for d in DENSITIES:
            jp = data_dir / f"{name}_conv_k{d}" / "energy.json"
            if jp.exists():
                r = json.loads(jp.read_text())
                energies[d] = r["energy"]
                natoms = r.get("natoms", natoms)

        if not energies:
            print(f"{name:<22s}  {'—':>6s}  {'—':>6s}")
            continue

        device = "GPU" if int(natoms) >= GPU_ATOM_THRESHOLD else "CPU"
        ref_e = energies[max(energies)] / int(natoms)

        print(f"{name:<22s} {str(natoms):>6s} {device:>6s}", end="")
        for d in DENSITIES:
            if d in energies:
                de = (energies[d] / int(natoms) - ref_e) * 1000
                marker = "*" if abs(de) < 1.0 else " "
                print(f" {marker}{de:>5.1f}", end="")
            else:
                print(f" {'—':>6s}", end="")
        print()

    print(f"\n  unit: meV/atom vs k{DENSITIES[-1]} reference  |  * = converged < 1 meV/atom\n")


# ── 10. CLI ───────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="K-point convergence test")
    parser.add_argument("--tick", type=str, help="Structure name")
    parser.add_argument("--density", type=int, help="Reciprocal density")
    parser.add_argument("--device", type=str, default="cpu", choices=["cpu", "gpu"])
    parser.add_argument("--slurm", action="store_true", help="Submit all 24 jobs")
    parser.add_argument("--rerun", action="store_true")
    parser.add_argument("--summary", action="store_true", help="Print convergence table")
    args = parser.parse_args()

    if args.slurm:
        submit_jobs(rerun=args.rerun)
    elif args.tick and args.density:
        run_tick(args.tick, args.density, device=args.device, rerun=args.rerun)
    elif args.summary:
        print_summary()
    else:
        parser.print_help()
