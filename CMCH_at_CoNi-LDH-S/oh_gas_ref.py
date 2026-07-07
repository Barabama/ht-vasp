import json
import time
import logging
import argparse
from pathlib import Path

import numpy as np
from pymatgen.core import Structure, Lattice, Molecule
from htvasp.workflows import StaticWorker
from htvasp.slurm import SlurmJobManager, SlurmConfig

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_gam'",
}

# 气相分子：大盒子 + Gamma-only + 无 U + 无溶剂
GLOBAL_INCAR = {
    "ENCUT": 520,
    "PREC": "Accurate",
    "ALGO": "Fast",
    "NELM": 200,
    "EDIFF": 1e-6,
    "ISPIN": 2,
    "ISMEAR": 0,
    "SIGMA": 0.01,
    "LREAL": False,  # 气相分子用精确 FFT
    "GGA": "PE",
    "IVDW": 0,  # 气相分子不需要 DFT-D3
    "LDAU": False,  # 气相分子不需要 U
    "LWAVE": False,
    "LCHARG": False,
}

RELAX_INCAR = {
    "IBRION": 2,
    "ISIF": 0,  # 只弛豫原子，不弛豫晶格
    "NSW": 50,
    "EDIFFG": -0.02,
}

STATIC_INCAR = {
    "IBRION": -1,
    "NSW": 0,
    "ALGO": "Normal",
    "NELM": 200,
}

root_dir = Path(__file__).parent.resolve()
flow_dir = Path("/nfs_ssd/tmp")
conda_env = str(root_dir.parent / ".conda")
cpu_config = SlurmConfig(
    job_name="cpu_job",
    output_log="cpu_job.log",
    nodes=1,
    ntasks=4,
    memory="20G",
    partition="partCPU",
    conda_env=conda_env,
    module_name="vasp-cpu",
)

molecules = {
    # # H₂O 几何: O-H = 0.957 Å, H-O-H = 104.5°
    # "h2o": Molecule(
    #     ["O", "H", "H"],
    #     [[0, 0, 0], [0, 0.757, 0.587], [0, -0.757, 0.587]],
    # ),
    # # H₂ 几何: H-H = 0.74 Å
    # "h2": Molecule(
    #     ["H", "H"],
    #     [[0, 0, 0], [0, 0, 0.74]],
    # ),
    # OH 自由基: O-H = 0.97 Å, 自旋极化 + 偶极修正
    "oh_gas": Molecule(
        ["O", "H"],
        [[0, 0, 0], [0, 0, 0.97]],
    ),
}

# Per-molecule INCAR overrides（叠加到 GLOBAL_INCAR 之上）
MOLECULE_INCAR = {
    "h2o": {},
    "h2": {},
    "oh_gas": {
        "NELECT": 7,          # O (6) + H (1) = 7, 中性自由基
        "MAGMOM": {"O":1, "H":1},          # 初猜磁矩 2 μB（O 有 2 个未配对电子）
        "LDIPOL": True,       # 偶极修正（OH 有 ~1.7 D 的偶极矩）
        "IDIPOL": 3,          # 沿 z 轴
        "LREAL": "Auto",      # 大盒子 + 高 ENCUT 用 Auto 避免 FEXCF 溢出
    },
}


def build_gas_molecule(mol: Molecule, box_size: float = 15.0) -> Structure:
    """将气相分子放入大盒子中心。"""
    # 将分子放在盒子中心
    mol.translate_sites(range(len(mol)), [box_size / 2, box_size / 2, box_size / 2])
    lattice = Lattice.cubic(box_size)
    return Structure(lattice, mol.species, mol.cart_coords, coords_are_cartesian=True)


def run_gas_calc(name: str, rerun: bool = False):
    store_dir = root_dir / "data" / name
    json_path = store_dir / "static_out.json"

    if not rerun and json_path.exists():
        log.info(f"{name} already done")
        return

    log.info(f"Building {name}...")
    structure = build_gas_molecule(molecules[name])

    try:
        mol_incar = MOLECULE_INCAR.get(name, {})
        worker = StaticWorker(
            vasp_args=VASP_ARGS,
            global_incar={**GLOBAL_INCAR, **mol_incar},
            relax_incar=RELAX_INCAR,
            static_incar=STATIC_INCAR,
        )
        worker.run_flow(
            name=name,
            structure=structure,
            flow_dir=flow_dir,
            store_dir=store_dir,
            resume=not rerun,
        )
        output = worker.get_result("static") or {}
        worker.write_result(output, json_path)

        energy = output.get("energy", "N/A")
        log.info(f"{name}: E = {energy} eV")

    except Exception as e:
        log.error(f"{name} failed: {e}")


def submit_jobs(rerun: bool = False):
    manager = SlurmJobManager()

    for name in molecules:
        job_name = f"{name}_gas"
        job_cmd = f"python {__file__} --static {name} {'--rerun' if rerun else ''}"
        job_log = str(root_dir / "logs" / f"{job_name}.log")
        cpu_config.job_name = job_name
        cpu_config.output_log = job_log
        jid = manager.submit_command(
            command=job_cmd,
            config=cpu_config,
            work_dir=root_dir,
        )
        time.sleep(1)
        if not jid:
            log.error(f"Failed to submit job {job_name}")
            continue
        log.info(f"Submitted job {job_name} with jid {jid}")


def check_jobs():
    for name in molecules:
        path = root_dir / "data" / name / "static_out.json"
        if path.exists():
            with open(path) as f:
                d = json.load(f)
            print(f"{name}: E = {d['output']['energy']:.6f} eV")
        else:
            print(f"{name}: ❌ 未计算")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="OH gas reference")
    parser.add_argument("--slurm", action="store_true", help="Submit to Slurm")
    parser.add_argument("--static", type=str, help="Run static calculation")
    parser.add_argument("--check", action="store_true", help="Check if calculations are done")
    parser.add_argument("--rerun", action="store_true", help="Rerun calculations")
    args = parser.parse_args()
    if args.slurm:
        submit_jobs(args.rerun)
    if args.check:
        check_jobs()
    if args.static:
        run_gas_calc(args.static, rerun=args.rerun)
    else:
        parser.print_help()
