# main_dos.py

import json
import shutil
import logging
from pathlib import Path
from datetime import datetime

from custodian.vasp.handlers import VaspErrorHandler
from pymatgen.core import Structure

from src.workflow.flow_dos import DosWorker

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


def main():
    vasp_args = {
        "handlers": [VaspErrorHandler()],
        "vasp_cmd": "/bin/bash -c 'module load vasp-cpu && mpirun -np 40 -bootstrap=ssh vasp_std'",
        "vasp_gamma_cmd": "/bin/bash -c 'module load vasp-cpu && mpirun -np 40 -bootstrap=ssh vasp_gam'",
    }
    # vasp_args = {
    #     "handlers": [VaspErrorHandler()],
    #     "vasp_cmd": "/bin/bash -c 'module load vasp-gpu && mpirun -np 1 vasp_gpu_std'",
    #     "vasp_gamma_cmd": "/bin/bash -c 'module load vasp-gpu && mpirun -np 1 vasp_gpu_gam'",
    # }
    potcar_functional = "PBE_64"
    incar_settings = {
        "ENCUT": 520,
        "NELM": 200,
        "LREAL": False,  # for H, C, O
        "SYMPREC": 1e-7,
        "LDAU": True,  # 强烈建议启用 Hubbard U（Co 的 d 电子局域化）
        "LDAUTYPE": 2,  # Dudarev 方法（最常用）
        # 改为使用pymatgen推荐的字典形式设置LDAU参数
        "LDAUL": {"Co": 2},  # Co: l=2 (d轨道)
        "LDAUU": {"Co": 3.5},  # Co 的 U 值（典型范围 3.0–5.0 eV，建议用 3.5 ）
        "LDAUJ": {"Co": 0.0},
        # Parallel
        "NCORE": 2,
        "KPAR": 4,
    }
    worker = DosWorker(vasp_args, potcar_functional, incar_settings)

    name = "Co2H2CO5"
    dosdir = Path(f"data/{name}/dos")
    if dosdir.exists():
        shutil.rmtree(dosdir)
    dosdir.mkdir(parents=True, exist_ok=True)

    json_path = Path(f"data/{name}/dos_result.json")

    # Skip if already done
    if json_path.exists():
        with open(json_path, "r", encoding="utf-8") as jf:
            result = json.load(jf)
        if result["state"] == "successful":
            log.info(f"{name} already done")
            return

    # Run dos
    poscar = Path("data/poscars/POSCAR-Co2H2CO5")
    struct = Structure.from_file(poscar)
    dos_data = worker.run_dos(name, struct, dosdir)
    result = (
        {"name": name, "state": "successful", **dos_data}
        if dos_data
        else {"name": name, "state": "failed"}
    )
    with open(json_path, "w", encoding="utf-8") as jf:
        json.dump(result, jf, indent=2, cls=DateTimeEncoder)
    log.info(f"{name} done")


if __name__ == "__main__":
    main()
