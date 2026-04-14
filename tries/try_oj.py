#!/usr/bin/env python3
"""
try_oj.py - Validate OstravaJ workflow for binary metal alloys

This script validates the OstravaJ workflow for calculating magnetic exchange
interactions in binary metal alloys without relying on atomate2/pymatgen.
Uses vaspkit for VASP input generation where possible.

Stage 2.1: OJ配置dict/JSON化支持
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Literal, Optional, Union

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
log = logging.getLogger(__name__)


@dataclass
class OJConfig:
    """OJ配置类，支持dict/JSON序列化"""
    j_count: int | None = None
    dist_cutoff: float | None = None
    magnetic_ion_types: list[str] | None = None
    noncollinear: bool = False
    base_spin: float | list[float] = 1.0
    extend_poscar: tuple[int, int, int] = (1, 1, 1)

    def to_dict(self) -> dict:
        """转换为字典"""
        d = asdict(self)
        return d

    def to_json(self, fp: str | Path) -> None:
        """保存为JSON文件"""
        with open(fp, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_dict(cls, d: dict) -> "OJConfig":
        """从字典创建"""
        if "extend_poscar" in d and isinstance(d["extend_poscar"], list):
            d["extend_poscar"] = tuple(d["extend_poscar"])
        return cls(**d)

    @classmethod
    def from_json(cls, fp: str | Path) -> "OJConfig":
        """从JSON文件加载"""
        with open(fp) as f:
            d = json.load(f)
        return cls.from_dict(d)

    def validate(self) -> None:
        """验证配置有效性"""
        if self.j_count is None and self.dist_cutoff is None:
            raise ValueError("Either j_count or dist_cutoff must be specified")
        if self.j_count is not None and self.dist_cutoff is not None:
            raise ValueError("Cannot specify both j_count and dist_cutoff")

OJ_CONFIGS = {
    "J_count": "J_count",
    "dist_cutoff": "dist_cutoff",
    "magnetic_ion_types": "magnetic_ion_types",
    "noncollinear": "noncollinear",
    "base_spin": "base_spin",
    "extend_poscar": "extend_poscar",
}


class OstravaJWorkflow:
    def __init__(
        self,
        workdir: Path | str,
        system_name: str = "test",
        potcar_path: str = "",
    ):
        self.workdir = Path(workdir)
        self.system_name = system_name
        self.potcar_path = potcar_path
        self.input_dir = self.workdir / "input"
        self.run_dir = self.workdir / "run_oj"

    def prepare_input_from_config(
        self,
        oj_config: OJConfig,
        structure_type: Literal["bcc", "fcc", "hcp"] = "bcc",
        elem1: str = "Fe",
        elem2: str = "Co",
        encut: int = 400,
    ) -> OJConfig:
        """从OJConfig对象准备输入文件

        Args:
            oj_config: OJ配置对象
            structure_type: 晶体结构类型
            elem1: 第一元素
            elem2: 第二元素
            encut: 截断能

        Returns:
            验证后的OJConfig对象
        """
        oj_config.validate()
        log.info(f"Preparing input for {self.system_name} with OJConfig")
        self.input_dir.mkdir(parents=True, exist_ok=True)

        if structure_type == "bcc":
            poscar_content = self.create_poscar_bcc(elem1, elem2)
        elif structure_type == "fcc":
            poscar_content = self.create_poscar_fcc(elem1, elem2)
        elif structure_type == "hcp":
            poscar_content = self.create_poscar_hcp(elem1, elem2)
        else:
            raise ValueError(f"Unknown structure type: {structure_type}")

        with open(self.input_dir / "POSCAR", "w") as f:
            f.write(poscar_content)

        with open(self.input_dir / "INCAR", "w") as f:
            f.write(self.create_incar(encut=encut))

        if self.potcar_path and os.path.exists(self.potcar_path):
            log.info(f"Copying POTCAR from {self.potcar_path}")
            shutil.copy(self.potcar_path, self.input_dir / "POTCAR")
        else:
            log.warning("POTCAR not found. Please provide valid POTCAR path.")

        oj_conf_str = self.create_oj_conf_from_config(oj_config)
        with open(self.input_dir / "OJ.conf", "w") as f:
            f.write(oj_conf_str)

        log.info(f"Input files prepared in {self.input_dir}")
        log.info(f"OJConfig: {oj_config.to_dict()}")
        return oj_config

    def create_oj_conf_from_config(self, oj_config: OJConfig) -> str:
        """从OJConfig对象生成OJ.conf字符串"""
        lines = []

        if oj_config.j_count is not None:
            lines.append(f"J_count {oj_config.j_count}")
        elif oj_config.dist_cutoff is not None:
            lines.append(f"dist_cutoff {oj_config.dist_cutoff}")

        if oj_config.magnetic_ion_types:
            lines.append(f"magnetic_ion_types {' '.join(oj_config.magnetic_ion_types)}")

        lines.append("noncollinear 1" if oj_config.noncollinear else "noncollinear 0")

        if isinstance(oj_config.base_spin, list) and len(oj_config.base_spin) == 3:
            lines.append(f"base_spin {' '.join(str(x) for x in oj_config.base_spin)}")
        else:
            lines.append(f"base_spin {oj_config.base_spin}")

        lines.append(f"extend_poscar {oj_config.extend_poscar[0]} {oj_config.extend_poscar[1]} {oj_config.extend_poscar[2]}")

        return "\n".join(lines) + "\n"

    def prepare_input_from_structure(
        self,
        structure: "Structure",
        oj_config: OJConfig,
        encut: int = 400,
    ) -> OJConfig:
        """从pymatgen Structure对象准备输入文件

        Args:
            structure: pymatgen Structure对象
            oj_config: OJ配置对象
            encut: 截断能

        Returns:
            验证后的OJConfig对象
        """
        try:
            from pymatgen.core import Structure
        except ImportError:
            log.error("pymatgen is required for prepare_input_from_structure")
            raise

        oj_config.validate()
        log.info(f"Preparing input from pymatgen Structure: {self.system_name}")
        self.input_dir.mkdir(parents=True, exist_ok=True)

        structure.to(filename=str(self.input_dir / "POSCAR"), fmt="poscar")
        log.info(f"POSCAR written from Structure")

        with open(self.input_dir / "INCAR", "w") as f:
            f.write(self.create_incar(encut=encut))

        if self.potcar_path and os.path.exists(self.potcar_path):
            log.info(f"Copying POTCAR from {self.potcar_path}")
            shutil.copy(self.potcar_path, self.input_dir / "POTCAR")
        else:
            log.warning("POTCAR not found. Please provide valid POTCAR path.")

        oj_conf_str = self.create_oj_conf_from_config(oj_config)
        with open(self.input_dir / "OJ.conf", "w") as f:
            f.write(oj_conf_str)

        log.info(f"Input files prepared in {self.input_dir}")
        log.info(f"Structure: {len(structure)} sites, {structure.composition.reduced_formula}")
        log.info(f"OJConfig: {oj_config.to_dict()}")
        return oj_config

    def create_poscar_bcc(self, elem1: str, elem2: str) -> str:
        poscar = f"""B2 {elem1}-{elem2} structure
1.0
  2.850000 0.000000 0.000000
  0.000000 2.850000 0.000000
  0.000000 0.000000 2.850000
  {elem1} {elem2}
  1 1
Direct
  0.000000 0.000000 0.000000
  0.500000 0.500000 0.500000
"""
        return poscar

    def create_poscar_fcc(self, elem1: str, elem2: str) -> str:
        poscar = f"""L12 {elem1}-{elem2} structure
1.0
  3.520000 0.000000 0.000000
  0.000000 3.520000 0.000000
  0.000000 0.000000 3.520000
  {elem1} {elem2}
  1 3
Direct
  0.000000 0.000000 0.000000
  0.000000 0.500000 0.500000
  0.500000 0.000000 0.500000
  0.500000 0.500000 0.000000
"""
        return poscar

    def create_poscar_hcp(self, elem1: str, elem2: str) -> str:
        poscar = f"""{elem1}-{elem2} HCP structure
1.0
  2.500000 0.000000 0.000000
 -1.250000 2.165000 0.000000
  0.000000 0.000000 4.100000
  {elem1} {elem2}
  2 2
Direct
  0.000000 0.000000 0.000000
  0.000000 0.000000 0.500000
  0.333333 0.666667 0.250000
  0.333333 0.666667 0.750000
"""
        return poscar

    def create_incar(self, encut: int = 400, ispin: int = 2) -> str:
        incar = f"""System = {self.system_name}
ENCUT = {encut}
ISPIN = {ispin}
ISMEAR = 1
SIGMA = 0.2
ALGO = Fast
NELM = 100
EDIFF = 1E-5
EDIFFG = -0.05
IBRION = 2
ISIF = 3
NSW = 50
LREAL = Auto
PREC = Normal
LWAVE = False
LCHARG = False
KPAR = 2
NCORE = 4
"""
        return incar

    def create_oj_conf(
        self,
        j_count: int | None = None,
        dist_cutoff: float | None = None,
        magnetic_ion_types: list[str] | None = None,
        noncollinear: bool = False,
        base_spin: float | list[float] = 1.0,
        extend_poscar: tuple[int, int, int] = (1, 1, 1),
    ) -> str:
        lines = []

        if j_count is not None:
            lines.append(f"J_count {j_count}")
        elif dist_cutoff is not None:
            lines.append(f"dist_cutoff {dist_cutoff}")
        else:
            raise ValueError("Either j_count or dist_cutoff must be specified")

        if magnetic_ion_types:
            lines.append(f"magnetic_ion_types {' '.join(magnetic_ion_types)}")

        if noncollinear:
            lines.append("noncollinear 1")
        else:
            lines.append("noncollinear 0")

        if isinstance(base_spin, list) and len(base_spin) == 3:
            lines.append(f"base_spin {' '.join(str(x) for x in base_spin)}")
        else:
            lines.append(f"base_spin {base_spin}")

        lines.append(f"extend_poscar {extend_poscar[0]} {extend_poscar[1]} {extend_poscar[2]}")

        return "\n".join(lines) + "\n"

    def prepare_input(
        self,
        structure_type: Literal["bcc", "fcc", "hcp"] = "bcc",
        elem1: str = "Fe",
        elem2: str = "Co",
        encut: int = 400,
        j_count: int | None = None,
        dist_cutoff: float | None = None,
        magnetic_ion_types: list[str] | None = None,
        noncollinear: bool = False,
        base_spin: float | list[float] = 1.0,
        extend_poscar: tuple[int, int, int] = (2, 2, 2),
    ):
        log.info(f"Preparing input for {self.system_name}")
        self.input_dir.mkdir(parents=True, exist_ok=True)

        if structure_type == "bcc":
            poscar_content = self.create_poscar_bcc(elem1, elem2)
        elif structure_type == "fcc":
            poscar_content = self.create_poscar_fcc(elem1, elem2)
        elif structure_type == "hcp":
            poscar_content = self.create_poscar_hcp(elem1, elem2)
        else:
            raise ValueError(f"Unknown structure type: {structure_type}")

        with open(self.input_dir / "POSCAR", "w") as f:
            f.write(poscar_content)

        with open(self.input_dir / "INCAR", "w") as f:
            f.write(self.create_incar(encut=encut))

        if self.potcar_path and os.path.exists(self.potcar_path):
            log.info(f"Copying POTCAR from {self.potcar_path}")
            shutil.copy(self.potcar_path, self.input_dir / "POTCAR")
        else:
            log.warning("POTCAR not found. Please provide valid POTCAR path.")

        oj_conf = self.create_oj_conf(
            j_count=j_count,
            dist_cutoff=dist_cutoff,
            magnetic_ion_types=magnetic_ion_types or [elem1, elem2],
            noncollinear=noncollinear,
            base_spin=base_spin,
            extend_poscar=extend_poscar,
        )
        with open(self.input_dir / "OJ.conf", "w") as f:
            f.write(oj_conf)

        log.info(f"Input files prepared in {self.input_dir}")

    def generate_configs(self, dry_run: bool = False, all_configs: bool = False) -> bool:
        log.info("Running ostravaj generate...")
        ostravaj_sh = Path(__file__).parent / "ostravaj" / "ostravaj.sh"
        cmd = [str(ostravaj_sh), "generate", "-i", str(self.input_dir)]
        if dry_run:
            cmd.append("--dry-run")
        if all_configs:
            cmd.append("--all")
        if not dry_run:
            cmd.extend(["-r", str(self.run_dir)])

        result = subprocess.run(cmd, capture_output=True, text=True)
        log.info(f"Generate stdout:\n{result.stdout}")
        if result.stderr:
            log.warning(f"Generate stderr:\n{result.stderr}")

        if result.returncode != 0:
            log.error(f"Generate failed with return code {result.returncode}")
            return False

        if dry_run:
            if "Success" in result.stdout or "suitable magnetic configurations" in result.stdout:
                log.info("Generate (dry-run) completed successfully")
                return True
            else:
                log.error("Generate (dry-run) may have failed")
                return False

        flip_dirs = list(self.run_dir.glob("flip*"))
        log.info(f"Generated {len(flip_dirs)} magnetic configurations")
        return len(flip_dirs) > 0

    def submit_vasp_jobs(
        self,
        slurm_config: Optional["SlurmConfig"] = None,
        vasp_cmd: str = "srun vasp_std",
        conda_env: str = "",
        module_name: str = "",
    ) -> list[str]:
        """提交所有flip*目录的VASP任务到Slurm

        Args:
            slurm_config: Slurm配置，如果为None则使用默认CPU配置
            vasp_cmd: VASP命令
            conda_env: Conda环境名称
            module_name: 模块名称

        Returns:
            提交的job_id列表
        """
        from htvasp.slurm import SlurmJobManager, SlurmConfig

        manager = SlurmJobManager()
        if slurm_config is None:
            slurm_config = manager.get_cpu_config(
                job_name=f"oj-{self.system_name}",
                ntasks=48,
                memory="20G",
            )

        flip_dirs = sorted(self.run_dir.glob("flip*"))
        if not flip_dirs:
            log.warning("No flip* directories found to submit")
            return []

        job_ids = []
        for flip_dir in flip_dirs:
            job_id = manager.submit_command(
                command=vasp_cmd,
                config=slurm_config,
                conda_env=conda_env,
                module_name=module_name,
                workdir=flip_dir,
            )
            if job_id:
                job_ids.append(job_id)
                log.info(f"Submitted {flip_dir.name}: job_id={job_id}")
            else:
                log.error(f"Failed to submit {flip_dir.name}")

        log.info(f"Total submitted: {len(job_ids)} jobs")
        return job_ids

    def solve(self) -> dict | None:
        log.info("Running ostravaj solve...")
        ostravaj_sh = Path(__file__).parent / "ostravaj" / "ostravaj.sh"
        cmd = [str(ostravaj_sh), "solve", "-r", str(self.run_dir)]

        result = subprocess.run(cmd, capture_output=True, text=True)
        log.info(f"Solve stdout:\n{result.stdout}")
        if result.stderr:
            log.warning(f"Solve stderr:\n{result.stderr}")

        solution_file = self.run_dir / "OJ_solution.json"
        if solution_file.exists():
            with open(solution_file) as f:
                solution = json.load(f)
            log.info(f"Solution loaded from {solution_file}")
            return solution
        else:
            log.error("Solution file not found")
            return None

    def solve_with_mock(self) -> dict | None:
        log.info("Running solve with mock data...")
        import sys
        ostravaj_path = str(Path(__file__).parent / "ostravaj")
        if ostravaj_path not in sys.path:
            sys.path.insert(0, ostravaj_path)
        import jmixer.config
        import jmixer.solver
        import jmixer.vasprun

        input_config = jmixer.config.InputConfig(str(self.run_dir / "input"))
        s = jmixer.solver.ExchangeSystemSolver(input_config)

        flip_dirs = sorted(self.run_dir.glob("flip*"))
        for flip_dir in flip_dirs:
            poscar = jmixer.poscar.Poscar(str(flip_dir / "POSCAR"))
            with open(flip_dir / "INCAR") as f:
                for line in f:
                    if line.startswith("MAGMOM"):
                        spins = [float(x) for x in line.split("=")[1].split()]
                        break

            class MockVaspResult:
                def __init__(self, fname, moments_list, total_energy=-400.0):
                    self.fname = fname
                    self._moments = moments_list
                    self.total_energy = total_energy
                def getMoments(self):
                    return self._moments

            mock_result = MockVaspResult(str(flip_dir / "vasprun.xml"), spins)
            s.addVaspRun(mock_result, poscar)

        sol = s.solve()
        log.info(f"Solution:\n{sol}")

        if sol:
            solution_data = {
                "J_reprs": sol.J_reprs,
                "Js": sol.Js.tolist() if hasattr(sol.Js, 'tolist') else list(sol.Js),
                "Tc_MFA": sol.Tc_MFA,
                "Tc_RPA": sol.Tc_RPA,
                "resid": sol.resid,
                "is_complete": sol.is_complete,
            }
            solution_file = self.run_dir / "OJ_solution.json"
            with open(solution_file, "w") as f:
                json.dump(solution_data, f, indent=2)
            return solution_data
        return None


def run_mock_vasp(flip_dir: Path):
    """Mock VASP calculation for testing purposes."""
    poscar_file = flip_dir / "POSCAR"
    incar_file = flip_dir / "INCAR"

    if not poscar_file.exists() or not incar_file.exists():
        return False

    with open(incar_file) as f:
        incar_content = f.read()

    magmom_line = ""
    for line in incar_content.split("\n"):
        if line.strip().startswith("MAGMOM"):
            magmom_line = line
            break

    mock_energy = -400.0 + len(magmom_line) * 0.01

    vasprun_xml = f"""<?xml version="1.0" encoding="iso-8859-1"?>
<modeling>
 <generator name="vasp" version="5.4.4">
  <i name="date">2024_01_01</i>
  <i name="time">12:00:00</i>
 </generator>
 <parameters>
  <separator name="electronic">
   <i type="int" name="NELECT">0</i>
  </separator>
 </parameters>
 <atominfo>
  <atoms type="int">{len(magmom_line.split()) - 2}</atoms>
  <types type="int">1</types>
 </atominfo>
 <structure>
  <crystal>
   <varray name="basis">
    <v>5.0 0.0 0.0</v>
    <v>0.0 5.0 0.0</v>
    <v>0.0 0.0 5.0</v>
   </varray>
  </crystal>
  <varray name="positions">
   <v>0.0 0.0 0.0</v>
   <v>0.5 0.5 0.5</v>
  </varray>
 </structure>
 <calculation>
  <scstep>
   <e name="efrree">{mock_energy}</e>
  </scstep>
  <energy>
   <i name="e_fren">0.0</i>
   <i name="e_0">0.0</i>
  </energy>
  <varray name="meandrift">
   <v>0.0 0.0 0.0</v>
  </varray>
  <varray name="totalmag">
   <v>{len(magmom_line.split()[2:])}</v>
  </varray>
 </calculation>
</modeling>
"""
    with open(flip_dir / "vasprun.xml", "w") as f:
        f.write(vasprun_xml)

    return True


def main():
    parser = argparse.ArgumentParser(description="Validate OstravaJ workflow")
    parser.add_argument("--workdir", type=str, default="./oj_test",
                        help="Working directory")
    parser.add_argument("--system", type=str, default="FeCo",
                        help="System name")
    parser.add_argument("--elem1", type=str, default="Fe",
                        help="First element")
    parser.add_argument("--elem2", type=str, default="Co",
                        help="Second element")
    parser.add_argument("--structure", type=str, default="bcc",
                        choices=["bcc", "fcc", "hcp"],
                        help="Crystal structure")
    parser.add_argument("--encut", type=int, default=400,
                        help="ENCUT value")
    parser.add_argument("--j-count", type=int, default=None,
                        help="Number of J interactions")
    parser.add_argument("--dist-cutoff", type=float, default=None,
                        help="Distance cutoff for J interactions")
    parser.add_argument("--potcar", type=str, default="",
                        help="Path to POTCAR file")
    parser.add_argument("--extend-poscar", type=str, default="2,2,2",
                        help="Extend POSCAR (nx,ny,nz)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Dry run (don't create files)")
    parser.add_argument("--skip-vasp", action="store_true",
                        help="Skip VASP calculations (use mock)")
    parser.add_argument("--config", type=str, default=None,
                        help="Path to JSON config file for OJ parameters")
    parser.add_argument("--slurm", action="store_true",
                        help="Submit VASP jobs to Slurm")
    parser.add_argument("--ntasks", type=int, default=48,
                        help="Number of tasks for Slurm")
    parser.add_argument("--partition", type=str, default="partCPU",
                        choices=["partCPU", "partGPU"],
                        help="Slurm partition")
    parser.add_argument("--memory", type=str, default="20G",
                        help="Memory per node")
    parser.add_argument("--module", type=str, default="",
                        help="Module to load")
    parser.add_argument("--monitor", action="store_true",
                        help="Monitor Slurm jobs until completion")

    args = parser.parse_args()

    extend_parts = args.extend_poscar.split(",")
    extend_poscar = tuple(int(x) for x in extend_parts)

    workflow = OstravaJWorkflow(
        workdir=args.workdir,
        system_name=args.system,
        potcar_path=args.potcar,
    )

    if args.config:
        log.info(f"Loading OJ config from {args.config}")
        oj_config = OJConfig.from_json(args.config)
        workflow.prepare_input_from_config(
            oj_config=oj_config,
            structure_type=args.structure,
            elem1=args.elem1,
            elem2=args.elem2,
            encut=args.encut,
        )
    else:
        workflow.prepare_input(
            structure_type=args.structure,
            elem1=args.elem1,
            elem2=args.elem2,
            encut=args.encut,
            j_count=args.j_count,
            dist_cutoff=args.dist_cutoff,
            extend_poscar=extend_poscar,
        )

    if not workflow.generate_configs(dry_run=args.dry_run):
        log.error("Config generation failed")
        return 1

    if args.dry_run:
        log.info("Dry run completed")
        return 0

    flip_dirs = list(workflow.run_dir.glob("flip*"))
    if not flip_dirs:
        log.error("No flip directories found")
        return 1

    if args.skip_vasp:
        log.info("Running mock VASP calculations...")
        for flip_dir in flip_dirs:
            run_mock_vasp(flip_dir)
            log.info(f"Mock calculated: {flip_dir}")
        solution = workflow.solve_with_mock()
    elif args.slurm:
        from htvasp.slurm import SlurmJobManager, SlurmConfig
        manager = SlurmJobManager()
        if args.partition == "partGPU":
            slurm_config = manager.get_gpu_config(
                job_name=f"oj-{args.system}",
                ntasks=args.ntasks,
                memory=args.memory,
            )
        else:
            slurm_config = manager.get_cpu_config(
                job_name=f"oj-{args.system}",
                ntasks=args.ntasks,
                memory=args.memory,
            )
        job_ids = workflow.submit_vasp_jobs(
            slurm_config=slurm_config,
            vasp_cmd="srun vasp_std",
            module_name=args.module,
        )
        log.info(f"Submitted {len(job_ids)} jobs: {job_ids}")

        if args.monitor and job_ids:
            log.info("Monitoring jobs until completion...")
            while True:
                pending = []
                for jid in job_ids:
                    status_dict = manager.get_job_status(jid)
                    if status_dict:
                        job_state = status_dict.get("JobState", "UNKNOWN")
                        if job_state in ["PENDING", "RUNNING", "CONFIGURING"]:
                            pending.append(jid)
                if not pending:
                    log.info("All jobs completed")
                    break
                log.info(f"Still running: {pending}")
                import time
                time.sleep(60)
        return 0
    else:
        log.info("=" * 60)
        log.info("VASP calculations need to be run manually")
        log.info(f"Run directory: {workflow.run_dir}")
        log.info("Submit VASP jobs for all flip* directories")
        log.info("=" * 60)
        return 0
    if solution:
        log.info("=" * 60)
        log.info("SOLUTION SUMMARY:")
        log.info(f"  J pairs: {solution.get('J_reprs', [])}")
        log.info(f"  J values (meV): {[f'{j*1000:.4f}' for j in solution.get('Js', [])]}")
        log.info(f"  Tc MFA (K): {solution.get('Tc_MFA', 'N/A')}")
        log.info(f"  Tc RPA (K): {solution.get('Tc_RPA', 'N/A')}")
        log.info(f"  Fitting error: {solution.get('resid', 'N/A')}")
        log.info("=" * 60)
        return 0
    else:
        log.error("Solve failed")
        return 1


if __name__ == "__main__":
    sys.exit(main() or 0)