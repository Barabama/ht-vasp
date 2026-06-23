"""
偶极修正补跑 — 对已有 static 结果重新跑静态计算，仅加 LDIPOL/IDIPOL.

全部文件操作使用 Python (gzip/shutil)，不产生中间 shell 脚本文件.
用法:
  # 本地准备 + 提交 Slurm (含 prepare → vasp → save 全流程)
  python patch_dipole.py --submit

  # 或在计算节点手动执行其中一步:
  python patch_dipole.py --prepare CoNiOH2S-noH-slab-flip   # 准备输入文件
  vasp_std                                                   # 运行 VASP
  python patch_dipole.py --save   CoNiOH2S-noH-slab-flip   # 压缩结果
"""

import argparse
import gzip
import logging
import shutil
import subprocess
import sys
from pathlib import Path

from htvasp.slurm import SlurmJobManager

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

DATA_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/data")
SSD_TMP = Path("/nfs_ssd/tmp")

SYSTEMS = [
    # "CoNiOH2S-noH-slab-flip",  # ✅ done (verified)
    "CoNiOH2-slab",
    "CoNiOH2S-noH-slab",
    "CoMnH2CO5-slab",
    "hetero_intrinsic",
    "hetero_s_doped",
    "hetero_s_exposed",
]

VASP_INIT = ". /etc/profile.d/modules.sh && module load vasp-gpu"


# ═══════════════════════════════════════════════
# Prepare
# ═══════════════════════════════════════════════

def prepare(name: str) -> Path:
    """解压输入文件, 写入带偶极修正的 INCAR.

    Returns: 工作目录 Path.
    """
    src = DATA_DIR / name / "3-static"
    work = SSD_TMP / f"{name}-dipole"

    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    _gunzip(src / "CONTCAR.gz",  work / "POSCAR")
    _gunzip(src / "POTCAR.gz",   work / "POTCAR")
    _gunzip(src / "KPOINTS.gz",  work / "KPOINTS")
    _gunzip(src / "WAVECAR.gz",  work / "WAVECAR")
    _gunzip(src / "CHGCAR.gz",   work / "CHGCAR")

    with gzip.open(src / "INCAR.gz", "rt") as f:
        incar = f.read()
    incar += "\n# Dipole correction (patch_dipole.py)\n"
    incar += "LDIPOL = .TRUE.\n"
    incar += "IDIPOL = 3\n"
    (work / "INCAR").write_text(incar)

    # WAVECAR from non-dipole calc is incompatible with LDIPOL boundary condition
    # → delete it so VASP starts fresh wavefunctions (ISTART=1→auto 0)
    wf = work / "WAVECAR"
    if wf.exists():
        wf.unlink()
    # Keep CHGCAR as initial charge density guess (accelerates convergence)

    n = len(list(work.iterdir()))
    log.info("prepare %s → %s  (%d files)", name, work, n)
    return work


# ═══════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════

def save(name: str):
    """压缩 LOCPOT/OUTCAR 等结果写入持久目录后删除临时目录."""
    work = SSD_TMP / f"{name}-dipole"
    if not work.exists():
        log.warning("save %s: workdir not found %s", name, work)
        return

    dst = DATA_DIR / name / "3-static-dipole"
    dst.mkdir(parents=True, exist_ok=True)

    for f in ("LOCPOT", "OUTCAR", "OSZICAR", "INCAR", "vasprun.xml"):
        fp = work / f
        if fp.exists():
            with open(fp, "rb") as f_in:
                with gzip.open(dst / f"{f}.gz", "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

    shutil.rmtree(work, ignore_errors=True)
    log.info("save %s → %s  (%d files)", name, dst, len(list(dst.iterdir())))


# ═══════════════════════════════════════════════
# Submit (prepare + vasp + save 全流程)
# ═══════════════════════════════════════════════

def submit_all():
    """每个体系独立提交一个 GPU job，并行运行."""
    manager = SlurmJobManager()

    for name in SYSTEMS:
        # ── Python 准备输入文件 ──
        work = prepare(name)

        # ── sbatch 命令（仅有的一次 shell 调用就是跑 VASP） ──
        save_cmd = (
            f"cd /nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S && "
            f"/nfs_hdd/2025/gaominliang/ht-vasp/.conda/bin/python patch_dipole.py "
            f"--save {name}"
        )
        run_cmd = (
            f"{VASP_INIT} && cd {work} && srun vasp_std && {save_cmd}"
        )

        output_log = str(DATA_DIR.parent / "logs" / f"{name}-dipole.log")
        config = manager.get_gpu_config(
            job_name=f"{name}-dipole",
            output_log=output_log,
            nodes=1, ntasks=1, memory="20G",
        )

        jid = manager.submit_command(command=run_cmd, config=config)
        if jid:
            log.info("submit %s → job %s", name, jid)
            prev_job_id = jid
        else:
            log.error("submit FAILED for %s", name)


# ═══════════════════════════════════════════════
# Utils
# ═══════════════════════════════════════════════

def _gunzip(src: Path, dst: Path):
    with gzip.open(src, "rb") as f_in:
        with open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


# ═══════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Dipole correction patch")
    parser.add_argument("--submit", action="store_true", help="准备 + 提交 Slurm")
    parser.add_argument("--prepare", type=str, help="仅准备: 体系名", metavar="NAME")
    parser.add_argument("--save", type=str, help="仅保存: 体系名", metavar="NAME")
    args = parser.parse_args()

    if args.prepare:
        prepare(args.prepare)
    elif args.save:
        save(args.save)
    elif args.submit:
        submit_all()
    else:
        parser.print_help()
