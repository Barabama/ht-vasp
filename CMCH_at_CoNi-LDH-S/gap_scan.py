#!/usr/bin/env python3
"""
gap 扫描 — 构建不同初始界面间距的本征异质结并提交弛豫

用法:
  # 查看下一步需运行的 gap
  python gap_scan.py --status

  # 构建 POSCAR + 提交 Slurm (逐个提交, 建议一次一个)
  python gap_scan.py --gap 0.5 --submit
  python gap_scan.py --gap 1.5 --submit

  # 仅构建 POSCAR (不提交)
  python gap_scan.py --gap 0.5 --build-only

  # 对比所有 gap 的结果
  python gap_scan.py --compare

说明:
  gap=1.0 已做 (data/gap_scan/intrinsic_gap_1.0/)
  gap=2.0 已做 (data/hetero_intrinsic/)
  本脚本负责: gap=0.5 和 gap=1.5
"""

import argparse
import json
import logging
import shutil
import sys
import time
from pathlib import Path

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.interface import Interface

from htvasp.slurm import SlurmJobManager

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

ROOT = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S")
DATA = ROOT / "data"
FLOW_DIR = Path("/nfs_ssd/tmp")
CONDA = str(ROOT.parent / ".conda")
POSCAR_DIR = DATA / "heterostructures"
SCAN_DIR = DATA / "gap_scan"

# VASP 参数 (与 nscf_hetero.py 一致)
VASP_ARGS_GPU = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-gpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-gpu && srun vasp_gam'",
}

GLOBAL_INCAR = {
    "ENCUT": 520, "PREC": "Accurate", "ALGO": "Fast",
    "NELM": 100, "EDIFF": 1e-6, "ISPIN": 2,
    "MAGMOM": {"Co": 3.0, "Mn": 5.0, "Ni": 2.0, "C": 0.6, "H": 0.6, "O": 0.6},
    "AMIX": 0.1, "BMIX": 1e-4, "AMIX_MAG": 0.4, "BMIX_MAG": 1e-4,
    "LREAL": "Auto", "KPAR": 4, "NCORE": 4,
    "GGA": "PE", "IVDW": 12,
    "LDAU": True, "LDAUTYPE": 2, "LDAUPRINT": 1, "LASPH": True, "LMAXMIX": 4,
    "LDAUL": {"Co": 2, "Mn": 2, "Ni": 2, "C": -1, "H": -1, "O": -1},
    "LDAUU": {"Co": 3.32, "Mn": 5.00, "Ni": 6.20},
    "LDAUJ": {"Co": 0.0, "Mn": 0.0, "Ni": 0.0},
}

RELAX_INCAR = {
    "ISMEAR": 0, "SIGMA": 0.05, "IBRION": 2, "ISIF": 2,
    "NELM": 100, "NSW": 150, "EDIFFG": -0.05,
}

# 已有数据的 gap
EXISTING = {
    0.5: None,
    1.0: SCAN_DIR / "intrinsic_gap_1.0",
    1.5: None,
    2.0: DATA / "hetero_intrinsic",
}

# 晶格匹配参数 (与 build_heterostructure.py 一致)
VAC_BOTTOM = 2.0
VAC_TOP = 12.0
TOTAL_C = 28.0


# ═══════════════════════════════════════════════════════════════
# 构建
# ═══════════════════════════════════════════════════════════════

def load_slabs():
    """读取优化后 slab (与 build_heterostructure.py 一致)"""
    def load(name):
        with open(DATA / name / "static_out.json") as f:
            d = json.load(f)
        return Structure.from_dict(d["output"]["structure"])

    cmch = load("CoMnH2CO5-slab")
    ldh = load("CoNiOH2-slab")

    # 平均晶格
    avg_a = (cmch.lattice.a + ldh.lattice.a) / 2
    avg_b = (cmch.lattice.b + ldh.lattice.b) / 2

    cmch_m = _strain_to(cmch, avg_a, avg_b)
    ldh_m = _strain_to(ldh, avg_a, avg_b)
    return cmch_m, ldh_m, avg_a, avg_b


def _strain_to(slab, target_a, target_b):
    old = slab.lattice
    new_lat = Lattice.from_parameters(
        target_a, target_b, old.c,
        old.alpha, old.beta, old.gamma,
    )
    return Structure(new_lat, slab.species, slab.frac_coords,
                     site_properties=slab.site_properties,
                     coords_are_cartesian=False)


def build_hetero(gap_value: float) -> Structure:
    """用给定 gap 构建本征异质结，返回 Structure (含 selective_dynamics)."""
    cmch_m, ldh_m, avg_a, avg_b = load_slabs()

    sub_thick = cmch_m.cart_coords[:, 2].max() - cmch_m.cart_coords[:, 2].min()
    extra_vac = VAC_BOTTOM + VAC_TOP

    iface = Interface.from_slabs(
        substrate_slab=cmch_m,
        film_slab=ldh_m,
        gap=gap_value,
        vacuum_over_film=extra_vac,
        center_slab=True,
    )

    # shift 到底部真空 = VAC_BOTTOM, c = TOTAL_C
    zs = iface.cart_coords[:, 2]
    labels = np.array(iface.site_properties["interface_label"])
    sub_zs = zs[labels == "substrate"]
    shift = VAC_BOTTOM - sub_zs.min()

    coords_new = iface.cart_coords.copy()
    coords_new[:, 2] += shift
    new_lat = Lattice.from_parameters(
        iface.lattice.a, iface.lattice.b, TOTAL_C,
        iface.lattice.alpha, iface.lattice.beta, iface.lattice.gamma,
    )
    iface_final = Structure(new_lat, iface.species, coords_new,
                            site_properties=iface.site_properties,
                            coords_are_cartesian=True)

    # selective_dynamics: CMCH 底部下半 FFF, 其余 TTT
    zs_final = iface_final.cart_coords[:, 2]
    labels_final = np.array(iface_final.site_properties["interface_label"])
    sub_mask = labels_final == "substrate"
    sub_zs = zs_final[sub_mask]
    sub_z_mid = (sub_zs.min() + sub_zs.max()) / 2

    sd = []
    for i in range(len(iface_final)):
        if sub_mask[i] and zs_final[i] < sub_z_mid:
            sd.append([True, True, False])     # TTF = z固定, xy可动
        else:
            sd.append([True, True, True])      # TTT

    iface_final.add_site_property("selective_dynamics", sd)

    # 报告
    film_zs = zs_final[~sub_mask]
    true_gap = film_zs.min() - sub_zs.max()
    n_fff = sum(1 for d in sd if d == [False, False, False])
    n_ttf = sum(1 for d in sd if d == [True, True, False])
    n_free = sum(1 for d in sd if d == [True, True, True])

    log.info("构建 gap=%.1f: 原子=%d, CMCH z=[%.2f,%.2f], LDH z=[%.2f,%.2f], 界面间隙=%.2f, %d TTF + %d TTT",
             gap_value, len(iface_final), sub_zs.min(), sub_zs.max(),
             film_zs.min(), film_zs.max(), true_gap, n_ttf, n_free)

    return iface_final


# ═══════════════════════════════════════════════════════════════
# 提交
# ═══════════════════════════════════════════════════════════════

def submit_gap(gap_value: float):
    """构建 POSCAR + 准备 INCAR + 提交 Slurm."""
    store_dir = SCAN_DIR / f"intrinsic_gap_{gap_value}"
    if store_dir.exists():
        log.warning("目录已存在: %s, 跳过", store_dir)
        return

    # 1. 构建异质结
    struct = build_hetero(gap_value)

    # 2. 准备输入文件
    store_dir.mkdir(parents=True)

    # POSCAR
    struct.to(str(store_dir / "POSCAR"), fmt="poscar")

    # INCAR (relax)
    incar_lines = [
        "SYSTEM = CMCH@CoNi-LDH intrinsic gap scan",
        "ENCUT = 520",
        "PREC = Accurate",
        "ALGO = Fast",
        "NELM = 100",
        "EDIFF = 1E-6",
        "EDIFFG = -0.05",
        "ISPIN = 2",
        "MAGMOM = 6*5.0 12*3.0 6*2.0 36*0.6 6*0.6 54*0.6",
        "AMIX = 0.1",
        "BMIX = 1E-4",
        "AMIX_MAG = 0.4",
        "BMIX_MAG = 1E-4",
        "LREAL = Auto",
        "KPAR = 4",
        "NCORE = 4",
        "GGA = PE",
        "IVDW = 12",
        "LDAU = .TRUE.",
        "LDAUTYPE = 2",
        "LDAUPRINT = 1",
        "LASPH = .TRUE.",
        "LMAXMIX = 4",
        "LDAUL = 2 2 2 -1 -1 -1",
        "LDAUU = 3.32 5.00 6.20 0.0 0.0 0.0",
        "LDAUJ = 0.0 0.0 0.0 0.0 0.0 0.0",
        "ISMEAR = 0",
        "SIGMA = 0.05",
        "IBRION = 2",
        "ISIF = 2",
        "NSW = 150",
        "LWAVE = .TRUE.",
        "LCHARG = .FALSE.",
    ]
    (store_dir / "INCAR").write_text("\n".join(incar_lines) + "\n")

    # KPOINTS: gamma-centered, 1x1x1 (超胞较大)
    (store_dir / "KPOINTS").write_text("Auto\n0\nGamma\n1 1 1\n")

    # POTCAR: 从现有异质结复制 (需解压, VASP 不读 .gz)
    src_potcar = DATA / "hetero_intrinsic" / "1-relax_1" / "POTCAR.gz"
    if src_potcar.exists():
        import gzip
        with gzip.open(src_potcar, "rb") as f_in:
            (store_dir / "POTCAR").write_bytes(f_in.read())
        log.info("POTCAR → %s/POTCAR (%d bytes)", store_dir, (store_dir / "POTCAR").stat().st_size)
    else:
        log.warning("POTCAR 不可用, 需手动复制")

    log.info("输入文件就绪: %s  (%d files)", store_dir, len(list(store_dir.iterdir())))

    # 3. 提交 Slurm
    job_name = f"gap_scan_{gap_value}"
    command = (
        f"cd {store_dir} && "
        f". /etc/profile.d/modules.sh && module load vasp-gpu && "
        f"srun vasp_std"
    )
    output_log = str(ROOT / "logs" / f"{job_name}.log")

    manager = SlurmJobManager()
    config = manager.get_gpu_config(
        job_name=job_name,
        output_log=output_log,
        nodes=1, ntasks=1, memory="20G",
        conda_env=CONDA,
    )

    jid = manager.submit_command(command=command, config=config, workdir=str(store_dir))
    if jid:
        log.info("提交 gap=%.1f → job %s", gap_value, jid)
    else:
        log.error("提交失败 gap=%.1f", gap_value)

    return jid


# ═══════════════════════════════════════════════════════════════
# 状态 / 对比
# ═══════════════════════════════════════════════════════════════

def print_status():
    """显示各 gap 的状态"""
    print(f"\n{'='*70}")
    print(f"  Gap Scan Status")
    print(f"{'='*70}")
    print(f"  {'Gap':>6s}  {'状态':>20s}  {'Dir':>30s}")
    print(f"  {'─'*6}  {'─'*20}  {'─'*30}")

    for g in [0.5, 1.0, 1.5, 2.0]:
        if g == 2.0:
            d = DATA / "hetero_intrinsic"
        else:
            d = SCAN_DIR / f"intrinsic_gap_{g}"

        if d.exists():
            contcar = d / "CONTCAR"
            is_running = False
            if contcar.exists() and (time.time() - contcar.stat().st_mtime) < 3600:
                is_running = True

            oszicar = d / "OSZICAR"
            n_steps = 0
            if oszicar.exists():
                with open(oszicar) as f:
                    for line in f:
                        if line.strip().split()[0].isdigit() and 'F=' in line:
                            n_steps += 1

            if is_running:
                status = f"运行中 ({n_steps}步)"
            elif contcar.exists() and not is_running:
                status = f"已完成 ({n_steps}步)"
            else:
                status = f"准备中"
        else:
            status = "未开始"

        print(f"  {g:>6.1f}  {status:>20s}  {str(d):>30s}")

    print()


def print_compare():
    """对比所有 gap 的结果"""
    print(f"\n{'='*70}")
    print(f"  Gap Scan Comparison")
    print(f"{'='*70}")
    print(f"  {'Gap':>6s}  {'界面间距(Å)':>12s}  {'能量(eV)':>14s}  {'步数':>6s}  {'maxF(eV/Å)':>12s}")
    print(f"  {'─'*6}  {'─'*12}  {'─'*14}  {'─'*6}  {'─'*12}")

    for g in [0.5, 1.0, 1.5, 2.0]:
        if g == 2.0:
            d = DATA / "hetero_intrinsic" / "2-relax_2"
        else:
            d = SCAN_DIR / f"intrinsic_gap_{g}"

        if not d.exists():
            print(f"  {g:>6.1f}  {'N/A':>12s}")
            continue

        # 读 CONTCAR
        contcar = d / "CONTCAR"
        if not contcar.exists():
            contcar = d / "CONTCAR.gz"
        if not contcar.exists():
            print(f"  {g:>6.1f}  {'无CONTCAR':>12s}")
            continue

        # 界面间距
        import gzip
        try:
            if str(contcar).endswith(".gz"):
                with gzip.open(contcar, "rt") as f:
                    s = Structure.from_str(f.read(), fmt="poscar")
            else:
                s = Structure.from_file(contcar)

            sites = sorted([(site.coords[2], site.species_string) for site in s], key=lambda x: x[0])
            zs = [x[0] for x in sites]
            z_min_mat = min(zs) + (max(zs)-min(zs))*0.05
            z_max_mat = max(zs) - (max(zs)-min(zs))*0.05

            gap_final = 0
            for i in range(len(zs)-1):
                dz = zs[i+1]-zs[i]
                z_mid = (zs[i]+zs[i+1])/2
                if dz > gap_final and 0.3 < dz and z_mid > z_min_mat and z_mid < z_max_mat:
                    gap_final = dz
        except Exception:
            gap_final = -1

        # 能量
        oszicar = d / "OSZICAR"
        energy = None
        n_steps = 0
        if oszicar.exists():
            with open(oszicar) as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 8 and parts[0].isdigit() and parts[1] == 'F=' and parts[3] == 'E0=':
                        energy = float(parts[4])
                        n_steps += 1

        # max force (从 OUTCAR)
        max_force = None
        outcar = d / "OUTCAR"
        if outcar.exists():
            with open(outcar) as f:
                for line in f:
                    if "FORCES" in line and "max" in line:
                        try:
                            max_force = float(line.split()[3])
                        except:
                            pass

        gap_s = f"{gap_final:.4f}" if gap_final > 0 else "N/A"
        e_s = f"{energy:.4f}" if energy else "N/A"
        f_s = f"{max_force:.4f}" if max_force else "N/A"
        print(f"  {g:>6.1f}  {gap_s:>12s}  {e_s:>14s}  {n_steps:>6d}  {f_s:>12s}")

    print()


# ═══════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gap scan for intrinsic heterojunction")
    parser.add_argument("--status", action="store_true", help="显示所有 gap 的状态")
    parser.add_argument("--compare", action="store_true", help="对比所有 gap 的结果")
    parser.add_argument("--gap", type=float, help="初始间隙 (0.5 或 1.5)")
    parser.add_argument("--submit", action="store_true", help="构建 + 提交 Slurm")
    parser.add_argument("--build-only", action="store_true", help="仅构建 POSCAR, 不提交")
    args = parser.parse_args()

    if args.status:
        print_status()
    elif args.compare:
        print_compare()
    elif args.gap:
        if args.gap not in (0.5, 1.5):
            log.error("gap 必须是 0.5 或 1.5 (1.0和2.0已做)")
            sys.exit(1)
        if args.build_only:
            struct = build_hetero(args.gap)
            out_dir = SCAN_DIR / f"intrinsic_gap_{args.gap}"
            out_dir.mkdir(parents=True, exist_ok=True)
            struct.to(str(out_dir / "POSCAR"), fmt="poscar")
            log.info("POSCAR → %s/POSCAR", out_dir)
        elif args.submit:
            submit_gap(args.gap)
        else:
            parser.print_help()
    else:
        parser.print_help()
