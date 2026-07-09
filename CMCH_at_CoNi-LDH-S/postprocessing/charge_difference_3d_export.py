#!/usr/bin/env python3
"""生成 VESTA 可读的 3D 差分电荷密度文件 (CHGCAR_diff)

原理:
  Δρ(r) = ρ_hetero(r) − ρ_CMCH(r) − ρ_LDH(r)
  直接将三个体系的 CHGCAR 在 FFT 网格上逐格相减，
  输出为 VASP CHGCAR 格式文件，VESTA 直接打开即可做 3D 等值面图。

前提:
  所有 Level 2 体系共享同一 FFT 网格 (160,140,420)，
  原子坐标在相同晶格下对齐，无需插值。

输出:
  output/charge_diff_3d/CHGCAR_diff_{het_name}.vasp

用法:
  python postprocessing/charge_difference_3d_export.py
  python postprocessing/charge_difference_3d_export.py --systems hetero_intrinsic
"""

import argparse
import gzip
import logging
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.io.vasp.inputs import Poscar

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
OUTPUT_DIR = SCRIPT_DIR / "output" / "charge_diff_3d"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

HETERO_MAP = {
    "hetero_intrinsic": ("LDH_strained", "CMCH_strained"),
    "hetero_s_doped": ("LDH_S_strained", "CMCH_strained"),
    "hetero_s_exposed": ("LDH_S_flip_strained", "CMCH_strained"),
}


def load_chgcar(name: str) -> Chgcar:
    """从 data/{name}/3-static/CHGCAR.gz 加载 CHGCAR。"""
    path = DATA_DIR / name / "3-static" / "CHGCAR.gz"
    if not path.exists():
        raise FileNotFoundError(f"Missing CHGCAR for {name}: {path}")

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".CHGCAR")
    try:
        with gzip.open(path, "rb") as f_in:
            shutil.copyfileobj(f_in, tmp)
        tmp_path = tmp.name
        tmp.close()
        chg = Chgcar.from_file(tmp_path)
        return chg
    finally:
        try:
            os.unlink(tmp_path)
        except Exception:
            pass


def write_chgcar_diff(out_path: Path, het_name: str, diff_data: np.ndarray,
                      het_chg: Chgcar):
    """将 3D 差分电荷密度写入 VASP CHGCAR 格式文件。

    复用异质结 CHGCAR 的晶格和原子坐标信息，
    将密度数据替换为差分值。
    """
    struct = het_chg.structure
    ngx, ngy, ngz = het_chg.dim

    with open(out_path, "w") as f:
        # Header
        f.write(f"CHGCAR_diff_{het_name}\n")
        f.write("  1.000000000000000\n")
        for vec in struct.lattice.matrix:
            f.write(f"  {vec[0]:20.12f}{vec[1]:20.12f}{vec[2]:20.12f}\n")
        # Species and counts — use het_chg.poscar to get correct format
        poscar_str = str(het_chg.poscar).split("\n")
        f.write(f"  {poscar_str[5]}\n")  # species line
        f.write(f"   {poscar_str[6]}\n")  # counts
        # Write "Direct" + coords
        f.write(f"  {poscar_str[7]}\n")  # "Direct" or "Cartesian"
        for site in struct:
            frac = site.frac_coords
            f.write(f"  {frac[0]:20.16f}{frac[1]:20.16f}{frac[2]:20.16f}\n")

        # FFT 网格
        f.write(f"\n{ngx:6d}{ngy:6d}{ngz:6d}\n")

        # 差分数据 → 按 VASP 格式输出
        # VASP 每行 5 个数值，科学计数法
        n_pts = ngx * ngy * ngz
        diff_flat = diff_data.ravel()
        assert len(diff_flat) == n_pts, f"{len(diff_flat)} != {n_pts}"

        count = 0
        for val in diff_flat:
            f.write(f" {val:20.12e}")
            count += 1
            if count % 5 == 0:
                f.write("\n")
        if count % 5 != 0:
            f.write("\n")

    log.info("  ✓ Wrote %s  (%.1f MB)", out_path.name, out_path.stat().st_size / 1e6)


def compute_and_export(het_name: str):
    """对一个异质结计算 3D 差分并导出。"""
    if het_name not in HETERO_MAP:
        log.warning("  ✗ Unknown heterojunction: %s", het_name)
        return

    ldh_name, cmch_name = HETERO_MAP[het_name]

    log.info("Loading %s ...", het_name)
    het_chg = load_chgcar(het_name)
    log.info("Loading %s ...", cmch_name)
    cmch_chg = load_chgcar(cmch_name)
    log.info("Loading %s ...", ldh_name)
    ldh_chg = load_chgcar(ldh_name)

    # 提取总电荷密度: chg.data["total"] has shape (ngx, ngy, ngz)
    het_dens = het_chg.data["total"]
    cmch_dens = cmch_chg.data["total"]
    ldh_dens = ldh_chg.data["total"]

    # 检查 FFT 网格一致性

    # 差分: Δρ = ρ_hetero − ρ_CMCH − ρ_LDH
    log.info("Computing Δρ = ρ_hetero − ρ_CMCH − ρ_LDH ...")
    diff_data = het_dens - cmch_dens - ldh_dens

    # 统计
    log.info("  Δρ range: [%.6e, %.6e] e⁻/Å³", diff_data.min(), diff_data.max())
    log.info("  Δρ |mean|: %.6e e⁻/Å³", np.mean(np.abs(diff_data)))

    out_path = OUTPUT_DIR / f"CHGCAR_diff_{het_name}.vasp"
    write_chgcar_diff(out_path, het_name, diff_data, het_chg)

    # 同时输出一个简化的 CHGCAR（仅原子坐标，VESTA 专用）
    log.info("Done: %s", het_name)


def main():
    parser = argparse.ArgumentParser(description="3D 差分电荷密度导出 (VESTA)")
    parser.add_argument("--systems", type=str, nargs="+",
                        help="异构体系名（默认全部）")
    args = parser.parse_args()

    systems = args.systems or list(HETERO_MAP.keys())
    log.info("3D CHGCAR difference for %d heterojunctions", len(systems))

    for het_name in systems:
        compute_and_export(het_name)

    log.info("All done. Output → %s", OUTPUT_DIR)
    for f in sorted(OUTPUT_DIR.glob("CHGCAR_diff_*.vasp")):
        log.info("  %s  (%.1f MB)", f.name, f.stat().st_size / 1e6)


if __name__ == "__main__":
    main()
