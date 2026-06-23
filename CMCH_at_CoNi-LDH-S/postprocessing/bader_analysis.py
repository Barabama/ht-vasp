"""
Bader 电荷分析 — 10 体系批量分析

在临时目录解压 AECCAR0/AECCAR2 → 调用 bader CLI → 提取原子电荷/电荷转移/磁矩

用法:
  python bader_analysis.py                      # 全部 10 体系
  python bader_analysis.py --systems CoNiOH2     # 指定体系

输出:
  data/{system}/bader_results.json  (单体系详细结果)
"""

import argparse
import json
import logging
import shutil
import tempfile
from pathlib import Path

import numpy as np
from pymatgen.command_line.bader_caller import BaderAnalysis
from monty.json import MontyEncoder

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

DATA_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/data")
BADER_BIN = "/nfs_ssd/softwares/bader/bin/bader"

SYSTEMS_ALL = [
    "CoNiOH2", "CoNiOH2S-noH", "CoMnH2CO5",
    "CoNiOH2-slab", "CoNiOH2S-noH-slab", "CoNiOH2S-noH-slab-flip", "CoMnH2CO5-slab",
    "hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed",
]

ZVAL = {"Co": 9.0, "Ni": 10.0, "O": 6.0, "H": 1.0, "C": 4.0, "Mn": 7.0, "S": 6.0}


def analyze_bader(name: str) -> dict:
    """对指定体系运行 Bader 分析。

    流程:
      1. 在临时目录解压 AECCAR0/AECCAR2/POTCAR
      2. BaderAnalysis.from_path() 调用 chgsum + bader CLI
      3. 提取原子电荷、电荷转移、磁矩
      4. 清理临时文件
    """
    static_dir = DATA_DIR / name / "3-static"
    if not static_dir.exists():
        return {"name": name, "error": f"3-static not found: {static_dir}"}

    required = ["AECCAR0.gz", "AECCAR2.gz", "POTCAR.gz"]
    for f in required:
        if not (static_dir / f).exists():
            return {"name": name, "error": f"Missing {f}"}

    try:
        with tempfile.TemporaryDirectory(prefix=f"bader_{name}_") as tmpdir:
            tmp = Path(tmpdir)

            # 解压输入文件（bader CLI 不能读 .gz）
            for f in required + ["CHGCAR.gz"]:
                _gunzip(static_dir / f, tmp / f.replace(".gz", ""))

            # 运行 Bader 分析
            # from_path 不传 bader_path，通过 __init__ 后的内部属性设置
            bader = BaderAnalysis.from_path(str(tmp))

            # 设置 bader binary 路径（from_path 不支持此参数）
            bader.bader_path = BADER_BIN

            result = {"name": name}
            result["charge"] = []
            result["charge_transfer"] = []
            result["atomic_volume"] = []

            # 逐原子提取电荷
            n_atoms = len(bader.chgcar.structure)
            for i in range(n_atoms):
                chg = bader.get_charge(i)
                nelect = ZVAL.get(bader.chgcar.structure[i].species_string, 6.0)
                # get_charge_transfer returns charge - ZVAL
                # 学术惯例: charge_transfer = ZVAL - charge, 正值=失电子
                ct = nelect - chg

                result["charge"].append(float(chg))
                result["charge_transfer"].append(float(ct))

            # atom volume from bader output (ACF.dat column 4)
            # 解析 ACF.dat 获取原子体积
            acf_path = tmp / "ACF.dat"
            if acf_path.exists():
                with open(acf_path) as f:
                    lines = f.readlines()
                # ACF.dat format: # X Y Z CHARGE MIN_DIST ATOMIC_VOL
                for line in lines[2:]:  # skip headers
                    parts = line.strip().split()
                    if len(parts) >= 6:
                        result["atomic_volume"].append(float(parts[5]))
            else:
                result["atomic_volume"] = []

            # 磁矩需从 OUTCAR 获取（bader CLI 不支持自旋极化 CFG 格式）
            # 已在 comp_hetero.py 中通过 Outcar.total_mag 提取

            # 保存到持久目录
            out_dir = DATA_DIR / name / "bader_analysis"
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "bader_results.json").write_text(
                json.dumps(result, indent=2, cls=MontyEncoder)
            )

            n_elements = len(result["charge"])
            log.info("  ✅ %s: %d atoms, avg |ct| = %.4f e",
                      name, n_elements,
                      np.mean(np.abs(result["charge_transfer"])))

            return result

    except Exception as e:
        log.error("  ❌ %s: %s", name, e)
        return {"name": name, "error": str(e)}


def _gunzip(src: Path, dst: Path):
    import gzip
    with gzip.open(src, "rb") as f_in:
        with open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def main():
    parser = argparse.ArgumentParser(description="Bader charge analysis")
    parser.add_argument("--systems", type=str, nargs="+",
                        help="体系名（默认全部）")
    parser.add_argument("--parallel", action="store_true",
                        help="并行运行")
    args = parser.parse_args()

    systems = args.systems or SYSTEMS_ALL

    log.info("Bader analysis for %d systems", len(systems))

    results = {}
    for name in systems:
        log.info("Processing %s ...", name)
        results[name] = analyze_bader(name)

    # 汇总
    log.info("\n" + "=" * 60)
    log.info("SUMMARY")
    log.info("=" * 60)
    for name, r in results.items():
        if "error" in r:
            log.info("  ✗ %s: %s", name, r["error"])
        elif "charge" in r:
            ct = r["charge_transfer"]
            s_idx = next((i for i in range(len(ct)) if _element_at(results, name, i) == "S"), None)
            s_ct = ct[s_idx] if s_idx is not None else None
            s_str = f", S ct={s_ct:+.3f}" if s_ct is not None else ""
            log.info("  ✓ %s: %d atoms, ct range [%+.3f, %+.3f]%s",
                      name, len(ct), min(ct), max(ct), s_str)

    log.info("Done.")


def _element_at(results, name, idx):
    """辅助: 获取指定原子索引的元素符号."""
    try:
        static_dir = DATA_DIR / name / "3-static"
        import gzip
        from pymatgen.core import Structure
        with gzip.open(static_dir / "CONTCAR.gz", "rt") as f:
            struct = Structure.from_str(f.read(), fmt="poscar")
        return struct[idx].species_string
    except Exception:
        return None


if __name__ == "__main__":
    main()
