"""
Bader 电荷分析 — 10 体系批量分析

在临时目录解压 AECCAR0/AECCAR2 → 调用 bader CLI → 提取原子电荷/电荷转移

输出:
  output/bader_summary.csv                  逐体系汇总表
  output/bader_atoms_{system}.csv           逐原子详细数据
  output/bader.json                         完整数值存档
  data/{system}/bader_analysis/bader_results.json  (单体系结果，支持快速恢复)

用法:
  python bader_analysis.py                      # 全部 10 体系
  python bader_analysis.py --systems CoNiOH2    # 指定体系
  python bader_analysis.py --force              # 强制重新计算
"""

import argparse
import gzip
import json
import logging
import shutil
import tempfile
from pathlib import Path

import numpy as np
from pymatgen.command_line.bader_caller import BaderAnalysis
from pymatgen.core import Structure

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
OUTPUT_DIR = SCRIPT_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

BADER_BIN = "bader"  # assume on PATH

# 系统命名 → 数据目录映射（解耦命名与目录结构）
DATA_DIR_MAP = {
    "LDH_bulk": "CoNiOH2",
    "LDH_S_bulk": "CoNiOH2S-noH",
    "CMCH_bulk": "CMCH_bulk",
    "CMCH_strained": "CMCH_strained",
    "LDH_strained": "LDH_strained",
    "LDH_S_strained": "LDH_S_strained",
    "LDH_S_flip_strained": "LDH_S_flip_strained",
    "LDH2_strained": "LDH2_strained",
    "LDH_strained_Co+1": "LDH_strained_Co+1",
    "LDH_strained_Co-1": "LDH_strained_Co-1",
    "LDH_S_strained_Co+1": "LDH_S_strained_Co+1",
    "LDH_S_strained_Co-1": "LDH_S_strained_Co-1",
    "hetero_intrinsic": "hetero_intrinsic",
    "hetero_s_doped": "hetero_s_doped",
    "hetero_s_exposed": "hetero_s_exposed",
}

SYSTEMS_ALL = [
    "LDH_bulk", "LDH_S_bulk", "CMCH_bulk",
    "CMCH_strained", "LDH_strained", "LDH_S_strained", "LDH_S_flip_strained", "LDH2_strained",
    "hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed",
]

SYSTEM_LABELS = {
    "LDH_bulk": "LDH (Bulk)", "LDH_S_bulk": "LDH+S (Bulk)", "CMCH_bulk": "CMCH (Bulk)",
    "CMCH_strained": "CMCH (Strained)", "LDH_strained": "LDH (Strained)",
    "LDH_S_strained": "LDH+S (Strained)", "LDH_S_flip_strained": "LDH+S (Strained-flip)",
    "LDH2_strained": "LDH (Bilayer)",
    "hetero_intrinsic": "Intrinsic Het.", "hetero_s_doped": "S-Doped Het.", "hetero_s_exposed": "S-Exposed Het.",
}

# ═══════════════════════════════════════════════
# Bader 分析
# ═══════════════════════════════════════════════

def _parse_potcar_zvals(potcar_gz_path: Path) -> dict[str, float]:
    """从 POTCAR.gz 解析 ZVAL，避免硬编码错误（如 Mn_pv 的 ZVAL=13 而非 7）。"""
    import gzip
    zvals = {}
    with gzip.open(potcar_gz_path, "rt") as f:
        text = f.read()
    for block in text.split("End of Dataset"):
        titel = None
        zval = None
        for line in block.split("\n"):
            if "TITEL" in line:
                titel = line.strip().split()[-2]  # e.g. Co, Mn_pv
            if "ZVAL" in line:
                zval = float(line.strip().split("=")[-1].split()[0])
        if titel and zval:
            # Mn_pv → Mn (用元素符号前半部分匹配)
            elem = titel.split("_")[0]
            zvals[elem] = zval
    return zvals


def _dir(name: str) -> Path:
    """返回系统名对应的数据目录（通过 DATA_DIR_MAP 解耦命名与目录路径）。"""
    return DATA_DIR / DATA_DIR_MAP.get(name, name)


def analyze_bader(name: str) -> dict:
    """对指定体系运行 Bader 分析，返回结果字典."""
    static_dir = _dir(name) / "3-static"
    if not static_dir.exists():
        return {"name": name, "error": f"3-static not found: {static_dir}"}

    required = ["AECCAR0.gz", "AECCAR2.gz", "POTCAR.gz"]
    for f in required:
        if not (static_dir / f).exists():
            return {"name": name, "error": f"Missing {f}"}

    # 从 POTCAR 解析 ZVAL（取代硬编码的 ZVAL 字典）
    potcar_zvals = _parse_potcar_zvals(static_dir / "POTCAR.gz")

    try:
        with tempfile.TemporaryDirectory(prefix=f"bader_{name}_") as tmpdir:
            tmp = Path(tmpdir)

            for f in required + ["CHGCAR.gz"]:
                _gunzip(static_dir / f, tmp / f.replace(".gz", ""))

            bader = BaderAnalysis.from_path(str(tmp))
            bader.bader_path = BADER_BIN

            result = {"name": name, "charge": [], "charge_transfer": [], "atomic_volume": []}

            n_atoms = len(bader.chgcar.structure)
            for i in range(n_atoms):
                chg = bader.get_charge(i)
                elem = bader.chgcar.structure[i].species_string
                if elem not in potcar_zvals:
                    log.warning("ZVAL not found for element %s in POTCAR, falling back to 6.0", elem)
                nelect = potcar_zvals.get(elem, 6.0)
                result["charge"].append(float(chg))
                result["charge_transfer"].append(float(nelect - chg))

            acf_path = tmp / "ACF.dat"
            if acf_path.exists():
                with open(acf_path) as f:
                    lines = f.readlines()
                for line in lines[2:]:
                    parts = line.strip().split()
                    if len(parts) >= 6:
                        result["atomic_volume"].append(float(parts[5]))

            return result

    except Exception as e:
        log.error("  ❌ %s: %s", name, e)
        return {"name": name, "error": str(e)}


def _gunzip(src: Path, dst: Path):
    with gzip.open(src, "rb") as f_in:
        with open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def _load_cached_result(name: str) -> dict | None:
    """从 data/{name}/bader_analysis/bader_results.json 加载已有结果."""
    path = _dir(name) / "bader_analysis" / "bader_results.json"
    if not path.exists():
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return None


def _save_cached_result(name: str, result: dict):
    """保存结果到 data/{name}/bader_analysis/bader_results.json."""
    out_dir = _dir(name) / "bader_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "bader_results.json").write_text(json.dumps(result, indent=2))


def _get_structure(name: str) -> Structure | None:
    contcar = _dir(name) / "3-static" / "CONTCAR.gz"
    if not contcar.exists():
        return None
    try:
        with gzip.open(contcar, "rt") as f:
            return Structure.from_str(f.read(), fmt="poscar")
    except Exception:
        return None


# ═══════════════════════════════════════════════
# 数据导出 (CSV / JSON)
# ═══════════════════════════════════════════════

def export_summary_csv(results: dict):
    """汇总表 → output/bader_summary.csv"""
    rows = []
    for name in SYSTEMS_ALL:
        r = results.get(name, {})
        if "error" in r:
            continue
        ct = r.get("charge_transfer", [])
        avg_ct = float(np.mean(np.abs(ct))) if ct else None
        rows.append({
            "system": name,
            "label": SYSTEM_LABELS.get(name, name),
            "n_atoms": len(ct),
            "charge_transfer_mean (e)": float(np.mean(ct)) if ct else "",
            "charge_transfer_abs_mean (e)": avg_ct if avg_ct is not None else "",
            "charge_transfer_min (e)": min(ct) if ct else "",
            "charge_transfer_max (e)": max(ct) if ct else "",
        })

    if not rows:
        return
    header = list(rows[0].keys())
    with open(OUTPUT_DIR / "bader_summary.csv", "w") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in header) + "\n")
    log.info("Data  → output/bader_summary.csv")


def export_atoms_csv(results: dict):
    """逐原子数据 → output/bader_atoms_{system}.csv"""
    for name in SYSTEMS_ALL:
        r = results.get(name, {})
        if "error" in r or "charge" not in r:
            continue

        struct = _get_structure(name)
        charges = r["charge"]
        cts = r["charge_transfer"]
        vols = r.get("atomic_volume", [])

        with open(OUTPUT_DIR / f"bader_atoms_{name}.csv", "w") as f:
            f.write("index,element,charge (e),charge_transfer (e),atomic_volume (A^3)\n")
            for i in range(len(charges)):
                elem = struct[i].species_string if struct and i < len(struct) else "?"
                vol = vols[i] if i < len(vols) else ""
                ct = cts[i] if i < len(cts) else ""
                f.write(f"{i+1},{elem},{charges[i]:.6f},{ct:.6f},{vol}\n")
    log.info("Data  → output/bader_atoms_*.csv")


def export_results_json(results: dict):
    """完整数值存档 → output/bader.json"""
    with open(OUTPUT_DIR / "bader.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    log.info("Data  → output/bader.json")


def load_results() -> dict | None:
    """从 output/bader.json 加载上次的汇总结果."""
    path = OUTPUT_DIR / "bader.json"
    if not path.exists():
        return None
    try:
        with open(path) as f:
            data = json.load(f)
        log.info("Loaded from %s", path.name)
        return data
    except Exception as e:
        log.warning("Cannot load %s: %s", path.name, e)
        return None


# ═══════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Bader charge analysis")
    parser.add_argument("--systems", type=str, nargs="+", help="体系名（默认全部）")
    parser.add_argument("--force", action="store_true", help="强制重新计算（忽略已有结果）")
    args = parser.parse_args()

    systems = args.systems or SYSTEMS_ALL

    # 优先从 output JSON 恢复
    results = None if args.force else load_results()

    if results is None:
        log.info("Bader analysis for %d systems", len(systems))
        results = {}
        for name in systems:
            # 快速路径：直接加载 data/{name}/bader_analysis/bader_results.json
            cached = None if args.force else _load_cached_result(name)

            if cached is not None:
                results[name] = cached
                log.info("  ⏩ %s: loaded from cache", name)
            else:
                log.info("Processing %s ...", name)
                r = analyze_bader(name)
                results[name] = r
                if "error" not in r:
                    _save_cached_result(name, r)
    else:
        log.info("Loaded %d system results from output/bader.json", len(results))

    # Print summary table
    print("\n" + "=" * 100)
    print("BADER CHARGE ANALYSIS SUMMARY")
    print("=" * 100)
    print(f"  {'System':<28s}  {'Atoms':>6s}  {'CT mean':>10s}  {'CT |mean|':>10s}  {'CT min':>10s}  {'CT max':>10s}")
    print(f"  {'─'*28}  {'─'*6}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*10}")
    for name in SYSTEMS_ALL:
        r = results.get(name, {})
        if "error" in r:
            print(f"  {SYSTEM_LABELS.get(name, name):<28s}  ERROR: {r['error']}")
            continue
        ct = r.get("charge_transfer", [])
        if not ct:
            print(f"  {SYSTEM_LABELS.get(name, name):<28s}  (no data)")
            continue
        avg = np.mean(ct)
        abs_avg = np.mean(np.abs(ct))
        print(f"  {SYSTEM_LABELS.get(name, name):<28s}  {len(ct):>6d}  "
              f"{avg:>+10.4f}  {abs_avg:>10.4f}  {min(ct):>+10.4f}  {max(ct):>+10.4f}")

    # Data export
    log.info("Exporting data...")
    export_summary_csv(results)
    export_atoms_csv(results)
    export_results_json(results)

    log.info("Done. Output → %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
