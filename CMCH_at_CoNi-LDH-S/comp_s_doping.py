"""
S 掺杂 CoNi 氢氧化物 — 6 体系多维度对比分析

体系 (3×3×2 超胞):
  Group 1 — 本征基准:
    CoNiHO              : 本征态 (Co9 Ni9 H36 O36, 90 atoms)

  Group 2 — SH⁻ 掺杂 (S-H 保留, 90 atoms):
    CoNiHOS-Co3         : S 与 3 Co 配位
    CoNiHOS-Ni3         : S 与 3 Ni 配位
    CoNiHOS-Co1Ni2      : S 与 1 Co + 2 Ni 配位
    CoNiHOS-Co2Ni1      : S 与 2 Co + 1 Ni 配位

  Group 3 — S²⁻ 掺杂 (S 无 H, 89 atoms):
    CoNiHOS-Co3-noH     : S 与 3 Co 配位, 脱氢态

分析维度:
  A. 6 体系总览对比
  B. 本征 → SH-Co3 掺杂效应
  C. SH-Co3 → S-Co3-noH 脱氢效应
  D. 4 种 SH 配位环境对比
  E. Bader 电荷 + 磁性综合分析
"""

import argparse
import json
import logging
import warnings
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

from monty.json import MontyDecoder
from pymatgen.core import Structure
from pymatgen.electronic_structure.dos import CompleteDos, Spin
from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine
from pymatgen.electronic_structure.core import OrbitalType
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.io.vasp.inputs import Incar
from pymatgen.io.vasp.outputs import Outcar

warnings.filterwarnings("ignore")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s[%(levelname)s]%(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["font.size"] = 11
plt.rcParams["axes.linewidth"] = 1.2
plt.rcParams["figure.dpi"] = 100

# ═══════════════════════════════════════════════════════════════════════════
# 常量
# ═══════════════════════════════════════════════════════════════════════════

SYSTEMS_ALL = [
    "CoNiHO",
    "CoNiHOS-Co3",
    "CoNiHOS-Ni3",
    "CoNiHOS-Co1Ni2",
    "CoNiHOS-Co2Ni1",
    "CoNiHOS-Co3-noH",
]

SYSTEMS_DOPED_SH = ["CoNiHOS-Co3", "CoNiHOS-Ni3", "CoNiHOS-Co1Ni2", "CoNiHOS-Co2Ni1"]
SYSTEMS_DOPED_ALL = SYSTEMS_DOPED_SH + ["CoNiHOS-Co3-noH"]

SYSTEM_LABELS = {
    "CoNiHO":              "Pristine",
    "CoNiHOS-Co3":         "SH-3Co",
    "CoNiHOS-Ni3":         "SH-3Ni",
    "CoNiHOS-Co1Ni2":      "SH-Co1Ni2",
    "CoNiHOS-Co2Ni1":      "SH-Co2Ni1",
    "CoNiHOS-Co3-noH":     "S-3Co (no H)",
}

SYSTEM_COLORS = {
    "CoNiHO":              "#1f77b4",
    "CoNiHOS-Co3":         "#2ca02c",
    "CoNiHOS-Ni3":         "#d62728",
    "CoNiHOS-Co1Ni2":      "#ff7f0e",
    "CoNiHOS-Co2Ni1":      "#9467bd",
    "CoNiHOS-Co3-noH":     "#8c564b",
}

# ZVAL reference for charge transfer
ZVAL = {"Co": 9.0, "Ni": 10.0, "O": 6.0, "H": 1.0, "S": 6.0}

# Output directory (separate from input data/)
DATA_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/data")
OUTPUT_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/doping-output")


# ═══════════════════════════════════════════════════════════════════════════
# 数据加载
# ═══════════════════════════════════════════════════════════════════════════

def load_dos(dos_path: Path) -> CompleteDos | None:
    if not dos_path.exists():
        log.warning("DOS file not found: %s", dos_path)
        return None
    try:
        with open(dos_path, encoding="utf-8") as f:
            return json.load(f, cls=MontyDecoder).vasp_objects["dos"]
    except Exception as e:
        log.warning("Failed to load DOS from %s: %s", dos_path, e)
        return None


def load_band(band_path: Path) -> BandStructure | BandStructureSymmLine | None:
    if not band_path.exists():
        log.warning("Band file not found: %s", band_path)
        return None
    try:
        with open(band_path, encoding="utf-8") as f:
            return json.load(f, cls=MontyDecoder).vasp_objects["bandstructure"]
    except Exception as e:
        log.warning("Failed to load band from %s: %s", band_path, e)
        return None


def load_store(store_path: Path) -> list[dict]:
    with open(store_path, encoding="utf-8") as f:
        return json.load(f)


def load_structure(contcar_path: Path) -> Structure | None:
    if not contcar_path.exists():
        log.warning("CONTCAR not found: %s", contcar_path)
        return None
    try:
        return Structure.from_file(contcar_path)
    except Exception as e:
        log.warning("Failed to load structure from %s: %s", contcar_path, e)
        return None


def load_bader(bader_path: Path) -> dict | None:
    if not bader_path.exists():
        log.warning("Bader file not found: %s", bader_path)
        return None
    try:
        with open(bader_path, encoding="utf-8") as f:
            return json.load(f)
    except Exception as e:
        log.warning("Failed to load Bader from %s: %s", bader_path, e)
        return None


# ═══════════════════════════════════════════════════════════════════════════
# 信息提取
# ═══════════════════════════════════════════════════════════════════════════

def extract_static_info(store_data: list[dict]) -> dict:
    """从 store.json 提取 static 步骤的基本信息."""
    for job in store_data:
        if job.get("name") == "static":
            output = job.get("output", {})
            inner = output.get("output", {})
            return {
                "energy": inner.get("energy"),
                "energy_per_atom": inner.get("energy_per_atom"),
                "bandgap": inner.get("bandgap"),
                "nsites": output.get("nsites"),
                "formula": output.get("formula_pretty"),
                "state": output.get("state"),
            }
    return {}


def extract_mag_from_outcar(outcar_path: Path) -> tuple[float | None, list[float]]:
    try:
        outcar = Outcar(outcar_path)
        total = outcar.total_mag
        atomic = [m["tot"] for m in outcar.magnetization] if outcar.magnetization else []
        return total, atomic
    except Exception:
        return None, []


def s_metal_bonds(structure: Structure) -> dict:
    """S 与最近邻金属的键长信息."""
    s_idx = next((i for i, s in enumerate(structure) if s.species_string == "S"), None)
    if s_idx is None:
        return {"has_S": False}

    s_site = structure[s_idx]
    neighbors = structure.get_neighbors(s_site, r=3.5)
    metal_nn = [
        (n.index, n.species_string, n.nn_distance)
        for n in neighbors if n.species_string in ("Co", "Ni")
    ]
    metal_nn.sort(key=lambda x: x[2])
    top3 = metal_nn[:3]

    # Also find nearest H for SH vs S distinction
    h_nn = [
        (n.index, n.nn_distance)
        for n in neighbors if n.species_string == "H"
    ]
    h_nn.sort(key=lambda x: x[1])

    return {
        "has_S": True,
        "s_index": s_idx,
        "nearest_metals": [(x[0], x[1]) for x in top3],
        "metal_types": [x[1] for x in top3],
        "bond_lengths": [x[2] for x in top3],
        "avg_bond_length": np.mean([x[2] for x in top3]) if top3 else None,
        "nearest_H_dist": h_nn[0][1] if h_nn else None,
        "has_SH_bond": h_nn and h_nn[0][1] < 1.5,
    }


def d_band_center(dos: CompleteDos, element: str, erange: tuple = (-10, 5)) -> float | None:
    try:
        return dos.get_band_center(band=OrbitalType.d, elements=[element], erange=erange)
    except Exception:
        return None


def pdos_at_fermi(dos: CompleteDos, efermi: float, window: float = 0.3) -> dict:
    energies = dos.energies
    mask = np.abs(energies - efermi) <= window
    result = {"S_p": 0.0, "Co_d": 0.0, "Ni_d": 0.0, "total": 0.0}

    for site in dos.structure:
        symbol = site.species_string
        if symbol not in ("S", "Co", "Ni"):
            continue
        spd = dos.get_site_spd_dos(site)
        if symbol == "S" and OrbitalType.p in spd:
            for spin in (Spin.up, Spin.down):
                if spin in spd[OrbitalType.p].densities:
                    result["S_p"] += np.mean(spd[OrbitalType.p].densities[spin][mask])
        elif symbol in ("Co", "Ni") and OrbitalType.d in spd:
            key = f"{symbol}_d"
            for spin in (Spin.up, Spin.down):
                if spin in spd[OrbitalType.d].densities:
                    result[key] += np.mean(spd[OrbitalType.d].densities[spin][mask])

    total = dos.get_densities(Spin.up)
    if Spin.down in dos.densities:
        total += dos.get_densities(Spin.down)
    result["total"] = np.mean(total[mask])
    return result


def spin_polarization(dos: CompleteDos, e_window: float = 0.05, threshold: float = 1e-3) -> float | None:
    """计算费米能级处自旋极化率, 带阈值保护.

    仅当 E_F ± e_window 窗口内总 DOS 均值 > threshold states/eV 时计算;
    半导体/绝缘体返回 None (输出中显示为 "—").

    Args:
        dos: CompleteDos 对象
        e_window: 费米能级附近能量窗口 (eV), 默认 ±0.05
        threshold: 总 DOS 阈值 (states/eV), 默认 1e-3

    Returns:
        float | None: 自旋极化率 [0, 1], 或 None 表示不适用 (绝缘体/半导体).
    """
    try:
        energies = dos.energies - dos.efermi
        mask = np.abs(energies) <= e_window

        dos_up = dos.get_densities(Spin.up)[mask]
        if Spin.down in dos.densities:
            dos_down = dos.get_densities(Spin.down)[mask]
        else:
            return None  # 非自旋极化计算

        total_dos = np.mean(dos_up + dos_down)

        if total_dos < threshold:
            return None  # 半导体/绝缘体: DOS 可忽略

        # 用插值点值计算自旋极化率
        n_F = dos.get_interpolated_value(dos.efermi)
        n_F_up = n_F[Spin.up]
        if Spin.down not in n_F:
            return None
        n_F_down = n_F[Spin.down]
        if (n_F_up + n_F_down) == 0:
            return None
        return abs((n_F_up - n_F_down) / (n_F_up + n_F_down))
    except Exception:
        return None


def vbm_cbm(band) -> dict:
    result = {"vbm": None, "cbm": None}
    try:
        result["vbm"] = band.get_vbm()["energy"]
        result["cbm"] = band.get_cbm()["energy"]
    except Exception:
        pass
    return result


def lattice_info(structure: Structure) -> dict:
    lat = structure.lattice
    return {"a": lat.a, "b": lat.b, "c": lat.c, "volume": lat.volume}


# ═══════════════════════════════════════════════════════════════════════════
# Bader 分析
# ═══════════════════════════════════════════════════════════════════════════

def _normalize_bader(bader_data: dict, structure: Structure) -> dict:
    """Handle both flat (bader_analysis_from_path) and atoms (bader_charge_analysis.py) formats.

    统一电荷转移符号约定: charge_transfer = ZVAL - bader_charge
    正值表示原子失去电子, 负值表示获得电子 (学术惯例).
    """
    if "charge" in bader_data:
        # Flat format — negate pre-computed charge_transfer to academic convention
        if "charge_transfer" in bader_data:
            bader_data = dict(bader_data)  # shallow copy to avoid mutating original
            bader_data["charge_transfer"] = [-ct for ct in bader_data["charge_transfer"]]
        return bader_data
    if "atoms" in bader_data:
        atoms = bader_data["atoms"]
        n = len(structure)
        return {
            "charge": [a["bader_charge"] for a in atoms],
            # Negate: original stored as bader_charge - ZVAL, we want ZVAL - bader_charge
            "charge_transfer": [-a["charge_transfer"] for a in atoms],
            "magmom": [a["magmom"] for a in atoms],
            "atomic_volume": [a["atomic_vol"] for a in atoms],
            "min_dist": [a["min_dist"] for a in atoms],
            "vacuum_charge": bader_data.get("vacuum_charge", 0),
            "vacuum_volume": bader_data.get("vacuum_volume", 0),
            "reference_used": bader_data.get("reference_used", True),
            "bader_version": bader_data.get("bader_version", -1),
        }
    return bader_data


def analyze_bader_local(bader_data: dict, structure: Structure, name: str) -> dict:
    """
    提取 S 局域 Bader 信息: S 电荷、S 周围金属的电荷和磁矩.
    """
    if not bader_data:
        return {}

    bader_data = _normalize_bader(bader_data, structure)
    charges = bader_data.get("charge", [])
    ct = bader_data.get("charge_transfer", [])
    magmoms = bader_data.get("magmom", [])
    vols = bader_data.get("atomic_volume", [])

    s_idx = next((i for i, s in enumerate(structure) if s.species_string == "S"), None)
    if s_idx is None:
        return {"has_S": False}

    s_bond = s_metal_bonds(structure)
    nn_sites = [x[0] for x in s_bond.get("nearest_metals", [])]

    def _safe(lst, i):
        return lst[i] if i < len(lst) else None

    result = {
        "has_S": True,
        "s_index": s_idx,
        "S_bader_charge": _safe(charges, s_idx),
        "S_charge_transfer": _safe(ct, s_idx),
        "S_magmom": _safe(magmoms, s_idx),
        "S_atomic_vol": _safe(vols, s_idx),
        "nn_metal_bader": [],
    }

    for idx in nn_sites:
        elem = str(structure[idx].species)
        result["nn_metal_bader"].append({
            "index": idx,
            "element": elem,
            "bader_charge": _safe(charges, idx),
            "charge_transfer": _safe(ct, idx),
            "magmom": _safe(magmoms, idx),
            "atomic_vol": _safe(vols, idx),
        })

    # Aggregates
    nn_ct = [x["charge_transfer"] for x in result["nn_metal_bader"] if x["charge_transfer"] is not None]
    nn_mag = [x["magmom"] for x in result["nn_metal_bader"] if x["magmom"] is not None]
    result["nn_avg_transfer"] = np.mean(nn_ct) if nn_ct else None
    result["nn_avg_magmom"] = np.mean(nn_mag) if nn_mag else None

    # All Co Bader charges and magmoms
    co_indices = [i for i, s in enumerate(structure) if s.species_string == "Co"]
    ni_indices = [i for i, s in enumerate(structure) if s.species_string == "Ni"]
    result["Co_bader_charges"] = [_safe(charges, i) for i in co_indices]
    result["Co_magmoms"] = [_safe(magmoms, i) for i in co_indices]
    result["Ni_bader_charges"] = [_safe(charges, i) for i in ni_indices]
    result["Ni_magmoms"] = [_safe(magmoms, i) for i in ni_indices]
    result["Co_avg_bader_charge"] = np.mean([c for c in result["Co_bader_charges"] if c is not None])
    result["Co_avg_magmom"] = np.mean([m for m in result["Co_magmoms"] if m is not None])
    result["Ni_avg_bader_charge"] = np.mean([c for c in result["Ni_bader_charges"] if c is not None])
    result["Ni_avg_magmom"] = np.mean([m for m in result["Ni_magmoms"] if m is not None])

    return result


# ═══════════════════════════════════════════════════════════════════════════
# 系统分析主函数
# ═══════════════════════════════════════════════════════════════════════════

def analyze_system(data_dir: Path, name: str) -> dict:
    log.info("Analyzing %s ...", name)
    result = {"name": name}

    # Pre-check: store.json must exist
    sp = data_dir / "store.json"
    if not sp.exists():
        log.warning("store.json not found for %s, skipping", name)
        return result
    result.update(extract_static_info(load_store(sp)))

    # DOS
    dp = data_dir / "dos_out.json"
    if dp.exists():
        dos = load_dos(dp)
        if dos:
            result["dos"] = dos
            result["efermi"] = dos.efermi
            result["bandgap_dos"] = dos.get_gap()
            result["d_center_Co"] = d_band_center(dos, "Co")
            result["d_center_Ni"] = d_band_center(dos, "Ni")
            result["pdos_ef"] = pdos_at_fermi(dos, dos.efermi)
            result["spin_pol"] = spin_polarization(dos)

    # Band
    bp = data_dir / "band_out.json"
    if bp.exists():
        band = load_band(bp)
        if band:
            result["band"] = band
            result["is_line_mode"] = isinstance(band, BandStructureSymmLine)
            result.update(vbm_cbm(band))

    # Structure
    contcar = data_dir / "3-static" / "CONTCAR.gz"
    if contcar.exists():
        structure = load_structure(contcar)
        if structure:
            result["structure"] = structure
            result["lattice"] = lattice_info(structure)
            result["S_bonds"] = s_metal_bonds(structure)

    # Magnetic moments (OUTCAR)
    outcar = data_dir / "3-static" / "OUTCAR.gz"
    if outcar.exists():
        total_mag, atomic_mags = extract_mag_from_outcar(outcar)
        result["total_mag"] = total_mag
        result["atomic_mags"] = atomic_mags
        if atomic_mags and "structure" in result:
            struct = result["structure"]
            co_i = [i for i, s in enumerate(struct) if s.species_string == "Co"]
            ni_i = [i for i, s in enumerate(struct) if s.species_string == "Ni"]
            result["Co_mags"] = [atomic_mags[i] for i in co_i if i < len(atomic_mags)]
            result["Ni_mags"] = [atomic_mags[i] for i in ni_i if i < len(atomic_mags)]
            result["Co_avg_mag"] = np.mean(result["Co_mags"]) if result["Co_mags"] else None
            result["Ni_avg_mag"] = np.mean(result["Ni_mags"]) if result["Ni_mags"] else None

    # Bader
    bf = data_dir / "bader_analysis" / "bader_results.json"
    if bf.exists():
        bader_data = load_bader(bf)
        result["bader"] = bader_data
        if bader_data and "structure" in result:
            result["bader_local"] = analyze_bader_local(bader_data, result["structure"], name)

    return result


# ═══════════════════════════════════════════════════════════════════════════
# A ─ 6 体系总览对比表格
# ═══════════════════════════════════════════════════════════════════════════

def overview_table(results: dict) -> pd.DataFrame:
    rows = []

    specs = [
        ("Formula",                    "formula",           "{}"),
        ("N atoms",                    "nsites",            "{}"),
        ("Total Energy (eV)",          "energy",            "{:.4f}"),
        ("Energy/atom (eV)",           "energy_per_atom",   "{:.4f}"),
        ("Band Gap (eV)",              "bandgap_dos",       "{:.4f}"),
        ("VBM (eV)",                   "vbm",               "{:.4f}"),
        ("CBM (eV)",                   "cbm",               "{:.4f}"),
        ("E_Fermi (eV)",               "efermi",            "{:.4f}"),
        ("Total Mag (μB)",             "total_mag",         "{:.4f}"),
        ("Co avg Mag (μB)",            "Co_avg_mag",        "{:.4f}"),
        ("Ni avg Mag (μB)",            "Ni_avg_mag",        "{:.4f}"),
        ("Spin Polarization",          "spin_pol",          "{:.4f}"),
        ("d-center Co (eV)",           "d_center_Co",       "{:.4f}"),
        ("d-center Ni (eV)",           "d_center_Ni",       "{:.4f}"),
        ("Volume (Å³)",                None,                None),
        ("S-M avg bond (Å)",           None,                None),
        ("S Coordination",             None,                None),
        ("S-H bond?",                  None,                None),
        ("S Bader charge (e)",         None,                None),
        ("S charge transfer (e)",      None,                None),
        ("S magmom (μB)",              None,                None),
        ("NN metals avg transfer (e)",  None,                None),
        ("NN metals avg magmom (μB)",   None,                None),
    ]

    def _fmt(val, fmt_str):
        if val is None:
            return "—"
        try:
            return fmt_str.format(val)
        except (ValueError, TypeError):
            return str(val)

    for label, key, fmt in specs:
        row = {"Property": label}
        for name in SYSTEMS_ALL:
            r = results[name]
            if label == "S Coordination":
                sb = r.get("S_bonds", {})
                row[name] = "/".join(sb.get("metal_types", [])) if sb.get("has_S") else "—"
            elif label == "S-M avg bond (Å)":
                sb = r.get("S_bonds", {})
                row[name] = _fmt(sb.get("avg_bond_length"), "{:.4f}")
            elif label == "S-H bond?":
                sb = r.get("S_bonds", {})
                dist = sb.get("nearest_H_dist")
                row[name] = f"✓ ({dist:.2f} Å)" if sb.get("has_SH_bond") else f"✗ (nearest {dist:.2f} Å)" if dist else "—"
            elif label == "Volume (Å³)":
                lat = r.get("lattice", {})
                row[name] = _fmt(lat.get("volume"), "{:.2f}")
            elif label in ("S Bader charge (e)", "S charge transfer (e)", "S magmom (μB)",
                           "NN metals avg transfer (e)", "NN metals avg magmom (μB)"):
                bl = r.get("bader_local", {})
                field_map = {
                    "S Bader charge (e)": "S_bader_charge",
                    "S charge transfer (e)": "S_charge_transfer",
                    "S magmom (μB)": "S_magmom",
                    "NN metals avg transfer (e)": "nn_avg_transfer",
                    "NN metals avg magmom (μB)": "nn_avg_magmom",
                }
                fk = field_map.get(label, "")
                val = bl.get(fk)
                row[name] = _fmt(val, "{:.4f}")
            elif key and fmt:
                val = r.get(key)
                row[name] = _fmt(val, fmt)
            else:
                row[name] = "—"
        rows.append(row)

    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════════════
# B ─ 本征 → SH-Co3 掺杂效应 (Pristine vs SH-3Co)
# ═══════════════════════════════════════════════════════════════════════════

def analyze_doping_effect(results: dict):
    """本征 vs SH-Co3 对比: S 替代 OH 的影响."""
    p = results["CoNiHO"]
    d = results["CoNiHOS-Co3"]

    print("\n" + "=" * 80)
    print("B. DOPING EFFECT: Pristine → SH-3Co")
    print("=" * 80)

    # Energy
    e_p, e_d = p.get("energy"), d.get("energy")
    if e_p and e_d:
        delta_e = e_d - e_p
        print(f"\n  Substitution energy (no chemical potentials):")
        print(f"    ΔE = E(SH-3Co) - E(Pristine) = {delta_e:+.4f} eV  ({delta_e/90*1000:+.2f} meV/atom)")

    # Structure
    lat_p, lat_d = p.get("lattice", {}), d.get("lattice", {})
    if lat_p and lat_d:
        print(f"\n  Lattice distortion:")
        print(f"    {'':>10s}  {'Pristine':>12s}  {'SH-3Co':>12s}  {'Δ':>10s}")
        for axis in ("a", "c", "volume"):
            vp, vd = lat_p.get(axis), lat_d.get(axis)
            if vp and vd:
                print(f"    {axis:>10s}  {vp:12.4f}  {vd:12.4f}  {vd-vp:+10.4f}")

    # Band gap
    gp, gd = p.get("bandgap_dos"), d.get("bandgap_dos")
    if gp is not None and gd is not None:
        print(f"\n  Band gap: {gp:.4f} → {gd:.4f} eV  (Δ = {gd-gp:+.4f} eV)")

    # Magnetism
    print(f"\n  Magnetic moments:")
    print(f"    {'':>10s}  {'Pristine':>12s}  {'SH-3Co':>12s}")
    print(f"    {'Total':>10s}  {p.get('total_mag', 0):12.4f}  {d.get('total_mag', 0):12.4f}")
    print(f"    {'Co avg':>10s}  {p.get('Co_avg_mag', 0):12.4f}  {d.get('Co_avg_mag', 0):12.4f}")
    print(f"    {'Ni avg':>10s}  {p.get('Ni_avg_mag', 0):12.4f}  {d.get('Ni_avg_mag', 0):12.4f}")

    # S-M bond lengths
    sb = d.get("S_bonds", {})
    if sb.get("has_S"):
        print(f"\n  S-M bond lengths in SH-3Co:")
        for idx, elem, bl in zip(sb["nearest_metals"], sb["metal_types"], sb["bond_lengths"]):
            print(f"    {elem} site {idx[0]}: {bl:.4f} Å")

    # Bader — 3 Co around S in doped vs corresponding Co in pristine
    bl_d = d.get("bader_local", {})
    print(f"\n  Bader charge comparison (S-local Co atoms):")
    nn_list = sb.get("nearest_metals", [])
    if bl_d and nn_list:
        print(f"    {'Site':>10s}  {'Element':>8s}  {'SH-3Co CT':>12s}  {'Pristine mag':>14s}")
        p_mags = p.get("atomic_mags", [])
        for nn in nn_list:
            idx, elem = nn
            d_ct = next((x["charge_transfer"] for x in bl_d["nn_metal_bader"] if x["index"] == idx), None)
            p_mag = p_mags[idx] if idx < len(p_mags) else None
            ct_s = f"{d_ct:+.4f} e" if d_ct is not None else "N/A"
            pm_s = f"{p_mag:+.4f} μB" if p_mag is not None else "N/A"
            print(f"    {idx:>10d}  {elem:>8s}  {ct_s:>12s}  {pm_s:>14s}")


# ═══════════════════════════════════════════════════════════════════════════
# C ─ SH-Co3 → S-Co3-noH 脱氢效应
# ═══════════════════════════════════════════════════════════════════════════

def analyze_dehydrogenation(results: dict):
    """SH-3Co vs S-3Co-noH 对比: H 脱除的影响."""
    sh = results["CoNiHOS-Co3"]
    s0 = results["CoNiHOS-Co3-noH"]

    print("\n" + "=" * 80)
    print("C. DEHYDROGENATION EFFECT: SH-3Co → S-3Co (no H)")
    print("=" * 80)

    # Dehydrogenation energy
    e_sh, e_s = sh.get("energy"), s0.get("energy")
    if e_sh and e_s:
        # E(S-Co3) - E(SH-Co3) = energy cost to remove H
        de = e_s - e_sh
        print(f"\n  Raw ΔE (S-Co3 − SH-Co3) = {de:+.4f} eV")

    # Structure
    sb_sh = sh.get("S_bonds", {})
    sb_s = s0.get("S_bonds", {})
    if sb_sh.get("has_S") and sb_s.get("has_S"):
        print(f"\n  S-M bond length change after H removal:")
        bl_sh = sb_sh["bond_lengths"]
        bl_s = sb_s["bond_lengths"]
        metals = sb_sh["metal_types"]
        for i in range(3):
            print(f"    {metals[i]}-S: {bl_sh[i]:.4f} → {bl_s[i]:.4f} Å  (Δ = {bl_s[i]-bl_sh[i]:+.4f} Å)")

    # Band gap
    gp_sh, gp_s = sh.get("bandgap_dos"), s0.get("bandgap_dos")
    if gp_sh is not None and gp_s is not None:
        print(f"\n  Band gap: {gp_sh:.4f} → {gp_s:.4f} eV  (Δ = {gp_s-gp_sh:+.4f} eV)")

    # PDOS at Ef
    print(f"\n  PDOS at E_F (±0.3 eV):")
    pef_sh = sh.get("pdos_ef", {})
    pef_s = s0.get("pdos_ef", {})
    if pef_sh and pef_s:
        total_sh = pef_sh.get("total", 1)
        total_s = pef_s.get("total", 1)
        print(f"    {'':>12s}  {'SH-3Co':>12s}  {'S-3Co':>12s}")
        for k in ("S_p", "Co_d", "Ni_d", "total"):
            v_sh = pef_sh.get(k, 0) * 100 / total_sh if total_sh else 0
            v_s = pef_s.get(k, 0) * 100 / total_s if total_s else 0
            print(f"    {k:>12s}  {pef_sh.get(k,0):12.4f}  {pef_s.get(k,0):12.4f}")

    # Bader local: S and 3 Co
    bl_sh = sh.get("bader_local", {})
    bl_s = s0.get("bader_local", {})
    if bl_sh and bl_s:
        print(f"\n  Bader charge transfer change:")
        print(f"    {'Atom':>12s}  {'SH-3Co':>10s}  {'S-3Co':>10s}  {'Δ':>10s}")
        s_ct_sh = bl_sh.get("S_charge_transfer")
        s_ct_s = bl_s.get("S_charge_transfer")
        if s_ct_sh and s_ct_s:
            print(f"    {'S':>12s}  {s_ct_sh:+10.4f}  {s_ct_s:+10.4f}  {s_ct_s-s_ct_sh:+10.4f}")
        for nb_sh, nb_s in zip(bl_sh.get("nn_metal_bader", []), bl_s.get("nn_metal_bader", [])):
            ct_sh, ct_s = nb_sh.get("charge_transfer"), nb_s.get("charge_transfer")
            if ct_sh is not None and ct_s is not None:
                label = f"{nb_sh['element']}({nb_sh['index']})"
                print(f"    {label:>12s}  {ct_sh:+10.4f}  {ct_s:+10.4f}  {ct_s-ct_sh:+10.4f}")

        print(f"\n  Magnetic moment change (S-local):")
        print(f"    {'Atom':>12s}  {'SH-3Co':>10s}  {'S-3Co':>10s}  {'Δ':>10s}")
        s_mm_sh = bl_sh.get("S_magmom")
        s_mm_s = bl_s.get("S_magmom")
        if s_mm_sh is not None and s_mm_s is not None:
            print(f"    {'S':>12s}  {s_mm_sh:+10.4f}  {s_mm_s:+10.4f}  {s_mm_s-s_mm_sh:+10.4f}")
        for nb_sh, nb_s in zip(bl_sh.get("nn_metal_bader", []), bl_s.get("nn_metal_bader", [])):
            mm_sh, mm_s = nb_sh.get("magmom"), nb_s.get("magmom")
            if mm_sh is not None and mm_s is not None:
                label = f"{nb_sh['element']}({nb_sh['index']})"
                print(f"    {label:>12s}  {mm_sh:+10.4f}  {mm_s:+10.4f}  {mm_s-mm_sh:+10.4f}")

    # Spin polarization
    sp_sh, sp_s = sh.get("spin_pol"), s0.get("spin_pol")
    if sp_sh is not None and sp_s is not None:
        print(f"\n  Spin polarization at E_F: {sp_sh:.4f} → {sp_s:.4f}")

    # d-band centers
    print(f"\n  d-band center shift:")
    dco_sh, dco_s = sh.get("d_center_Co"), s0.get("d_center_Co")
    dni_sh, dni_s = sh.get("d_center_Ni"), s0.get("d_center_Ni")
    if dco_sh and dco_s:
        print(f"    Co d-center: {dco_sh:.4f} → {dco_s:.4f} eV  (Δ = {dco_s-dco_sh:+.4f})")
    if dni_sh and dni_s:
        print(f"    Ni d-center: {dni_sh:.4f} → {dni_s:.4f} eV  (Δ = {dni_s-dni_sh:+.4f})")


# ═══════════════════════════════════════════════════════════════════════════
# D ─ 4 种 SH 配位环境对比
# ═══════════════════════════════════════════════════════════════════════════

def analyze_coordination(results: dict):
    """4 种 SH 掺杂配位环境全面对比."""
    print("\n" + "=" * 80)
    print("D. COORDINATION EFFECT: 4 SH-Doped Systems")
    print("=" * 80)

    names = SYSTEMS_DOPED_SH

    # Formation energies (raw, no μ)
    e_pristine = results["CoNiHO"].get("energy")
    print(f"\n  Formation energies (ΔE = E_doped - E_pristine, same atom count):")
    print(f"    {'System':<20s}  {'ΔE (eV)':>12s}  {'ΔE (meV/atom)':>15s}")
    for name in names:
        e_d = results[name].get("energy")
        if e_pristine and e_d:
            de = e_d - e_pristine
            print(f"    {name:<20s}  {de:+12.4f}  {de/90*1000:+15.2f}")

    # S-M bond lengths
    print(f"\n  S-M bond lengths:")
    print(f"    {'System':<20s}  {'Coordination':<18s}  {'Bond1 (Å)':>10s}  {'Bond2 (Å)':>10s}  {'Bond3 (Å)':>10s}  {'Avg (Å)':>10s}")
    for name in names:
        sb = results[name].get("S_bonds", {})
        if sb.get("has_S"):
            bl = sb["bond_lengths"]
            mt = sb["metal_types"]
            print(f"    {name:<20s}  {str(mt):<18s}  {bl[0]:10.4f}  {bl[1]:10.4f}  {bl[2]:10.4f}  {sb['avg_bond_length']:10.4f}")

    # Band gaps
    print(f"\n  Band gaps:")
    for name in names:
        gap = results[name].get("bandgap_dos")
        print(f"    {name:<20s}  {gap:.4f} eV" if gap is not None else f"    {name:<20s}  N/A")

    # Magnetism — S-local 3 metal atoms
    print(f"\n  Magnetic moments — S-neighboring metal atoms:")
    print(f"    {'System':<20s}  {'M1 (μB)':>10s}  {'M2 (μB)':>10s}  {'M3 (μB)':>10s}  {'Avg (μB)':>10s}  {'S (μB)':>10s}")
    for name in names:
        bl = results[name].get("bader_local", {})
        nn = bl.get("nn_metal_bader", [])
        mags = [x["magmom"] for x in nn if x["magmom"] is not None]
        s_mag = bl.get("S_magmom")
        if len(mags) >= 3:
            s_str = f"{s_mag:10.4f}" if s_mag is not None else "       N/A"
            print(f"    {name:<20s}  {mags[0]:10.4f}  {mags[1]:10.4f}  {mags[2]:10.4f}  {np.mean(mags):10.4f}  {s_str}")

    # Bader charge transfer — S and neighbors
    print(f"\n  Bader charge transfer — S and nearest metals:")
    print(f"    {'System':<20s}  {'S (e)':>10s}  {'M1 (e)':>10s}  {'M2 (e)':>10s}  {'M3 (e)':>10s}  {'M-avg (e)':>10s}")
    for name in names:
        bl = results[name].get("bader_local", {})
        nn = bl.get("nn_metal_bader", [])
        cts = [x["charge_transfer"] for x in nn if x["charge_transfer"] is not None]
        s_ct = bl.get("S_charge_transfer")
        if len(cts) >= 3 and s_ct is not None:
            print(f"    {name:<20s}  {s_ct:+10.4f}  {cts[0]:+10.4f}  {cts[1]:+10.4f}  {cts[2]:+10.4f}  {np.mean(cts):+10.4f}")

    # d-band centers
    print(f"\n  d-band centers (relative to E_F):")
    print(f"    {'System':<20s}  {'Co d-center':>14s}  {'Ni d-center':>14s}")
    for name in names:
        d_co = results[name].get("d_center_Co")
        d_ni = results[name].get("d_center_Ni")
        co_s = f"{d_co:.4f} eV" if d_co is not None else "N/A"
        ni_s = f"{d_ni:.4f} eV" if d_ni is not None else "N/A"
        print(f"    {name:<20s}  {co_s:>14s}  {ni_s:>14s}")

    # PDOS at Ef
    print(f"\n  PDOS contribution at E_F (±0.3 eV):")
    print(f"    {'System':<20s}  {'S-3p':>10s}  {'Co-3d':>10s}  {'Ni-3d':>10s}  {'Total':>10s}")
    for name in names:
        pef = results[name].get("pdos_ef", {})
        total = pef.get("total", 0) or 1
        s_p = pef.get("S_p", 0) / total * 100
        co_d = pef.get("Co_d", 0) / total * 100
        ni_d = pef.get("Ni_d", 0) / total * 100
        print(f"    {name:<20s}  {s_p:9.1f}%  {co_d:9.1f}%  {ni_d:9.1f}%  {total:10.4f}")


# ═══════════════════════════════════════════════════════════════════════════
# E ─ Bader 电荷 + 磁性综合分析 (全 6 体系)
# ═══════════════════════════════════════════════════════════════════════════

def analyze_bader_comprehensive(results: dict):
    """全 6 体系 Bader 电荷和磁性总结."""
    print("\n" + "=" * 80)
    print("E. BADER CHARGE + MAGNETIC COMPREHENSIVE ANALYSIS")
    print("=" * 80)
    print("  注: 电荷转移 = ZVAL - bader_charge, 正值 = 失去电子, 负值 = 获得电子")

    # S charge state across all doped systems
    print(f"\n  S atom Bader analysis:")
    print(f"    {'System':<22s}  {'Bader(e)':>10s}  {'Transfer(e)':>12s}  {'Vol(Å³)':>10s}  {'Mag(μB)':>10s}  {'Nearest H':>12s}")
    for name in SYSTEMS_DOPED_ALL:
        bl = results[name].get("bader_local", {})
        sb = results[name].get("S_bonds", {})
        if bl.get("has_S"):
            bc = bl.get("S_bader_charge")
            ct = bl.get("S_charge_transfer")
            vol = bl.get("S_atomic_vol")
            mm = bl.get("S_magmom")
            nh = sb.get("nearest_H_dist", None)
            nh_s = f"{nh:.2f} Å" if nh else "—"
            bc_s = f"{bc:10.4f}" if bc is not None else "        —"
            ct_s = f"{ct:+12.4f}" if ct is not None else "          —"
            vol_s = f"{vol:10.2f}" if vol is not None else "        —"
            mm_s = f"{mm:10.4f}" if mm is not None else "        —"
            print(f"    {name:<22s}  {bc_s}  {ct_s}  {vol_s}  {mm_s}  {nh_s:>12s}")

    # Co valence comparison: For the 3 Co nearest to S
    print(f"\n  Co neighbors around S — detailed Bader + magnetic:")
    for name in SYSTEMS_DOPED_ALL:
        bl = results[name].get("bader_local", {})
        if not bl.get("has_S"):
            continue
        print(f"\n  ── {name} ──")
        print(f"    {'Index':>6s}  {'Elem':>6s}  {'Bader(e)':>10s}  {'Transfer(e)':>12s}  {'Mag(μB)':>10s}  {'Vol(Å³)':>10s}")
        for nb in bl.get("nn_metal_bader", []):
            print(f"    {nb['index']:6d}  {nb['element']:>6s}  {nb['bader_charge']:10.4f}  "
                  f"{nb['charge_transfer']:+12.4f}  {nb['magmom']:10.4f}  {nb['atomic_vol']:10.2f}")

    # All Co average Bader charge across systems
    print(f"\n  Element-average Bader charge transfer:")
    print(f"    {'System':<22s}  {'Co avg CT':>12s}  {'Ni avg CT':>12s}  {'Co avg Mag':>12s}  {'Ni avg Mag':>12s}")
    for name in SYSTEMS_ALL:
        bl = results[name].get("bader_local", {})
        co_ct = bl.get("Co_avg_bader_charge")
        ni_ct = bl.get("Ni_avg_bader_charge")
        co_mm = bl.get("Co_avg_magmom")
        ni_mm = bl.get("Ni_avg_magmom")
        # Convert bader charge to transfer (学术惯例: ZVAL - bader_charge, 正值 = 失电子)
        co_transfer = (ZVAL["Co"] - co_ct) if co_ct is not None else None
        ni_transfer = (ZVAL["Ni"] - ni_ct) if ni_ct is not None else None
        co_s = f"{co_transfer:+12.4f}" if co_transfer is not None else "        N/A"
        ni_s = f"{ni_transfer:+12.4f}" if ni_transfer is not None else "        N/A"
        co_m = f"{co_mm:12.4f}" if co_mm is not None else "        N/A"
        ni_m = f"{ni_mm:12.4f}" if ni_mm is not None else "        N/A"
        print(f"    {name:<22s}  {co_s}  {ni_s}  {co_m}  {ni_m}")


# ═══════════════════════════════════════════════════════════════════════════
# 数据导出工具
# ═══════════════════════════════════════════════════════════════════════════

def save_plot_data(data: dict, output: Path):
    """Save companion JSON data alongside a plot figure.

    Args:
        data: Dict with keys 'plot_type', 'parameters', 'systems' (list), 'data'.
        output: Plot output path (without extension). Companion saved as {output}_data.json.
    """
    json_path = output.with_name(output.name + "_data.json").with_suffix("")
    if json_path.suffix != ".json":
        json_path = json_path.with_suffix(".json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=_json_default)
    log.info("Plot data → %s", json_path)


def _json_default(obj):
    """JSON serializer for numpy arrays and other non-standard types."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


# ═══════════════════════════════════════════════════════════════════════════
# 绘图工具
# ═══════════════════════════════════════════════════════════════════════════

def _save_figure(fig, output: Path, label: str):
    """Save figure as PNG + PDF, then close."""
    for fmt in ("png", "pdf"):
        out_path = output.with_suffix(f".{fmt}")
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        log.info("%s → %s", label, out_path)
    plt.close(fig)


def _annotate_bandgap(ax, gap, orientation: str = "v"):
    """Shade bandgap region on DOS (v) or band (h) plot."""
    if gap is None or gap <= 0.05:
        return
    if orientation == "v":
        ax.axvspan(0, gap, color="gray", alpha=0.1, lw=0)
        ax.axvline(0, color="gray", ls=":", lw=0.5)
        ax.axvline(gap, color="gray", ls=":", lw=0.5)
    else:
        ax.axhspan(0, gap, color="gray", alpha=0.1, lw=0)
        ax.axhline(0, color="gray", ls=":", lw=0.5)
        ax.axhline(gap, color="gray", ls=":", lw=0.5)


def _plot_scatter_panel(ax, x_vals, y_vals, labels, colors, xlabel, ylabel, title,
                         fit=None, note=None):
    """Plot a single scatter panel with annotations, optional fit line, and R²."""
    ax.scatter(x_vals, y_vals, c=colors, s=60, zorder=5, edgecolors="k", linewidth=0.5)
    for xi, yi, nm in zip(x_vals, y_vals, labels):
        ax.annotate(nm, (xi, yi), textcoords="offset points", xytext=(5, 5),
                    fontsize=7, alpha=0.9,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8, ec="none"))
    if fit and fit.get("slope") is not None and fit.get("r2") is not None:
        xl = np.linspace(min(x_vals) * 0.95, max(x_vals) * 1.05, 50)
        ax.plot(xl, fit["slope"] * xl + fit["intercept"], "k--", lw=0.8, alpha=0.4)
        r2_text = f"R² = {fit['r2']:.3f}"
        if note:
            r2_text += f"  ({note})"
        ax.text(0.05, 0.92, r2_text, transform=ax.transAxes, fontsize=8, va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8, ec="none"))
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.25)


# ═══════════════════════════════════════════════════════════════════════════
# 可视化
# ═══════════════════════════════════════════════════════════════════════════

def plot_dos_overview(results: dict, output: Path):
    """6 体系 TDOS 总览 (2×3 布局)."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 9))
    axes = axes.flatten()
    e_range = (-8, 4)

    for i, name in enumerate(SYSTEMS_ALL):
        ax = axes[i]
        dos = results[name].get("dos")
        if not dos:
            ax.text(0.5, 0.5, "No DOS", ha="center", va="center", transform=ax.transAxes)
            continue

        energies = dos.energies - dos.efermi
        mask = (energies >= e_range[0]) & (energies <= e_range[1])
        e_plot = energies[mask]
        dos_up = dos.get_densities(Spin.up)[mask]
        dos_down = dos.get_densities(Spin.down)[mask] if Spin.down in dos.densities else np.zeros_like(dos_up)

        color = SYSTEM_COLORS[name]
        ax.plot(e_plot, dos_up, color=color, lw=1.0, alpha=0.9)
        ax.plot(e_plot, -dos_down, color=color, lw=1.0, alpha=0.4)
        ax.fill_between(e_plot, 0, dos_up, color=color, alpha=0.2)
        ax.fill_between(e_plot, 0, -dos_down, color=color, alpha=0.1)
        ax.axvline(0, color="gray", ls="--", lw=0.6)
        _annotate_bandgap(ax, dos.get_gap())
        ax.set_xlim(e_range)
        ax.set_title(SYSTEM_LABELS[name], fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.25)

    # ── Companion data export ──
    plot_data = {"plot_type": "dos_overview", "parameters": {"e_range": list(e_range)}, "systems": []}
    for sys_name in SYSTEMS_ALL:
        dos = results[sys_name].get("dos")
        if dos:
            energies = dos.energies - dos.efermi
            mask = (energies >= e_range[0]) & (energies <= e_range[1])
            up = dos.get_densities(Spin.up)[mask]
            down = dos.get_densities(Spin.down)[mask] if Spin.down in dos.densities else np.zeros_like(up)
            plot_data["systems"].append({
                "name": sys_name, "label": SYSTEM_LABELS[sys_name],
                "energies": energies[mask], "dos_up": up, "dos_down": down,
                "efermi": dos.efermi, "gap": dos.get_gap(),
            })
    save_plot_data(plot_data, output)

    fig.suptitle("TDOS Overview — 6 Systems", fontsize=14, fontweight="bold")
    fig.text(0.5, 0.01, "E − E$_F$ (eV)", ha="center", fontsize=12)
    fig.text(0.01, 0.5, "DOS (states/eV)", va="center", rotation="vertical", fontsize=12)
    plt.tight_layout(rect=[0.02, 0.03, 1, 0.95])
    _save_figure(fig, output, "DOS overview")


def plot_dehydrogenation_pdos(results: dict, output: Path):
    """SH-3Co vs S-3Co PDOS 元素投影对比."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    e_range = (-8, 4)

    for idx, name in enumerate(["CoNiHOS-Co3", "CoNiHOS-Co3-noH"]):
        dos = results[name].get("dos")
        if not dos:
            continue
        energies = dos.energies - dos.efermi
        mask = (energies >= e_range[0]) & (energies <= e_range[1])
        e_plot = energies[mask]
        color = SYSTEM_COLORS[name]

        # Total DOS
        ax = axes[idx][0]
        dos_up = dos.get_densities(Spin.up)[mask]
        dos_down = dos.get_densities(Spin.down)[mask] if Spin.down in dos.densities else np.zeros_like(dos_up)
        ax.plot(e_plot, dos_up, color=color, lw=1.0)
        ax.plot(e_plot, -dos_down, color=color, lw=1.0, alpha=0.5)
        ax.fill_between(e_plot, 0, dos_up, color=color, alpha=0.15)
        ax.fill_between(e_plot, 0, -dos_down, color=color, alpha=0.08)
        ax.axvline(0, color="gray", ls="--", lw=0.6)
        _annotate_bandgap(ax, dos.get_gap())
        ax.set_xlim(e_range)
        ax.set_title(f"{SYSTEM_LABELS[name]} — TDOS", fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.25)

        # Site-projected: S
        ax = axes[idx][1]
        try:
            s_idx = next(i for i, s in enumerate(dos.structure) if s.species_string == "S")
            spd = dos.get_site_spd_dos(dos.structure[s_idx])
            if OrbitalType.p in spd:
                spd_p = spd[OrbitalType.p]
                for spin_label, spin, ls in [("up", Spin.up, "-"), ("down", Spin.down, "--")]:
                    if spin in spd_p.densities:
                        d = spd_p.densities[spin][mask]
                        ax.plot(e_plot, d if spin == Spin.up else -d,
                                color="#d62728" if spin == Spin.up else "#ff9896",
                                ls=ls, lw=1.2,
                                label=f"S-3p {spin_label}")
            ax.axvline(0, color="gray", ls="--", lw=0.6)
            ax.set_xlim(e_range)
            ax.set_title(f"{SYSTEM_LABELS[name]} — S-3p PDOS", fontsize=11, fontweight="bold")
            ax.legend(fontsize=8, loc="upper right")
            ax.grid(True, alpha=0.25)
        except StopIteration:
            ax.text(0.5, 0.5, "No S", ha="center", va="center", transform=ax.transAxes)

    # ── Companion data export ──
    plot_data = {"plot_type": "dehydrogenation_pdos", "parameters": {"e_range": list(e_range)}, "systems": []}
    for name in ["CoNiHOS-Co3", "CoNiHOS-Co3-noH"]:
        dos = results[name].get("dos")
        if dos:
            energies = dos.energies - dos.efermi
            mask = (energies >= e_range[0]) & (energies <= e_range[1])
            up = dos.get_densities(Spin.up)[mask]
            down = dos.get_densities(Spin.down)[mask] if Spin.down in dos.densities else np.zeros_like(up)
            s_p_up, s_p_down = [], []
            try:
                s_idx = next(i for i, s in enumerate(dos.structure) if s.species_string == "S")
                spd = dos.get_site_spd_dos(dos.structure[s_idx])
                if OrbitalType.p in spd:
                    spd_p = spd[OrbitalType.p]
                    s_p_up = spd_p.densities[Spin.up][mask] if Spin.up in spd_p.densities else []
                    s_p_down = spd_p.densities[Spin.down][mask] if Spin.down in spd_p.densities else []
            except Exception:
                pass
            plot_data["systems"].append({
                "name": name, "label": SYSTEM_LABELS[name],
                "energies": energies[mask], "tdos_up": up, "tdos_down": down,
                "s_pdos_up": s_p_up, "s_pdos_down": s_p_down,
                "efermi": dos.efermi, "gap": dos.get_gap(),
            })
    save_plot_data(plot_data, output)

    fig.suptitle("Dehydrogenation Effect: SH-3Co vs S-3Co", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    _save_figure(fig, output, "Dehydrogenation PDOS")


def plot_s_local_magmom(results: dict, output: Path):
    """S 周围 3 金属磁矩条形图 (6 体系)."""
    doped = [n for n in SYSTEMS_ALL if n != "CoNiHO"]  # All doped, including noH
    n_sys = len(doped)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Bader charge transfer
    ax = axes[0]
    x = np.arange(n_sys)
    width = 0.25
    for j, pos in enumerate([0, 1, 2]):  # 3 nearest metals
        vals = []
        for name in doped:
            bl = results[name].get("bader_local", {})
            nn = bl.get("nn_metal_bader", [])
            vals.append(nn[j]["charge_transfer"] if j < len(nn) and nn[j].get("charge_transfer") is not None else np.nan)
        bars = ax.bar(x + j * width, vals, width, label=f"M{j+1}", alpha=0.85)
    ax.set_xticks(x + width)
    ax.set_xticklabels([SYSTEM_LABELS[n] for n in doped], rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Charge Transfer (e)", fontsize=12)
    ax.set_title("S-Neighbor Metals: Bader Charge Transfer", fontsize=12, fontweight="bold")
    ax.axhline(0, color="gray", ls="-", lw=0.8)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25, axis="y")

    # Magnetic moments
    ax = axes[1]
    for j, pos in enumerate([0, 1, 2]):
        vals = []
        for name in doped:
            bl = results[name].get("bader_local", {})
            nn = bl.get("nn_metal_bader", [])
            vals.append(nn[j]["magmom"] if j < len(nn) and nn[j].get("magmom") is not None else np.nan)
        ax.bar(x + j * width, vals, width, label=f"M{j+1}", alpha=0.85)
    ax.set_xticks(x + width)
    ax.set_xticklabels([SYSTEM_LABELS[n] for n in doped], rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Magnetic Moment (μB)", fontsize=12)
    ax.set_title("S-Neighbor Metals: Local Magnetic Moment", fontsize=12, fontweight="bold")
    ax.axhline(0, color="gray", ls="-", lw=0.8)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25, axis="y")

    # ── Companion data export ──
    plot_data = {"plot_type": "s_local_charge_mag", "parameters": {}, "systems": []}
    for name in doped:
        bl = results[name].get("bader_local", {})
        nn = bl.get("nn_metal_bader", [])
        cts = [nn[j]["charge_transfer"] if j < len(nn) and nn[j].get("charge_transfer") is not None else None for j in range(3)]
        mags = [nn[j]["magmom"] if j < len(nn) and nn[j].get("magmom") is not None else None for j in range(3)]
        plot_data["systems"].append({
            "name": name, "label": SYSTEM_LABELS[name],
            "nn_charge_transfers": cts, "nn_magmoms": mags,
            "S_charge_transfer": bl.get("S_charge_transfer"),
            "S_magmom": bl.get("S_magmom"),
        })
    save_plot_data(plot_data, output)

    fig.suptitle("S-Local Environment: Charge & Magnetism", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    _save_figure(fig, output, "S-local magmom")


def plot_band_overview(results: dict, output: Path):
    """6 体系能带总览 (2×3)."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    for i, name in enumerate(SYSTEMS_ALL):
        ax = axes[i]
        band = results[name].get("band")
        if band and results[name].get("is_line_mode"):
            plotter = BSPlotter(band)
            pd_data = plotter.bs_plot_data(zero_to_efermi=True)
            all_e = []
            for bidx, dists in enumerate(pd_data["distances"]):
                for spin_key in ("1", "-1"):
                    if spin_key in pd_data["energy"] and bidx < len(pd_data["energy"][spin_key]):
                        for ib in range(len(pd_data["energy"][spin_key][bidx])):
                            e_vals = pd_data["energy"][spin_key][bidx][ib]
                            color = "#2171b5" if spin_key == "1" else "#cb181d"
                            alpha = 0.85 if spin_key == "1" else 0.5
                            ax.plot(dists, e_vals, color=color, lw=0.6, alpha=alpha)
                            all_e.extend(e_vals)
            # Bandgap annotation
            vbm = results[name].get("vbm")
            cbm = results[name].get("cbm")
            if vbm is not None and cbm is not None:
                band_gap = cbm - vbm
                _annotate_bandgap(ax, band_gap, orientation="h")
            ax.axhline(0, color="gray", ls="--", lw=0.6)
            if all_e:
                m = (max(all_e) - min(all_e)) * 0.1
                ax.set_ylim(min(all_e) - m, max(all_e) + m)
            ticks = pd_data["ticks"]
            for dist, label in zip(ticks["distance"], ticks["label"]):
                ax.axvline(dist, color="k", ls="-", lw=0.4, alpha=0.25)
            ax.set_xticks(ticks["distance"])
            ax.set_xticklabels([l.replace("$\\Gamma$", "Γ").replace("$", "") for l in ticks["label"]], fontsize=7)
            ax.set_xlim(ticks["distance"][0], ticks["distance"][-1])
        else:
            ax.text(0.5, 0.5, "No band data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(SYSTEM_LABELS[name], fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.25)

    # ── Companion data export ──
    plot_data = {"plot_type": "band_overview", "parameters": {}, "systems": []}
    for name in SYSTEMS_ALL:
        band = results[name].get("band")
        if band and results[name].get("is_line_mode"):
            plotter = BSPlotter(band)
            pd_data = plotter.bs_plot_data(zero_to_efermi=True)
            bands_up, bands_down = [], []
            for bidx, dists in enumerate(pd_data["distances"]):
                for spin_key, spin_data in [("1", bands_up), ("-1", bands_down)]:
                    if spin_key in pd_data["energy"] and bidx < len(pd_data["energy"][spin_key]):
                        for ib in range(len(pd_data["energy"][spin_key][bidx])):
                            spin_data.append({"band_idx": ib, "distances": dists,
                                              "energies": pd_data["energy"][spin_key][bidx][ib]})
            vbm = results[name].get("vbm")
            cbm = results[name].get("cbm")
            plot_data["systems"].append({
                "name": name, "label": SYSTEM_LABELS[name],
                "efermi": band.efermi,
                "vbm": vbm, "cbm": cbm,
                "k_labels": pd_data["ticks"]["label"],
                "k_distances": pd_data["ticks"]["distance"],
                "bands_up": bands_up, "bands_down": bands_down,
            })
    save_plot_data(plot_data, output)

    fig.suptitle("Band Structure Overview — 6 Systems", fontsize=14, fontweight="bold")
    fig.text(0.5, 0.005, "k-path", ha="center", fontsize=11)
    fig.text(0.005, 0.5, "E − E$_F$ (eV)", va="center", rotation="vertical", fontsize=11)
    plt.tight_layout(rect=[0.02, 0.02, 1, 0.96])
    _save_figure(fig, output, "Band overview")


# ═══════════════════════════════════════════════════════════════════════════
# F ─ 脱氢态磁矩异常诊断 (CoNiHOS-Co3-noH: 45→36 μB)
# ═══════════════════════════════════════════════════════════════════════════

def diagnose_magnetic_anomaly(data_dir: Path, ref_dir: Path):
    """诊断 CoNiHOS-Co3-noH 脱氢态磁矩从 45 μB 骤降至 36 μB 的原因.

    检查项:
      1. INCAR NELECT / MAGMOM 设置 (via pymatgen Incar)
      2. 理论电子数 vs Outcar.nelect
      3. OUTCAR 磁矩收敛轨迹 (手动解析 — pymatgen Outcar 不暴露电子步历史)
    """
    import gzip  # 仅用于 OUTCAR 电子步磁矩收敛解析 (无 pymatgen API)

    print("\n" + "=" * 80)
    print("F. MAGNETIC CONVERGENCE CHECK: SH-3Co vs S-3Co-noH (NELECT=422 fix)")
    print("=" * 80)

    noh_dir = data_dir / "CoNiHOS-Co3-noH" / "3-static"
    sh_dir = ref_dir / "CoNiHOS-Co3" / "3-static"

    # ── 1. INCAR 检查 (pymatgen Incar.from_file) ──
    print("\n── 1. INCAR NELECT / MAGMOM 参数检查 ──")
    for label, d in [("SH-3Co", sh_dir), ("S-3Co-noH", noh_dir)]:
        incar_path = d / "INCAR.gz"
        if incar_path.exists():
            incar = Incar.from_file(str(incar_path))
            nelect_val = incar.get("NELECT", None)
            magmom_val = incar.get("MAGMOM", None)
            if nelect_val is not None:
                print(f"  {label} INCAR NELECT = {nelect_val}")
            else:
                print(f"  {label} INCAR: NELECT 未显式设置 (使用默认值 = sum of ZVAL)")
            if magmom_val is not None:
                print(f"  {label} INCAR MAGMOM = {magmom_val}")
            else:
                print(f"  {label} INCAR: MAGMOM 未设置")

    # ── 2. 理论电子数验证 (Structure.from_file + Outcar.nelect) ──
    print("\n── 2. 理论电子数验证 ──")

    for label, d in [("SH-3Co", sh_dir), ("S-3Co-noH", noh_dir)]:
        poscar = d / "CONTCAR.gz"
        outcar_path = d / "OUTCAR.gz"
        if poscar.exists():
            struct = Structure.from_file(poscar)
            n_atoms = len(struct)
            zval_sum = sum(ZVAL.get(s.species_string, 0) for s in struct)
            counts = {}
            for s in struct:
                counts[s.species_string] = counts.get(s.species_string, 0) + 1
            # 交叉验证: Outcar.nelect vs Σ ZVAL
            actual_nelect = None
            if outcar_path.exists():
                try:
                    actual_nelect = Outcar(outcar_path).nelect
                except Exception:
                    pass
            nelect_str = f"NELECT(实际) = {actual_nelect:.1f}" if actual_nelect is not None else "N/A"
            print(f"  {label}: {n_atoms} atoms, 组成: {counts}")
            print(f"         Σ ZVAL(理论) = {zval_sum:.0f}, {nelect_str}")

    # ── 3. OUTCAR 磁矩收敛过程 ──
    print("\n── 3. CoNiHOS-Co3-noH 磁矩收敛轨迹 (OUTCAR) ──")
    outcar_path = noh_dir / "OUTCAR.gz"
    if outcar_path.exists():
        with gzip.open(outcar_path, "rt") as f:
            lines = f.readlines()

        mag_hist = []
        for l in lines:
            if "number of electron" in l:
                parts = l.split()
                try:
                    mag_idx = parts.index("magnetization")
                    mag = float(parts[mag_idx + 1])
                    nelec = float(parts[parts.index("electron") + 1])
                    mag_hist.append((nelec, mag))
                except (ValueError, IndexError):
                    pass

        if mag_hist:
            print(f"  总步数: {len(mag_hist)}")
            print(f"  初始磁矩: {mag_hist[0][1]:.2f} μB  (NELECT = {mag_hist[0][0]:.1f})")
            print(f"  最终磁矩: {mag_hist[-1][1]:.2f} μB  (NELECT = {mag_hist[-1][0]:.1f})")

            # 收敛过程概览 (选取代表性步骤)
            n = len(mag_hist)
            indices = [0]
            for i in range(1, n):
                if abs(mag_hist[i][1] - mag_hist[i - 1][1]) > 5:
                    indices.append(i)
            if n - 1 not in indices:
                indices.append(n - 1)
            for frac in [0.25, 0.5, 0.75]:
                idx = int(n * frac)
                if idx not in indices:
                    indices.append(idx)
            indices = sorted(set(indices))[:12]

            print(f"\n  磁矩演化 (selective steps):")
            print(f"    {'Step':>6s}  {'NELECT':>10s}  {'Mag (μB)':>10s}  {'ΔMag':>10s}")
            prev_mag = mag_hist[0][1]
            for idx in indices:
                ne, mag = mag_hist[idx]
                delta = mag - prev_mag
                marker = " ← 突变" if abs(delta) > 10 else ""
                print(f"    {idx:>6d}  {ne:10.1f}  {mag:10.4f}  {delta:+10.4f}{marker}")
                prev_mag = mag

            # 最后10步
            print(f"\n  最后10步磁矩收敛细节:")
            print(f"    {'Step':>6s}  {'NELECT':>10s}  {'Mag (μB)':>12s}")
            for i in range(max(0, n - 10), n):
                ne, mag = mag_hist[i]
                print(f"    [{i:>4d}]  {ne:10.1f}  {mag:12.6f}")



# ═══════════════════════════════════════════════════════════════════════════
# G ─ 缺陷形成能的化学势修正
# ═══════════════════════════════════════════════════════════════════════════

# Reference chemical potentials from molecular DFT (user-provided)
MU_H = -3.375       # ½ H₂
MU_O = -4.925       # ½ O₂
MU_S = -2.9375      # ⅛ S₈


def analyze_formation_energy(results: dict):
    """计算 S 掺杂的绝对缺陷形成能 (含化学势修正).

    反应方程式:
      SH⁻ 掺杂: OH⁻ → SH⁻,  ΔE_f = E(doped) - E(pristine) - μ_S + μ_O
      S²⁻ 掺杂: OH⁻ + H⁺ → S²⁻ + □,  ΔE_f = E(noH) + μ_H - E(pristine) - μ_S + μ_O
    """
    print("\n" + "=" * 80)
    print("G. DEFECT FORMATION ENERGY WITH CHEMICAL POTENTIAL CORRECTION")
    print("=" * 80)
    print(f"  Reference chemical potentials:")
    print(f"    μ_H  = {MU_H:+.4f} eV  (½ H₂ = -6.75 eV)")
    print(f"    μ_O  = {MU_O:+.4f} eV  (½ O₂ = -9.85 eV)")
    print(f"    μ_S  = {MU_S:+.4f} eV  (⅛ S₈ = -23.50 eV)")

    e_pristine = results["CoNiHO"].get("energy")
    if e_pristine is None:
        print("  ERROR: Pristine energy not available")
        return

    print(f"\n  Formation energies (relative to pristine CoNiHO + reservoirs):")
    print(f"  {'System':<22s}  {'Raw ΔE/eV':>12s}  {'Correction':>12s}  {'ΔE_f/eV':>12s}  {'ΔE_f/atom':>14s}")
    print(f"  {'─'*22}  {'─'*12}  {'─'*12}  {'─'*12}  {'─'*14}")

    rows = []
    for name in SYSTEMS_DOPED_ALL:
        e_d = results[name].get("energy")
        if e_d is None:
            continue
        raw_de = e_d - e_pristine

        n_atoms = results[name].get("nsites", "?")
        if name == "CoNiHOS-Co3-noH":
            # S²⁻ doping: remove OH, add S, lose H to reservoir
            correction = -MU_S + MU_O + MU_H  # +μ_H because H goes back to reservoir
        else:
            # SH⁻ doping: replace OH with SH
            correction = -MU_S + MU_O

        de_f = raw_de + correction
        de_f_per_atom = de_f / n_atoms * 1000 if isinstance(n_atoms, (int, float)) else None
        row = {
            "System": SYSTEM_LABELS.get(name, name),
            "Raw ΔE (eV)": round(raw_de, 4),
            "Correction (eV)": round(correction, 4),
            "ΔE_f (eV)": round(de_f, 4),
            "ΔE_f (meV/atom)": round(de_f_per_atom, 1) if de_f_per_atom else None,
        }
        rows.append(row)
        corr_s = f"{correction:+12.4f}"
        de_s = f"{de_f:+12.4f}"
        pa_s = f"{de_f_per_atom:+14.1f}" if de_f_per_atom else "N/A"
        print(f"  {name:<22s}  {raw_de:+12.4f}  {corr_s}  {de_s}  {pa_s}")

    if rows:
        print(f"\n  Correction formula:")
        print(f"    SH⁻ doping: correction = −μ_S + μ_O = {(-MU_S + MU_O):+.4f} eV")
        print(f"    S²⁻ doping: correction = −μ_S + μ_O + μ_H = {(-MU_S + MU_O + MU_H):+.4f} eV")
        print(f"    ΔE_f = Raw ΔE + Correction")

    return rows


# ═══════════════════════════════════════════════════════════════════════════
# H ─ 元素分辨 PDOS 总览图 (2×3)
# ═══════════════════════════════════════════════════════════════════════════

def plot_pdos_overview(results: dict, output: Path):
    """6 体系元素分辨 PDOS: Co-3d, Ni-3d, O-2p, S-3p (自旋分辨)."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    e_range = (-8, 4)

    element_config = {
        "Co": ("d", "#2171b5", "#6baed6"),     # 3d: blue
        "Ni": ("d", "#238b45", "#74c476"),     # 3d: green
        "O":  ("p", "#cb181d", "#fb6a4a"),     # 2p: red
        "S":  ("p", "#d94801", "#fd8d3c"),     # 3p: orange (×5 scaled)
    }

    for i, name in enumerate(SYSTEMS_ALL):
        ax = axes[i]
        dos = results[name].get("dos")
        if not dos:
            ax.text(0.5, 0.5, "No DOS", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(SYSTEM_LABELS[name], fontsize=11, fontweight="bold")
            continue

        energies = dos.energies - dos.efermi
        mask = (energies >= e_range[0]) & (energies <= e_range[1])
        e_plot = energies[mask]

        struct = dos.structure
        has_S = any(s.species_string == "S" for s in struct)

        for elem, (orb_str, color_up, color_down) in element_config.items():
            if elem == "S" and not has_S:
                continue
            try:
                spd = dos.get_element_spd_dos(elem)
                orb = OrbitalType.d if orb_str == "d" else OrbitalType.p
                if orb not in spd:
                    continue
                elem_dos = spd[orb]
                scale = 5.0 if elem == "S" else 1.0
                for spin, spin_label, color, ls in [
                    (Spin.up, "up", color_up, "-"),
                    (Spin.down, "down", color_down, "--"),
                ]:
                    if spin in elem_dos.densities:
                        d = elem_dos.densities[spin][mask] * scale
                        ax.plot(e_plot, d if spin == Spin.up else -d,
                                color=color, ls=ls, lw=1.0, alpha=0.85)
            except Exception:
                pass

        _annotate_bandgap(ax, dos.get_gap())
        ax.axvline(0, color="gray", ls="--", lw=0.6)
        ax.set_xlim(e_range)
        ylim = ax.get_ylim()
        ax.set_ylim(ylim[0], ylim[1] * 1.15)
        ax.set_title(SYSTEM_LABELS[name], fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.25)

    # Legend on last panel
    legend_elements = [
        Line2D([0], [0], color="#2171b5", lw=1.5, label="Co-3d (up)"),
        Line2D([0], [0], color="#2171b5", lw=1.5, ls="--", alpha=0.5, label="Co-3d (down)"),
        Line2D([0], [0], color="#238b45", lw=1.5, label="Ni-3d (up)"),
        Line2D([0], [0], color="#238b45", lw=1.5, ls="--", alpha=0.5, label="Ni-3d (down)"),
        Line2D([0], [0], color="#cb181d", lw=1.5, label="O-2p (up)"),
        Line2D([0], [0], color="#cb181d", lw=1.5, ls="--", alpha=0.5, label="O-2p (down)"),
        Line2D([0], [0], color="#d94801", lw=1.5, label="S-3p (×5, up)"),
        Line2D([0], [0], color="#d94801", lw=1.5, ls="--", alpha=0.5, label="S-3p (×5, down)"),
    ]
    axes[-1].legend(handles=legend_elements, fontsize=7, loc="upper right", ncol=2,
                    framealpha=0.9, columnspacing=0.5, handlelength=1.2)

    # ── Companion data export ──
    plot_data = {"plot_type": "pdos_overview", "parameters": {"e_range": list(e_range), "S_scale": 5.0}, "systems": []}
    for name in SYSTEMS_ALL:
        dos = results[name].get("dos")
        if dos:
            energies = dos.energies - dos.efermi
            mask = (energies >= e_range[0]) & (energies <= e_range[1])
            sys_data = {"name": name, "label": SYSTEM_LABELS[name],
                       "energies": energies[mask], "efermi": dos.efermi, "gap": dos.get_gap()}
            for elem in ["Co", "Ni", "O", "S"]:
                try:
                    spd = dos.get_element_spd_dos(elem)
                    orb = OrbitalType.d if elem in ("Co", "Ni") else OrbitalType.p
                    if orb in spd:
                        for spin_label, spin in [("up", Spin.up), ("down", Spin.down)]:
                            if spin in spd[orb].densities:
                                key = f"{elem}_{'d' if elem in ('Co','Ni') else 'p'}_{spin_label}"
                                sys_data[key] = spd[orb].densities[spin][mask]
                except Exception:
                    pass
            plot_data["systems"].append(sys_data)
    save_plot_data(plot_data, output)

    fig.suptitle("Element-Resolved PDOS — 6 Systems", fontsize=14, fontweight="bold")
    fig.text(0.5, 0.01, "E − E$_F$ (eV)", ha="center", fontsize=12)
    fig.text(0.01, 0.5, "PDOS (states/eV)", va="center", rotation="vertical", fontsize=12)
    plt.tight_layout(rect=[0.02, 0.03, 1, 0.95])
    _save_figure(fig, output, "PDOS overview")


# ═══════════════════════════════════════════════════════════════════════════
# I ─ 构效关系散点图 (3-panel)
# ═══════════════════════════════════════════════════════════════════════════

def _linear_fit(x: np.ndarray, y: np.ndarray):
    """Return (slope, intercept, r_squared) or (None, None, None) if invalid."""
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    if len(x) < 3:
        return None, None, None
    slope, intercept, rval, _pval, _err = stats.linregress(x, y)
    return slope, intercept, rval ** 2


def _collect_scatter_data(results, names, x_key, y_key, x_sub=None, y_sub=None):
    """Collect (x, y, label, color) for scatter plots from results dict."""
    xv, yv, lb, cl = [], [], [], []
    for name in names:
        r = results[name]
        x = r.get(x_key) if x_sub is None else r.get(x_sub, {}).get(x_key)
        y = r.get(y_key) if y_sub is None else r.get(y_sub, {}).get(y_key)
        if x is not None and y is not None:
            xv.append(x); yv.append(y)
            lb.append(SYSTEM_LABELS[name]); cl.append(SYSTEM_COLORS[name])
    return xv, yv, lb, cl


def plot_structure_activity(results: dict, output: Path):
    """3-panel structure-activity scatter plots."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    doped = SYSTEMS_DOPED_ALL

    # ── Panel 1: Co d-band center vs Bandgap ──
    xv, yv, lb, cl = _collect_scatter_data(results, doped, "d_center_Co", "bandgap_dos")
    rp = results["CoNiHO"]
    if rp.get("d_center_Co") is not None and rp.get("bandgap_dos") is not None:
        xv.append(rp["d_center_Co"]); yv.append(rp["bandgap_dos"])
        lb.append("Pristine"); cl.append(SYSTEM_COLORS["CoNiHO"])
    fit1 = dict(zip(("slope", "intercept", "r2"), _linear_fit(np.array(xv), np.array(yv))))
    _plot_scatter_panel(axes[0], xv, yv, lb, cl, "Co d-band center (eV)", "Band Gap (eV)",
                         "d-band Center vs Bandgap", fit1, "weak expected")
    axes[0].axhline(0, color="gray", ls=":", lw=0.6)
    p1_points = [{"name": n, "x": float(x), "y": float(y)} for n, x, y in zip(lb, xv, yv)]

    # ── Panel 2: S-M avg bond length vs S charge transfer ──
    xv2, yv2, lb2, cl2 = _collect_scatter_data(results, doped, "avg_bond_length", "S_charge_transfer",
                                                  x_sub="S_bonds", y_sub="bader_local")
    fit2 = dict(zip(("slope", "intercept", "r2"), _linear_fit(np.array(xv2), np.array(yv2))))
    _plot_scatter_panel(axes[1], xv2, yv2, lb2, cl2, "S–M Avg Bond Length (Å)", "S Charge Transfer (e)",
                         "Bond Length vs Charge Transfer", fit2)
    p2_points = [{"name": n, "x": float(x), "y": float(y)} for n, x, y in zip(lb2, xv2, yv2)]

    # ── Panel 3: NN metal avg magmom vs NN metal avg charge transfer ──
    xv3, yv3, lb3, cl3 = _collect_scatter_data(results, doped, "nn_avg_transfer", "nn_avg_magmom",
                                                  y_sub="bader_local", x_sub="bader_local")
    fit3 = dict(zip(("slope", "intercept", "r2"), _linear_fit(np.array(xv3), np.array(yv3))))
    _plot_scatter_panel(axes[2], xv3, yv3, lb3, cl3, "NN Metal Avg Charge Transfer (e)",
                         "NN Metal Avg Magmom (μB)", "Charge Transfer vs Local Magnetism", fit3)
    p3_points = [{"name": n, "x": float(x), "y": float(y)} for n, x, y in zip(lb3, xv3, yv3)]

    # ── Companion data export ──
    panels_data = [
        {"title": "d-band Center vs Bandgap", "x_label": "Co d-band center (eV)",
         "y_label": "Band Gap (eV)", "points": p1_points, "fit": fit1 if fit1["slope"] else None},
        {"title": "Bond Length vs Charge Transfer", "x_label": "S-M Avg Bond Length (Å)",
         "y_label": "S Charge Transfer (e)", "points": p2_points, "fit": fit2 if fit2["slope"] else None},
        {"title": "Charge Transfer vs Local Magnetism", "x_label": "NN Metal Avg Charge Transfer (e)",
         "y_label": "NN Metal Avg Magmom (μB)", "points": p3_points, "fit": fit3 if fit3["slope"] else None},
    ]
    save_plot_data({"plot_type": "structure_activity", "parameters": {}, "panels": panels_data}, output)

    fig.suptitle("Structure–Activity Relationships: S-Doped CoNi Hydroxide", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    _save_figure(fig, output, "Structure-activity")


# ═══════════════════════════════════════════════════════════════════════════
# J ─ 主数据导出 (AI-readable structured JSON)
# ═══════════════════════════════════════════════════════════════════════════

def export_all_data(results: dict, output_dir: Path):
    """Export all extracted data as a single structured JSON file."""
    export = {
        "meta": {
            "description": "S-doped CoNi hydroxide 6-system analysis — pure data export",
            "n_systems": len(SYSTEMS_ALL),
            "zval_reference": ZVAL,
            "chemical_potentials": {"mu_H": MU_H, "mu_O": MU_O, "mu_S": MU_S},
        },
        "overview": {},
        "doping_effect": {},
        "dehydrogenation": {},
        "coordination": {},
        "bader_comprehensive": {},
        "magnetic_check": {},
        "formation_energy": {},
    }

    # ── Overview ──
    df = overview_table(results)
    export["overview"] = {
        "columns": list(df.columns),
        "rows": df.to_dict(orient="records"),
    }

    # ── Doping effect (Pristine vs SH-3Co) ──
    p = results.get("CoNiHO", {})
    d = results.get("CoNiHOS-Co3", {})
    export["doping_effect"] = {
        "pristine_energy": p.get("energy"),
        "sh3co_energy": d.get("energy"),
        "delta_e": (d.get("energy", 0) - p.get("energy", 0)) if p.get("energy") and d.get("energy") else None,
        "pristine_lattice": p.get("lattice", {}),
        "sh3co_lattice": d.get("lattice", {}),
        "pristine_gap": p.get("bandgap_dos"),
        "sh3co_gap": d.get("bandgap_dos"),
        "pristine_total_mag": p.get("total_mag"),
        "sh3co_total_mag": d.get("total_mag"),
        "pristine_co_avg_mag": p.get("Co_avg_mag"),
        "pristine_ni_avg_mag": p.get("Ni_avg_mag"),
        "sh3co_co_avg_mag": d.get("Co_avg_mag"),
        "sh3co_ni_avg_mag": d.get("Ni_avg_mag"),
        "s_bond_lengths": d.get("S_bonds", {}).get("bond_lengths", []),
        "s_metal_types": d.get("S_bonds", {}).get("metal_types", []),
    }

    # ── Dehydrogenation (SH-3Co vs S-3Co-noH) ──
    sh = results.get("CoNiHOS-Co3", {})
    s0 = results.get("CoNiHOS-Co3-noH", {})
    bl_sh = sh.get("bader_local", {})
    bl_s = s0.get("bader_local", {})
    export["dehydrogenation"] = {
        "sh3co_energy": sh.get("energy"),
        "s3co_noh_energy": s0.get("energy"),
        "raw_delta_e": (s0.get("energy", 0) - sh.get("energy", 0)) if sh.get("energy") and s0.get("energy") else None,
        "sh3co_gap": sh.get("bandgap_dos"),
        "s3co_noh_gap": s0.get("bandgap_dos"),
        "sh3co_s_bonds": sh.get("S_bonds", {}).get("bond_lengths", []),
        "s3co_noh_s_bonds": s0.get("S_bonds", {}).get("bond_lengths", []),
        "sh3co_S_charge_transfer": bl_sh.get("S_charge_transfer"),
        "sh3co_S_magmom": bl_sh.get("S_magmom"),
        "s3co_noh_S_charge_transfer": bl_s.get("S_charge_transfer"),
        "s3co_noh_S_magmom": bl_s.get("S_magmom"),
        "sh3co_nn_metals": bl_sh.get("nn_metal_bader", []),
        "s3co_noh_nn_metals": bl_s.get("nn_metal_bader", []),
        "sh3co_co_dcenter": sh.get("d_center_Co"),
        "sh3co_ni_dcenter": sh.get("d_center_Ni"),
        "s3co_noh_co_dcenter": s0.get("d_center_Co"),
        "s3co_noh_ni_dcenter": s0.get("d_center_Ni"),
        "sh3co_pdos_ef": sh.get("pdos_ef", {}),
        "s3co_noh_pdos_ef": s0.get("pdos_ef", {}),
    }

    # ── Coordination (4 SH-doped) ──
    coord_data = {}
    for name in SYSTEMS_DOPED_SH:
        r = results.get(name, {})
        sb = r.get("S_bonds", {})
        bl = r.get("bader_local", {})
        coord_data[name] = {
            "label": SYSTEM_LABELS[name],
            "energy": r.get("energy"),
            "gap": r.get("bandgap_dos"),
            "s_bond_lengths": sb.get("bond_lengths", []),
            "s_metal_types": sb.get("metal_types", []),
            "s_avg_bond": sb.get("avg_bond_length"),
            "nn_magmoms": [x.get("magmom") for x in bl.get("nn_metal_bader", [])],
            "nn_charge_transfers": [x.get("charge_transfer") for x in bl.get("nn_metal_bader", [])],
            "S_charge_transfer": bl.get("S_charge_transfer"),
            "S_magmom": bl.get("S_magmom"),
            "co_dcenter": r.get("d_center_Co"),
            "ni_dcenter": r.get("d_center_Ni"),
            "pdos_ef": r.get("pdos_ef", {}),
        }
    export["coordination"] = coord_data

    # ── Bader comprehensive ──
    bader_data = {}
    for name in SYSTEMS_DOPED_ALL:
        bl = results.get(name, {}).get("bader_local", {})
        bader_data[name] = {
            "label": SYSTEM_LABELS[name],
            "S_bader_charge": bl.get("S_bader_charge"),
            "S_charge_transfer": bl.get("S_charge_transfer"),
            "S_magmom": bl.get("S_magmom"),
            "S_atomic_vol": bl.get("S_atomic_vol"),
            "nn_metal_bader": bl.get("nn_metal_bader", []),
            "Co_avg_bader_charge": bl.get("Co_avg_bader_charge"),
            "Ni_avg_bader_charge": bl.get("Ni_avg_bader_charge"),
            "Co_avg_magmom": bl.get("Co_avg_magmom"),
            "Ni_avg_magmom": bl.get("Ni_avg_magmom"),
        }
    export["bader_comprehensive"] = bader_data

    # ── Magnetic check (raw INCAR/OUTCAR values) ──
    mag_check = {}
    for name in ["CoNiHOS-Co3", "CoNiHOS-Co3-noH"]:
        d = DATA_DIR / name / "3-static"
        entry = {"name": name, "label": SYSTEM_LABELS[name]}
        incar_path = d / "INCAR.gz"
        if incar_path.exists():
            incar = Incar.from_file(str(incar_path))
            entry["INCAR_NELECT"] = incar.get("NELECT", None)
            entry["INCAR_MAGMOM"] = incar.get("MAGMOM", None)
        outcar_path = d / "OUTCAR.gz"
        if outcar_path.exists():
            outcar = Outcar(outcar_path)
            entry["OUTCAR_nelect"] = outcar.nelect
            entry["OUTCAR_total_mag"] = outcar.total_mag
        poscar = d / "CONTCAR.gz"
        if poscar.exists():
            struct = Structure.from_file(poscar)
            entry["n_atoms"] = len(struct)
            entry["composition"] = dict(Counter(s.species_string for s in struct))
            entry["zval_sum"] = sum(ZVAL.get(s.species_string, 0) for s in struct)
        mag_check[name] = entry
    export["magnetic_check"] = mag_check

    # ── Formation energy ──
    e_pristine = results.get("CoNiHO", {}).get("energy")
    fe_data = []
    if e_pristine:
        for name in SYSTEMS_DOPED_ALL:
            e_d = results.get(name, {}).get("energy")
            n_atoms = results.get(name, {}).get("nsites", "?")
            if e_d:
                raw_de = e_d - e_pristine
                if name == "CoNiHOS-Co3-noH":
                    correction = -MU_S + MU_O + MU_H
                else:
                    correction = -MU_S + MU_O
                de_f = raw_de + correction
                fe_data.append({
                    "system": name, "label": SYSTEM_LABELS[name],
                    "raw_delta_e": raw_de, "correction": correction, "delta_e_f": de_f,
                    "n_atoms": n_atoms,
                    "delta_e_f_mev_per_atom": de_f / n_atoms * 1000 if isinstance(n_atoms, (int, float)) else None,
                })
    export["formation_energy"] = fe_data

    json_path = output_dir / "all_analysis_data.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(export, f, indent=2, default=_json_default)
    log.info("Master data export → %s", json_path)
    return export


# ═══════════════════════════════════════════════════════════════════════════
# 主函数
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="S-doped CoNi hydroxide 6-system analysis")
    parser.add_argument("--parallel", action="store_true",
                        help="Use parallel execution (default: serial)")
    parser.add_argument("--skip-plots", action="store_true",
                        help="Skip all figure generation")
    args = parser.parse_args()

    base_dir = DATA_DIR
    dirs = {name: base_dir / name for name in SYSTEMS_ALL}

    print("=" * 80)
    print("COMPREHENSIVE ANALYSIS: S-DOPED CoNi HYDROXIDE (6 Systems)")
    print("=" * 80)
    print()
    for name in SYSTEMS_ALL:
        print(f"  {name:<25s} → {SYSTEM_LABELS[name]}")
    print()

    # ── Analysis (serial by default) ──
    use_parallel = args.parallel
    mode_str = "parallel" if use_parallel else "serial"
    log.info("Analyzing 6 systems (%s)...", mode_str)
    results = {}

    if use_parallel:
        with ProcessPoolExecutor(max_workers=len(SYSTEMS_ALL)) as executor:
            futures = {executor.submit(analyze_system, dir_path, name): name
                       for name, dir_path in dirs.items()}
            for future in as_completed(futures):
                name = futures[future]
                try:
                    results[name] = future.result()
                    log.info("  ✓ %s", name)
                except Exception as e:
                    log.error("  ✗ %s: %s", name, e)
                    results[name] = {"name": name}
    else:
        for name, dir_path in dirs.items():
            try:
                results[name] = analyze_system(dir_path, name)
                log.info("  ✓ %s", name)
            except Exception as e:
                log.error("  ✗ %s: %s", name, e)
                results[name] = {"name": name}

    # ── Master data export ──
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    export_all_data(results, OUTPUT_DIR)

    # ── A: Overview table ──
    print("\n" + "=" * 80)
    print("A. 6-SYSTEM OVERVIEW TABLE")
    print("=" * 80)
    print("  注: Bader 电荷转移 = ZVAL - bader_charge, 正值 = 失去电子, 负值 = 获得电子")
    df = overview_table(results)
    print(df.to_string(index=False))
    csv_path = OUTPUT_DIR / "overview_6_systems.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # ── B: Doping effect (Pristine → SH-3Co) ──
    analyze_doping_effect(results)

    # ── C: Dehydrogenation (SH-3Co → S-3Co-noH) ──
    analyze_dehydrogenation(results)

    # ── D: Coordination effect (4 SH-doped) ──
    analyze_coordination(results)

    # ── E: Bader comprehensive ──
    analyze_bader_comprehensive(results)

    # ── F: Magnetic anomaly diagnosis ──
    diagnose_magnetic_anomaly(base_dir, base_dir)

    # ── G: Formation energy with chemical potentials ──
    analyze_formation_energy(results)

    # ── Plots ──
    if not args.skip_plots:
        print("\n" + "=" * 80)
        print("GENERATING FIGURES")
        print("=" * 80)
        plot_dos_overview(results, OUTPUT_DIR / "dos_overview_6")
        plot_dehydrogenation_pdos(results, OUTPUT_DIR / "dehydrogenation_pdos")
        plot_s_local_magmom(results, OUTPUT_DIR / "s_local_charge_mag")
        plot_band_overview(results, OUTPUT_DIR / "band_overview_6")
        plot_pdos_overview(results, OUTPUT_DIR / "pdos_overview_6")
        plot_structure_activity(results, OUTPUT_DIR / "structure_activity")
    else:
        log.info("Skipping plots (--skip-plots)")

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
