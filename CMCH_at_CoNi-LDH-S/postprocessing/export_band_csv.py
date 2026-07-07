#!/usr/bin/env python3
"""导出能带结构逐点数据 → output/band_curves/

每体系输出一个 CSV: k_distance, band1_up, band2_up, ..., band1_down, band2_down, ...

用法:
  python postprocessing/export_band_csv.py
"""

import json
import logging
import os
from pathlib import Path

import numpy as np
from monty.json import MontyDecoder
from pymatgen.electronic_structure.bandstructure import Spin

from pymatgen.electronic_structure.bandstructure import Spin

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

DATA_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/data")
OUTPUT_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/postprocessing/output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SYSTEMS = [
    "CoNiOH2", "CoNiOH2S-noH", "CoMnH2CO5",
    "CMCH_strained", "LDH_strained", "LDH_S_strained", "LDH_S_flip_strained", "LDH2_strained",
    "hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed",
]


def export_band_structure(name: str):
    path = DATA_DIR / name / "band_out.json"
    if not path.exists():
        log.warning("  ✗ %s: no band_out.json", name)
        return False

    try:
        with open(path) as f:
            raw = json.load(f)
        bd_dict = raw.get("vasp_objects", {}).get("bandstructure")
        if bd_dict is None:
            log.warning("  ✗ %s: no bandstructure object", name)
            return False

        bs = MontyDecoder().process_decoded(bd_dict)
        n_kpts = len(bs.kpoints)

        # ── k 点路径标注 ───────────────────────────────────────────
        k_distances = list(bs.distance)
        branches_arr = [""] * n_kpts
        labels_arr = [""] * n_kpts
        branch_meta = []

        for b in bs.branches:
            start, end = b["start_index"], b["end_index"]
            # 部分体系的 end_index 可能 = n_kpts-1 (末点)
            actual_end = min(end, n_kpts - 1)
            bname = b.get("name", "?")

            for ik in range(start, actual_end + 1):
                branches_arr[ik] = bname

            # 端点标签（从 kpoints 对象自带 label 取，保底用分支名解析）
            kp_s = bs.kpoints[start]
            kp_e = bs.kpoints[actual_end]
            sl = getattr(kp_s, "label", "") or bname.split("-")[0]
            el = getattr(kp_e, "label", "") or bname.split("-")[-1]
            labels_arr[start] = sl
            labels_arr[actual_end] = el

            branch_meta.append({
                "name": bname,
                "start_index": start,
                "end_index": actual_end,
                "start_label": sl,
                "end_label": el,
            })

        # ── 能带本征值（保持与 k_distances 一致的顺序） ───────────
        eigenvalues = []
        spin_channels = list(bs.bands.keys())

        for spin in spin_channels:
            for band_idx in range(bs.nb_bands):
                eigvals = []
                for b in bs.branches:
                    start = b["start_index"]
                    end = min(b["end_index"], n_kpts - 1)
                    for ik in range(start, end + 1):
                        eigvals.append(bs.bands[spin][band_idx][ik])
                eigenvalues.append(eigvals)

        # ── 写出 CSV ──────────────────────────────────────────────
        os.makedirs(OUTPUT_DIR / "band_curves", exist_ok=True)
        fpath = OUTPUT_DIR / "band_curves" / f"band_{name}.csv"

        n_spin = len(spin_channels)
        n_band = bs.nb_bands
        spin_labels = {Spin.up: "up", Spin.down: "down"}

        header_parts = ["k_index", "k_distance", "branch", "label"]
        for s in spin_channels:
            for b in range(n_band):
                header_parts.append(f"band_{b+1}_{spin_labels.get(s, s)}")

        kpt_indices = list(range(n_kpts))
        data = np.column_stack([kpt_indices, k_distances,
                                branches_arr, labels_arr] + eigenvalues)

        # np.savetxt 混合 str/float 需 fmt="%s"
        np.savetxt(fpath, data, delimiter=",", header=",".join(header_parts),
                   comments="", fmt="%s")
        log.info("  ✓ %s: %s kpts × %s bands × %s spins → %s",
                 name, n_kpts, n_band, len(spin_channels), fpath.name)

        # ── 写出分支结构 JSON（画图用） ──────────────────────────
        branch_path = OUTPUT_DIR / "band_curves" / f"band_{name}_branches.json"
        with open(branch_path, "w") as f:
            json.dump({"branches": branch_meta}, f, indent=2)
        log.info("  ✓ %s: branches info → %s", name, branch_path.name)

        return True

    except Exception as e:
        log.warning("  ✗ %s: %s", name, e)
        return False


def main():
    log.info("Exporting band structure curves for %d systems ...", len(SYSTEMS))
    ok = sum(1 for name in SYSTEMS if export_band_structure(name))
    log.info("Done: %d/%d systems exported", ok, len(SYSTEMS))
    print(f"\n  Output: {OUTPUT_DIR}/band_curves/*.csv")


if __name__ == "__main__":
    main()
