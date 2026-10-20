#!/usr/bin/env python
"""Emit the recipe + validation pair for all 84 endpoints.

This is report 8.3 step 1: the recipe is a REQUIRED input, and the coverage is a
separate record.  Everything here is derived from decisions that were actually
made, and every field that the module cannot fill from a decided value is left
absent so `Recipe.from_dict` raises -- rather than being papered over with a
plausible default.

Where each field comes from
---------------------------
structure_class      from the pool inventory (report 8.0d)
site_occupancy       from the cell itself, for L1_2 (the class alone does not
                     say which direction -- report 4.12.3)
representation       "conventional" for all 84 (report 8.0d)
sym_reduce           True for the ordered classes; SER_HCP is recorded but the
                     validation says its symmetry was never audited
supercell            2x2x2 for all (user decision)
kpoints              mesh_for(V_cell, 16), even + shift 0.5; supercell mesh
                     = cell/2 with the Gamma-parity rule (report 4.8.1)
phonon_density       16 (report 4.11.10, "Stage 0 baseline" -- the validation
                     says it was only measured on three BCC/B2 representatives)
smearing / displacement / potcar   frozen (reports 4.11.3, 4.11.5, 5.1)
magnetic_policy      from a measurement (MEASURED_MAGNETIC): the three
                     original Al-Cr results plus the 18-endpoint spin screen
                     of 2026-10-20; "unresolved" if a Cr-bearing endpoint has
                     none; "nonmagnetic" otherwise
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from pymatgen.core import Structure

import build_stage0 as S0
from recipe import Recipe, Validation, load_pair

DENSITY = 16
SUPERCELL = 2
RECORDED_AT = "2026-10-13"

# Magnetic policy.  The value is "confirmed_branch" only where a calculation
# actually showed a magnetic ground state, per report 7.8.3's rule about not
# letting an untested claim look like a setting.
MEASURED_MAGNETIC = {
    "BCC-Al-Cr-2": ("confirmed_branch",
                    "report 4.9.3: magnetic ground state 283.7 meV/atom below "
                    "non-magnetic"),
    "FCC-Al-Cr-4": ("confirmed_branch",
                    "report 8.0a(j): Cr=3.0 branch 22.96 meV/atom below the "
                    "zero-moment one"),
    "FCC-Cr-Al-4": ("confirmed_branch",
                    "report 8.0a(j): Cr=3.0 branch 183.89 meV/atom below the "
                    "zero-moment one"),
    # Spin screen 2026-10-20: free ISIF=3 relaxations from nm / fm / afm
    # starts, 5 meV/atom threshold, AFM states excluded where the cell changed
    # by more than 2%.  data/endpoints/<name>/spin-screen.json holds the
    # energies, moments and per-state lattice.
    "BCC-Cr-Cr-2": ("confirmed_branch",
                    "spin screen 2026-10-20: afm 14.5 meV/atom below "
                    "non-magnetic (|m|=1.13, cell preserved)"),
    "SER-Cr-1": ("confirmed_branch",
                 "spin screen 2026-10-20: same 2-atom bcc cell as "
                 "BCC-Cr-Cr-2, afm 14.5 meV/atom below non-magnetic"),
    "BCC-Cr-Hf-2": ("confirmed_branch",
                    "spin screen 2026-10-20: fm 77.5 meV/atom below "
                    "non-magnetic"),
    "BCC-Cr-Ti-2": ("confirmed_branch",
                    "spin screen 2026-10-20: fm 21.4 meV/atom below "
                    "non-magnetic"),
    "BCC-Cr-Zr-2": ("confirmed_branch",
                    "spin screen 2026-10-20: fm 91.5 meV/atom below "
                    "non-magnetic"),
    "FCC-Cr-Hf-4": ("confirmed_branch",
                    "spin screen 2026-10-20: fm 102.9 meV/atom below "
                    "non-magnetic"),
    "FCC-Cr-Nb-4": ("confirmed_branch",
                    "spin screen 2026-10-20: fm 63.8 meV/atom below "
                    "non-magnetic"),
    "FCC-Cr-Zr-4": ("confirmed_branch",
                    "spin screen 2026-10-20: fm 128.4 meV/atom below "
                    "non-magnetic"),
    "FCC-Hf-Cr-4": ("confirmed_branch",
                    "spin screen 2026-10-20: fm 15.7 meV/atom below "
                    "non-magnetic (free afm collapsed the cell)"),
    "FCC-Nb-Cr-4": ("confirmed_branch",
                    "spin screen 2026-10-20: fm 6.5 meV/atom below "
                    "non-magnetic; constrained afm at the fm geometry is "
                    "0.3 meV/atom above fm and the nm start itself converges "
                    "magnetic"),
    "FCC-Ta-Cr-4": ("confirmed_branch",
                    "spin screen 2026-10-20: fm 5.8 meV/atom below "
                    "non-magnetic; constrained afm at the fm geometry is "
                    "1.0 meV/atom above fm and nm is 0.2 above fm"),
    "FCC-Zr-Cr-4": ("confirmed_branch",
                    "spin screen 2026-10-20: fm 41.1 meV/atom below "
                    "non-magnetic (free afm collapsed the cell)"),
    "BCC-Cr-Nb-2": ("nonmagnetic",
                    "spin screen 2026-10-20: the fm start converges back to "
                    "the non-magnetic state (gain 0.0 meV/atom)"),
    "BCC-Cr-Ta-2": ("nonmagnetic",
                    "spin screen 2026-10-20: fm is 7.0 meV/atom above "
                    "non-magnetic"),
    "FCC-Cr-Cr-4": ("nonmagnetic",
                    "spin screen 2026-10-20: fm and afm both converge to zero "
                    "moment (gain 0.0 meV/atom)"),
    "FCC-Cr-Ta-4": ("nonmagnetic",
                    "spin screen 2026-10-20: the fm start converges to zero "
                    "moment (gain 0.0 meV/atom)"),
    "FCC-Cr-Ti-4": ("nonmagnetic",
                    "spin screen 2026-10-20: fm is 2.9 meV/atom below "
                    "non-magnetic, inside the 5 meV/atom threshold"),
    "FCC-Ti-Cr-4": ("nonmagnetic",
                    "spin screen 2026-10-20: fm is 2.9 meV/atom below "
                    "non-magnetic, inside the 5 meV/atom threshold (the free "
                    "afm start collapsed the cell)"),
}

# Per-class validation.  `measured_on` names what was ACTUALLY run; anything the
# recipe contains but the runs did not cover goes in `not_measured`.
CLASS_VALIDATION = {
    "BCC_B2": {
        "stage0_coverage": "measured",
        "measured_on": ["BCC-Ta-Ta-2", "BCC-Al-Ta-2", "BCC-Al-Cr-2",
                        "BCC-Hf-V-2"],
        "not_measured": ["含 Cr 端基的自旋筛选（只有 BCC-Al-Cr-2 测过）"],
        "qha_applicability": "not_yet_assessed",
        "recipe_robustness": "phonon_recipe_robust",
    },
    "FCC_L12": {
        "stage0_coverage": "partially_measured",
        "measured_on": ["FCC-Zr-Al-4", "FCC-Al-Cr-4", "FCC-Cr-Al-4"],
        "not_measured": ["声子层（Stage 0C 只跑了单体积、非 d=16）",
                         "含 Cr 之外元素组合的 E(V)",
                         "d=16 声子密度"],
        "qha_applicability": "not_yet_assessed",
        "recipe_robustness": "recipe_untested",
    },
    "SER_FCC": {
        "stage0_coverage": "partially_measured",
        "measured_on": ["SER-Al-1"],
        "not_measured": ["声子层", "conventional 路径的声子超胞（Stage 0C 走 primitive）"],
        "qha_applicability": "not_yet_assessed",
        "recipe_robustness": "recipe_untested",
    },
    "SER_BCC": {
        "stage0_coverage": "partially_measured",
        "measured_on": ["SER-Ta-1"],
        "not_measured": ["声子层", "SER-Nb-1 完全没跑", "SER-Cr-1 的自旋态"],
        "qha_applicability": "not_yet_assessed",
        "recipe_robustness": "recipe_untested",
    },
    "SER_HCP": {
        "stage0_coverage": "partially_measured",
        "measured_on": ["SER-Zr-1"],
        "not_measured": ["声子层", "SER-Ti-1 / SER-Hf-1 完全没跑",
                         "HCP 的对称性约化（c/a 与内坐标）从未审计"],
        "qha_applicability": "not_yet_assessed",
        "recipe_robustness": "recipe_untested",
    },
    "LOW_SYMMETRY": {
        "stage0_coverage": "untested",
        "measured_on": [],
        "not_measured": ["全部：池内没有低对称端基，这一类是为将来预留的"],
        "qha_applicability": "not_yet_assessed",
        "recipe_robustness": "recipe_untested",
    },
}

# The 5 BCC aluminides that report 4.2 directly audited as soft, plus the
# AlTa endpoint that the mode-following work carried further.
SOFT_MODE_ENDPOINTS = {
    "BCC-Al-Al-2": "report 4.2: mesh minimum -3.3983 THz",
    "BCC-Al-Ta-2": "report 4.2: mesh minimum -3.6234 THz; report 4.13.13d: "
                   "ISIF=2 relaxation leaves B2 (-130.4 meV/atom)",
    "BCC-Al-Ti-2": "report 4.2: mesh minimum -3.6543 THz",
    "BCC-Ta-Ti-2": "report 4.2: mesh minimum -1.9407 THz",
    "BCC-Ti-Ti-2": "report 4.2: mesh minimum -4.4343 THz",
}


def class_of(rec: dict, structure: Structure) -> str:
    if rec["class"] == "BCC_B2":
        return "BCC_B2"
    if rec["class"] == "FCC_L12":
        return "FCC_L12"
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    cs = SpacegroupAnalyzer(structure, symprec=1e-3).get_crystal_system()
    return {"cubic": ("SER_FCC" if _is_fcc(structure) else "SER_BCC"),
            "hexagonal": "SER_HCP"}.get(cs, "LOW_SYMMETRY")


def _is_fcc(structure: Structure) -> bool:
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    return SpacegroupAnalyzer(structure, symprec=1e-3).get_space_group_symbol() \
        in ("Fm-3m", "F-43m")


def occupancy_for(rec: dict, structure: Structure) -> dict:
    """Record which element sits on which Wyckoff site.

    L1_2 NEEDS this (the two directions are different compounds, report
    4.12.3).  B2 does not need it to tell endpoints apart -- AB == BA -- but
    recording it costs nothing and lets the checker catch a chemistry swap.
    """
    if rec["class"] not in ("FCC_L12", "BCC_B2"):
        return {}
    import numpy as np
    at_1a, rest = None, []
    for s in structure:
        f = np.array(s.frac_coords) % 1.0
        if np.allclose(f, 0.0, atol=1e-4) or np.allclose(f, 1.0, atol=1e-4):
            at_1a = str(s.specie)
        else:
            rest.append(str(s.specie))
    if rec["class"] == "FCC_L12":
        return {"1a": at_1a, "3c": rest}
    return {"1a": at_1a, "1b": rest[0] if rest else None}


def kpoints_for(structure: Structure) -> dict:
    n_cell = S0.mesh_for(float(structure.volume), DENSITY)
    assert n_cell % SUPERCELL == 0, (n_cell, SUPERCELL)
    n_sup = n_cell // SUPERCELL
    cell_shift = 0.5 if n_cell % 2 == 0 else 0.0
    sup_shift = 0.5 if n_sup % 2 == 0 else 0.0
    return {"grid": [n_cell] * 3, "shift": [cell_shift] * 3,
            "supercell_grid": [n_sup] * 3, "supercell_shift": [sup_shift] * 3,
            "supercell_factor": SUPERCELL, "cell_atoms": len(structure),
            "folding_verified": True,
            "n3V": n_cell ** 3 * float(structure.volume)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True, type=Path,
                    help="pool_structures.json from build_pool_structures.py")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    pool = json.loads(args.pool.read_text())
    eps = pool["endpoints"]
    out: dict = {}
    problems: list[str] = []

    for name, rec in eps.items():
        st = Structure.from_dict(rec["structure"])
        cls = class_of(rec, st)
        mag, mag_why = "nonmagnetic", None
        # `any(list.count(...))` is a TypeError -- count returns an int, and
        # any() needs an iterable.  Membership is what was meant.
        if "Cr" in rec.get("elements", []):
            mag, mag_why = MEASURED_MAGNETIC.get(
                name, ("unresolved", "含 Cr 但未做自旋筛选"))
        elif rec["class"] in ("BCC_B2", "FCC_L12"):
            mag_why = ("无 Cr；记录为 nonmagnetic 但未做自旋对照 —— "
                       "这是默认预期，不是验证过的结论")
        recipe = {
            "structure_class": cls,
            "site_occupancy": occupancy_for(rec, st),
            "representation": "conventional",
            "sym_reduce": True,
            "supercell": [SUPERCELL] * 3,
            "kpoints": kpoints_for(st),
            "phonon_density": DENSITY,
            "smearing": {"ISMEAR": 1, "SIGMA": 0.20},
            "displacement_angstrom": 0.01,
            "potcar": {"functional": "PBE",
                       "variants": {"Nb": "Nb_pv", "default": "standard"}},
            "magnetic_policy": mag,
        }
        v = CLASS_VALIDATION[cls]
        validation = dict(v, recorded_at=RECORDED_AT)
        # per-endpoint override: an endpoint that was directly audited as soft
        # says so, regardless of what its class default is.
        if name in SOFT_MODE_ENDPOINTS:
            validation = dict(validation,
                              qha_applicability="finite_q_soft_mode_diagnostic")
            validation["not_measured"] = list(validation["not_measured"]) + [
                f"软模直接证据：{SOFT_MODE_ENDPOINTS[name]}"]
        validation["measured_on"] = list(validation["measured_on"])
        validation["not_measured"] = list(validation["not_measured"])

        try:
            r, vv = load_pair({"recipe": recipe, "validation": validation})
        except Exception as exc:                              # noqa: BLE001
            problems.append(f"{name}: {exc}")
            continue
        # recipe vs the actual structure
        probs = r.check_against(st)
        if probs:
            problems.append(f"{name}: " + "; ".join(probs))
        out[name] = {"physical_state_id": rec.get("physical_state_id"),
                     "n_atoms": len(st),
                     "recipe": r.to_dict(), "validation": vv.to_dict(),
                     "magnetic_note": mag_why,
                     "vpa_is_proxy": rec.get("is_proxy")}

    summary = {"schema": "legacy-qha-pool/pool-recipes/v1",
               "n_endpoints": len(out),
               "n_problems": len(problems),
               "problems": problems,
               "endpoints": out}

    print(f"{'class':>14} {'n':>4} {'coverage':>20} {'qha_applicability':>30}")
    print("-" * 74)
    by = {}
    for name, r in out.items():
        by.setdefault(r["recipe"]["structure_class"], []).append(r)
    for cls, rs in sorted(by.items()):
        cov = sorted({r["validation"]["stage0_coverage"] for r in rs})
        app = sorted({r["validation"]["qha_applicability"] for r in rs})
        print(f"{cls:>14} {len(rs):>4} {','.join(cov):>20} {','.join(app):>30}")
    print()
    mag_counts = {}
    for r in out.values():
        mag_counts[r["recipe"]["magnetic_policy"]] = \
            mag_counts.get(r["recipe"]["magnetic_policy"], 0) + 1
    print("magnetic_policy 分布:", mag_counts)
    print()
    if problems:
        print(f"!! {len(problems)} 个端基没通过：")
        for p in problems[:10]:
            print("   ", p)
    else:
        print(f"全部 {len(out)} 个端基：recipe 必填齐全 + 结构一致性检查通过")

    if args.out:
        args.out.write_text(json.dumps(summary, indent=2) + "\n")
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
