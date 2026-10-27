#!/usr/bin/env python
"""Regression tests for the recipe gate (report 8.3 step 1).

The gate's whole value is that it FAILS when a field is absent or contradicts
the structure.  A test that only checks "a good recipe loads" would pass against
a module that ignores every field, so most of what follows is negative cases.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure

from recipe import (MissingRecipeField, Recipe, RecipeInconsistent, Validation,
                    load_pair)

PASS, FAIL = [], []


def check(name: str, fn) -> None:
    try:
        fn()
        PASS.append(name)
    except Exception as exc:                                  # noqa: BLE001
        FAIL.append((name, f"{type(exc).__name__}: {exc}"))


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------

B2 = Structure(Lattice.cubic(3.2246), ["Al", "Ta"],
               [[0, 0, 0], [0.5, 0.5, 0.5]])
L12 = Structure(Lattice.cubic(4.0389), ["Al", "Cr", "Cr", "Cr"],
                [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])

GOOD_RECIPE = {
    "structure_class": "BCC_B2",
    "site_occupancy": {"1a": "Al", "1b": "Ta"},
    "representation": "conventional",
    "sym_reduce": True,
    "supercell": [2, 2, 2],
    "kpoints": {"grid": [16, 16, 16], "shift": [0.5, 0.5, 0.5],
                "supercell_grid": [8, 8, 8],
                "supercell_shift": [0.5, 0.5, 0.5],
                "supercell_factor": 2, "cell_atoms": 2,
                "folding_verified": True,
                "n3V": 16 ** 3 * float(B2.volume)},
    "phonon_density": 16,
    "smearing": {"ISMEAR": 1, "SIGMA": 0.20},
    "displacement_angstrom": 0.01,
    "potcar": {"functional": "PBE",
               "variants": {"Nb": "Nb_pv", "default": "standard"}},
    "magnetic_policy": "nonmagnetic",
}
GOOD_VALIDATION = {
    "stage0_coverage": "measured",
    "measured_on": ["BCC-Al-Ta-2"],
    "not_measured": [],
    "qha_applicability": "not_yet_assessed",
    "recipe_robustness": "phonon_recipe_robust",
    "recorded_at": "2026-10-13",
}


def variant(recipe=None, validation=None):
    """Deep copy every time.

    `dict(GOOD_RECIPE, **overrides)` is a SHALLOW copy, so the nested dicts
    (`kpoints`, `smearing`, `potcar`) stay shared with the module-level fixture.
    The negative tests delete and overwrite keys inside those nested dicts, so
    with a shallow copy the FIRST test to touch `kpoints` empties it for every
    later test -- which is exactly what happened, and it looked like a module
    bug rather than a test bug.
    """
    import copy
    d = copy.deepcopy(GOOD_RECIPE)
    d.update(copy.deepcopy(recipe or {}))
    v = copy.deepcopy(GOOD_VALIDATION)
    v.update(copy.deepcopy(validation or {}))
    return {"recipe": d, "validation": v}


# --------------------------------------------------------------------------
# negative cases: missing fields must raise and NAME the field
# --------------------------------------------------------------------------

def _missing_raises(key: str):
    def fn():
        d = variant()
        del d["recipe"][key]
        try:
            Recipe.from_dict(d["recipe"])
        except MissingRecipeField as exc:
            assert key in exc.missing, f"raised but did not name {key}: {exc}"
            return
        raise AssertionError(f"missing {key} did NOT raise")
    return fn


for _k in ("structure_class", "site_occupancy", "representation", "sym_reduce",
           "supercell", "kpoints", "phonon_density", "smearing",
           "displacement_angstrom", "potcar", "magnetic_policy"):
    check(f"缺 recipe.{_k} → 报错并点名", _missing_raises(_k))

for _k in ("grid", "shift", "cell_atoms", "supercell_grid",
           "supercell_factor", "folding_verified", "n3V"):
    def _fn(k=_k):
        d = variant()
        del d["recipe"]["kpoints"][k]
        try:
            Recipe.from_dict(d["recipe"])
        except MissingRecipeField as exc:
            assert k in exc.missing, f"did not name kpoints.{k}: {exc}"
            return
        raise AssertionError(f"missing kpoints.{k} did NOT raise")
    check(f"缺 kpoints.{_k} → 报错并点名", _fn)

for _k in ("stage0_coverage", "measured_on", "not_measured",
           "qha_applicability", "recipe_robustness", "recorded_at"):
    def _fn(k=_k):
        d = dict(GOOD_VALIDATION)
        del d[k]
        try:
            Validation.from_dict(d)
        except MissingRecipeField as exc:
            assert k in exc.missing
            return
        raise AssertionError(f"missing validation.{k} did NOT raise")
    check(f"缺 validation.{_k} → 报错并点名", _fn)


# --------------------------------------------------------------------------
# negative cases: a recipe that parses but contradicts itself / the structure
# --------------------------------------------------------------------------

def _expect_inconsistent(recipe=None, validation=None, structure=B2,
                         needle=""):
    def fn():
        d = variant(recipe, validation)
        try:
            r, v = load_pair(d)
        except RecipeInconsistent as exc:
            assert needle in str(exc), f"raised but message lacks {needle!r}: {exc}"
            return
        probs = r.check_against(structure) + v.check_internal()
        assert probs, "no problem reported"
        assert any(needle in p for p in probs), f"problems lack {needle!r}: {probs}"
    return fn


check("supercell 与 kpoints.supercell_factor 不一致 → 报错",
      _expect_inconsistent({"supercell": [2, 3, 2]}, needle="supercell_factor"))
check("折叠配对破坏（supercell_grid × m ≠ grid）→ 报错",
      _expect_inconsistent({"kpoints": dict(GOOD_RECIPE["kpoints"],
                                            supercell_grid=[7, 8, 8])},
                           needle="fold-matched"))
check("grid 不能被 supercell_factor 整除 → 报错",
      _expect_inconsistent({"kpoints": dict(GOOD_RECIPE["kpoints"],
                                            grid=[15, 16, 16],
                                            n3V=15 ** 3 * float(B2.volume))},
                           needle="not divisible"))
check("网格排除 Γ → 报错",
      _expect_inconsistent({"kpoints": dict(GOOD_RECIPE["kpoints"],
                                            shift=[0.0, 0.5, 0.5])},
                           needle="excludes Gamma"))
check("超胞网格排除 Γ → 报错",
      _expect_inconsistent({"kpoints": dict(GOOD_RECIPE["kpoints"],
                                            supercell_grid=[8, 8, 8],
                                            supercell_shift=[0.0, 0.5, 0.5])},
                           needle="kpoints[supercell]"))
check("smearing 不是冻结值 → 报错",
      _expect_inconsistent({"smearing": {"ISMEAR": -5, "SIGMA": 0.2}},
                           needle="frozen"))
check("位移不是冻结值 → 报错",
      _expect_inconsistent({"displacement_angstrom": 0.05}, needle="0.01"))
check("POTCAR 变体与冻结选择不符 → 报错",
      _expect_inconsistent({"potcar": {"functional": "PBE",
                                       "variants": {"Nb": "Nb_sv"}}},
                           needle="Nb_pv"))
check("POTCAR 编造一个没有记录的例外 → 报错",
      _expect_inconsistent({"potcar": {"functional": "PBE",
                                       "variants": {"Ta": "Ta_pv"}}},
                           needle="no recorded exception"))
check("sym_reduce 写成字符串 → 报错",
      _expect_inconsistent({"sym_reduce": "true"}, needle="must be a bool"))
check("structure_class 不在词表 → 报错",
      _expect_inconsistent({"structure_class": "B2"},
                           needle="not in"))
check("structure_class=BCC_B2 但结构不是立方 → 报错",
      _expect_inconsistent({"structure_class": "BCC_B2",
                            "kpoints": dict(GOOD_RECIPE["kpoints"],
                                            cell_atoms=2)},
                           structure=Structure(
                               Lattice.hexagonal(3.2, 5.1), ["Al", "Ta"],
                               [[0, 0, 0], [1 / 3, 2 / 3, 0.5]]),
                           needle="expects a cubic"))
check("kpoints.n3V 与实际不符 → 报错",
      _expect_inconsistent({"kpoints": dict(GOOD_RECIPE["kpoints"],
                                            n3V=1.0)}, needle="n3V"))
check("kpoints.cell_atoms 与实际不符 → 报错",
      _expect_inconsistent({"kpoints": dict(GOOD_RECIPE["kpoints"],
                                            cell_atoms=4)},
                           needle="cell_atoms"))
check("L1₂ 缺 site_occupancy → 报错",
      _expect_inconsistent({"structure_class": "FCC_L12", "site_occupancy": {},
                            "kpoints": dict(GOOD_RECIPE["kpoints"],
                                            cell_atoms=4,
                                            n3V=16 ** 3 * float(L12.volume))},
                           structure=L12, needle="requires site_occupancy"))
check("L1₂ 占位写反（1a/3c 对调）→ 报错",
      _expect_inconsistent({"structure_class": "FCC_L12",
                            "site_occupancy": {"1a": "Cr",
                                               "3c": ["Al", "Al", "Al"]},
                            "kpoints": dict(GOOD_RECIPE["kpoints"],
                                            cell_atoms=4,
                                            n3V=16 ** 3 * float(L12.volume))},
                           structure=L12, needle="site_occupancy"))
check("validation: coverage=measured 但 measured_on 为空 → 报错",
      _expect_inconsistent({"kpoints": dict(GOOD_RECIPE["kpoints"])},
                           {"measured_on": [], "recipe_robustness": "recipe_untested"},
                           needle="empty measured_on"))
check("validation: coverage=untested 却有 measured_on → 报错",
      _expect_inconsistent({"kpoints": dict(GOOD_RECIPE["kpoints"])},
                           {"stage0_coverage": "untested",
                            "measured_on": ["X"], "recipe_robustness": "recipe_untested"},
                           needle="untested"))


# --------------------------------------------------------------------------
# positive cases
# --------------------------------------------------------------------------

def _good_loads():
    r, v = load_pair(variant())
    assert r.structure_class == "BCC_B2"
    assert v.stage0_coverage == "measured"
    assert not r.check_against(B2), r.check_against(B2)
check("完整 recipe + validation → 通过", _good_loads)


def _l12_ok():
    r, v = load_pair(variant(
        {"structure_class": "FCC_L12",
         "site_occupancy": {"1a": "Al", "3c": ["Cr", "Cr", "Cr"]},
         "kpoints": dict(GOOD_RECIPE["kpoints"], cell_atoms=4,
                         n3V=16 ** 3 * float(L12.volume))}))
    assert not r.check_against(L12), r.check_against(L12)
check("L1₂ 占位正确（Al@1a）→ 通过", _l12_ok)


def _l12_other_direction_ok():
    other = Structure(Lattice.cubic(4.0389), ["Cr", "Al", "Al", "Al"],
                      [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
    r, v = load_pair(variant(
        {"structure_class": "FCC_L12",
         "site_occupancy": {"1a": "Cr", "3c": ["Al", "Al", "Al"]},
         "kpoints": dict(GOOD_RECIPE["kpoints"], cell_atoms=4,
                         n3V=16 ** 3 * float(other.volume))}))
    assert not r.check_against(other), r.check_against(other)
check("L1₂ 另一个方向（Cr@1a）→ 也通过（两个方向都合法）", _l12_other_direction_ok)


def _two_axes_are_separate():
    """The whole point of report 7.8.3's second revision: an endpoint can be
    recipe-robust AND QHA-inapplicable at the same time."""
    v = Validation.from_dict(dict(
        GOOD_VALIDATION,
        qha_applicability="finite_q_soft_mode_diagnostic",
        recipe_robustness="phonon_recipe_robust"))
    assert v.qha_applicability != "harmonic_QHA_applicable"
    assert v.recipe_robustness == "phonon_recipe_robust"
    assert not v.check_internal(), v.check_internal()
check("qha_applicability 与 recipe_robustness 可同时取不同值",
      _two_axes_are_separate)


# --------------------------------------------------------------------------
# the real pool
# --------------------------------------------------------------------------

POOL = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json")


def _real_pool():
    if not POOL.is_file():
        print("  (跳过：pool_recipes.json 不存在)")
        return
    d = json.loads(POOL.read_text())
    assert d["n_endpoints"] == 84, d["n_endpoints"]
    assert d["n_problems"] == 0, d["problems"][:3]
    seen = {}
    for name, r in d["endpoints"].items():
        rr = Recipe.from_dict(r["recipe"])
        vv = Validation.from_dict(r["validation"])
        assert not rr.check_internal(), (name, rr.check_internal())
        assert not vv.check_internal(), (name, vv.check_internal())
        seen[rr.structure_class] = seen.get(rr.structure_class, 0) + 1
    assert sum(seen.values()) == 84, seen
    print(f"  真实池：84 端基全部通过，类别分布 {seen}")
check("真实池的 84 条 recipe 全部可加载且自洽", _real_pool)


# --------------------------------------------------------------------------

print()
for n in PASS:
    print(f"  PASS  {n}")
for n, e in FAIL:
    print(f"  FAIL  {n}\n        {e}")
print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
sys.exit(1 if FAIL else 0)
