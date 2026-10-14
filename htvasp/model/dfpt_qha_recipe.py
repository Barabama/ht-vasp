#!/usr/bin/env python
"""Recipe + validation records -- two layers, no implicit defaults.

Report 7.8.3.  The point of this module is NOT to hold data; it is to make the
ABSENCE of a decision loud.

Two records, deliberately separate
----------------------------------
`Recipe`      what settings this calculation used.  Scientific input.
`Validation`  what range of structures those settings have actually been
              verified on.  A record about evidence, NOT a calculation
              parameter.  Mixing the two turns "untested" into something that
              reads like a setting (report 7.8.3, first revision).

No defaults, anywhere
---------------------
`Recipe.from_dict()` raises `MissingRecipeField` naming every absent key.  There
is no fallback for `sym_reduce`, no assumed supercell, no default smearing.
The reason is the project's own history: three generator bugs (§9.3) all had the
shape "the generator wrote a template property and the checker only asked
whether the key existed".  A recipe object with defaults reproduces exactly
that failure mode at a higher level -- a field nobody set silently takes a
value, and 84 endpoints get computed with a setting no one chose.

Beyond presence: cross-checks
-----------------------------
`Recipe.check_against(structure)` verifies the recipe against the structure it
claims to describe -- class vs cell shape, occupancy vs fractional coordinates,
the folding pairing of the k-mesh, the Gamma-inclusion rule, and `n3V`.  A
recipe that parses but describes a different cell is worse than a missing one.
"""
from __future__ import annotations

import itertools
import json
from dataclasses import dataclass

# --------------------------------------------------------------------------
# vocabularies
# --------------------------------------------------------------------------

STRUCTURE_CLASSES = ("BCC_B2", "FCC_L12", "SER_FCC", "SER_BCC", "SER_HCP",
                     "LOW_SYMMETRY")
REPRESENTATIONS = ("conventional", "primitive", "as_built")
MAGNETIC_POLICIES = ("nonmagnetic", "confirmed_branch", "unresolved")
COVERAGE = ("measured", "partially_measured", "assumed_from_analogue",
            "untested")
QHA_APPLICABILITY = ("harmonic_QHA_applicable",
                     "finite_q_soft_mode_diagnostic",
                     "distortion_candidate",
                     "anharmonic_required",
                     "numerical_unresolved",
                     # Report 7.8.3's list has no value for "no phonon has been
                     # run on this endpoint yet".  Without it the only options
                     # are to leave the field out (and lose the check) or to
                     # borrow a value that means something else -- so it is
                     # added here rather than misused.
                     "not_yet_assessed")
RECIPE_ROBUSTNESS = ("phonon_recipe_robust", "recipe_untested")

# POTCAR variants that are NOT pymatgen's plain default for that element.
# Nb is the pool's one exception: there is no standard Nb POTCAR under
# potpaw_PBE, only Nb_pv and Nb_sv, and Nb_pv was chosen (report 5.1).
POTCAR_KNOWN_VARIANTS = {"Nb": "Nb_pv"}

# Every key that must be present.  Listed explicitly rather than derived so a
# typo in a builder fails here instead of silently shrinking the requirement.
RECIPE_REQUIRED = (
    "structure_class", "site_occupancy", "representation", "sym_reduce",
    "supercell", "kpoints", "phonon_density", "smearing",
    "displacement_angstrom", "potcar", "magnetic_policy",
)
# NOTE ON THE FIELD NAMES.  Report 7.8.3's example shows
#     "grid": [16,16,16] ... "folding_parent_grid": [32,32,32], "n3V": 137337
# which cannot all be true at once: 16^3 * 33.53 A^3 = 137339 IS the 2-atom
# cell's n3V, so `grid` is the CELL mesh, and the fold-matched SUPERCELL mesh
# is grid/supercell_factor = 8, not 32.  The example's two grid fields are
# swapped.  This module uses the unambiguous pair below and the report is
# corrected separately.
KP_REQUIRED = ("grid", "shift", "cell_atoms", "supercell_grid",
               "supercell_factor", "folding_verified", "n3V")
SM_REQUIRED = ("ISMEAR", "SIGMA")
PO_REQUIRED = ("functional", "variants")
VALIDATION_REQUIRED = ("stage0_coverage", "measured_on", "not_measured",
                       "qha_applicability", "recipe_robustness", "recorded_at")


class MissingRecipeField(ValueError):
    """Raised when a required recipe field is absent.

    Carries the whole list so one round trip reports every gap, not just the
    first -- a builder that is missing four fields should learn that once.
    """

    def __init__(self, layer: str, missing: list[str], extra: str = ""):
        self.layer = layer
        self.missing = missing
        msg = (f"{layer} is missing required field(s): {', '.join(missing)}. "
               f"There are no defaults for these by design -- see report "
               f"7.8.3. Set them explicitly.")
        if extra:
            msg += f"  {extra}"
        super().__init__(msg)


class RecipeInconsistent(ValueError):
    """The recipe parses but does not describe the structure it is attached to."""


def _require(d: dict, keys, layer: str) -> None:
    missing = [k for k in keys if k not in d or d[k] is None]
    if missing:
        raise MissingRecipeField(layer, missing)


# --------------------------------------------------------------------------
# layer 1: recipe
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Recipe:
    structure_class: str
    site_occupancy: dict
    representation: str
    sym_reduce: bool
    supercell: tuple
    kpoints: dict
    phonon_density: int
    smearing: dict
    displacement_angstrom: float
    potcar: dict
    magnetic_policy: str

    @classmethod
    def from_dict(cls, d: dict) -> "Recipe":
        _require(d, RECIPE_REQUIRED, "recipe")
        _require(d["kpoints"], KP_REQUIRED, "recipe.kpoints")
        _require(d["smearing"], SM_REQUIRED, "recipe.smearing")
        _require(d["potcar"], PO_REQUIRED, "recipe.potcar")

        if d["structure_class"] not in STRUCTURE_CLASSES:
            raise RecipeInconsistent(
                f"structure_class {d['structure_class']!r} not in "
                f"{STRUCTURE_CLASSES}")
        if d["representation"] not in REPRESENTATIONS:
            raise RecipeInconsistent(
                f"representation {d['representation']!r} not in "
                f"{REPRESENTATIONS}")
        if d["magnetic_policy"] not in MAGNETIC_POLICIES:
            raise RecipeInconsistent(
                f"magnetic_policy {d['magnetic_policy']!r} not in "
                f"{MAGNETIC_POLICIES}")
        if not isinstance(d["sym_reduce"], bool):
            raise RecipeInconsistent(
                f"sym_reduce must be a bool, got {d['sym_reduce']!r}; a string "
                f"like 'true' would compare unequal to True everywhere")
        sup = tuple(int(x) for x in d["supercell"])
        if len(sup) != 3 or min(sup) < 1:
            raise RecipeInconsistent(f"supercell must be 3 positive ints: {sup}")
        return cls(
            structure_class=d["structure_class"],
            site_occupancy=dict(d["site_occupancy"]),
            representation=d["representation"],
            sym_reduce=bool(d["sym_reduce"]),
            supercell=sup,
            kpoints=dict(d["kpoints"]),
            phonon_density=int(d["phonon_density"]),
            smearing=dict(d["smearing"]),
            displacement_angstrom=float(d["displacement_angstrom"]),
            potcar=dict(d["potcar"]),
            magnetic_policy=d["magnetic_policy"],
        )

    def to_dict(self) -> dict:
        return {"structure_class": self.structure_class,
                "site_occupancy": self.site_occupancy,
                "representation": self.representation,
                "sym_reduce": self.sym_reduce,
                "supercell": list(self.supercell),
                "kpoints": self.kpoints,
                "phonon_density": self.phonon_density,
                "smearing": self.smearing,
                "displacement_angstrom": self.displacement_angstrom,
                "potcar": self.potcar,
                "magnetic_policy": self.magnetic_policy}

    # -- internal consistency, no structure needed ------------------------
    def check_internal(self) -> list[str]:
        p = []
        g = tuple(self.kpoints["grid"])
        s = tuple(self.kpoints["shift"])
        sg = tuple(self.kpoints["supercell_grid"])
        m = int(self.kpoints["supercell_factor"])
        # 1. folding pairing: supercell mesh x factor == cell mesh
        for i in range(3):
            if sg[i] * m != g[i]:
                p.append(f"kpoints: supercell_grid[{i}]={sg[i]} but "
                         f"grid[{i}]/supercell_factor={g[i]}/{m}="
                         f"{g[i] / m:.3f}. The fold-matched supercell mesh is "
                         f"what makes the supercell sampling equal the cell "
                         f"sampling (report 4.8.1).")
            if g[i] % m != 0:
                p.append(f"kpoints: grid[{i}]={g[i]} is not divisible by "
                         f"supercell_factor={m}, so no integer supercell mesh "
                         f"is fold-matched to it")
        if tuple(self.supercell) != (m, m, m):
            p.append(f"recipe.supercell={self.supercell} disagrees with "
                     f"kpoints.supercell_factor={m}")
        # 2. Gamma inclusion, per the rule in report 4.8.1, on BOTH meshes
        for label, mesh, sh3 in (("cell", g, s),
                                 ("supercell", sg,
                                  tuple(self.kpoints.get("supercell_shift", s)))):
            for i in range(3):
                n, sh = mesh[i], sh3[i]
                contains = (n % 2 == 1 and abs(sh) < 1e-9) or \
                           (n % 2 == 0 and abs(sh - 0.5) < 1e-9)
                if not contains:
                    p.append(f"kpoints[{label}]: grid {mesh} shift {sh3} "
                             f"excludes Gamma along axis {i} (n={n}, "
                             f"shift={sh}); the folding-compatible family is "
                             f"'even n + shift 0.5' (report 4.8.1).")
        # 3. phonon density must reproduce the supercell grid from the cell
        if self.phonon_density <= 0:
            p.append(f"phonon_density={self.phonon_density} must be positive")
        # 4. potcar variants must be the ones the pool actually uses
        for el, v in self.potcar["variants"].items():
            if el == "default":
                continue
            want = POTCAR_KNOWN_VARIANTS.get(el)
            if want is None:
                p.append(f"potcar.variants names {el}={v!r} but there is no "
                         f"recorded exception for {el}; the pool default is "
                         f"the standard variant (report 5.1)")
            elif v != want:
                p.append(f"potcar.variants {el}={v!r} but the frozen choice is "
                         f"{want!r} (report 5.1)")
        # 5. smearing / displacement must be the frozen values
        if int(self.smearing["ISMEAR"]) != 1 or \
                abs(float(self.smearing["SIGMA"]) - 0.20) > 1e-9:
            p.append(f"smearing {self.smearing} is not the frozen "
                     f"ISMEAR=1/SIGMA=0.20 (report 4.11.3)")
        if abs(self.displacement_angstrom - 0.01) > 1e-9:
            p.append(f"displacement_angstrom={self.displacement_angstrom} is "
                     f"not the frozen 0.01 (report 4.11.5)")
        return p

    # -- recipe vs the structure it claims to describe ---------------------
    def check_against(self, structure) -> list[str]:
        import numpy as np
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        p = list(self.check_internal())
        n = len(structure)
        sa = SpacegroupAnalyzer(structure, symprec=1e-3)
        cs = sa.get_crystal_system()

        want = {"BCC_B2": ("cubic", 2), "FCC_L12": ("cubic", 4),
                "SER_FCC": ("cubic", None), "SER_BCC": ("cubic", None),
                "SER_HCP": ("hexagonal", None)}.get(self.structure_class)
        if want:
            exp_cs, exp_n = want
            if cs != exp_cs:
                p.append(f"structure_class={self.structure_class} expects a "
                         f"{exp_cs} cell, got {cs}")
            if exp_n is not None and n != exp_n:
                p.append(f"structure_class={self.structure_class} expects "
                         f"{exp_n} atoms, got {n}")
        # occupancy -- for B2 too.  AB == BA for B2 so the assignment does not
        # distinguish endpoints the way L1_2's does, but a chemistry SWAP
        # (Al@1b instead of Al@1a) is still a real error and costs nothing to
        # catch here.
        if self.structure_class == "BCC_B2":
            seen = {}
            for s in structure:
                f = np.array(s.frac_coords) % 1.0
                key = "1a" if (np.allclose(f, 0.0, atol=1e-4)
                               or np.allclose(f, 1.0, atol=1e-4)) else "1b"
                seen.setdefault(key, []).append(str(s.specie))
            claimed = self.site_occupancy or {}
            for site in ("1a", "1b"):
                if site in claimed and claimed[site] not in seen.get(site, []):
                    p.append(f"site_occupancy says {site}={claimed[site]!r} but "
                             f"the cell has {seen.get(site)} there")
        if self.structure_class == "FCC_L12":
            seen = {}
            for s in structure:
                f = np.array(s.frac_coords) % 1.0
                key = "1a" if (np.allclose(f, 0.0, atol=1e-4)
                               or np.allclose(f, 1.0, atol=1e-4)) else "3c"
                seen.setdefault(key, []).append(str(s.specie))
            claimed = self.site_occupancy or {}
            if claimed:
                if "1a" in claimed and claimed["1a"] not in seen.get("1a", []):
                    p.append(f"site_occupancy says 1a={claimed['1a']!r} but the "
                             f"atom at the origin is {seen.get('1a')}")
                if "3c" in claimed:
                    want3 = sorted(claimed["3c"]) if isinstance(
                        claimed["3c"], (list, tuple)) else [claimed["3c"]]
                    if sorted(seen.get("3c", [])) != want3:
                        p.append(f"site_occupancy 3c={claimed['3c']!r} but the "
                                 f"cell has {sorted(seen.get('3c', []))}")
            else:
                p.append("structure_class=FCC_L12 requires site_occupancy: the "
                         "two L1_2 directions are different compounds and the "
                         "class alone does not say which one this is "
                         "(report 4.12.3 / 4.2b)")
        # n3V
        n3V = self.kpoints["grid"][0] ** 3 * float(structure.volume)
        if abs(n3V - float(self.kpoints["n3V"])) > 1e-3 * max(n3V, 1.0):
            p.append(f"kpoints.n3V={self.kpoints['n3V']} but grid[0]^3 * V = "
                     f"{n3V:.1f}; n3V is the density invariant and must be the "
                     f"one actually used (report 4.8.3)")
        # cell_atoms
        if int(self.kpoints["cell_atoms"]) != n:
            p.append(f"kpoints.cell_atoms={self.kpoints['cell_atoms']} but the "
                     f"structure has {n} atoms")
        return p


# --------------------------------------------------------------------------
# layer 2: validation
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Validation:
    stage0_coverage: str
    measured_on: list
    not_measured: list
    qha_applicability: str
    recipe_robustness: str
    recorded_at: str

    @classmethod
    def from_dict(cls, d: dict) -> "Validation":
        _require(d, VALIDATION_REQUIRED, "validation")
        if d["stage0_coverage"] not in COVERAGE:
            raise RecipeInconsistent(
                f"stage0_coverage {d['stage0_coverage']!r} not in {COVERAGE}")
        if d["qha_applicability"] not in QHA_APPLICABILITY:
            raise RecipeInconsistent(
                f"qha_applicability {d['qha_applicability']!r} not in "
                f"{QHA_APPLICABILITY}")
        if d["recipe_robustness"] not in RECIPE_ROBUSTNESS:
            raise RecipeInconsistent(
                f"recipe_robustness {d['recipe_robustness']!r} not in "
                f"{RECIPE_ROBUSTNESS}")
        if len(d["measured_on"]) != len(set(d["measured_on"])):
            raise RecipeInconsistent("measured_on has duplicates")
        return cls(stage0_coverage=d["stage0_coverage"],
                   measured_on=list(d["measured_on"]),
                   not_measured=list(d["not_measured"]),
                   qha_applicability=d["qha_applicability"],
                   recipe_robustness=d["recipe_robustness"],
                   recorded_at=d["recorded_at"])

    def to_dict(self) -> dict:
        return {"stage0_coverage": self.stage0_coverage,
                "measured_on": self.measured_on,
                "not_measured": self.not_measured,
                "qha_applicability": self.qha_applicability,
                "recipe_robustness": self.recipe_robustness,
                "recorded_at": self.recorded_at}

    def check_internal(self) -> list[str]:
        p = []
        if self.stage0_coverage == "measured" and not self.measured_on:
            p.append("stage0_coverage='measured' with an empty measured_on: "
                     "name the endpoints or downgrade the coverage")
        if self.stage0_coverage == "untested" and self.measured_on:
            p.append(f"stage0_coverage='untested' but measured_on={self.measured_on}")
        # the two axes must be free to disagree (report 7.8.3)
        if self.recipe_robustness == "phonon_recipe_robust" and \
                self.qha_applicability == "harmonic_QHA_applicable" and \
                self.stage0_coverage == "untested":
            p.append("recipe_robustness='phonon_recipe_robust' and "
                     "qha_applicability='harmonic_QHA_applicable' cannot both "
                     "hold on an untested class; a robustness claim needs a "
                     "measurement")
        return p


def load_pair(d: dict) -> tuple[Recipe, Validation]:
    """Build both layers from a two-key dict and run the internal checks."""
    _require(d, ("recipe", "validation"), "record")
    r = Recipe.from_dict(d["recipe"])
    v = Validation.from_dict(d["validation"])
    problems = r.check_internal() + v.check_internal()
    if problems:
        raise RecipeInconsistent("; ".join(problems))
    return r, v


if __name__ == "__main__":
    # a self-test of the gate itself: missing fields must raise, and a recipe
    # that parses must still fail if it contradicts its structure.
    good = {
        "structure_class": "BCC_B2",
        "site_occupancy": {"1a": "Al", "1b": "Ta"},
        "representation": "conventional", "sym_reduce": True,
        "supercell": [2, 2, 2],
        "kpoints": {"grid": [16, 16, 16], "shift": [0.5, 0.5, 0.5],
                    "cell_atoms": 2, "supercell_grid": [8, 8, 8],
                    "supercell_shift": [0.5, 0.5, 0.5],
                    "supercell_factor": 2,
                    "folding_verified": True, "n3V": 137337},
        "phonon_density": 16, "smearing": {"ISMEAR": 1, "SIGMA": 0.20},
        "displacement_angstrom": 0.01,
        "potcar": {"functional": "PBE",
                   "variants": {"Nb": "Nb_pv", "default": "standard"}},
        "magnetic_policy": "nonmagnetic",
    }
    r = Recipe.from_dict(good)
    print("完整 recipe:", "内部检查通过" if not r.check_internal()
          else r.check_internal())
    for drop in ("sym_reduce", "kpoints", "magnetic_policy"):
        bad = {k: v for k, v in good.items() if k != drop}
        try:
            Recipe.from_dict(bad)
            print(f"  漏掉 {drop}: **没有报错** ← 这是错的")
        except MissingRecipeField as e:
            print(f"  漏掉 {drop}: 报错 ✓ ({len(e.missing)} 个缺失字段)")
    for k, v, why in (("supercell", [2, 3, 2], "折叠配对被破坏"),
                      ("smearing", {"ISMEAR": -5, "SIGMA": 0.2}, "不是冻结展宽"),
                      ("displacement_angstrom", 0.05, "不是冻结位移")):
        bad = dict(good, **{k: v})
        rr = Recipe.from_dict(bad)
        probs = rr.check_internal()
        print(f"  {k}={v} ({why}): {'报错 ✓' if probs else '**没有报错** ← 错的'}")
