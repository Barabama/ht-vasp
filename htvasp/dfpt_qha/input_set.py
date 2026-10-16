"""Recipe-backed VASP input generators for the DFPT-QHA stages.

Two stages, two generators:

``DfptRelaxSetGenerator``   ionic relaxation (``IBRION=2``) at a free (ISIF=3)
                            or fixed-volume (ISIF=2) cell, used by the EOS stage.
``DfptInputSetGenerator``   DFPT (``IBRION=8``) on the recipe's phonon
                            supercell, used by the phonon stage.

The recipe's settings are published through ``incar_updates``, the documented
extension point, so ``user_incar_settings`` stays what it says it is: user
overrides, applied on top.

POTCAR variants are deliberately NOT set.  ``atomate2``/``pymatgen`` resolve
the symbols from the selected functional; hard-coding a variant table was how a
dead ``"default": "standard"`` key silently produced ``Ta_pv`` for every
Ta-bearing endpoint.  The resolved symbols are what they are -- record them in
provenance rather than assume them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

from atomate2.vasp.sets.base import VaspInputGenerator
from pymatgen.io.vasp import Kpoints

from htvasp.model.dfpt_qha_recipe import Recipe

MAGMOM_BY_ELEMENT = {
    "Al": 1.0, "Cr": 3.0, "Hf": 1.0, "Nb": 3.0, "Ta": 3.0, "Ti": 3.0, "Zr": 1.0,
}

SHARED_INCAR: dict[str, Any] = {
    "ENCUT": 520,
    "PREC": "Normal",
    "ALGO": "Normal",
    "ISYM": 2,
    "LREAL": False,
    "LWAVE": False,
    "LCHARG": False,
    "ISTART": 0,
    "ICHARG": 2,
    "LORBIT": 11,
}

#: relaxation and static electronic settings (Stage 0A recipe)
RELAX_ELECTRONIC: dict[str, Any] = {
    "EDIFF": 1e-6, "NELM": 120, "NELMIN": 6, "NCORE": 8,
}

#: DFPT electronic settings, transcribed from the INCAR that actually ran.
#: NCORE and KPAR are mutually exclusive, so the DFPT family carries KPAR only.
DFPT_ELECTRONIC: dict[str, Any] = {
    "EDIFF": 1e-5, "NELM": 120, "NELMIN": 10, "KPAR": 2,
}

RELAX_MODE: dict[str, Any] = {"IBRION": 2, "NSW": 80, "POTIM": 0.2, "EDIFFG": -0.02}
DFPT_MODE: dict[str, Any] = {"IBRION": 8, "ISIF": 2, "NSW": 1}


class _RecipeGenerator(VaspInputGenerator):
    """Common plumbing: recipe-derived INCAR values and k-mesh selection."""

    recipe: Recipe | None = None
    use_supercell_kpoints: bool = False

    def _recipe_defaults(self) -> dict[str, Any]:
        recipe = self.recipe
        if recipe is None:
            raise ValueError(f"{type(self).__name__} requires an explicit Recipe")
        if recipe.magnetic_policy == "unresolved":
            raise ValueError(
                "magnetic_policy is 'unresolved': generating inputs would pick a "
                "spin state nobody decided"
            )
        return {
            **SHARED_INCAR,
            "ISMEAR": recipe.smearing.ISMEAR,
            "SIGMA": recipe.smearing.SIGMA,
            "ISPIN": 2,
            "MAGMOM": MAGMOM_BY_ELEMENT,
        }

    @property
    def kpoints_updates(self) -> Kpoints:
        recipe = self.recipe
        if recipe is None:
            raise ValueError(f"{type(self).__name__} requires an explicit Recipe")
        if self.use_supercell_kpoints:
            grid, shift = recipe.kpoints.supercell_grid, recipe.kpoints.supercell_shift
        else:
            grid, shift = recipe.kpoints.grid, recipe.kpoints.shift
        return Kpoints.monkhorst_automatic(kpts=tuple(grid), shift=tuple(shift))


@dataclass
class DfptRelaxSetGenerator(_RecipeGenerator):
    """Ionic relaxation at a free (ISIF=3) or fixed-volume (ISIF=2) cell."""

    recipe: Recipe | None = None
    isif: Literal[2, 3] = 2
    use_supercell_kpoints: bool = False
    user_incar_settings: dict[str, Any] = field(default_factory=dict)
    user_kpoints_settings: dict[str, Any] = field(default_factory=dict)
    user_potcar_settings: dict[str, str] = field(default_factory=dict)
    user_potcar_functional: str = "PBE"
    auto_ispin: bool = False
    auto_ismear: bool = False
    auto_metal_kpoints: bool = False

    @property
    def incar_updates(self) -> dict[str, Any]:
        return {
            **self._recipe_defaults(),
            **RELAX_ELECTRONIC,
            **RELAX_MODE,
            "ISIF": self.isif,
        }


@dataclass
class DfptInputSetGenerator(_RecipeGenerator):
    """VASP DFPT (``IBRION=8``) on the recipe's phonon supercell."""

    recipe: Recipe | None = None
    use_supercell_kpoints: bool = True
    user_incar_settings: dict[str, Any] = field(default_factory=dict)
    user_kpoints_settings: dict[str, Any] = field(default_factory=dict)
    user_potcar_settings: dict[str, str] = field(default_factory=dict)
    user_potcar_functional: str = "PBE"
    auto_ispin: bool = False
    auto_ismear: bool = False
    auto_metal_kpoints: bool = False

    @property
    def incar_updates(self) -> dict[str, Any]:
        return {**self._recipe_defaults(), **DFPT_ELECTRONIC, **DFPT_MODE}


__all__ = [
    "DFPT_ELECTRONIC",
    "DFPT_MODE",
    "DfptInputSetGenerator",
    "DfptRelaxSetGenerator",
    "MAGMOM_BY_ELEMENT",
    "RELAX_ELECTRONIC",
    "RELAX_MODE",
    "SHARED_INCAR",
]
