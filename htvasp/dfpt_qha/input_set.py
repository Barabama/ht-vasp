"""Recipe-backed VASP input generator for one DFPT job."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from atomate2.vasp.sets.base import VaspInputGenerator
from pymatgen.io.vasp import Kpoints

from htvasp.model.dfpt_qha_recipe import Recipe

MAGMOM_BY_ELEMENT = {"Al": 1.0, "Cr": 3.0, "Hf": 1.0, "Nb": 3.0, "Ta": 3.0, "Ti": 3.0, "Zr": 1.0}


@dataclass
class DfptInputSetGenerator(VaspInputGenerator):
    recipe: Recipe | None = None
    user_incar_settings: dict[str, Any] = field(default_factory=dict)
    user_kpoints_settings: dict[str, Any] = field(default_factory=dict)
    user_potcar_settings: dict[str, str] = field(default_factory=dict)
    user_potcar_functional: str = "PBE"
    stage: str = "dfpt"
    auto_ispin: bool = False
    auto_ismear: bool = False
    auto_metal_kpoints: bool = False

    def __post_init__(self) -> None:
        if self.stage != "dfpt":
            raise ValueError("DfptInputSetGenerator is only for the DFPT stage")
        if self.recipe is None:
            raise ValueError("DfptInputSetGenerator requires an explicit Recipe")
        if self.recipe.magnetic_policy == "unresolved":
            raise ValueError("unresolved magnetic policy cannot generate VASP inputs")
        self.user_potcar_settings.update(self.recipe.potcar.variants)
        self.user_incar_settings.update(self._recipe_incar)
        super().__post_init__()

    @property
    def _recipe_incar(self) -> dict[str, Any]:
        return {
            "ENCUT": 520, "PREC": "Normal", "ALGO": "Normal", "ISYM": 2,
            "LREAL": False, "LWAVE": False, "LCHARG": False,
            "ISTART": 0, "ICHARG": 2, "LORBIT": 11,
            "IBRION": 8, "ISIF": 2, "NSW": 1,
            "EDIFF": 1e-5, "NELM": 120, "NELMIN": 10, "KPAR": 2,
            "ISMEAR": self.recipe.smearing.ISMEAR,
            "SIGMA": self.recipe.smearing.SIGMA,
            "ISPIN": 2,
            "MAGMOM": MAGMOM_BY_ELEMENT,
        }

    @property
    def incar_updates(self) -> dict[str, Any]:
        return dict(self._recipe_incar)

    @property
    def kpoints_updates(self) -> Kpoints:
        return Kpoints.monkhorst_automatic(
            kpts=tuple(self.recipe.kpoints.grid),
            shift=tuple(self.recipe.kpoints.shift),
        )


__all__ = ["DfptInputSetGenerator", "MAGMOM_BY_ELEMENT"]
