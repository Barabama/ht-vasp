"""DFPT-QHA VASP input generators.

These classes deliberately follow atomate2's ``VaspInputGenerator`` contract:
they describe an input set; ``BaseVaspMaker`` writes and runs it.  They never
write files, invoke VASP, or submit to a scheduler themselves.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from atomate2.vasp.sets.base import VaspInputGenerator

from htvasp.model.dfpt_qha_recipe import Recipe


DFPT_SHARED_INCAR = {
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

STATIC_ELECTRONIC = {"EDIFF": 1e-6, "NELM": 120, "NELMIN": 6, "NCORE": 8}
DFPT_ELECTRONIC = {"EDIFF": 1e-5, "NELM": 120, "NELMIN": 10, "KPAR": 2}


@dataclass
class DfptQhaInputSetGenerator(VaspInputGenerator):
    """Base generator for one recipe-backed DFPT-QHA VASP job."""

    recipe: Recipe | None = None
    stage: str = "static"
    density: int | None = None
    user_incar_settings: dict[str, Any] = field(default_factory=dict)
    user_kpoints_settings: dict[str, Any] = field(default_factory=dict)
    user_potcar_functional: str | None = "PBE"
    auto_ispin: bool = False
    auto_ismear: bool = False
    auto_metal_kpoints: bool = False

    def __post_init__(self) -> None:
        if self.stage not in {"relax", "static", "dfpt"}:
            raise ValueError("stage must be 'relax', 'static', or 'dfpt'")
        if self.recipe is not None and self.stage == "dfpt":
            if self.recipe.displacement_angstrom != 0.01:
                raise ValueError("DFPT recipe displacement must be the frozen 0.01 A")
        super().__post_init__()

    @property
    def incar_updates(self) -> dict[str, Any]:
        if self.stage == "relax":
            stage = {"IBRION": 2, "ISIF": 3, "NSW": 80, "POTIM": 0.2, "EDIFFG": -0.02}
            electronic = STATIC_ELECTRONIC
        elif self.stage == "dfpt":
            # This is the production DFPT path.  Do not replace with IBRION=2.
            stage = {"IBRION": 8, "ISIF": 2, "NSW": 1}
            electronic = DFPT_ELECTRONIC
        else:
            stage = {"IBRION": -1, "ISIF": 2, "NSW": 0}
            electronic = STATIC_ELECTRONIC
        policy = self.recipe.magnetic_policy if self.recipe is not None else "confirmed_branch"
        updates = {**DFPT_SHARED_INCAR, **electronic, **stage, "ISMEAR": 1, "SIGMA": 0.20,
                    "ISPIN": 2 if policy != "unresolved" else 2}
        if policy == "unresolved":
            raise ValueError("cannot generate inputs for unresolved magnetic_policy")
        return {**updates, **self.user_incar_settings}

    @property
    def kpoints_updates(self) -> dict[str, Any]:
        if self.recipe is not None:
            grid = tuple(self.recipe.kpoints.grid)
            shift = tuple(self.recipe.kpoints.shift)
            if len(set(grid)) != 1 or len(set(shift)) != 1:
                raise ValueError("DFPT-QHA currently requires isotropic recipe kpoints")
            return {"grid_density": grid[0], "shift": shift}
        return dict(self.user_kpoints_settings)

    @property
    def potcar_updates(self) -> dict[str, Any]:
        return {}


@dataclass
class DfptRelaxSetGenerator(DfptQhaInputSetGenerator):
    """Recipe-backed initial or fixed-volume relaxation generator."""

    stage: str = "relax"


@dataclass
class DfptStaticSetGenerator(DfptQhaInputSetGenerator):
    """Recipe-backed static energy generator."""

    stage: str = "static"


@dataclass
class DfptSetGenerator(DfptQhaInputSetGenerator):
    """Recipe-backed VASP DFPT generator using ``IBRION=8``."""

    stage: str = "dfpt"


__all__ = [
    "DFPT_ELECTRONIC",
    "DFPT_SHARED_INCAR",
    "DfptQhaInputSetGenerator",
    "DfptRelaxSetGenerator",
    "DfptSetGenerator",
    "DfptStaticSetGenerator",
    "STATIC_ELECTRONIC",
]
