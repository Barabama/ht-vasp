"""Explicit DFPT-QHA recipe and model contracts.

This module contains data and validation only.  It does not import jobflow,
launch VASP, or write files.  Calculations consume these models through the
DFPT input generator and preserve the serialized recipe in TaskDoc metadata.
"""

from __future__ import annotations

from typing import Any, Literal, Mapping

from pydantic import BaseModel, ConfigDict, Field, model_validator

RECIPE_SCHEMA = "htvasp-dfpt-qha-recipe/v1"


class RecipeError(ValueError):
    """Base class for invalid or incomplete recipe data."""


class MissingRecipeField(RecipeError):
    def __init__(self, layer: str, missing: list[str]):
        self.layer, self.missing = layer, missing
        super().__init__(f"{layer} missing required field(s): {', '.join(missing)}")


class RecipeInconsistent(RecipeError):
    pass


class KpointsRecipe(BaseModel):
    model_config = ConfigDict(extra="forbid")

    grid: tuple[int, int, int]
    shift: tuple[float, float, float]
    cell_atoms: int = Field(gt=0)
    supercell_grid: tuple[int, int, int]
    supercell_shift: tuple[float, float, float]
    supercell_factor: int = Field(gt=0)
    folding_verified: bool
    n3V: float = Field(gt=0)

    @model_validator(mode="after")
    def validate_folding(self) -> "KpointsRecipe":
        for a, b in zip(self.grid, self.supercell_grid):
            if a != b * self.supercell_factor:
                raise ValueError("supercell_grid is not fold-matched to grid")
        return self


class SmearingRecipe(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ISMEAR: int
    SIGMA: float


class PotcarRecipe(BaseModel):
    model_config = ConfigDict(extra="forbid")

    functional: str
    variants: dict[str, str]


class Recipe(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal[RECIPE_SCHEMA] = RECIPE_SCHEMA
    structure_class: str
    site_occupancy: dict[str, Any]
    representation: Literal["conventional", "primitive", "as_built"]
    sym_reduce: bool
    supercell: tuple[int, int, int]
    kpoints: KpointsRecipe
    phonon_density: int = Field(gt=0)
    smearing: SmearingRecipe
    displacement_angstrom: float = 0.01
    potcar: PotcarRecipe
    magnetic_policy: Literal["nonmagnetic", "confirmed_branch", "unresolved"]

    @model_validator(mode="after")
    def validate_recipe(self) -> "Recipe":
        if self.supercell != (self.kpoints.supercell_factor,) * 3:
            raise ValueError("recipe supercell and kpoints factor disagree")
        if self.smearing.ISMEAR != 1 or abs(self.smearing.SIGMA - 0.20) > 1e-12:
            raise ValueError("recipe smearing is not the frozen ISMEAR=1/SIGMA=0.20")
        if abs(self.displacement_angstrom - 0.01) > 1e-12:
            raise ValueError("recipe displacement is not the frozen 0.01 Angstrom")
        return self

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Recipe":
        try:
            return cls.model_validate(data)
        except Exception as exc:
            if hasattr(exc, "errors"):
                missing = [
                    ".".join(str(x) for x in e["loc"])
                    for e in exc.errors() if e["type"] == "missing"
                ]
                if missing:
                    raise MissingRecipeField("recipe", missing) from exc
            raise RecipeInconsistent(str(exc)) from exc


class Validation(BaseModel):
    model_config = ConfigDict(extra="forbid")

    stage0_coverage: Literal["measured", "partially_measured", "assumed_from_analogue", "untested"]
    measured_on: list[str]
    not_measured: list[str]
    qha_applicability: str
    recipe_robustness: Literal["phonon_recipe_robust", "recipe_untested"]
    recorded_at: str


def load_pair(data: Mapping[str, Any]) -> tuple[Recipe, Validation]:
    missing = [key for key in ("recipe", "validation") if key not in data]
    if missing:
        raise MissingRecipeField("record", missing)
    return Recipe.from_dict(data["recipe"]), Validation.model_validate(data["validation"])


__all__ = [
    "KpointsRecipe", "MissingRecipeField", "PotcarRecipe", "Recipe", "RecipeError",
    "RecipeInconsistent", "RECIPE_SCHEMA", "SmearingRecipe", "Validation", "load_pair",
]
