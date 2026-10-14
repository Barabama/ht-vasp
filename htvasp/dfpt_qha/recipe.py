"""The endpoint recipe: scientific input that must be stated, never assumed.

Two records, deliberately separate
----------------------------------
``Recipe``      the settings an endpoint calculation uses.  Scientific input.
``Validation``  the range of structures those settings have actually been
                verified on.  A record about EVIDENCE, not a calculation
                parameter.  Merging the two is what turns "untested" into
                something that reads like a setting, so they are two models and
                neither implies anything about the other.

No defaults
-----------
Every field is required and both models forbid extra keys.  There is no
fallback ``sym_reduce``, no assumed supercell, no default smearing.  The
project's own history is the reason: three generator incidents all had the
shape "the builder wrote a template value and the checker only asked whether
the key existed".  A recipe with defaults reproduces that failure one level
up -- a field nobody chose silently takes a value, and the whole pool is
computed with a setting no one picked.

A recipe that cannot be constructed cannot be run
-------------------------------------------------
Self-consistency (the folding pairing, the Gamma rule, the frozen values) is
enforced when the model is built, so a ``Recipe`` that exists is one that can
be handed to a builder.  What a model cannot know is the structure it is
supposed to describe, so that half is the separate ``check_against`` method.

The comparisons here are exact on purpose.  ``supercell_factor`` is an int and
``supercell_grid`` is an int, so "the supercell mesh is fold-matched to the
cell mesh" is an equality, not a tolerance -- a mesh that is *nearly* matched
is a mesh that samples a different set of q-points.
"""

from __future__ import annotations

from datetime import date
from typing import Any, Literal, Mapping, get_args

from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

SCHEMA_VERSION = "htvasp-dfpt-qha-recipe/v1"

# --------------------------------------------------------------------------
# vocabularies
#
# Each is a Literal so the type checker and the runtime agree, with the tuple
# derived from it -- a second hand-written list is a second thing to forget to
# update, and the gap between the two is exactly what a closed vocabulary is
# supposed to close.
# --------------------------------------------------------------------------

StructureClass = Literal[
    "BCC_B2", "FCC_L12", "SER_FCC", "SER_BCC", "SER_HCP", "LOW_SYMMETRY"
]
Representation = Literal["conventional", "primitive", "as_built"]
MagneticPolicy = Literal["nonmagnetic", "confirmed_branch", "unresolved"]
Stage0Coverage = Literal[
    "measured", "partially_measured", "assumed_from_analogue", "untested"
]
QhaApplicability = Literal[
    "harmonic_QHA_applicable",
    "finite_q_soft_mode_diagnostic",
    "distortion_candidate",
    "anharmonic_required",
    "numerical_unresolved",
    # Without this value the only ways to describe an endpoint whose phonons
    # have not been run are to leave the field out -- losing the check -- or to
    # borrow a value that means something else.  It is added rather than
    # misused.
    "not_yet_assessed",
]
RecipeRobustness = Literal["phonon_recipe_robust", "recipe_untested"]

STRUCTURE_CLASSES = get_args(StructureClass)
REPRESENTATIONS = get_args(Representation)
MAGNETIC_POLICIES = get_args(MagneticPolicy)
STAGE0_COVERAGES = get_args(Stage0Coverage)
QHA_APPLICABILITIES = get_args(QhaApplicability)
RECIPE_ROBUSTNESSES = get_args(RecipeRobustness)

# --------------------------------------------------------------------------
# frozen values
#
# These are not preferences.  Each was fixed by a Stage 0 measurement, and the
# validators below reject anything else rather than accepting a "variant" that
# would silently make an endpoint incomparable with its neighbours.
# --------------------------------------------------------------------------

# Stage 0B: the force constants are indistinguishable across four widths
# (0.2-4.3% spread), so one width is used for the whole pool.
FROZEN_SMEARING: dict[str, Any] = {"ISMEAR": 1, "SIGMA": 0.20}
# Stage 0B: 0.01 A sits inside the harmonic regime for every representative.
FROZEN_DISPLACEMENT_ANGSTROM = 0.01
FROZEN_POTCAR_FUNCTIONAL = "PBE"

# POTCAR variants that are NOT pymatgen's plain default for that element.  Nb
# is the pool's only exception: potpaw_PBE has no standard Nb POTCAR, only
# Nb_pv and Nb_sv, and Nb_pv was chosen.  An element absent from this table
# has no recorded exception, so naming one for it is an error, not a choice.
POTCAR_KNOWN_VARIANTS = {"Nb": "Nb_pv"}

# ``structure_class`` -> (crystal system, required atom count or None).
# Classes absent from this table are not constrained here.
CLASS_EXPECTATIONS: dict[str, tuple[str, int | None]] = {
    "BCC_B2": ("cubic", 2),
    "FCC_L12": ("cubic", 4),
    "SER_FCC": ("cubic", None),
    "SER_BCC": ("cubic", None),
    "SER_HCP": ("hexagonal", None),
}

# Gamma is in a Monkhorst-Pack mesh iff every n_i is odd with shift 0, or every
# n_i is even with shift 0.5.  A mixed-parity mesh with shift 0.5 is not
# Gamma-centred -- a common silent error, so it is checked on BOTH meshes.
GAMMA_SHIFT = 0.5
_PARITY_TOL = 1e-9


class RecipeError(ValueError):
    """Base class for a recipe that cannot be used."""


class MissingRecipeField(RecipeError):
    """A required field is absent.

    ``missing`` holds fully-qualified dotted paths (``"kpoints.grid"``) and
    carries every gap at once, so a builder that is missing four fields learns
    that in one round trip instead of four.
    """

    def __init__(self, layer: str, missing: list[str]) -> None:
        self.layer = layer
        self.missing = sorted(missing)
        super().__init__(
            f"{layer} is missing required field(s): {', '.join(self.missing)}. "
            "There are no defaults for these by design -- set them explicitly."
        )


class RecipeInconsistent(RecipeError):
    """The recipe parses but contradicts itself or the structure it describes."""

    def __init__(self, problems: str | list[str]) -> None:
        self.problems = [problems] if isinstance(problems, str) else list(problems)
        super().__init__("; ".join(self.problems))


def _translate(exc: ValidationError, layer: str) -> RecipeError:
    """Turn pydantic's report into the domain's failure language.

    Only an all-missing report becomes ``MissingRecipeField``; anything else
    means a value that was present and wrong, which is a different problem and
    gets a different type.  Callers that want to react to "nobody filled this
    in yet" must not have to string-match pydantic's formatting to do it.
    """

    errors = exc.errors()
    if all(error["type"] == "missing" for error in errors):
        return MissingRecipeField(
            layer, [".".join(str(part) for part in error["loc"]) for error in errors]
        )
    return RecipeInconsistent(
        [
            f"{'.'.join(str(part) for part in error['loc'])}: {error['msg']}"
            for error in errors
        ]
    )


def axis_contains_gamma(n: int, shift: float) -> bool:
    """Whether Gamma sits on this one axis of an automatic mesh.

    Odd ``n`` with shift 0, or even ``n`` with shift 0.5.  A mixed-parity mesh
    with shift 0.5 is NOT Gamma-centred -- a common silent error.
    """

    if abs(shift) < _PARITY_TOL:
        return n % 2 == 1
    if abs(shift - GAMMA_SHIFT) < _PARITY_TOL:
        return n % 2 == 0
    return False


def gamma_centred(mesh: tuple[int, ...], shift: tuple[float, ...]) -> bool:
    """Whether Gamma is in the mesh, i.e. on EVERY axis.

    This is the single definition of the rule; the input gate imports it
    rather than restating it.
    """

    return all(axis_contains_gamma(n, s) for n, s in zip(mesh, shift))


def _gamma_problems(
    label: str, mesh: tuple[int, ...], shift: tuple[float, ...]
) -> list[str]:
    problems: list[str] = []
    for axis, (n, s) in enumerate(zip(mesh, shift)):
        if not axis_contains_gamma(n, s):
            problems.append(
                f"kpoints[{label}]: axis {axis} has grid {n} shift {s}, which "
                f"excludes Gamma (n={n} is {'even' if n % 2 == 0 else 'odd'}); "
                f"the folding-compatible family is 'even n + shift 0.5'"
            )
    return problems


class Kpoints(BaseModel):
    """The cell mesh and the phonon-supercell mesh, and the pairing between them.

    The pairing is what makes the supercell sampling equal the cell sampling:
    ``supercell_grid[i] * supercell_factor == grid[i]``.  Both meshes must
    independently contain Gamma, and ``n3V = grid[0] ** 3 * V`` is the density
    invariant that has to match the cell actually used -- it is carried
    explicitly rather than recomputed so that a mismatch is visible instead of
    being absorbed.

    Names are deliberately ``grid`` (the CELL mesh) and ``supercell_grid``.  An
    earlier draft of this contract called the cell mesh's fold-parent
    ``folding_parent_grid`` and gave it a value that contradicted ``n3V``; the
    two spellings here cannot be swapped without the arithmetic failing.
    """

    model_config = ConfigDict(extra="forbid")

    # cell
    grid: tuple[int, int, int]
    shift: tuple[float, float, float]
    cell_atoms: int = Field(gt=0)
    # supercell
    supercell_grid: tuple[int, int, int]
    supercell_shift: tuple[float, float, float]
    supercell_factor: int = Field(gt=0)
    folding_verified: bool = Field(strict=True)
    # density invariant
    n3V: float = Field(gt=0)

    @model_validator(mode="after")
    def _folding_pairing(self) -> "Kpoints":
        problems: list[str] = []
        for label, mesh in (("grid", self.grid), ("supercell_grid", self.supercell_grid)):
            if any(n < 1 for n in mesh):
                problems.append(f"kpoints.{label} {list(mesh)} must be all positive")
        for axis in range(3):
            cell, supercell = self.grid[axis], self.supercell_grid[axis]
            if supercell * self.supercell_factor != cell:
                quotient = cell / self.supercell_factor
                problems.append(
                    f"kpoints.supercell_grid[{axis}]={supercell} but "
                    f"grid[{axis}]/supercell_factor={cell}/{self.supercell_factor}"
                    f"={quotient:.3f}. The fold-matched supercell mesh is what "
                    f"makes the supercell sampling equal the cell sampling"
                )
        if problems:
            raise ValueError("; ".join(problems))
        return self

    @model_validator(mode="after")
    def _gamma_on_both_meshes(self) -> "Kpoints":
        problems: list[str] = []
        for label, mesh, shift in (
            ("cell", self.grid, self.shift),
            ("supercell", self.supercell_grid, self.supercell_shift),
        ):
            problems.extend(_gamma_problems(label, mesh, shift))
        if problems:
            raise ValueError("; ".join(problems))
        return self


class Smearing(BaseModel):
    """Electronic smearing, frozen for the whole pool.

    Not a per-endpoint choice: a Stage 0B calibration found the force
    constants indistinguishable across widths spanning 0.2-4.3%, so one width
    is used everywhere.  Any other value is a mistake rather than a variant,
    because it would make the endpoint's energies incomparable with its
    neighbours.
    """

    model_config = ConfigDict(extra="forbid")

    ISMEAR: int
    SIGMA: float

    @model_validator(mode="after")
    def _frozen(self) -> "Smearing":
        if self.ISMEAR != FROZEN_SMEARING["ISMEAR"] or abs(
            self.SIGMA - FROZEN_SMEARING["SIGMA"]
        ) > _PARITY_TOL:
            raise ValueError(
                f"smearing ISMEAR={self.ISMEAR} SIGMA={self.SIGMA} is not the "
                f"frozen {FROZEN_SMEARING}"
            )
        return self


class Potcar(BaseModel):
    """The potential set: one functional, plus the pool's recorded exceptions."""

    model_config = ConfigDict(extra="forbid")

    functional: str
    variants: dict[str, str]

    @model_validator(mode="after")
    def _frozen_functional(self) -> "Potcar":
        if self.functional != FROZEN_POTCAR_FUNCTIONAL:
            raise ValueError(
                f"potcar.functional={self.functional!r} is not the frozen "
                f"{FROZEN_POTCAR_FUNCTIONAL!r}"
            )
        return self

    @model_validator(mode="after")
    def _recorded_variants_only(self) -> "Potcar":
        problems: list[str] = []
        for element, variant in self.variants.items():
            if element == "default":
                continue
            recorded = POTCAR_KNOWN_VARIANTS.get(element)
            if recorded is None:
                problems.append(
                    f"potcar.variants names {element}={variant!r} but there is "
                    f"no recorded exception for {element}; the default is the "
                    f"standard variant"
                )
            elif variant != recorded:
                problems.append(
                    f"potcar.variants {element}={variant!r} but the frozen "
                    f"choice is {recorded!r}"
                )
        if problems:
            raise ValueError("; ".join(problems))
        return self


class Recipe(BaseModel):
    """What settings one endpoint calculation uses.

    Constructing one proves it is internally consistent.  ``check_against``
    proves it describes the structure it is attached to, which is the half no
    amount of self-consistency can establish -- a recipe that parses but
    describes a different cell is worse than a missing one.
    """

    model_config = ConfigDict(extra="forbid")

    structure_class: StructureClass
    site_occupancy: dict[str, Any]
    representation: Representation
    # strict: a recipe written with the string "true" would compare unequal to
    # True everywhere downstream while looking correct in the file.
    sym_reduce: bool = Field(strict=True)
    supercell: tuple[int, int, int]
    kpoints: Kpoints
    phonon_density: int = Field(gt=0)
    smearing: Smearing
    displacement_angstrom: float
    potcar: Potcar
    magnetic_policy: MagneticPolicy

    @model_validator(mode="after")
    def _supercell_matches_factor(self) -> "Recipe":
        if self.supercell != (self.kpoints.supercell_factor,) * 3:
            raise ValueError(
                f"recipe.supercell={list(self.supercell)} disagrees with "
                f"kpoints.supercell_factor={self.kpoints.supercell_factor}"
            )
        return self

    @model_validator(mode="after")
    def _frozen_displacement(self) -> "Recipe":
        if abs(self.displacement_angstrom - FROZEN_DISPLACEMENT_ANGSTROM) > _PARITY_TOL:
            raise ValueError(
                f"displacement_angstrom={self.displacement_angstrom} is not the "
                f"frozen {FROZEN_DISPLACEMENT_ANGSTROM}"
            )
        return self

    @model_validator(mode="after")
    def _l12_needs_occupancy(self) -> "Recipe":
        # The two L1_2 directions are different compounds and the class alone
        # does not say which one this is, so the occupancy is not optional.
        if self.structure_class == "FCC_L12" and not self.site_occupancy:
            raise ValueError(
                "structure_class=FCC_L12 requires site_occupancy: the two L1_2 "
                "directions are different compounds (A@1a/B@3c vs B@1a/A@3c)"
            )
        return self

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Recipe":
        """Build a recipe, reporting gaps in the domain's own terms."""

        try:
            return cls.model_validate(data)
        except ValidationError as exc:
            raise _translate(exc, "recipe") from exc

    def check_against(self, structure: Any) -> list[str]:
        """Problems with this recipe as a description of ``structure``.

        Returns a list rather than raising so a batch builder can report every
        endpoint's problems in one pass instead of stopping at the first.
        """

        import numpy as np
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        problems: list[str] = []
        n_atoms = len(structure)
        crystal_system = SpacegroupAnalyzer(structure, symprec=1e-3).get_crystal_system()

        expectation = CLASS_EXPECTATIONS.get(self.structure_class)
        if expectation is not None:
            expected_system, expected_atoms = expectation
            if crystal_system != expected_system:
                problems.append(
                    f"structure_class={self.structure_class} expects a "
                    f"{expected_system} cell, got {crystal_system}"
                )
            if expected_atoms is not None and n_atoms != expected_atoms:
                problems.append(
                    f"structure_class={self.structure_class} expects "
                    f"{expected_atoms} atoms, got {n_atoms}"
                )

        if self.structure_class in ("BCC_B2", "FCC_L12"):
            problems.extend(self._occupancy_problems(structure, np))

        n3V = self.kpoints.grid[0] ** 3 * float(structure.volume)
        if abs(n3V - float(self.kpoints.n3V)) > 1e-3 * max(n3V, 1.0):
            problems.append(
                f"kpoints.n3V={self.kpoints.n3V} but grid[0]^3 * V = {n3V:.1f}; "
                f"n3V is the density invariant and must be the one actually used"
            )
        if self.kpoints.cell_atoms != n_atoms:
            problems.append(
                f"kpoints.cell_atoms={self.kpoints.cell_atoms} but the structure "
                f"has {n_atoms} atoms"
            )
        return problems

    def _occupancy_problems(self, structure: Any, np: Any) -> list[str]:
        """Check the declared 1a/1b or 1a/3c assignment against the cell.

        Recorded for B2 as well as L1_2.  AB == BA for B2, so the assignment
        does not distinguish endpoints the way L1_2's does, but a chemistry
        swap is still a real error and costs nothing to catch here.
        """

        seen: dict[str, list[str]] = {}
        for site in structure:
            frac = np.array(site.frac_coords) % 1.0
            at_origin = np.allclose(frac, 0.0, atol=1e-4) or np.allclose(
                frac, 1.0, atol=1e-4
            )
            key = "1a" if at_origin else ("1b" if self.structure_class == "BCC_B2" else "3c")
            seen.setdefault(key, []).append(str(site.specie))

        claimed = self.site_occupancy or {}
        if self.structure_class == "BCC_B2":
            problems = []
            for site in ("1a", "1b"):
                if site in claimed and claimed[site] not in seen.get(site, []):
                    problems.append(
                        f"site_occupancy says {site}={claimed[site]!r} but the "
                        f"cell has {seen.get(site)} there"
                    )
            return problems

        problems = []
        if "1a" in claimed and claimed["1a"] not in seen.get("1a", []):
            problems.append(
                f"site_occupancy says 1a={claimed['1a']!r} but the atom at the "
                f"origin is {seen.get('1a')}"
            )
        if "3c" in claimed:
            declared = claimed["3c"]
            wanted = sorted(declared) if isinstance(declared, (list, tuple)) else [declared]
            if sorted(seen.get("3c", [])) != wanted:
                problems.append(
                    f"site_occupancy 3c={declared!r} but the cell has "
                    f"{sorted(seen.get('3c', []))}"
                )
        return problems


class Validation(BaseModel):
    """What range of structures a recipe has actually been verified on.

    Two independent axes on purpose.  ``qha_applicability`` says whether the
    harmonic approximation holds for this endpoint; ``recipe_robustness`` says
    whether the settings were tested at all.  An endpoint can be
    recipe-robust AND QHA-inapplicable at the same time -- collapsing the two
    into one field is what made "untested" read like a setting.
    """

    model_config = ConfigDict(extra="forbid")

    stage0_coverage: Stage0Coverage
    measured_on: list[str]
    not_measured: list[str]
    qha_applicability: QhaApplicability
    recipe_robustness: RecipeRobustness
    recorded_at: date

    @model_validator(mode="after")
    def _coverage_matches_evidence(self) -> "Validation":
        problems: list[str] = []
        if len(self.measured_on) != len(set(self.measured_on)):
            problems.append(f"measured_on has duplicates: {self.measured_on}")
        if self.stage0_coverage == "measured" and not self.measured_on:
            problems.append(
                "stage0_coverage='measured' with an empty measured_on: name the "
                "endpoints or downgrade the coverage"
            )
        if self.stage0_coverage == "untested" and self.measured_on:
            problems.append(
                f"stage0_coverage='untested' but measured_on={self.measured_on}"
            )
        if (
            self.recipe_robustness == "phonon_recipe_robust"
            and self.qha_applicability == "harmonic_QHA_applicable"
            and self.stage0_coverage == "untested"
        ):
            problems.append(
                "recipe_robustness='phonon_recipe_robust' and "
                "qha_applicability='harmonic_QHA_applicable' cannot both hold "
                "on an untested class; a robustness claim needs a measurement"
            )
        if problems:
            raise ValueError("; ".join(problems))
        return self

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Validation":
        try:
            return cls.model_validate(data)
        except ValidationError as exc:
            raise _translate(exc, "validation") from exc


class EndpointRecord(BaseModel):
    """The two layers as one serialized unit, plus what the endpoint is."""

    model_config = ConfigDict(extra="forbid")

    schema_version: Literal[SCHEMA_VERSION] = SCHEMA_VERSION
    physical_state_id: str | None = None
    n_atoms: int = Field(gt=0)
    recipe: Recipe
    validation: Validation
    magnetic_note: str | None = None
    vpa_is_proxy: bool | None = None

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "EndpointRecord":
        try:
            return cls.model_validate(data)
        except ValidationError as exc:
            raise _translate(exc, "record") from exc


def load_pair(data: Mapping[str, Any]) -> tuple[Recipe, Validation]:
    """Build both layers from a two-key mapping."""

    try:
        recipe_data = data["recipe"]
        validation_data = data["validation"]
    except KeyError as exc:
        raise MissingRecipeField("record", [str(exc.args[0])]) from exc
    return Recipe.from_dict(recipe_data), Validation.from_dict(validation_data)


__all__ = [
    "CLASS_EXPECTATIONS",
    "EndpointRecord",
    "FROZEN_DISPLACEMENT_ANGSTROM",
    "FROZEN_POTCAR_FUNCTIONAL",
    "FROZEN_SMEARING",
    "GAMMA_SHIFT",
    "Kpoints",
    "MAGNETIC_POLICIES",
    "MissingRecipeField",
    "POTCAR_KNOWN_VARIANTS",
    "Potcar",
    "QHA_APPLICABILITIES",
    "QhaApplicability",
    "RECIPE_ROBUSTNESSES",
    "REPRESENTATIONS",
    "Recipe",
    "RecipeError",
    "RecipeInconsistent",
    "RecipeRobustness",
    "SCHEMA_VERSION",
    "STAGE0_COVERAGES",
    "STRUCTURE_CLASSES",
    "Smearing",
    "Stage0Coverage",
    "StructureClass",
    "Validation",
    "axis_contains_gamma",
    "gamma_centred",
    "load_pair",
]
