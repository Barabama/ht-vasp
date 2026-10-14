"""Turn a `Recipe` into the VASP inputs for one stage of a DFPT-QHA endpoint.

The stages
----------
``relax``   ISIF=3 ionic + volume relaxation, to locate the endpoint's V0.
``static``  a single point at a fixed cell -- one E(V) sample.
``phonon``  DFPT (``IBRION=8``), which is what the production line uses.
``shape``   constant-volume shape relaxation (``ISIF=4``) for a cell that has
            a genuine shape degree of freedom at fixed volume (hcp's c/a).

Why ``shape`` exists as its own stage
-------------------------------------
The discriminator must be the CRYSTAL SYSTEM, not the lattice angles.  An
earlier version tested ``angles[0] != 90``, which is wrong in both directions:
the fcc primitive cell has 60 deg angles and the bcc primitive cell 109.47 deg,
and both are cubic (shape fixed by symmetry), while hcp has 90/90/120 and its
c/a is genuinely free.  That test handed ``ISIF=4`` to the two cells that did
not need it and ``ISIF=2`` to the one that did.

Submitting is gated
-------------------
``write_stage_inputs`` writes the inputs and then runs the input gate on the
directory it just wrote, returning the verdict.  There is no way to obtain a
written input set without a verdict attached, which is the point: the project
has three incidents that all came from inputs nobody checked.

The magnetic policy is not guessed
----------------------------------
``Recipe.magnetic_policy == "unresolved"`` makes ``magmom_for`` raise.  A
recipe that records its magnetism as undecided must not silently acquire the
table default and produce 84 endpoints of unexamined spin states.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from htvasp.dfpt_qha.gate import GateResult, validate_dir
from htvasp.dfpt_qha.recipe import FROZEN_SMEARING, Recipe

Stage = Literal["relax", "static", "phonon", "shape"]

STAGES: tuple[Stage, ...] = ("relax", "static", "phonon", "shape")

# --------------------------------------------------------------------------
# frozen INCAR
#
# The two families have DIFFERENT electronic settings and must not be merged.
# The values below are transcribed from INCARs that actually ran, not chosen
# here; the file each came from is named with it.
# --------------------------------------------------------------------------

#: settings shared by every stage
SHARED_INCAR: dict[str, Any] = {
    "ENCUT": 520,
    "PREC": "Normal",
    "ALGO": "Normal",
    "ISYM": 2,
    "LORBIT": 11,          # output-only: it is what makes per-site moments readable
    "LREAL": False,
    "LWAVE": False,
    "LCHARG": False,
    "ISTART": 0,
    "ICHARG": 2,
}

#: electronic convergence and parallel layout, per family
STATIC_ELECTRONIC: dict[str, Any] = {
    "EDIFF": 1e-6, "NELM": 120, "NELMIN": 6, "NCORE": 8,
}
#: transcribed from /nfs_ssd/tmp/legacy-qha-pool/e1-dfpt-al-ta/INCAR and
#: .../d0-dfpt-3x/INCAR, which agree.  DFPT is looser (1e-5), takes one ionic
#: step, and parallelises over k (KPAR=2) rather than bands -- NCORE and KPAR
#: are mutually exclusive, so NCORE is absent rather than overridden.
DFPT_ELECTRONIC: dict[str, Any] = {
    "EDIFF": 1e-5, "NELM": 120, "NELMIN": 10, "KPAR": 2,
}

#: ionic/electronic mode per stage
STAGE_MODE: dict[Stage, dict[str, Any]] = {
    "relax": {"IBRION": 2, "ISIF": 3, "NSW": 80, "POTIM": 0.2, "EDIFFG": -0.02},
    "static": {"IBRION": -1, "ISIF": 2, "NSW": 0},
    "phonon": {"IBRION": 8, "ISIF": 2, "NSW": 1},
    "shape": {"IBRION": 2, "ISIF": 4, "NSW": 80, "POTIM": 0.2, "EDIFFG": -0.02},
}

#: which electronic family each stage draws its convergence settings from
STAGE_ELECTRONIC: dict[Stage, dict[str, Any]] = {
    "relax": STATIC_ELECTRONIC,
    "static": STATIC_ELECTRONIC,
    "phonon": DFPT_ELECTRONIC,
    "shape": STATIC_ELECTRONIC,
}

#: kept for callers that want the whole frozen incar without a stage
FROZEN_INCAR: dict[str, Any] = {
    **SHARED_INCAR, **STATIC_ELECTRONIC, **STAGE_MODE["static"],
}

# Starting moments per element, per magnetic policy.  `confirmed_branch` is
# only ever set on an endpoint where a calculation actually showed a magnetic
# ground state -- see the recipe's `magnetic_policy`.
BASE_MOMENTS: dict[str, float] = {
    "Al": 1.0, "Cr": 3.0, "Hf": 1.0, "Nb": 3.0, "Ta": 3.0, "Ti": 3.0, "Zr": 1.0,
}

# Reference cell for the k-density definition: the legacy-relaxed B2 AlTa.
DENSITY_REFERENCE_A = 3.2246015525754168
DENSITY_DEFAULT = 16
GAMMA_SHIFT = 0.5


class InputSetError(ValueError):
    """A recipe cannot be turned into inputs without a decision being made."""


def round_even(value: float) -> int:
    """Nearest even integer, minimum 2.

    Halve, round, double -- rather than rounding and then nudging the odd
    result.  The two differ exactly on odd integers: 17.0 gives 16 here (17/2
    rounds to 8 under round-half-to-even) and 18 under the nudging form.  The
    halving form is the one the Stage 0 batches actually used, so keeping it
    is what makes a recomputed mesh equal a recorded one.

    Even is not cosmetic: only an even mesh takes shift 0.5 and still contains
    Gamma, which is what makes it fold-matchable to a supercell mesh.
    """

    return max(int(round(value / 2.0)) * 2, 2)


def mesh_for(volume_angstrom3: float, density: int = DENSITY_DEFAULT) -> int:
    """Even ``n`` whose k-sampling density matches ``density`` on the reference.

    The invariant is the number of k-points per unit reciprocal volume::

        n^3 * V_cell = const      (the BZ volume is (2*pi)^3 / V_cell)

    so ``n = (a_ref * density) / V_cell^(1/3)``.  For a cubic cell
    ``V^(1/3) == a`` and this reduces to ``n = a_ref*density/a``, but for a
    NON-cubic cell -- e.g. the primitive fcc cell -- using ``a`` instead of
    ``V^(1/3)`` under-samples it by 4^(1/3) = 1.59x.  That is why this takes
    the volume.
    """

    return round_even((DENSITY_REFERENCE_A * density) / volume_angstrom3 ** (1.0 / 3.0))


def potcar_symbols(structure: Any, recipe: Recipe) -> list[str]:
    """POTCAR symbols in POSCAR order, with the pool's recorded exceptions.

    This is not cosmetic.  ``Potcar(["Nb"], functional="PBE")`` raises
    FileNotFoundError: potpaw_PBE has no plain Nb, only Nb_pv and Nb_sv.  A
    builder that ignores ``recipe.potcar.variants`` cannot produce an
    Nb-bearing endpoint at all.
    """

    variants = recipe.potcar.variants
    default = variants.get("default", "standard")
    symbols: list[str] = []
    for species in dict.fromkeys(str(site.specie) for site in structure):
        # A variant entry holds the FULL POTCAR symbol ("Nb_pv"), not a suffix,
        # so it is used as-is.  "standard" means the plain element symbol.
        variant = variants.get(species, default)
        symbols.append(species if variant == "standard" else variant)
    return symbols


#: starting moments, by policy.  The policy is a CLAIM about the ground state,
#: not an input switch: the Stage 0 batches ran every endpoint spin-polarised
#: (ISPIN=2) with a starting moment, so that a "nonmagnetic" expectation is
#: something the calculation can falsify.  ISPIN=1 appeared exactly once, for
#: the dedicated magnetic screening case.
POLICY_ISPIN: dict[str, int] = {
    "nonmagnetic": 2,
    "confirmed_branch": 2,
    # "unresolved" is absent on purpose -- it raises rather than defaulting.
}


def magmom_for(structure: Any, recipe: Recipe) -> dict[str, float]:
    """Starting moments per element, for the recipe's magnetic policy.

    Raises for ``unresolved`` rather than falling back to the table: an
    endpoint whose spin state nobody has decided must not quietly acquire one.

    ``nonmagnetic`` and ``confirmed_branch`` produce the SAME starting moments.
    That is deliberate -- both run spin-polarised from the same table, and the
    policy records what is believed about the answer, not what is fed in.  The
    difference between them lives in ``Validation``, where ``measured_on`` has
    to be non-empty for a claim of a confirmed branch.
    """

    if recipe.magnetic_policy == "unresolved":
        raise InputSetError(
            "magnetic_policy is 'unresolved': the spin state of this endpoint "
            "has not been decided, so generating inputs would pick one by "
            "default.  Run the spin screening first, or set the policy "
            "explicitly."
        )
    return {
        str(site.specie): BASE_MOMENTS.get(str(site.specie), 1.0)
        for site in structure
    }


def incar_for(stage: Stage, recipe: Recipe, *, is_cubic: bool = True) -> dict[str, Any]:
    """The INCAR for one stage, using the recipe's frozen values."""

    if stage not in STAGE_MODE:
        raise InputSetError(f"unknown stage {stage!r}; expected one of {STAGES}")
    if stage == "shape" and is_cubic:
        # A cubic cell has no shape degree of freedom at fixed volume, so
        # ISIF=4 would be a no-op that costs an ionic loop.
        raise InputSetError(
            "stage='shape' is only meaningful for a non-cubic cell; a cubic "
            "cell's shape is fixed by symmetry"
        )
    incar = {**SHARED_INCAR, **STAGE_ELECTRONIC[stage], **STAGE_MODE[stage]}
    incar["ISMEAR"] = FROZEN_SMEARING["ISMEAR"]
    incar["SIGMA"] = FROZEN_SMEARING["SIGMA"]
    # ISPIN is written explicitly rather than left to VASP's default of 1.
    # Left out, MAGMOM is ignored entirely and an endpoint whose recipe claims
    # a magnetic branch would run non-magnetic while still looking correct.
    try:
        incar["ISPIN"] = POLICY_ISPIN[recipe.magnetic_policy]
    except KeyError as exc:
        raise InputSetError(
            f"magnetic_policy {recipe.magnetic_policy!r} has no ISPIN mapping; "
            f"'unresolved' must be screened before inputs can be generated"
        ) from exc
    return incar


def kpoints_for(recipe: Recipe) -> tuple[int, float]:
    """The cell mesh the recipe records, as ``(n, shift)``.

    Read from the recipe rather than recomputed, so the inputs and the record
    cannot disagree -- ``Recipe.check_against`` is what proves the recorded
    mesh matches the cell.
    """

    grid, shift = recipe.kpoints.grid, recipe.kpoints.shift
    if len(set(grid)) != 1 or len(set(shift)) != 1:
        raise InputSetError(
            f"the pool writes isotropic meshes; got grid={list(grid)} "
            f"shift={list(shift)}"
        )
    return grid[0], shift[0]


def supercell_kpoints_for(recipe: Recipe) -> tuple[int, float]:
    """The fold-matched supercell mesh, which is what the phonon cell uses."""

    grid, shift = recipe.kpoints.supercell_grid, recipe.kpoints.supercell_shift
    if len(set(grid)) != 1 or len(set(shift)) != 1:
        raise InputSetError(
            f"the pool writes isotropic meshes; got supercell_grid={list(grid)} "
            f"supercell_shift={list(shift)}"
        )
    return grid[0], shift[0]


def write_stage_inputs(
    directory: Path,
    structure: Any,
    recipe: Recipe,
    stage: Stage,
    *,
    is_cubic: bool = True,
    magmom_override: dict[str, float] | None = None,
) -> GateResult:
    """Write POSCAR/INCAR/KPOINTS/POTCAR and return the gate's verdict.

    The verdict is returned rather than merely logged, and there is no code
    path that writes the inputs without producing one -- a caller that ignores
    the result has to do so explicitly.
    """

    from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    structure.to(filename=str(directory / "POSCAR"), fmt="poscar")

    incar = incar_for(stage, recipe, is_cubic=is_cubic)
    # MAGMOM is written only when ISPIN=2.  With ISPIN=1 VASP ignores it, and a
    # file that carries a moment it does not use invites the reader to believe
    # the run was spin-polarised.
    if incar["ISPIN"] == 2:
        magmom = magmom_override or magmom_for(structure, recipe)
        incar["MAGMOM"] = [magmom[str(site.specie)] for site in structure]
    Incar(incar).write_file(str(directory / "INCAR"))

    mesh, shift = kpoints_for(recipe)
    Kpoints.from_str(
        "Automatic kpoint scheme\n0\nMonkhorst-Pack\n"
        f"{mesh} {mesh} {mesh}\n{shift} {shift} {shift}\n"
    ).write_file(str(directory / "KPOINTS"))

    Potcar(potcar_symbols(structure, recipe),
           functional=recipe.potcar.functional).write_file(str(directory / "POTCAR"))

    return validate_dir(directory)


__all__ = [
    "BASE_MOMENTS",
    "DENSITY_DEFAULT",
    "DENSITY_REFERENCE_A",
    "DFPT_ELECTRONIC",
    "FROZEN_INCAR",
    "InputSetError",
    "POLICY_ISPIN",
    "SHARED_INCAR",
    "STAGES",
    "STAGE_ELECTRONIC",
    "STAGE_MODE",
    "STATIC_ELECTRONIC",
    "Stage",
    "incar_for",
    "kpoints_for",
    "magmom_for",
    "mesh_for",
    "potcar_symbols",
    "round_even",
    "supercell_kpoints_for",
    "write_stage_inputs",
]
