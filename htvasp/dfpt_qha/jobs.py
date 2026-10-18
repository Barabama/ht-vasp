"""DFPT jobflow boundaries.

The collector consumes the output of an atomate2 ``BaseVaspMaker`` job.  It does
not manufacture a finite-displacement dataset: the dynamical matrix is read from
that job's own ``DYNMAT``.

Flow shape, per endpoint::

    CommonEosMaker                 E(V) and the relaxed cell at each volume
      -> spawn_dfpt_volume_jobs    one DFPT run per volume, on the recipe supercell
           -> collect_dfpt_task_doc    VASP DFPT output -> force constants
           -> validate_dfpt_task_doc   fail-closed: cell, shape, provenance
      -> collect_validated_volumes  gather the survivors
      -> classify_dfpt_stability    one verdict per volume, endpoint gate
      -> dfpt_thermal_properties    F(V,T) from the dynamical matrices
      -> analyze_dfpt_qha           G(T), V(T) and the EOS fit

The validation job is a pass-through, so a stopped volume stops the thermal and
QHA jobs behind it.  Stability is a single job after collection rather than one
per volume: the decision belongs to the endpoint (a soft volume fails the
endpoint instead of being dropped from the fit), and the VASP time is already
spent by then, so a per-volume fail-fast would buy nothing.

What VASP actually writes
-------------------------
``IBRION=8`` writes **DYNMAT**, not ``FORCE_CONSTANTS``.  The verified reference
directory has both only because ``FORCE_CONSTANTS`` was produced 90 minutes
after the run by an offline conversion; assuming VASP writes it means the
collector finds nothing after the allocation has already been spent.

``DYNMAT`` holds, for each perturbation ``(atom, direction)``, the force
response on every atom, in eV/Angstrom^2 because the perturbation is a unit
displacement.  That is the same quantity phonopy calls force constants, so the
matrix is recovered by replaying the perturbations through phonopy with
``distance=1.0`` and attaching the responses.  No reordering is needed here
because the DFPT job runs on phonopy's OWN supercell (see
``htvasp.dfpt_qha.phonons``), which is what makes VASP's atom order and
phonopy's index order the same list.
"""

from __future__ import annotations

import gzip
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from phonopy.harmonic.force_constants import (
    set_permutation_symmetry,
    set_translational_invariance,
)
from phonopy.interface.vasp import Vasprun as VaspVasprun
from atomate2.vasp.jobs.base import BaseVaspMaker
from atomate2.utils.path import strip_hostname
from jobflow import Response, job
from monty.io import zopen
from monty.json import jsanitize
from pymatgen.core import Structure

from htvasp.dfpt_qha.input_audit import audit_vasp_inputs
from htvasp.dfpt_qha.input_set import DfptInputSetGenerator, DfptRelaxSetGenerator
from htvasp.dfpt_qha.phonons import (
    assert_same_ordering,
    ordering_fingerprint,
    phonopy_for,
    supercell_structure,
)
from htvasp.dfpt_qha.task_doc import DfptForceConstantsProvenance, DfptTaskDoc
from htvasp.model.dfpt_qha_recipe import Recipe

#: 1 eV/Angstrom^2 expressed in GPa (for the bulk modulus, if needed).
DYNMAT_UNIT = "eV/Angstrom^2"

#: Above this relative acoustic-sum-rule residual the raw matrix is not a
#: force-constant matrix and must not be "repaired" into one.  VASP's own DFPT
#: Hessian is not ASR-consistent -- the archived Co-Ti runs sit at 8e-3 to
#: 2e-2, and one broken volume at 0.36 (Gamma acoustic -14 THz, optical
#: 26.8 THz where the healthy runs are 5-7 THz).  Enforcing translational
#: invariance is standard practice and changes the optical frequencies of the
#: healthy runs by 0.1-0.2%; on the broken one it would hide a 22% error.
ASR_CORRECTION_MAX_REL = 0.1
#: Gate applied to the matrix AFTER the correction.
ASR_TOL_REL = 1e-2

#: VASP's own Hessian.  This is the authoritative source: it is complete,
#: exactly symmetric, and independent of phonopy's symmetry-reduced finite
#: displacement list.  ``DYNMAT`` is the fallback.
VASPRUN_FC_SOURCE = "vasprun.xml dynmat hessian"
VASPRUN_FC_PARSER = "phonopy VASP vasprun parser"
DYNMAT_FC_SOURCE = "DYNMAT"
DYNMAT_FC_PARSER = "phonopy replay at unit displacement"


def _find_file(directory: Path, name: str) -> Path | None:
    """``name`` or ``name.gz`` inside ``directory``, whichever exists.

    BaseVaspMaker gzips a fixed file list at the end of every job, so a
    collector that only looks for the plain name breaks the moment that list
    changes.  Both spellings are accepted.
    """

    for candidate in (directory / name, directory / f"{name}.gz"):
        if candidate.is_file():
            return candidate
    return None


def parse_dynmat(path: Path | str) -> tuple[list[dict[str, Any]], int]:
    """Read VASP's ``DYNMAT``: the DFPT force response per perturbation.

    Layout::

        line 0   n_species  n_atoms  n_perturbations
        line 1   masses
        then, per perturbation:
            atom direction  displacement(3)
            n_atoms lines of (force_x force_y force_z)

    Returns ``(perturbations, n_atoms)`` with each perturbation as
    ``{"atom": 1-based index, "direction": 1-based, "displacement": [3],
    "forces": (n_atoms, 3) array}``.
    """

    with zopen(path, mode="rt") as handle:
        lines = handle.read().splitlines()
    if len(lines) < 3:
        raise ValueError(f"DYNMAT is too short to parse: {path}")

    header = lines[0].split()
    n_atoms = int(header[1])
    n_perturbations = int(header[2])
    if n_atoms < 1 or n_perturbations < 1:
        raise ValueError(f"DYNMAT header is not usable: {lines[0]!r}")

    perturbations: list[dict[str, Any]] = []
    index = 2
    for _ in range(n_perturbations):
        if index >= len(lines):
            raise ValueError(f"DYNMAT ended after {len(perturbations)} perturbations")
        head = lines[index].split()
        atom, direction = int(head[0]), int(head[1])
        displacement = [float(x) for x in head[2:5]] if len(head) >= 5 else [0.0, 0.0, 0.0]
        index += 1
        block = lines[index:index + n_atoms]
        if len(block) != n_atoms:
            raise ValueError(
                f"DYNMAT perturbation {len(perturbations) + 1} has {len(block)} "
                f"force lines, expected {n_atoms}"
            )
        forces = np.array([[float(x) for x in row.split()[:3]] for row in block], dtype=float)
        index += n_atoms
        perturbations.append(
            {
                "atom": atom,
                "direction": direction,
                "displacement": displacement,
                "forces": forces,
            }
        )
    return perturbations, n_atoms


def force_constants_from_vasprun(
    path: Path | str, supercell: Structure
) -> np.ndarray:
    """DFPT force constants from the ``<dynmat><hessian>`` block of vasprun.xml.

    VASP's Hessian is the authoritative DFPT output.  Unlike the ``DYNMAT``
    text file it is complete, and unlike a phonopy finite-displacement fit it
    does not depend on VASP's perturbation list matching phonopy's
    symmetry-reduced one -- which it does not for FCC/L1_2 and HCP.

    phonopy's parser handles the two conventions VASP writes (VASP 5 in
    eV/Angstrom^2 with the sqrt(m_i m_j) mass factor, VASP 6 in THz^2), and
    applies the same mass normalisation, so the result is directly comparable
    with the ``DYNMAT`` replay: on the archived BCC-Al-Ta-2 run the two agree
    to 2.4e-9 eV/Angstrom^2 after symmetrisation.
    """

    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as handle:
        force_constants, elements = VaspVasprun(handle).read_force_constants()

    force_constants = np.asarray(force_constants, dtype=float)
    expected = (len(supercell), len(supercell), 3, 3)
    if force_constants.shape != expected:
        raise ValueError(
            f"{path.name} carries force constants of shape "
            f"{force_constants.shape}, expected {expected}"
        )

    symbols = [str(site.specie) for site in supercell]
    atomtypes: list[str] = []
    for symbol in symbols:
        if not atomtypes or atomtypes[-1] != symbol:
            atomtypes.append(symbol)
    if list(elements) != atomtypes:
        raise ValueError(
            f"{path.name} declares atom types {list(elements)} but the cell is "
            f"{atomtypes}; the force-constant matrix is indexed by the other "
            f"order"
        )
    return force_constants


def force_constants_from_dynmat(
    cell: Structure,
    supercell: Structure,
    perturbations: list[dict[str, Any]],
    *,
    supercell_matrix: Any,
) -> np.ndarray:
    """Force constants from a DFPT ``DYNMAT``, in the supercell's atom order.

    phonopy is asked to generate its own perturbation list at unit displacement
    and then handed the DFPT responses in that same order.  The VASP order and
    phonopy's order are the same because the job ran on
    ``phonopy.supercell``; that is asserted rather than assumed.
    """

    phonon = phonopy_for(cell, supercell_matrix)
    phonopy_supercell = supercell_structure(cell, supercell_matrix)
    assert_same_ordering(supercell, phonopy_supercell)

    n_supercell = len(phonopy_supercell)
    phonon.generate_displacements(distance=1.0)
    expected = len(phonon.supercells_with_displacements)
    if expected != len(perturbations):
        raise ValueError(
            f"DYNMAT has {len(perturbations)} perturbations but phonopy expects "
            f"{expected} at unit displacement; the two are not the same "
            f"calculation and the responses cannot be attached"
        )
    forces = []
    # phonopy reports displacements as [atom_index, dx, dy, dz] with a
    # zero-based atom index; VASP's DYNMAT numbers atoms and directions from 1.
    for index, (record, phonopy_disp) in enumerate(zip(perturbations, phonon.displacements)):
        atom, vector = int(phonopy_disp[0]), [float(x) for x in phonopy_disp[1:4]]
        if record["atom"] != atom + 1:
            raise ValueError(
                f"perturbation {index + 1} displaces atom {record['atom']} but "
                f"phonopy expects atom {atom + 1}"
            )
        if not np.allclose(record["displacement"], vector, atol=1e-6):
            raise ValueError(
                f"perturbation {index + 1} displaces along "
                f"{record['displacement']} but phonopy expects {vector}; the two "
                f"are not the same calculation"
            )
        array = np.asarray(record["forces"], dtype=float)
        if array.shape != (n_supercell, 3):
            raise ValueError(
                f"perturbation {index + 1} carries forces of shape {array.shape}, "
                f"expected {(n_supercell, 3)}"
            )
        forces.append(array)
    phonon.forces = forces
    phonon.produce_force_constants()
    return np.array(phonon.force_constants, dtype=float)


@dataclass
class DfptRelaxMaker(BaseVaspMaker):
    """Ionic relaxation for the EOS stage: ISIF=3 (free) or ISIF=2 (fixed volume).

    Deliberately NOT the DFPT maker.  Wiring ``IBRION=8`` into the EOS stage
    produces an E(V) curve of unrelaxed strained geometries, because VASP
    ignores ISIF/NSW for a dynamical-matrix run.
    """

    name: str = "DFPT-QHA relax"
    recipe: Recipe | None = None
    isif: int = 2
    use_supercell_kpoints: bool = False
    input_set_generator: Any = None
    run_vasp_kwargs: dict[str, Any] = field(default_factory=dict)
    stop_children_kwargs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.input_set_generator is None:
            if self.recipe is None:
                raise ValueError("DfptRelaxMaker requires an explicit Recipe")
            self.input_set_generator = DfptRelaxSetGenerator(
                recipe=self.recipe,
                isif=self.isif,
                use_supercell_kpoints=self.use_supercell_kpoints,
            )
        # BaseVaspMaker is a plain dataclass in this atomate2 release and
        # defines no post-init hook.


@dataclass
class DfptVaspMaker(BaseVaspMaker):
    """A standard atomate2 BaseVaspMaker configured for IBRION=8.

    A dataclass (not a hand-written ``__init__``) because jobflow serializes
    makers through monty, which requires every constructor argument to exist
    as a dataclass field on the instance.
    """

    name: str = "DFPT job"
    recipe: Recipe | None = None
    use_supercell_kpoints: bool = True
    input_set_generator: Any = None
    run_vasp_kwargs: dict[str, Any] = field(default_factory=dict)
    stop_children_kwargs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.input_set_generator is None:
            if self.recipe is None:
                raise ValueError("DfptVaspMaker requires an explicit Recipe")
            self.input_set_generator = DfptInputSetGenerator(
                recipe=self.recipe, use_supercell_kpoints=self.use_supercell_kpoints
            )


def _task_payload(task_doc: Any) -> dict[str, Any]:
    """Reduce whatever a jobflow boundary handed us to a JSON-safe mapping.

    Three shapes arrive here, and which one depends on the runner rather than on
    the code:

    * a plain ``dict`` -- a job that wrote one;
    * a ``TaskDoc`` pydantic model -- ``BaseVaspMaker``'s output schema;
    * a dict monty's decoder rebuilt from JSON, whose nested ``input.kpoints``,
      ``structure`` and ``orig_inputs.incar`` are live pymatgen objects again.

    The third is the one that bites: it survives until something serializes the
    document a second time, and then raises ``Unable to serialize unknown type:
    Kpoints`` -- which happened three jobs downstream, inside the validation
    gate, where the cause is hardest to read.  Sanitizing at this boundary makes
    "this document is JSON" true for everything after it.

    ``jsanitize`` is asked for ``strict=True`` because ``strict=False`` falls
    back to ``str(obj)`` for anything it does not recognise, which silently
    turns a pydantic model into its repr.  ``strict=True`` raises instead, and
    the repr of a whole TaskDoc is exactly the kind of quiet corruption this
    project keeps finding.
    """

    payload = jsanitize(
        task_doc, strict=True, allow_bson=True, recursive_msonable=True
    )
    if not isinstance(payload, dict):
        raise ValueError(
            f"DFPT TaskDoc did not reduce to a mapping: {type(task_doc).__name__}"
        )
    # jsanitize keeps datetime objects, which the job store accepts but a plain
    # JSON dump does not; the document is read as data downstream, so they are
    # rendered as ISO strings.
    return json.loads(json.dumps(payload, default=str))


#: A DFPT run recomputes the ground state on its own (larger) cell.  Its energy
#: per atom must therefore agree with the static energy of the volume it was
#: built from.  Healthy agreement on this pool is a few meV/atom; a diverged
#: SCF shows up as hundreds or thousands.  See the BCC-Co-Ti-2 case where a
#: ground state that ran away to +4.4e3 eV still carried no "unconverged"
#: marker anywhere in the OUTCAR.
STATIC_ENERGY_TOL_MEV_PER_ATOM = 100.0


@job
def collect_dfpt_task_doc(
    task_doc: Any,
    *,
    recipe: dict[str, Any],
    cell: dict[str, Any],
    supercell_matrix: Any = None,
    vasprun_path: str | Path | None = None,
    expected_static_energy_per_atom: float | None = None,
) -> dict[str, Any]:
    """Convert one completed atomate2 TaskDoc into a DFPT TaskDoc.

    The collector is a pure jobflow analysis boundary: it reads the job's own
    DFPT output.  It never launches VASP and never scans a calculation tree it
    was not handed.

    Force constants come from VASP, in this order:

    1. ``<dynmat><hessian>`` in ``vasprun.xml`` -- complete and independent of
       phonopy's symmetry-reduced displacement list;
    2. the ``DYNMAT`` text file replayed through phonopy -- the fallback for a
       VASP that did not write the Hessian.

    ``cell`` is the primitive cell the supercell was built from.  It is passed
    explicitly rather than recovered from the supercell so that the ordering
    assertion compares two independently-known structures.
    """

    task_payload = _task_payload(task_doc)
    directory = task_payload.get("dir_name")

    # The TaskDoc state is the job-level verdict: emmet derives it from
    # ``vasprun.converged``, which is electronic AND ionic convergence of the
    # final step.  Custodian's ``VaspErrorHandler`` has 45 text signatures and
    # none of them matches "aborting loop EDIFF was not reached", so a
    # non-converged SCF reaches this point looking like a finished run.  The
    # collector is the last place where the run's own verdict is still attached,
    # so it is checked here rather than three jobs downstream.
    state = task_payload.get("state")
    if state != "successful":
        raise ValueError(
            f"DFPT TaskDoc is not successful: state={state!r}.  The calculation "
            f"did not converge (or was otherwise marked failed); its dynamical "
            f"matrix must not be used."
        )

    # TaskDoc.dir_name is an emmet URI ("host:/abs/path"); atomate2 strips the
    # host for the same value in copy_vasp_outputs, so the collector must too.
    job_dir = Path(strip_hostname(str(vasprun_path or directory)))

    structure = task_payload.get("structure")
    if structure is None:
        raise ValueError("DFPT TaskDoc carries no final structure to fingerprint")
    supercell = Structure.from_dict(structure)

    parsed = Recipe.from_dict(recipe)
    factor = int(parsed.kpoints.supercell_factor)
    if supercell_matrix is None:
        supercell_matrix = [[factor, 0, 0], [0, factor, 0], [0, 0, factor]]

    vasprun_file = _find_file(job_dir, "vasprun.xml")
    dynmat_file = _find_file(job_dir, "DYNMAT")
    if vasprun_file is None and dynmat_file is None:
        raise ValueError(
            f"DFPT job directory has neither vasprun.xml nor DYNMAT: {job_dir}.  "
            f"VASP writes at least one of them for IBRION=8; their absence means "
            f"the run did not reach the dynamical-matrix step."
        )

    # DYNMAT is always parsed when present: even when the Hessian is used, the
    # perturbation list it records is what lets a later reader reconstruct the
    # run, and it is the only independent cross-check on the Hessian.
    dynmat_provenance: dict[str, Any] | None = None
    perturbations: list[dict[str, Any]] = []
    if dynmat_file is not None:
        perturbations, n_atoms_dynmat = parse_dynmat(dynmat_file)
        if n_atoms_dynmat != len(supercell):
            raise ValueError(
                f"DYNMAT describes {n_atoms_dynmat} atoms but the TaskDoc "
                f"structure has {len(supercell)}"
            )
        dynmat_provenance = {
            "path": str(dynmat_file),
            "parser": "VASP DYNMAT text parser",
            "unit": DYNMAT_UNIT,
            "n_perturbations": len(perturbations),
            "perturbation_atoms": [int(p["atom"]) for p in perturbations],
            "perturbation_directions": [
                [float(x) for x in p["displacement"]] for p in perturbations
            ],
            "role": "fallback and audit",
        }

    force_constants: np.ndarray | None = None
    source = parser = None
    force_constants_path: Path | None = None
    vasprun_error: str | None = None
    if vasprun_file is not None:
        try:
            force_constants = force_constants_from_vasprun(vasprun_file, supercell)
            source, parser = VASPRUN_FC_SOURCE, VASPRUN_FC_PARSER
            force_constants_path = vasprun_file
        except Exception as exc:  # noqa: BLE001 - recorded, then the fallback runs
            vasprun_error = f"{type(exc).__name__}: {exc}"
    if force_constants is None:
        if dynmat_file is None:
            raise ValueError(
                f"no usable DFPT force constants in {job_dir}: vasprun.xml "
                f"carries no Hessian ({vasprun_error}) and there is no DYNMAT"
            )
        force_constants = force_constants_from_dynmat(
            Structure.from_dict(cell), supercell, perturbations,
            supercell_matrix=supercell_matrix,
        )
        source, parser = DYNMAT_FC_SOURCE, DYNMAT_FC_PARSER
        force_constants_path = dynmat_file
    shape = tuple(int(x) for x in force_constants.shape)
    if len(shape) != 4 or shape[0] != shape[1] or shape[2:] != (3, 3):
        raise ValueError(f"unexpected DFPT force_constants shape: {shape}")

    # Acoustic sum rule: sum_j Phi_ij must vanish for every i.  VASP does NOT
    # enforce it on the raw DFPT Hessian, and the residual is not small: the
    # archived VASP 5 Co-Ti volumes sit at 8e-3 to 2e-2 relative.  Rejecting
    # them outright would fail every endpoint; accepting them uncorrected puts
    # a spurious ~0.1 THz on the Gamma acoustic branches and contaminates the
    # thermal sums.  So the standard correction is applied here, the residual
    # before and after is recorded, and a matrix the correction cannot fix is
    # passed through untouched for the classification job to reject with a
    # stored verdict (rather than raising here where nothing is recorded).
    scale = float(np.abs(force_constants).max()) or 1.0
    residual_before = float(np.abs(force_constants.sum(axis=1)).max())
    asr_rel_before = residual_before / scale
    asr_correction_applied = False
    residual = residual_before
    if asr_rel_before <= ASR_CORRECTION_MAX_REL:
        corrected = np.array(force_constants, dtype=float)
        set_translational_invariance(corrected)
        set_permutation_symmetry(corrected)
        residual_after = float(np.abs(corrected.sum(axis=1)).max())
        if residual_after / scale <= ASR_TOL_REL:
            force_constants = corrected
            residual = residual_after
            asr_correction_applied = True

    # Physical convergence gate.  The task state above is the job-level verdict,
    # and it is derived from "did the last electronic step stop before NELM" --
    # which a run that diverged and then happened to stop early would also
    # satisfy.  This compares the DFPT run's own ground-state energy with the
    # static energy of the volume it was launched from: the two are the same
    # crystal at the same volume, so they must agree to a few meV/atom.
    energy_delta_mev_per_atom: float | None = None
    if expected_static_energy_per_atom is not None:
        output = task_payload.get("output") or {}
        energy = output.get("energy")
        n_atoms = len(supercell)
        if energy is None:
            raise ValueError("DFPT TaskDoc carries no total energy to compare")
        actual_per_atom = float(energy) / n_atoms
        energy_delta_mev_per_atom = 1e3 * (
            actual_per_atom - float(expected_static_energy_per_atom)
        )
        if abs(energy_delta_mev_per_atom) > STATIC_ENERGY_TOL_MEV_PER_ATOM:
            raise ValueError(
                f"DFPT ground state disagrees with the static energy of its "
                f"volume by {energy_delta_mev_per_atom:.1f} meV/atom "
                f"(DFPT {actual_per_atom:.6f} vs static "
                f"{float(expected_static_energy_per_atom):.6f} eV/atom, tolerance "
                f"{STATIC_ENERGY_TOL_MEV_PER_ATOM:g}).  The SCF did not converge to "
                f"the same state the EOS used; the dynamical matrix is not usable."
            )
    # custodian is allowed to correct inputs; what it changed is recorded here
    # rather than used to reject the run.  A POSCAR correction is the one that
    # matters most, because the matrix is indexed by that atom order.
    input_audit = audit_vasp_inputs(job_dir)

    doc = DfptTaskDoc(
        task_doc=task_payload,
        recipe=recipe,
        force_constants=force_constants.tolist(),
        force_constants_provenance=DfptForceConstantsProvenance(
            source=source,
            parser=parser,
            n_atoms=shape[0],
            shape=shape,
            unit=DYNMAT_UNIT,
            ordering_fingerprint=ordering_fingerprint(supercell),
            path=str(force_constants_path),
            max_acoustic_sum_residual=residual,
        ),
        provenance={
            "task_doc_dir": str(job_dir),
            "task_state": task_payload.get("state"),
            "n_perturbations": len(perturbations),
            "vasprun_hessian_error": vasprun_error,
            "asr_rel_before": asr_rel_before,
            "asr_rel_after": residual / scale,
            "asr_correction_applied": asr_correction_applied,
            "dynmat": dynmat_provenance,
            "input_audit": input_audit,
            "energy_vs_static_meV_per_atom": energy_delta_mev_per_atom,
            "energy_vs_static_tolerance_meV_per_atom": STATIC_ENERGY_TOL_MEV_PER_ATOM,
        },
    )
    return doc.model_dump(mode="json")


@job
def validate_dfpt_task_doc(
    dfpt_doc: dict[str, Any],
    *,
    recipe: dict[str, Any] | None = None,
) -> Response:
    """Stop descendants unless the DFPT result matches the recipe it claims.

    The shape check is not decoration.  A DFPT run on the wrong cell still
    returns a square ``(n, n, 3, 3)`` force-constant array, so a Gamma-only
    result from the 2-atom cell would otherwise be recorded as a valid volume
    and silently lose the finite-q instability the recipe exists to find.

    Stability is *not* decided here.  It needs every volume at once (the
    decision is per endpoint, and a single soft volume fails the endpoint
    rather than being silently dropped), so it runs as its own job after
    ``collect_validated_volumes`` -- see ``classify_dfpt_stability``.
    """

    try:
        parsed = DfptTaskDoc.model_validate(dfpt_doc)
        if not parsed.force_constants:
            raise ValueError("empty force constants")
        if recipe is not None:
            expected = Recipe.from_dict(recipe)
            expected_atoms = expected.kpoints.cell_atoms * int(
                expected.kpoints.supercell_factor
            ) ** 3
            if parsed.force_constants_provenance.n_atoms != expected_atoms:
                raise ValueError(
                    f"DFPT ran on {parsed.force_constants_provenance.n_atoms} atoms "
                    f"but the recipe's phonon supercell is {expected_atoms} atoms "
                    f"(cell {expected.kpoints.cell_atoms} x factor "
                    f"{expected.kpoints.supercell_factor})"
                )
    except Exception as exc:
        return Response(output={"ok": False, "error": str(exc)}, stop_children=True)
    try:
        document = parsed.model_dump(mode="json")
    except Exception as exc:  # noqa: BLE001 - an unserializable document is a finding
        return Response(
            output={
                "ok": False,
                "error": f"validated DFPT document could not be serialized: {exc}",
            },
            stop_children=True,
        )
    return Response(
        output={
            "ok": True,
            "schema_version": parsed.schema_version,
            "n_atoms": parsed.force_constants_provenance.n_atoms,
            "dfpt_doc": document,
        }
    )


__all__ = [
    "DYNMAT_UNIT",
    "DfptRelaxMaker",
    "DfptVaspMaker",
    "collect_dfpt_task_doc",
    "force_constants_from_dynmat",
    "parse_dynmat",
    "validate_dfpt_task_doc",
]
