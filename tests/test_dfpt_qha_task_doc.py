"""No-VASP tests for the DFPT collector.

The collector reads VASP's ``DYNMAT`` -- not ``FORCE_CONSTANTS``.  ``IBRION=8``
writes only ``DYNMAT``; the verified reference directory has both only because
``FORCE_CONSTANTS`` was produced 90 minutes later by an offline conversion.
Assuming VASP writes it means the collector finds nothing after the allocation
has already been spent, which is exactly what the first real Slurm run showed.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from jobflow import JobStore
from maggma.stores import MemoryStore
from pymatgen.core import Structure

from htvasp.dfpt_qha.jobs import (
    collect_dfpt_task_doc,
    force_constants_from_dynmat,
    parse_dynmat,
    validate_dfpt_task_doc,
)
from htvasp.dfpt_qha.phonons import phonopy_for, supercell_structure
from htvasp.dfpt_qha.task_doc import DfptTaskDoc
from htvasp.model.dfpt_qha_recipe import Recipe

RECIPE_FILE = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json")
STRUCTURE_FILE = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_structures.json")
ENDPOINT = "BCC-Al-Ta-2"
MATRIX = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]

#: A real DFPT job directory, if one happens to be on disk.
LIVE_CANDIDATES = sorted(Path("/nfs_ssd/tmp").glob("BCC-Al-Ta-2-*/1-DFPT_volume_1/5"))


def _store() -> JobStore:
    docs = MemoryStore()
    docs.connect()
    return JobStore(docs)


def recipe_dict() -> dict:
    return json.loads(RECIPE_FILE.read_text())["endpoints"][ENDPOINT]["recipe"]


def cell() -> Structure:
    return Structure.from_dict(
        json.loads(STRUCTURE_FILE.read_text())["endpoints"][ENDPOINT]["structure"]
    )


def synthetic_dynmat(cell: Structure, supercell: Structure, fc: np.ndarray) -> str:
    """Write a DYNMAT in VASP's layout for a known force-constant matrix.

    VASP's file holds, per perturbation ``(atom, direction)``, the force on every
    atom produced by a UNIT displacement along that direction.  With unit
    displacement the response IS the force-constant column, so
    ``F_i = -sum_j Phi[i, j] . u_j`` with ``u`` a single unit vector.
    """

    phonon = phonopy_for(cell, MATRIX)
    phonon.generate_displacements(distance=1.0)
    lines = [f"{len(set(str(s.specie) for s in supercell))} {len(supercell)} "
             f"{len(phonon.displacements)}", " ".join(["1.0"] * len(supercell))]
    for displacement in phonon.displacements:
        atom, vector = int(displacement[0]), np.array(displacement[1:4], dtype=float)
        u = np.zeros((len(supercell), 3))
        u[atom] = vector
        forces = -np.einsum("ijkl,jl->ik", fc, u)
        lines.append(
            f"{atom + 1} {int(np.argmax(np.abs(vector))) + 1} "
            f"{vector[0]:.4f} {vector[1]:.4f} {vector[2]:.4f}"
        )
        lines.extend(
            f" {f[0]:.12f} {f[1]:.12f} {f[2]:.12f}" for f in forces
        )
    return "\n".join(lines) + "\n"


def test_parse_dynmat_reads_every_perturbation(tmp_path):
    text = "    1    3    2\n 1.0 1.0 1.0\n 1 1 1.0 0.0 0.0\n" + "".join(
        " 0.1 0.2 0.3\n" for _ in range(3)
    ) + " 2 1 1.0 0.0 0.0\n" + "".join(" 0.4 0.5 0.6\n" for _ in range(3))
    path = tmp_path / "DYNMAT"
    path.write_text(text)
    perturbations, n_atoms = parse_dynmat(path)
    assert n_atoms == 3
    assert [p["atom"] for p in perturbations] == [1, 2]
    assert np.allclose(perturbations[0]["forces"][0], [0.1, 0.2, 0.3])
    assert np.allclose(perturbations[1]["forces"][2], [0.4, 0.5, 0.6])


def test_parse_dynmat_rejects_a_truncated_file(tmp_path):
    path = tmp_path / "DYNMAT"
    path.write_text("    1    3    2\n 1.0 1.0 1.0\n 1 1 1.0 0.0 0.0\n 0.1 0.2 0.3\n")
    with pytest.raises(ValueError, match="force lines"):
        parse_dynmat(path)


def test_force_constants_replay_rejects_a_perturbation_mismatch():
    """The replay must refuse forces that do not belong to the perturbations.

    With only the symmetry-irreducible perturbations in the file, phonopy has
    to expand them with the crystal symmetry.  Attaching a force block to the
    wrong perturbation would still produce a square, non-zero matrix -- so the
    atom index and displacement vector are both checked.
    """

    base = cell()
    supercell = supercell_structure(base, MATRIX)
    n = len(supercell)
    reference = phonopy_for(base, MATRIX)
    reference.generate_displacements(distance=1.0)
    reference.forces = [np.zeros((n, 3)) for _ in reference.displacements]
    reference.produce_force_constants()
    text = synthetic_dynmat(base, supercell, np.array(reference.force_constants))

    perturbations = parse_dynmat_text(text)
    # a well-formed file with the right count is accepted
    recovered = force_constants_from_dynmat(
        base, supercell, perturbations, supercell_matrix=MATRIX
    )
    assert recovered.shape == (n, n, 3, 3)

    # displacing the wrong atom is refused
    tampered = [dict(p) for p in perturbations]
    tampered[0]["atom"] = perturbations[0]["atom"] + 1
    with pytest.raises(ValueError, match="phonopy expects atom"):
        force_constants_from_dynmat(base, supercell, tampered, supercell_matrix=MATRIX)

    # a wrong displacement vector is refused
    tampered = [dict(p) for p in perturbations]
    tampered[0]["displacement"] = [0.0, 0.0, 1.0]
    with pytest.raises(ValueError, match="not the same calculation"):
        force_constants_from_dynmat(base, supercell, tampered, supercell_matrix=MATRIX)

    # a missing perturbation is refused
    with pytest.raises(ValueError, match="phonopy expects"):
        force_constants_from_dynmat(
            base, supercell, perturbations[:1], supercell_matrix=MATRIX
        )


def parse_dynmat_text(text: str):
    """parse_dynmat for an in-memory file (keeps the test free of tmp plumbing)."""

    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix="DYNMAT", delete=False) as handle:
        handle.write(text)
        name = handle.name
    try:
        return parse_dynmat(name)[0]
    finally:
        Path(name).unlink()


def test_force_constants_reject_a_mismatched_perturbation_count():
    base = cell()
    supercell = supercell_structure(base, MATRIX)
    n = len(supercell)
    fc = np.zeros((n, n, 3, 3))
    text = synthetic_dynmat(base, supercell, fc)
    perturbations = parse_dynmat_text(text)[:1]
    with pytest.raises(ValueError, match="phonopy expects"):
        force_constants_from_dynmat(base, supercell, perturbations, supercell_matrix=MATRIX)


def test_collector_requires_the_dynmat_file(tmp_path):
    with pytest.raises(ValueError, match="DYNMAT"):
        collect_dfpt_task_doc(
            {"dir_name": str(tmp_path), "state": "successful", "structure": {}},
            recipe=recipe_dict(),
            cell=cell().as_dict(),
            supercell_matrix=MATRIX,
        ).run(_store())


def test_validation_stops_on_a_malformed_document():
    result = validate_dfpt_task_doc({"bad": True}, recipe=recipe_dict()).run(_store())
    assert result.output["ok"] is False
    assert result.stop_children is True


@pytest.mark.skipif(not LIVE_CANDIDATES, reason="no live DFPT job directory on disk")
def test_collector_on_a_real_dynmat_from_a_real_run():
    """End to end on VASP output that this project actually produced.

    The job ran on phonopy's own supercell, so the ordering assertion passes;
    that is the property the whole design turns on.
    """

    import gzip

    from pymatgen.io.vasp import Poscar

    job_dir = LIVE_CANDIDATES[0]
    with gzip.open(job_dir / "POSCAR.gz", "rt") as handle:
        supercell = Poscar.from_str(handle.read()).structure
    base = Structure(
        np.array(supercell.lattice.matrix) / 2,
        ["Al", "Ta"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    payload = collect_dfpt_task_doc(
        {
            "dir_name": str(job_dir),
            "state": "successful",
            "structure": supercell.as_dict(),
        },
        recipe=recipe_dict(),
        cell=base.as_dict(),
        supercell_matrix=MATRIX,
    ).run(_store()).output
    parsed = DfptTaskDoc.model_validate(payload)
    provenance = parsed.force_constants_provenance
    assert provenance.source == "DYNMAT"
    assert provenance.n_atoms == 16
    assert provenance.shape == (16, 16, 3, 3)
    assert provenance.ordering_fingerprint
    scale = max(abs(x) for row in parsed.force_constants for m in row for v in m for x in v)
    assert provenance.max_acoustic_sum_residual / scale < 1e-6


@pytest.mark.skipif(not LIVE_CANDIDATES, reason="no live DFPT job directory on disk")
def test_validation_accepts_the_real_collected_document():
    import gzip

    from pymatgen.io.vasp import Poscar

    job_dir = LIVE_CANDIDATES[0]
    with gzip.open(job_dir / "POSCAR.gz", "rt") as handle:
        supercell = Poscar.from_str(handle.read()).structure
    base = Structure(
        np.array(supercell.lattice.matrix) / 2,
        ["Al", "Ta"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    collected = collect_dfpt_task_doc(
        {"dir_name": str(job_dir), "state": "successful",
         "structure": supercell.as_dict()},
        recipe=recipe_dict(),
        cell=base.as_dict(),
        supercell_matrix=MATRIX,
    ).run(_store()).output
    result = validate_dfpt_task_doc(collected, recipe=recipe_dict()).run(_store())
    assert result.output["ok"] is True
    # the validated document travels with the verdict, so the thermal job can
    # consume it; a gate returning only {ok: bool} would strand the data
    assert result.output["dfpt_doc"]["force_constants"]


def test_validation_survives_a_task_doc_holding_pymatgen_objects(tmp_path):
    """jobflow stores output as JSON and monty rebuilds what it recognises.

    So the document the validation job receives can carry a real ``Kpoints``
    inside ``input.kpoints``.  The gate re-serializes its verdict, and without
    sanitizing at the collector boundary that fails with "Unable to serialize
    unknown type: Kpoints" -- three jobs downstream of the cause, which is
    exactly where it is hardest to read.
    """

    from pymatgen.io.vasp import Kpoints

    # A phonopy-ordered supercell with a zero dynamical matrix: the point here
    # is the Kpoints object surviving serialization, not the physics.
    base = cell()
    supercell = supercell_structure(base, MATRIX)
    fc = np.zeros((len(supercell), len(supercell), 3, 3))
    job_dir = tmp_path / "job"
    job_dir.mkdir()
    (job_dir / "DYNMAT").write_text(synthetic_dynmat(base, supercell, fc))

    task_doc = {
        "dir_name": str(job_dir),
        "state": "successful",
        "structure": supercell.as_dict(),
        # a live pymatgen object, as monty's decoder would hand it back
        "input": {"kpoints": Kpoints.monkhorst_automatic(kpts=(8, 8, 8), shift=(0.5, 0.5, 0.5))},
    }
    collected = collect_dfpt_task_doc(
        task_doc, recipe=recipe_dict(), cell=base.as_dict(), supercell_matrix=MATRIX
    ).run(_store()).output
    result = validate_dfpt_task_doc(collected, recipe=recipe_dict()).run(_store())
    assert result.output["ok"] is True, result.output
    # the verdict must itself be serializable, since jobflow will store it
    import json

    json.dumps(result.output)


def test_validation_accepts_a_task_doc_that_came_back_as_a_model(tmp_path):
    """The document crosses three job boundaries; each one may rebuild it.

    jobflow stores output as JSON and monty's decoder rebuilds recognised
    objects, so by the time the gate sees ``task_doc`` it can be a live emmet
    ``TaskDoc`` again.  A field typed ``dict[str, Any]`` rejects that outright:
    "Input should be a valid dictionary ... input_type=TaskDoc".  This is the
    second time the same boundary has bitten -- the first was a Kpoints nested
    inside it.
    """

    from atomate2.vasp.jobs.base import get_vasp_task_document

    # A phonopy-ordered supercell with a zero dynamical matrix: the point is the
    # TaskDoc surviving the boundary, not the physics.
    base = cell()
    supercell = supercell_structure(base, MATRIX)
    fc = np.zeros((len(supercell), len(supercell), 3, 3))
    job_dir = tmp_path / "job"
    job_dir.mkdir()
    (job_dir / "DYNMAT").write_text(synthetic_dynmat(base, supercell, fc))

    # a document as the collector would emit it
    payload = collect_dfpt_task_doc(
        {"dir_name": str(job_dir), "state": "successful",
         "structure": supercell.as_dict()},
        recipe=recipe_dict(),
        cell=base.as_dict(),
        supercell_matrix=MATRIX,
    ).run(_store()).output

    # ...and the same document after a job boundary rebuilt the TaskDoc.
    # get_vasp_task_document needs a real directory, so the rebuilt form is
    # emulated with the one monty would construct from the stored JSON.
    from emmet.core.tasks import TaskDoc

    rebuilt = dict(payload)
    rebuilt["task_doc"] = TaskDoc.model_validate(payload["task_doc"])
    assert not isinstance(rebuilt["task_doc"], dict)

    result = validate_dfpt_task_doc(rebuilt, recipe=recipe_dict()).run(_store())
    assert result.output["ok"] is True, result.output
    import json

    json.dumps(result.output)


def test_validation_rejects_a_matrix_from_the_wrong_cell():
    """A Gamma-only result from the 2-atom cell is square and non-empty.

    Without the atom-count check against the recipe it would be recorded as a
    valid volume and the finite-q instability would be lost silently.
    """

    document = {
        "schema_version": "htvasp-dfpt-task/v1",
        "task_doc": {},
        "recipe": recipe_dict(),
        "force_constants": [[[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]]],
        "force_constants_provenance": {
            "source": "DYNMAT",
            "parser": "phonopy replay at unit displacement",
            "n_atoms": 2,
            "shape": [2, 2, 3, 3],
            "ordering_fingerprint": "x" * 64,
        },
    }
    result = validate_dfpt_task_doc(document, recipe=recipe_dict()).run(_store())
    assert result.output["ok"] is False
    assert result.stop_children is True
    assert "phonon supercell" in result.output["error"]
