"""One VASP DFPT job and its dynamical-matrix collector."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dataclasses import dataclass, field

from atomate2.vasp.jobs.base import BaseVaspMaker
from atomate2.utils.path import strip_hostname
from jobflow import Response, job
from monty.io import zopen
from pymatgen.io.vasp.outputs import Vasprun

from htvasp.dfpt_qha.task_doc import DfptForceConstantsProvenance, DfptTaskDoc
from htvasp.dfpt_qha.input_set import DfptInputSetGenerator, DfptRelaxSetGenerator
from htvasp.model.dfpt_qha_recipe import Recipe


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
    stop_children_kwargs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.input_set_generator is None:
            if self.recipe is None:
                raise ValueError("DfptVaspMaker requires an explicit Recipe")
            self.input_set_generator = DfptInputSetGenerator(
                recipe=self.recipe, use_supercell_kpoints=self.use_supercell_kpoints
            )


@job
def collect_dfpt_task_doc(
    task_doc: Any,
    *,
    recipe: dict[str, Any],
    vasprun_path: str | Path | None = None,
) -> dict[str, Any]:
    """Convert one completed atomate2 TaskDoc into a DFPT TaskDoc.

    ``vasprun_path`` is optional because the normal BaseVaspMaker output keeps
    the calculation directory in ``TaskDoc.dir_name``.  The collector remains
    a pure jobflow analysis boundary: it parses ``vasprun.xml`` with pymatgen,
    never launches VASP and never scans an arbitrary calculation tree.
    """

    if hasattr(task_doc, "model_dump"):
        task_payload = task_doc.model_dump(mode="json")
        directory = task_doc.dir_name
    else:
        task_payload = dict(task_doc)
        directory = task_payload.get("dir_name")
    raw_directory = str(vasprun_path or directory)
    directory_path = Path(strip_hostname(raw_directory))
    plain = directory_path / "vasprun.xml"
    compressed = directory_path / "vasprun.xml.gz"
    path = plain if plain.is_file() else compressed if compressed.is_file() else None
    if path is None:
        raise ValueError(f"DFPT TaskDoc directory has no vasprun.xml(.gz): {directory_path}")
    vasprun = Vasprun(str(path), parse_potcar_file=False)
    force_constants = getattr(vasprun, "force_constants", None)
    if force_constants is None:
        raise ValueError(f"DFPT vasprun has no <dynmat> force constants: {path}")
    shape = tuple(int(x) for x in force_constants.shape)
    if len(shape) != 4 or shape[0] != shape[1] or shape[2:] != (3, 3):
        raise ValueError(f"unexpected DFPT force_constants shape: {shape}")
    doc = DfptTaskDoc(
        task_doc=task_payload,
        recipe=recipe,
        force_constants=force_constants.tolist(),
        force_constants_provenance=DfptForceConstantsProvenance(
            source="vasprun.xml",
            parser="pymatgen.Vasprun.force_constants",
            n_atoms=shape[0],
            shape=shape,
            vasprun_path=str(path),
        ),
        provenance={"task_doc_dir": str(directory), "task_state": task_payload.get("state")},
    )
    return doc.model_dump(mode="json")


@job
def validate_dfpt_task_doc(
    dfpt_doc: dict[str, Any], *, recipe: dict[str, Any] | None = None
) -> Response:
    """Stop descendants unless the DFPT result matches the recipe it claims.

    The shape check is not decoration.  A DFPT run on the wrong cell still
    returns a square ``(n, n, 3, 3)`` force-constant array, so a Gamma-only
    result from the 2-atom cell would otherwise be recorded as a valid volume
    and silently lose the finite-q instability the recipe exists to find.
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
    return Response(
        output={
            "ok": True,
            "schema_version": parsed.schema_version,
            "n_atoms": parsed.force_constants_provenance.n_atoms,
            "dfpt_doc": parsed.model_dump(mode="json"),
        }
    )


__all__ = [
    "DfptRelaxMaker",
    "DfptVaspMaker",
    "collect_dfpt_task_doc",
    "validate_dfpt_task_doc",
]
