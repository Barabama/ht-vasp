"""OJ Input Set Generator - VASP input files for OstravaJ magnetic exchange calculations."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from atomate2.vasp.sets.base import VaspInputGenerator

if TYPE_CHECKING:
    from pymatgen.core import Structure
    from pymatgen.io.vasp import VaspInput

log = logging.getLogger(__name__)


@dataclass
class OJInputSetGenerator(VaspInputGenerator):
    """Generator for OstravaJ magnetic exchange calculation input sets.

    Inherits from VaspInputGenerator for atomate2 compatibility.

    Args:
        j_count: Number of J pairs to consider (mutually exclusive with dist_cutoff)
        dist_cutoff: Distance cutoff for J pairs in Angstrom
        magnetic_ion_types: List of magnetic ion types (e.g., ["Fe", "Co"])
        noncollinear: Whether to use noncollinear magnetism
        base_spin: Base spin value
        extend_poscar: Supercell extension factors (nx, ny, nz)
    """

    j_count: int | None = 4
    dist_cutoff: float | None = None
    magnetic_ion_types: list[str] = field(default_factory=list)
    noncollinear: bool = False
    base_spin: float | list[float] = 1.0
    extend_poscar: tuple[int, int, int] = (2, 2, 2)

    def __post_init__(self) -> None:
        if self.j_count is None and self.dist_cutoff is None:
            raise ValueError("Either j_count or dist_cutoff must be specified")
        super().__post_init__()

    @property
    def incar_updates(self) -> dict:
        """Get INCAR updates for OJ calculations."""
        return {
            "ENCUT": 500,
            "NELM": 100,
            "IBRION": 2,
            "ISIF": 2,
            "NSW": 50,
            "EDIFF": 1e-6,
            "EDIFFG": -0.02,
            "ISYM": 0,
            "LREAL": "Auto",
            "LORBIT": 10,
            "LWAVE": False,
            "LCHARG": False,
        }

    def get_oj_conf_string(self) -> str:
        """Generate OJ.conf file content."""
        lines = []

        if self.j_count is not None:
            lines.append(f"J_count {self.j_count}")
        elif self.dist_cutoff is not None:
            lines.append(f"dist_cutoff {self.dist_cutoff}")

        if self.magnetic_ion_types:
            lines.append(f"magnetic_ion_types {' '.join(self.magnetic_ion_types)}")

        lines.append("noncollinear 1" if self.noncollinear else "noncollinear 0")

        if isinstance(self.base_spin, list):
            lines.append(f"base_spin {' '.join(str(x) for x in self.base_spin)}")
        else:
            lines.append(f"base_spin {self.base_spin}")

        lines.append(f"extend_poscar {self.extend_poscar[0]} {self.extend_poscar[1]} {self.extend_poscar[2]}")

        return "\n".join(lines) + "\n"

    def get_input_set(
        self,
        structure: Structure | None = None,
        prev_dir: str | Path | None = None,
        potcar_spec: bool = False,
    ) -> VaspInput:
        """Get the VASP input set.

        OstravaJ requires INCAR without MAGMOM — it generates its own magnetic
        configurations. We strip any magmom site property and remove MAGMOM from INCAR.
        """
        if structure is not None and structure.site_properties.get("magmom"):
            structure = structure.copy()
            structure.remove_site_property("magmom")
            log.debug("Stripped 'magmom' site property for OJ input set")

        input_set = super().get_input_set(
            structure=structure,
            prev_dir=prev_dir,
            potcar_spec=potcar_spec,
        )

        if "MAGMOM" in input_set.incar:
            del input_set.incar["MAGMOM"]
            log.debug("Removed MAGMOM from INCAR for OJ input set")

        return input_set


def write_oj_input_set(
    structure: Structure,
    directory: Path | str,
    input_set_generator: OJInputSetGenerator,
    include_oj_conf: bool = True,
    **kwargs,
) -> None:
    """Write OJ input files to a directory.

    Args:
        structure: Input structure
        directory: Target directory
        input_set_generator: OJ input set generator
        include_oj_conf: Whether to include OJ.conf file
        **kwargs: Additional arguments passed to write_input
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    vis = input_set_generator.get_input_set(structure)

    if include_oj_conf:
        (directory / "OJ.conf").write_text(input_set_generator.get_oj_conf_string())

    clean_prev = kwargs.pop("clean_prev", True)
    if clean_prev:
        for filename in ("POSCAR", "KPOINTS", "POTCAR", "INCAR"):
            filepath = directory / filename
            if filepath.exists():
                filepath.unlink()

    vis.write_input(directory, **kwargs)
    log.debug(f"Written OJ input set to {directory}")
