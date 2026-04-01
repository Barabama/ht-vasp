"""
OJ Input Set Generator

Generates VASP input files for OstravaJ magnetic exchange calculations.
Compatible with atomate2 VaspInputGenerator for seamless integration.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

from atomate2.vasp.sets.base import VaspInputGenerator

if TYPE_CHECKING:
    from pymatgen.core import Structure
    from pymatgen.io.vasp import Kpoints
    from pymatgen.io.vasp.sets import UserPotcarFunctional

log = logging.getLogger(__name__)


@dataclass
class OJInputSetGenerator(VaspInputGenerator):
    """Generator for OstravaJ magnetic exchange calculation input sets.

    This class inherits from VaspInputGenerator to ensure compatibility with
    atomate2 workflows, allowing seamless integration with other Makers like
    StaticMaker for combined workflows.

    OJ-specific parameters:
        j_count: Number of J pairs to consider (mutually exclusive with dist_cutoff)
        dist_cutoff: Distance cutoff for J pairs in Angstrom
        magnetic_ion_types: List of magnetic ion types (e.g., ["Fe", "Co"])
        noncollinear: Whether to use noncollinear magnetism
        base_spin: Base spin value
        extend_poscar: Supercell extension factors (nx, ny, nz)
        oj_conf_settings: OJ.conf configuration parameters

    Inherited parameters from VaspInputGenerator:
        See VaspInputGenerator documentation for full list.
        Key ones relevant to OJ:
        - user_incar_settings: Override INCAR settings
        - user_kpoints_settings: Override KPOINTS settings
        - reciprocal_density: K-point density (default: 100)
        - reciprocal_density_metal: K-point density for metals (default: 400)
        - auto_metal_kpoints: Auto-detect metals (default: True)
        - force_gamma: Force gamma-centered mesh (default: True)

    Example:
        >>> from pymatgen.core import Structure
        >>> from htvasp.oj import OJInputSetGenerator
        >>>
        >>> structure = Structure.from_file("POSCAR")
        >>> generator = OJInputSetGenerator(
        ...     j_count=3,
        ...     magnetic_ion_types=["Fe", "Co"],
        ...     reciprocal_density=100,
        ... )
        >>> input_set = generator.get_input_set(structure)
        >>> input_set.write_input("./flip000")
    """

    name: str = "oj_input"
    j_count: int | None = 2
    dist_cutoff: float | None = None
    magnetic_ion_types: list[str] = field(default_factory=list)
    noncollinear: bool = False
    base_spin: float | list[float] = 1.0
    extend_poscar: tuple[int, int, int] = (2, 2, 2)
    reciprocal_density: float = 100
    reciprocal_density_metal: float = 400
    auto_metal_kpoints: bool = True
    force_gamma: bool = True
    user_potcar_functional: UserPotcarFunctional = "PBE_64"
    user_incar_settings: dict[str, Any] = field(default_factory=dict)
    user_kpoints_settings: dict[str, Any] | Kpoints | None = None
    user_potcar_settings: dict[str, Any] = field(default_factory=dict)
    constrain_total_magmom: bool = False
    sort_structure: bool = False
    auto_ismear: bool = False
    auto_ispin: bool = False

    def __post_init__(self) -> None:
        """Initialize the input set generator."""
        super().__post_init__()

    @property
    def incar_updates(self) -> dict:
        """Get INCAR updates for OJ calculations.

        Returns:
            Dictionary of INCAR parameter updates.
        """
        return {
            "ENCUT": 500,
            "ISTART": 0,
            "ICHARG": 2,
            "ISMEAR": 1,
            "SIGMA": 0.1,
            "ALGO": "Normal",
            "NELM": 100,
            "NELMIN": 6,
            "NELMDL": -6,
            "IBRION": 2,
            "ISIF": 3,
            "NSW": 50,
            "POTIM": 0.2,
            "EDIFF": 1e-6,
            "EDIFFG": -0.01,
            "ISPIN": 2,
            "KPAR": 2,
            "NCORE": 2,
            "ISYM": 0,
            "LREAL": "Auto",
            "PREC": "Accurate",
            "SYMPREC": 1e-7,
            "LORBIT": 11,
            "LWAVE": False,
            "LCHARG": False,
        }

    def get_oj_conf_string(self) -> str:
        """Generate OJ.conf file content.

        Returns:
            OJ.conf configuration string.
        """
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

        lines.append(
            f"extend_poscar {self.extend_poscar[0]} {self.extend_poscar[1]} {self.extend_poscar[2]}"
        )

        return "\n".join(lines) + "\n"

    def get_input_set(
        self,
        structure: Structure | None = None,
        prev_dir: str | Path | None = None,
        potcar_spec: bool = False,
    ) -> "VaspInputSet":
        """Get the VASP input set.

        Args:
            structure: A structure to generate the input set for. If None, the structure
                must have been set previously.
            prev_dir: A previous calculation directory to copy output files from.
            potcar_spec: Instead of generating a POTCAR, use a list of POTCAR symbols.

        Returns:
            A VaspInputSet object.
        """
        # Call parent get_input_set
        input_set = super().get_input_set(
            structure=structure,
            prev_dir=prev_dir,
            potcar_spec=potcar_spec,
        )
        return input_set


def write_oj_input_set(
    structure: Structure,
    directory: Path | str,
    input_set_generator: OJInputSetGenerator,
    include_oj_conf: bool = True,
    **kwargs,
) -> None:
    """Write OJ input files to a directory.

    This function writes VASP input files using the OJInputSetGenerator,
    ensuring compatibility with atomate2 workflows.

    Args:
        structure: Input structure
        directory: Target directory
        input_set_generator: OJ input set generator
        include_oj_conf: Whether to include OJ.conf file
        **kwargs: Additional arguments passed to write_input
    """
    from pymatgen.io.vasp.inputs import VaspInput

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    vis = input_set_generator.get_input_set(structure)

    if include_oj_conf:
        oj_conf = input_set_generator.get_oj_conf_string()
        (directory / "OJ.conf").write_text(oj_conf)

    clean_prev = kwargs.pop("clean_prev", True)
    if clean_prev:
        for filename in ("POSCAR", "KPOINTS", "POTCAR", "INCAR"):
            filepath = directory / filename
            if filepath.exists():
                filepath.unlink()

    vis.write_input(directory, **kwargs)
    log.info(f"Written OJ input set to {directory}")
