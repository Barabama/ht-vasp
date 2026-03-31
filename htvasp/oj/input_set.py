"""
OJ Input Set

Generates input files for OstravaJ magnetic exchange calculations.
"""

import logging
from pathlib import Path
from typing import Literal

from pymatgen.core import Structure
from pymatgen.io.vasp import Kpoints, Potcar

from htvasp.oj.config import OJConfig

log = logging.getLogger(__name__)


def write_oj_inputs(
    structure: Structure,
    directory: Path | str,
    config: OJConfig,
    potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64",
    kppa: int = 1000,
) -> list[Path]:
    """Write OJ input files (POSCAR, INCAR, POTCAR, KPOINTS, OJ.conf)

    Args:
        structure: pymatgen Structure
        directory: Target directory
        config: OJ configuration
        potcar_functional: POTCAR functional type
        kppa: K-points per atom for automatic k-point mesh generation

    Returns:
        List of written file paths
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    written = []

    poscar_path = directory.joinpath("POSCAR")
    structure.to(filename=str(poscar_path), fmt="poscar")
    written.append(poscar_path)

    incar_path = directory.joinpath("INCAR")
    incar_path.write_text(config.to_incar_string())
    written.append(incar_path)

    kpoints_path = directory.joinpath("KPOINTS")
    try:
        kpoints = Kpoints.automatic_density(structure, kppa=kppa)
        kpoints.write_file(str(kpoints_path))
        written.append(kpoints_path)
        log.debug(f"Written KPOINTS with {kppa} k-points per atom")
    except Exception as e:
        log.warning(f"Failed to write KPOINTS: {e}")
        log.warning("KPOINTS generation skipped.")

    potcar_path = directory.joinpath("POTCAR")
    try:
        symbols = [str(el) for el in structure.species]
        potcar = Potcar(symbols=symbols, functional=potcar_functional)
        potcar.write_file(str(potcar_path))
        written.append(potcar_path)
        log.debug(f"Written POTCAR with functional {potcar_functional}")
    except Exception as e:
        log.warning(f"Failed to write POTCAR: {e}")
        log.warning("POTCAR generation skipped. Please ensure POTCAR is available.")

    oj_conf_path = directory.joinpath("OJ.conf")
    oj_conf_path.write_text(config.to_oj_conf_string())
    written.append(oj_conf_path)

    log.debug(f"Written {len(written)} input files to {directory}")
    return written
