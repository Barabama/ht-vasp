"""
OJ Configuration

Configuration for OstravaJ magnetic exchange calculations using pydantic.
"""

from typing import Any

from pydantic import BaseModel, Field, field_validator, model_validator


class OJConfig(BaseModel):
    """OJ configuration for magnetic exchange calculations

    Attributes:
        j_count: Number of J pairs to consider (mutually exclusive with dist_cutoff)
        dist_cutoff: Distance cutoff for J pairs in Angstrom
        magnetic_ion_types: List of magnetic ion types (e.g., ["Fe", "Co"])
        noncollinear: Whether to use noncollinear magnetism
        base_spin: Base spin value
        extend_poscar: Supercell extension factors (nx, ny, nz)
        reciprocal_density: K-point density by reciprocal volume (same as atomate2)
        reciprocal_density_metal: K-point density for metallic systems
        auto_metal_kpoints: Automatically use higher density for metallic systems
        force_gamma: Force gamma-centered k-point mesh
        incar: INCAR settings for VASP calculations
    """

    j_count: int | None = 2
    dist_cutoff: float | None = None
    magnetic_ion_types: list[str] = Field(default_factory=list)
    noncollinear: bool = False
    base_spin: float | list[float] = 1.0
    extend_poscar: tuple[int, int, int] = (2, 2, 2)
    reciprocal_density: float = 100
    reciprocal_density_metal: float = 400
    auto_metal_kpoints: bool = True
    force_gamma: bool = True

    incar: dict[str, Any] = Field(
        default_factory=lambda: {
            "ENCUT": 500,
            "ISTART": 0,
            "ICHARG": 2,
            # Electronic
            "ISMEAR": 1,
            "SIGMA": 0.1,
            "ALGO": "Normal",
            "NELM": 100,
            "NELMIN": 6,
            "NELMDL": -6,
            # Ionic
            "IBRION": 2,
            "ISIF": 3,
            "NSW": 50,
            "POTIM": 0.2,
            "EDIFF": 1e-6,
            "EDIFFG": -0.01,
            # Magnetic
            "ISPIN": 2,
            # Precision
            "KPAR": 2,
            "NCORE": 2,
            "ISYM": 0,
            "LREAL": "Auto",
            "PREC": "Accurate",
            "SYMPREC": 1e-7,
            # Output
            "LORBIT": 11,
            "LWAVE": False,
            "LCHARG": False,
        }
    )

    @model_validator(mode="after")
    def validate_j_params(self):
        if self.j_count is None and self.dist_cutoff is None:
            raise ValueError("Either j_count or dist_cutoff must be specified")
        return self

    @field_validator("extend_poscar", mode="before")
    @classmethod
    def validate_extend_poscar(cls, v):
        if isinstance(v, list):
            return tuple(v)
        return v

    def to_oj_conf_string(self) -> str:
        """Generate OJ.conf file content"""
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

    def to_incar_string(self) -> str:
        """Generate INCAR file content"""
        lines = []
        for key, value in self.incar.items():
            if isinstance(value, bool):
                lines.append(f"{key} = .{str(value).upper()}.")
            else:
                lines.append(f"{key} = {value}")
        return "\n".join(lines) + "\n"

    def update_incar(self, **kwargs) -> "OJConfig":
        """Return a new instance with updated INCAR settings

        Supports:
        - Updating existing parameters
        - Adding new parameters
        - Removing parameters (set value to None)

        Args:
            **kwargs: INCAR parameters to update

        Returns:
            New OJConfig instance with updated INCAR

        Example:
            >>> config = OJConfig()
            >>> config2 = config.update_incar(ENCUT=520, ISPIN=None)  # Update ENCUT, remove ISPIN
        """
        new_incar = {**self.incar}

        for key, value in kwargs.items():
            if value is None:
                new_incar.pop(key, None)
            else:
                new_incar[key] = value

        return OJConfig(**{**self.model_dump(), "incar": new_incar})

    def merge_incar(self, incar_dict: dict[str, Any]) -> "OJConfig":
        """Merge a dictionary of INCAR settings

        Args:
            incar_dict: Dictionary of INCAR settings to merge

        Returns:
            New OJConfig instance with merged INCAR
        """
        return self.update_incar(**incar_dict)
