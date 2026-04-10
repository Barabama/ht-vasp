"""
OJ Task Document

Output schema for OstravaJ magnetic exchange calculations.
"""

from datetime import datetime
from functools import cached_property
from typing import Any

from pydantic import BaseModel, Field, field_validator


def format_j_repr(j_repr: list | str) -> str:
    """Format J representation as string

    Args:
        j_repr: J representation, either as list ["Fe", "Fe", 1] or string "Fe-Fe-1"

    Returns:
        Formatted string like "Fe-Fe-1"

    Example:
        >>> format_j_repr(["Fe", "Fe", 1])
        'Fe-Fe-1'
        >>> format_j_repr("Fe-Fe-1")
        'Fe-Fe-1'
    """
    if isinstance(j_repr, str):
        return j_repr
    elif isinstance(j_repr, list):
        return "-".join(str(x) for x in j_repr)
    else:
        return str(j_repr)


class OJResult(BaseModel):
    """Result from OstravaJ magnetic exchange calculation

    A comprehensive output model containing exchange parameters J, Curie temperatures,
    magnetic moment information, equation system diagnostics, and per-config detailed data.
    """

    # Original fields
    J_reprs: list[Any] = Field(default_factory=list, description="J pair representations")
    Js: list[float] = Field(default_factory=list, description="Exchange parameters (eV)")
    Tc_MFA: float | None = Field(None, description="Curie temperature MFA (K)")
    Tc_RPA: float | None = Field(None, description="Curie temperature RPA (K)")
    num_configs: int = Field(0, description="Number of magnetic configurations")
    magnetic_ion_types: list[str] = Field(default_factory=list)
    vasp_success_rate: float = Field(1.0, description="VASP success rate")

    # New fields
    E_DLM: float | None = Field(None, description="DLM energy in eV/atom")
    D_stiff: float | None = Field(None, description="Spin-wave stiffness in eV*Ang^2")
    is_complete: bool = Field(True, description="Whether solution is complete")
    rank: int = Field(0, description="Matrix rank")
    resid: float = Field(0.0, description="Fitting residual")
    condition_number: float = Field(0.0, description="Equation system condition number")
    singular_values: list[float] = Field(default_factory=list, description="Singular values of equation system")

    # megnetic information
    avg_magnetic_moment: float = Field(0.0, description="Average absolute moment per magnetic ion in mu_B")
    total_magnetic_moment: float = Field(0.0, description="Total magnetic moment in mu_B")
    ground_state_moments: list[Any] = Field(default_factory=list, description="Per-atom moments of ground state")
    atom_types: list[str] = Field(default_factory=list, description="Atom type labels")

    # data of each config
    configs_data: list[dict[str, Any]] = Field(default_factory=list, description="Per-config detailed data")
    energies: list[float] = Field(default_factory=list, description="Total energies for all configs")
    J_mat: list[list[float]] = Field(default_factory=list, description="J coefficient matrix")

    @field_validator("J_reprs", mode="before")
    @classmethod
    def ensure_list(cls, v):
        """Ensure J_reprs list is a list of Any."""
        return v or []

    @classmethod
    def from_solution(cls, solution: dict[str, Any]) -> "OJResult":
        """Create OJResult from solution dict (步骤10)"""
        return cls(
            # Original fields
            J_reprs=solution.get("J_reprs", []),
            Js=solution.get("Js", []),
            Tc_MFA=solution.get("Tc_MFA"),
            Tc_RPA=solution.get("Tc_RPA"),
            num_configs=solution.get("num_configs", 0),
            magnetic_ion_types=solution.get("magnetic_ion_types", []),
            vasp_success_rate=solution.get("vasp_success_rate", 1.0),
            # New fields
            E_DLM=solution.get("E_DLM"),
            D_stiff=solution.get("D_stiff"),
            is_complete=solution.get("is_complete", True),
            rank=solution.get("rank", 0),
            resid=solution.get("resid", 0.0),
            condition_number=solution.get("condition_number", 0.0),
            singular_values=solution.get("singular_values", []),
            avg_magnetic_moment=solution.get("avg_magnetic_moment", 0.0),
            total_magnetic_moment=solution.get("total_magnetic_moment", 0.0),
            ground_state_moments=solution.get("ground_state_moments", []),
            atom_types=solution.get("atom_types", []),
            configs_data=solution.get("configs_data", []),
            energies=solution.get("energies", []),
            J_mat=solution.get("J_mat", []),
        )

    @cached_property
    def j_pairs_dict(self) -> dict[str, float]:
        """Get J pairs as a dictionary with formatted keys (cached).

        Returns:
            Dictionary mapping J representation to J value

        Example:
            >>> result = OJResult(J_reprs=[["Fe", "Fe", 1]], Js=[10.5])
            >>> result.j_pairs_dict
            {'Fe-Fe-1': 10.5}
        """
        if not self.J_reprs or not self.Js:
            return {}

        formatted_reprs = [format_j_repr(j) for j in self.J_reprs]
        return dict(zip(formatted_reprs, self.Js))

    def get_j_pairs_dict(self) -> dict[str, float]:
        """Get J pairs as a dictionary with formatted keys

        Returns:
            Dictionary mapping J representation to J value

        Example:
            >>> result = OJResult(J_reprs=[["Fe", "Fe", 1]], Js=[10.5])
            >>> result.get_j_pairs_dict()
            {'Fe-Fe-1': 10.5}
        """
        return self.j_pairs_dict

    def to_summary(self) -> dict[str, Any]:
        """Return summary dict"""
        summary = {
            "J_pairs": self.get_j_pairs_dict(),
            "Tc_MFA": self.Tc_MFA,
            "Tc_RPA": self.Tc_RPA,
            "E_DLM": self.E_DLM,
            "avg_magnetic_moment": self.avg_magnetic_moment,
            "total_magnetic_moment": self.total_magnetic_moment,
            "D_stiff": self.D_stiff,
            "condition_number": self.condition_number,
            "is_complete": self.is_complete,
            "num_configs": self.num_configs,
            "vasp_success_rate": self.vasp_success_rate,
        }

        warnings = []
        if self.Tc_RPA is not None and self.Tc_RPA < 0:
            warnings.append("Negative Tc_RPA detected")
        if self.condition_number > 1e10:
            warnings.append("High condition number - solution may be unreliable")
        if not self.is_complete:
            warnings.append("Solution is NOT complete")
        if warnings:
            summary["warnings"] = warnings

        return summary

    def __repr__(self) -> str:
        j_pairs = self.get_j_pairs_dict()
        return f"OJResult(J_pairs={j_pairs}, Tc_MFA={self.Tc_MFA}, Tc_RPA={self.Tc_RPA})"
