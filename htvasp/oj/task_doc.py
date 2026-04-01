"""
OJ Task Document

Output schema for OstravaJ magnetic exchange calculations.
"""

from datetime import datetime
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

    A lightweight output model containing exchange parameters J and Curie temperatures.
    """

    J_reprs: list[Any] = Field(default_factory=list, description="J pair representations")
    Js: list[float] = Field(default_factory=list, description="Exchange parameters (meV)")
    Tc_MFA: float | None = Field(None, description="Curie temperature MFA (K)")
    Tc_RPA: float | None = Field(None, description="Curie temperature RPA (K)")
    num_configs: int = Field(0, description="Number of magnetic configurations")
    magnetic_ion_types: list[str] = Field(default_factory=list)
    vasp_success_rate: float = Field(1.0, description="VASP success rate")

    @field_validator("J_reprs", mode="before")
    @classmethod
    def ensure_list(cls, v):
        """Ensure J_reprs list is a list of Any."""
        return v or []

    @classmethod
    def from_solution(cls, solution: dict[str, Any]) -> "OJResult":
        """Create OJResult from solution dict"""
        return cls(
            J_reprs=solution.get("J_reprs", []),
            Js=solution.get("Js", []),
            Tc_MFA=solution.get("Tc_MFA"),
            Tc_RPA=solution.get("Tc_RPA"),
            num_configs=solution.get("num_configs", 0),
            magnetic_ion_types=solution.get("magnetic_ion_types", []),
            vasp_success_rate=solution.get("vasp_success_rate", 1.0),
        )

    def get_j_pairs_dict(self) -> dict[str, float]:
        """Get J pairs as a dictionary with formatted keys

        Returns:
            Dictionary mapping J representation to J value

        Example:
            >>> result = OJResult(J_reprs=[["Fe", "Fe", 1]], Js=[10.5])
            >>> result.get_j_pairs_dict()
            {'Fe-Fe-1': 10.5}
        """
        if not self.J_reprs or not self.Js:
            return {}

        formatted_reprs = [format_j_repr(j) for j in self.J_reprs]
        return dict(zip(formatted_reprs, self.Js))

    def to_summary(self) -> dict[str, Any]:
        """Return summary dict"""
        has_negative_tc = False
        if self.Tc_RPA is not None and self.Tc_RPA < 0:
            has_negative_tc = True

        summary = {
            "J_pairs": self.get_j_pairs_dict(),
            "Tc_MFA": self.Tc_MFA,
            "Tc_RPA": self.Tc_RPA,
            "num_configs": self.num_configs,
            "vasp_success_rate": self.vasp_success_rate,
        }

        if has_negative_tc:
            summary["warnings"] = ["Negative Tc_RPA detected (physical anomaly)"]

        return summary

    def __repr__(self) -> str:
        j_pairs = self.get_j_pairs_dict()
        return f"OJResult(J_pairs={j_pairs}, Tc_MFA={self.Tc_MFA}, Tc_RPA={self.Tc_RPA})"
