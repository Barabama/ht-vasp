"""
OJ Task Document

Output schema for OstravaJ magnetic exchange calculations.
"""

from datetime import datetime
from typing import Any

from pydantic import BaseModel, Field, field_validator


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

    def to_summary(self) -> dict[str, Any]:
        """Return summary dict"""
        # Check for negative Tc_RPA (physical anomaly)
        tc_rpa = self.Tc_RPA
        has_negative_tc = False
        if tc_rpa < 0:
            has_negative_tc = True
        
        summary = {
            "J_pairs": dict(zip(self.J_reprs, self.Js)) if self.J_reprs else {},
            "Tc_MFA": self.Tc_MFA,
            "Tc_RPA": self.Tc_RPA,
            "num_configs": self.num_configs,
            "vasp_success_rate": self.vasp_success_rate,
        }
        
        if has_negative_tc:
            summary["warnings"] = ["Negative Tc_RPA detected (physical anomaly)"]
        
        return summary

    def __repr__(self) -> str:
        pairs = dict(zip(self.J_reprs, self.Js)) if self.J_reprs else {}
        return f"OJResult(J_pairs={pairs}, Tc_MFA={self.Tc_MFA}, Tc_RPA={self.Tc_RPA})"
