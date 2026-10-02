"""
TB2J Task Document

Output schema for TB2J magnetic exchange calculations (J_ij from Wannier90).
"""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field, field_validator

from htvasp.tb2j.magnetic_contract import (
    MAGNETIC_CONTRACT_SCHEMA_VERSION,
    empty_magnetic_contract,
)

# Keys required in every entry of ``j_pairs``.
_J_PAIR_KEYS = ("i", "j", "R", "J_iso_meV", "vector", "distance")
_ORBITALS_PER_ATOM = 9  # s + p + d, matches input_set.ORBITALS_PER_ATOM

TB2J_SCHEMA_VERSION = "emtdb-magnetic-dataset/v1"
TB2J_J_UNIT = "meV"
TB2J_MOMENT_UNIT = "mu_B"
TB2J_LENGTH_UNIT = "angstrom"
TB2J_HAMILTONIAN = "H=-sum(J_ij S_i dot S_j)"
TB2J_PAIR_COUNTING = "outgoing_per_site"
TB2J_REVERSE_PAIR_POLICY = "allow_single_outgoing_validate_reverse_if_present"


class Tb2jResult(BaseModel):
    """Result from a TB2J magnetic exchange calculation.

    Contains the parsed exchange parameters (J_iso and full 3x3 J tensors in
    meV), provenance (elements, num_atoms, efermi, num_wann, kmesh) and the
    path to the raw ``exchange.out``.

    ``j_pairs`` entries have the shape::

        {
            "i": "Co1",            # atom label of site i
            "j": "Ni1",            # atom label of site j
            "R": [0, 0, 0],        # lattice vector (fractional)
            "J_iso_meV": 4.8243,   # isotropic exchange in meV
            "vector": [1.406, 1.406, 1.406],   # displacement vector (Angstrom)
            "distance": 2.435,                 # interatomic distance (Angstrom)
        }
    """

    elements: list[str] = Field(default_factory=list, description="Magnetic elements")
    num_atoms: int = Field(0, description="Number of atoms in the cell")
    structure: dict[str, Any] = Field(
        default_factory=dict,
        description="Pymatgen Structure.as_dict() from the W90 static POSCAR"
    )
    efermi: float | None = Field(None, description="Fermi energy (eV)")
    num_wann: int = Field(0, description="Number of Wannier orbitals")
    kmesh: tuple[int, int, int] = Field((9, 9, 9), description="TB2J kmesh")
    j_pairs: list[dict[str, Any]] = Field(
        default_factory=list, description="Exchange pairs (see module docstring)"
    )
    j_tensors: list[list[list[float]]] = Field(
        default_factory=list, description="3x3 exchange tensors (meV), parallel to j_pairs"
    )
    atom_moments: list[dict[str, Any]] = Field(
        default_factory=list,
        description="Per-site Wannier moments, including w_magmom in mu_B",
    )
    exchange_out_path: str = Field("", description="Path to the raw exchange.out")
    warnings: list[str] = Field(default_factory=list, description="Non-fatal warnings")
    schema_version: str = Field(TB2J_SCHEMA_VERSION, description="Canonical magnetic dataset schema")
    source: dict[str, Any] = Field(
        default_factory=dict,
        description="Raw exchange.out source identity and SHA-256",
    )
    units: dict[str, str] = Field(
        default_factory=lambda: {
            "J": TB2J_J_UNIT,
            "J_iso": TB2J_J_UNIT,
            "distance": TB2J_LENGTH_UNIT,
            "vector": TB2J_LENGTH_UNIT,
            "w_magmom": TB2J_MOMENT_UNIT,
        },
        description="Units for numeric magnetic quantities",
    )
    exchange_convention: dict[str, Any] = Field(
        default_factory=lambda: {
            "hamiltonian": TB2J_HAMILTONIAN,
            "j_positive_is_ferromagnetic": True,
            "pair_counting": TB2J_PAIR_COUNTING,
            "spin_normalization": "unit",
            "reverse_pair_policy": TB2J_REVERSE_PAIR_POLICY,
        },
        description="Explicit downstream magnetic convention",
    )
    provenance: dict[str, Any] = Field(
        default_factory=dict,
        description="Traceability and producer contract metadata",
    )
    audit: dict[str, Any] = Field(
        default_factory=dict,
        description="Parser and production-contract audit status",
    )
    magnetic_contract: dict[str, Any] = Field(
        default_factory=empty_magnetic_contract,
        description="Versioned htvasp magnetic contract",
    )
    magnetic_observables: dict[str, Any] = Field(
        default_factory=dict,
        description="Derived TC/TN/BMAGN and MFA exchange observables",
    )

    @field_validator("magnetic_contract")
    @classmethod
    def _check_magnetic_contract(cls, value: dict[str, Any]) -> dict[str, Any]:
        if value and value.get("schema_version") != MAGNETIC_CONTRACT_SCHEMA_VERSION:
            raise ValueError(
                "magnetic_contract.schema_version must be "
                f"{MAGNETIC_CONTRACT_SCHEMA_VERSION!r}"
            )
        return value

    @field_validator("magnetic_observables")
    @classmethod
    def _check_magnetic_observables(cls, value: dict[str, Any]) -> dict[str, Any]:
        if value:
            status = value.get("status")
            if status not in {"computed", "nonmagnetic", "not_applicable", "unsupported", "needs_input"}:
                raise ValueError(f"unsupported magnetic_observables.status: {status!r}")
            for key in ("TC", "TN", "BMAGN", "exchange_matrix", "eigenvalues", "conversion_trace"):
                if key not in value:
                    raise ValueError(f"magnetic_observables missing key: {key}")
        return value

    @field_validator("j_pairs")
    @classmethod
    def _check_j_pairs(cls, v: list[dict[str, Any]]) -> list[dict[str, Any]]:
        for pair in v:
            missing = [k for k in _J_PAIR_KEYS if k not in pair]
            if missing:
                raise ValueError(f"j_pairs entry missing keys {missing}: {pair}")
        return v

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Tb2jResult":
        """Build a Tb2jResult from the dict returned by the tb2j_solve job.

        Args:
            data: Output dict of ``htvasp.tb2j.jobs.tb2j_solve``.

        Returns:
            Tb2jResult instance.
        """
        num_atoms = data.get("num_atoms") or 0
        num_wann = data.get("num_wann") or (num_atoms * _ORBITALS_PER_ATOM)
        return cls(
            elements=list(data.get("elements") or []),
            num_atoms=num_atoms,
            structure=dict(data.get("structure") or {}),
            efermi=data.get("efermi"),
            num_wann=num_wann,
            kmesh=tuple(data.get("kmesh") or (9, 9, 9)),
            j_pairs=list(data.get("j_pairs") or []),
            j_tensors=list(data.get("j_tensors") or []),
            atom_moments=list(data.get("atom_moments") or []),
            exchange_out_path=data.get("exchange_out_path") or "",
            warnings=list(data.get("warnings") or []),
            schema_version=data.get("schema_version", TB2J_SCHEMA_VERSION),
            source=dict(data.get("source") or {}),
            units=dict(data.get("units") or {
                "J": TB2J_J_UNIT,
                "J_iso": TB2J_J_UNIT,
                "distance": TB2J_LENGTH_UNIT,
                "vector": TB2J_LENGTH_UNIT,
                "w_magmom": TB2J_MOMENT_UNIT,
            }),
            exchange_convention=dict(data.get("exchange_convention") or {
                "hamiltonian": TB2J_HAMILTONIAN,
                "j_positive_is_ferromagnetic": True,
                "pair_counting": TB2J_PAIR_COUNTING,
                "spin_normalization": "unit",
                "reverse_pair_policy": TB2J_REVERSE_PAIR_POLICY,
            }),
            provenance=dict(data.get("provenance") or {}),
            audit=dict(data.get("audit") or {}),
            magnetic_contract=dict(data.get("magnetic_contract") or empty_magnetic_contract()),
            magnetic_observables=dict(data.get("magnetic_observables") or {}),
        )

    @classmethod
    def from_solution(cls, solution: dict[str, Any]) -> "Tb2jResult":
        """Alias of :meth:`from_dict` for parity with OJResult.from_solution."""
        return cls.from_dict(solution)

    def _shell_stats(self) -> list[dict[str, Any]]:
        """Group j_pairs into coordination shells by rounded distance."""
        shells: dict[float, list[float]] = {}
        for pair in self.j_pairs:
            key = round(float(pair["distance"]), 3)
            shells.setdefault(key, []).append(float(pair["J_iso_meV"]))
        stats = []
        for dist in sorted(shells):
            js = shells[dist]
            stats.append(
                {
                    "distance": dist,
                    "count": len(js),
                    "mean_J_meV": round(sum(js) / len(js), 4),
                    "max_abs_J_meV": round(max(abs(j) for j in js), 4),
                }
            )
        return stats

    def to_summary(self) -> dict[str, Any]:
        """Return a concise summary with J shell statistics.

        Returns:
            Dict with provenance, pair/tensor counts, per-shell stats, the
            strongest coupling and any warnings.
        """
        j_iso_values = [float(p["J_iso_meV"]) for p in self.j_pairs]
        shells = self._shell_stats()
        summary: dict[str, Any] = {
            "elements": self.elements,
            "num_atoms": self.num_atoms,
            "efermi": self.efermi,
            "num_wann": self.num_wann,
            "kmesh": list(self.kmesh),
            "num_j_pairs": len(self.j_pairs),
            "num_j_tensors": len(self.j_tensors),
            "num_shells": len(shells),
            "shells": shells[:10],
            "exchange_out_path": self.exchange_out_path,
            "magnetic_observables": self.magnetic_observables,
        }
        if j_iso_values:
            max_abs_idx = max(range(len(j_iso_values)), key=lambda k: abs(j_iso_values[k]))
            strongest = dict(self.j_pairs[max_abs_idx])
            summary["max_abs_J_meV"] = round(abs(j_iso_values[max_abs_idx]), 4)
            summary["strongest_pair"] = strongest
        if self.warnings:
            summary["warnings"] = self.warnings
        return summary

    def __repr__(self) -> str:
        n = len(self.j_pairs)
        strongest = "n/a"
        if self.j_pairs:
            js = [abs(float(p["J_iso_meV"])) for p in self.j_pairs]
            strongest = f"{max(js):.3f} meV"
        return (
            f"Tb2jResult(elements={self.elements}, num_atoms={self.num_atoms}, "
            f"num_wann={self.num_wann}, j_pairs={n}, max_abs_J={strongest})"
        )
