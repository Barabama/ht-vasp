"""SQS structure-generation adapters.

The orchestration layer depends on a small, typed generator boundary.  The
optional ``sqsgenerator`` package is imported only when the concrete adapter
is used, so planning and validation remain usable without the engine.
"""

from __future__ import annotations

import importlib.metadata
from dataclasses import dataclass, field
from numbers import Real
from typing import Any, Mapping, Protocol

from pymatgen.core import Structure

from htvasp.model.sqs import SqsCompositionPoint
from htvasp.model.templates import structure_checksum


class SqsGenerationError(RuntimeError):
    """Raised when an SQS engine cannot produce a valid structure."""


@dataclass(frozen=True)
class GeneratedSqs:
    """An SQS structure and engine provenance returned by an adapter."""

    structure: Structure
    provenance: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.structure, Structure):
            raise TypeError("GeneratedSqs.structure must be a pymatgen Structure")
        object.__setattr__(self, "provenance", dict(self.provenance))


class SqsGenerator(Protocol):
    """Protocol implemented by concrete SQS generation engines."""

    engine_name: str

    def generate(
        self,
        template: Structure,
        point: SqsCompositionPoint,
        *,
        element_a: str,
        element_b: str,
        iterations: int,
        seed: int,
        engine_config: Mapping[str, Any] | None = None,
    ) -> GeneratedSqs:
        """Generate one structure for an integer composition point."""


class SqsgeneratorAdapter:
    """Strict adapter for the public ``sqsgenerator`` Python API.

    The input structure must already be the expanded supercell.  The adapter
    records the requested seed but does not claim reproducibility because the
    installed engine has no documented seed configuration key.
    """

    engine_name = "sqsgenerator"

    def generate(
        self,
        template: Structure,
        point: SqsCompositionPoint,
        *,
        element_a: str,
        element_b: str,
        iterations: int,
        seed: int,
        engine_config: Mapping[str, Any] | None = None,
    ) -> GeneratedSqs:
        if len(template) != point.n_sites:
            raise SqsGenerationError(
                f"Expanded template has {len(template)} sites, expected {point.n_sites}"
            )
        if iterations <= 0:
            raise ValueError("SQS iterations must be positive")

        try:
            from sqsgenerator import StructureFormat, optimize
        except ImportError as exc:  # pragma: no cover - depends on environment
            raise SqsGenerationError(
                "sqsgenerator is required for SQS generation; install it explicitly"
            ) from exc

        config: dict[str, Any] = dict(engine_config or {})
        config.update(
            {
                "iterations": int(iterations),
                "iteration_mode": config.get("iteration_mode", "random"),
                "sublattice_mode": config.get("sublattice_mode", "interact"),
                "composition": {
                    str(element_a): point.element_a_count,
                    str(element_b): point.element_b_count,
                },
                "structure": {
                    "lattice": template.lattice.matrix.tolist(),
                    "coords": template.frac_coords.tolist(),
                    "species": [str(element_a)] * len(template),
                },
            }
        )

        try:
            result_pack = optimize(config)
            best = result_pack.best()
            generated = Structure.from_str(
                best.structure().dump(StructureFormat.poscar), fmt="poscar"
            )
        except Exception as exc:
            raise SqsGenerationError(f"sqsgenerator failed for {point.label}: {exc}") from exc

        expected = {
            str(element_a): point.element_a_count,
            str(element_b): point.element_b_count,
        }
        actual = {
            str(element): int(amount)
            for element, amount in generated.composition.get_el_amt_dict().items()
        }
        if len(generated) != point.n_sites or actual != expected:
            raise SqsGenerationError(
                f"sqsgenerator returned composition {actual} for {point.label}; "
                f"expected {expected}"
            )

        objective = getattr(best, "objective", None)
        if not isinstance(objective, Real):
            objective = None
        try:
            version = importlib.metadata.version("sqsgenerator")
        except importlib.metadata.PackageNotFoundError:  # pragma: no cover
            version = "unknown"

        provenance = {
            "engine": self.engine_name,
            "engine_version": version,
            "iterations": int(iterations),
            "seed": int(seed),
            "seed_applied": False,
            "seed_note": (
                "sqsgenerator has no documented seed configuration key; "
                "seed is recorded for provenance only"
            ),
            "config": config,
            "objective": float(objective) if objective is not None else None,
            "structure_checksum": structure_checksum(generated),
        }
        return GeneratedSqs(generated, provenance)


__all__ = [
    "GeneratedSqs",
    "SqsGenerationError",
    "SqsGenerator",
    "SqsgeneratorAdapter",
]
